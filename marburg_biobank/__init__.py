import zipfile
import os
import numpy as np
import pandas as pd
try:
    from functools import lru_cache
except (ImportError, AttributeError):
    from functools32 import lru_cache
try:
    import cPickle as pickle
except ImportError:
    import pickle

datasets_to_cache = 32


def lazy_member(field):
    """Evaluate a function once and store the result in the member (an object specific in-memory cache)
    Beware of using the same name in subclasses!
    """
    def decorate(func):
        if field == func.__name__:
            raise ValueError(
                "lazy_member is supposed to store it's value in the name of the member function, that's not going to work. Please choose another name (prepend an underscore...")

        def doTheThing(*args, **kw):
            if not hasattr(args[0], field):
                setattr(args[0], field, func(*args, **kw))
            return getattr(args[0], field)
        return doTheThing
    return decorate


class OvcaBiobank(object):
    """An interface to a dump of our Biobank.
    Also used internally by the biobank website to access the data.

    In essence, a souped up dict of pandas dataframes stored
    as pickles in a zip file with memory caching"""

    def __init__(self, filename):
        self.filename = filename
        self.zf = zipfile.ZipFile(filename)
        self._cached_datasets = {}

    def get_all_patients(self):
        df = self.get_dataset('_meta/patient_compartment_dataset')
        return set(df['patient'].unique())

    def number_of_patients(self):
        """How many patients/indivuums are in all datasets?"""
        return len(self.get_all_patients())

    def number_of_datasets(self):
        """How many different datasets do we have"""
        return len(self.list_datasets())

    def get_compartments(self):
        """Get all compartments we have data for"""
        pcd = self.get_dataset('_meta/patient_compartment_dataset')
        return pcd['compartment'].unique()

    @lru_cache(datasets_to_cache)
    def get_dataset_compartments(self, dataset):
        """Get available compartments in dataset @dataset"""
        pcd = self.get_dataset('_meta/patient_compartment_dataset')
        pcd = pcd[pcd.dataset == dataset]
        return pcd['compartment'].unique()

    @lru_cache(datasets_to_cache)
    def get_variables_and_units(self, dataset):
        """What variables are availabe in a dataset?"""
        df = self.get_dataset(dataset)
        if len(df['unit'].cat.categories) == 1:
            vars = df['variable'].unique()
            unit = df['unit'].iloc[0]
            return set([(v, unit) for v in vars])
        else:
            x = df[['variable', 'unit']].drop_duplicates(['variable', 'unit'])
            return set(zip(x['variable'], x['unit']))

    def get_possible_values(self, variable, unit):
        pass

    @lazy_member('_cache_list_datasets')
    def list_datasets(self):
        """What datasets to we have"""
        return sorted([name for name in self.zf.namelist() if
                       not name.startswith('_') and not os.path.basename(name).startswith('_')])

    @lazy_member('_cache_list_datasets_incl_meta')
    def list_datasets_including_meta(self):
        """What datasets to we have"""
        return sorted(self.zf.namelist())

    @lazy_member('_datasets_with_name_lookup')
    def datasets_with_name_lookup(self):
        return [ds for (ds, df) in self.iter_datasets() if 'name' in df.columns]

    def name_lookup(self, dataset, variable):
        df = self.get_dataset(dataset)
        # todo: optimize using where?
        return df[df.variable == variable]['name'].iloc[0]

    @lru_cache(maxsize=datasets_to_cache)
    def get_wide(self, dataset, standardized=False):
        """Return dataset in row=variable, column=patient format.
        if @standardized is True Index is always (variable, unit) or (variable, unit, name), and columns always (patient, compartment)
        Otherwise, unit and compartment will be left of if there is only a single value for them in the dataset"""
        df = self.get_dataset(dataset)
        columns = ['patient']
        index = ['variable']
        if standardized or len(df.compartment.cat.categories) > 1:
            columns.append('compartment')
        if standardized or len(df.unit.cat.categories) > 1:
            index.append('unit')
        if 'name' in df.columns:
            index.append('name')
        return self.to_wide(df, index, columns)

    def to_wide(self, df, index=['variable', ], columns=['patient', 'compartment'], values='value', sort_on_first_level=False):
        """Convert a dataset (or filtered dataset) to a wide DataFrame.
        Preferred to pd.pivot_table manually because it is
           a) faster and
           b) avoids a bunch of pitfalls when working with categorical data and
           c) makes sure the columns are dtype=float if they contain nothing but floats

        index = variable,unit
        columns = (patient, compartment)
        """
        df = df[['value'] + index + columns]
        set_index_on = index + columns
        columns_pos = tuple(range(len(index), len(index) + len(columns)))
        res = df.set_index(set_index_on).unstack(columns_pos)
        c = res.columns
        c = c.droplevel(0)
        # this removes categories from the levels of the index. Absolutly
        # necessar.
        if isinstance(c, pd.MultiIndex):
            c = pd.MultiIndex([list(x) for x in c.levels],
                              labels=c.labels, names=c.names)
        else:
            c = list(c)
        res.columns = c
        if sort_on_first_level:
            # sort on first level - ie. patient, not compartment - slow though
            res = res[sorted(list(res.columns))]
        for c in res.columns:
            try:
                res[c] = res[c].astype(float)
            except ValueError:
                pass
        return res

    @lru_cache(maxsize=datasets_to_cache)
    def get_excluded_patients(self, dataset):
        """Which patients are excluded from this particular dataset (or globally?"""
        global_exclusion_df = self.get_dataset('clinical/_other_exclusion')
        excluded = set(global_exclusion_df['patient'].unique())
        # local exclusion from this dataset
        try:
            exclusion_df = self.get_dataset(os.path.dirname(
                dataset) + '/' + '_' + os.path.basename(dataset) + '_exclusion')
            excluded.update(exclusion_df['patient'].unique())
        except KeyError:
            pass
        return excluded

    @lru_cache(maxsize=1)
    def get_exclusion_reasons(self):
        """Get exclusion information for all the datasets + globally"""
        result = {}
        global_exclusion_df = self.get_dataset('clinical/_other_exclusion')
        for tup in global_exclusion_df.itertuples():
            if tup.patient not in result:
                result[tup.patient] = {}
            result[tup.patient]['global'] = tup.reason
        for dataset in self.list_datasets():
            try:
                exclusion_df = self.get_dataset(os.path.dirname(
                    dataset) + '/' + '_' + os.path.basename(dataset) + '_exclusion')
                for tup in exclusion_df.itertuples():
                    if tup.patient not in result:
                        result[tup.patient] = {}
                    result[tup.patient][dataset] = tup.reason
            except KeyError:
                pass
        return result

    def iter_datasets(self, yield_meta=False):
        if yield_meta:
            l = self.list_datasets_including_meta()
        else:
            l = self.list_datasets()
        for name in l:
            yield name, self.get_dataset(name)

    @lru_cache(datasets_to_cache)
    def get_dataset(self, name):
        """Retrieve a dataset"""
        if name not in self.list_datasets_including_meta():
            raise KeyError("No such dataset: %s.\nAvailable: %s" %
                           (name, self.list_datasets_including_meta()))
        else:
            fh = self.zf.open(name)
            try:
                return pickle.load(fh, encoding='latin1')
            except TypeError:  # older pickle.load has no encoding, only needed in python3 anyhow
                return pickle.load(fh)
