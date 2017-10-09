import zipfile
import os
import numpy as np
import pandas as pd
try:
    from functools import lru_cache
except (ImportError, AttributeError):
    # don't know how to tell setup.py that we only need functools32 when under 2.7.
    # so we'll just include a copy (*bergh*)
    import sys
    sys.path.append(os.path.join(os.path.dirname(__file__), 'functools32'))
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
    def get_wide(self, dataset, apply_exclusion=True, standardized=False):
        """Return dataset in row=variable, column=patient format.
        if @standardized is True Index is always (variable, unit) or (variable, unit, name), and columns always (patient, compartment)
        Otherwise, unit and compartment will be left of if there is only a single value for them in the dataset
         if @apply_exclusion is True, excluded patients will be filtered from DataFrame

        """
        df = self.get_dataset(dataset)
        columns = ['patient']
        index = ['variable']
        if standardized or len(df.compartment.cat.categories) > 1:
            columns.append('compartment')
        if standardized or len(df.unit.cat.categories) > 1:
            index.append('unit')
        if 'name' in df.columns:
            index.append('name')
        dfw = self.to_wide(df, index, columns)
        if apply_exclusion:
            return self.apply_exclusion(dataset, dfw)
        else:
            return dfw

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
        # necessary, or you can't add columns later otherwise
        if isinstance(c, pd.MultiIndex):
            c = pd.MultiIndex([list(x) for x in c.levels],
                              labels=c.labels, names=c.names)
        else:
            c = list(c)
        res.columns = c
        if isinstance(c, list):
            res.columns.names = columns
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
        """Which patients are excluded from this particular dataset (or globally)?.

        May return a set of patient_id, or (patient_id, compartment) tuples if only
        certain compartments where excluded.

        """
        global_exclusion_df = self.get_dataset('clinical/_other_exclusion')
        excluded = set(global_exclusion_df['patient'].unique())
        # local exclusion from this dataset
        try:
            exclusion_df = self.get_dataset(os.path.dirname(
                dataset) + '/' + '_' + os.path.basename(dataset) + '_exclusion')
            if 'compartment' in exclusion_df.columns:
                excluded.update(zip(exclusion_df['patient'], exclusion_df['compartment']))
            else:
                excluded.update(exclusion_df['patient'].unique())
        except KeyError:
            pass
        return excluded

    def apply_exclusion(self, dataset_name, df):
        if not dataset_name in self.list_datasets():
            raise KeyError(dataset_name)
        excluded = self.get_excluded_patients(dataset_name)
        if df.columns.names == ('patient', 'compartment'):
            to_remove = [x for x in df.columns if x in excluded]
            return df.drop(to_remove, axis=1)
        if df.columns.names == ('patient', ):
            #single compartment dataset, ignore compartment in excluded
            excluded = set([x[0] if isinstance(x, tuple) else x for x in excluded])
            to_remove = [x for x in df.columns if x in excluded]
            return df.drop(to_remove, axis=1)
        elif 'patient' in df.columns and 'compartment' in df.columns:
            keep = np.ones((len(df),), np.bool)

            for x in excluded:
                if isinstance(x, tuple):
                    keep = keep & ~((df['patient'] == x[0]) & (df['compartment'] == x[1]))
                else:
                    keep = keep & ~(df['patient'] == x)
            return df[keep]
        else:
            raise ValueError("Sorry, not a tall or wide DataFrame that I know how to handle.")

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
    def get_dataset(self, name, apply_exclusion=False):
        """Retrieve a dataset"""
        if name not in self.list_datasets_including_meta():
            raise KeyError("No such dataset: %s.\nAvailable: %s" %
                           (name, self.list_datasets_including_meta()))
        else:
            with self.zf.open(name) as op:
                try:
                    df = pd.read_msgpack(op.read())
                    if apply_exclusion:
                        df = self.apply_exclusion(name, df)
                    return df
                except KeyError as e:
                    if "KeyError: u'category'" in str(e):
                        raise ValueError("Your pandas is too old. You need at least version 0.18")

    def get_comment(self, name):
        comments = self.get_dataset('_meta/comments')
        match = comments.path == name
        if match.any():
            return comments[match].iloc[0]['comment']
        else:
            return ''
