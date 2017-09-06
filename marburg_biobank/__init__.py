import zipfile
import os
import numpy as np
import pandas as pd
try:
    import cPickle as pickle
except ImportError:
    import pickle



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

    def get_dataset_compartments(self, dataset):
        """Get available compartments in dataset @dataset"""
        pcd = self.get_dataset('_meta/patient_compartment_dataset')
        pcd = pcd[pcd.dataset == dataset]
        return pcd['compartment'].unique()

    def get_variables_and_units(self, dataset):
        """What variables are availabe in a dataset?"""
        df = self.get_dataset(dataset)
        return df[['variable','unit']].groupby(['variable','unit']).groups.keys()

    def get_possible_values(self, variable, unit):
        pass

    @lazy_member('_cache_list_datasets')
    def list_datasets(self):
        """What datasets to we have"""
        return sorted([ name for name in self.zf.namelist() if
                not name.startswith('_') and not os.path.basename(name).startswith('_')])

    @lazy_member('_cache_list_datasets_incl_meta')
    def list_datasets_including_meta(self):
        """What datasets to we have"""
        return sorted(self.zf.namelist())

    @lazy_member('_datasets_with_name_lookup')
    def datasets_with_name_lookup(self):
        return  [ds for (ds, df) in self.iter_datasets() if 'name' in df.columns]

    def name_lookup(self, dataset, variable):
        df = self.get_dataset(dataset)
        return df[df.variable == variable]['name'].iloc[0]  # todo: optimize using where?

    def get_excluded_patients(self, dataset):
        global_exclusion_df = self.get_dataset('clinical/_other_exclusion')
        excluded = set(global_exclusion_df['patient'].unique())
        #local exclusion from this dataset
        try:
            exclusion_df = self.get_dataset(os.path.dirname(dataset) + '/' + '_' + os.path.basename(dataset) + '_exclusion')
            excluded.update(exclusion_df['patient'].unique())
        except KeyError:
            pass
        return excluded



    def iter_datasets(self, yield_meta=False):
        if yield_meta:
            l = self.list_datasets_including_meta()
        else:
            l = self.list_datasets()
        for name in l:
                yield name, self.get_dataset(name)

    def get_dataset(self, name):
        """Retrieve a dataset"""
        if name not in self.list_datasets_including_meta():
            raise KeyError("No such dataset: %s.\nAvailable: %s" %
                           (name, self.list_datasets_including_meta()))
        else:
            if name not in self._cached_datasets:
                self._cached_datasets[name] = self._get_dataset(name)
        return self._cached_datasets[name].copy()

    def _get_dataset(self, name):
        fh = self.zf.open(name)
        try:
            return pickle.load(fh, encoding='latin1')
        except TypeError:  # older pickle.load has no encoding, only needed in python3 anyhow
            return pickle.load(fh)
