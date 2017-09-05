import zipfile
import pandas as pd
try:
    import cPickle as pickle
except ImportError:
    import pickle


class OvcaBiobank(object):
    """An interface to a dump of our Biobank.
    Also used internally by the biobank website to access the data.

    In essence, a souped up dict of pandas dataframes stored
    as pickles in a zip file with memory caching"""

    def __init__(self, filename):
        self.filename = filename
        self.zf = zipfile.open(filename)
        self._cached_datasets = {}

    def number_of_patients(self):
        """How many patients/indivuums are in all datasets?"""
        df = self.get_dataset('_meta/patient_compartment_dataset')
        return set(df['patient'].unique())

    def number_of_datasets(self):
        """How many different datasets do we have"""
        return len(self.list_datasets())

    @lazy_member('_chache_list_datasets')
    def list_datasets(self):
        """What datasets to we have"""
        return self.zf.namelist()

    def get_dataset(self, name):
        """Retrieve a dataset"""
        if name not in self.list_datasets():
            raise KeyError("No such dataset: %s" % name)
        else:
            if name not in self._cached_datasets:
                self._cached_datasets[name] = self._get_dataset(name)
        return self._cached_datasets[name].copy()

    def _get_dataset(self, name):
        fh = self.zipfile.open(name)
        return pickle.load(fh)


def create_biobank(
        dict_of_dataframes, name, revision, filename):
    """Create a file suitable for biobank consumption.
    Assumes all dataframes have at least variable, unit, patient, compartment and value columns
    """
    dict_of_dataframes['meta/biobank'] = pd.DataFrame([
        {'variable': 'biobank', 'value': name, },
        {'variable': 'revision', 'value': revision, },
    ])
    patient_compartment_dataset = {
        'patient': [], 'compartment': [], 'dataset': []}
    for name, df in dict_of_dataframes:
        if not name.startswith('meta'):
            must_have = {'variable', 'unit', 'patient', 'compartment'}
            missing = must_have.difference(df.columns)
            if missing:
                raise ValueError("%s is missing columns: %s" % missing)
            here = set(
                [df[['patient', 'compartment']].itertuples(index=False, name=None)])
            for p, c in here:
                patient_compartment_dataset['patient'].append(p)
                patient_compartment_dataset['compartment'].append(c)
                patient_compartment_dataset['dataset'].append(name)
    dict_of_dataframes[
        'meta/patient_compartment_dataset'] = pd.DataFrame(patient_compartment_dataset)
    with zipfile.open(filename, 'w') as op:
        for name, df in dict_of_dataframes:
            op.writestr(name, pickle.dumps(df, pickle.HIGHEST_PROTOCOL))


def lazy_member(field):
    """Evaluate a function once and store the result in the member (an object specific in-memory cache)
    Beware of using the same name in subclasses!
    """
    def decorate(func):
        if field == func.func_name:
            raise ValueError(
                "lazy_member is supposed to store it's value in the name of the member function, that's not going to work. Please choose another name (prepend an underscore...")

        def doTheThing(*args, **kw):
            if not hasattr(args[0], field):
                setattr(args[0], field, func(*args, **kw))
            return getattr(args[0], field)
        return doTheThing
    return decorate
