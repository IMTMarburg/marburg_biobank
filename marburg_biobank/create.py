import pandas as pd
import re
import pickle
import zipfile
import os


must_have_columns = ['variable', 'unit', 'patient',  'compartment', 'value', ]
allowed_compartments = {'n.a.', 'ascites',
                        'TAM', 'TAT', 'TU',
                        'TU_L','TU_m', 'TU_G', 'TU_s', 'TU_sc',
                        'pMPH',
                        'PBMC',
                        'buffy coat',
                        'plasma',
                        }


def check_patient_id(patient_id):
    if patient_id.startswith("OVCA"):
        if not re.match("^OVCA\d+$", patient_id):
            raise ValueError("Patient id must be OVCA\\d if it starts with OVCA")


def check_dataframe(name, df):
    basename = os.path.basename(name)
    if not basename.startswith('_') and not name.startswith('_'):  # no fixed requirements on _meta dfs
        missing = set(must_have_columns).difference(df.columns)
        if missing:
            raise ValueError("%s is missing columns: %s, had %s" % (name, missing, df.columns))
    elif name.endswith('_exclusion'):
        mhc = ['patient', 'reason']
        missing = set(mhc).difference(df.columns)
        if missing:
            raise ValueError("%s is missing columns: %s, had %s" % (name, missing, df.columns))

    if 'compartment' in df.columns:
        x = set(df['compartment'].unique()).difference(allowed_compartments)
        if x:
            raise ValueError("invalid compartment(s) found in %s: %s" % (name, x,))
    if 'patient' in df.columns:
        [check_patient_id(x) for x in df['patient']]
    for x in 'variable', 'unit':
        if x in df.columns:
            if pd.isnull(df[x]).any():
                raise ValueError("%s must not be nan in %s" % (x, name))

def fix_the_darn_string(x):
    if isinstance(x, str):
        x = x.decode('utf-8')
    try:
        return unicode(x)
    except:
        print(repr(x))
        print(type(x))
        print(x)
        import pickle
        with open("debug.dat", 'w') as op:
            pickle.dump(x, op)
        raise

def categorical_where_appropriate(df):
    """make sure numerical columns are numeric
    and string columns that have less than 10% unique values are categorical
    and everything is unicode!

    """
    for c in df.columns:
        if df.dtypes[c] == object:
            try:
                df[c] = pd.to_numeric(df[c], errors='raise')
            except ValueError:
                if len(df[c].unique()) <= len(df) * 0.1:
                    df[c] = pd.Categorical(df[c])
                    new_cats = [fix_the_darn_string(x) for x in df[c].cat.categories]
                    df[c].cat.categories = new_cats
                else:
                    df[c] = [fix_the_darn_string(x) for x in df[c]]
    df.columns = [fix_the_darn_string(x) for x in df.columns]
    df.index.names = [ fix_the_darn_string(x) for x in df.index.names]


def create_biobank(
        dict_of_dataframes, name, revision, filename):
    """Create a file suitable for biobank consumption.
    Assumes all dataframes have at least variable, unit, patient, compartment and value columns
    """
    dict_of_dataframes['_meta/biobank'] = pd.DataFrame([
        {'variable': 'biobank', 'value': name, },
        {'variable': 'revision', 'value': revision, },
    ])
    patient_compartment_dataset = {
        'patient': [], 'compartment': [], 'dataset': []}
    for name, df in dict_of_dataframes.items():
        basename = os.path.basename(name)
        check_dataframe(name, df)
        categorical_where_appropriate(df)
        if not name.startswith('_') and not basename.startswith('_'):
            here = set(
                df[['patient', 'compartment']].itertuples(index=False, name=None))
            for p, c in here:
                patient_compartment_dataset['patient'].append(p)
                patient_compartment_dataset['compartment'].append(c)
                patient_compartment_dataset['dataset'].append(name)
            print ('converting types', name)
            # enforce alphabetical column order after default columns
            df = df[must_have_columns +
                    sorted([x for x in df.columns if x not in must_have_columns])]
    patient_compartment_dataset = pd.DataFrame(patient_compartment_dataset)
    categorical_where_appropriate(patient_compartment_dataset)
    dict_of_dataframes['_meta/patient_compartment_dataset'] = patient_compartment_dataset
    zfs = zipfile.ZipFile(filename, 'w')
    for name, df in dict_of_dataframes.items():
        zfs.writestr(name, df.to_msgpack())


def split_seperate_me(out_df):
    """Helper for creating biobank compatible dataframes.
    splits a column 'seperate_me' with OVCA12-compartment
    into seperate patient and compartment columns"""
    return out_df.assign(
        patient=[
            x[:x.find('-')] if '-' in x else x for x in out_df['seperate_me']],
        compartment=[x[x.find('-') + 1:] if '-' in x else 'n.a.' for x in out_df['seperate_me']]).drop('seperate_me', axis=1)


def write_dfs(dict_of_dfs):
    """Helper used by the notebooks to dump the dataframes for import"""
    for name, df_and_comment in dict_of_dfs.items():
        df, comment = df_and_comment
        check_dataframe(name, df)
        d = os.path.dirname(name)
        target_path = os.path.join('../../processed', d)
        if not os.path.exists(target_path):
            os.makedirs(target_path)
        fn = os.path.join(target_path, os.path.basename(name))
        df.to_pickle(fn)
        with open(fn, 'a') as op:
            pickle.dump(comment, op, pickle.HIGHEST_PROTOCOL)
