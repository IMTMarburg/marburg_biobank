import pandas as pd
import time
import re
import pickle
import zipfile
import os
import collections



# for the primary data
must_have_columns = ['variable', 'unit', 'value', 'patient', ]
# for 'secondary' datasets
must_have_columns_secondary = ['variable', 'unit', 'value']
allowed_cells = {
    'n.a',
    'T',
    'macrophage',
    'tumor',
    'tumor_s',
    'tumor_sc',
    'tumor_m',
    'tumor_L',
    'tumor_G',
    'MDSC',
}
allowed_tissues = {'blood',
                   'ascites',
                   'n.a.'
                   }
allowed_disease_states = {
    'cancer',
    'healthy',
    'benign',
    'n.a.',


}


def check_patient_id(patient_id):
    if patient_id.startswith("OVCA"):
        if not re.match("^OVCA\d+$", patient_id):
            raise ValueError(
                "Patient id must be OVCA\\d if it starts with OVCA")
        return 'cancer'
    else:
        return 'non-cancer'


def check_dataframe(name, df):
    if 'variable' in df.columns:
        df = df.assign(
            variable=[x.encode('utf-8') if isinstance(x, str) else x for x in df.variable])
    for c in 'compartment', 'seperate_me':
        if c in df.columns:
            raise ValueError(
                "%s must no longer be a df column - %s " % (c, name))
    basename = os.path.basename(name)
    # no fixed requirements on _meta dfs
    if not basename.startswith('_') and not name.startswith('_'):
        if name.startswith('secondary'):
            mh = set(must_have_columns_secondary)
        else:
            mh = set(must_have_columns)
        missing = mh.difference(df.columns)
        if missing:
            raise ValueError("%s is missing columns: %s, had %s" %
                             (name, missing, df.columns))
    elif name.endswith('_exclusion'):
        mhc = ['patient', 'reason']
        missing = set(mhc).difference(df.columns)
        if missing:
            raise ValueError("%s is missing columns: %s, had %s" %
                             (name, missing, df.columns))

    for column, allowed_values in [
        ('cell', allowed_cells),
        ('tissue', allowed_tissues),
        ('disease_state', allowed_disease_states),
    ]:
        if column in df.columns and not name.startswith('secondary/'):
            x = set(df[column].unique()).difference(allowed_values)
            if x:
                raise ValueError("invalid %s found in %s: %s" % (column, name, x,))

    if 'patient' in df.columns and not name.endswith('_exclusion'):
        states = set([check_patient_id(x) for x in df['patient']])
        if len(states) > 1:
            if 'disease_state' not in df.columns:
                raise ValueError(
                    "Datasets mixing cancer and non cancer data need a disease_state column:%s"  % (name,))

    for x in 'variable', 'unit':
        if x in df.columns:
            if pd.isnull(df[x]).any():
                raise ValueError("%s must not be nan in %s" % (x, name))
    if not basename.startswith('_') and not name.startswith('_'):
        for vu, group in df.groupby(['variable', 'unit']):
            variable, unit = vu
            if unit == 'string':
                pass
            elif unit == 'timestamp':
                for v in group.value:
                    if not isinstance(v, pd.Timestamp):
                        raise ValueError(
                            "Not timestamp data in %s %s" % vu)
            elif unit == 'bool':
                if set(group.value.unique()) != set([True, False]):
                    raise ValueError(
                        "Unexpected values for bool variables in %s %s" % vu)
            else:
                if not ((group.value.dtype == int) & (group.value.dtype == float)): #might not be floaty enough
                    for v in group.value:
                        if not isinstance(v, float) and not isinstance(v, int):
                            raise ValueError("Non float in %s, %s" % vu)



def fix_the_darn_string(x):
    if isinstance(x, bool):
        return x
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
    to_assign = {}
    for c in df.columns:
        if df.dtypes[c] == object:
            try:
                to_assign[c] = pd.to_numeric(df[c], errors='raise')
            except (ValueError, TypeError):
                if len(df[c].unique()) <= len(df) * 0.1:
                    to_assign[c] = pd.Categorical(df[c])
                    new_cats = [fix_the_darn_string(x) for x in to_assign[c].categories]
                    to_assign[c].categories = new_cats
                else:
                    to_assign[c] = [fix_the_darn_string(x) for x in df[c]]
    df = df.assign(**to_assign)
    df.columns = [fix_the_darn_string(x) for x in df.columns]
    df.index.names = [fix_the_darn_string(x) for x in df.index.names]
    return df

def extract_patient_compartment_meta(dict_of_dfs):
    output = []
    from . import known_compartment_columns
    columns = ['patient'] + known_compartment_columns
    for name in dict_of_dfs:
        if (not name.startswith('secondary/') and 
            not name.startswith('_') and not 
            os.path.basename(name).startswith('_')):
            df = dict_of_dfs[name]
            subset = df[[x for x in columns if x in df.columns]]
            subset = subset[~subset.duplicated()]
            for idx, row in subset.iterrows():
                row[u'dataset'] = unicode(name)
                output.append(row)
    return pd.DataFrame(output)


def create_biobank(
        dict_of_dataframes, name, revision, filename):
    """Create a file suitable for biobank consumption.
    Assumes all dataframes have at least variable, unit, patient, compartment and value columns
    """
    dict_of_dataframes['_meta/biobank'] = pd.DataFrame([
        {'variable': 'biobank', 'value': name, },
        {'variable': 'revision', 'value': revision, },
    ])
    for name, df in dict_of_dataframes.items():
        print("handling", name)
        basename = os.path.basename(name)
        s = time.time()
        check_dataframe(name, df)
        print('check time', time.time() - s)
        s = time.time()
        df = categorical_where_appropriate(df)
        print('cat time', time.time() - s)
        s = time.time()
        # enforce alphabetical column order after default columns
        df = df[[x for x in must_have_columns if x in df.columns] +
                sorted([x for x in df.columns if x not in must_have_columns])]
        print('column order time', time.time() - s)
        dict_of_dataframes[name] = df
    s = time.time()
    dict_of_dataframes["_meta/patient_compartment_dataset"] = extract_patient_compartment_meta(dict_of_dataframes)
    print("patient_compartment_dataset_time", time.time() - s)
    print("now writing zip file")
    zfs = zipfile.ZipFile(filename, 'w')
    for name, df in dict_of_dataframes.items():
        zfs.writestr(name, df.to_msgpack())
    zfs.close()
    #one last check it's all numbers...
    print("checking float")
    from . import OvcaBiobank
    bb = OvcaBiobank(filename)
    for ds in bb.list_datasets():
        if ds.startswith('secondary/'):
            continue
        print ds
        df = bb.get_wide(ds, filter_func=lambda df: df[~df.unit.isin(['timestamp','string', 'bool'])])
        #df = bb.get_wide(ds)    
        for idx, row in df.iterrows():
	    if row.dtype != float:
	        print("Error in %s %s, dtype was %s" % (ds, idx, row.dtype))    



def split_seperate_me(out_df, in_order=['patient', 'tissue']):
    """Helper for creating biobank compatible dataframes.
    splits a column 'seperate_me' with OVCA12-compartment
    into seperate patient and compartment columns"""
    split = [x.split("-") for x in out_df['seperate_me']]
    return out_df.assign(**{
        x: [y[ii] for y in split]
        for (ii, x) in enumerate(in_order)
    }).drop('seperate_me', axis=1)


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
