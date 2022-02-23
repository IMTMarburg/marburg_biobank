import zipfile
import os
import numpy as np
import pandas as pd
from pathlib import Path

__version__ = '0.156'

try:
    from functools import lru_cache
except (ImportError, AttributeError):
    # don't know how to tell setup.py that we only need functools32 when under 2.7.
    # so we'll just include a copy (*bergh*)
    import sys

    sys.path.append(os.path.join(os.path.dirname(__file__), "functools32"))
    from functools32 import lru_cache


class WideNotSupported(ValueError):
    def __init__(self):
        self.message = (
            ".get_wide() is not supported for this dataset. Use .get_dataset() instead"
        )


class CantApplyExclusion(ValueError):
    pass

datasets_to_cache = 32

known_compartment_columns = [
    "compartment",
    "cell_type",
    "disease",
    "culture_method",  # for those cells we can't take into sequencing ex vivo
    # these are only for backward compability
    "tissue",
    "disease-state",
]  # tissue


def lazy_member(field):
    """Evaluate a function once and store the result in the member (an object specific in-memory cache)
    Beware of using the same name in subclasses!
    """

    def decorate(func):
        if field == func.__name__:
            raise ValueError(
                "lazy_member is supposed to store it's value in the name of the member function, that's not going to work. Please choose another name (prepend an underscore..."
            )

        def doTheThing(*args, **kw):
            if not hasattr(args[0], field):
                setattr(args[0], field, func(*args, **kw))
            return getattr(args[0], field)

        return doTheThing

    return decorate


class Biobank(object):
    """An interface to a dump of our Biobank.
    Also used internally by the biobank website to access the data.

    In essence, a souped up dict of pandas dataframes stored
    as pickles in a zip file with memory caching"""

    def __init__(self, filename):
        self.filename = filename
        self.zf = zipfile.ZipFile(filename)
        if not "_meta/_data_format" in self.zf.namelist():
            self.data_format = "msg_pack"
        else:
            with self.zf.open("_meta/_data_format") as op:
                self.data_format = op.read().decode("utf-8")
        if self.data_format not in ("msg_pack", "parquet"):
            raise ValueError(
                "Unexpected data format (%s). Do you need to update marburg_biobank"
                % (self.data_format)
            )
        self._cached_datasets = {}

    @property
    def tall(self):
        return _BiobankItemAccessor(self.list_datasets, lambda dataset: self.get_dataset(dataset, apply_exclusion=True))

    @property
    def wide(self):
        return _BiobankItemAccessor(self.list_datasets, lambda dataset: self.get_wide(dataset, apply_exclusion=True))

    def get_all_patients(self):
        df = self.get_dataset("_meta/patient_compartment_dataset")
        return set(df["patient"].unique())

    def number_of_patients(self):
        """How many patients/indivuums are in all datasets?"""
        return len(self.get_all_patients())

    def number_of_datasets(self):
        """How many different datasets do we have"""
        return len(self.list_datasets())

    def get_compartments(self):
        """Get all compartments we have data for"""
        pcd = self.get_dataset("_meta/patient_compartment_dataset")
        return pcd

    @lru_cache(datasets_to_cache)
    def get_dataset_compartments(self, dataset):
        """Get available compartments in dataset @dataset"""
        ds = self.get_dataset(dataset)
        columns = self.get_dataset_compartment_columns(dataset)
        if not columns:
            return []
        else:
            sub_ds = ds[columns]
            sub_ds = sub_ds[~sub_ds.duplicated()]
            result = []
            for dummy_idx, row in sub_ds.iterrows():
                result.append(tuple([row[x] for x in columns]))
            return set(result)

    @lru_cache(datasets_to_cache)
    def get_dataset_compartment_columns(self, dataset):
        """Get available compartments columns in dataset @dataset"""
        ds = self.get_dataset(dataset)
        columns = [
            x for x in known_compartment_columns if x in ds.columns
        ]  # compartment included for older datasets
        return columns

    @lru_cache(datasets_to_cache)
    def get_variables_and_units(self, dataset):
        """What variables are availabe in a dataset?"""
        df = self.get_dataset(dataset)
        if len(df["unit"].cat.categories) == 1:
            vars = df["variable"].unique()
            unit = df["unit"].iloc[0]
            return set([(v, unit) for v in vars])
        else:
            x = df[["variable", "unit"]].drop_duplicates(["variable", "unit"])
            return set(zip(x["variable"], x["unit"]))

    def get_possible_values(self, dataset, variable, unit):
        df = self.get_dataset(dataset)
        return df["value"][(df["variable"] == variable) & (df["unit"] == unit)].unique()

    @lazy_member("_cache_list_datasets")
    def list_datasets(self):
        """What datasets to we have"""
        if self.data_format == "msg_pack":
            return sorted(
                [
                    name
                    for name in self.zf.namelist()
                    if not name.startswith("_")
                    and not os.path.basename(name).startswith("_")
                ]
            )
        elif self.data_format == "parquet":
            return sorted(
                [
                    name[: name.rfind("/")]
                    for name in self.zf.namelist()
                    if not name.startswith("_")
                    and not os.path.basename(name[: name.rfind("/")]).startswith("_")
                    and name.endswith("/0")
                ]
            )

    @lazy_member("_cache_list_datasets_incl_meta")
    def list_datasets_including_meta(self):
        """What datasets to we have"""
        if self.data_format == "msg_pack":
            return sorted(self.zf.namelist())
        elif self.data_format == "parquet":
            import re

            raw = self.zf.namelist()
            without_numbers = [
                x if not re.search("/[0-9]+$", x) else x[: x.rfind("/")] for x in raw
            ]
            return sorted(set(without_numbers))

    @lazy_member("_datasets_with_name_lookup")
    def datasets_with_name_lookup(self):
        return [ds for (ds, df) in self.iter_datasets() if "name" in df.columns]

    def name_lookup(self, dataset, variable):
        df = self.get_dataset(dataset)
        # todo: optimize using where?
        return df[df.variable == variable]["name"].iloc[0]

    def variable_or_name_to_variable_and_unit(self, dataset, variable_or_name):
        df = self.get_dataset(dataset)[["variable", "name", "unit"]]
        rows = df[(df.variable == variable_or_name) | (df.name == variable_or_name)]
        if len(rows["variable"].unique()) > 1:
            raise ValueError(
                "variable_or_name_to_variable led to multiple variables (%i): %s"
                % (len(rows["variable"].unique()), rows["variable"].unique())
            )
        try:
            r = rows.iloc[0]
        except IndexError:
            raise KeyError("Not found: %s" % variable_or_name)
        return r["variable"], r["unit"]

    def _get_dataset_columns_meta(self):
        import json

        with self.zf.open("_meta/_to_wide_columns") as op:
            return json.loads(op.read().decode("utf-8"))

    def has_wide(self, dataset):
        if dataset.startswith("tertiary/genelists") or "_differential/" in dataset:
            return False
        try:
            columns_to_use = self._get_dataset_columns_meta()
        except KeyError:
            return True
        if dataset in columns_to_use and not columns_to_use[dataset]:
            return False
        return True

    @lru_cache(maxsize=datasets_to_cache)
    def get_wide(
        self,
        dataset,
        apply_exclusion=True,
        standardized=False,
        filter_func=None,
        column="value",
    ):
        """Return dataset in row=variable, column=patient format.
        if @standardized is True Index is always (variable, unit) or (variable, unit, name), 
        and columns always (patient, [compartment, cell_type, disease])

        Otherwise, unit and compartment will be left off if there is only a 
        single value for them in the dataset
        if @apply_exclusion is True, excluded patients will be filtered from DataFrame

         @filter_func is run on the dataset before converting to wide, it
         takes a df, returns a modified df

        """
        dataset = self.dataset_exists(dataset)
        if not self.has_wide(dataset):
            raise WideNotSupported()
        df = self.get_dataset(dataset)
        if filter_func:
            df = filter_func(df)

        index = ["variable"]
        columns = self._get_wide_columns(dataset, df, standardized)
        if standardized or len(df.unit.cat.categories) > 1:
            index.append("unit")
        if "name" in df.columns:
            index.append("name")
        # if 'somascan' in dataset:
        # raise ValueError(dataset, df.columns, index ,columns)
        dfw = self.to_wide(df, index, columns, column=column)
        if apply_exclusion:
            try:
                return self.apply_exclusion(dataset, dfw)
            except CantApplyExclusion:
                return dfw
        else:
            return dfw

    def _get_wide_columns(self, dataset, tall_df, standardized):
        try:
            columns_to_use = self._get_dataset_columns_meta()
        except KeyError:
            columns_to_use = {}
        if dataset in columns_to_use:
            columns = columns_to_use[dataset]
            if standardized:
                for x in known_compartment_columns:
                    if not x in columns:
                        columns.append(x)
                    if x in tall_df.columns and (
                        (
                            hasattr(tall_df[x], "cat")
                            and (len(tall_df[x].cat.categories) > 1)
                        )
                        or (len(tall_df[x].unique()) > 1)
                    ):
                        pass
                    else:
                        if standardized and x not in tall_df.columns:
                            tall_df = tall_df.assign(**{x: np.nan})
        else:
            if "vid" in tall_df.columns and not "patient" in tall_df.columns:
                columns = ["vid"]
            elif "patient" in tall_df.columns:
                columns = ["patient"]
            else:
                raise ValueError(
                    "Do not know how to convert this dataset to wide format."
                    " Retrieve it get_dataset() and call to_wide() manually with appropriate parameters."
                )
            for x in known_compartment_columns:
                if x in tall_df.columns or (standardized and x != "compartment"):
                    if not x in columns:
                        columns.append(x)
                    if x in tall_df.columns and (
                        (
                            hasattr(tall_df[x], "cat")
                            and (len(tall_df[x].cat.categories) > 1)
                        )
                        or (len(tall_df[x].unique()) > 1)
                    ):
                        pass
                    else:
                        if standardized and x not in tall_df.columns:
                            tall_df = tall_df.assign(**{x: np.nan})
                        elif not standardized:
                            if (
                                hasattr(tall_df[x], "cat")
                                and (len(tall_df[x].cat.categories) == 1)
                            ) or (len(tall_df[x].unique()) == 1):
                                if x in columns:
                                    columns.remove(x)
        return columns

    def to_wide(
        self,
        df,
        index=["variable"],
        columns=known_compartment_columns,
        sort_on_first_level=False,
        column='value',
    ):
        """Convert a dataset (or filtered dataset) to a wide DataFrame.
        Preferred to pd.pivot_table manually because it is
           a) faster and
           b) avoids a bunch of pitfalls when working with categorical data and
           c) makes sure the columns are dtype=float if they contain nothing but floats

        index = variable,unit
        columns = (patient, compartment, cell_type)
        """
        if columns == known_compartment_columns:
            columns = [x for x in columns if x in df.columns]
        # raise ValueError(df.columns,index,columns)
        chosen = [column] + index + columns
        df = df.loc[:, [x for x in chosen if x in df.columns]]
        for x in chosen:
            if x not in df.columns:
                df = df.assign(**{x: np.nan})
        set_index_on = index + columns
        columns_pos = tuple(range(len(index), len(index) + len(columns)))
        res = df.set_index(set_index_on).unstack(columns_pos)
        c = res.columns
        c = c.droplevel(0)
        # this removes categories from the levels of the index. Absolutly
        # necessary, or you can't add columns later otherwise
        if isinstance(c, pd.MultiIndex):
            try:
                c = pd.MultiIndex(
                    [list(x) for x in c.levels], codes=c.codes, names=c.names
                )
            except AttributeError:
                c = pd.MultiIndex(
                    [list(x) for x in c.levels], labels=c.labels, names=c.names
                )
        else:
            c = list(c)
        res.columns = c
        single_unit = not 'unit' in df.columns or len(df['unit'].unique()) == 1
        if isinstance(c, list):
            res.columns.names = columns
        if sort_on_first_level:
            # sort on first level - ie. patient, not compartment - slow though
            res = res[sorted(list(res.columns))]
        for c in res.columns:
            x = res[c].fillna(value=np.nan, inplace=False)
            if (x == None).any():  # noqa: E711
                raise ValueError("here")
            if single_unit: # don't do this for multiple units -> might have multiple dtypes
                try:
                    res[c] = pd.to_numeric(x, errors="raise")
                except (ValueError, TypeError):  # leaving the Nones as Nones
                    pass
        return res

    @lru_cache(maxsize=datasets_to_cache)
    def get_excluded_patients(self, dataset):
        """Which patients are excluded from this particular dataset (or globally)?.

        May return a set of patient_id, or tuples of (('patient', 'x'y'), ('compartment1', 'xyz'),...) tuples if only
        certain compartments where excluded.

        """
        try:
            global_exclusion_df = self.get_dataset("clinical/_other_exclusion")
            excluded = set(global_exclusion_df["patient"].unique())
        except KeyError:
            excluded = set()
        # local exclusion from this dataset
        try:
            exclusion_path = (
                os.path.dirname(dataset)
                + "/"
                + "_"
                + os.path.basename(dataset)
                + "_exclusion"
            )
            exclusion_df = self.get_dataset(exclusion_path)
        except KeyError:
            return excluded
        columns = ["patient"] + self.get_dataset_compartment_columns(dataset)
        columns = [x for x in columns if x in exclusion_df.columns]
        res = exclusion_df[columns]
        if set(res.columns) == set(["patient"]):
            excluded.update(exclusion_df["patient"].unique())
        else:
            for idx, row in exclusion_df.iterrows():
                d = []
                for c in columns:
                    d.append((c, row[c]))
                excluded.add(tuple(d))
        return excluded

    def apply_exclusion(self, dataset_name, df):
        dataset_name = self.dataset_exists(dataset_name)
        excluded = self.get_excluded_patients(dataset_name)
        # columns = ["patient"] + self.get_dataset_compartment_columns(dataset_name)
        if "patient" in df.columns:  # a tall dataset
            keep = np.ones((len(df),), np.bool)
            for x in excluded:
                if isinstance(x, tuple):
                    matching = np.ones((len(df),), np.bool)
                    for column, value in x:
                        matching &= df[column] == value
                    keep = keep & ~matching
                else:
                    keep = keep & ~(df["patient"] == x)
            return df[keep]
        elif df.index.names[0] == "variable":  # a wide dataset...
            to_remove = []
            for c in df.columns:
                if isinstance(c, tuple):
                    if c[0] in excluded:  # patient totaly excluded
                        to_remove.append(c)
                    else:
                        key = tuple(zip(df.columns.names, c))
                        if key in excluded:
                            to_remove.append(c)
                else:
                    if c in excluded:
                        to_remove.append(c)
            return df.drop(to_remove, axis=1)
        else:
            raise CantApplyExclusion(
                "Sorry, not a tall or wide DataFrame that I know how to handle."
            )

    @lru_cache(maxsize=1)
    def get_exclusion_reasons(self):
        """Get exclusion information for all the datasets + globally"""
        result = {}
        global_exclusion_df = self.get_dataset("clinical/_other_exclusion")
        for tup in global_exclusion_df.itertuples():
            if tup.patient not in result:
                result[tup.patient] = {}
            result[tup.patient]["global"] = tup.reason
        for dataset in self.list_datasets():
            try:
                exclusion_df = self.get_dataset(
                    os.path.dirname(dataset)
                    + "/"
                    + "_"
                    + os.path.basename(dataset)
                    + "_exclusion"
                )
                for tup in exclusion_df.itertuples():
                    if tup.patient not in result:
                        result[tup.patient] = {}
                    result[tup.patient][dataset] = tup.reason
            except KeyError:
                pass
        return result

    def iter_datasets(self, yield_meta=False):
        if yield_meta:
            lst = self.list_datasets_including_meta()
        else:
            lst = self.list_datasets()
        for name in lst:
            yield name, self.get_dataset(name)

    def dataset_exists(self, name):
        datasets = self.list_datasets_including_meta()
        out = False
        if name in datasets:
            out = name
        else:
            next = "primary/" + name
            if next in datasets:
                out = next
            else:
                if name.startswith('secondary/'):
                    next = 'tertiary' + name[name.find('/'):]
                    if next in datasets:
                        out = next
        if not out:
            msg = "No such dataset: %s." % name
            import difflib

            msg += "Suggestions: "
            for x in difflib.get_close_matches(name, datasets):
                msg += " " + x + " "
            msg += ". Use .list_datasets() to view all datasets"
            raise KeyError(msg)
        return out

    def __load_df_from_parquet(self, name):
        try:
            import pyarrow
        except ImportError:
            try:
                import fastparquet
            except ImportError:
                raise ValueError("marburg_biobank needs either pyarrow or fastparquet")

        try:
            with self.zf.open(name) as op:
                return pd.read_parquet(op)
        except Exception as e:
            if (
                "UnsupportedOperation" in str(e)
                or "has no attribute" in str(e)
                or "UnsupportedOperation" in repr(e)
            ):  # python prior 3.7 has no seek on zipfiles
                import io

                with self.zf.open(name) as op:
                    b = io.BytesIO()
                    b.write(op.read())
                    b.seek(0)
                    return pd.read_parquet(b)
            elif 'not a path-like object' in str(e):
                import tempfile
                with tempfile.NamedTemporaryFile(suffix=".biobank.parquet") as tf:
                    with self.zf.open(name) as op:
                        tf.write(op.read())
                    tf.flush()
                    return pd.read_parquet(tf.name)
            else:
                raise
        raise NotImplementedError()

    @lru_cache(datasets_to_cache)
    def get_dataset(self, name, apply_exclusion=False):
        """Retrieve a dataset"""
        name = self.dataset_exists(name)
        if self.data_format == "msg_pack":
            try:
                import mbf_pandas_msgpack
            except (ImportError, AttributeError):
                raise ImportError("Please install mbf-pandas-msgpack to read this old school biobank file")
            with self.zf.open(name) as op:
                try:
                    df = mbf_pandas_msgpack.read_msgpack(op.read())
                except KeyError as e:
                    if "KeyError: u'category'" in str(e):
                        raise ValueError(
                            "Your pandas is too old. You need at least version 0.18"
                        )
        elif self.data_format == "parquet":
            try:
                import pyarrow
            except ImportError:
                try:
                    import fastparquet
                except ImportError:
                    raise ValueError("marburg_biobank needs either pyarrow or fastparquet")

            ds = self.zf.namelist()
            ii = 0
            dfs = []
            sub_name = name + "/" + str(ii)
            while sub_name in ds:
                dfs.append(self.__load_df_from_parquet(sub_name))
                ii += 1
                sub_name = name + "/" + str(ii)
            if not dfs:  # not actually a unit splitted dataframe - meta?
                df = self.__load_df_from_parquet(name)
            elif len(dfs) == 1:
                df = dfs[0]
            else:
                categoricals = set()
                for df in dfs:
                    for c, dt in df.dtypes.items():
                        if dt.name == "category":
                            categoricals.add(c)
                df = pd.concat(dfs)
                reps = {c: pd.Categorical(df[c]) for c in categoricals}
                if reps:
                    df = df.assign(**reps)
        else:
            raise ValueError(
                "Unexpected data format. Do you need to upgrade marburg_biobank?"
            )
        if apply_exclusion:
            try:
                df = self.apply_exclusion(name, df)
            except CantApplyExclusion:
                return df
        return df

    def get_comment(self, name):
        comments = self.get_dataset("_meta/comments")
        if len(comments) == 0:
            return ""
        match = comments.path == name
        if match.any():
            return comments[match].iloc[0]["comment"]
        else:
            return ""

    def get_changelog(self):
        try:
            return self.get_dataset("_meta/_changelog").sort_values("revision")
        except KeyError:
            raise ValueError(
                "This revision of the biobank did not include a change log."
            )

def biobank_to_url(biobank):
    if biobank.lower() == 'ovca':
        return "https://mbf.imt.uni-marburg.de/biobank"
    elif biobank.lower() == 'paad':
        return "https://mbf.imt.uni-marburg.de/biobank_paad"
    else:
        raise ValueError(f"Don't know how to download {biobank}")


def _find_newest_revision(username, password, revision, biobank):
    import requests
    url = biobank_to_url(biobank) + '/download/find_newest_revision'
    if revision: # find teh newest sub release (eg. find 20.3 from 20)
        url += "?revision=%s" % revision
    r = requests.get(
        url, stream=True, auth=requests.auth.HTTPBasicAuth(username, password)
    )
    if r.status_code != 200:
        raise ValueError("Non 200 OK Return - was %s" % r.status_code)
    return r.text


passwd_file = Path("~/.ovca_biobank_password").expanduser()
def query_user():
    import sys
    from getpass import getpass
    if passwd_file.exists():
        print("Reading password from " + str(passwd_file))
        username, password = passwd_file.read_text().split("\n")[:2]
        store = False
    else:
        username = input("please enter your username")
        password = getpass("please enter your password")
        store = True
    return username, password, store



def download_and_open(username=False, password=False, revision=None, biobank='ovca'):
    from pathlib import Path
    import requests
    import shutil

    store = False
    if username is False and password is False:
        print("query")
        username,password, store = query_user()

    newest = _find_newest_revision(username, password, revision, biobank)
    if store:
        print("Storing password to " + str(passwd_file))
        passwd_file.write_text("%s\n%s\n" % (username, password))
    if revision is None:
        print("newest revision is", newest)
    else:
        print("newest revision for %s is %s" % (revision, newest))
    fn = "marburg_%s_biobank_%s.zip" % (biobank, newest)
    if not Path(fn).exists():
        print("downloading biobank revision %s" % newest)
        url = biobank_to_url(biobank) + "/download/marburg_biobank?revision=%s" % newest
        r = requests.get(
            url, stream=True, auth=requests.auth.HTTPBasicAuth(username, password)
        )
        if r.status_code != 200:
            raise ValueError("Non 200 OK Return - was %s" % r.status_code)
        r.raw.decode_content = True
        fh = open(fn, "wb")
        shutil.copyfileobj(r.raw, fh)
        fh.close()
    else:
        print("using local copy %s" % fn)
    return Biobank(fn)


class _BiobankItemAccessor:
    def __init__(self, list_callback, get_callback):
        self.list_callback = list_callback
        self.get_callback = get_callback

    def __getitem__(self, key):
        return self.get_callback(key)

    def _ipython_key_completions_(self):
        return self.list_callback()

OvcaBiobank = Biobank # old school code support
