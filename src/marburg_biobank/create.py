import pandas as pd
import tempfile
import inspect
import pypipegraph as ppg
from pathlib import Path
import time
import re
import pickle
import zipfile
import os
import json
import base64
from . import WideNotSupported


settings = None


def apply_ovca_settings():
    global settings
    if settings is not None and settings["what"] != "OVCA":
        raise ValueError("different apply_*_settings being called")

    def check_patient_id(patient_id):
        if patient_id.startswith("OVCA"):
            if not re.match(r"^OVCA\d+(R[0-9]*)?$", patient_id):
                raise ValueError("Patient id must follow OVCA\\d(R[0-9]*)? if it starts with OVCA")
            return "cancer"
        elif patient_id.startswith("OC"):
            raise ValueError("OVCA patients must not start with OC")
        else:
            return "non-cancer"

    settings = {
        "what": "OVCA",
        # for the primary data
        "must_have_columns": ["variable", "unit", "value", "patient"],
        # for 'secondary' datasets
        "must_have_columns_secondary": ["variable", "unit", "value"],
        # for gene lists
        "must_have_columns_tertiary_genelists": ["stable_id", "gene"],
        "allowed_cells": {
            "T",
            "macrophage",
            "tumor",
            "tumor_s",
            "tumor_sc",
            "tumor_m",
            "tumor_L",
            "tumor_G",
            "MDSC",
            "NK",
            "n.a.",
            "adipocyte",
            "HPMC",
            "CAF",
        },
        "allowed_compartments": {"blood", "ascites", "n.a.", "omentum"},
        "allowed_disease_states": {"cancer", "healthy", "benign", "n.a."},
        "check_patient_id": check_patient_id,
        'database_filename_template': 'marburg_ovca_revision_%s.zip'
    }


def apply_paad_settings():
    "for the pancreas biobank"
    global settings
    if settings is not None and settings["what"] != "PAAD":
        raise ValueError("different apply_*_settings being called")

    def check_patient_id(patient_id):
        if patient_id.startswith("ACH"):
            if not re.match(r"^ACH-\d+$", patient_id):
                raise ValueError("Patient id must be ACH\\d if it starts with ACH")
            return "PAAD"
        else:
            raise ValueError(
                "PAAD patients must start with ACH (non-cancer samples yet to be suported in apply_paad_settings"
            )

    settings = {
        "what": "PAAD",
        # for the primary data
        "must_have_columns": ["variable", "unit", "value", "patient"],
        # for 'secondary' datasets
        "must_have_columns_secondary": ["variable", "unit", "value"],
        # for gene lists
        "must_have_columns_tertiary_genelists": ["stable_id", "gene"],
        "allowed_cells": {"solid_tumor_mix",},
        "allowed_compartments": {"tumor"},  # -
        "allowed_disease_states": {"PAAD",},
        "check_patient_id": check_patient_id,
        'database_filename_template': 'marburg_paad_biobank_revision_%s.zip'
    }


def check_dataframe(name, df):
    # why was this done?
    # if "variable" in df.columns:
    # df = df.assign(
    # variable=[
    # x.encode("utf-8") if isinstance(x, str) else x for x in df.variable
    # ]
    # )
    if settings is None:
        raise ValueError("Must call apply_*_settings (eg. apply_ovca_settings) first")
    for c in "seperate_me":
        if c in df.columns:
            raise ValueError("%s must no longer be a df column - %s " % (c, name))
    if "compartment" in df.columns and not "disease" in df.columns:
        raise ValueError("Columns must now be cell_type/disease/compartment split")
    if "patient" in df.columns:
        for patient in df["patient"]:
            settings["check_patient_id"](patient)
        #
    # dataframes ofter now are _actual_name/0-9+,
    # but possibly only after writing it out...
    if re.search("/[0-9]+$", name):
        name = name[: name.rfind("/")]
    basename = os.path.basename(name)
    # no fixed requirements on _meta dfs
    if not basename.startswith("_") and not name.startswith("_"):
        if (
            "_differential/" in name
            or "/genomics/"  # no special requirements for differential datasets for now
            in name  # mutation data is weird enough.
        ):
            mh = set()
        elif name.startswith("secondary") or name.startswith('tertiary/transcriptomics'):
            mh = set(settings["must_have_columns_secondary"])
        elif name.startswith("tertiary/genelists"):
            mh = set(settings["must_have_columns_tertiary_genelists"])
        elif name.startswith("tertiary/survival"):
            mh = set()
        else:
            mh = set(settings["must_have_columns"])
            for c in "cell", "disease_state", "tissue":
                if c in df.columns:
                    raise ValueError(
                        "%s must no longer be a df column - %s " % (c, name)
                    )

        missing = mh.difference(df.columns)
        if missing:
            raise ValueError(
                "%s is missing columns: %s, had %s" % (name, missing, df.columns)
            )
    elif name.endswith("_exclusion"):
        mhc = ["patient", "reason"]
        missing = set(mhc).difference(df.columns)
        if missing:
            raise ValueError(
                "%s is missing columns: %s, had %s" % (name, missing, df.columns)
            )

    for column, allowed_values in [
        ("cell_type", settings["allowed_cells"]),
        ("compartment", settings["allowed_compartments"]),
        ("disease", settings["allowed_disease_states"]),
    ]:
        if column in df.columns and not name.startswith("secondary/") and not name.startswith('tertiary/'):
            x = set(df[column].unique()).difference(allowed_values)
            if x:
                raise ValueError(
                    "invalid %s found in %s: %s - check marburg_biobank/create.py, allowed_* if you want to extend it"
                    % (column, name, x)
                )

    if "patient" in df.columns and not name.endswith("_exclusion"):
        states = set([settings["check_patient_id"](x) for x in df["patient"]])
        if len(states) > 1:
            if "disease" not in df.columns:
                raise ValueError(
                    "Datasets mixing cancer and non cancer data need a disease column:%s"
                    % (name,)
                )

    for x in "variable", "unit":
        if x in df.columns:
            try:
                if pd.isnull(df[x]).any():
                    raise ValueError("%s must not be nan in %s" % (x, name))
                if df[x].str.startswith(" ").any():
                    raise ValueError("At least one %s started with a space" % x)
                if df[x].str.endswith(" ").any():
                    raise ValueError("At least one %s ended with a space" % x)
            except:
                print("column", x)
                raise

    if (
        not basename.startswith("_")
        and not name.startswith("_")
        and not name.startswith("tertiary")
        and mh  # was not '_differential/' in name
    ):
        for vu, group in df.groupby(["variable", "unit"]):
            variable, unit = vu
            if unit == "string":
                pass
            elif unit == "timestamp":
                for v in group.value:
                    if not isinstance(v, pd.Timestamp):
                        raise ValueError("Not timestamp data in %s %s" % vu)
            elif unit == "bool":
                if set(group.value.unique()) != set([True, False]):
                    raise ValueError(
                        "Unexpected values for bool variables in %s %s" % vu
                    )
            else:
                if not (
                    (group.value.dtype == int) & (group.value.dtype == float)
                ):  # might not be floaty enough
                    for v in group.value:
                        if not isinstance(v, float) and not isinstance(v, int):
                            raise ValueError("Non float in %s, %s" % vu)


def fix_the_darn_string(x):
    if isinstance(x, bool):
        return x
    if isinstance(x, bytes):
        x = x.decode("utf-8")
    try:
        return str(x)
    except:  # noqa:E722
        print(repr(x))
        print(type(x))
        print(x)
        import pickle

        with open("debug.dat", "w") as op:
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
                to_assign[c] = pd.to_numeric(df[c], errors="raise")
            except (ValueError, TypeError):
                if len(df[c].unique()) <= len(df) * 0.3 or c == "patient":
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

    columns = ["patient"] + known_compartment_columns
    for name in dict_of_dfs:
        if (
            not name.startswith("secondary/")
            and not name.startswith("tertiary/")
            and not name.startswith("_")
            and not os.path.basename(name).startswith("_")
        ):
            df = dict_of_dfs[name]
            subset = df[[x for x in columns if x in df.columns]]
            subset = subset[~subset.duplicated()]
            for idx, row in subset.iterrows():
                row[u"dataset"] = str(name)
                output.append(row)
    return pd.DataFrame(output)


def create_biobank(dict_of_dataframes, name, revision, filename, to_wide_columns):
    """Create a file suitable for biobank consumption.
    Assumes all dataframes pass check_dataframe
    """
    if settings is None:
        raise ValueError("Must call apply_*_settings (eg. apply_ovca_settings) first")
    dict_of_dataframes["_meta/biobank"] = pd.DataFrame(
        [
            {"variable": "biobank", "value": name},
            {"variable": "revision", "value": revision},
        ]
    )
    for name, df in dict_of_dataframes.items():
        print("handling", name)
        # basename = os.path.basename(name)
        s = time.time()
        check_dataframe(name, df)
        print("check time", time.time() - s)
        s = time.time()
        df = categorical_where_appropriate(df)
        print("cat time", time.time() - s)
        s = time.time()
        # enforce alphabetical column order after default columns
        df = df[
            [x for x in settings["must_have_columns"] if x in df.columns]
            + sorted([x for x in df.columns if x not in settings["must_have_columns"]])
        ]
        print("column order time", time.time() - s)
        dict_of_dataframes[name] = df
    s = time.time()
    dict_of_dataframes[
        "_meta/patient_compartment_dataset"
    ] = extract_patient_compartment_meta(dict_of_dataframes)
    print("patient_compartment_dataset_time", time.time() - s)
    print("now writing zip file")
    zfs = zipfile.ZipFile(filename, "w")
    for name, df in dict_of_dataframes.items():
        tf = tempfile.NamedTemporaryFile(mode="w+b", suffix=".pq")
        df.to_parquet(tf)
        tf.flush()
        tf.seek(0, 0)
        zfs.writestr(name, tf.read())
    zfs.writestr("_meta/_to_wide_columns", json.dumps(to_wide_columns))
    zfs.writestr("_meta/_data_format", "parquet")
    zfs.close()
    # one last check it's all numbers...
    print("checking float")
    from . import OvcaBiobank

    # check that we can do the get_wide on all of them
    bb = OvcaBiobank(filename)
    for ds in bb.list_datasets():
        try:
            df = bb.get_wide(
                ds,
                filter_func=lambda df: df[
                    ~df.unit.isin(["timestamp", "string", "bool"])
                ],
            )
        except WideNotSupported:
            continue
        except:
            print("issue is in", ds)
            raise
        # df = bb.get_wide(ds)
        for idx, row in df.iterrows():
            if row.dtype != float:
                print("Error in %s %s, dtype was %s" % (ds, idx, row.dtype))


def split_seperate_me(out_df, in_order=["patient", "compartment"]):
    """Helper for creating biobank compatible dataframes.
    splits a column 'seperate_me' with OVCA12-compartment
    into seperate patient and compartment columns"""
    split = [x.split("-") for x in out_df["seperate_me"]]
    return out_df.assign(
        **{x: [y[ii] for y in split] for (ii, x) in enumerate(in_order)}
    ).drop("seperate_me", axis=1)


def write_dfs(dict_of_dfs):
    """Helper used by the notebooks to dump the dataframes for import"""
    for name, df_and_comment in dict_of_dfs.items():
        if isinstance(df_and_comment, tuple):
            df, _comment = df_and_comment
        else:
            df = df_and_comment
        check_dataframe(name, df)
        d = os.path.dirname(name)
        target_path = os.path.join("/project/processed", d)
        if not os.path.exists(target_path):
            os.makedirs(target_path)
        fn = os.path.join(target_path, os.path.basename(name))
        df.to_pickle(fn)
        # with open(fn, "a") as op:
        # pickle.dump(comment, op, pickle.HIGHEST_PROTOCOL)


exporting_classes = []


def exporting_class(cls):
    exporting_classes.append(cls)
    return cls


def prep_desc(x):
    x = x.strip()
    import re

    x = re.sub("\n[ ]+", "\n", x)
    return x


def exporting_method(output_name, description, input_files, deps, wide_columns):
    def inner(func):
        frame = inspect.stack()[1]
        filename = frame.filename
        cwd = os.getcwd()
        os.chdir(Path(filename).parent)
        func._output_name = output_name
        func._description = prep_desc(description)
        func._input_files = [Path(x).absolute() for x in input_files]
        func._deps = deps
        func._wide_columns = wide_columns
        if not isinstance(wide_columns, list):
            raise ValueError("wide_columns must be a list of str")
        os.chdir(cwd)
        func._abs_filename = str(Path(func.__code__.co_filename).absolute())
        return func

    return inner


def run_exports(gen_additional_jobs=None, handle_ppg=True, settings='ovca'):
    if settings == 'ovca':
        apply_ovca_settings()
    else:
        raise ValueError("unknow setting value", settings)

    old = Path(os.getcwd()).absolute()
    os.chdir("/project")
    if handle_ppg:
        ppg.new_pipegraph()
    # os.chdir(old)
    to_wide_columns = {}
    jobs = []
    for cls in exporting_classes:
        instance = cls()
        if hasattr(instance, "exports"):
            instance.exports()
        
        out_prefix = getattr(instance, "out_prefix", "")
        for method_name in dir(instance):
            method = getattr(instance, method_name)
            if hasattr(method, "_output_name"):
                print(cls.__name__, method.__name__)
                output_filename = (
                    "/project/processed/" + out_prefix + method._output_name + ".units"
                )
                cwd = str(Path(method._abs_filename).parent)

                def write(output_filename=output_filename, method=method, cwd=cwd):
                    os.chdir(cwd)
                    df = method()
                    os.chdir("/project")
                    check_dataframe(out_prefix + method._output_name, df)
                    Path(output_filename).parent.mkdir(exist_ok=True, parents=True)
                    if "unit" in df:
                        for ii, (unit, sub_df) in enumerate(
                            df.groupby("unit", sort=True)
                        ):
                            try:
                                sub_df.to_parquet(
                                    output_filename[: output_filename.rfind(".")]
                                    + "."
                                    + str(ii)
                                    + ".parquet"
                                )
                            except:
                                sub_df.to_pickle("debug.pickle")
                                raise

                        Path(output_filename).write_text(
                            json.dumps(sorted(df.unit.unique()))
                        )
                    else:
                        df.to_parquet(
                            output_filename[: output_filename.rfind(".")] + ".0.parquet"
                        )
                        Path(output_filename).write_text(json.dumps(["nounit"]))
                    Path(output_filename + ".desc").write_text(method._description)

                job = ppg.MultiFileGeneratingJob(
                    [output_filename, output_filename + ".desc"], write
                )
                job.depends_on(
                    ppg.FunctionInvariant(output_filename + "_inner_func", method)
                )
                if method._input_files:
                    job.depends_on(ppg.MultiFileInvariant(method._input_files))
                if method._deps:
                    if hasattr(method._deps, "__call__"):
                        deps = method._deps(method.__self__)
                    else:
                        deps = method._deps
                    job.depends_on(deps)

                print(output_filename)
                print("")
                os.chdir("/project")
                jobs.append(job)
                to_wide_columns[out_prefix + method._output_name] = method._wide_columns

    def dump_to_wide_columns(output_filename):
        Path(output_filename).write_text(json.dumps(to_wide_columns))

    jobs.append(
        ppg.FileGeneratingJob(
            "/project/processed/_to_wide_columns.json", dump_to_wide_columns
        ).depends_on(
            ppg.ParameterInvariant(
                "/project/processed/_to_wide_columns.json",
                ppg.util.freeze(to_wide_columns),
            )
        )
    )

    old = Path(os.getcwd()).absolute()
    if handle_ppg:
        os.chdir("/project")
        ppg.run_pipegraph()
    os.chdir(old)
    return jobs


def PseudoNotebookRun(notebook_python_file, target_object, chdir=False):
    notebook_python_file = str(notebook_python_file)
    inv = ppg.FileInvariant(notebook_python_file)

    def run():
        import marburg_biobank.create

        source = Path(notebook_python_file).read_text()
        collector = {}

        def write_dfs(d):
            res = {}
            for k, v in d.items():
                if isinstance(v, tuple):
                    collector[k] = v[0]  # throw away description
                else:
                    collector[k] = v
            return res

        def get_dummy_ipython():
            class DummyIpython:
                def run_line_magic(self, *args, **kwargs):
                    pass

            return DummyIpython()

        marburg_biobank.create.write_dfs = write_dfs
        g = globals().copy()
        g["get_ipython"] = get_dummy_ipython
        g['here'] = Path(notebook_python_file).parent.absolute()
        ppg.util.global_pipegraph = None
        if chdir:
            os.chdir(Path(notebook_python_file).parent)
        exec(source, g)
        os.chdir("/project")
        return collector

    return ppg.CachedAttributeLoadingJob(
        notebook_python_file + ".result", target_object, "data", run
    ).depends_on(inv)
