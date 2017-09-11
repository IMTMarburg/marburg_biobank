#marburg_biobank
##Introduction

The marburg_biobank python module offers a high level interface to the data sets
stored in the [Ovarian Cancer Effusion Biobank and Database])(https://www.ovara.net/biobank).

The basic usage is as follows:
```python
import marburg_biobank
db = marburg_biobank.OvcaBiobank("marburg_ovca_revision_5.zip") #  you need to download that file from your biobank.
print(db.list_datasets())
df_wide = db.get_wide('transcriptomics/rnaseq')  # to retrieve the data in a one sample per column / one row per measured variable format
df_tall = db.get_dataset('transcriptomics/rnaseq') # to retrieve the data in one row per data point format
```


##Data formats available
###wide
Using ```db.get_wide(dataset)```:

A pandas DataFrame that looks like this

| Index | Patient12, TAM | Patient12, TU | PatientX, Compartment
| ----- | --------------- | -------------- | ------------------------
| **VariableA, unitA** | 23.23 | 112.2 | nan |
| **VariableB, unitB** | 3.23 | 12.2 | 12.7 |


Caveats: If a dataset has only one compartment, the compartment information is ommited by get_wide(), unless .get_wide(standardized=True) is used.
The same applies for the unit in the index.
If there is a 'name' column in dataset, it get's added to the index, regardless of the value of standardized.

###tall
Using: ```db.get_dataset(dataset)```):

A pandas DataFrame that looks like this

|variable | unit | patient | compartment | value | optional columns...|
| ------- | ---- | ------- | ----------- | ----- | -----|
| variableA | unitA | Patient12 | TAM | 23.23|
| variableA | unitA | Patient12 | TU | 112.2|
| variableB | unitB | Patient13 | TAM | 3.23|
| variableB | unitB | Patient13 | TU | 12.2|

This is the internal storage format.


##compartments
 Compartments are an abstraction on top of 'cells' and 'bio-liquid'. Examples are Tumor associated macrophages (TAMs), Tumor cells (TU), ascites, blood...
 ```db.get_compartments()``` provides a list

##Datasets

Datasets are organized two levels deep. The first one defines the
*omics being measured (transcriptomics, proteomics, ... or 'clinical'), while
the second levels defines the actual method (RNaseq, FACS,...)

Survival data is in clinical/survival. Please remember: if using [https://pypi.python.org/pypi/lifelines](lifelines), censored and event are negations of each other.

##Excluded patients:
Patients are excluded from our studies on two levels.

* On global level (for example because their malignancy was not high grade serous ovarian carcinoma)
* On a per dataset level.

To query what patients are excluded use ```db.get_excluded_patients(dataset)```. Dataset may be an empty string, in which case you will receive
only the globally excluded patients.

```db.get_exclusion_reasons()``` Lists for each patient (and datasets) why they were excluded.
