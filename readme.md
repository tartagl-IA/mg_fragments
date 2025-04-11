# 🧪 MG Fragments  
[![DOI](https://zenodo.org/badge/955623938.svg)](https://doi.org/10.5281/zenodo.15199144)

## Introduction

### 1. Environment Setup

Execute the following command to setup the Python environment (considering a `python 3.6+` version already installed)

```console
$ python -m venv .venv
$ source .venv/bin/activate
$ python -m pip install -r requirements.txt
```

### 2. Database Setup

1. Download ChEMBLdb SQL Lite from [here](https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/) (file is named like `chembl_35_sqlite.tar.gz`)
2. Extract the db into the folder `db_chembl`


### 3. Run the App

```console
$ streamlit run app/home.py
```

## Applications

### 1. 📊 Database Explorer
The app permits to import from ChEMBLdb into a browsable DB the desidered molecules filtered by Target ID.
The app gives also the possibility to remove the previously imported molecules.


### 2. 🔬 Molecule Explorer
The app permits the generation of the fragments and the subsequentially creation of the clusters.
The final centroids of the clusters can be stored into SDF files.


### 3. 🖼️ Molecule Viewer
The app allows to view the molecules stored into the SDF files.