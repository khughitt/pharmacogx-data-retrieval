PharmacoGx Data Retrieval Pipeline
==================================

## Overview

This repository contains a simple
[Snakemake](https://snakemake.readthedocs.io/en/stable/) pipeline for retrieving and processing molecular
profiling and drug screen data available on
[PharmacoGx](https://pharmacodb.ca/pharmacogx).

Datasets specified in a config file are downloaded, cleaned, and stored as tabular data
files, for ease of incorporation into external pipelines and workflows.

Basic filtering, transformation, and imputation logic is also provided.

## Usage

To begin, create a copy of the repo and cd to the directory:

```
git clone https://github.com/khughitt/pharmacogx-data-retrieval
cd pharmacogx-data-retrieval
```

Next, create and activate a conda environment with the neccessary dependencies using:

```sh
conda create -n pharmacogx --file requirements.txt
conda activate pharmacogx
```

Finally, after editing the YAML configuration file to make any desired adjustments, run
the pipeline using the following command:

```
snakemake --configfile config/config-v1.0.yml -j<NUM THREADS>
```

Where `<NUM THREADS>` is the number of threads to allow the pipeline to use.


