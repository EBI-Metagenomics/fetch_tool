[![Testing](https://github.com/EBI-Metagenomics/fetch_tool/actions/workflows/test.yml/badge.svg)](https://github.com/EBI-Metagenomics/fetch_tool/actions/workflows/test.yml)
[![PyPI version](https://badge.fury.io/py/fetch-tool.svg)](https://badge.fury.io/py/fetch-tool)
[![Docker Repository on Quay](https://quay.io/repository/microbiome-informatics/fetch-tool/status "Docker Repository on Quay")](https://quay.io/repository/microbiome-informatics/fetch-tool)

# Microbiome Informatics ENA fetch tool

Set of tools which allows you to fetch RAW read and assembly files from the European Nucleotide Archive (ENA).

## Install fetch tool

Install from Pypi

```bash
$ pip install fetch-tool
```

Install from the git repo

```bash
$ pip install https://github.com/EBI-Metagenomics/fetch_tool/archive/master.zip
```

#### Configuration options

The tool has a number of options, with sensible defaults for the most common use cases.

Setup the configuration file, the template [fetchdata-config-template.json](config/fetchdata-config-template.json) for the configuration file.

The required fields are:
-
  - ena_api_user
  - ena_api_password

## Fetch read files (amplicon and WGS data)

### Usage

```bash
$ fetch-read-tool -h
usage: fetch-read-tool [-h] [-p PROJECTS [PROJECTS ...] | -l PROJECT_LIST] [-d DIR] [-v] [--version] [-f] [--ignore-errors] [--private] [-i] [-c CONFIG_FILE] [--fix-desc-file] [-ru RUNS [RUNS ...]
                       | --run-list RUN_LIST]

optional arguments:
  -h, --help            show this help message and exit
  -p PROJECTS [PROJECTS ...], --projects PROJECTS [PROJECTS ...]
                        Whitespace separated list of project accession(s)
  -l PROJECT_LIST, --project-list PROJECT_LIST
                        File containing line-separated project list
  -d DIR, --dir DIR     Base directory for downloads
  -v, --verbose         Verbose
  --version             Version
  -f, --force           Ignore download errors and force re-download all files
  --ignore-errors       Ignore download errors and continue
  --private             Use when fetching private data
  -i, --interactive     interactive mode - allows you to skip failed downloads.
  -c CONFIG_FILE, --config-file CONFIG_FILE
                        Alternative config file
  --fix-desc-file       Fixed runs in project description file
  -ru RUNS [RUNS ...], --runs RUNS [RUNS ...]
                        Run accession(s), whitespace separated. Use to download only certain project runs
  --run-list RUN_LIST   File containing line-separated run accessions
```

### Example

Download amplicon study:

```bash
$ fetch-read-tool -p SRP062869 -v -d /home/<user>/temp/
```

## Fetch assembly files

### Usage

```
fetch-assembly-tool -h
usage: fetch-assembly-tool [-h] [-p PROJECTS [PROJECTS ...] | -l PROJECT_LIST] [-d DIR] [-v] [--version] [-f] [--ignore-errors] [--private] [-i] [-c CONFIG_FILE] [--fix-desc-file]
                           [-as ASSEMBLIES [ASSEMBLIES ...]] [--assembly-type {primary metagenome,binned metagenome,metatranscriptome}] [--assembly-list ASSEMBLY_LIST]

optional arguments:
  -h, --help            show this help message and exit
  -p PROJECTS [PROJECTS ...], --projects PROJECTS [PROJECTS ...]
                        Whitespace separated list of project accession(s)
  -l PROJECT_LIST, --project-list PROJECT_LIST
                        File containing line-separated project list
  -d DIR, --dir DIR     Base directory for downloads
  -v, --verbose         Verbose
  --version             Version
  -f, --force           Ignore download errors and force re-download all files
  --ignore-errors       Ignore download errors and continue
  --private             Use when fetching private data
  -i, --interactive     interactive mode - allows you to skip failed downloads.
  -c CONFIG_FILE, --config-file CONFIG_FILE
                        Alternative config file
  --fix-desc-file       Fixed runs in project description file
  -as ASSEMBLIES [ASSEMBLIES ...], --assemblies ASSEMBLIES [ASSEMBLIES ...]
                        Assembly ERZ accession(s), whitespace separated. Use to download only certain project assemblies
  --assembly-type {primary metagenome,binned metagenome,metatranscriptome}
                        Assembly type
  --assembly-list ASSEMBLY_LIST
                        File containing line-separated assembly accessions
```

### Example

Download assembly study:

```bash
$ fetch-assembly-tool -p ERP111288 -v -d /home/<user>/temp/
```

# How to set up your development environment

We recommend you to use [miniconda|conda](https://docs.conda.io/en/latest/miniconda.html) to manage the environment.

Clone the repo and install the requirements.

```
$ git clone git@github.com:EBI-Metagenomics/fetch_tool.git
$ cd fetch_tool
$ # activate anv (conda activate xxx)
$ pip install .[dev]
```

## Pre-commit hooks

Setup the git [pre-commit hook](https://pre-commit.com/):

```bash
pre-commit install
```

*Why?*

pre-commit will run a set of pre-configured tools before allowing you to commit files. You can find the currently configure hooks and configurations in [.pre-commit-config.yaml](./.pre-commit-config.yaml)

## Tests

This repo uses [pytest](https://docs.pytest.org).
