[![Testing](https://github.com/EBI-Metagenomics/fetch_tool/actions/workflows/test.yml/badge.svg)](https://github.com/EBI-Metagenomics/fetch_tool/actions/workflows/test.yml)

# Microbiome Informatics ENA fetch tool

Set of tools which allows you to fetch RAW read and assembly files from the European Nucleotide Archive (ENA).

## How to set up your development environment

We recommend you to use [miniconda|conda](https://docs.conda.io/en/latest/miniconda.html) to manage the environment.

Clone the repo and install the requirements.

```
$ git clone git@github.com:EBI-Metagenomics/fetch_tool.git
$ cd fetch_tool
$ # activate anv (conda activate xxx)
$ pip install -r requirements-dev.txt
```

### Pre-commit hooks

Setup the git [pre-commit hook](https://pre-commit.com/):

```bash
pre-commit install
```

*Why?*

pre-commit will run a set of pre-configured tools before allowing you to commit files. You can find the currently configure hooks and configurations in [.pre-commit-config.yaml](./.pre-commit-config.yaml)

### Tests

This repo uses [pytest](https://docs.pytest.org).

It requires the aspera cli installed in the default location (`install-aspera.sh` with no parameters).

To run the test suite:
```bash
pytest
```

## Install fetch tool

### Using Conda

```bash
$ conda create -q -n fetch_tool python=3.8
$ conda activate fetch_tool
```

Install from requirements file

```bash
$ git clone git@github.com:EBI-Metagenomics/fetch_tool.git
$ pip install -r requirements.txt
$ pip install -U fetch_tool/
```

Install from Git repo

```bash
$ pip install git+ssh://git@github.com/EBI-Metagenomics/fetch_tool.git
```

Install from private Git repo with access token (access token can be found in centralised password file)

```bash
$ pip install -U git+https://{access_token}@github.com/EBI-Metagenomics/fetch_tool@master
```

#### Configuration file

Setup the configuration file, the template [config/fetchdata-config-template.json](fetchdata-config-template.json) for the configuration file.

The required fields are:
- For Aspera
  - aspera_bin (the path to ascp, usually in the aspera installation under /cli/bin)
  - aspera_cert (the path to the ascp provided cert, usually in the aspera installation under /cli/etc/asperaweb_id_dsa.openssh)
- To pull private ENA data
  - ena_api_user
  - ena_api_password

### Install Aspera

Taken from: http://docs.transfer.sdo.ebi.ac.uk/protocols/aspera/#endpoints

## Install

Run the `install-aspera.sh` command here, it has only one optional parameter (the installation folder).

```bash
./install path/to/installation-i-want
```

Otherwise it will install it in $PWD/aspera-cli

## Fetch read files (amplicon and WGS data)

### Usage

```bash
$ fetch-read-tool -h
usage: fetch-read-tool [-h] [-p PROJECTS]
                    [-ru PROJECT_RUNS [PROJECT_RUNS ...]] [-c CONFIG_FILE]
                    -d DDIR [-f] [-w] [-r] [-l PLIST] [-v] [-o OUTPUT_FILE]
                    [-i]

Tool to fetch project sequence data from ENA

optional arguments:
    -h, --help            show this help message and exit
    -p PROJECTS, --project PROJECTS
                        Project accession(s)
    -ru PROJECT_RUNS [PROJECT_RUNS ...], --runs PROJECT_RUNS [PROJECT_RUNS ...]
                        Run accession(s) whitespace separated. That option is
                        useful if you want to download only certain project
                        runs
    -c CONFIG_FILE, --config CONFIG_FILE
                        Configuration file [json]
    -d DDIR, --dir DDIR   Base directory for downloads
    -f, --force           Force the download, even if it hits cases of submitted
                        files only. But it will not download the submitted
                        files.
    -w, --use_view        Use DB views rather than FTP
    -r, --private         Data is private
    -l PLIST, --project_list PLIST
                        Project list
    -v, --verbose         Verbose
    -o OUTPUT_FILE, --output OUTPUT_FILE
                        Output summary [json]
    -i, --interactive     interactive mode - allows you to skip failed
                        downloads.
```

### Example


Download amplicon study:

    $ fetch-read-tool -p SRP062869 -c fetchdata-config-local.json -v -d /home/<user>/temp/

## Fetch assembly files


### Usage

    $ fetch-assembly-tool -h
    usage: fetch-assembly-tool [-h] -p PROJECT
                           [-as PROJECT_ASSEMBLIES [PROJECT_ASSEMBLIES ...]]
                           [-c CONFIG_FILE] [-d DDIR] [-s {ftp,filesystem}]
                           [-v] [-pr {1.0,2.0,3.0,4.0,4.1}] [-i]

    Tool to fetch assemblies from ENA

    optional arguments:
      -h, --help            show this help message and exit
      -p PROJECT, --project PROJECT
                            Project accession
      -as PROJECT_ASSEMBLIES [PROJECT_ASSEMBLIES ...], --assemblies PROJECT_ASSEMBLIES [PROJECT_ASSEMBLIES ...]
                            Analysis accession(s) (e.g. ERZ773283) whitespace
                            separated. That option is useful if you want to
                            download only certain project analyses
      -c CONFIG_FILE, --config CONFIG_FILE
                            Configuration file [json]
      -d DDIR, --dir DDIR   Base directory for downloads
      -s {ftp,filesystem}, --source {ftp,filesystem}
                            Source of the RAW files.
      -v, --verbose         Verbose
      -pr {1.0,2.0,3.0,4.0,4.1}, --pipeline-version {1.0,2.0,3.0,4.0,4.1}
                            Specify pipeline version e.g. 4.1
      -i, --interactive     interactive mode - allows you to skip failed
                            downloads.

### Example

Download assembly study:

```bash
$ fetch-assembly-tool -p ERP111288 -c fetchdata-config-local.json -v -d /home/<user>/temp/
```
