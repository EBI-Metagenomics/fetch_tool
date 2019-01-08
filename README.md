Set of tools which allows you to fetch RAW read and assembly files from the European Nucleotide Archive (ENA).


Install fetch tool
============================

Using Conda
-----------

    $ conda create -q -n fetch_tool python=3.6.4
    $ source activate fetch_tool

Install from requirements file

    $ git clone git@github.com:EBI-Metagenomics/fetch_tool.git
    $ pip install -r requirements.txt

Install from Git repo

    $ pip install git+ssh://git@github.com/EBI-Metagenomics/fetch_tool.git


Fetch read files (amplicon and WGS data)
=====

Usage
-----


    $ fetch_reads.py -h
    usage: fetch_reads.py [-h] [-p PROJECTS] [-ru PROJECT_RUNS [PROJECT_RUNS ...]]
                      [-c CONFIG_FILE] -d DDIR [-f] [-w] [-r] [-s] [-l PLIST]
                      [-v] [-z] [-o OUTPUT_FILE] [-i]

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
      -s, --skip            Skip runs with missing/incorrect data
      -l PLIST, --project_list PLIST
                        Project list
      -v, --verbose         Verbose
      -z, --unzip           Unzip .gz files
      -o OUTPUT_FILE, --output OUTPUT_FILE
                        Output summary [json]
      -i, --interactive     interactive mode - allows you to skip failed
                        downloads.


Example
--------

Download amplicon study:

    $ fetch-read-tool -p SRP062869 -c fetchdata-config-local.json -v -d /home/<user>/temp/

Fetch assembly files
=====

Usage
-----


    $ fetch_reads.py -h
    usage: fetch-read-tool [-h] -p PROJECT [-c CONFIG_FILE] [-d DDIR]
                           [-s {ftp,filesystem}] [-v]
                           [-pr {1.0,2.0,3.0,4.0,4.1}] [-i]

    Tool to fetch assemblies from ENA

    optional arguments:
      -h, --help            show this help message and exit
      -p PROJECT, --project PROJECT
                        Project accession
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

Example
--------

Download assembly study:

    $ fetch-assembly-tool -p ERP111288 -c fetchdata-config-local.json -v -d /home/<user>/temp/