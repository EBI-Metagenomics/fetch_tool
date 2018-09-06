Tool which allows you to fetch RAW read files from the European Nucleotide Archive (ENA).


Install fetch tool
============================

    pip install -U ??


Usage
=====

    $ fetch_data.py -h
    usage: fetch_data.py [-h] [-p PROJECTS] [-ru PROJECT_RUNS [PROJECT_RUNS ...]]
                     [-c CONFIG_FILE] -d DDIR [-f] [-w] [-r] [-s] [-l PLIST]
                     [-v] [-z] [-o OUTPUT_FILE] [-a] [-i]

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
      -a, --include_assembly
                            Data for project including assembly. FALSE means do
                            not include assemblies (default).
      -i, --interactive     interactive mode - allows you to skip failed
                            downloads.


Examples
========

Download metadata:

    $ fetch_data.py -v -p ERP000554 -d {output dir}
