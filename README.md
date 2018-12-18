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

    $ pip install git+ssh://github.com/EBI-Metagenomics/fetch_tool.git


Usage
=====

    $ fetch_assemblies.py -h
    usage: fetch_assemblies.py [-h] -p PROJECT [-c CONFIG_FILE] [-d DDIR]
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


Examples
========

Download assembly study:

    $ fetch_assemblies.py -p ERP022256 -c fetchdata-config-local.json -v
