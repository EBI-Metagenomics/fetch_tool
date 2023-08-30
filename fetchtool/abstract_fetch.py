#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright 2018-2022 EMBL - European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


import argparse
import copy
import ftplib
import hashlib
import json
import logging
import os
import re
import subprocess
import sys
from abc import ABC, abstractmethod
from importlib.metadata import version

import pandas as pd
import requests
from flufl.lock import Lock
from pandas.errors import EmptyDataError

from fetchtool.exceptions import ENAFetch204, ENAFetch401, ENAFetchFail

PRIVATE_ENA_FTP = "ftp.dcc-private.ebi.ac.uk"
PUBLIC_ENA_FTP = "ftp.ebi.ac.uk"


class AbstractDataFetcher(ABC):
    DEFAULT_HEADERS = [
        "study_id",
        "sample_id",
        "run_id",
        "analysis_id",
        "library_layout",
        "library_strategy",
        "library_source",
        "file",
        "file_path",
    ]
    ACCESSION_FIELD = None
    ACCESSION_REGEX = r"([EDS]R[RZS]\d+)"
    PROGRAM_EXIT_MSG = "Program will exit now!"
    NO_DATA_MSG = "No entries found!"

    def __init__(self, argv=sys.argv[1:]):
        self.args = self._parse_args(argv)
        self._validate_args()

        self.set_logging(self.args.verbose)
        self.create_output_dir(self.args.dir)
        self.base_dir = self.args.dir

        config_file = os.getenv("FETCH_TOOL_CONFIG", None)

        if not self.args.config_file and not config_file:
            raise ValueError("Missing configuration file. It shoud be provided using -c or setting the env variable $FETCH_TOOL_CONFIG")

        with open(self.args.config_file or config_file) as f:
            self.config = json.load(f)

        self.ENA_API_USER = self.config["ena_api_username"]
        self.ENA_API_PASSWORD = self.config["ena_api_password"]

        self.interactive_mode = self.args.interactive
        self.private_mode = self.args.private
        self.force_mode = self.args.force
        self.desc_file_only = self.args.fix_desc_file
        self.ignore_errors = self.args.ignore_errors

        self.prod_user = os.environ.get("USER") == "emgpr"

        self._process_additional_args()
        if self.args.projects or self.args.project_list:
            self.projects = self._get_project_accessions(self.args)
            self.sanity_check_project_accessions()

    @abstractmethod
    def _validate_args(self):
        pass

    @abstractmethod
    def _process_additional_args(self):
        pass

    def _get_project_accessions(self, args):
        projects = []
        if args.projects:
            projects = args.projects
        elif args.project_list:
            projects = self._read_line_sep_file(args.project_list)
        logging.debug("Found projects " + ", ".join(projects))
        return projects

    @staticmethod
    def _read_line_sep_file(filename):
        with open(filename) as fhandler:
            data = [line.strip() for line in fhandler.readlines()]
        return data

    def _parse_args(self, argv):
        parser = argparse.ArgumentParser()
        project_args = parser.add_mutually_exclusive_group()
        project_args.add_argument(
            "-p",
            "--projects",
            help="Whitespace separated list of project accession(s)",
            nargs="+",
        )
        project_args.add_argument("-l", "--project-list", help="File containing line-separated project list")
        parser.add_argument("-d", "--dir", help="Base directory for downloads", default=os.getcwd())
        parser.add_argument("-v", "--verbose", help="Verbose", action="count")
        parser.add_argument("--version", help="Version", action="version", version=version("fetch-tool"))
        parser.add_argument(
            "-f",
            "--force",
            help="Ignore download errors and force re-download all files",
            action="store_true",
        )
        parser.add_argument(
            "--ignore-errors",
            help="Ignore download errors and continue",
            action="store_true",
        )
        parser.add_argument("--private", help="Use when fetching private data", action="store_true")
        parser.add_argument(
            "-i",
            "--interactive",
            help="interactive mode - allows you to skip failed downloads.",
            action="store_true",
        )
        parser.add_argument("-c", "--config-file", required=False, help="Alternative config file")
        parser.add_argument(
            "--fix-desc-file",
            help="Fixed runs in project description file",
            action="store_true",
        )
        parser = self.add_arguments(parser)
        return parser.parse_args(argv)

    @staticmethod
    def add_arguments(parser):
        return parser

    @staticmethod
    def set_logging(level):
        if level:
            if level > 1:
                loglevel = logging.DEBUG
            else:
                loglevel = logging.INFO
        else:
            loglevel = logging.WARN

        logging.basicConfig(
            format="%(levelname)s %(asctime)s - %(message)s",
            datefmt="%Y/%m/%d %I:%M:%S %p",
            level=loglevel,
        )

    @staticmethod
    def create_output_dir(dirname):
        os.makedirs(dirname, exist_ok=True)

    @abstractmethod
    def _retrieve_project_info_from_api(self, project_accession):
        pass

    @abstractmethod
    def _filter_accessions_from_args(self, data, fieldname):
        pass

    def fetch(self):
        for project_accession in self.projects:
            self.fetch_project(project_accession)

    def filter_by_accessions(self, new_data):
        if not self.force_mode:
            new_data = self._filter_accessions_from_args(new_data, self.ACCESSION_FIELD)
        return new_data

    def fetch_project(self, project_accession):
        new_data = self.retrieve_project(project_accession)
        if not new_data:  # exit function if there is no data and skip to the next study
            return
        if not self.desc_file_only and not self.force_mode:
            logging.info("Filtering study entries...")
            logging.info("Number of entries before filtering: {}".format(len(new_data)))
            new_data = self.filter_by_accessions(new_data)
            logging.info("Number of entries after filtering: {}.".format(len(new_data)))
        if len(new_data) == 0:
            logging.warning(self.NO_DATA_MSG)
            return
        secondary_project_accession = project_accession

        os.makedirs(self.get_project_workdir(secondary_project_accession), exist_ok=True)

        self.write_project_files(secondary_project_accession, new_data)

        if not self.desc_file_only:
            self.download_raw_files(project_accession, new_data)

    def retrieve_project(self, project_accession):
        new_runs = self._retrieve_project_info_from_api(project_accession)
        return new_runs

    def download_raw_files(self, project_accession, new_runs):
        raw_dir = self.get_project_rawdir(project_accession)
        os.makedirs(raw_dir, exist_ok=True)
        for run in new_runs:
            download_sources = run["DATA_FILE_PATH"]
            filenames = run["file"]
            file_md5s = run["MD5"]
            for dl_file, dl_name in zip(download_sources, filenames):
                dest = os.path.join(raw_dir, dl_name)
                self.download_raw_file(dl_file, dest, file_md5s, self.private_mode)

    def download_raw_file(self, dl_file, dest, dl_md5s, is_public):
        """
        Returns true if file was re-downloaded
        """
        filename = os.path.basename(dest)
        file_downloaded = False
        if not self._is_file_valid(dest, dl_md5s) or self.force_mode:
            silentremove(dest)
            try:
                file_downloaded = self.download_aspera(dest, dl_file)
                if not file_downloaded:
                    logging.info("Aspera didn't work.. trying FTP with lftp")
                    file_downloaded = self.download_lftp(dest, dl_file)
                if not file_downloaded:
                    logging.info("FTP didn't worked... trying wget now...")
                    file_downloaded = self.download_wget(dest, dl_file)
            except Exception as e:
                if self.ignore_errors:
                    logging.warning(e)
                else:
                    raise e
        else:
            logging.info("File {} already exists and MD5 matches, skipping download".format(filename))

        if not self._is_file_valid(dest, dl_md5s):
            msg = "MD5 of downloaded file {} does not match expected MD5".format(filename)
            if self.ignore_errors:
                logging.error(msg)
            else:
                raise EnvironmentError(msg)

        return file_downloaded

    def get_project_workdir(self, project_accession):
        return os.path.join(self.base_dir, project_accession)

    def get_project_rawdir(self, project_accession):
        return os.path.join(self.base_dir, project_accession, "raw")

    def get_project_download_file(self, project_accession):
        return os.path.join(self.get_project_workdir(project_accession), "download")

    def get_project_insdc_txt_file(self, project_accession):
        return os.path.join(self.get_project_workdir(project_accession), project_accession + "insdc.txt")

    def read_download_data(self, project_accession):
        filepath = self.get_project_download_file(project_accession)
        with open(filepath) as f:
            return f.readlines()

    @staticmethod
    def create_empty_file(filepath):
        open(filepath, "a").close()

    def write_project_download_file(self, project_accession, new_rows):
        new_download_rows = []
        for run in new_rows:
            for file_path, file in zip(run["file_path"], run["file"]):
                row = file_path + "\t" + file + "\n"
                new_download_rows.append(row)

        download_file = self.get_project_download_file(project_accession)
        if not os.path.isfile(download_file):
            self.create_empty_file(download_file)

        lock_file = download_file + ".lock"
        # Lock for 1 minute
        with Lock(lock_file, lifetime=60):
            existing_rows = set(self.read_download_data(project_accession))
            existing_rows = existing_rows.union(set(new_download_rows))
            with open(download_file, "w+") as f:
                f.writelines(sorted(existing_rows))

    def get_project_filepath(self, project_accession):
        return os.path.join(self.get_project_workdir(project_accession), project_accession + ".txt")

    def read_project_description_file(self, project_accession):
        filepath = self.get_project_filepath(project_accession)
        return pd.read_csv(filepath, sep="\t")

    @staticmethod
    def clean_data_row(data):
        clean_data = copy.deepcopy(data)
        for field in ["file", "file_path"]:
            clean_data[field] = ";".join(clean_data[field])
        return clean_data

    def get_downloaded_raw_file_accessions(self, project_accession):
        raw_dir = self.get_project_rawdir(project_accession)
        try:
            files = filter(
                len,
                map(lambda r: re.findall(self.ACCESSION_REGEX, r), os.listdir(raw_dir)),
            )
            accessions = {f[0] for f in list(files)}
        except FileNotFoundError:
            accessions = set()
        return accessions

    def generate_expected_desc_data(self, project_accession, existing_data, project_data):
        accessions = self.get_downloaded_raw_file_accessions(project_accession)
        if "run_id" in existing_data:
            accessions = accessions.union(existing_data["run_id"].tolist())
        if "analysis_id" in existing_data:
            accessions = accessions.union(existing_data["analysis_id"].tolist())

        project_data = list(
            filter(
                lambda r: (r.get("run_id") or r["analysis_id"]) in accessions,
                project_data,
            )
        )
        return project_data

    @staticmethod
    def remove_project_desc_duplicates(df):
        return (
            df.assign(counts=df.count(axis=1))
            .sort_values(by=["run_id", "analysis_id"])
            .drop_duplicates(subset=["run_id", "analysis_id"], keep="last")
        )

    def add_missing_headers(self, df):
        for h in self.DEFAULT_HEADERS:
            if h not in df:
                df[h] = None
        return df

    def write_project_description_file(self, project_accession, new_rows):
        project_data = list(map(self.clean_data_row, new_rows))

        project_file = self.get_project_filepath(project_accession)

        lock_file = project_file + ".lock"
        # We lock the file for 2 minutes, should be enough
        with Lock(lock_file, lifetime=120):
            # Fallback in case empty file exists
            try:
                project_runs = self.read_project_description_file(project_accession)
            except (EmptyDataError, FileNotFoundError):
                project_runs = pd.DataFrame(project_data)

            headers = self.DEFAULT_HEADERS

            if self.desc_file_only:
                project_data = self.generate_expected_desc_data(project_accession, project_runs, project_data)

            project_runs = pd.concat([project_runs, pd.DataFrame(project_data)], sort=True)
            project_runs = self.add_missing_headers(project_runs)
            project_runs = self.remove_project_desc_duplicates(project_runs)

            project_runs = project_runs.fillna("n/a").sort_values(by=["run_id", "analysis_id"])
            project_runs.to_csv(project_file, sep="\t", index=False, columns=headers)

    def get_api_credentials(self):
        return self.config["ena_api_username"] + ":" + self.config["ena_api_password"]

    @staticmethod
    def _is_rawdata_filetype(filename):
        return any(x in filename for x in [".fa", ".fna", ".fasta", ".fq", "fastq"])

    def _filter_secondary_files(self, joined_file_names, md5s):
        file_names = joined_file_names.split(";")
        md5s = md5s.split(";")
        filename_md5s = zip(file_names, md5s)
        filtered_filename_md5s = [(f, md5) for f, md5 in filename_md5s if self._is_rawdata_filetype(f)]
        filtered_file_names, filtered_md5s = zip(*filtered_filename_md5s)
        return filtered_file_names, filtered_md5s

    def _get_raw_filenames(self, filepaths, md5s, run_id, is_submitted_file):
        """
            Rename file names if submitted files or if generated assemblies
        :param filepaths:
        :param md5s:
        :param run_id:
        :param is_submitted_file:
        :return:
        """
        filepaths, md5s = self._filter_secondary_files(filepaths, md5s)
        if is_submitted_file or (not is_submitted_file and run_id.startswith("ERZ")):
            file_names = self._rename_raw_files(filepaths, run_id)
        else:
            file_names = [os.path.basename(f) for f in filepaths]
        # print("file path:{}".format(filepaths))
        # print("file names:{}".format(file_names))
        return filepaths, file_names, md5s

    @staticmethod
    def _rename_raw_files(file_names, run_id):
        file_names = [f.lower() for f in file_names]
        if any([".fastq" in fn for fn in file_names]):
            filetype = ".fastq.gz"
        elif any(x in ";".join(file_names) for x in [".fasta", ".fna", ".fa"]):
            filetype = ".fasta.gz"
        else:
            raise ValueError("Unknown sequence file format: " + ",".join(file_names))
        if len(file_names) == 1:
            return [run_id + filetype]
        else:
            return [run_id + "_" + str(i + 1) + filetype for i, _ in enumerate(file_names)]

    def _retrieve_ena_url(self, url, raise_on_204=True):
        """Request json from ENA
        raise_on_204: raise ENAFetch204 if the response status code i 204
        """
        attempt = 0
        request_params = {"url": url}
        if self.private_mode:
            request_params["auth"] = (self.ENA_API_USER, self.ENA_API_PASSWORD)
        while attempt <= self.config["url_max_attempts"]:
            try:
                response = requests.get(**request_params)
                if response.status_code == 200:
                    return response.json()
                if response.status_code == 204:
                    if raise_on_204:
                        raise ENAFetch204("No Runs/Assemblies found. Check if study is metagenomic")
                    logging.info("Run/Assembly may not be metagenomic. Skipping...")
                    return
                elif response.status_code == 401:
                    raise ENAFetch401("Invalid Username or Password!")
                else:
                    logging.warning(
                        "Received the following unknown response code from the " "Portal API server:\n{}".format(response.status_code)
                    )
            except requests.exceptions.RequestException as e:
                logging.warning("Request exception. " "Exception:\n {}".format(e))
            attempt += 1

        error_message = "Failed to open url " + url + " after " + str(attempt) + " attempts. "

        raise ENAFetchFail(error_message)

    @staticmethod
    def _is_file_valid(dest, file_md5):
        if os.path.exists(dest):
            basename = os.path.basename(dest)
            if md5(dest) in file_md5:
                return True
            else:
                logging.info("File {} exists, but MD5 does not match".format(basename))
        return False

    def download_wget(self, dest, url, auth=True):
        """Download the files on the url using wget."""
        if url[:4] == "ftp.":
            url = "ftp://" + url
        attempt = 0
        while True:
            try:
                logging.info("Downloading file from FTP server..." + url)
                download_command = [
                    "wget",
                    "-v",
                    "--user={}".format(self.ENA_API_USER),
                    "--password={}".format(self.ENA_API_PASSWORD) if auth else "-q",
                    "-t",
                    "5",
                    "-O",
                    dest,
                    url,
                ]
                retcode = subprocess.call(download_command)
                if retcode:
                    logging.error("Error downloading the file from " + url)
                else:
                    logging.info("Done.")
                break
            except IOError as err:
                logging.error("Error downloading the file from " + url)
                logging.error(err)
                attempt += 1
            if attempt >= self.config["url_max_attempts"]:
                logging.critical("Failed to retrieve" + url + " after " + str(attempt) + " attempts")
                if self.interactive_mode:
                    var = input("Please type C to continue to fetch the next sequence file or anything else to exit: ")
                    if not var.upper().startswith("C"):
                        logging.info("Exiting now")
                        sys.exit(0)
                    else:
                        break
                else:
                    if self.force_mode:
                        logging.warning("Force mode is activated. Will skip the download of this run and move onto the next sequence!")
                        break
                    else:
                        logging.warning(
                            "Too many failed attempts. Program will exit now. "
                            + "Try again to fetch the data in interactive mode (-i option)!"
                        )
                        sys.exit(1)

    def download_lftp(self, dest, url):
        """Download from ENA FTP server.
        Usage example, to get file path and names from full FTP URL
        - url = ftp.sra.ebi.ac.uk/vol1/sequence/ERZ166/ERZ1669403/contig.fa.gz
        - path list = ['vol1', 'sequence', 'ERZ166', 'ERZ1669403']
        - path = vol1/sequence/ERZ166/ERZ1669403
        - filename = contig.fasta.gz
        """
        server = PRIVATE_ENA_FTP if self.private_mode else PUBLIC_ENA_FTP
        path_list = url.split("ebi.ac.uk/")[-1].split("/")[:-1]
        path = "/".join(path_list)
        file_name = os.path.basename(url)
        attempt = 0
        while attempt <= 3:
            try:
                with ftplib.FTP(server, timeout=300) as ftp:
                    logging.info("Downloading file from FTP server..." + url)
                    if self.private_mode:
                        logging.info("Logging in...")
                        ftp.login(self.ENA_API_USER, self.ENA_API_PASSWORD)
                    else:
                        logging.info("Logging as anonymous")
                        ftp.login()
                    ftp.cwd(path)
                    logging.info("Getting the file...")
                    # store with the same name
                    with open(dest, "wb") as output_file:
                        ftp.retrbinary("RETR " + file_name, output_file.write)
                    logging.info("File " + dest + " downloaded.")
                    return True
            except ftplib.all_errors as e:
                logging.error(e)
                attempt += 1
        else:
            return False

    def download_aspera(self, dest: str, url: str) -> bool:
        """Download using the aspera cli.
        Usage example, to get file path and names from full FTP URL
        - url = ftp.sra.ebi.ac.uk/vol1/sequence/ERZ166/ERZ1669403/contig.fa.gz
        """
        ASPERA_SERVER = self.config.get("aspera_server", "fasp.ebi.ac.uk")
        ASPERA_BIN = os.environ.get("ASPERA_BIN") or self.config.get("aspera_bin")
        # The cert is needed by the aspera cli tool (asperaweb_id_dsa.openssh) - which usually is in <installation>/cli/etc/"
        ASPERA_CERT = os.environ.get("ASPERA_CERT") or self.config.get("aspera_cert")
        if ASPERA_BIN is None or ASPERA_CERT is None:
            logging.error("Aspera needs the binary ('aspera_bin') and the cert ('aspera_cert') config values")
            return False

        ASPERA_PORT = self.config.get("aspera_port", 33001)
        ASPERA_ENA_PUBLIC_USER = self.config.get("aspera_ena_public_user", "era-fasp")

        path = "/".join(url.split("ebi.ac.uk/")[-1].split("/")[:-1])
        file_name = os.path.basename(url)

        aspera_user_host = f"{ASPERA_ENA_PUBLIC_USER}@{ASPERA_SERVER}:{path}/{file_name}"
        ascp_command = [
            ASPERA_BIN,
            "-l",
            "300m",
            "-P",
            str(ASPERA_PORT),
            "-i",
            ASPERA_CERT,
        ]
        if self.private_mode:
            # For private ones we need to remove the certificate
            os.environ["ASPERA_SCP_PASS"] = self.ENA_API_PASSWORD
            aspera_user_host = f"{self.ENA_API_USER}@{ASPERA_SERVER}:{path}/{file_name}"
            del ascp_command[-2:]

        ascp_command.extend([aspera_user_host, dest])

        # only for debug logging level
        if logging.DEBUG >= logging.root.level:
            # "-L-", # print logging info, useful in case it fails
            index = 5 if self.private_mode else 3
            ascp_command.insert(index, "-L-")

        try:
            logging.info("Downloading with Aspera")
            logging.info(" ".join(ascp_command))
            result = subprocess.run(ascp_command, capture_output=True, text=True, check=True)
            logging.info(result.stdout)
            logging.debug(result.stderr)
        except Exception as error:
            logging.exception(error)
            logging.error("Failed to download the files with aspera")
            return False

        return True

    @staticmethod
    def get_md5_file(filename):
        return filename + ".md5"

    def read_md5_file(self, filename):
        with open(self.get_md5_file(filename)) as f:
            try:
                return f.readlines()[0]
            except IndexError:
                return None

    def write_md5(self, filename):
        md5_dest = self.get_md5_file(filename)
        md5_val = md5(filename)
        with open(md5_dest, "w+") as f:
            f.write(md5_val)

    @abstractmethod
    def map_project_info_to_row(self, data):
        pass

    def write_project_files(self, project_accession, new_runs):
        new_run_rows = list(map(self.map_project_info_to_row, new_runs))
        self.write_project_description_file(project_accession, new_run_rows)
        if not self.desc_file_only:
            self.write_project_download_file(project_accession, new_run_rows)

    @staticmethod
    def is_study_accession(accession):
        study_accssion_re = r"([ESD]RP\d{6,})"
        match = re.match(study_accssion_re, accession)
        if match and len(match.group(0)) == len(accession):
            return True
        return False

    def sanity_check_project_accessions(self):
        for study_acc in self.projects:
            if not self.is_study_accession(study_acc):
                logging.error("Encountered an invalid study accession: {}".format(study_acc))
                logging.info("Program will exit now!")
                sys.exit(1)


def silentremove(filename):
    try:
        os.remove(filename)
    except OSError as e:
        if type(e) != FileNotFoundError:
            raise


def md5(fname):
    hash_md5 = hashlib.md5()
    with open(fname, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()
