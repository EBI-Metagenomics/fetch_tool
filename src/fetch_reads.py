#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright 2018 EMBL - European Bioinformatics Institute
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
import json
import logging
import os
import re
import subprocess
import sys
import urllib.request
from subprocess import call

from src import ENAAPIUtils as apiutils

__author__ = "Hubert Denise, Simon Potter, Maxim Scheremetjew"
__copyright__ = "Copyright (c) 2018 EMBL - European Bioinformatics Institute"
__version__ = "1.0-rc1"
__status__ = "Development"

"""
    Represents a Python 2 version of the raw data fetch script. This version does not support assemblies and private
    data.

    Program exit codes:
        0: Succeeded to fetch raw data
        1: Something failed. The warning message should tell you what.
        5: No raw data available yet in ENA.
"""

script_dir = os.path.dirname(os.path.abspath(__file__))
default_configfile_basename = os.path.join(script_dir, os.pardir, "fetchdata-config-default.json")

LATEST_PIPELINE_VERSION = '4.1'


class ENADataFetcher(object):
    def __init__(self, ssh_max_attempts, url_max_attempts, ena_root_path,
                 ena_login_hosts):
        self.ssh_max_attempts = ssh_max_attempts
        self.url_max_attempts = url_max_attempts
        self.ena_root_path = ena_root_path
        self.ena_login_hosts = ena_login_hosts
        self.ena_project_url = 'http://www.ebi.ac.uk/ena/data/warehouse/filereport?accession={0}&result=read_run&fields=study_accession,secondary_study_accession,sample_accession,secondary_sample_accession,experiment_accession,run_accession,tax_id,scientific_name,instrument_model,library_layout,fastq_ftp,submitted_ftp,library_strategy,library_source&download=txt'

    @staticmethod
    def _combine_analyses(metadata, wgs):
        combine_analyses = []
        wgs_dict = {}
        for assembly in wgs:
            analysis_id = assembly['ASSEMBLY_ID']
            wgs_dict[analysis_id] = {'CONTIG_CNT': assembly['CONTIG_CNT'],
                                     'FILENAME': assembly['WGS_ACC'],
                                     'GC_ID': assembly['GC_ID']}
        for analysis in metadata:
            analysis_id = analysis['ANALYSIS_ID']
            if analysis_id in wgs_dict:
                new_analysis = analysis.copy()
                new_analysis.update(wgs_dict[analysis_id])
                combine_analyses.append(new_analysis)
            else:
                logging.error(
                    'No ENA contig name for ' + analysis_id + ' yet. Exiting')
                sys.exit()
        return combine_analyses

    def _trigger_ftp_download(self, url, fn, interactive):
        attempt = 0
        while True:
            try:
                logging.info("Downloading file from FTP server..." + url)
                download_command = ["wget", "-q", "-t", "5", "-O", fn, url]
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
            if attempt >= self.url_max_attempts:
                logging.critical("Failed to retrieve" + url + " after " + str(
                    attempt) + " attempts")
                if interactive:
                    var = input(
                        "Please type C to continue to fetch the next sequence file or anything else to exit: ")
                    if not var.upper().startswith('C'):
                        logging.info("Exiting now")
                        sys.exit(0)
                    else:
                        break
                else:
                    if force:
                        logging.warning(
                            "Force mode is activated. Will skip the download of this run and move onto the next sequence!")
                        break
                    else:
                        logging.warning(
                            "Too many failed attempts. Program will exit now. " +
                            "Try again to fetch the data in interactive mode (-i option)!")
                        sys.exit(1)

    def _retrieve_url(self, url):
        attempt = 0
        response = None
        while True:
            try:
                response = urllib.request.urlopen(url)
                break
            except urllib.request.URLError:
                logging.warning("Error opening url " + url)
                attempt += 1
            if attempt >= self.url_max_attempts:
                logging.critical("Failed to open url " + url + " after " + str(
                    attempt) + " attempts")
                sys.exit(1)
        return response

    @staticmethod
    def _get_column_headers_list():
        return ['study_id', 'sample_id', 'run_id', 'library_layout', 'file',
                'file_path', 'tax_id',
                'scientific_name', 'library_strategy', 'library_source',
                'LATEST_PIPELINE_VERSION',
                'analysis_status', 'sample_biome', 'opt:assembly_id',
                'opt:analysis_id']

    @staticmethod
    def _get_broker_name(api_result, accession_id, broker_name_field):
        if accession_id in api_result:
            metadata_dict = api_result.get(accession_id)
            if broker_name_field in metadata_dict:
                return metadata_dict.get(broker_name_field)
            else:
                logging.warning(
                    "Could not retrieve the broker name for accession - " + accession_id + " - from the ENA api.")
        logging.warning(
            "Could not retrieve any metadata for accession - " + accession_id + " - from the ENA api.")
        return None

    @staticmethod
    def _get_file_name(file_path, run_id, counter, num_of_files):
        if "fasta" in file_path:
            return run_id + ".fasta.gz"
        elif "fastq" in file_path:
            paired_end_number_label = "_" + str(
                counter) if num_of_files > 1 else ""
            return run_id + paired_end_number_label + ".fastq.gz"
        else:
            logging.warning("Unknown sequence file format: " + file_path)
            sys.exit(1)

    def populate_download_file(self, downloadable_files, path_regex,
                               file_names, run_id, existing_download_rows,
                               new_entries, is_submitted_file=False,
                               file_dir='raw'):
        """

        :param is_submitted_file:
        :param run_id:
        :param downloadable_files:
        :param path_regex:
        :param dl:
        :param file_names:
        :param file_dir:
        :return:
        """
        file_path = ''
        for counter, downloadable_file in enumerate(downloadable_files, 1):
            path_match = path_regex.match(downloadable_file)
            if not path_match:
                logging.error(
                    "Failed to match file to obtain path " + downloadable_file + "\n")
                sys.exit(1)

            file_path = path_match.group(1)
            if is_submitted_file:
                file_name = self._get_file_name(file_path.lower(), run_id,
                                                counter,
                                                len(downloadable_files))
            else:
                file_name = path_match.group(2)
            file_names.append(file_name)
            row = '\t'.join([downloadable_file,
                             os.sep.join([file_dir, file_names[-1]])]) + '\n'
            if row not in existing_download_rows:
                new_entries.append(row)
        return file_path

    def _retrieve_project_info_ftp(self, acc, run_id_list, user_pass, api_url):
        data = self._retrieve_url(self.ena_project_url.format(acc))
        data = [s.decode().rstrip() for s in data]
        # check line count - must be min 2, inc header
        if len(data) < 2:
            return False
        path_re = re.compile('.*?/(.*)/(.*)')

        project_file = acc + '.txt'
        download_file = 'download'
        project_data, download_entries = [], []
        if os.path.isfile(project_file):
            with open(project_file, 'r') as fh:
                project_data = fh.readlines()
        if os.path.isfile(download_file):
            with open(download_file, 'r') as dl:
                download_entries = dl.readlines()
        new_download_entries = []
        new_project_rows = []

        headers_line = '\t'.join(self._get_column_headers_list()) + '\n'
        # Write header line if not present
        if len(project_data) == 0:
            new_project_rows.append(headers_line)

        for line in data[1:]:
            processed_files = []
            submitted_files = []
            fields = line.split('\t')
            if len(fields) > 10:
                fastq_ftp = fields[10]
                if fastq_ftp:
                    processed_files = fastq_ftp.split(';')
                else:
                    run_id = fields[5]
                    sample_id = fields[3]
                    submitted_ftp = fields[11]
                    if submitted_ftp:
                        logging.warning(
                            "Found a run - " + run_id + " - with submitted files only!")
                        if self._trusted_broker_check(sample_id, api_url,
                                                      user_pass,
                                                      result_type='sample'):
                            submitted_files = submitted_ftp.split(';')
                    else:
                        logging.warning(
                            "Found a run - " + run_id + " - without any processed or submitted files!")

                    if force:
                        logging.warning(
                            "Force mode is activated. Will skip the download of this run and move onto the next sequence file!")
            else:
                logging.error(
                    "Unexpected number of fields in downloaded TXT file!")
            run_id = fields[5]
            sample_id = fields[3]
            if not run_id_list or run_id in run_id_list:
                file_dir = 'raw'
                file_names = []
                if not os.path.isdir(file_dir):
                    os.mkdir(file_dir)
                file_path = ''
                if len(processed_files) > 0:
                    file_path = self.populate_download_file(processed_files,
                                                            path_re,
                                                            file_names, run_id,
                                                            download_entries,
                                                            new_download_entries)

                if len(submitted_files) > 0:
                    file_path = self.populate_download_file(submitted_files,
                                                            path_re,
                                                            file_names, run_id,
                                                            download_entries,
                                                            new_download_entries,
                                                            is_submitted_file=True)
                new_row = '\t'.join(
                    [fields[1], sample_id, fields[5], fields[9],
                     ';'.join(file_names), file_path,
                     fields[6],
                     fields[7], fields[12], fields[13], LATEST_PIPELINE_VERSION,
                     'COMPLETED',
                     'biome_placeholder', 'opt:assembly_id',
                     'opt:analysis_id']) + '\n'
                if new_row not in project_data:
                    new_project_rows.append(new_row)
            else:
                logging.debug(run_id + " is filtered out.")
        with open(project_file, 'a') as fh:
            fh.writelines(new_project_rows)
        with open('download', 'a') as dl:
            dl.writelines(new_download_entries)
        return True

    def _trusted_broker_check(self, accession, api_url, user_pass,
                              result_type='study'):
        """

        :type accession: str
        :type api_url: str
        :type user_pass: str
        :return:
        """
        if result_type not in ['study', 'sample']:
            logging.warning(
                "Unsupported result type specified: " + result_type)
        logging.info(
            "Checking if accession - " + accession + " - belongs to a trusted broker...")
        api_query = apiutils.create_study_search_batch_query([accession])
        accession_field = 'SECONDARY_STUDY_ACCESSION'
        if result_type == 'sample':
            api_query = apiutils.create_sample_search_batch_query([accession])
            accession_field = 'SECONDARY_SAMPLE_ACCESSION'

        broker_name_field = 'BROKER_NAME'
        search_fields = [broker_name_field, accession_field]
        response_str = apiutils.retrieve_metadata(result_type, api_query,
                                                  search_fields, user_pass,
                                                  api_url)
        processed_api_result = apiutils.parse_response_str(response_str,
                                                           accession_field)
        broker_name = self._get_broker_name(processed_api_result, accession,
                                            broker_name_field)
        if broker_name and broker_name in trusted_brokers:
            logging.info(
                "Trusted broker. Continue fetching submitted files...")
            return True
        else:
            logging.warning(
                "Untrusted broker for accession - " + accession)
            return False

    @staticmethod
    def _get_run_accessions(runs_list_dict):
        """

        :param runs_list_dict: List of dictionaries, containing run metadata.
        :return: List of run accessions
        """
        result = []
        for run_entry in runs_list_dict:
            run_acc = run_entry['RUN_ID']
            if run_acc and len(run_acc) > 0:
                result.append(run_acc)
        return result

    def _retrieve_project_info_db(self, acc, file, run_id_list, mode, api_url,
                                  user_pass, eradao):
        from ERADAO import ERADAO
        runs = ERADAO(eradao).retrieve_generated_files(acc)
        run_accession_list = self._get_run_accessions(runs)
        if self._trusted_broker_check(acc, api_url, user_pass):
            submitted_files = ERADAO(eradao).retrieve_submitted_files(acc)
            for submitted_file in submitted_files:
                run_id = submitted_file['RUN_ID']
                if run_id not in run_accession_list:
                    runs.append(submitted_file)
        if len(runs) < 1:
            return False
        with open(acc + '.txt', mode) as fh, open('download', "w") as dl:
            run_data = {}
            if mode == 'w':
                fh.write("\t".join(self._get_column_headers_list()) + "\n")
            path_re = re.compile('(.*)/(.*)')
            for run in runs:
                run_id = run['RUN_ID']
                sample_id = run['SAMPLE_ID']
                is_submitted_file = True if run[
                                                'DATA_FILE_ROLE'] == 'SUBMISSION_FILE' else False
                if not run_id_list or run_id in run_id_list:
                    path_match = path_re.match(run['DATA_FILE_PATH'])
                    if not path_match:
                        logging.error(
                            "Failed to match file to obtain path " + file + "\n")
                        sys.exit(1)
                    file_path = path_match.group(1)
                    file_dir = 'raw'
                    if not os.path.isdir(file_dir):
                        os.mkdir(file_dir)

                    if is_submitted_file:
                        if run_id in run_data:
                            file_name_1 = self._get_file_name(
                                file_path.lower(), run_id, 1, 2)
                            file_name_2 = self._get_file_name(
                                file_path.lower(), run_id, 2, 2)
                            file_name = file_name_1 + ';' + file_name_2
                        else:
                            file_name = self._get_file_name(file_path.lower(),
                                                            run_id, 1, 1)
                    else:
                        file_name = path_match.group(2)

                    if run_id in run_data:
                        if is_submitted_file:
                            run_data[run_id][3] = file_name
                        else:
                            run_data[run_id][3] += ';' + file_name

                    else:
                        run_data[run_id] = [run['STUDY_ID'], sample_id,
                                            run['LIBRARY_LAYOUT'], file_name,
                                            file_path,
                                            run['TAX_ID'], 'n/a',
                                            run['LIBRARY_STRATEGY'],
                                            run['LIBRARY_SOURCE'],
                                            LATEST_PIPELINE_VERSION, 'COMPLETED',
                                            'biome_placeholder',
                                            'opt:assembly_id',
                                            'opt: analysis_id']
                    dl.write('\t'.join([run['DATA_FILE_PATH'], os.sep.join(
                        [file_dir, file_name])]) + '\n')
                else:
                    logging.debug(run_id + " is filtered out.")
            for run_id in run_data:
                fh.write('\t'.join(
                    [run_data[run_id][0], run_data[run_id][1], run_id,
                     run_data[run_id][2], run_data[run_id][3],
                     run_data[run_id][4], run_data[run_id][5],
                     run_data[run_id][6], run_data[run_id][7],
                     run_data[run_id][8], run_data[run_id][9],
                     run_data[run_id][10], run_data[run_id][11], 'n/a',
                     'n/a']) + '\n')
        return True

    def download_project(self, acc, use_view, is_private, eradao, prod_user, run_id_list, interactive, user_pass,
                         api_url):
        summary_file = acc + ".txt"
        use_dcc_metagenome = False
        exit_tag = 0
        header = "Y"
        no_run_data_msg = "No run data for project " + acc
        prod_user_msg = "Using the production user for fetching data from dcc metagenome!"
        if use_view:
            use_dcc_metagenome = True
            if not self._retrieve_project_info_db(acc, summary_file,
                                                  run_id_list, 'w', api_url,
                                                  user_pass, eradao):
                logging.warning(no_run_data_msg)
                exit_tag += 1  # sys.exit(1)
            else:
                header = "N"
        else:
            if not self._retrieve_project_info_ftp(acc, run_id_list, user_pass,
                                                   api_url):
                logging.error(no_run_data_msg)
                exit_tag += 2
        if exit_tag >= 2:
            logging.warning("No data available for download!")
            sys.exit(5)

        logging.info(
            "Next step: Changing permission of the current working dir...")
        current_working_dir = os.getcwd()
        chmod_command = ['chmod', '-R', '775', current_working_dir]
        rv = subprocess.call(chmod_command)
        if rv:
            logging.error(
                "Could not change permission of the current working dir: " + current_working_dir)
            sys.exit(1)

        if os.path.isfile('download'):
            with open('download', 'r') as f:
                for line in f:
                    url, local_file = line.rstrip().split('\t')
                    if run_id_list:
                        run_id = os.path.basename(local_file).split('.')[0]
                        if run_id not in run_id_list:
                            continue
                    if use_dcc_metagenome and 'wgs' not in line:
                        self.download_fastq_dcc(url, local_file, prod_user,
                                                interactive)
                    else:
                        self.download_fastq_ftp(url, local_file, interactive)
        else:
            logging.warning("Nothing to download!")

    def download_fastq_dcc(self, path, fn, prod_user, interactive):
        logging.info("Retrieving data from /nfs/dcc_metagenome...")
        if not path[0] == '/' and 'wgs' not in path:
            path = os.sep.join([self.ena_root_path, path])
        current_working_dir = os.getcwd()

        fn = os.sep.join([current_working_dir, fn])
        attempt = 1
        n_hosts = len(self.ena_login_hosts)
        command = []
        while True:
            if not prod_user:
                command = ["sudo", "-H", "-u", "emgpr"]
            # command.extend(['ssh', '-o', 'ConnectTimeout=3', ena_login_hosts[attempt % n_hosts], scp_call])
            command.extend(['scp', '-o', 'StrictHostKeyChecking=no', '-o',
                            'ConnectTimeout=3',
                            self.ena_login_hosts[
                                attempt % n_hosts] + ":{}".format(path), fn])
            print(command)
            rv = call(command)
            if not rv:
                return
            if attempt >= self.ssh_max_attempts:
                logging.error('Failed to run ' + ' '.join(command))
                if interactive:
                    var = input(
                        "Please type C to continue to fetch the next sequence file or anything else to exit: ")
                    if not var.upper().startswith('C'):
                        logging.warning("Exiting now...")
                        sys.exit(0)
                    else:
                        break
                else:
                    if force:
                        logging.warning(
                            "Force mode is activated. Will skip the download of this run and move onto the next sequence!")
                        break
                    else:
                        logging.warning(
                            "Too many failed attempts. Program will exit now." +
                            " Try again to fetch the data in interactive mode (-i option)!")
                        sys.exit(1)
            attempt += 1

    def download_fastq_ftp(self, url, fn, interactive):
        if url[:3] == 'ftp':
            url = 'ftp://' + url
        self._trigger_ftp_download(url, fn, interactive)


def read_project_list(fn):
    plist = []
    try:
        with open(fn, 'r') as f:
            for line in f:
                plist.append(line.rstrip())
    except FileNotFoundError:
        logging.error("File \"" + fn + "\" not found ")
        sys.exit(1)
    return plist

def main():
    parser = argparse.ArgumentParser(
        description="Tool to fetch project sequence data from ENA")
    parser.add_argument("-p", "--project", help="Project accession(s)",
                        dest='projects', required=False,
                        action='append')
    parser.add_argument("-ru", "--runs",
                        help="Run accession(s) whitespace separated. That option is useful if you want to download only "
                             "certain project runs",
                        dest='project_runs',
                        nargs='+',
                        required=False)
    parser.add_argument("-c", "--config",
                        help="Configuration file [json]",
                        dest='config_file',
                        required=False)
    parser.add_argument("-d", "--dir", help="Base directory for downloads",
                        dest='ddir', required=True)
    parser.add_argument("-f", "--force",
                        help="Force the download, even if it hits cases of submitted files only. But it "
                             "will not download the submitted files.",
                        dest='force', required=False,
                        action='store_true')
    parser.add_argument("-w", "--use_view",
                        help="Use DB views rather than FTP", dest='use_view',
                        required=False,
                        action='store_true')
    parser.add_argument("-r", "--private", help="Data is private",
                        dest='is_private', required=False,
                        action='store_true')
    parser.add_argument("-s", "--skip",
                        help="Skip runs with missing/incorrect data",
                        dest='skip', required=False,
                        action='store_true')
    parser.add_argument("-l", "--project_list", help="Project list",
                        dest='plist', required=False)
    parser.add_argument("-v", "--verbose", help="Verbose", dest='verbosity',
                        required=False, action='count')
    parser.add_argument("-z", "--unzip", help="Unzip .gz files", dest='unzip',
                        required=False, action='store_true')
    parser.add_argument("-o", "--output", help="Output summary [json]",
                        dest='output_file', required=False)
    parser.add_argument("-i", "--interactive",
                        help="interactive mode - allows you to skip failed downloads.",
                        dest="interactive",
                        action="store_true",
                        required=False,
                        default=False)
    args = parser.parse_args()

    # Identifying user name
    user_name = os.environ.get('USER')
    prod_user = False
    if user_name == 'emgpr':
        prod_user = True

    projects = args.projects
    run_id_list = args.project_runs
    project_list = args.plist
    is_private = args.is_private
    use_view = args.use_view
    interactive = args.interactive
    if not use_view and is_private:
        use_view = True
    ddir = os.path.abspath(args.ddir)
    force = args.force
    verbosity = args.verbosity
    config_file = args.config_file
    skip = args.skip
    unzip = args.unzip
    output_file = args.output_file
    if output_file:
        if output_file[0] != '/':
            output_file = os.sep.join([os.getcwd(), output_file])

    if verbosity:
        if verbosity > 1:
            loglevel = logging.DEBUG
        else:
            loglevel = logging.INFO
    else:
        loglevel = logging.WARN

    logging.basicConfig(format='%(levelname)s %(asctime)s - %(message)s',
                        datefmt='%Y/%m/%d %I:%M:%S %p',
                        level=loglevel)

    if projects and project_list:
        logging.error("Only one argument -p or -l")
        sys.exit(1)
    if not (projects or project_list):
        logging.error("Need either argument -p or -l")
        sys.exit(1)

    if project_list:
        logging.info("Reading project accessions from " + project_list)
        projects = read_project_list(project_list)
    else:
        logging.info("Taking project accessions from command line")

    if not config_file:
        # try to load default config file, which is in the same location as the script itself
        config_file = default_configfile_basename
        if not os.path.exists(config_file):
            logging.error(
                "Configuration file with database parameters required")
            sys.exit(1)
        else:
            pass  # Default config files does exist and can be loaded

    eradao = None
    if use_view:
        from oracle_db_access_object import OracleDataAccessObject
        from oracle_db_connection import OracleDBConnection

        with open(config_file) as fh:
            config = json.load(fh)
            eradao = OracleDataAccessObject(
                OracleDBConnection(config["eraUser"], config["eraPassword"],
                                   config["eraHost"], config["eraPort"],
                                   config["eraInstance"]))


    with open(config_file) as fh:
        config = json.load(fh)
        ssh_max_attempts = config["ssh_max_attempts"]
        url_max_attempts = config["url_max_attempts"]
        ena_root_path = config["ena_root_path"]
        ena_login_hosts = config["ena_login_hosts"]
        program = ENADataFetcher(ssh_max_attempts, url_max_attempts,
                                 ena_root_path, ena_login_hosts)

        api_url, trusted_brokers = config["enaAPIUrl"], config[
            "trustedBrokers"]
        api_username, api_password = config["enaAPIUsername"], config[
            "enaAPIPassword"]
        if api_username and api_password:
            api_credentials = api_username + ":" + api_password
        else:
            api_credentials = ""

    for pacc in projects:
        # Create project directory
        pdir = os.sep.join([ddir, pacc])
        if not os.path.isdir(pdir):
            os.mkdir(pdir)
        os.chdir(pdir)
        logging.info("Handling project " + pacc)

        # Make a copy of the web uploader config file (a template version sleeps in the template sub folder)
        program.download_project(pacc, use_view, is_private, eradao, prod_user, run_id_list, interactive,
                                 api_credentials, api_url)

    if output_file:
        with open(output_file, 'w') as of:
            data = {}
            data['download_path'] = ddir
            data['config_file'] = os.path.abspath(config_file)
            of.write(json.dumps(data, sort_keys=True, indent=4) + '\n')



if __name__ == '__main__':
    main()
