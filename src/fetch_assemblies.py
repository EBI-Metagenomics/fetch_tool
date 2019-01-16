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
import collections
import json
import logging
import os
import sys
from pathlib import Path
from subprocess import call

import sh

from src.rename_fasta_header_util import rename_raw_file, \
    change_fasta_headers
from src.fetch_reads import LATEST_PIPELINE_VERSION

__author__ = "Maxim Scheremetjew"
__copyright__ = "Copyright (c) 2018 EMBL - European Bioinformatics Institute"

"""
    Represents a Python 2 version of the raw data fetch script.
    This version does not support assemblies and private data.

    Program exit codes:
        0: Succeeded to fetch raw data
        1: Something failed. The warning message should tell you what.
        5: No raw data available yet in ENA.
"""

ENAAPIConfig = collections.namedtuple('ENAAPIConfig',
                                      'enaAPIUrl, api_credentials, '
                                      'master_password, swagger_url')

ERAPRODatabaseConfig = collections.namedtuple('ERAPRODatabaseConfig',
                                              'host, port, instance, user, password')

ENADatabaseConfig = collections.namedtuple('ENADatabaseConfig',
                                           'host, port, instance, user, password')

DCCMetagenomeConfig = collections.namedtuple('DCCMetagenomeConfig',
                                             'ssh_max_attempts, url_max_attempts,ena_root_path, ena_login_hosts')


class ENADataFetcher(object):
    def __init__(self, project_acc, pipeline_version, output_dir,
                 dcc_meta_confg, prod_user, era_db_config, ena_db_config,
                 interactive, api_config, assembly_id_list):
        self.project_acc = project_acc
        self.assembly_id_list = assembly_id_list
        self.pipeline_version = pipeline_version
        self.output_dir = output_dir
        self.interactive = interactive
        self.dcc_meta_confg = dcc_meta_confg
        self.prod_user = prod_user
        self.era_db_config = era_db_config
        self.ena_db_config = ena_db_config
        self.ena_project_url = 'http://www.ebi.ac.uk/ena/data/warehouse/filereport?accession={0}&result=read_run&fields=study_accession,secondary_study_accession,sample_accession,secondary_sample_accession,experiment_accession,run_accession,tax_id,scientific_name,instrument_model,library_layout,fastq_ftp,submitted_ftp,library_strategy,library_source&download=txt'
        self.filename_to_analysis_id_lookup = {}
        self.api_config = api_config
        self.project_txt_file = os.path.join(output_dir, project_acc + '.txt')
        self.download_file = os.path.join(output_dir, 'download')

    @staticmethod
    def _uncompress_file(file_path, check_integrity=False):
        cmd = []
        if check_integrity:
            cmd.append("-t")
        cmd.append(file_path)
        sh.gunzip(cmd)

    @staticmethod
    def compress_file(file_path):
        sh.gzip([file_path])

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

    def _filter_analyses(self, analyses):
        if not self.assembly_id_list:
            return analyses
        result = []
        for analysis_dict in analyses:
            assembly_id = analysis_dict.get('ANALYSIS_ID')
            if assembly_id in self.assembly_id_list:
                result.append(analysis_dict)
        return result

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

    def populate_download_file(self, downloadable_files, path_regex, dl,
                               file_names, run_id, is_submitted_file=False,
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

            dl.write('\t'.join([downloadable_file, os.sep.join(
                [file_dir, file_names[-1]])]) + '\n')
        return file_path

    def _uncompress_raw_files(self):
        flag_file = os.path.join(self.output_dir, 'decompress-success')
        if not os.path.exists(flag_file):
            with open(self.download_file, 'r') as f:
                logging.info("Uncompressing RAW files...")
                for line in f:
                    url, local_file = line.rstrip().split('\t')
                    absolute_local_file_path = os.path.join(self.output_dir,
                                                            local_file)
                    success_flag_file = ''.join(
                        [absolute_local_file_path, '.success'])
                    if not os.path.exists(absolute_local_file_path):
                        if os.path.exists(success_flag_file):
                            logging.info(
                                "File {0} already decompressed successfully.".format(
                                    absolute_local_file_path))
                        else:
                            logging.error(
                                "The following assembly file does not exist "
                                "and has not be decompressed successfully:\n{0}\n"
                                "Please investigate! Program will exit now!".format(
                                    absolute_local_file_path))
                            sys.exit(1)
                    else:
                        decompressed_file = absolute_local_file_path.replace(
                            '.gz', '')
                        if os.path.exists(decompressed_file):
                            os.remove(decompressed_file)
                        self._uncompress_file(
                            file_path=absolute_local_file_path,
                            check_integrity=True)
                        self._uncompress_file(absolute_local_file_path)
                        Path(success_flag_file).touch()
                logging.info("Finished uncompressing RAW files.")
                Path(flag_file).touch()
        else:
            logging.info("Files already decompressed!")

    def _retrieve_assembly_info_db(self, eradao, enadao):
        logging.info("Reading assembly meta data from the ENA...")
        acc = self.project_acc
        from src.ERADAO import ERADAO
        metadata_analyses = ERADAO(eradao).retrieve_assembly_metadata(acc)
        if len(metadata_analyses) < 1:
            return False
        project_acc = metadata_analyses[0]['PROJECT_ID']
        # to deal wit the fact than some NCBI study do not have project_id
        if len(project_acc) < 3:
            project_acc = acc
        from src.ENADAO import ENADAO
        wgs_analyses = ENADAO(enadao).retrieve_assembly_data(project_acc)
        print(json.dumps(wgs_analyses, indent=4))
        if len(wgs_analyses) < 1:
            logging.error(
                "Failed to retrieve contigs and assembly data for study " + acc + "\n")
            return False
        analyses = self._combine_analyses(metadata_analyses, wgs_analyses)
        filtered_analyses = self._filter_analyses(analyses)

        existing_project_entries, existing_download_entries = [], []
        if os.path.isfile(self.project_txt_file):
            with open(self.project_txt_file, 'r') as fh:
                existing_project_entries = fh.readlines()
        if os.path.isfile(self.download_file):
            with open(self.download_file, 'r') as dl:
                existing_download_entries = dl.readlines()
        new_download_entries = []
        new_project_rows = []

        headers_line = "\t".join(self._get_column_headers_list()) + '\n'
        if len(existing_project_entries) == 0:
            new_project_rows.append(headers_line)

        for analysis in filtered_analyses:
            analysis_id = analysis['ANALYSIS_ID']
            assembly_id = analysis['GC_ID']
            fh_path = analysis['DATA_FILE_PATH']
            basename = os.path.basename(fh_path)
            self.filename_to_analysis_id_lookup[basename] = analysis_id
            path = os.path.join(self.dcc_meta_confg.ena_root_path, fh_path)
            if len(path) < 26:
                logging.error(
                    "Failed to get contig information to obtain path ")
                sys.exit(1)
            local_file_name = os.path.basename(fh_path)
            future_file_name = analysis_id + '.fasta.gz'
            # Deal with the download file entries
            new_download_entry = '\t'.join(
                [path, os.sep.join(['raw', local_file_name])]) + '\n'
            if new_download_entry not in existing_download_entries:
                new_download_entries.append(new_download_entry)

            # Deal with the project file entries
            new_project_data_entry = '\t'.join(
                [analysis['STUDY_ID'], analysis['SAMPLE_ID'], analysis_id,
                 'FASTA',
                 future_file_name, fh_path, analysis['TAX_ID'], 'n/a',
                 'ASSEMBLY', 'METAGENOMIC',
                 self.pipeline_version, 'COMPLETED', 'biome_placeholder',
                 assembly_id, analysis_id]) + '\n'
            if new_project_data_entry not in existing_project_entries:
                new_project_rows.append(new_project_data_entry)

        with open(self.project_txt_file, 'a') as fh:
            fh.writelines(new_project_rows)
        with open(self.download_file, 'a') as dl:
            dl.writelines(new_download_entries)

        return True

    def fetch_data_from_ftp(self):
        # TODO: Implement!
        sys.exit(1)

    def fetch_data_from_fs(self):

        era_db_config = self.era_db_config
        ena_db_config = self.ena_db_config
        acc = self.project_acc

        from src.oracle_db_access_object import OracleDataAccessObject
        from src.oracle_db_connection import OracleDBConnection

        eradao = OracleDataAccessObject(
            OracleDBConnection(era_db_config.user, era_db_config.password,
                               era_db_config.host, era_db_config.port,
                               era_db_config.instance))

        enadao = OracleDataAccessObject(
            OracleDBConnection(ena_db_config.user, ena_db_config.password,
                               ena_db_config.host, ena_db_config.port,
                               ena_db_config.instance))

        exit_tag = 0
        no_assembly_data_msg = "No assembly data for project " + acc
        found_assembly_data_msg = "Found assembly data for project " + acc
        prod_user_msg = "Using the production user for fetching data from dcc metagenome!"

        if not self._retrieve_assembly_info_db(eradao, enadao):
            logging.warning(no_assembly_data_msg)
            exit_tag += 1  # sys.exit(1)
        else:
            logging.info(found_assembly_data_msg)
            logging.info(prod_user_msg)
        if exit_tag >= 1:
            logging.warning("No assembly data available for download!")
            sys.exit(5)

        self._check_if_file_exists(self.download_file)

        self._download_raw_files_from_dcc()

        self._uncompress_raw_files()

        webin_account_id = self._retrieve_webin_account_id(eradao)
        self._change_fasta_headers(webin_account_id)

        self._delete_original_files()

        self._compress_new_files()

        self._chmod_dir(self.output_dir)

    def _download_raw_files_from_dcc(self):
        flag_file = os.path.join(self.output_dir, 'download-success')
        if not os.path.exists(flag_file):
            with open(self.download_file, 'r') as f:
                logging.info("Copying RAW data from /nfs/dcc_metagenome...")
                for line in f:
                    url, local_file = line.rstrip().split('\t')
                    self.download_fastq_dcc(url, local_file)
                logging.info("Finished RAW data download.")
                Path(flag_file).touch()
        else:
            logging.info("Files already downloaded!")

    def download_fastq_dcc(self, path, fn):
        dcc_meta_confg = self.dcc_meta_confg

        fn_absolute_path = os.sep.join([self.output_dir, fn])
        attempt = 1
        n_hosts = len(dcc_meta_confg.ena_login_hosts)
        command = []
        if os.path.exists(fn_absolute_path):
            try:
                logging.info("Checking integrity of file:\n{0}".format(
                    fn_absolute_path))
                self._uncompress_file(fn_absolute_path, True)
                logging.info("File already downloaded successfully")
                return
            except:
                logging.warning(
                    "File is corrupted! Re-downloading file again...")
                os.remove(fn_absolute_path)

        while True:
            if not self.prod_user:
                command = ["sudo", "-H", "-u", "emgpr"]
            command.extend(['scp', '-o', 'StrictHostKeyChecking=no', '-o',
                            'ConnectTimeout=3', dcc_meta_confg.ena_login_hosts[
                                attempt % n_hosts] + ":{}".format(path),
                            fn_absolute_path])
            print(command)
            rv = call(command)
            if not rv:
                return
            if attempt >= dcc_meta_confg.ssh_max_attempts:
                logging.error('Failed to run ' + ' '.join(command))
                if self.interactive:
                    var = input(
                        "Please type C to continue to fetch the next sequence file or anything else to exit: ")
                    if not var.upper().startswith('C'):
                        logging.warning("Exiting now...")
                        sys.exit(0)
                    else:
                        break
                else:
                    logging.warning(
                        "Too many failed attempts. Program will exit now."
                        " Try again to fetch the data in interactive mode "
                        "(-i option)!")
                    sys.exit(1)
            attempt += 1

    @staticmethod
    def _check_if_file_exists(file_path):
        if not os.path.isfile(file_path):
            logging.warning(
                "The following file does not exist:\n" + file_path)
            logging.warning("Program will exit now!")
            sys.exit(1)

    def _rename_raw_files(self, file_path, output_dir):
        with open(file_path, 'r') as f:
            logging.info("Renaming RAW files...")
            for line in f:
                url, local_file = line.rstrip().split('\t')
                basename = os.path.basename(local_file)
                analysis_id = self.filename_to_analysis_id_lookup.get(basename)
                source = os.path.join(output_dir, local_file)
                destination = os.path.join(output_dir, 'raw',
                                           analysis_id + '.fasta')
                rename_raw_file(source=source, destination=destination)
            logging.info("Finished renaming RAW files.")

    def _change_fasta_headers(self, webin_account_id):
        flag_file = os.path.join(self.output_dir, 'fasta-header-success')
        if not os.path.exists(flag_file):
            with open(self.download_file, 'r') as f:
                logging.info("Changing FASTA headers...")
                for line in f:
                    url, local_file = line.rstrip().split('\t')
                    basename = os.path.basename(local_file)
                    analysis_id = self.filename_to_analysis_id_lookup.get(
                        basename)
                    source = os.path.join(self.output_dir, local_file)
                    source = source.replace('.gz', '')
                    destination = os.path.join(self.output_dir, 'raw',
                                               analysis_id + '.fasta')
                    swagger_credentials = (
                        webin_account_id, self.api_config.master_password)
                    change_fasta_headers(source,
                                         destination,
                                         analysis_id,
                                         self.api_config.swagger_url,
                                         swagger_credentials,
                                         self.api_config.enaAPIUrl,
                                         self.api_config.api_credentials)
                logging.info("Finished changing FASTA headers.")
                Path(flag_file).touch()
        else:
            logging.info("FASTA headers already changed!")

    def _retrieve_webin_account_id(self, eradao):
        from src.ERADAO import ERADAO
        results = ERADAO(eradao).retrieve_webin_accound_id(
            self.project_acc)
        return results[0]['SUBMISSION_ACCOUNT_ID']

    def _delete_original_files(self):
        with open(self.download_file, 'r') as f:
            logging.info("Deleting original files...")
            for line in f:
                url, local_file = line.rstrip().split('\t')
                local_file = local_file.replace('.gz', '')
                absolute_local_file_path = os.path.join(self.output_dir,
                                                        local_file)
                if os.path.exists(absolute_local_file_path):
                    os.remove(absolute_local_file_path)
            logging.info("Finished file deletion.")

    def _compress_new_files(self):
        flag_file = os.path.join(self.output_dir, 'compress-success')
        if not os.path.exists(flag_file):
            with open(self.project_txt_file, 'r') as f:
                logging.info("Compressing new files...")
                for line in f:
                    if not line.startswith('study_id'):
                        row = line.rstrip().split('\t')
                        local_file_name = ''.join([row[2], '.fasta'])
                        absolute_local_file_path = os.path.join(
                            self.output_dir,
                            'raw',
                            local_file_name)
                        self.compress_file(file_path=absolute_local_file_path)
                logging.info("Finished file compression.")
                Path(flag_file).touch()
        else:
            logging.info("Files already compressed!")

    @staticmethod
    def _chmod_dir(directory):
        logging.info("Changing permission of {0}".format(directory))
        sh.chmod(['774', '-R', directory])


def create_folder(directory):
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError as err:
        print("OS error: {0}".format(err))
    except:
        print("Unexpected error:", sys.exc_info()[0])
        raise


def load_configuration(config_file):
    try:
        with open(config_file) as fh:
            config = json.load(fh)
            # Load ENA ssh details
            ssh_max_attempts = config["ssh_max_attempts"]
            url_max_attempts = config["url_max_attempts"]
            ena_root_path = config["ena_root_path"]
            ena_login_hosts = config["ena_login_hosts"]
            dcc_meta_confg = DCCMetagenomeConfig(
                ssh_max_attempts=ssh_max_attempts,
                url_max_attempts=url_max_attempts,
                ena_root_path=ena_root_path,
                ena_login_hosts=ena_login_hosts)
            # Load ENA API details
            api_url, trusted_brokers = config["enaAPIUrl"], config[
                "trustedBrokers"]
            api_username, api_password, \
            master_password, swagger_url = config["enaAPIUsername"], \
                                           config["enaAPIPassword"], \
                                           config["enaMasterPassword"], \
                                           config["swagger_url"]
            if api_username and api_password:
                api_credentials = (api_username, api_password)
            else:
                api_credentials = ()

            api_config = ENAAPIConfig(enaAPIUrl=api_url,
                                      api_credentials=api_credentials,
                                      master_password=master_password,
                                      swagger_url=swagger_url)

            # Load ENA databases connection details
            eraUser, eraPassword, eraHost, = config["eraUser"], \
                                             config["eraPassword"], \
                                             config["eraHost"]
            eraPort, eraInstance = config["eraPort"], config["eraInstance"]
            era_db_config = ERAPRODatabaseConfig(host=eraHost,
                                                 port=eraPort,
                                                 instance=eraInstance,
                                                 user=eraUser,
                                                 password=eraPassword)
            enaUser, enaPassword, enaHost, = config["enaUser"], \
                                             config["enaPassword"], \
                                             config["enaHost"]
            enaPort, enaInstance = config["enaPort"], config["enaInstance"]
            ena_db_config = ENADatabaseConfig(host=enaHost,
                                              port=enaPort,
                                              instance=enaInstance,
                                              user=enaUser,
                                              password=enaPassword)

            return api_config, era_db_config, ena_db_config, dcc_meta_confg

    except OSError as err:
        print("OS error: {0}".format(err))
    except:
        print("Unexpected error:", sys.exc_info()[0])
        raise


def create_output_folders(ddir, project_acc):
    output_dir = os.path.join(ddir, project_acc)
    if not os.path.exists(output_dir):
        sh.mkdir(['-m', '774', output_dir])
    raw_dir = os.path.join(output_dir, 'raw')
    if not os.path.exists(raw_dir):
        sh.mkdir(['-m', '774', raw_dir])
    return output_dir


def main():
    parser = argparse.ArgumentParser(
        description="Tool to fetch assemblies from ENA")
    parser.add_argument("-p", "--project", help="Project accession",
                        dest='project', required=True)
    parser.add_argument("-as", "--assemblies",
                        help="Analysis accession(s) (e.g. ERZ773283) whitespace separated. That option is useful if you want to download only "
                             "certain project analyses",
                        dest='project_assemblies',
                        nargs='+',
                        required=False)
    parser.add_argument("-c", "--config",
                        help="Configuration file [json]",
                        dest='config_file',
                        default="fetchdata-config-default.json",
                        required=False)
    parser.add_argument("-d", "--dir", help="Base directory for downloads",
                        dest='ddir', default=os.getcwd())
    parser.add_argument("-s", "--source", help="Source of the RAW files.",
                        dest='source', required=False,
                        choices=['ftp', 'filesystem'], default='filesystem')
    parser.add_argument("-v", "--verbose", help="Verbose",
                        dest='verbosity',
                        required=False, action='count')
    parser.add_argument("-pr", "--pipeline-version",
                        help="Specify pipeline version e.g. 4.1",
                        dest='pipeline_version',
                        choices=['1.0', '2.0', '3.0', '4.0', '4.1'],
                        required=False, default=LATEST_PIPELINE_VERSION)
    parser.add_argument("-i", "--interactive",
                        help="interactive mode - allows you to skip failed downloads.",
                        dest="interactive",
                        action="store_true",
                        required=False,
                        default=False)
    args = parser.parse_args()

    # Setting up log levels
    verbosity = args.verbosity
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

    logging.info("Initialising the program...")
    # Identifying user name
    user_name = os.environ.get('USER')
    prod_user = False
    if user_name == 'emgpr':
        prod_user = True
    logging.debug(
        "The program is executed by the production user emgpr: " + str(
            prod_user))

    project_acc = args.project
    assembly_id_list = args.project_assemblies
    source = args.source
    config_file = args.config_file
    pipeline_version = args.pipeline_version
    interactive = args.interactive

    # Creating output folders
    output_dir = create_output_folders(args.ddir, project_acc)

    logging.debug("Loading configurations...")
    api_config, era_db_config, ena_db_config, dcc_meta_confg = load_configuration(
        config_file)

    program = ENADataFetcher(project_acc, pipeline_version, output_dir,
                             dcc_meta_confg, prod_user, era_db_config,
                             ena_db_config, interactive, api_config,
                             assembly_id_list)

    logging.info("Starting the program...")
    if source == 'filesystem':
        program.fetch_data_from_fs()
    else:
        program.fetch_data_from_ftp()
    logging.info("Program finished.")


if __name__ == '__main__':
    main()
