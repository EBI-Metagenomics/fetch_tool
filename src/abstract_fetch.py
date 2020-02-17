import hashlib
from abc import ABC, abstractmethod
import json
import os
import logging
import sys
import argparse
import re

import copy
from filelock import UnixFileLock
from pandas.errors import EmptyDataError
import pandas as pd
import requests
import urllib.request
from subprocess import call

script_dir = os.path.dirname(os.path.abspath(__file__))

config_file = os.getenv('FETCH_TOOL_CONFIG',
                        os.path.realpath(os.path.join(script_dir, os.pardir,
                                                      "fetchdata-config-default.json")))


class AbstractDataFetcher(ABC):
    DEFAULT_HEADERS = ['study_id', 'sample_id', 'run_id', 'analysis_id', 'library_layout', 'library_strategy',
                       'library_source', 'file', 'file_path']
    ACCESSION_FIELD = None
    ACCESSION_REGEX = r'([EDS]R[RZS]\d+)'

    def __init__(self, argv=sys.argv[1:]):
        self.args = self._parse_args(argv)
        self._validate_args()

        self.set_logging(self.args.verbose)
        self.create_output_dir(self.args.dir)
        self.base_dir = self.args.dir

        self.config = load_config(self.args.config_file)

        self.interactive_mode = self.args.interactive
        self.private_mode = self.args.private
        self.force_mode = self.args.force
        self.desc_file_only = self.args.fix_desc_file
        self.ignore_errors = self.args.ignore_errors

        self.prod_user = os.environ.get('USER') == 'emgpr'

        self._process_additional_args()
        self.projects = self.args.projects
        self.sanity_check_project_accessions()
        self.project_accessions = self._get_project_accessions(self.args)

        if self.args.private:
            self.init_era_dao()
            self.init_ena_dao()
            self.study_brokers = self._get_studies_brokers(self.project_accessions)
        else:
            self.eradao = None
            self.enadao = None

    def init_era_dao(self):
        self.eradao = self.load_oracle_connection(self.config['eraUser'],
                                                  self.config['eraPassword'],
                                                  self.config['eraHost'],
                                                  self.config['eraPort'],
                                                  self.config['eraInstance'])

    def init_ena_dao(self):
        self.enadao = self.load_oracle_connection(self.config['enaUser'],
                                                  self.config['enaPassword'],
                                                  self.config['enaHost'],
                                                  self.config['enaPort'],
                                                  self.config['enaInstance'])

    @abstractmethod
    def _process_additional_args(self):
        pass

    def _get_project_accessions(self, args):
        projects = []
        if args.projects:
            projects = args.projects
        elif args.project_list:
            projects = self._read_line_sep_file(args.project_list)
        logging.debug('Found projects ' + ", ".join(projects))
        return projects

    @staticmethod
    def _read_line_sep_file(filename):
        with open(filename) as f:
            data = [l.strip() for l in f.readlines()]
        return data

    def _parse_args(self, argv):
        parser = argparse.ArgumentParser()
        project_args = parser.add_mutually_exclusive_group()
        project_args.add_argument('-p', '--projects', help='Whitespace separated list of project accession(s)',
                                  nargs='+')
        project_args.add_argument("-l", "--project-list", help="File containing line-separated project list")
        parser.add_argument('-d', '--dir', help='Base directory for downloads', default=os.getcwd())
        parser.add_argument('-v', '--verbose', help='Verbose', action='count')
        parser.add_argument('-f', '--force', help='Ignore download errors and force re-download all files',
                            action='store_true')
        parser.add_argument('--ignore-errors', help='Ignore download errors and continue', action='store_true')
        parser.add_argument('--private', help='Use when fetching private data', action='store_true')
        parser.add_argument('-i', '--interactive', help='interactive mode - allows you to skip failed downloads.',
                            action='store_true')
        parser.add_argument('-c', '--config-file', help='Alternative config file', default=config_file)
        parser.add_argument('--fix-desc-file', help='Fixed runs in project description file', action='store_true')
        parser = self.add_arguments(parser)
        return parser.parse_args(argv)

    def _validate_args(self):
        if not self.args.projects and not self.args.project_list:
            raise ValueError(
                'No projects specified, please set --projects <ERP... SRP...> or --project_list <projects.txt>')

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

        logging.basicConfig(format='%(levelname)s %(asctime)s - %(message)s',
                            datefmt='%Y/%m/%d %I:%M:%S %p',
                            level=loglevel)

    @staticmethod
    def create_output_dir(dirname):
        os.makedirs(dirname, exist_ok=True)

    @abstractmethod
    def _retrieve_project_info_db(self, project_accession):
        pass

    @abstractmethod
    def _retrieve_project_info_ftp(self, project_accession):
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
        if not self.desc_file_only and not self.force_mode:
            new_data = self.filter_by_accessions(new_data)
        if len(new_data) == 0:
            logging.warning('No new data found')
            return
        project_accession = new_data[0]['STUDY_ID']

        os.makedirs(self.get_project_workdir(project_accession), exist_ok=True)

        self.write_project_files(project_accession, new_data)

        if not self.desc_file_only:
            self.download_raw_files(project_accession, new_data)

    def retrieve_project(self, project_accession):
        if self.private_mode:
            new_runs = self._retrieve_project_info_db(project_accession)
        else:
            new_runs = self._retrieve_project_info_ftp(project_accession)
        return new_runs

    def download_raw_files(self, project_accession, new_runs):
        raw_dir = self.get_project_rawdir(project_accession)
        os.makedirs(raw_dir, exist_ok=True)
        for run in new_runs:
            download_sources = run['DATA_FILE_PATH']
            filenames = run['file']
            file_md5s = run['MD5']
            for dl_file, dl_name in zip(download_sources, filenames):
                dest = os.path.join(raw_dir, dl_name)
                self.download_raw_file(dl_file, dest, file_md5s)

    def download_raw_file(self, dl_file, dest, dl_md5s):
        """
            Returns true if file was re-downloaded
        """
        filename = os.path.basename(dest)
        file_downloaded = False
        if not self._is_file_valid(dest, dl_md5s) or self.force_mode:
            silentremove(dest)
            try:
                if self.private_mode:
                    self.download_dcc(dest, dl_file)
                else:
                    self.download_ftp(dest, dl_file)
                file_downloaded = True
            except Exception as e:
                if self.ignore_errors:
                    logging.warning(e)
                else:
                    raise e
        else:
            logging.info('File {} already exists, skipping download'.format(filename))

        if not self._is_file_valid(dest, dl_md5s):
            msg = 'MD5 of downloaded file {} does not match expected MD5'.format(filename)
            if self.ignore_errors:
                logging.error(msg)
            else:
                raise EnvironmentError(msg)

        return file_downloaded

    def get_project_workdir(self, project_accession):
        return os.path.join(self.base_dir, project_accession)

    def get_project_rawdir(self, project_accession):
        return os.path.join(self.base_dir, project_accession, 'raw')

    def get_project_download_file(self, project_accession):
        return os.path.join(self.get_project_workdir(project_accession), 'download')

    def read_download_data(self, project_accession):
        filepath = self.get_project_download_file(project_accession)
        with open(filepath) as f:
            return f.readlines()

    @staticmethod
    def create_empty_file(filepath):
        open(filepath, 'a').close()

    def write_project_download_file(self, project_accession, new_rows):
        new_download_rows = []
        for run in new_rows:
            for file_path, file in zip(run['file_path'], run['file']):
                row = file_path + '\t' + file + '\n'
                new_download_rows.append(row)

        download_file = self.get_project_download_file(project_accession)
        if not os.path.isfile(download_file):
            self.create_empty_file(download_file)

        lock_file = download_file + '.lock'
        with UnixFileLock(lock_file):
            existing_rows = set(self.read_download_data(project_accession))
            existing_rows = existing_rows.union(set(new_download_rows))
            with open(download_file, 'w+') as f:
                f.writelines(sorted(existing_rows))

    def get_project_filepath(self, project_accession):
        return os.path.join(self.get_project_workdir(project_accession), project_accession + '.txt')

    def read_project_description_file(self, project_accession):
        filepath = self.get_project_filepath(project_accession)
        return pd.read_csv(filepath, sep='\t')

    @staticmethod
    def check_files_downloaded(raw_dir, filenames):
        # Filter NaN values from pandas
        if not isinstance(filenames, str):
            return False
        filepaths = [os.path.join(raw_dir, f) for f in filenames.split(';')]
        return all(list(map(os.path.isfile, filepaths)))

    @staticmethod
    def clean_data_row(data):
        clean_data = copy.deepcopy(data)
        for field in ['file', 'file_path']:
            clean_data[field] = ";".join(clean_data[field])
        return clean_data

    def get_downloaded_raw_file_accessions(self, project_accession):
        raw_dir = self.get_project_rawdir(project_accession)
        try:
            files = filter(len, map(lambda r: re.findall(self.ACCESSION_REGEX, r), os.listdir(raw_dir)))
            accessions = {f[0] for f in list(files)}
        except FileNotFoundError:
            accessions = set()
        return accessions

    def generate_expected_desc_data(self, project_accession, existing_data, project_data):
        accessions = self.get_downloaded_raw_file_accessions(project_accession)
        if 'run_id' in existing_data:
            accessions = accessions.union(existing_data['run_id'].tolist())
        if 'analysis_id' in existing_data:
            accessions = accessions.union(existing_data['analysis_id'].tolist())

        project_data = list(filter(lambda r: (r.get('run_id') or r['analysis_id']) in accessions, project_data))
        return project_data

    @staticmethod
    def remove_project_desc_duplicates(df):
        return df.assign(counts=df.count(axis=1)) \
            .sort_values(by=['run_id', 'analysis_id']) \
            .drop_duplicates(subset=['run_id', 'analysis_id'], keep='last')

    def add_missing_headers(self, df):
        for h in self.DEFAULT_HEADERS:
            if h not in df:
                df[h] = None
        return df

    def write_project_description_file(self, project_accession, new_rows):
        project_data = list(map(self.clean_data_row, new_rows))

        project_file = self.get_project_filepath(project_accession)

        lock_file = project_file + '.lock'
        with UnixFileLock(lock_file):
            # Fallback in case empty file exists
            try:
                project_runs = self.read_project_description_file(project_accession)
            except (EmptyDataError, FileNotFoundError):
                project_runs = pd.DataFrame(project_data)
            headers = self.DEFAULT_HEADERS

            if self.desc_file_only:
                project_data = self.generate_expected_desc_data(project_accession, project_runs, project_data)
            project_runs = project_runs.append(project_data, sort=True)
            project_runs = self.add_missing_headers(project_runs)
            project_runs = self.remove_project_desc_duplicates(project_runs)

            project_runs = project_runs.fillna('n/a').sort_values(by=['run_id', 'analysis_id'])
            project_runs.to_csv(project_file, sep='\t', index=False, columns=headers)

    def get_api_credentials(self):
        return self.config['enaAPIUsername'] + ':' + self.config['enaAPIPassword']

    def get_trusted_brokers(self):
        return self.config['trustedBrokers']

    def get_ena_api_url(self):
        return

    def _get_studies_brokers(self, study_accessions):
        logging.info("Retrieving study broker name for study {}...".format(study_accessions))
        headers = {'Accept': '*/*',
                   'Content-Type': 'application/x-www-form-urlencoded'}
        data = {
            'dataPortal': 'metagenome',
            'result': 'study',
            'query': "secondary_study_accession%3D" + "%20OR%20secondary_study_accession%3D".join(study_accessions),
            'fields': 'secondary_study_accession,broker_name',
            'format': 'json'
        }

        r = requests.post(self.config['enaAPIUrl'] + 'search', headers=headers, data=data,
                          auth=(self.config['enaAPIUsername'], self.config['enaAPIPassword']))

        if r.status_code != 200:
            logging.warning(
                "Could NOT retrieve study broker name using the ENA Portal API! Got response code {}".format(
                    r.status_code))
            logging.warning("Possible reasons for that could be:\n - private study is"
                            " not indexed yet\n - study is not labelled as metagenome study ")
            return {acc: '' for acc in study_accessions}
        else:
            json_data = r.json()
            return {d['secondary_study_accession']: d['broker_name'] for d in json_data}

    def _study_has_permitted_broker(self, study_accession):
        broker = None if study_accession not in self.study_brokers else self.study_brokers.get(study_accession)
        return broker and broker != '' and broker in self.config['trustedBrokers']

    def _is_rawdata_filetype(self, filename):
        return any(x in filename for x in ['.fa', '.fna', '.fasta', '.fq', 'fastq'])

    def _filter_secondary_files(self, joined_file_names, md5s):
        file_names = joined_file_names.split(';')
        md5s = md5s.split(';')
        filename_md5s = zip(file_names, md5s)
        filtered_filename_md5s = [(f, md5) for f, md5 in filename_md5s if self._is_rawdata_filetype(f)]
        filtered_file_names, filtered_md5s = zip(*filtered_filename_md5s)
        return filtered_file_names, filtered_md5s

    def _get_raw_filenames(self, filepaths, md5s, run_id, is_submitted_file):
        filepaths, md5s = self._filter_secondary_files(filepaths, md5s)
        if is_submitted_file:
            file_names = self._rename_raw_files(filepaths, run_id)
        else:
            file_names = [os.path.basename(f) for f in filepaths]
        return filepaths, file_names, md5s

    @staticmethod
    def _rename_raw_files(file_names, run_id):
        file_names = [f.lower() for f in file_names]
        if any([".fastq" in fn for fn in file_names]):
            filetype = ".fastq.gz"
        elif any(x in ";".join(file_names) for x in ['.fasta', '.fna', '.fa']):
            filetype = ".fasta.gz"
        else:
            raise ValueError("Unknown sequence file format: " + ",".join(file_names))
        if len(file_names) == 1:
            return [run_id + filetype]
        else:
            return [run_id + '_' + str(i + 1) + filetype for i, _ in enumerate(file_names)]

    @staticmethod
    def load_oracle_connection(user, password, host, port, instance):
        from src.oracle_db_access_object import OracleDataAccessObject
        from src.oracle_db_connection import OracleDBConnection
        return OracleDataAccessObject(OracleDBConnection(user, password, host, port, instance))

    def _retrieve_ena_url(self, url):
        attempt = 0
        response = None
        while True:
            try:
                response = urllib.request.urlopen(url)
                break
            except urllib.request.URLError as e:
                logging.error(e)
                logging.warning("Error opening url " + url)
                attempt += 1
            if attempt >= self.config['url_max_attempts']:
                logging.critical("Failed to open url " + url + " after " + str(
                    attempt) + " attempts")
                sys.exit(1)
        data = [s.decode().rstrip() for s in response]
        headers = data[0].split('\t')
        data = data[1:]
        data = [{k: v for k, v in zip(headers, d.split('\t'))} for d in data]  # Convert rows to list of dictionaries
        return data

    @staticmethod
    def _is_file_valid(dest, file_md5):
        if os.path.exists(dest):
            basename = os.path.basename(dest)
            if md5(dest) in file_md5:
                return True
            else:
                logging.info('File {} exists, but MD5 does not match'.format(basename))
        return False

    def download_dcc(self, dest, download_file):
        attempt = 1
        n_hosts = len(self.config['ena_login_hosts'])
        command = []
        while True:
            if not self.prod_user:
                command = ['sudo', '-H', '-u', 'emgpr']
            host = self.config['ena_login_hosts'][attempt % n_hosts]
            host_path = '{}:/nfs/dcc_metagenomics/{}'.format(host, download_file)
            logging.info(host_path)
            command.extend(['scp', '-o', 'StrictHostKeyChecking=no', '-o', 'ConnectTimeout=3', host_path, dest])
            rv = call(command)
            if not rv:
                if not self.prod_user:
                    command = ["sudo", "-H", "-u", "emgpr"]
                    command.extend(['chmod', '775', dest])
                    rv = call(command)
                    if not rv:
                        return
                    else:
                        logging.warning('Failed to rename ' + str(dest))
                break
            if attempt >= self.config['ssh_max_attempts']:
                logging.error('Failed to run ' + ' '.join(command))
                if self.interactive_mode:
                    var = input(
                        'Please type C to continue to fetch the next sequence file or anything else to exit: ')
                    if not var.upper().startswith('C'):
                        logging.warning('Exiting now...')
                        sys.exit(0)
                    else:
                        break
                else:
                    if self.force_mode:
                        logging.warning('Force mode is activated. Will skip the download of this run '
                                        'and move onto the next sequence!')
                        break
                    else:
                        raise EnvironmentError('Too many failed attempts. Program will exit now.')
            attempt += 1

    def download_ftp(self, dest, url):
        if url[:4] == 'ftp.':
            url = 'ftp://' + url
        attempt = 0
        while True:
            try:
                logging.info("Downloading file from FTP server..." + url)
                download_command = ["wget", "-v" if self.args.verbose else "-q", "-t", "5", "-O", dest, url]
                retcode = call(download_command)
                if retcode:
                    logging.error("Error downloading the file from " + url)
                else:
                    logging.info("Done.")
                break
            except IOError as err:
                logging.error("Error downloading the file from " + url)
                logging.error(err)
                attempt += 1
            if attempt >= self.config['url_max_attempts']:
                logging.critical("Failed to retrieve" + url + " after " + str(
                    attempt) + " attempts")
                if self.interactive_mode:
                    var = input(
                        "Please type C to continue to fetch the next sequence file or anything else to exit: ")
                    if not var.upper().startswith('C'):
                        logging.info("Exiting now")
                        sys.exit(0)
                    else:
                        break
                else:
                    if self.force_mode:
                        logging.warning(
                            "Force mode is activated. Will skip the download of this run and move onto the next sequence!")
                        break
                    else:
                        logging.warning(
                            "Too many failed attempts. Program will exit now. " +
                            "Try again to fetch the data in interactive mode (-i option)!")
                        sys.exit(1)

    @staticmethod
    def get_md5_file(filename):
        return filename + '.md5'

    def read_md5_file(self, filename):
        with open(self.get_md5_file(filename)) as f:
            try:
                return f.readlines()[0]
            except IndexError:
                return None

    def write_md5(self, filename):
        md5_dest = self.get_md5_file(filename)
        md5_val = md5(filename)
        with open(md5_dest, 'w+') as f:
            f.write(md5_val)

    @abstractmethod
    def map_project_info_db_row(self, data):
        pass

    def write_project_files(self, project_accession, new_runs):
        new_run_rows = list(map(self.map_project_info_db_row, new_runs))
        self.write_project_description_file(project_accession, new_run_rows)
        if not self.desc_file_only:
            self.write_project_download_file(project_accession, new_run_rows)

    @staticmethod
    def is_study_accession(accession):
        study_accssion_re = r'([ESD]RP\d{6,})'
        matches = re.match(study_accssion_re, accession)
        if len(matches.groups()) == 1 and matches.end() == len(accession):
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
    except OSError as e:  # this would be "except OSError, e:" before Python 2.6
        if type(e) != FileNotFoundError:  # errno.ENOENT = no such file or directory
            raise  # re-raise exception if a different error occurred


def md5(fname):
    hash_md5 = hashlib.md5()
    with open(fname, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()


def load_config(config_file):
    with open(config_file) as f:
        return json.load(f)
