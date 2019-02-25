import hashlib
from abc import ABC, abstractmethod
import json
import os
import logging
import sys
import argparse

from filelock import SoftFileLock
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
    DEFAULT_HEADERS = None

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

        self.prod_user = os.environ.get('USER') == 'emgpr'

        self._process_additional_args()
        self.projects = self.args.projects
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
        project_args.add_argument('-p', '--projects', help='Project accession(s)', nargs='+')
        project_args.add_argument("-l", "--project-list", help="File containing line-separated project list")
        parser.add_argument('-d', '--dir', help='Base directory for downloads', default=os.getcwd())
        parser.add_argument('-v', '--verbose', help='Verbose', action='count')
        parser.add_argument('-f', '--force', help='Force', action='store_true')
        parser.add_argument('--private', help='Use when fetching private data', action='store_true')
        parser.add_argument('-i', '--interactive', help='interactive mode - allows you to skip failed downloads.',
                            action='store_true')
        parser.add_argument('-c', '--config-file', help='Alternative config file', default=config_file)
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
    def write_project_files(self, project_accession, new_runs):
        pass

    @abstractmethod
    def _retrieve_project_info_db(self, project_accession):
        pass

    @abstractmethod
    def _retrieve_project_info_ftp(self, project_accession):
        pass

    def fetch(self):
        for project_accession in self.projects:
            self.fetch_project(project_accession)

    def fetch_project(self, project_accession):
        new_runs = self.retrieve_project(project_accession)
        if len(new_runs) == 0 and not self.force_mode:
            logging.warning('No new data found')
            return
        project_accession = new_runs[0]['STUDY_ID']

        os.makedirs(self.get_project_workdir(project_accession), exist_ok=True)

        self.write_project_files(project_accession, new_runs)

        self.download_raw_files(project_accession, new_runs)

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
            download_sources = run['DATA_FILE_PATH'].split(';')
            file_md5s = run['MD5'].split(';')
            for dl_file, dl_name, dl_md5 in zip(download_sources, run['files'], file_md5s):
                dest = os.path.join(raw_dir, dl_name)
                self.download_raw_file(dl_file, dest, dl_md5)

    def download_raw_file(self, dl_file, dest, dl_md5):
        filename = os.path.basename(dest)

        if not self._is_file_valid(dest, dl_md5) or self.force_mode:
            silentremove(dest)
            if self.private_mode:
                self.download_dcc(dest, dl_file)
            else:
                self.download_ftp(dest, dl_file)
        else:
            logging.info('File {} already exists, skipping download'.format(filename))

        if not self._is_file_valid(dest, dl_md5):
            logging.error('MD5 of downloaded file {} does not match expected MD5'.format(filename))

    def get_project_workdir(self, project_accession):
        return os.path.join(self.base_dir, project_accession[0:7], project_accession)

    def get_project_rawdir(self, project_accession):
        return os.path.join(self.base_dir, project_accession[0:7], project_accession, 'raw')

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
            for file in run['file'].split(';'):
                row = '\t'.join([run['file_path'], os.path.join('raw', os.path.basename(file))]) + '\n'
                new_download_rows.append(row)

        download_file = self.get_project_download_file(project_accession)
        if not os.path.isfile(download_file):
            self.create_empty_file(download_file)

        download_file_lock = download_file + '.lock'
        dh_lock = SoftFileLock(download_file_lock)
        with dh_lock:
            existing_rows = set(self.read_download_data(project_accession))
            existing_rows = existing_rows.union(set(new_download_rows))
            with open(download_file, 'w+') as f:
                f.writelines(sorted(existing_rows))

    def get_project_filepath(self, project_accession):
        return os.path.join(self.get_project_workdir(project_accession), project_accession + '.txt')

    def read_project_description_file(self, project_accession):
        filepath = self.get_project_filepath(project_accession)
        return pd.read_csv(filepath, sep='\t')

    def write_project_description_file(self, project_accession, new_rows):
        project_file = self.get_project_filepath(project_accession)
        project_file_lock = project_file + '.lock'
        new_file = not os.path.isfile(project_file)
        if new_file:
            self.create_empty_file(project_file)
        fh_lock = SoftFileLock(project_file_lock)
        with fh_lock:
            # Fallback in case empty file exists
            if not new_file:
                try:
                    project_runs = self.read_project_description_file(project_accession)
                except EmptyDataError:
                    new_file = True

            if new_file:
                project_runs = pd.DataFrame(new_rows)
                headers = self.DEFAULT_HEADERS
            else:
                headers = list(project_runs.columns.values)
                project_runs = project_runs.append(new_rows, sort=True)
                project_runs = project_runs.drop_duplicates(subset='run_id')

            project_runs = project_runs.fillna('n/a').sort_values(by='run_id')
            project_runs.to_csv(project_file, sep='\t', index=False, columns=headers)

    def get_api_credentials(self):
        return self.config['enaAPIUsername'] + ':' + self.config['enaAPIPassword']

    def get_trusted_brokers(self):
        return self.config['trustedBrokers']

    def get_ena_api_url(self):
        return

    def _get_studies_brokers(self, study_accessions):
        headers = {'Accept': 'text/plain',
                   'Content-Type': 'application/x-www-form-urlencoded'}
        data = {
            'dataPortal': 'metagenome',
            'result': 'study',
            'query': "secondary_study_accession%3D" + "%20OR%20secondary_study_accession%3D".join(study_accessions),
            'fields': 'secondary_study_accession,broker_name',
            'format': 'json'
        }

        r = requests.post(self.config['enaAPIUrl'], headers=headers, data=data,
                          auth=(self.config['enaAPIUsername'], self.config['enaAPIPassword']))
        response_str = r.text

        if r.status_code != 200:
            raise ValueError(r.text)
        json_data = json.loads(response_str)
        return {d['secondary_study_accession']: d['broker_name'] for d in json_data}

    def _study_has_permitted_broker(self, study_accession):
        broker = self.study_brokers[study_accession]
        return broker == '' or broker in self.config['trustedBrokers']

    def _is_rawdata_filetype(self, filename):
        return any(x in filename for x in ['.fa', '.fna', '.fasta', '.fq', 'fastq'])

    def _filter_secondary_files(self, joined_file_names):
        return ";".join([f for f in joined_file_names.split(';') if self._is_rawdata_filetype(f)])

    def _get_raw_filenames(self, file_names, run_id, is_submitted_file):
        file_names = self._filter_secondary_files(file_names)
        if is_submitted_file:
            file_names = self._rename_raw_files(file_names, run_id)
        else:
            file_names = [os.path.basename(f) for f in file_names.split(';')]
        return file_names

    @staticmethod
    def _rename_raw_files(file_names, run_id):
        file_names = file_names.lower()
        if any(x in file_names for x in ['.fasta', '.fna']):
            filetype = ".fasta.gz"
        elif ".fastq" in file_names:
            filetype = ".fastq"
        else:
            raise ValueError("Unknown sequence file format: " + file_names)
        files = file_names.split(';')
        if len(files) == 1:
            return [run_id + filetype]
        else:
            return [run_id + '_' + str(i) + filetype for i, _ in enumerate(files)]

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
            if md5(dest) == file_md5:
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
                        logging.warning(
                            'Too many failed attempts. Program will exit now.' +
                            ' Try again to fetch the data in interactive mode (-i option)!')
                        sys.exit(1)
            attempt += 1

    def download_ftp(self, dest, url):
        if url[:4] == 'ftp.':
            url = 'ftp://' + url
        attempt = 0
        while True:
            try:
                logging.info("Downloading file from FTP server..." + url)
                download_command = ["wget", "-q", "-t", "5", "-O", dest, url]
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
