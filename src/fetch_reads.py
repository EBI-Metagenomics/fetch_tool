import re
import logging
from operator import itemgetter
from src.ERADAO import ERADAO

from src.abstract_fetch import AbstractDataFetcher

path_re = re.compile(r'(.*)/(.*)')


class FetchReads(AbstractDataFetcher):
    ENA_PROJECT_URL = 'https://www.ebi.ac.uk/ena/portal/api/filereport?accession={0}&result=read_run&' \
                      'fields=study_accession,secondary_study_accession,sample_accession,secondary_sample_accession,' \
                      'experiment_accession,run_accession,instrument_model,library_layout,' \
                      'fastq_ftp,fastq_md5,submitted_ftp,submitted_md5,library_strategy,broker_name,library_source&' \
                      'download=true'

    def __init__(self, argv=None):
        self.runs = None
        self.ACCESSION_FIELD = 'RUN_ID'
        super().__init__(argv)

    @staticmethod
    def add_arguments(parser):
        runs_group = parser.add_mutually_exclusive_group()
        runs_group.add_argument("-ru", "--runs", nargs='+',
                                help="Run accession(s), whitespace separated. Use to download only certain project runs")
        runs_group.add_argument("--run-list", help='File containing line-separated run accessions')
        return parser

    def _validate_args(self):
        if not any([self.args.runs, self.args.run_list, self.args.projects, self.args.project_list]):
            raise ValueError('No data specified, please use -ru, --run-list, -p or --project-list')
        elif not self.args.private and not (self.args.projects or self.args.project_list):
            raise NotImplementedError('Fetching studies from runs via FTP is not supported due '
                                      'to performance issues, please use --private mode')

    def _process_additional_args(self):
        if self.args.run_list:
            self.runs = self._read_line_sep_file(self.args.run_list)
        else:
            self.runs = self.args.runs

        if not self.args.projects and not self.args.project_list:
            logging.info('Fetching projects from list of runs')
            self.args.projects = self._get_project_accessions_from_runs(self.runs)

    def _get_project_accessions_from_runs(self, runs):
        self.init_era_dao()
        # Get all generated study data from ENA
        study_records = ERADAO(self.eradao).retrieve_study_accessions_from_runs(runs)
        return [s['STUDY_ID'] for s in study_records]

    @staticmethod
    def _get_study_run_accessions(project_data):
        return list(filter(bool, [entry['RUN_ID'] for entry in project_data]))

    def _filter_accessions_from_args(self, run_data, run_accession_field):
        if self.runs:
            run_data = list(filter(lambda r: r[run_accession_field] in self.runs, run_data))
        return run_data

    def map_project_info_db_row(self, run):
        return {
            'study_id': run['STUDY_ID'],
            'sample_id': run['SAMPLE_ID'],
            'run_id': run['RUN_ID'],
            'library_layout': run['LIBRARY_LAYOUT'],
            'file': run['file'],
            'file_path': run['DATA_FILE_PATH'],
            'library_strategy': run['LIBRARY_STRATEGY'],
            'library_source': run['LIBRARY_SOURCE']
        }

    def _retrieve_era_submitted_data(self, project_accession):
        return ERADAO(self.eradao).retrieve_submitted_files(project_accession)

    def _retrieve_era_generated_data(self, project_accession):
        return ERADAO(self.eradao).retrieve_generated_files(project_accession)

    def _retrieve_project_info_db(self, project_accession):
        # Get all generated study data from ENA
        study_run_data = self._retrieve_era_generated_data(project_accession)
        # Add all runs which don't have generated data
#        if self._study_has_permitted_broker(project_accession):
        submitted_files = self._retrieve_era_submitted_data(project_accession)
        existing_accessions = self._get_study_run_accessions(study_run_data)
        study_run_data.extend([f for f in submitted_files if f['RUN_ID'] not in existing_accessions])
#        else:
#            broker = None if project_accession not in self.study_brokers else self.study_brokers.get(project_accession)
#            logging.debug(
#                'Study {} does not come from a trusted broker ({}). Submitted_ftp files will be ignored'.format(
#                    project_accession, broker))
        for data in study_run_data:
            is_submitted_file = data['DATA_FILE_ROLE'] == 'SUBMITTED_FILE'
            data['DATA_FILE_PATH'], data['file'], data['MD5'] = self._get_raw_filenames(data['DATA_FILE_PATH'],
                                                                                        data['MD5'],
                                                                                        data['RUN_ID'],
                                                                                        is_submitted_file)
        return study_run_data

#    @staticmethod
#    def is_trusted_ftp_data(data, trusted_brokers):
#        return data['fastq_ftp'] != '' or (data['submitted_ftp'] != '' and data['broker_name'] in trusted_brokers)

#    def _filter_ftp_broker_names(self, data):
#        trusted_brokers = self.config['trustedBrokers']
#        return [d for d in data if self.is_trusted_ftp_data(d, trusted_brokers)]

    def map_datafields_ftp_2_db(self, rundata):
        is_submitted_file = rundata['submitted_ftp'] is not ''
        rundata['STUDY_ID'] = rundata.pop('secondary_study_accession')
        rundata['SAMPLE_ID'] = rundata.pop('secondary_sample_accession')
        rundata['RUN_ID'] = rundata.pop('run_accession')
        rundata['DATA_FILE_ROLE'] = 'SUBMISSION_FILE' if is_submitted_file else 'GENERATED_FILE'
        file_paths = rundata.get('fastq_ftp') or rundata.get('submitted_ftp')
        md5s = rundata.get('fastq_md5') or rundata.get('submitted_md5')
        rundata['DATA_FILE_PATH'], rundata['file'], rundata['MD5'] = self._get_raw_filenames(file_paths,
                                                                                             md5s,
                                                                                             rundata['RUN_ID'],
                                                                                             is_submitted_file)
        for key in ('fastq_ftp', 'submitted_ftp', 'fastq_md5', 'submitted_md5'):
            del rundata[key]
        rundata['LIBRARY_STRATEGY'] = rundata.pop('library_strategy')
        rundata['LIBRARY_SOURCE'] = rundata.pop('library_source')
        rundata['LIBRARY_LAYOUT'] = rundata.pop('library_layout')
        return rundata

    def _filter_secondary_files(self, joined_file_names, md5s):
        filtered_file_names, filtered_md5s = super()._filter_secondary_files(joined_file_names, md5s)
        # Case to remove runs with 3 fastq files => keep only _1 and _2
        if len(filtered_file_names) == 3:
            keep_files = ['_' in f for f in filtered_file_names]
            indexes = [i for i, x in enumerate(keep_files) if x]
            filtered_file_names = itemgetter(*indexes)(filtered_file_names)
            filtered_md5s = itemgetter(*indexes)(filtered_md5s)

        return filtered_file_names, filtered_md5s

    def _retrieve_project_info_ftp(self, project_accession):
        data = self._retrieve_ena_url(self.ENA_PROJECT_URL.format(project_accession))
#        trusted_data = self._filter_ftp_broker_names(data)
        return list(map(self.map_datafields_ftp_2_db, data))


def main():
    data_fetcher = FetchReads()
    data_fetcher.fetch()


if __name__ == '__main__':
    main()
