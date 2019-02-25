import re
import logging
import os
import sys
from src.ERADAO import ERADAO

from src.abstract_fetch import AbstractDataFetcher

path_re = re.compile(r'(.*)/(.*)')


class FetchAssemblies(AbstractDataFetcher):
    DEFAULT_HEADERS = ['study_id', 'sample_id', 'run_id', 'library_layout', 'file', 'file_path', 'tax_id',
                       'scientific_name', 'library_strategy', 'library_source']

    ENA_PROJECT_URL = 'http://www.ebi.ac.uk/ena/data/warehouse/filereport?accession={0}&result=analysis&fields=analysis_accession,study_accession,secondary_study_accession,sample_accession,secondary_sample_accession,analysis_title,analysis_type,center_name,first_public,last_updated,study_title,tax_id,scientific_name,analysis_alias,study_alias,submitted_md5,submitted_ftp,sample_alias,broker_name,sample_title&download=txt'

    def __init__(self, argv=None):
        self.assemblies = None
        super().__init__(argv)

    @staticmethod
    def add_arguments(parser):
        assembly_group = parser.add_mutually_exclusive_group()
        assembly_group.add_argument("-as", "--assemblies", nargs='+',
                                    help="Assembly ERZ accession(s), whitespace separated. "
                                         "Use to download only certain project assemblies")
        assembly_group.add_argument("--assembly-list", help='File containing line-separated assembly accessions')
        return parser

    def _validate_args(self):
        if not any([self.args.assemblies, self.args.assembly_list, self.args.projects, self.args.project_list]):
            raise ValueError('No data specified, please use -as, --assembly-list, -p or --project-list')
        elif not self.args.private and not (self.args.projects or self.args.project_list):
            raise NotImplementedError('Fetching studies from runs via FTP is not supported due '
                                      'to performance issues, please use --private mode')

    def _process_additional_args(self):
        if self.args.assembly_list:
            self.assemblies = self._read_line_sep_file(self.args.assemblies_list)
        else:
            self.assemblies = self.args.assemblies

        if not self.args.projects and not self.args.project_list:
            logging.info('Fetching projects from list of runs')
            self.args.projects = self._get_studies_from_assemblies(self.assemblies)

    @staticmethod
    def _get_study_run_accessions(project_data):
        return list(filter(bool, [entry['RUN_ID'] for entry in project_data]))

    def _filter_assemblies_from_args(self, assembly_data, assembly_accession_field):
        if self.assemblies:
            assembly_data = list(filter(lambda r: r[assembly_accession_field] in self.assemblies, assembly_data))
        return assembly_data

    @staticmethod
    def _filter_assembly_accessions(assembly_accessions, existing_assembly_accessions):
        return list(filter(lambda r: r not in existing_assembly_accessions, assembly_accessions))

    def _filter_assembly_from_existing_downloads(self, project_accession, insertable_assemblies, run_accession_field):
        try:
            existing_assemblies = list(self.read_project_description_file(project_accession)['assembly_id'])
            return list(filter(lambda r: r[run_accession_field] not in existing_assemblies, insertable_assemblies))
        except FileNotFoundError:
            return insertable_assemblies

    def map_project_info_db_row(self, assembly):
        is_submitted_file = assembly['DATA_FILE_ROLE'] == 'SUBMISSION_FILE'
        assembly['DATA_FILE_PATH'] = self._filter_secondary_files(assembly['DATA_FILE_PATH'])
        assembly['files'] = self._get_raw_filenames(assembly['DATA_FILE_PATH'], assembly['RUN_ID'], is_submitted_file)
        return {
            'study_id': assembly['STUDY_ID'],
            'sample_id': assembly['SAMPLE_ID'],
            'run_id': assembly['RUN_ID'],
            'file': ";".join(assembly['files']),
            'file_path': assembly['DATA_FILE_PATH'],
            'tax_id': assembly['TAX_ID'],
            'scientific_name': 'n/a',
        }

    def _retrieve_project_info_db(self, project_accession):
        # Get all generated study data from ENA
        study_run_data = ERADAO(self.eradao).retrieve_generated_files(project_accession)
        # Add all runs which don't have generated data
        if self._study_has_permitted_broker(project_accession):
            submitted_files = ERADAO(self.eradao).retrieve_submitted_files(project_accession)
            existing_accessions = self._get_study_run_accessions(study_run_data)
            study_run_data.extend([f for f in submitted_files if f['RUN_ID'] not in existing_accessions])
        else:
            broker = self.study_brokers[project_accession]
            logging.debug('Study {} does not have a trusted broker ({}), submitted_ftp files will be ignored'.format(
                project_accession, broker))

        # Allow force mode to bypass filtering
        if not self.force_mode:
            study_run_data = self._filter_assemblies_from_args(study_run_data, 'RUN_ID')
            study_run_data = self._filter_assembly_from_existing_downloads(project_accession, study_run_data, 'RUN_ID')
        return study_run_data

    @staticmethod
    def is_trusted_ftp_data(data, trusted_brokers):
        return data['submitted_ftp'] != '' and data['broker_name'] in trusted_brokers

    def _filter_ftp_broker_names(self, data):
        trusted_brokers = self.config['trustedBrokers']
        return [d for d in data if self.is_trusted_ftp_data(d, trusted_brokers)]

    def map_datafields_ftp_2_db(self, rundata):
        is_submitted_file = rundata['submitted_ftp'] is not ''
        rundata['STUDY_ID'] = rundata.pop('secondary_study_accession')
        rundata['SAMPLE_ID'] = rundata.pop('secondary_sample_accession')
        rundata['DATA_FILE_ROLE'] = 'SUBMISSION_FILE' if is_submitted_file else 'GENERATED_FILE'
        rundata['DATA_FILE_PATH'] = rundata.pop('submitted_ftp')
        rundata['MD5'] = rundata.pop('submitted_md5')
        rundata['RUN_ID'] = rundata.pop('analysis_accession')
        rundata['TAX_ID'] = rundata.pop('tax_id')
        rundata['files'] = self._get_raw_filenames(rundata['DATA_FILE_PATH'], rundata['RUN_ID'], is_submitted_file)
        return rundata

    def _retrieve_project_info_ftp(self, project_accession):
        data = self._retrieve_ena_url(self.ENA_PROJECT_URL.format(project_accession))
        trusted_data = self._filter_ftp_broker_names(data)
        insertable_runs = self._filter_assemblies_from_args(trusted_data, 'analysis_accession')
        insertable_runs = self._filter_assembly_from_existing_downloads(project_accession, insertable_runs,
                                                                        'analysis_accession')
        return list(map(self.map_datafields_ftp_2_db, insertable_runs))

    def write_project_files(self, project_accession, new_assemblies):
        new_assembly_rows = list(map(self.map_project_info_db_row, new_assemblies))
        self.write_project_description_file(project_accession, new_assembly_rows)

        new_download_rows = []
        for run in new_assembly_rows:
            for file in run['file'].split(';'):
                row = '\t'.join([run['file_path'], os.path.join('raw', os.path.basename(file))]) + '\n'
                new_download_rows.append(row)
                # self.write_project_download_file(project_accession, new_download_rows)


def main():
    data_fetcher = FetchAssemblies()
    data_fetcher.fetch()


if __name__ == '__main__':
    main()
