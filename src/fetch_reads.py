import re
import sys
import logging
from operator import itemgetter
from src.ERADAO import ERADAO

from src.abstract_fetch import AbstractDataFetcher

path_re = re.compile(r'(.*)/(.*)')


class FetchReads(AbstractDataFetcher):
    ENA_PORTAL_API_URL = 'https://www.ebi.ac.uk/ena/portal/api/search?dataPortal=metagenome&dccDataOnly=false&result=' \
                          'read_run&format=json&query=secondary_study_accession=%22{0}%22&fields=study_accession,' \
                          'secondary_study_accession,sample_accession,secondary_sample_accession,experiment_accession,' \
                          'run_accession,instrument_model,library_layout,fastq_ftp,fastq_md5,submitted_ftp,submitted_md5,' \
                          'library_strategy,broker_name,library_source&download=true'
    ENA_PORTAL_API_BY_RUN = 'https://www.ebi.ac.uk/ena/portal/api/search?dataPortal=metagenome&dccDataOnly=false&result' \
                            '=read_run&format=json&query=run_accession=%22{0}%22&fields=secondary_study_accession&' \
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
        #elif not self.args.private and not (self.args.projects or self.args.project_list):
        #    raise NotImplementedError('Fetching studies from runs via FTP is not supported due '
        #                              'to performance issues, please use --private mode')

    def _process_additional_args(self):
        if self.args.run_list:
            self.runs = self._read_line_sep_file(self.args.run_list)
        else:
            self.runs = self.args.runs

        if not self.args.projects and not self.args.project_list:
            logging.info('Fetching projects from list of assemblies')
            self.args.projects = self._get_project_accessions_from_runs(self.runs)

    def _retrieve_project_info_from_api(self, project_accession):
        data = self._retrieve_ena_url(self.ENA_PORTAL_API_URL.format(project_accession))
        logging.info("Retrieved {count} runs for study {project_accession} from "
                     "the ENA Portal API.".format(count=len(data), project_accession=project_accession))
        for d in data:
            if d['fastq_ftp'] == '':
                data.remove(d)
                logging.info("The ftp location for run {} is not available yet".format(d['run_accession']))
        return list(map(self.map_datafields_ftp_2_data, data))

    def map_datafields_ftp_2_data(self, rundata):
        is_submitted_file = rundata['submitted_ftp'] is not ''
        rundata['STUDY_ID'] = rundata.pop('secondary_study_accession')
        rundata['SAMPLE_ID'] = rundata.pop('secondary_sample_accession')
        rundata['RUN_ID'] = rundata.pop('run_accession')
        rundata['DATA_FILE_ROLE'] = 'SUBMISSION_FILE' if is_submitted_file else 'GENERATED_FILE'
        file_paths = rundata.get('fastq_ftp') #or rundata.get('submitted_ftp')
        is_valid_filetype = [self._is_rawdata_filetype(f) for f in file_paths.split(';')]
        if not False in is_valid_filetype:
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

    def _get_project_accessions_from_runs(self, runs):
        project_list = set()
        for run in runs:
            data = self._retrieve_ena_url(self.ENA_PORTAL_API_BY_RUN.format(run))
            [project_list.add(d['secondary_study_accession']) for d in data]
        return project_list


def main():
    data_fetcher = FetchReads()
    data_fetcher.fetch()


if __name__ == '__main__':
    main()
