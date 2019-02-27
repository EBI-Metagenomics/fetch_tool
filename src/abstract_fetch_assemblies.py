import re
import logging
import shutil
import sys
import os
import requests
import gzip

from src.ENADAO import ENADAO
from src.ERADAO import ERADAO

from src.abstract_fetch import AbstractDataFetcher

path_re = re.compile(r'(.*)/(.*)')


class FetchAssemblies(AbstractDataFetcher):
    DEFAULT_HEADERS = ['study_id', 'sample_id', 'analysis_id', 'files', 'file_path']

    ENA_PROJECT_URL = 'http://www.ebi.ac.uk/ena/data/warehouse/filereport?accession={0}&result=analysis&fields=analysis_accession,study_accession,secondary_study_accession,sample_accession,secondary_sample_accession,analysis_title,analysis_type,center_name,first_public,last_updated,study_title,tax_id,scientific_name,analysis_alias,study_alias,submitted_md5,submitted_ftp,sample_alias,broker_name,sample_title&download=txt'

    def __init__(self, argv=None):
        self.assemblies = None
        super().__init__(argv)
        self.init_era_dao()
        self.init_ena_dao()

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
            self.args.projects = self._get_project_accessions_from_assemblies(self.assemblies)

    def _get_project_accessions_from_assemblies(self, assemblies):
        self.init_era_dao()
        # Get all generated study data from ENA
        study_records = ERADAO(self.eradao).retrieve_study_accessions(assemblies)
        return [s['STUDY_ID'] for s in study_records]

    @staticmethod
    def _get_study_run_accessions(project_data):
        return list(filter(bool, [entry['ANALYSIS_ID'] for entry in project_data]))

    def _filter_assemblies_from_args(self, assembly_data, assembly_accession_field):
        if self.assemblies:
            assembly_data = list(filter(lambda r: r[assembly_accession_field] in self.assemblies, assembly_data))
        return assembly_data

    @staticmethod
    def _filter_assembly_accessions(assembly_accessions, existing_assembly_accessions):
        return list(filter(lambda r: r not in existing_assembly_accessions, assembly_accessions))

    def _filter_assembly_from_existing_downloads(self, project_accession, insertable_assemblies, run_accession_field):
        try:
            existing_assemblies = list(self.read_project_description_file(project_accession)['analysis_id'])
            return list(filter(lambda r: r[run_accession_field] not in existing_assemblies, insertable_assemblies))
        except FileNotFoundError:
            return insertable_assemblies

    def map_project_info_db_row(self, assembly):
        return {
            'study_id': assembly['STUDY_ID'],
            'sample_id': assembly['SAMPLE_ID'],
            'analysis_id': assembly['ANALYSIS_ID'],
            'files': assembly['files'],
            'file_path': assembly['DATA_FILE_PATH'],
            'tax_id': assembly['TAX_ID'],
            'scientific_name': 'n/a',
            'md5': assembly['MD5']
        }

    def _retrieve_project_info_db(self, project_accession):
        # Get all generated study data from ENA
        metadata_analyses = ERADAO(self.eradao).retrieve_assembly_metadata(project_accession)
        project_acc = metadata_analyses[0]['PROJECT_ID']

        wgs_analyses = ENADAO(self.enadao).retrieve_assembly_data(project_acc)

        study_analyses = self._combine_analyses(metadata_analyses, wgs_analyses)

        # Allow force mode to bypass filtering
        return study_analyses

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

    @staticmethod
    def is_trusted_ftp_data(data, trusted_brokers):
        is_trusted = data['submitted_ftp'] != '' and \
                     (data.get('broker_name') in trusted_brokers
                      or data.get('center_name') in trusted_brokers)
        if not is_trusted:
            logging.warning('Study contains untrusted broker {} and centername {}'.format(data.get('broker_name'),
                                                                                          data.get('center_name')))
        return is_trusted

    def _filter_ftp_broker_names(self, data):
        trusted_brokers = self.config['trustedBrokers']
        return [d for d in data if self.is_trusted_ftp_data(d, trusted_brokers)]

    def map_datafields_ftp_2_db(self, rundata):
        is_submitted_file = rundata['submitted_ftp'] is not ''
        rundata['STUDY_ID'] = rundata.pop('secondary_study_accession')
        rundata['SAMPLE_ID'] = rundata.pop('secondary_sample_accession')
        rundata['DATA_FILE_PATH'] = rundata.pop('submitted_ftp')
        rundata['ANALYSIS_ID'] = rundata.pop('analysis_accession')
        rundata['TAX_ID'] = rundata.pop('tax_id')
        rundata['DATA_FILE_PATH'], rundata['files'], rundata['MD5'] = self._get_raw_filenames(rundata['DATA_FILE_PATH'],
                                                                                              rundata.pop(
                                                                                                  'submitted_md5'),
                                                                                              rundata['ANALYSIS_ID'],
                                                                                              is_submitted_file)
        return rundata

    def filter_by_accessions(self, project_accession, data):
        if not self.force_mode:
            data = self._filter_assemblies_from_args(data, 'ANALYSIS_ID')
            data = self._filter_assembly_from_existing_downloads(project_accession, data, 'ANALYSIS_ID')
        return data

    def _retrieve_project_info_ftp(self, project_accession):
        data = self._retrieve_ena_url(self.ENA_PROJECT_URL.format(project_accession))
        insertable_runs = self._filter_ftp_broker_names(data)
        return list(map(self.map_datafields_ftp_2_db, insertable_runs))

    def write_project_files(self, project_accession, new_assemblies):
        new_assembly_rows = list(map(self.map_project_info_db_row, new_assemblies))
        self.write_project_description_file(project_accession, new_assembly_rows, 'analysis_id')

        self.write_project_download_file(project_accession, new_assembly_rows)

    project_webin_accounts = {}

    def _retrieve_webin_account_id(self, project_accession):
        if project_accession not in self.project_webin_accounts:
            results = ERADAO(self.eradao).retrieve_webin_accound_id(project_accession)
            self.project_webin_accounts[project_accession] = results[0]['SUBMISSION_ACCOUNT_ID']
        return self.project_webin_accounts[project_accession]

    def get_scientific_name(self, analysis_accession):
        '''
        Function to get analysis_accession (ERZ_id), analysis_alias (file_name) and scientific_name for all assemblies of a project
        input: study_id
        input: dictionary generated by the extract_data function (dic[ERZ_id]:[contig_name, GCA id])
        input API_URL: url for the ENA API
        output: api_dic: dictionary with ERZ_id as keys and values are list of file_name and scientific_name
        '''
        payload = {
            "result": "analysis",
            "query": 'analysis_accession="{0}"'.format(analysis_accession),
            "fields": "analysis_accession,scientific_name",
            "dataPortal": "metagenome",
            "format": "json",
        }
        headers = {'content-type': 'application/x-www-form-urlencoded'}
        try:
            res = requests.post(
                self.config['enaAPIUrl'],
                headers=headers,
                data=payload,
                auth=(self.config['enaAPIUsername'], self.config['enaAPIPassword'])
            )
        except:
            raise
        results = res.json()
        if len(results) == 0:
            logging.error("No results received for analysis object: {0}".format(
                analysis_accession))
            logging.error("Shutting down the program now!")
        elif len(results) > 1:
            logging.error(
                "Unexpected number of results received for analysis object: {0}".format(
                    analysis_accession))
        else:
            scientific_name = results[0]['scientific_name']
            return scientific_name

    def download_raw_files(self, project_accession, new_runs):
        raw_dir = self.get_project_rawdir(project_accession)
        os.makedirs(raw_dir, exist_ok=True)
        for run in new_runs:
            download_sources = run['DATA_FILE_PATH']
            filenames = run['files']
            file_md5s = run['MD5']
            for dl_file, dl_name, dl_md5 in zip(download_sources, filenames, file_md5s):
                dest = os.path.join(raw_dir, dl_name)
                self.download_raw_file(dl_file, dest, dl_md5)
                self.rename_fasta_headers(dest, project_accession, run['ANALYSIS_ID'])

    def get_contig_range_from_api(self, project_accession, analysis_id):
        '''
        Function to get the contig range from the ENA swagger
        input: ERZ_id: individual analysis_id
        output: string corresponding to the names of the first and last contigs separated by hyphen (ex: OFEO01000001-OFEO01173648)
        '''
        webin = self._retrieve_webin_account_id(project_accession)

        url = f'https://www.ebi.ac.uk/ena/submit/report/analysis-process/{analysis_id}/'
        payload = {"format": "json"}
        headers = {"Accept": "*/*"}
        creds = (webin, self.config['enaMasterPassword'])
        res = requests.get(
            url,
            headers=headers,
            data=payload,
            auth=creds
        )
        data = res.json()[0]
        data_contig = data['report']['acc'].split(",")[0]
        return data_contig

    @staticmethod
    def parse_wgs_seq_acc_range(line):
        """
            This method parses the following line:
            scaffolds:FWWM01000001-FWWM01391746
            <assembly-type>:<wgs seq set acc><number>-<wgs seq set acc><number>

            This method parses the WGS sequence set accession, FWWM01 and
            the ordinal number range, 000001-391746

            The line we have to parse is retrieved from an API call.
            The API is maintained by the ENA. Further down below you will find
            an example call and an example response:

            Example call:
            curl -X GET "https://www.ebi.ac.uk/ena/submit/report/analysis-process/ERZ404940?format=json&max-results=100" -H  "accept: */*" -H  "authorization: Basic V2ViaW4tNDI5NzA6M25AIUAyLTEyOA=="

            Example response body (05/10/2018):
            [
              {
                "report": {
                  "id": "ERZ404940",
                  "analysisType": "SEQUENCE_ASSEMBLY",
                  "acc": "scaffolds:FWWM01000001-FWWM01391746",
                  "processingStatus": "COMPLETED",
                  "processingError": null
                },
                "links": []
              }
            ]
        :param amount:
        :return:
        """
        wgs_seq_set_acc, start, end = None, None, None
        if not line:
            logging.warning("Value for line is undefined!")
            return None, None, None
        warn_msg = "Unexpected format of line: {0}".format(line)
        # Split case like  contigs:ODOD01000001-ODOD01010452,contigs(set):ODOD01,genome:GCA_900206415.1
        first_split = line.split(',')[0]
        # Split scaffolds:FWPU01000001-FWPU01204773
        second_split = first_split.split(':')
        if len(second_split) == 2:
            # Split FWPU01000001-FWPU01204773
            third_split = second_split[1].split('-')
            if len(third_split) == 2:
                start = third_split[0][6:]
                end = third_split[1][6:]
                # wgs seq set acc
                wgs_seq_set_acc = third_split[0][:6]
            else:
                logging.warning(warn_msg)
        else:
            logging.warning(warn_msg)
        return wgs_seq_set_acc, int(start), int(end)

    def rename_fasta_headers(self, fasta_file, project_accession, analysis_id):
        scientific_name = self.get_scientific_name(analysis_id)
        contig_names = self.get_contig_range_from_api(project_accession, analysis_id)
        wgs_seq_set_acc, first_contig_number, last_contig_number = self.parse_wgs_seq_acc_range(contig_names)
        counter = 0
        bak_fasta_file = fasta_file + '.bak'
        shutil.move(fasta_file, bak_fasta_file)
        with gzip.open(bak_fasta_file, 'rt', encoding='utf-8') as original_f:
            with gzip.open(fasta_file, 'wt', encoding='utf-8') as new_f:
                for line in original_f:
                    if line.startswith('>'):
                        contig_name = line[1:]
                        contig_number = (first_contig_number + counter)
                        contig_number = f'{contig_number:08}'
                        line = f'>ENA|{wgs_seq_set_acc}{contig_number}|{wgs_seq_set_acc}{contig_number}.1 ' \
                               f'{scientific_name} genome assembly, contig: {contig_name}'
                        counter += 1
                    new_f.write(line)


def main():
    data_fetcher = FetchAssemblies()
    data_fetcher.fetch()


if __name__ == '__main__':
    main()
