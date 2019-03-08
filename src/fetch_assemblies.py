import re
import logging
import os
import requests
import gzip

from src.ENADAO import ENADAO
from src.ERADAO import ERADAO

from src.abstract_fetch import AbstractDataFetcher

path_re = re.compile(r'(.*)/(.*)')


class FetchAssemblies(AbstractDataFetcher):
    DEFAULT_HEADERS = ['study_id', 'sample_id', 'analysis_id', 'file', 'file_path']

    ENA_PROJECT_URL = 'http://www.ebi.ac.uk/ena/data/warehouse/filereport?accession={0}&result=analysis&fields=analysis_accession,study_accession,secondary_study_accession,sample_accession,secondary_sample_accession,analysis_title,analysis_type,center_name,first_public,last_updated,study_title,analysis_alias,study_alias,submitted_md5,submitted_ftp,sample_alias,broker_name,sample_title&download=txt'

    def __init__(self, argv=None):
        self.ACCESSION_FIELD = 'ANALYSIS_ID'
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
            raise NotImplementedError('Fetching studies from assemblies via FTP is not supported due '
                                      'to performance issues, please use --private mode')

    def _process_additional_args(self):
        if self.args.assembly_list:
            self.assemblies = self._read_line_sep_file(self.args.assembly_list)
        else:
            self.assemblies = self.args.assemblies

        if not self.args.projects and not self.args.project_list:
            logging.info('Fetching projects from list of assemblies')
            self.args.projects = self._get_project_accessions_from_assemblies(self.assemblies)

    def _get_project_accessions_from_assemblies(self, assemblies):
        self.init_era_dao()
        # Get all generated study data from ENA
        study_records = ERADAO(self.eradao).retrieve_study_accessions_from_analyses(assemblies)
        return [s['STUDY_ID'] for s in study_records]

    @staticmethod
    def _get_study_assembly_accessions(project_data):
        return list(filter(bool, [entry['ANALYSIS_ID'] for entry in project_data]))

    def _filter_accessions_from_args(self, assembly_data, assembly_accession_field):
        if self.assemblies:
            assembly_data = list(filter(lambda r: r[assembly_accession_field] in self.assemblies, assembly_data))
        return assembly_data

    @staticmethod
    def _filter_assembly_accessions(assembly_accessions, existing_assembly_accessions):
        return list(filter(lambda r: r not in existing_assembly_accessions, assembly_accessions))

    def _filter_accessions_from_existing_downloads(self, project_accession, insertable_assemblies,
                                                   assembly_accession_field):
        try:
            existing_runs = self.read_project_description_file(project_accession).to_dict('records')
        except FileNotFoundError:
            return insertable_assemblies

        raw_dir = self.get_project_rawdir(project_accession)
        existing_runs = list(filter(lambda r: self.check_files_downloaded(raw_dir, r['analysis_id']+'.fasta.gz'), existing_runs))
        existing_run_ids = [r['analysis_id'] for r in existing_runs]
        return list(filter(lambda r: r[assembly_accession_field] not in existing_run_ids, insertable_assemblies))

    def map_project_info_db_row(self, assembly):
        return {
            'study_id': assembly['STUDY_ID'],
            'sample_id': assembly['SAMPLE_ID'],
            'analysis_id': assembly['ANALYSIS_ID'],
            'file': assembly['file'],
            'file_path': assembly['DATA_FILE_PATH'],
            'scientific_name': 'n/a',
            'md5': assembly['MD5']
        }

    def _get_assembly_metadata(self, secondary_project_accession):
        return ERADAO(self.eradao).retrieve_assembly_metadata(secondary_project_accession)

    def _get_study_wgs_analyses(self, primary_project_accession):
        return ENADAO(self.enadao).retrieve_assembly_data(primary_project_accession)

    def _retrieve_project_info_db(self, project_accession):
        # Get all generated study data from ENA
        metadata_analyses = self._get_assembly_metadata(project_accession)
        project_acc = metadata_analyses[0]['PROJECT_ID']

        wgs_analyses = self._get_study_wgs_analyses(project_acc)

        study_analyses = self._combine_analyses(metadata_analyses, wgs_analyses)
        # Allow force mode to bypass filtering
        for data in study_analyses:
            data['DATA_FILE_PATH'], data['file'], data['MD5'] = self._get_raw_filenames(
                data['DATA_FILE_PATH'],
                data['MD5'],
                data['ANALYSIS_ID'],
                True)
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

    def map_datafields_ftp_2_db(self, assemblydata):
        is_submitted_file = assemblydata['submitted_ftp'] is not ''
        assemblydata['STUDY_ID'] = assemblydata.pop('secondary_study_accession')
        assemblydata['SAMPLE_ID'] = assemblydata.pop('secondary_sample_accession')
        assemblydata['DATA_FILE_PATH'] = assemblydata.pop('submitted_ftp')
        assemblydata['ANALYSIS_ID'] = assemblydata.pop('analysis_accession')
        assemblydata['DATA_FILE_PATH'], assemblydata['file'], assemblydata['MD5'] = self._get_raw_filenames(
            assemblydata['DATA_FILE_PATH'],
            assemblydata.pop('submitted_md5'),
            assemblydata['ANALYSIS_ID'],
            is_submitted_file)
        return assemblydata

    def _retrieve_project_info_ftp(self, project_accession):
        data = self._retrieve_ena_url(self.ENA_PROJECT_URL.format(project_accession))
        insertable_assemblies = self._filter_ftp_broker_names(data)
        return list(map(self.map_datafields_ftp_2_db, insertable_assemblies))

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
        res = requests.post(
            self.config['enaAPIUrl'] + 'search/',
            headers=headers,
            data=payload,
            auth=(self.config['enaAPIUsername'], self.config['enaAPIPassword'])
        )
        if res.status_code == 401:
            raise EnvironmentError('Error 401: Authentication credentials missing for ENA API')
        
        results = res.json()
        if len(results) == 0:
            logging.error("No results received for analysis object: {0}".format(
                analysis_accession))
            logging.error("Shutting down the program now!")
        elif len(results) > 1:
            raise ValueError("Unexpected number of results received for analysis object: {0}".format(
                    analysis_accession))
        else:
            scientific_name = results[0]['scientific_name']
            return scientific_name

    def download_raw_files(self, project_accession, new_assemblies):
        raw_dir = self.get_project_rawdir(project_accession)
        os.makedirs(raw_dir, exist_ok=True)
        for assembly in new_assemblies:
            download_sources = assembly['DATA_FILE_PATH']
            filenames = assembly['file']
            file_md5s = assembly['MD5']
            for dl_file, dl_name, dl_md5 in zip(download_sources, filenames, file_md5s):
                dest = os.path.join(raw_dir, dl_name)
                dir, basename = os.path.split(dest)
                name, ext = os.path.splitext(basename)
                unmapped_dest = os.path.join(dir, name + '_unmapped' + ext)
                was_downloaded = self.download_raw_file(dl_file, unmapped_dest, dl_md5)

                mapped_md5 = self.get_md5_file(dl_file)
                if (not self._is_file_valid(dest, mapped_md5)) or was_downloaded or self.force_mode:
                    logging.debug('Remapping ' + dl_file)
                    self.rename_fasta_headers(unmapped_dest, dest, project_accession, assembly['ANALYSIS_ID'])

    def get_contig_range_from_api(self, project_accession, analysis_id):
        '''
        Function to get the contig range from the ENA swagger
        input: ERZ_id: individual analysis_id
        output: string corresponding to the names of the first and last contigs separated by hyphen (ex: OFEO01000001-OFEO01173648)
        '''
        logging.debug('Getting contig names')

        webin = self._retrieve_webin_account_id(project_accession)

        url = 'https://www.ebi.ac.uk/ena/submit/report/analysis-process/{}/'.format(analysis_id)
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
        logging.debug('Got contig names')

        return data_contig

    @staticmethod
    def parse_wgs_seq_acc_range(line):
        """
            Parses a line in format scaffolds:FWWM01000001-FWWM01391746
            in format <assembly-type>:<wgs seq set acc><number>-<wgs seq set acc><number>
        :param line:
        :return:
        """
        wgs_seq_set_acc, start, end = re.findall(r'.+:(.{6})(\d+)-\1(\d+)', line)[0]
        return wgs_seq_set_acc, int(start), int(end)

    def rename_fasta_headers(self, unmapped_fasta_file, fasta_file, project_accession, analysis_id):
        scientific_name = self.get_scientific_name(analysis_id)
        contig_names = self.get_contig_range_from_api(project_accession, analysis_id)
        wgs_seq_set_acc, first_contig_number, last_contig_number = self.parse_wgs_seq_acc_range(contig_names)
        counter = 0
        logging.debug('Iterating through ' + fasta_file)
        with gzip.open(unmapped_fasta_file, 'rt', encoding='utf-8') as original_f:
            with gzip.open(fasta_file, 'wt', encoding='utf-8') as new_f:
                for line in original_f:
                    if line.startswith('>'):
                        contig_name = line[1:]
                        contig_number = (first_contig_number + counter)
                        contig_number = '{:08}'.format(contig_number)
                        line = '>ENA|{wgs_seq_set_acc}{contig_number}|{wgs_seq_set_acc}{contig_number}.1 ' \
                               '{scientific_name} genome assembly, ' \
                               'contig: {contig_name}'.format(wgs_seq_set_acc=wgs_seq_set_acc,
                                                              contig_number=contig_number,
                                                              scientific_name=scientific_name,
                                                              contig_name=contig_name)
                        counter += 1
                    new_f.write(line)
        logging.debug('Finished re-mapping fasta file')
        self.write_md5(fasta_file)


def main():
    data_fetcher = FetchAssemblies()
    data_fetcher.fetch()


if __name__ == '__main__':
    main()
