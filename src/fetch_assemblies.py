import re
import logging
import os
import requests
import gzip
import json

import sys

from src.ENADAO import ENADAO
from src.ERADAO import ERADAO

from src.abstract_fetch import AbstractDataFetcher

path_re = re.compile(r'(.*)/(.*)')




class Analysis(object):
    def __init__(self, analysis_accession, assembly_type, status_id):
        self.analysis_accession = analysis_accession
        self.assembly_type = assembly_type
        self.status_id = status_id


class FetchAssemblies(AbstractDataFetcher):
    ENA_PORTAL_API_URL = 'https://www.ebi.ac.uk/ena/portal/api/search?dataPortal=metagenome&dccDataOnly=false&result=' \
                         'analysis&format=json&query=secondary_study_accession=%22{0}%22%20AND%20assembly_type=%22{1}%22&' \
                         'fields=analysis_accession,study_accession,secondary_study_accession,sample_accession,' \
                         'secondary_sample_accession,analysis_title,analysis_type,center_name,first_public,last_updated' \
                         ',study_title,analysis_alias,study_alias,submitted_md5,submitted_ftp,generated_md5,' \
                         'generated_ftp,sample_alias,broker_name,sample_title,assembly_type&download=true'
    ENA_PORTAL_API_BY_RUN = 'https://www.ebi.ac.uk/ena/portal/api/search?dataPortal=metagenome&dccDataOnly=false&result' \
                            '=analysis&format=json&query=analysis_accession=%22{0}%22%20AND%20assembly_type=%22{1}%22&' \
                            'fields=secondary_study_accession&download=true'
    ENA_WGS_SET_API_URL = 'https://www.ebi.ac.uk/ena/portal/api/search?result=wgs_set&query=study_accession=' \
                          '%22{0}%22&fields=study_accession,assembly_type,scientific_name,fasta_file&format=tsv'

    #http: // www.ebi.ac.uk / ena / data / warehouse / filereport?accession = {0} & result = read_run & fields =

    PROGRAM_EXIT_MSG = "Program will exit now!"

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
        # TODO: Valid what's the correct value for binned metagenomes!
        # TODO: Also the older INSDC style metagenome assemblies produced by MGnify are not support that way at moment
        # TODO: The returned assembly type for those is an empty string, this needs to be handle differently
        # TODO: I imagine you could retrieve all assembly types and filter these empty stringed assembly types
        parser.add_argument('--assembly-type', help="Assembly type",
                                    choices=["primary metagenome", "Metagenome-Assembled Genome (MAG)",
                                             "binned metagenome", "metatranscriptome"],
                                    default="primary metagenome")
        assembly_group.add_argument("--assembly-list", help="File containing line-separated assembly accessions")

        return parser

    def _validate_args(self):
        if not any([self.args.assemblies, self.args.assembly_list, self.args.projects, self.args.project_list]):
            raise ValueError('No data specified, please use -as, --assembly-list, -p or --project-list')
        #elif not self.args.private and not (self.args.projects or self.args.project_list):
        #    raise NotImplementedError('Fetching studies from assemblies via FTP is not supported due '
        #                              'to performance issues, please use --private mode')

    def _process_additional_args(self):
        self.assembly_type = self.args.assembly_type

        if self.args.assembly_list:
            self.assemblies = self._read_line_sep_file(self.args.assembly_list)
        else:
            self.assemblies = self.args.assemblies

        if not self.args.projects and not self.args.project_list:
            logging.info('Fetching projects from list of assemblies')
            self.args.projects = self._get_project_accessions_from_assemblies(self.assemblies)

            #logging.error('Please specify the secondary project ID')
            #logging.warning(self.PROGRAM_EXIT_MSG)
            #sys.exit(1)

    def _retrieve_project_info_from_api(self, project_accession):
        #self._retrieve_insdc_assemblies(project_accession)
        data = self._retrieve_ena_url(self.ENA_PORTAL_API_URL.format(project_accession, self.assembly_type))
        logging.info("Retrieved {count} assemblies for study {project_accession} from "
                     "the ENA Portal API.".format(count=len(data), project_accession=project_accession))
        for d in data:
            if d['analysis_type'] != 'SEQUENCE_ASSEMBLY':
                data.remove(d)
            if not d['generated_ftp']:
                data.remove(d)
                logging.info("The generated ftp location for run {} is not available yet".format(d['analysis_accession']))

        #data_filtered = [d for d in data if d['analysis_type'] == 'SEQUENCE_ASSEMBLY' and len(d['generated_ftp'])]
        return list(map(self.map_datafields_ftp_2_data, data))


    def map_datafields_ftp_2_data(self, assemblydata):
        is_submitted_file = assemblydata['submitted_ftp'] is not ''
        assemblydata['STUDY_ID'] = assemblydata.pop('secondary_study_accession')
        assemblydata['SAMPLE_ID'] = assemblydata.pop('secondary_sample_accession')
        assemblydata['DATA_FILE_PATH'] = assemblydata.pop('generated_ftp')
        assemblydata['ANALYSIS_ID'] = assemblydata.pop('analysis_accession')
        is_valid_filetype = self._is_rawdata_filetype(assemblydata['DATA_FILE_PATH']) #filter empties and non fasta
        if is_valid_filetype:
            assemblydata['DATA_FILE_PATH'], assemblydata['file'], assemblydata['MD5'] = self._get_raw_filenames(
                assemblydata['DATA_FILE_PATH'],
                assemblydata.pop('generated_md5'),
                assemblydata['ANALYSIS_ID'],
                is_submitted_file)
            return assemblydata

    def _filter_accessions_from_args(self, assembly_data, assembly_accession_field):
        if self.assemblies:
            data = list(filter(lambda r: (r[assembly_accession_field] in self.assemblies), assembly_data))
            return data
        else:
            return assembly_data


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

    def _get_study_wgs_analyses(self, primary_project_accession):
        return ENADAO(self.enadao).retrieve_assembly_data(primary_project_accession)

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
            analysis_id = analysis['accession']
            if analysis_id in wgs_dict:
                new_analysis = analysis.copy()
                new_analysis.update(wgs_dict[analysis_id])
                new_analysis['DATA_FILE_PATH'] = new_analysis.pop('fasta_file')
                combine_analyses.append(new_analysis)
        return combine_analyses

    def _retrieve_insdc_assemblies(self, study_accession):
        # Step 1: Retrieve analysis results including assembly type from ENA Portal API
        mapped_data = []
        json_data = self._retrieve_ena_url(self.ENA_WGS_SET_API_URL.format('PRJEB22493'))
        logging.info("Retrieved {count} assemblies of type wgs_set from ENA Portal API.".format(count=len(json_data)))
        if len(json_data):
            study_assembly_data = self._get_study_wgs_analyses(study_accession)
            study_analyses = self._combine_analyses(json_data, study_assembly_data)
            self._download_insdc_assemblies(study_analyses, study_accession)
            for data in study_analyses:
                mapped_out = self.map_insdc_data_for_files(data, study_accession)
                mapped_data.append(mapped_out)
            self.add_insdc_to_file(mapped_data, study_accession)
        else:
            logging.info("No INSDC style assemblies. Skipping...")

    def _download_insdc_assemblies(self, study_analyses, study_accession):
        raw_dir = self.get_project_rawdir(study_accession)
        for study in study_analyses:
            filename = study['ASSEMBLY_ID'] + '.fasta.gz'
            dest = os.path.join(raw_dir, filename)
            unmapped_dest = dest + '.unmapped'
            try:
                logging.info("Downloading INSDC style assembly {}".format(filename))
                self.download_ftp(dest, study['DATA_FILE_PATH'], auth=False)
                self.rename_fasta_headers(unmapped_dest, dest, study['ASSEMBLY_ID'])
            except Exception as e:
                if self.ignore_errors:
                    logging.warning(e)
                else:
                    raise e

    def rename_fasta_headers(self, unmapped_fasta_file, fasta_file, analysis_id):
        counter = 1
        logging.debug('Iterating through ' + fasta_file)
        with gzip.open(unmapped_fasta_file, 'rt', encoding='utf-8') as original_f:
            with gzip.open(fasta_file, 'wt', encoding='utf-8') as new_f:
                for line in original_f:
                    if line.startswith('>'):
                        contig_name = line[1:]
                        line = '>{analysis_id}.{counter} {contig_name}'.format(analysis_id=analysis_id,
                                                                               counter=counter,
                                                                               contig_name=contig_name)
                        counter += 1
                    new_f.write(line)
        logging.debug('Finished re-mapping fasta file')
        self.write_md5(fasta_file)

    @staticmethod
    def map_insdc_data_for_files(self, insdc_data, study_accession,):
        return {
            'study_id': study_accession,
            'sample_id': insdc_data['SAMPLE_ID'],
            'analysis_id': insdc_data['accession'],
            'file': insdc_data['ASSEMBLY_ID'] + '.fasta.gz',
            'file_path': insdc_data['DATA_FILE_PATH'],
        }

    def add_insdc_to_file(self, data, study_accession):
        headers = ['study_id', 'sample_id', 'analysis_id', 'file', 'file_path']
        with open(self.get_project_insdc_txt_file(study_accession), 'w') as txt_file:
            txt_file.write('\t'.join(headers) +'/n')
            for d in data:
                json.dump(d, txt_file)
                txt_file.write('/n')

    def _get_project_accessions_from_assemblies(self, assemblies):
        project_list = set()
        for assembly in assemblies:
            data = self._retrieve_ena_url(self.ENA_PORTAL_API_BY_RUN.format(assembly, self.assembly_type))
            [project_list.add(d['secondary_study_accession']) for d in data]
        return project_list


def main():
    data_fetcher = FetchAssemblies()
    data_fetcher.fetch()


if __name__ == '__main__':
    main()
