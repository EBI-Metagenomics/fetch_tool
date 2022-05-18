import re
import logging
import os
import gzip
import json

from src.ENADAO import ENADAO

from src.abstract_fetch import AbstractDataFetcher
from src.exceptions import NoDataError

path_re = re.compile(r'(.*)/(.*)')


class Analysis(object):
    def __init__(self, analysis_accession, assembly_type, status_id):
        self.analysis_accession = analysis_accession
        self.assembly_type = assembly_type
        self.status_id = status_id


class FetchAssemblies(AbstractDataFetcher):
    ENA_PORTAL_BASE_API_URL = 'https://www.ebi.ac.uk/ena/portal/api/search?'

    ENA_PORTAL_FIELDS = [
        'analysis_accession',
        'study_accession',
        'secondary_study_accession',
        'sample_accession',
        'secondary_sample_accession',
        'analysis_title',
        'analysis_type',
        'center_name',
        'first_public',
        'last_updated',
        'study_title',
        'analysis_alias',
        'study_alias',
        'submitted_md5',
        'submitted_ftp',
        'generated_md5',
        'generated_ftp',
        'sample_alias',
        'broker_name',
        'sample_title',
        'assembly_type',
    ]

    ENA_PORTAL_RUN_FIELDS = 'secondary_study_accession'

    ENA_PORTAL_PARAMS = [
        'dataPortal=metagenome',
        'dccDataOnly=false',
        'result=analysis',
        'format=json',
        'download=true',
        'fields='
    ]

    # query
    ENA_PORTAL_QUERY = "query=secondary_study_accession=%22{0}%22%20AND%20assembly_type=%22{1}%22"
    ENA_PORTAL_RUN_QUERY = "query=analysis_accession=%22{0}%22%20AND%20assembly_type=%22{1}%22"

    ENA_PORTAL_API_URL = ENA_PORTAL_BASE_API_URL + "&".join(ENA_PORTAL_PARAMS) + ','.join(ENA_PORTAL_FIELDS) + "&" + ENA_PORTAL_QUERY
    ENA_PORTAL_API_BY_RUN = ENA_PORTAL_BASE_API_URL + "&".join(ENA_PORTAL_PARAMS) + ENA_PORTAL_RUN_FIELDS + "&" + ENA_PORTAL_RUN_QUERY

    #not in use
    ENA_WGS_SET_API_URL = 'https://www.ebi.ac.uk/ena/portal/api/search?result=wgs_set&query=study_accession=' \
                          '%22{0}%22&fields=study_accession,assembly_type,scientific_name,fasta_file&format=tsv'

    def __init__(self, argv=None):
        self.ACCESSION_FIELD = 'ANALYSIS_ID'
        self.assemblies = None
        super().__init__(argv)
        self.init_ena_dao()

    @staticmethod
    def add_arguments(parser):
        assembly_group = parser.add_mutually_exclusive_group()
        assembly_group.add_argument("-as", "--assemblies", nargs='+',
                                    help="Assembly ERZ accession(s), whitespace separated. "
                                         "Use to download only certain project assemblies")
        # TODO: Also the older INSDC style metagenome assemblies produced by MGnify are not support that way at moment
        parser.add_argument('--assembly-type', help="Assembly type",
                                    choices=["primary metagenome", "Metagenome-Assembled Genome (MAG)",
                                             "binned metagenome", "metatranscriptome"],
                                    default="primary metagenome")
        assembly_group.add_argument("--assembly-list", help="File containing line-separated assembly accessions")

        return parser

    def _validate_args(self):
        if not any([self.args.assemblies, self.args.assembly_list, self.args.projects, self.args.project_list]):
            raise ValueError('No data specified, please use -as, --assembly-list, -p or --project-list')

    def _process_additional_args(self):
        self.assembly_type = self.args.assembly_type

        if self.args.assembly_list:
            self.assemblies = self._read_line_sep_file(self.args.assembly_list)
        else:
            self.assemblies = self.args.assemblies

        if not self.args.projects and not self.args.project_list:
            logging.info('Fetching projects from list of assemblies')
            self.args.projects = self._get_project_accessions_from_assemblies(self.assemblies)

    def _retrieve_project_info_from_api(self, project_accession):
        #self._retrieve_insdc_assemblies(primary_project_accession)
        raise_error = False if len(self.projects) > 1 else True     #allows script to continue to next project if one fails
        data = self._retrieve_ena_url(self.ENA_PORTAL_API_URL.format(project_accession, self.assembly_type), raise_on_204=raise_error)
        if not data:
            return
        logging.info("Retrieved {count} assemblies for study {project_accession} from "
                     "the ENA Portal API.".format(count=len(data), project_accession=project_accession))
        mapped_data = []
        for d in data:
            if not d['generated_ftp']:
                logging.info(
                    "The generated ftp location for assembly {} is not available yet".format(d['analysis_accession']))
            if d['analysis_type'] == 'SEQUENCE_ASSEMBLY' and d['generated_ftp']:
                if self._is_rawdata_filetype(os.path.basename(d['generated_ftp'])):  #filter filenames not fasta
                    raw_data_file_path, file_, md5_ = self._get_raw_filenames(
                        d.get('generated_ftp'),
                        d.get('generated_md5'),
                        d.get('analysis_accession'),
                        bool(d.get('submitted_ftp'))
                    )
                    mapped_data.append({
                        'STUDY_ID': d.get('secondary_study_accession'),
                        'SAMPLE_ID': d.get('secondary_sample_accession'),
                        'ANALYSIS_ID': d.get('analysis_accession'),
                        'DATA_FILE_PATH': raw_data_file_path,
                        'file': file_,
                        'MD5': md5_
                    })
        return mapped_data

    def _filter_accessions_from_args(self, assembly_data, assembly_accession_field):
        if self.assemblies:
            data = list(filter(lambda r: (r[assembly_accession_field] in self.assemblies), assembly_data))
            return data
        else:
            return assembly_data

    def map_project_info_to_row(self, assembly):
        return {
            'study_id': assembly['STUDY_ID'],
            'sample_id': assembly['SAMPLE_ID'],
            'analysis_id': assembly['ANALYSIS_ID'],
            'file': assembly['file'],
            'file_path': assembly['DATA_FILE_PATH'],
            'scientific_name': 'n/a',
            'md5': assembly['MD5']
        }

    def _get_project_accessions_from_assemblies(self, assemblies):
        project_list = set()
        for assembly in assemblies:
            data = self._retrieve_ena_url(self.ENA_PORTAL_API_BY_RUN.format(assembly, self.assembly_type),
                                          raise_on_204=False)
            if data:
                [project_list.add(d['secondary_study_accession']) for d in data]
        if not len(project_list):
            raise NoDataError(self.NO_DATA_MSG)
        return project_list

    """
    The following functions are for old-style assemblies and are not in use.
    Remove after resubmission of old assemblies and replace with simple fetch of INSDC submitted files
    (with no renaming or file reformatting).
    """

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
        # not in use yet. needs primary accession.
        mapped_data = []
        json_data = self._retrieve_ena_url(self.ENA_WGS_SET_API_URL.format(study_accession))
        logging.info("Retrieved {count} assemblies of type wgs_set from ENA Portal API.".format(count=len(json_data)))
        if len(json_data):
            study_assembly_data = self._get_study_wgs_analyses(study_accession)
            study_analyses = self._combine_analyses(json_data, study_assembly_data)
            self._download_insdc_assemblies(study_analyses, study_accession)
            for data in study_analyses:
                mapped_data.append({
                    'study_id': study_accession,
                    'sample_id': data['SAMPLE_ID'],
                    'analysis_id': data['accession'],
                    'file': data['ASSEMBLY_ID'] + '.fasta.gz',
                    'file_path': data['DATA_FILE_PATH']
                })
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

    def add_insdc_to_file(self, data, study_accession):
        headers = ['study_id', 'sample_id', 'analysis_id', 'file', 'file_path']
        with open(self.get_project_insdc_txt_file(study_accession), 'w') as txt_file:
            txt_file.write('\t'.join(headers) +'/n')
            for d in data:
                json.dump(d, txt_file)
                txt_file.write('/n')


def main():
    data_fetcher = FetchAssemblies()
    data_fetcher.fetch()


if __name__ == '__main__':
    main()
