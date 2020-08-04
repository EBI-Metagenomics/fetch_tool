import os
import pytest
import shutil
from unittest.mock import patch

from copy import deepcopy
from src import fetch_assemblies as afa

FIXTURES_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, 'fixtures'))


class TestFetchAssemblies:
    def test_argparse_should_include_additional_args(self):
        fetch = afa.FetchAssemblies(argv=['-p', 'ERP104225'])

        args = fetch.args
        accepted_args = {'projects', 'project_list', 'dir', 'verbose', 'force',
                         'private', 'interactive', 'config_file', 'assemblies', 'assembly_list', 'fix_desc_file',
                         'ignore_errors'}
        assert set(vars(args)) == accepted_args

    def test_validate_args_should_raise_exception_as_no_data_specified(self):
        with pytest.raises(ValueError):
            afa.FetchAssemblies(argv=['-f'])

    def test_process_additional_args_should_set_assemblies_from_arglist(self):
        assemblies = ['ERZ477685']
        fetch = afa.FetchAssemblies(argv=['-p', 'ERP104225', '-as'] + assemblies)
        assert fetch.assemblies == assemblies

    def test_process_additional_args_args_should_set_assemblies_from_file(self, tmpdir):
        assemblies = ['ERZ477685']
        tmpdir = str(tmpdir)
        assembly_file = os.path.join(tmpdir, 'assemblies.txt')
        with open(assembly_file, 'w') as f:
            f.write('\n'.join(assemblies))
        fetch = afa.FetchAssemblies(argv=['-p', 'ERP104225', '--assembly-list', assembly_file])
        assert fetch.assemblies == assemblies

    def test_validate_args_should_raise_error_if_fetching_studies_without_project_in_public_mode(self, tmpdir):
        assemblies = ['ERZ477685']
        tmpdir = str(tmpdir)
        runfile = os.path.join(tmpdir, 'assemblies.txt')
        with open(runfile, 'w') as f:
            f.write('\n'.join(assemblies))

        with pytest.raises(NotImplementedError):
            afa.FetchAssemblies(argv=['--assembly-list', runfile])

    def test_get_project_accessions_from_assemblies_should_return_empty(self):
        assert [] == afa.FetchAssemblies._get_study_assembly_accessions([])

    def test_get_project_accessions_from_assemblies_should_raise_keyerror_on_invalid_project_data(self):
        with pytest.raises(KeyError):
            afa.FetchAssemblies._get_study_assembly_accessions([{'PROJECT_ID': 'ERP104225'}])

    def test_get_project_accessions_from_assemblies_should_return_analysis_id(self):
        assembly_id = 'ERZ477685'
        assert [assembly_id] == afa.FetchAssemblies._get_study_assembly_accessions([{'ANALYSIS_ID': assembly_id}])

    def test_filter_accessions_from_args_should_return_empty(self):
        fetch = afa.FetchAssemblies(argv=['-p', 'ERP104225'])
        fetch.assemblies = 'ERZ477685'
        assert [] == fetch._filter_accessions_from_args([], 'analysis_id')

    def test_filter_accessions_from_args_should_return_filtered_assemblies(self):
        fetch = afa.FetchAssemblies(argv=['-p', 'ERP104225'])
        assemblies = ['ERZ477685', 'ERZ477686']
        fetch.assemblies = assemblies
        run_data = [
            {'analysis_id': 'ERZ477685'},
            {'analysis_id': 'ERZ477686'},
            {'analysis_id': 'ERZ477687'},
        ]
        assert run_data[0:2] == fetch._filter_accessions_from_args(run_data, 'analysis_id')

    def test_filter_assembly_accessions_should_return_empty(self):
        fetch = afa.FetchAssemblies(argv=['-p', 'ERP104225'])
        assert [] == fetch._filter_assembly_accessions([], [])

    def test_filter_assembly_accessions_should_not_return_accession(self):
        new_acc = 'ERZ477685'
        fetch = afa.FetchAssemblies(argv=['-p', 'ERP104225'])
        assert [] == fetch._filter_assembly_accessions([new_acc], ['ERZ477684', new_acc])

    def test_filter_assembly_accessions_should_return_accession(self):
        new_acc = 'ERZ477685'
        fetch = afa.FetchAssemblies(argv=['-p', 'ERP104225'])
        assert [new_acc] == fetch._filter_assembly_accessions([new_acc], ['ERZ477684'])

    def test_map_project_info_db_row_should_copy_fields(self):
        raw_data = {
            'STUDY_ID': 'ERP104225',
            'SAMPLE_ID': 'ERS599830',
            'ANALYSIS_ID': 'ERZ477685',
            'file': 'ERZ477685.fasta.gz',
            'DATA_FILE_PATH': '/tmp/ERP104225/ERZ477685.fasta.gz',
            'MD5': ('md51',),
        }
        fetch = afa.FetchAssemblies(argv=['-p', 'ERP003634'])
        transform = fetch.map_project_info_db_row(raw_data)
        equivalent_fields = (
            ('STUDY_ID', 'study_id'),
            ('SAMPLE_ID', 'sample_id'),
            ('ANALYSIS_ID', 'analysis_id'),
        )
        for f1, f2 in equivalent_fields:
            assert raw_data[f1] == transform[f2]

    def test_map_project_info_db_row_should_generate_file_list(self):
        raw_data = {
            'STUDY_ID': 'ERP104225',
            'SAMPLE_ID': 'ERS599830',
            'ANALYSIS_ID': 'ERZ477685',
            'file': 'ERZ477685.fasta.gz',
            'DATA_FILE_PATH': '/tmp/ERP104225/ERZ477685.fasta.gz',
            'DATA_FILE_ROLE': 'SUBMITTED',
            'MD5': ('md51',),
        }
        fetch = afa.FetchAssemblies(argv=['-p', 'ERP003634'])
        transform = fetch.map_project_info_db_row(raw_data)
        assert transform['file'] == 'ERZ477685.fasta.gz'

#    def test_is_trusted_ftp_data_should_be_trusted_from_broker(self):
#        fetch = afa.FetchAssemblies(argv=['-p', 'ERP003634'])
#        trusted_brokers = ['EMG']
#        assert fetch.is_trusted_ftp_data({'fasta_ftp': '',
#                                          'submitted_ftp': 'filepath',
#                                          'broker_name': 'EMG'},
#                                         trusted_brokers)#

#    def test_is_trusted_ftp_data_should_be_trusted_as_generated_data(self):
#        fetch = afa.FetchAssemblies(argv=['-p', 'ERP003634'])
#        trusted_brokers = ['EMG']
#        assert fetch.is_trusted_ftp_data({'fasta_ftp': 'filepath', 'submitted_ftp': 'filepath', 'broker_name': 'EMG'},
#                                         trusted_brokers)

#    def test_filter_ftp_broker_names_should_return_empty(self):
#        fetch = afa.FetchAssemblies(argv=['-p', 'ERP003634'])
#        assert [] == fetch._filter_ftp_broker_names([])

#    def test_filter_ftp_broker_names_should_filter_brokers(self):
#        fetch = afa.FetchAssemblies(argv=['-p', 'ERP003634'])
#        fetch.config['trustedBrokers'] = ['EMG']
#        run_data = [
#            {'fasta_ftp': '', 'submitted_ftp': 'datafile', 'broker_name': 'EMG'},
#            {'fasta_ftp': '', 'submitted_ftp': 'datafile', 'broker_name': 'NOT_EMG'},
#        ]
#        assert 'EMG' == fetch._filter_ftp_broker_names(run_data)[0]['broker_name']

#    def test_filter_ftp_broker_names_should_allow_generated_file(self):
#        fetch = afa.FetchAssemblies(argv=['-p', 'ERP003634'])
#        fetch.config['trustedBrokers'] = ['EMG']
#        run_data = [
#            {'fasta_ftp': '', 'submitted_ftp': 'datafile', 'broker_name': 'EMG'},
#            {'fasta_ftp': 'datafile', 'submitted_ftp': 'datafile', 'broker_name': 'NOT_EMG'},
#        ]
#        assert 1 == len(fetch._filter_ftp_broker_names(run_data))

    def test_map_datafields_ftp_2_db_should_map_all_fields(self):
        raw_data = {
            'secondary_study_accession': 'ERP104225',
            'secondary_sample_accession': 'ERS599830',
            'analysis_accession': 'ERZ477685',
            'submitted_ftp': '/tmp/ERP104225/ERZ477685.fasta.gz',
            'submitted_md5': 'md51'
        }
        fetch = afa.FetchAssemblies(argv=['-p', 'ERP003634'])
        assert fetch.map_datafields_ftp_2_db(deepcopy(raw_data)) == {
            'STUDY_ID': raw_data['secondary_study_accession'],
            'SAMPLE_ID': raw_data['secondary_sample_accession'],
            'ANALYSIS_ID': raw_data['analysis_accession'],
            'DATA_FILE_PATH': ('/tmp/ERP104225/ERZ477685.fasta.gz',),
            'MD5': ('md51',),
            'file': ['ERZ477685.fasta.gz'],
        }

    def test_retrieve_project_info_ftp_should_fetch_all_assemblies(self, tmpdir):
        fetch = afa.FetchAssemblies(argv=['-p', 'ERP112670', '-d', str(tmpdir)])
        assemblies = fetch._retrieve_project_info_ftp('ERP112670')
        assert len(assemblies) == 3

    def mock_get_assembly_metadata(self, *args, **kwargs):
        return [
            {'STUDY_ID': 'ERP104225', 'PROJECT_ID': 'PRJEB22544', 'SAMPLE_ID': 'SRS591962', 'ANALYSIS_ID': 'ERZ477682',
             'DATA_FILE_FORMAT': 'FASTA', 'DATA_FILE_PATH': 'ERZ477/ERZ477682/SRR1238172.fasta.gz', 'TAX_ID': '410658',
             'BYTES': 9499, 'MD5': '9d8192b159ac14b80d68e0537f88cb88'},
            {'STUDY_ID': 'ERP104225', 'PROJECT_ID': 'PRJEB22544', 'SAMPLE_ID': 'SRS591979', 'ANALYSIS_ID': 'ERZ477684',
             'DATA_FILE_FORMAT': 'FASTA', 'DATA_FILE_PATH': 'ERZ477/ERZ477685/SRR1238384.fasta.gz', 'TAX_ID': '410658',
             'BYTES': 8149, 'MD5': 'f7355303727b9a4508e66b1e1d06bbb0'}]

    def mock_get_study_wgs_analyses(self, *args, **kwargs):
        return [{'PROJECT_ACC': 'PRJEB22544', 'SAMPLE_ID': 'SRS591976', 'ASSEMBLY_ID': 'ERZ477684', 'WGS_ACC': 'OEEN01',
                 'GC_ID': 'GCA_900216895', 'CONTIG_CNT': 220},
                {'PROJECT_ACC': 'PRJEB22544', 'SAMPLE_ID': 'SRS591962', 'ASSEMBLY_ID': 'ERZ477682', 'WGS_ACC': 'OEEM01',
                 'GC_ID': 'GCA_900216875', 'CONTIG_CNT': 223}]

    @patch('src.fetch_assemblies.FetchAssemblies._get_assembly_metadata')
    @patch('src.fetch_assemblies.FetchAssemblies._get_study_wgs_analyses')
    @patch('src.fetch_assemblies.FetchAssemblies._get_studies_brokers')
#    @patch('src.fetch_assemblies.FetchAssemblies._study_has_permitted_broker')
    def test_retrieve_project_info_db_should_merge_data_sources(self, mocked_class1, mocked_class2, mocked_class3,
                                                                    mocked_class_4, tmpdir):
        afa.FetchAssemblies._get_assembly_metadata = self.mock_get_assembly_metadata
        afa.FetchAssemblies._get_study_wgs_analyses = self.mock_get_study_wgs_analyses
#        afa.FetchAssemblies._study_has_permitted_broker = lambda *args, **kwargs: False
        afa.FetchAssemblies._get_studies_brokers = lambda *args, **kwargs: {'ERP003634': ''}
        fetch = afa.FetchAssemblies(argv=['-p', 'ERP003634', '-d', str(tmpdir), '--private'])
        assemblies = fetch._retrieve_project_info_db('ERP003634')
        assert len(assemblies) == 2

    @patch('src.fetch_assemblies.ERADAO.retrieve_study_accessions_from_analyses')
    @patch('src.fetch_assemblies.FetchAssemblies._get_assembly_metadata')
    @patch('src.fetch_assemblies.FetchAssemblies._get_study_wgs_analyses')
    @patch('src.fetch_assemblies.FetchAssemblies._get_studies_brokers')
    def test_process_additional_args_should_find_study_accessions_for_assemblies(self, mocked_class1, mocked_class2,
                                                                                 mocked_class3, mocked_class4, tmpdir):
        study_accession = 'ERP104225'
        analysis_id = 'ERZ477684'
        afa.ERADAO.retrieve_study_accessions_from_analyses = lambda *args, **kwargs: [{'STUDY_ID': study_accession}]
        afa.FetchAssemblies._get_assembly_metadata = self.mock_get_assembly_metadata
        afa.FetchAssemblies._get_study_wgs_analyses = self.mock_get_study_wgs_analyses
#        afa.FetchAssemblies._study_has_permitted_broker = lambda *args, **kwargs: True
        afa.FetchAssemblies._get_studies_brokers = lambda *args, **kwargs: {'ERP003634': ''}
        fetch = afa.FetchAssemblies(argv=['-as', analysis_id, '-d', str(tmpdir), '--private'])
        fetch._process_additional_args()
        assert fetch.assemblies == [analysis_id]
        assert fetch.args.projects == [study_accession]

    def test_process_additional_args_should_raise_exception_if_no_private_flag(self, tmpdir):
        with pytest.raises(NotImplementedError):
            afa.FetchAssemblies(argv=['-as', 'ERZ477684', '-d', str(tmpdir)])
