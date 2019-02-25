import os
import pytest
import shutil

from copy import deepcopy

from unittest.mock import patch

from src import abstract_fetch_reads as afr

FIXTURES_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, 'fixtures'))


class TestFetchReads:
    def test_argparse_should_include_additional_args(self):
        fetch = afr.FetchReads(argv=['-p', 'ERP001736'])

        args = fetch.args
        accepted_args = ['projects', 'project_list', 'dir', 'verbose', 'force',
                         'private', 'interactive', 'config_file', 'runs', 'run_list']
        assert all([hasattr(args, argname) for argname in accepted_args])

    def test_validate_args_should_raise_exception_as_no_data_specified(self):
        with pytest.raises(ValueError):
            afr.FetchReads(argv=['-f'])

    def test_process_additional_args_should_set_runs_from_arglist(self):
        runs = ['ERR599038', 'ERR599039']
        fetch = afr.FetchReads(argv=['-p', 'ERP001736', '-ru'] + runs)
        assert fetch.runs == runs

    def test_process_additional_args_args_should_set_runs_from_file(self, tmpdir):
        runs = ['ERR599038', 'ERR599039']
        tmpdir = str(tmpdir)
        runfile = os.path.join(tmpdir, 'runs.txt')
        with open(runfile, 'w') as f:
            f.write('\n'.join(runs))
        fetch = afr.FetchReads(argv=['-p', 'ERP001736', '--run-list', runfile])
        assert fetch.runs == runs

    def test_validate_args_should_raise_error_if_fetching_studies_without_project_in_public_mode(self, tmpdir):
        runs = ['ERR599038', 'ERR599039']
        tmpdir = str(tmpdir)
        runfile = os.path.join(tmpdir, 'runs.txt')
        with open(runfile, 'w') as f:
            f.write('\n'.join(runs))

        with pytest.raises(NotImplementedError):
            afr.FetchReads(argv=['--run-list', runfile])

    def test_get_project_accessions_from_runs_should_return_empty(self):
        assert [] == afr.FetchReads._get_study_run_accessions([])

    def test_get_project_accessions_from_runs_should_raise_keyerror_on_invalid_project_data(self):
        with pytest.raises(KeyError):
            afr.FetchReads._get_study_run_accessions([{'PROJECT_ID': 'ERP001736'}])

    def test_get_project_accessions_from_runs_should_return_run_id(self):
        run_id = 'ERR599830'
        assert [run_id] == afr.FetchReads._get_study_run_accessions([{'RUN_ID': run_id}])

    def test_filter_runs_from_existing_downloads(self):
        pass

    def test_filter_runs_from_args_should_return_empty(self):
        fetch = afr.FetchReads(argv=['-p', 'ERP001736'])
        fetch.runs = 'ERR599830'
        assert [] == fetch._filter_runs_from_args([], 'run_id')

    def test_filter_runs_from_args_should_return_filtered_runs(self):
        fetch = afr.FetchReads(argv=['-p', 'ERP001736'])
        runs = ['ERR599830', 'ERR599831']
        fetch.runs = runs
        run_data = [
            {'run_id': 'ERR599830'},
            {'run_id': 'ERR599831'},
            {'run_id': 'ERR599832'},
        ]
        assert run_data[0:2] == fetch._filter_runs_from_args(run_data, 'run_id')

    def test_filter_runs_from_existing_downloads_should_not_filter_as_no_file_present(self):
        fetch = afr.FetchReads(argv=['-p', 'ERP001736'])
        # Assert description file does not exist
        with pytest.raises(FileNotFoundError):
            fetch.read_project_description_file('ERP001736')
        run_data = [
            {'run_id': 'ERR599830'},
            {'run_id': 'ERR599831'},
            {'run_id': 'ERR599832'},
        ]
        assert run_data == fetch._filter_runs_from_existing_downloads('ERP001736', run_data, 'ERR599083')

    def test_filter_runs_from_existing_downloads_should_filter_using_description_file(self, tmpdir):
        tmpdir = str(tmpdir)
        project_dir = os.path.join(FIXTURES_DIR, 'ERP0036')
        shutil.copytree(project_dir, tmpdir + '/ERP0036')
        fetch = afr.FetchReads(argv=['-p', 'ERP003634', '-d', tmpdir])
        new_runs = [{'run_id': 'ERR315856'}]
        assert new_runs == fetch._filter_runs_from_existing_downloads('ERP003634', new_runs, 'run_id')

    def test_map_project_info_db_row_should_copy_fields(self):
        raw_data = {
            'STUDY_ID': 'ERP001736',
            'SAMPLE_ID': 'ERS599830',
            'RUN_ID': 'ERR599830',
            'LIBRARY_LAYOUT': 'PAIRED',
            'TAX_ID': '1231',
            'LIBRARY_SOURCE': 'METAGENOMIC',
            'LIBRARY_STRATEGY': 'WGS',
            'DATA_FILE_PATH': '/tmp/ERP001736/ERR599383_1.fastq.gz;/tmp/ERP001736/ERR599383_2.fastq.gz',
            'DATA_FILE_ROLE': 'SUBMITTED',
        }
        fetch = afr.FetchReads(argv=['-p', 'ERP003634'])
        transform = fetch.map_project_info_db_row(raw_data)
        equivalent_fields = (
            ('STUDY_ID', 'study_id'),
            ('SAMPLE_ID', 'sample_id'),
            ('RUN_ID', 'run_id'),
            ('LIBRARY_LAYOUT', 'library_layout'),
            ('LIBRARY_STRATEGY', 'library_strategy'),
            ('LIBRARY_SOURCE', 'library_source'),
        )
        for f1, f2 in equivalent_fields:
            assert raw_data[f1] == transform[f2]

    def test_map_project_info_db_row_should_generate_file_list(self):
        raw_data = {
            'STUDY_ID': 'ERP001736',
            'SAMPLE_ID': 'ERS599830',
            'RUN_ID': 'ERR599830',
            'LIBRARY_LAYOUT': 'PAIRED',
            'TAX_ID': '1231',
            'LIBRARY_SOURCE': 'METAGENOMIC',
            'LIBRARY_STRATEGY': 'WGS',
            'DATA_FILE_PATH': '/tmp/ERP001736/ERR599383_1.fastq.gz;/tmp/ERP001736/ERR599383_2.fastq.gz',
            'DATA_FILE_ROLE': 'SUBMITTED',
        }
        fetch = afr.FetchReads(argv=['-p', 'ERP003634'])
        transform = fetch.map_project_info_db_row(raw_data)
        assert transform['file'] == 'ERR599383_1.fastq.gz;ERR599383_2.fastq.gz'

    def test_is_trusted_ftp_data_should_be_trusted_from_broker(self):
        fetch = afr.FetchReads(argv=['-p', 'ERP003634'])
        trusted_brokers = ['EMG']
        assert fetch.is_trusted_ftp_data({'fastq_ftp': '',
                                          'submitted_ftp': 'filepath',
                                          'broker_name': 'EMG'},
                                         trusted_brokers)

    def test_is_trusted_ftp_data_should_be_trusted_as_generated_data(self):
        fetch = afr.FetchReads(argv=['-p', 'ERP003634'])
        trusted_brokers = ['EMG']
        assert fetch.is_trusted_ftp_data({'fastq_ftp': 'filepath', 'submitted_ftp': 'filepath', 'broker_name': 'EMG'},
                                         trusted_brokers)

    def test_filter_ftp_broker_names_should_return_empty(self):
        fetch = afr.FetchReads(argv=['-p', 'ERP003634'])
        assert [] == fetch._filter_ftp_broker_names([])

    def test_filter_ftp_broker_names_should_filter_brokers(self):
        fetch = afr.FetchReads(argv=['-p', 'ERP003634'])
        fetch.config['trustedBrokers'] = ['EMG']
        run_data = [
            {'fastq_ftp': '', 'submitted_ftp': 'datafile', 'broker_name': 'EMG'},
            {'fastq_ftp': '', 'submitted_ftp': 'datafile', 'broker_name': 'NOT_EMG'},
        ]
        assert 'EMG' == fetch._filter_ftp_broker_names(run_data)[0]['broker_name']

    def test_filter_ftp_broker_names_should_allow_generated_files(self):
        fetch = afr.FetchReads(argv=['-p', 'ERP003634'])
        fetch.config['trustedBrokers'] = ['EMG']
        run_data = [
            {'fastq_ftp': '', 'submitted_ftp': 'datafile', 'broker_name': 'EMG'},
            {'fastq_ftp': 'datafile', 'submitted_ftp': 'datafile', 'broker_name': 'NOT_EMG'},
        ]
        assert 2 == len(fetch._filter_ftp_broker_names(run_data))

    def test_map_datafields_ftp_2_db_should_map_all_fields(self):
        raw_data = {
            'secondary_study_accession': 'ERP001736',
            'secondary_sample_accession': 'ERS599830',
            'run_accession': 'ERR599830',
            'tax_id': '1231',
            'library_source': 'METAGENOMIC',
            'library_strategy': 'WGS',
            'library_layout': 'PAIRED',
            'fastq_ftp': '/tmp/ERP001736/ERR599383_1.fastq.gz;/tmp/ERP001736/ERR599383_2.fastq.gz',
            'fastq_md5': 'md51;md52',
            'submitted_ftp': '',
            'submitted_md5': ''
        }
        fetch = afr.FetchReads(argv=['-p', 'ERP003634'])
        assert fetch.map_datafields_ftp_2_db(deepcopy(raw_data)) == {
            'STUDY_ID': raw_data['secondary_study_accession'],
            'SAMPLE_ID': raw_data['secondary_sample_accession'],
            'RUN_ID': raw_data['run_accession'],
            'TAX_ID': raw_data['tax_id'],
            'LIBRARY_SOURCE': raw_data['library_source'],
            'LIBRARY_STRATEGY': raw_data['library_strategy'],
            'LIBRARY_LAYOUT': raw_data['library_layout'],
            'DATA_FILE_ROLE': 'GENERATED_FILE',
            'DATA_FILE_PATH': raw_data['fastq_ftp'],
            'MD5': raw_data['fastq_md5'],
            'files': ['ERR599383_1.fastq.gz', 'ERR599383_2.fastq.gz'],
        }

    def test_retrieve_project_info_ftp_should_fetch_all_runs(self, tmpdir):
        tmpdir = str(tmpdir)
        fetch = afr.FetchReads(argv=['-p', 'ERP003634', '-d', tmpdir])
        runs = fetch._retrieve_project_info_ftp('ERP110686')
        assert len(runs) == 2

    def mock_db_response_generated_data(self, *args, **kwargs):
        return [{'STUDY_ID': 'ERP113309', 'SAMPLE_ID': 'ERS3033973', 'RUN_ID': 'ERR3063510', 'LIBRARY_LAYOUT': 'SINGLE',
                 'DATA_FILE_PATH': 'fastq/ERR306/000/ERR3063510/ERR3063510.fastq.gz', 'TAX_ID': '256318',
                 'LIBRARY_STRATEGY': 'WGS', 'LIBRARY_SOURCE': 'GENOMIC', 'DATA_FILE_ROLE': 'GENERATED_FILE',
                 'MD5': '2cb9441bbc157d76286b385fd0b7f0c4'},
                {'STUDY_ID': 'ERP113309', 'SAMPLE_ID': 'ERS3042514', 'RUN_ID': 'ERR3086151', 'LIBRARY_LAYOUT': 'SINGLE',
                 'DATA_FILE_PATH': 'fastq/ERR308/001/ERR3086151/ERR3086151.fastq.gz', 'TAX_ID': '256318',
                 'LIBRARY_STRATEGY': 'WGS', 'LIBRARY_SOURCE': 'GENOMIC', 'DATA_FILE_ROLE': 'GENERATED_FILE',
                 'MD5': '899eb5ab522ebc2c98bc567f3c3ad7d8'}
                ]

    def mock_db_response_submitted_data(self, *args, **kwargs):
        return [{'STUDY_ID': 'ERP113309', 'SAMPLE_ID': 'ERS3042515', 'RUN_ID': 'ERR3086152', 'LIBRARY_LAYOUT': 'SINGLE',
                 'DATA_FILE_PATH': 'fastq/ERR308/001/ERR3086151/ERR3086151.fastq.gz', 'TAX_ID': '256318',
                 'LIBRARY_STRATEGY': 'WGS', 'LIBRARY_SOURCE': 'GENOMIC', 'DATA_FILE_ROLE': 'SUBMITTED_FILE',
                 'MD5': '899eb5ab522ebc2c98bc567f3c3ad7d8'}]

    @patch('src.abstract_fetch_reads.FetchReads._retrieve_era_generated_data')
    @patch('src.abstract_fetch_reads.FetchReads._study_has_permitted_broker')
    def test_retrieve_project_info_db_should_not_use_generated_data(self, mocked_class1, mocked_class2, tmpdir):
        tmpdir = str(tmpdir)
        afr.FetchReads._retrieve_era_generated_data = self.mock_db_response_generated_data
        afr.FetchReads._study_has_permitted_broker = lambda *args, **kwargs: False
        fetch = afr.FetchReads(argv=['-p', 'ERP113309', '-d', tmpdir, '--private'])
        runs = fetch._retrieve_project_info_db('ERP113309')
        assert len(runs) == 2

    @patch('src.abstract_fetch_reads.FetchReads._retrieve_era_generated_data')
    @patch('src.abstract_fetch_reads.FetchReads._retrieve_era_submitted_data')
    @patch('src.abstract_fetch_reads.FetchReads._study_has_permitted_broker')
    def test_retrieve_project_info_db_should_add_submitted_data(self, mocked_class1, mocked_class2,
                                                                mocked_class3, tmpdir):
        tmpdir = str(tmpdir)
        afr.FetchReads._retrieve_era_generated_data = self.mock_db_response_generated_data
        afr.FetchReads._retrieve_era_submitted_data = self.mock_db_response_submitted_data
        afr.FetchReads._study_has_permitted_broker = lambda *args, **kwargs: True
        fetch = afr.FetchReads(argv=['-p', 'ERP113309', '-d', tmpdir, '--private'])
        runs = fetch._retrieve_project_info_db('ERP113309')
        assert len(runs) == 3
