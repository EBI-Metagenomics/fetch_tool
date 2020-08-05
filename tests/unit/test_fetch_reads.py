import json
import os
import pytest
import shutil
import sys

from copy import deepcopy

from unittest.mock import patch

from src import fetch_reads as afr, abstract_fetch

FIXTURES_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, 'fixtures'))


class TestFetchReads:
    def test_argparse_should_include_additional_args(self):
        fetch = afr.FetchReads(argv=['-p', 'ERP001736'])

        args = fetch.args
        accepted_args = {'projects', 'project_list', 'dir', 'verbose', 'force',
                         'private', 'interactive', 'config_file', 'runs', 'run_list', 'fix_desc_file', 'ignore_errors'}
        assert set(vars(args)) == accepted_args

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

    def test_filter_by_accessions_should_return_empty(self):
        fetch = afr.FetchReads(argv=['-p', 'ERP001736'])
        fetch.runs = 'ERR599830'
        assert [] == fetch.filter_by_accessions([])

    def test_filter_by_accessions_should_skip_with_force(self):
        fetch = afr.FetchReads(argv=['-p', 'ERP001736', '-f'])
        run_accession = 'ERR599830'
        run_data = [{'run_id': run_accession}]
        assert run_data == fetch.filter_by_accessions(run_data)

    def test_filter_by_accessions_should_not_filter_new_run(self):
        fetch = afr.FetchReads(argv=['-p', 'ERP001736'])
        fetch.runs = ['ERR599830', 'ERR599831']
        run_data = [{'RUN_ID': 'ERR599831'}]
        assert run_data == fetch.filter_by_accessions(run_data)

    def test_filter_accessions_from_args_should_return_empty(self):
        fetch = afr.FetchReads(argv=['-p', 'ERP001736'])
        fetch.runs = 'ERR599830'
        assert [] == fetch._filter_accessions_from_args([], 'run_id')

    def test_filter_accessions_from_args_should_return_filtered_runs(self):
        fetch = afr.FetchReads(argv=['-p', 'ERP001736'])
        runs = ['ERR599830', 'ERR599831']
        fetch.runs = runs
        run_data = [
            {'run_id': 'ERR599830'},
            {'run_id': 'ERR599831'},
            {'run_id': 'ERR599832'},
        ]
        assert run_data[0:2] == fetch._filter_accessions_from_args(run_data, 'run_id')

    def test_map_project_info_db_row_should_copy_fields(self):
        raw_data = {
            'STUDY_ID': 'ERP001736',
            'SAMPLE_ID': 'ERS599830',
            'RUN_ID': 'ERR599830',
            'LIBRARY_LAYOUT': 'PAIRED',
            'LIBRARY_SOURCE': 'METAGENOMIC',
            'LIBRARY_STRATEGY': 'WGS',
            'file': 'ERR599383_1.fastq.gz;ERR599383_2.fastq.gz',
            'DATA_FILE_PATH': '/tmp/ERP001736/ERR599383_1.fastq.gz;/tmp/ERP001736/ERR599383_2.fastq.gz',
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
            'LIBRARY_SOURCE': 'METAGENOMIC',
            'LIBRARY_STRATEGY': 'WGS',
            'file': 'ERR599383_1.fastq.gz;ERR599383_2.fastq.gz',
            'DATA_FILE_PATH': '/tmp/ERP001736/ERR599383_1.fastq.gz;/tmp/ERP001736/ERR599383_2.fastq.gz',
            'DATA_FILE_ROLE': 'SUBMITTED',
        }
        fetch = afr.FetchReads(argv=['-p', 'ERP003634'])
        transform = fetch.map_project_info_db_row(raw_data)
        assert transform['file'] == 'ERR599383_1.fastq.gz;ERR599383_2.fastq.gz'

#    def test_is_trusted_ftp_data_should_be_trusted_from_broker(self):
#        fetch = afr.FetchReads(argv=['-p', 'ERP003634'])
#        trusted_brokers = ['EMG']
#        assert fetch.is_trusted_ftp_data({'fastq_ftp': '',
#                                          'submitted_ftp': 'filepath',
#                                          'broker_name': 'EMG'},
#                                         trusted_brokers)

#    def test_is_trusted_ftp_data_should_be_trusted_as_generated_data(self):
#        fetch = afr.FetchReads(argv=['-p', 'ERP003634'])
#        trusted_brokers = ['EMG']
#        assert fetch.is_trusted_ftp_data({'fastq_ftp': 'filepath', 'submitted_ftp': 'filepath', 'broker_name': 'EMG'},
#                                         trusted_brokers)

#    def test_filter_ftp_broker_names_should_return_empty(self):
#        fetch = afr.FetchReads(argv=['-p', 'ERP003634'])
#        assert [] == fetch._filter_ftp_broker_names([])

#    def test_filter_ftp_broker_names_should_filter_brokers(self):
#        fetch = afr.FetchReads(argv=['-p', 'ERP003634'])
#        fetch.config['trustedBrokers'] = ['EMG']
#        run_data = [
#            {'fastq_ftp': '', 'submitted_ftp': 'datafile', 'broker_name': 'EMG'},
#            {'fastq_ftp': '', 'submitted_ftp': 'datafile', 'broker_name': 'NOT_EMG'},
#        ]
#        assert 1 == len(fetch._filter_ftp_broker_names(run_data))

#    def test_filter_ftp_broker_names_should_allow_generated_file(self):
#        fetch = afr.FetchReads(argv=['-p', 'ERP003634'])
#        fetch.config['trustedBrokers'] = ['EMG']
#        run_data = [
#            {'fastq_ftp': '', 'submitted_ftp': 'datafile', 'broker_name': 'EMG'},
#            {'fastq_ftp': 'datafile', 'submitted_ftp': 'datafile', 'broker_name': 'NOT_EMG'},
#        ]
#        assert 2 == len(fetch._filter_ftp_broker_names(run_data))

    def test_map_datafields_ftp_2_db_should_map_all_fields(self):
        raw_data = {
            'secondary_study_accession': 'ERP001736',
            'secondary_sample_accession': 'ERS599830',
            'run_accession': 'ERR599830',
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
            'LIBRARY_SOURCE': raw_data['library_source'],
            'LIBRARY_STRATEGY': raw_data['library_strategy'],
            'LIBRARY_LAYOUT': raw_data['library_layout'],
            'DATA_FILE_ROLE': 'GENERATED_FILE',
            'DATA_FILE_PATH': ('/tmp/ERP001736/ERR599383_1.fastq.gz', '/tmp/ERP001736/ERR599383_2.fastq.gz'),
            'MD5': ('md51', 'md52'),
            'file': ['ERR599383_1.fastq.gz', 'ERR599383_2.fastq.gz'],
        }

    def test_retrieve_project_info_ftp_should_fetch_all_runs(self, tmpdir):
        fetch = afr.FetchReads(argv=['-p', 'ERP003634', '-d', str(tmpdir)])
        runs = fetch._retrieve_project_info_ftp('ERP110686')
        assert len(runs) == 2

    def mock_db_response_generated_data(self, *args, **kwargs):
        return [{'STUDY_ID': 'ERP113309', 'SAMPLE_ID': 'ERS3033973', 'RUN_ID': 'ERR3063510', 'LIBRARY_LAYOUT': 'SINGLE',
                 'DATA_FILE_PATH': 'fastq/ERR306/000/ERR3063510/ERR3063510.fastq.gz',
                 'LIBRARY_STRATEGY': 'WGS', 'LIBRARY_SOURCE': 'GENOMIC', 'DATA_FILE_ROLE': 'GENERATED_FILE',
                 'MD5': '2cb9441bbc157d76286b385fd0b7f0c4'},
                {'STUDY_ID': 'ERP113309', 'SAMPLE_ID': 'ERS3042514', 'RUN_ID': 'ERR3086151', 'LIBRARY_LAYOUT': 'SINGLE',
                 'DATA_FILE_PATH': 'fastq/ERR308/001/ERR3086151/ERR3086151.fastq.gz',
                 'LIBRARY_STRATEGY': 'WGS', 'LIBRARY_SOURCE': 'GENOMIC', 'DATA_FILE_ROLE': 'GENERATED_FILE',
                 'MD5': '899eb5ab522ebc2c98bc567f3c3ad7d8'}
                ]

    def mock_db_response_submitted_data(self, *args, **kwargs):
        return [{'STUDY_ID': 'ERP113309', 'SAMPLE_ID': 'ERS3042515', 'RUN_ID': 'ERR3086152', 'LIBRARY_LAYOUT': 'SINGLE',
                 'DATA_FILE_PATH': 'fastq/ERR308/001/ERR3086151/ERR3086151.fastq.gz',
                 'LIBRARY_STRATEGY': 'WGS', 'LIBRARY_SOURCE': 'GENOMIC', 'DATA_FILE_ROLE': 'SUBMITTED_FILE',
                 'MD5': '899eb5ab522ebc2c98bc567f3c3ad7d8'}]

    @patch('src.fetch_reads.FetchReads._retrieve_era_generated_data')
#    @patch('src.fetch_reads.FetchReads._get_studies_brokers')
    @patch('src.fetch_reads.FetchReads._retrieve_era_submitted_data')
#    @patch('src.fetch_reads.FetchReads._study_has_permitted_broker')
    def test_retrieve_project_info_db_should_not_use_generated_data(self, mocked_class1, mocked_class2, tmpdir):
        afr.FetchReads._retrieve_era_generated_data = self.mock_db_response_generated_data
        afr.FetchReads._retrieve_era_submitted_data = self.mock_db_response_submitted_data
#        afr.FetchReads._study_has_permitted_broker = lambda *args, **kwargs: False
        afr.FetchReads._get_studies_brokers = lambda *args, **kwargs: {'ERP113309': ''}
        fetch = afr.FetchReads(argv=['-p', 'ERP113309', '-d', str(tmpdir), '--private'])
        runs = fetch._retrieve_project_info_db('ERP113309')
        assert len(runs) == 2

    @patch('src.fetch_reads.FetchReads._retrieve_era_generated_data')
    @patch('src.fetch_reads.FetchReads._retrieve_era_submitted_data')
#    @patch('src.fetch_reads.FetchReads._study_has_permitted_broker')
    def test_retrieve_project_info_db_should_add_submitted_data(self, mocked_class1, mocked_class2, tmpdir):
        afr.FetchReads._retrieve_era_generated_data = self.mock_db_response_generated_data
        afr.FetchReads._retrieve_era_submitted_data = self.mock_db_response_submitted_data
#        afr.FetchReads._study_has_permitted_broker = lambda *args, **kwargs: True
        afr.FetchReads._get_studies_brokers = lambda *args, **kwargs: {'ERP113309': ''}
        fetch = afr.FetchReads(argv=['-p', 'ERP113309', '-d', str(tmpdir), '--private'])
        runs = fetch._retrieve_project_info_db('ERP113309')
        assert len(runs) == 3

    @patch('src.fetch_reads.ERADAO.retrieve_study_accessions_from_runs')
    def test_process_additional_args_should_find_study_accessions_for_runs(self, mocked_class1, tmpdir):
        study_accession = 'ERP001736'
        run_id = 'ERR599083'
        afr.ERADAO.retrieve_study_accessions_from_runs = lambda *args, **kwargs: [{'STUDY_ID': study_accession}]
        afr.FetchReads._retrieve_era_generated_data = self.mock_db_response_generated_data
        afr.FetchReads._retrieve_era_submitted_data = self.mock_db_response_submitted_data
#        afr.FetchReads._study_has_permitted_broker = lambda *args, **kwargs: True
        fetch = afr.FetchReads(argv=['-ru', run_id, '-d', str(tmpdir), '--private'])
        fetch._process_additional_args()
        assert fetch.runs == [run_id]
        assert fetch.args.projects == [study_accession]

    def test_process_additional_args_should_raise_exception_if_no_private_flag(self, tmpdir):
        with pytest.raises(NotImplementedError):
            afr.FetchReads(argv=['-ru', 'ERR599083', '-d', str(tmpdir)])

    def test_write_project_files_should_create_both_file(self, tmpdir):
        tmpdir = str(tmpdir)
        study_accession = 'ERP110686'
        project_dir = os.path.join(tmpdir, study_accession)
        os.makedirs(project_dir)
        run_data = [
            {'study_accession': 'PRJEB28479', 'sample_accession': 'SAMEA4883561', 'experiment_accession': 'ERX2789866',
             'instrument_model': 'unspecified', 'broker_name': 'MGRAST',
             'STUDY_ID': 'ERP110686', 'SAMPLE_ID': 'ERS2702567', 'RUN_ID': 'ERR2777789',
             'DATA_FILE_ROLE': 'SUBMISSION_FILE',
             'DATA_FILE_PATH': ('ftp.sra.ebi.ac.uk/vol1/run/ERR277/ERR2777789/140210.050.upload.fna.trim.gz',),
             'file': ['ERR2777789.fasta.gz'], 'MD5': ('7935d13d964cc6bc5038f7706ec3e1c4',),
             'LIBRARY_STRATEGY': 'AMPLICON', 'LIBRARY_SOURCE': 'METAGENOMIC', 'LIBRARY_LAYOUT': 'SINGLE'}]
        fetch = afr.FetchReads(argv=['-p', study_accession, '-ru', 'ERR2777789', '-d', tmpdir])
        fetch.write_project_files(study_accession, run_data)
        project_desc_file = os.path.join(project_dir, study_accession + '.txt')
        with open(project_desc_file) as f:
            data = f.readlines()
        assert len(data) == 2
        assert data[0] == '\t'.join(afr.FetchReads.DEFAULT_HEADERS) + '\n'
        assert data[1] == 'ERP110686	ERS2702567	ERR2777789	n/a	SINGLE	AMPLICON	METAGENOMIC	' \
                          'ERR2777789.fasta.gz	' \
                          'ftp.sra.ebi.ac.uk/vol1/run/ERR277/ERR2777789/140210.050.upload.fna.trim.gz\n'
        download_file = os.path.join(project_dir, 'download')
        with open(download_file) as f:
            download_data = f.readlines()
        assert len(download_data) == 1

    @patch.object(afr.FetchReads, 'fetch')
    def test_main_should_call_fetch(self, mock):
        test_args = ['scriptname', '-p', 'ERP001736']
        with patch.object(sys, 'argv', test_args):
            afr.main()
        assert mock.called
