import os
import pytest
import sys
from pathlib import Path

from unittest.mock import patch

from src import fetch_assemblies

FIXTURES_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, 'fixtures'))


class TestFetchAssemblies:
    def test_argparse_should_include_additional_args(self):
        fetch = fetch_assemblies.FetchAssemblies(argv=['-p', 'ERP104225'])

        args = fetch.args
        accepted_args = {'projects', 'project_list', 'dir', 'verbose', 'force',
                         'private', 'interactive', 'config_file', 'assemblies', 'assembly_list', 'assembly_type',
                         'fix_desc_file', 'ignore_errors'}
        assert set(vars(args)) == accepted_args

    def test_validate_args_should_raise_exception_as_no_data_specified(self):
        with pytest.raises(ValueError):
            fetch_assemblies.FetchAssemblies(argv=['-f'])

    def test_process_additional_args_should_set_assemblies_from_arglist(self):
        assemblies = ['ERZ477685']
        fetch = fetch_assemblies.FetchAssemblies(argv=['-p', 'ERP104225', '-as'] + assemblies)
        assert fetch.assemblies == assemblies

    def test_process_additional_args_args_should_set_assemblies_from_file(self, tmpdir):
        assemblies = ['ERZ477685']
        tmpdir = str(tmpdir)
        assembly_file = os.path.join(tmpdir, 'assemblies.txt')
        with open(assembly_file, 'w') as f:
            f.write('\n'.join(assemblies))
        fetch = fetch_assemblies.FetchAssemblies(argv=['-p', 'ERP104225', '--assembly-list', assembly_file])
        assert fetch.assemblies == assemblies

    def test_process_additional_args_args_should_set_assemblies_from_type(self):
        mags = ['ERZ1069976']
        type = 'Metagenome-Assembled Genome (MAG)'
        fetch = fetch_assemblies.FetchAssemblies(argv=['-p', 'ERP116715', '--assembly-type', type, '-as'] + mags)
        assert fetch.assemblies == mags

    def test_filter_accessions_from_args_should_return_empty(self):
        fetch = fetch_assemblies.FetchAssemblies(argv=['-p', 'ERP104225'])
        fetch.assemblies = 'ERZ477685'
        assert [] == fetch._filter_accessions_from_args([], 'analysis_id')

    def test_filter_accessions_from_args_should_return_filtered_assemblies(self):
        fetch = fetch_assemblies.FetchAssemblies(argv=['-p', 'ERP104225'])
        assemblies = ['ERZ477685', 'ERZ477686']
        fetch.assemblies = assemblies
        run_data = [
            {'analysis_id': 'ERZ477685'},
            {'analysis_id': 'ERZ477686'},
            {'analysis_id': 'ERZ477687'},
        ]
        assert run_data[0:2] == fetch._filter_accessions_from_args(run_data, 'analysis_id')

    def test_map_project_info_to_row_should_copy_fields(self):
        raw_data = {
            'STUDY_ID': 'ERP104225',
            'SAMPLE_ID': 'ERS599830',
            'ANALYSIS_ID': 'ERZ477685',
            'file': 'ERZ477685.fasta.gz',
            'DATA_FILE_PATH': '/tmp/ERP104225/ERZ477685.fasta.gz',
            'MD5': ('md51',),
        }
        fetch = fetch_assemblies.FetchAssemblies(argv=['-p', 'ERP003634'])
        transform = fetch.map_project_info_to_row(raw_data)
        equivalent_fields = (
            ('STUDY_ID', 'study_id'),
            ('SAMPLE_ID', 'sample_id'),
            ('ANALYSIS_ID', 'analysis_id'),
        )
        for f1, f2 in equivalent_fields:
            assert raw_data[f1] == transform[f2]
        assert transform['file'] == 'ERZ477685.fasta.gz'

    @staticmethod
    def mock_get_study_from_assembly(self, *args, **kwargs):
        return [{'analysis_accession': 'ERZ1505406', 'secondary_study_accession': 'ERP123564'}]

    """
    1. INVALID = incorrect file format
    2. INVALID = REFERENCE ALIGNMENT should be SEQUENCE ASSEMBLY
    3. INVALID = no file paths
    4. VALID
    """

    @staticmethod
    def mock_get_assembly_metadata(self, *args, **kwargs):
        return [
            {'analysis_accession': 'ERZ1505403', 'study_accession': 'PRJEB39980',
             'analysis_type': 'SEQUENCE_ASSEMBLY', 'center_name': 'EMG', 'first_public': '2020-08-22',
             'submitted_md5': 'bbf202d0482c7bb359e5afe9157b578d;3b0213d90cc8902e29817d69bad524a8',
             'submitted_ftp': 'ftp.sra.ebi.ac.uk/vol1/ERZ150/ERZ1505403/ERR1654124.sub.gz;ftp.sra.ebi.ac.uk'
                              '/vol1/ERZ150/ERZ1505403/ERZ1505403.md5',
             'generated_md5': 'be8676c333d8d075c32db548273a9707',
             'generated_ftp': 'ftp.sra.ebi.ac.uk/vol1/sequence/ERZ150/ERZ1505404/contig.txt.gz',
             'assembly_type': 'primary metagenome', 'secondary_study_accession': 'ERP123564',
             'secondary_sample_accession': 'ERS1360485'},
            {'analysis_accession': 'ERZ1505404', 'study_accession': 'PRJEB39980',
             'analysis_type': 'REFERENCE_ALIGNMENT', 'center_name': 'EMG', 'first_public': '2020-08-22',
             'submitted_md5': 'bbf202d0482c7bb359e5afe9157b578d;3b0213d90cc8902e29817d69bad524a9',
             'submitted_ftp': 'ftp.sra.ebi.ac.uk/vol1/ERZ150/ERZ1505404/ERR1654124.fasta.gz;ftp.sra.ebi.ac.uk'
                              '/vol1/ERZ150/ERZ1505404/ERZ1505404.md5',
             'generated_md5': 'be8676c333d8d075c32db548273a9708',
             'generated_ftp': 'ftp.sra.ebi.ac.uk/vol1/sequence/ERZ150/ERZ1505404/contig.fa.gz',
             'assembly_type': 'primary metagenome', 'secondary_study_accession': 'ERP123564',
             'secondary_sample_accession': 'ERS1360485'},
            {'analysis_accession': 'ERZ1505405', 'study_accession': 'PRJEB39980',
             'analysis_type': 'SEQUENCE_ASSEMBLY', 'center_name': 'EMG', 'first_public': '2020-08-22',
             'submitted_md5': '99bac814522144e55c30369db7a6bf64;3fb2045485655f1faf735bc0a2919b92',
             'submitted_ftp': 'ftp.sra.ebi.ac.uk/vol1/ERZ150/ERZ1505405/ERR1654123.fasta.gz;ftp.sra.ebi.ac.uk'
                              '/vol1/ERZ150/ERZ1505405/ERZ1505405.md5',
             'generated_md5': '',
             'generated_ftp': '',
             'assembly_type': 'primary metagenome', 'secondary_study_accession': 'ERP123564',
             'secondary_sample_accession': 'ERS1360484'},
            {'analysis_accession': 'ERZ1505406', 'study_accession': 'PRJEB39980',
             'analysis_type': 'SEQUENCE_ASSEMBLY', 'center_name': 'EMG', 'first_public': '2020-08-22',
             'submitted_md5': '589552f30ffa7927c33cd121cab59de1;d63816a798413d5b8c361b87319c9b08',
             'submitted_ftp': 'ftp.sra.ebi.ac.uk/vol1/ERZ150/ERZ1505406/ERZ1505406.md5;ftp.sra.ebi.ac.uk/'
                              'vol1/ERZ150/ERZ1505406/ERR1654119.fasta.gz',
             'generated_md5': 'e0e3c99075ae482083b3e8ca35e401e2',
             'generated_ftp': 'ftp.sra.ebi.ac.uk/vol1/sequence/ERZ150/ERZ1505406/contig.fa.gz',
             'assembly_type': 'primary metagenome', 'secondary_study_accession': 'ERP123564',
             'secondary_sample_accession': 'ERS1360480'}
        ]

    @patch('src.fetch_assemblies.FetchAssemblies._retrieve_ena_url')
    def test_process_additional_args_should_find_study_accessions_for_assemblies(self, mocked_class1, tmpdir):
        study_accession = 'ERP123564'
        analysis_id = 'ERZ1505406'
        fetch_assemblies.FetchAssemblies._retrieve_ena_url = self.mock_get_study_from_assembly
        fetch = fetch_assemblies.FetchAssemblies(argv=['-as', analysis_id, '-d', str(tmpdir)])
        fetch._validate_args()
        fetch._process_additional_args()
        assert fetch.args.projects == {study_accession}

    @patch('src.fetch_assemblies.FetchAssemblies._retrieve_ena_url')
    def test_retrieve_project_should_return_only_valid_assemblies_and_check_md5(self, mocked_class1, tmpdir):
        study_accession = 'ERP123564'
        valid_assembly_file = 'ERZ1505406.fasta.gz'
        valid_assembly_md5 = 'e0e3c99075ae482083b3e8ca35e401e2'
        download_files = ['download', 'download.lock', 'ERP123564.txt', 'ERP123564.txt.lock']
        fetch_assemblies.FetchAssemblies._retrieve_ena_url = self.mock_get_assembly_metadata
        fetch_assemblies.FetchAssemblies.download_lftp = True
        fetch = fetch_assemblies.FetchAssemblies(argv=['-p', study_accession, '-d', str(tmpdir), '--private'])
        assemblies = fetch._retrieve_project_info_from_api(study_accession)
        for x in assemblies:
            for file in x['file']:
                assembly_path = tmpdir / file
                Path(str(assembly_path)).touch()
        assert len(assemblies) == 1
        assert os.listdir(str(tmpdir)) == [valid_assembly_file]
        assert not fetch._is_file_valid(str(assembly_path), valid_assembly_md5)
        project_dir = tmpdir / study_accession
        os.mkdir(str(project_dir))
        os.chdir(str(tmpdir))
        fetch.write_project_files(study_accession, assemblies)
        assert [os.path.exists(str(project_dir / x)) for x in download_files]
        with open(str(project_dir / 'download')) as f:
            download_data = f.readlines()
            assert len(download_data) == 1
        with open(str(project_dir / 'ERP123564.txt')) as t:
            txt_data = t.readlines()
            assert len(txt_data) == 2

    @patch.object(fetch_assemblies.FetchAssemblies, 'fetch')
    def test_main_should_call_fetch(self, mock):
        test_args = ['scriptname', '-p', 'ERP123564']
        with patch.object(sys, 'argv', test_args):
            fetch_assemblies.main()
        assert mock.called
