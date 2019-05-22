import os

import subprocess

import unittest

import pytest

FIXTURES_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, 'fixtures'))


def call_cmd(cmd):
    ret = subprocess.call(cmd, stdout=subprocess.PIPE, shell=True)
    assert ret == 0


class WorkingDir:
    def __init__(self, tmpdir):
        self.tmpdir = str(tmpdir)
        self.prev_dir = None

    def __enter__(self):
        self.prev_dir = os.getcwd()
        os.chdir(self.tmpdir)

    def __exit__(self, exc_type, exc_val, exc_tb):
        os.chdir(self.prev_dir)


study_id = 'ERP111374'


def validate_full_study(tmpdir):
    study_dir = os.path.join(tmpdir, study_id)
    assert os.path.exists(study_dir)
    study_file = os.path.join(study_dir, study_id + '.txt')
    with open(study_file) as f:
        study_data = f.readlines()

    assert study_data[0].strip().split('\t') == \
        ['study_id', 'sample_id', 'run_id', 'analysis_id', 'library_layout',
         'library_strategy', 'library_source', 'file', 'file_path']
    assert study_data[1].strip() == 'ERP111374\tSRS2137356\tn/a\tERZ773362\tn/a\tn/a\tn/a\tERZ773362.fasta.gz\t' \
                            'ftp.sra.ebi.ac.uk/vol1/ERZ773/ERZ773362/SRR5465818.fasta.gz'

    assert study_data[2].strip() == 'ERP111374\tSRS2137355\tn/a\tERZ773363\tn/a\tn/a\tn/a\tERZ773363.fasta.gz\t' \
                            'ftp.sra.ebi.ac.uk/vol1/ERZ773/ERZ773363/SRR5465817.fasta.gz'

    download_file = os.path.join(study_dir, 'download')
    with open(download_file) as f:
        download_data = f.readlines()
    assert len(download_data) == 2

    raw_dir = os.path.join(study_dir, 'raw')
    expected_file = ['ERZ773362.fasta.gz', 'ERZ773363.fasta.gz']
    for f in expected_file:
        f_path = os.path.join(raw_dir, f)
        assert os.path.getsize(f_path) > 0


@pytest.mark.skipif(os.environ.get('TRAVIS') == 'true', reason='Skipped as running on TravisCI')
class TestFetchCompleteStudyAssemblies:
    def test_fetch_all_study_data(self, tmpdir):
        with WorkingDir(tmpdir):
            call_cmd('fetch-assembly-tool -p {} -v -v -d {}'.format(study_id, str(tmpdir)))
            validate_full_study(tmpdir)
