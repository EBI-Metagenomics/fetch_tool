import os

import subprocess

import unittest

import pytest

FIXTURES_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, 'fixtures'))


def call_cmd(cmd):
    print(cmd)
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


study_id = 'ERP104225'


def validate_full_study(tmpdir):
    study_dir = os.path.join(tmpdir, study_id[0:7], study_id)
    assert os.path.exists(study_dir)
    study_file = os.path.join(study_dir, study_id + '.txt')
    with open(study_file) as f:
        study_data = f.readlines()
    assert study_data[0] == 'study_id\tsample_id\tanalysis_id\tfiles\tfile_path\n'
    assert study_data[1] == 'ERP104225\tSRS591962\tERZ477682\tERZ477682.fasta.gz\t' \
                            'ftp.sra.ebi.ac.uk/vol1/ERZ477/ERZ477682/SRR1238172.fasta.gz\n'
    assert study_data[2] == 'ERP104225\tSRS591965\tERZ477683\tERZ477683.fasta.gz\t' \
                            'ftp.sra.ebi.ac.uk/vol1/ERZ477/ERZ477683/SRR1238319.fasta.gz\n'
    assert study_data[3] == 'ERP104225\tSRS591976\tERZ477684\tERZ477684.fasta.gz\t' \
                            'ftp.sra.ebi.ac.uk/vol1/ERZ477/ERZ477684/SRR1238351.fasta.gz\n'
    assert study_data[4] == 'ERP104225\tSRS591979\tERZ477685\tERZ477685.fasta.gz\t' \
                            'ftp.sra.ebi.ac.uk/vol1/ERZ477/ERZ477685/SRR1238384.fasta.gz\n'

    download_file = os.path.join(study_dir, 'download')
    with open(download_file) as f:
        download_data = f.readlines()
    assert len(download_data) == 4

    raw_dir = os.path.join(study_dir, 'raw')
    expected_files = ['ERZ477682.fasta.gz', 'ERZ477683.fasta.gz', 'ERZ477684.fasta.gz', 'ERZ477685.fasta.gz', ]
    for f in expected_files:
        f_path = os.path.join(raw_dir, f)
        assert os.path.getsize(f_path) > 0


@pytest.mark.skipif(os.environ.get('TRAVIS') == 'true', reason='Skipped as running on TravisCI')
class TestFetchCompleteStudyAssemblies:
    def test_fetch_all_study_data(self, tmpdir):
        with WorkingDir(tmpdir):
            call_cmd('fetch-assembly-tool -p {} -v -v -d {}'.format(study_id, str(tmpdir)))
            validate_full_study(tmpdir)
