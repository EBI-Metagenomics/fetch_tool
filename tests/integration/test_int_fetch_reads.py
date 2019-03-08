import os

import subprocess


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


study_id = 'ERP110686'


def validate_full_study(tmpdir):
    study_dir = os.path.join(tmpdir, study_id)
    assert os.path.exists(study_dir)
    study_file = os.path.join(study_dir, study_id + '.txt')
    with open(study_file) as f:
        study_data = f.readlines()
    assert study_data[0] == 'study_id\tsample_id\trun_id\tlibrary_layout\tfile\tfile_path\ttax_id\t' \
                            'scientific_name\tlibrary_strategy\tlibrary_source\n'
    assert study_data[1] == 'ERP110686\tERS2702567\tERR2777789\tSINGLE\tERR2777789.fasta.gz\t' \
                            'ftp.sra.ebi.ac.uk/vol1/run/ERR277/ERR2777789/140210.050.upload.fna.trim.gz\t' \
                            '256318\tn/a\tAMPLICON\tMETAGENOMIC\n'
    assert study_data[2] == 'ERP110686\tERS2702568\tERR2777790\tSINGLE\tERR2777790.fasta.gz\t' \
                            'ftp.sra.ebi.ac.uk/vol1/run/ERR277/ERR2777790/140211.050.upload.fna.trim.gz\t' \
                            '256318\tn/a\tAMPLICON\tMETAGENOMIC\n'

    download_file = os.path.join(study_dir, 'download')
    with open(download_file) as f:
        download_data = f.readlines()
    assert len(download_data) == 2

    raw_dir = os.path.join(study_dir, 'raw')
    expected_file = ['ERR2777789.fasta.gz', 'ERR2777790.fasta.gz']
    for f in expected_file:
        f_path = os.path.join(raw_dir, f)
        assert os.path.getsize(f_path) > 0


def validate_single_study_run(tmpdir):
    tmpdir = str(tmpdir)
    study_dir = os.path.join(tmpdir, study_id)
    assert os.path.exists(study_dir)
    study_file = os.path.join(study_dir, study_id + '.txt')
    with open(study_file) as f:
        study_data = f.readlines()
    assert study_data[0] == 'study_id\tsample_id\trun_id\tlibrary_layout\tfile\tfile_path\ttax_id\t' \
                            'scientific_name\tlibrary_strategy\tlibrary_source\n'
    assert study_data[1] == 'ERP110686\tERS2702567\tERR2777789\tSINGLE\tERR2777789.fasta.gz\t' \
                            'ftp.sra.ebi.ac.uk/vol1/run/ERR277/ERR2777789/140210.050.upload.fna.trim.gz\t' \
                            '256318\tn/a\tAMPLICON\tMETAGENOMIC\n'

    download_file = os.path.join(study_dir, 'download')
    with open(download_file) as f:
        download_data = f.readlines()
    assert len(download_data) == 1

    raw_dir = os.path.join(study_dir, 'raw')
    expected_file = ['ERR2777789.fasta.gz']
    for f in expected_file:
        f_path = os.path.join(raw_dir, f)
        assert os.path.getsize(f_path) > 0

@pytest.mark.skipif(os.environ.get('TRAVIS') == 'true', reason='Skipped as running on TravisCI')
class TestFetchCompleteStudyReads:
    def test_fetch_all_study_data(self, tmpdir):
        with WorkingDir(tmpdir):
            call_cmd('fetch-read-tool -p {} -v -v -d {}'.format(study_id, str(tmpdir)))
            validate_full_study(tmpdir)

    def test_fetch_single_run(self, tmpdir):
        with WorkingDir(tmpdir):
            call_cmd('fetch-read-tool -p {} -ru {} -d {}'.format(study_id, 'ERR2777789', str(tmpdir)))
            validate_single_study_run(tmpdir)

    def test_fetch_sequential_runs(self, tmpdir):
        with WorkingDir(tmpdir):
            call_cmd('fetch-read-tool -p {} -ru {} -d {}'.format(study_id, 'ERR2777789', str(tmpdir)))
            validate_single_study_run(tmpdir)
            call_cmd('fetch-read-tool -p {} -ru {} -d {}'.format(study_id, 'ERR2777790', str(tmpdir)))
            validate_full_study(tmpdir)
