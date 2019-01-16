import pytest
import tempfile
from subprocess import Popen, PIPE
import logging
import os
import re
import unittest
import shutil

FIXTURES_DIR = os.path.join(os.path.dirname(__file__), 'fixtures')


def exec_fetch_cmd(study_id, dir, runs=None):
    cmd = 'fetch-read-tool -p {} -d {} -v '.format(study_id, str(dir))
    if runs:
        cmd += '-ru ' + " ".join(runs)
    p = Popen(cmd, stdout=PIPE, stderr=PIPE, shell=True)
    p.wait()
    stdout, stderr = p.communicate()
    logging.info(stdout)
    logging.error(stderr)
    assert p.returncode == 0


def check_all_files_downloaded(raw_dir, expected_files):
    for f in expected_files:
        f_path = os.path.join(raw_dir, f)
        assert os.path.getsize(f_path) > 0
    assert len(expected_files) == len(os.listdir(raw_dir))


class PrepDownloadTest(unittest.TestCase):
    runs = None

    @classmethod
    def setUpClass(cls):
        cls.tmpdir = str(tempfile.mkdtemp())
        cls.study_id = 'ERP003634'
        cls.study_dir = None
        cls.project_file = None
        print(cls.tmpdir)

        exec_fetch_cmd(cls.study_id, cls.tmpdir, cls.runs)
        cls.study_dir = os.path.join(cls.tmpdir, cls.study_id)
        cls.project_file = os.path.join(cls.study_dir, cls.study_id + '.txt')
        assert os.path.exists(cls.project_file)

    @classmethod
    def tearDownClass(cls):
        # shutil.rmtree(cls.tmpdir)
        pass


class TestFetchReads(PrepDownloadTest):
    # Fetch all runs in 1 go and assert that they were all fetched correctly

    def test_download_project_file_should_contain_correct_data(self):
        with open(self.project_file) as f:
            project_data = f.readlines()
        expected_project_file = os.path.join(FIXTURES_DIR, self.study_id + '.txt')
        with open(expected_project_file) as f:
            expected_project_data = f.readlines()

        assert set(project_data) == set(expected_project_data)

    def test_download_project_data_should_fetch_non_empty_files(self):
        with open(self.project_file) as f:
            project_data = f.read()
        raw_files = re.findall('(\wRR\d+(?:_1|_2)?\.fastq\.gz)', project_data)
        raw_dir = os.path.join(self.study_dir, 'raw')
        check_all_files_downloaded(raw_dir, raw_files)


class TestFetchReadsRuns(PrepDownloadTest):
    # Fetch 3 runs in 1 go and assert that they were all fetched correctly

    @classmethod
    def setUpClass(cls):
        cls.runs = ['ERR315850', 'ERR315851', 'ERR315852']
        super(TestFetchReadsRuns, cls).setUpClass()
        cls.project_file = os.path.join(FIXTURES_DIR, cls.study_id + '_runs.txt')

    def test_download_project_file_has_only_specified_runs(self):
        with open(self.project_file) as f:
            project_data = set(f.readlines())

        downloaded_project_data = os.path.join(self.study_dir, self.study_id + '.txt')
        with open(downloaded_project_data) as f:
            downloaded_project_data = set(f.readlines())
        assert project_data == downloaded_project_data


class TestFetchReadsRunsSequence(PrepDownloadTest):
    # Sequentially fetch 3 different runs and assert that they were all fetched correctly

    @classmethod
    def setUpClass(cls):
        cls.runs = ['ERR315850']
        super(TestFetchReadsRunsSequence, cls).setUpClass()
        cls.project_file = os.path.join(FIXTURES_DIR, cls.study_id + '_runs.txt')
        cls.runs = ['ERR315851']
        exec_fetch_cmd(cls.study_id, cls.tmpdir, cls.runs)
        cls.runs = ['ERR315852']
        exec_fetch_cmd(cls.study_id, cls.tmpdir, cls.runs)

    def test_download_project_file_has_only_specified_runs(self):
        with open(self.project_file) as f:
            project_data = set(f.readlines())

        downloaded_project_data = os.path.join(self.study_dir, self.study_id + '.txt')
        with open(downloaded_project_data) as f:
            downloaded_project_data = set(f.readlines())
        assert project_data == downloaded_project_data
