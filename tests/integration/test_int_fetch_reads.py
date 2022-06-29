#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright 2018-2022 EMBL - European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import csv
import os
import subprocess

import pytest

FIXTURES_DIR = os.path.abspath(
    os.path.join(os.path.dirname(__file__), os.pardir, "fixtures")
)


root_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))


def call_cmd(cmd):
    # FIXME: this pythonpath hack is _very_ ugly
    ret = subprocess.call(
        f"PYTHONPATH={root_path}:$PYTHONPATH {cmd}", stdout=subprocess.PIPE, shell=True
    )
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


study_id = "ERP110686"


def validate_full_study(tmpdir):
    study_dir = os.path.join(tmpdir, study_id)
    assert os.path.exists(study_dir)
    study_file = os.path.join(study_dir, study_id + ".txt")
    with open(study_file) as f:
        study_data = f.readlines()

    rows = list(csv.reader(study_data, delimiter="\t"))

    assert rows[0] == [
        "study_id",
        "sample_id",
        "run_id",
        "analysis_id",
        "library_layout",
        "library_strategy",
        "library_source",
        "file",
        "file_path",
    ]
    assert rows[1] == [
        "ERP110686",
        "ERS2702567",
        "ERR2777789",
        "n/a",
        "SINGLE",
        "AMPLICON",
        "METAGENOMIC",
        "ERR2777789.fasta.gz",
        "ftp.sra.ebi.ac.uk/vol1/run/ERR277/ERR2777789/140210.050.upload.fna.trim.gz",
    ]
    assert rows[2] == [
        "ERP110686",
        "ERS2702568",
        "ERR2777790",
        "n/a",
        "SINGLE",
        "AMPLICON",
        "METAGENOMIC",
        "ERR2777790.fasta.gz",
        "ftp.sra.ebi.ac.uk/vol1/run/ERR277/ERR2777790/140211.050.upload.fna.trim.gz",
    ]

    download_file = os.path.join(study_dir, "download")
    with open(download_file) as f:
        download_data = f.readlines()
    assert len(download_data) == 2

    raw_dir = os.path.join(study_dir, "raw")
    expected_file = ["ERR2777789.fasta.gz", "ERR2777790.fasta.gz"]
    for f in expected_file:
        f_path = os.path.join(raw_dir, f)
        assert os.path.getsize(f_path) > 0


def validate_single_study_run(tmpdir):
    tmpdir = str(tmpdir)
    study_dir = os.path.join(tmpdir, study_id)
    assert os.path.exists(study_dir)
    study_file = os.path.join(study_dir, study_id + ".txt")
    with open(study_file) as f:
        study_data = f.readlines()

    rows = list(csv.reader(study_data, delimiter="\t"))

    assert rows[0] == [
        "study_id",
        "sample_id",
        "run_id",
        "analysis_id",
        "library_layout",
        "library_strategy",
        "library_source",
        "file",
        "file_path",
    ]

    assert rows[1] == [
        "ERP110686",
        "ERS2702567",
        "ERR2777789",
        "n/a",
        "SINGLE",
        "AMPLICON",
        "METAGENOMIC",
        "ERR2777789.fasta.gz",
        "ftp.sra.ebi.ac.uk/vol1/run/ERR277/ERR2777789/140210.050.upload.fna.trim.gz",
    ]

    download_file = os.path.join(study_dir, "download")
    with open(download_file) as f:
        download_data = f.readlines()
    assert len(download_data) == 1

    raw_dir = os.path.join(study_dir, "raw")
    expected_file = ["ERR2777789.fasta.gz"]
    for f in expected_file:
        f_path = os.path.join(raw_dir, f)
        assert os.path.getsize(f_path) > 0


@pytest.mark.skip(reason="FIXME")
class TestFetchCompleteStudyReads:
    def test_fetch_all_study_data(self, tmpdir):
        with WorkingDir(tmpdir):
            call_cmd("fetch_reads.py -p {} -v -v -d {}".format(study_id, str(tmpdir)))
            validate_full_study(tmpdir)

    def test_fetch_single_run(self, tmpdir):
        with WorkingDir(tmpdir):
            call_cmd(f"fetch_reads.py -v -p {study_id} -ru ERR2777789 -d {tmpdir}")
            validate_single_study_run(tmpdir)

    def test_fetch_sequential_runs(self, tmpdir):
        with WorkingDir(tmpdir):
            call_cmd(f"fetch_reads.py -v -p {study_id} -ru ERR2777789 -d {tmpdir}")
            validate_single_study_run(tmpdir)
            call_cmd(f"fetch_reads.py -v -p {study_id} -ru ERR2777790 -d {tmpdir}")
            validate_full_study(tmpdir)
