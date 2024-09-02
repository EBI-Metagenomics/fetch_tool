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
from unittest.mock import patch

import pytest

from fetchtool import fetch_assemblies

FIXTURES_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, "fixtures"))


class WorkingDir:
    def __init__(self, tmpdir):
        self.tmpdir = str(tmpdir)
        self.prev_dir = None

    def __enter__(self):
        self.prev_dir = os.getcwd()
        os.chdir(self.tmpdir)

    def __exit__(self, exc_type, exc_val, exc_tb):
        os.chdir(self.prev_dir)


study_id = "ERP111374"


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
        'instrument_model',
        'instrument_platform',
        "file",
        "file_path",
    ]

    assert rows[1] == [
        "ERP111374",
        "SRS2137356",
        "n/a",
        "ERZ773362",
        "n/a",
        "n/a",
        "n/a",
        "n/a",
        "n/a",
        "ERZ773362.fasta.gz",
        # -> old structure, remove after checking with ENA
        # "ftp.sra.ebi.ac.uk/vol1/ERZ773/ERZ773362/SRR5465818.fasta.gz",
        "ftp.sra.ebi.ac.uk/vol1/sequence/ERZ773/ERZ773362/contig.fa.gz",
    ]

    assert rows[2] == [
        "ERP111374",
        "SRS2137355",
        "n/a",
        "ERZ773363",
        "n/a",
        "n/a",
        "n/a",
        "n/a",
        "n/a",
        "ERZ773363.fasta.gz",
        # -> old structure, remove after checking with ENA
        # "ftp.sra.ebi.ac.uk/vol1/ERZ773/ERZ773363/SRR5465817.fasta.gz",
        "ftp.sra.ebi.ac.uk/vol1/sequence/ERZ773/ERZ773363/contig.fa.gz",
    ]

    download_file = os.path.join(study_dir, "download")
    with open(download_file) as f:
        download_data = f.readlines()
    assert len(download_data) == 2

    raw_dir = os.path.join(study_dir, "raw")
    expected_file = ["ERZ773362.fasta.gz", "ERZ773363.fasta.gz"]
    for f in expected_file:
        f_path = os.path.join(raw_dir, f)
        assert os.path.getsize(f_path) > 0


@pytest.mark.flaky
class TestFetchCompleteStudyAssemblies:
    @patch("fetchtool.abstract_fetch.AbstractDataFetcher.download_lftp")
    @patch("fetchtool.abstract_fetch.AbstractDataFetcher.download_rsync")
    def test_fetch_all_study_data(self, lftp_mock, rsync_mock, tmpdir):
        lftp_mock.return_value = False
        rsync_mock.return_value = False

        with WorkingDir(tmpdir):
            fassemblies = fetch_assemblies.FetchAssemblies(["-p", study_id, "-v", "-d", str(tmpdir)])
            fassemblies.fetch()
            validate_full_study(tmpdir)
            lftp_mock.call_count = 2
            rsync_mock.call_count = 2

    @patch("fetchtool.fetch_assemblies.AbstractDataFetcher.download_lftp")
    @patch("fetchtool.fetch_assemblies.AbstractDataFetcher.download_wget")
    def test_fetch_sequential_runs_with_rsync(self, wget_mock, lftp_mock, tmpdir):
        def return_false(*args, **kwargs):
            return False

        lftp_mock.side_effect = return_false
        wget_mock.side_effect = return_false
        with WorkingDir(tmpdir):
            fassemblies = fetch_assemblies.FetchAssemblies(["-p", study_id, "-d", str(tmpdir)])
            fassemblies.fetch()
            validate_full_study(tmpdir)

            assert lftp_mock.called
            assert wget_mock.called is False  # this was it's only called after rsync fails
