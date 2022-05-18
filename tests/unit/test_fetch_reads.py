import os
import pytest
import sys
from pathlib import Path

from unittest.mock import patch

from src import fetch_reads

FIXTURES_DIR = os.path.abspath(
    os.path.join(os.path.dirname(__file__), os.pardir, "fixtures")
)


class TestFetchReads:
    def test_argparse_should_include_additional_args(self):
        fetch = fetch_reads.FetchReads(argv=["-p", "ERP001736"])

        args = fetch.args
        accepted_args = {
            "projects",
            "project_list",
            "dir",
            "verbose",
            "force",
            "private",
            "interactive",
            "config_file",
            "runs",
            "run_list",
            "fix_desc_file",
            "ignore_errors",
        }
        assert set(vars(args)) == accepted_args

    def test_validate_args_should_raise_exception_as_no_data_specified(self):
        with pytest.raises(ValueError):
            fetch_reads.FetchReads(argv=["-f"])

    def test_process_additional_args_should_set_runs_from_arglist(self):
        runs = ["ERR599038", "ERR599039"]
        fetch = fetch_reads.FetchReads(argv=["-p", "ERP001736", "-ru"] + runs)
        assert fetch.runs == runs

    def test_process_additional_args_args_should_set_runs_from_file(self, tmpdir):
        runs = ["ERR599038", "ERR599039"]
        tmpdir = str(tmpdir)
        runfile = os.path.join(tmpdir, "runs.txt")
        with open(runfile, "w") as f:
            f.write("\n".join(runs))
        fetch = fetch_reads.FetchReads(argv=["-p", "ERP001736", "--run-list", runfile])
        assert fetch.runs == runs

    def test_filter_by_accessions_should_return_empty(self):
        fetch = fetch_reads.FetchReads(argv=["-p", "ERP001736"])
        fetch.runs = "ERR599830"
        assert [] == fetch.filter_by_accessions([])

    def test_filter_accessions_from_args_should_return_filtered_runs(self):
        fetch = fetch_reads.FetchReads(argv=["-p", "ERP001736"])
        runs = ["ERR599830", "ERR599831"]
        fetch.runs = runs
        run_data = [
            {"run_id": "ERR599830"},
            {"run_id": "ERR599831"},
            {"run_id": "ERR599832"},
        ]
        assert run_data[0:2] == fetch._filter_accessions_from_args(run_data, "run_id")

    def test_map_project_info_db_row_should_copy_fields(self):
        raw_data = {
            "STUDY_ID": "ERP001736",
            "SAMPLE_ID": "ERS599830",
            "RUN_ID": "ERR599830",
            "LIBRARY_LAYOUT": "PAIRED",
            "LIBRARY_SOURCE": "METAGENOMIC",
            "LIBRARY_STRATEGY": "WGS",
            "file": "ERR599383_1.fastq.gz;ERR599383_2.fastq.gz",
            "DATA_FILE_PATH": "/tmp/ERP001736/ERR599383_1.fastq.gz;/tmp/ERP001736/ERR599383_2.fastq.gz",
        }
        fetch = fetch_reads.FetchReads(argv=["-p", "ERP003634"])
        transform = fetch.map_project_info_to_row(raw_data)
        equivalent_fields = (
            ("STUDY_ID", "study_id"),
            ("SAMPLE_ID", "sample_id"),
            ("RUN_ID", "run_id"),
            ("LIBRARY_LAYOUT", "library_layout"),
            ("LIBRARY_STRATEGY", "library_strategy"),
            ("LIBRARY_SOURCE", "library_source"),
        )
        for f1, f2 in equivalent_fields:
            assert raw_data[f1] == transform[f2]
        assert transform["file"] == "ERR599383_1.fastq.gz;ERR599383_2.fastq.gz"

    @staticmethod
    def mock_get_study_from_run(self, *args, **kwargs):
        return [
            {"run_accession": "ERR2777789", "secondary_study_accession": "ERP110686"}
        ]

    """
    1. INVALID = incorrect file format
    2. INVALID = no file paths
    3. VALID
    """

    @staticmethod
    def mock_get_run_metadata(self, *args, **kwargs):
        return [
            {
                "study_accession": "PRJEB28479",
                "secondary_study_accession": "ERP110686",
                "sample_accession": "SAMEA4883561",
                "secondary_sample_accession": "ERS2702567",
                "experiment_accession": "ERX2789866",
                "run_accession": "ERR2777788",
                "instrument_model": "unspecified",
                "library_layout": "PAIRED",
                "fastq_ftp": "ftp.sra.ebi.ac.uk/vol1/fastq/ERR277/009/ERR2777788/ERR2777788_1.txt.gz;"
                "ftp.sra.ebi.ac.uk/vol1/fastq/ERR277/009/ERR2777790/ERR2777788_2.txt.gz",
                "fastq_md5": "39f9956b66880e386d741eea2a0e54c2;9e6db19a2ef56383e8e426784ffff425",
                "submitted_ftp": "ftp.sra.ebi.ac.uk/vol1/run/ERR277/ERR2777788/140210.050.upload.txt.trim.gz",
                "submitted_md5": "7935d13d964cc6bc5038f7706ec3e1c5",
                "library_strategy": "AMPLICON",
                "broker_name": "MGRAST",
                "library_source": "METAGENOMIC",
            },
            {
                "study_accession": "PRJEB28479",
                "secondary_study_accession": "ERP110686",
                "sample_accession": "SAMEA4883561",
                "secondary_sample_accession": "ERS2702567",
                "experiment_accession": "ERX2789866",
                "run_accession": "ERR2777789",
                "instrument_model": "unspecified",
                "library_layout": "PAIRED",
                "fastq_ftp": "",
                "fastq_md5": "",
                "submitted_ftp": "ftp.sra.ebi.ac.uk/vol1/run/ERR277/ERR2777789/140210.050.upload.fna.trim.gz",
                "submitted_md5": "7935d13d964cc6bc5038f7706ec3e1c4",
                "library_strategy": "AMPLICON",
                "broker_name": "MGRAST",
                "library_source": "METAGENOMIC",
            },
            {
                "study_accession": "PRJEB28479",
                "secondary_study_accession": "ERP110686",
                "sample_accession": "SAMEA4883562",
                "secondary_sample_accession": "ERS2702568",
                "experiment_accession": "ERX2789867",
                "run_accession": "ERR2777790",
                "instrument_model": "unspecified",
                "library_layout": "PAIRED",
                "fastq_ftp": "ftp.sra.ebi.ac.uk/vol1/fastq/ERR277/009/ERR2777790/ERR2777790_1.fastq.gz;"
                "ftp.sra.ebi.ac.uk/vol1/fastq/ERR277/009/ERR2777790/ERR2777790_2.fastq.gz",
                "fastq_md5": "39f9956b66880e386d741eea2a0e54c1;9e6db19a2ef56383e8e426784ffff424",
                "submitted_ftp": "ftp.sra.ebi.ac.uk/vol1/run/ERR277/ERR2777790/140211.050.upload.fna.trim.gz:"
                "ftp.sra.ebi.ac.uk/vol1/run/ERR277/ERR2777790/140211.050.upload.fna.trim.gz",
                "submitted_md5": "39f9956b66880e386d741eea2a0e54c1;9e6db19a2ef56383e8e426784ffff424",
                "library_strategy": "AMPLICON",
                "broker_name": "MGRAST",
                "library_source": "METAGENOMIC",
            },
        ]

    @patch("src.fetch_reads.FetchReads._retrieve_ena_url")
    def test_process_additional_args_should_find_study_accessions_for_runs(
        self, mocked_class1, tmpdir
    ):
        study_accession = "ERP110686"
        run_id = "ERR2777789"
        fetch_reads.FetchReads._retrieve_ena_url = self.mock_get_study_from_run
        fetch = fetch_reads.FetchReads(argv=["-ru", run_id, "-d", str(tmpdir)])
        fetch._validate_args()
        fetch._process_additional_args()
        assert fetch.args.projects == {study_accession}

    @patch("src.fetch_reads.FetchReads._retrieve_ena_url")
    def test_retrieve_project_should_return_only_valid_reads_and_check_md5(
        self, mocked_class1, tmpdir
    ):
        study_accession = "ERP110686"
        valid_file_for = ("ERR2777790_1.fastq.gz", "39f9956b66880e386d741eea2a0e54c1")
        valid_file_rev = ("ERR2777790_2.fastq.gz", "9e6db19a2ef56383e8e426784ffff424")
        download_files = [
            "download",
            "download.lock",
            "ERP110686.txt",
            "ERP110686.txt.lock",
        ]
        fetch_reads.FetchReads._retrieve_ena_url = self.mock_get_run_metadata
        fetch_reads.FetchReads.download_lftp = True
        fetch = fetch_reads.FetchReads(
            argv=["-p", study_accession, "-d", str(tmpdir), "--private"]
        )
        runs = fetch._retrieve_project_info_from_api(study_accession)
        for x in runs:
            for file in x["file"]:
                run_path = tmpdir / file
                Path(str(run_path)).touch()
        assert len(runs) == 1
        assert (
            os.listdir(str(tmpdir)).sort()
            == ["ERR2777790_2.fastq.gz", "ERR2777790_1.fastq.gz"].sort()
        )
        for x, y in [valid_file_for, valid_file_rev]:
            assert not fetch._is_file_valid(str(tmpdir / x), y)
        project_dir = tmpdir / study_accession
        os.mkdir(str(project_dir))
        os.chdir(str(tmpdir))
        fetch.write_project_files(study_accession, runs)
        assert [os.path.exists(str(project_dir / x)) for x in download_files]
        with open(str(project_dir / "download")) as f:
            download_data = f.readlines()
            assert len(download_data) == 2
        with open(str(project_dir / "ERP110686.txt")) as t:
            txt_data = t.readlines()
            assert len(txt_data) == 2

    @patch.object(fetch_reads.FetchReads, "fetch")
    def test_main_should_call_fetch(self, mock):
        test_args = ["scriptname", "-p", "ERP110686"]
        with patch.object(sys, "argv", test_args):
            fetch_reads.main()
        assert mock.called
