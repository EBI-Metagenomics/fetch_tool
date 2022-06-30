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

import gzip
import json
import logging
import os
import re

from fetchtool.abstract_fetch import AbstractDataFetcher
from fetchtool.exceptions import ENAFetch204, NoDataError

path_re = re.compile(r"(.*)/(.*)")


class Analysis(object):
    def __init__(self, analysis_accession, assembly_type, status_id):
        self.analysis_accession = analysis_accession
        self.assembly_type = assembly_type
        self.status_id = status_id


class FetchAssemblies(AbstractDataFetcher):
    ENA_PORTAL_BASE_API_URL = "https://www.ebi.ac.uk/ena/portal/api/search?"

    ENA_PORTAL_FIELDS = [
        "analysis_accession",
        "study_accession",
        "secondary_study_accession",
        "sample_accession",
        "secondary_sample_accession",
        "analysis_title",
        "analysis_type",
        "center_name",
        "first_public",
        "last_updated",
        "study_title",
        "analysis_alias",
        "study_alias",
        "submitted_md5",
        "submitted_ftp",
        "generated_md5",
        "generated_ftp",
        "sample_alias",
        "broker_name",
        "sample_title",
        "assembly_type",
    ]

    ENA_PORTAL_RUN_FIELDS = "secondary_study_accession"

    ENA_PORTAL_PARAMS = [
        "dataPortal=metagenome",
        "dccDataOnly=false",
        "result=analysis",
        "format=json",
        "download=true",
        "fields=",
    ]

    # query
    ENA_PORTAL_QUERY = (
        "query=secondary_study_accession=%22{0}%22%20AND%20assembly_type=%22{1}%22"
    )
    ENA_PORTAL_RUN_QUERY = (
        "query=analysis_accession=%22{0}%22%20AND%20assembly_type=%22{1}%22"
    )

    ENA_PORTAL_API_URL = (
        ENA_PORTAL_BASE_API_URL
        + "&".join(ENA_PORTAL_PARAMS)
        + ",".join(ENA_PORTAL_FIELDS)
        + "&"
        + ENA_PORTAL_QUERY
    )
    ENA_PORTAL_API_BY_RUN = (
        ENA_PORTAL_BASE_API_URL
        + "&".join(ENA_PORTAL_PARAMS)
        + ENA_PORTAL_RUN_FIELDS
        + "&"
        + ENA_PORTAL_RUN_QUERY
    )

    ENA_FILEREPORT_URL = "https://www.ebi.ac.uk/ena/portal/api/filereport"

    # not in use
    ENA_WGS_SET_API_URL = (
        "https://www.ebi.ac.uk/ena/portal/api/search?result=wgs_set&query=study_accession="
        "%22{0}%22&fields=study_accession,assembly_type,scientific_name,fasta_file&format=tsv"
    )

    def __init__(self, argv=None):
        self.ACCESSION_FIELD = "ANALYSIS_ID"
        self.assemblies = None
        super().__init__(argv)

    @staticmethod
    def add_arguments(parser):
        assembly_group = parser.add_mutually_exclusive_group()
        assembly_group.add_argument(
            "-as",
            "--assemblies",
            nargs="+",
            help="Assembly ERZ accession(s), whitespace separated. "
            "Use to download only certain project assemblies",
        )
        # TODO: Also the older INSDC style metagenome assemblies produced by MGnify are not support that way at moment
        parser.add_argument(
            "--assembly-type",
            help="Assembly type",
            choices=[
                "primary metagenome",
                "binned metagenome",
                "metatranscriptome",
            ],
            default="primary metagenome",
        )
        assembly_group.add_argument(
            "--assembly-list", help="File containing line-separated assembly accessions"
        )

        return parser

    def _validate_args(self):
        if not any(
            [
                self.args.assemblies,
                self.args.assembly_list,
                self.args.projects,
                self.args.project_list,
            ]
        ):
            raise ValueError(
                "No data specified, please use -as, --assembly-list, -p or --project-list"
            )

    def _process_additional_args(self):
        self.assembly_type = self.args.assembly_type

        if self.args.assembly_list:
            self.assemblies = self._read_line_sep_file(self.args.assembly_list)
        else:
            self.assemblies = self.args.assemblies

        if not self.args.projects and not self.args.project_list:
            logging.info("Fetching projects from list of assemblies")
            self.args.projects = self._get_project_accessions_from_assemblies(
                self.assemblies
            )

    def _retrieve_project_info_from_api(self, project_accession):
        raise_error = (
            False if len(self.projects) > 1 else True
        )  # allows script to continue to next project if one fails

        data = None
        fetch_204_ex = None
        try:
            data = self._retrieve_ena_url(
                self.ENA_PORTAL_API_URL.format(project_accession, self.assembly_type),
                raise_on_204=raise_error,
            )
        except ENAFetch204 as ex:
            logging.info(
                "The data portal API returned a 204, trying the filereport API"
            )
            fetch_204_ex = ex

        # Alternative endpoint, the filereport which is the used on the ENA website
        if fetch_204_ex:
            file_report_url = (
                self.ENA_FILEREPORT_URL
                + "?"
                + "&".join(
                    [
                        "&".join(self.ENA_PORTAL_PARAMS)
                        + ",".join(self.ENA_PORTAL_FIELDS),
                        f"accession={project_accession}",
                    ]
                )
            )
            data = self._retrieve_ena_url(file_report_url)
            if not data:
                logging.error(
                    f"It was not possible to fetch data from the Portal API or the Filereport API for project {project_accession}"
                )
                return

        if not data:
            logging.error(
                f"It was not possible to fetch data from the Portal API for project {project_accession}"
            )
            return

        logging.info(
            "Retrieved {count} assemblies for study {project_accession} from "
            "the ENA Portal API.".format(
                count=len(data), project_accession=project_accession
            )
        )
        mapped_data = []
        for d in data:
            if not d["generated_ftp"]:
                logging.info(
                    "The generated ftp location for assembly {} is not available yet".format(
                        d["analysis_accession"]
                    )
                )
            if d["analysis_type"] == "SEQUENCE_ASSEMBLY" and d["generated_ftp"]:
                if self._is_rawdata_filetype(
                    os.path.basename(d["generated_ftp"])
                ):  # filter filenames not fasta
                    raw_data_file_path, file_, md5_ = self._get_raw_filenames(
                        d.get("generated_ftp"),
                        d.get("generated_md5"),
                        d.get("analysis_accession"),
                        bool(d.get("submitted_ftp")),
                    )
                    mapped_data.append(
                        {
                            "STUDY_ID": d.get("secondary_study_accession"),
                            "SAMPLE_ID": d.get("secondary_sample_accession"),
                            "ANALYSIS_ID": d.get("analysis_accession"),
                            "DATA_FILE_PATH": raw_data_file_path,
                            "file": file_,
                            "MD5": md5_,
                        }
                    )
        return mapped_data

    def _filter_accessions_from_args(self, assembly_data, assembly_accession_field):
        if self.assemblies:
            data = list(
                filter(
                    lambda r: (r[assembly_accession_field] in self.assemblies),
                    assembly_data,
                )
            )
            return data
        else:
            return assembly_data

    def map_project_info_to_row(self, assembly):
        return {
            "study_id": assembly["STUDY_ID"],
            "sample_id": assembly["SAMPLE_ID"],
            "analysis_id": assembly["ANALYSIS_ID"],
            "file": assembly["file"],
            "file_path": assembly["DATA_FILE_PATH"],
            "scientific_name": "n/a",
            "md5": assembly["MD5"],
        }

    def _get_project_accessions_from_assemblies(self, assemblies):
        project_list = set()
        for assembly in assemblies:
            data = self._retrieve_ena_url(
                self.ENA_PORTAL_API_BY_RUN.format(assembly, self.assembly_type),
                raise_on_204=False,
            )
            if data:
                [project_list.add(d["secondary_study_accession"]) for d in data]
        if not len(project_list):
            raise NoDataError(self.NO_DATA_MSG)
        return project_list


def main():
    data_fetcher = FetchAssemblies()
    data_fetcher.fetch()


if __name__ == "__main__":
    main()
