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


class ERADAO(object):
    """
    Data access object to ENA's ERA schema instances (predominantly ERAPRO)
    """

    def __init__(self, data_access_object):
        """
        Constructor
        """
        self.data_access_object = data_access_object

    def retrieve_study_accessions_from_runs(self, run_ids):
        """
            Returns a list of study_ids For run accessions in run_ids
        :param run_ids:
        :return:
        """
        runs = ",".join(["'" + run + "'" for run in run_ids])
        query = "SELECT study_id FROM v_mgp_run_file WHERE run_id IN ({}) GROUP BY study_id".format(
            runs
        )
        return self.data_access_object._runQuery(query)

    def retrieve_study_accessions_from_analyses(self, analysis_ids):
        """
            Returns a list of study_ids For run accessions in run_ids
        :param run_ids:
        :return:
        """
        analyses = ",".join(["'" + id + "'" for id in analysis_ids])
        query = "SELECT study_id FROM v_mgp_assembly_file WHERE analysis_id IN ({}) GROUP BY study_id".format(
            analyses
        )
        return self.data_access_object._runQuery(query)

    def retrieve_submitted_files(self, study_id):
        """
            Returns a list of submitted files.
        :type study_id: str
        :return: list[dict]
        """

        query = "SELECT study_id, sample_id, run_id, library_layout, LISTAGG(data_file_path, ';') WITHIN group( order by run_id) as data_file_path, tax_id, library_strategy, library_source, data_file_role, LISTAGG(md5, ';') WITHIN group( order by run_id) as md5 FROM v_mgp_run_file WHERE study_id = '{}' and DATA_FILE_ROLE = 'SUBMISSION_FILE' and (DATA_FILE_FORMAT='FASTA' OR DATA_FILE_FORMAT='FASTQ') and DATA_FILE_PATH not like 'err%' GROUP BY study_id, sample_id, run_id, library_layout, tax_id, library_strategy, library_source, data_file_role".format(
            study_id
        )
        return self.data_access_object._runQuery(query)

    def retrieve_generated_files(self, study_id):
        """
            Returns a list of generated files.
        :type study_id: str
        :return:
        """

        query = "SELECT study_id, sample_id, run_id, library_layout, LISTAGG(data_file_path, ';') WITHIN group( order by run_id) as data_file_path, tax_id, library_strategy, library_source, data_file_role, LISTAGG(md5, ';') WITHIN group( order by run_id) as md5 FROM v_mgp_run_file WHERE (study_id = '{}' and DATA_FILE_ROLE = 'GENERATED_FILE' and DATA_FILE_FORMAT='FASTQ' and DATA_FILE_PATH not like 'err%') GROUP BY study_id, sample_id, run_id, library_layout, tax_id, library_strategy, library_source, data_file_role".format(
            study_id
        )
        return self.data_access_object._runQuery(query)

    def retrieve_assembly_generated_files(self, study_id):
        """
        Returns assembly metadata
        :param study_id:
        :return:
        """

        query = "SELECT study_id, project_id, sample_id, analysis_id, data_file_format, data_file_role, LISTAGG(data_file_path, ';') WITHIN group( order by analysis_id) as data_file_path, tax_id, LISTAGG(md5, ';') WITHIN group( order by analysis_id) as md5 FROM v_mgp_assembly_file WHERE data_file_format = 'FASTA' and DATA_FILE_ROLE = 'GENERATED_FILE' and study_id = '{}' AND ANALYSIS_STATUS!='supressed' AND ANALYSIS_STATUS!='cancelled' GROUP BY study_id, project_id, sample_id, analysis_id, data_file_format, tax_id, data_file_role ".format(
            study_id
        )
        return self.data_access_object._runQuery(query)

    def retrieve_assembly_submitted_files(self, study_id):
        """
        Returns assembly metadata
        :param study_id:
        :return:
        """

        query = "SELECT study_id, project_id, sample_id, analysis_id, data_file_format, data_file_role, LISTAGG(data_file_path, ';') WITHIN group( order by analysis_id) as data_file_path, tax_id, LISTAGG(md5, ';') WITHIN group( order by analysis_id) as md5 FROM v_mgp_assembly_file WHERE data_file_format = 'FASTA' and DATA_FILE_ROLE = 'SUBMISSION_FILE' and study_id = '{}' AND ANALYSIS_STATUS!='supressed' AND ANALYSIS_STATUS!='cancelled' GROUP BY study_id, project_id, sample_id, analysis_id, data_file_format, tax_id, data_file_role ".format(
            study_id
        )
        return self.data_access_object._runQuery(query)

    def retrieve_webin_accound_id(self, study_id):
        """
        Returns ENA's submission account id
        :param study_id:
        :return:
        """

        query = "SELECT SUBMISSION_ACCOUNT_ID FROM ERA.V_MGP_ASSEMBLY_STUDY WHERE study_id = '{}'".format(
            study_id
        )
        return self.data_access_object._runQuery(query)


class EmgDBError(Exception):
    """EMG Data access exception"""

    def __init__(self, msg):
        self.msg = msg
