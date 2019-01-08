__author__ = 'scp'

"""
Created on 21 Nove 2016

@author: Simon Potter
"""


class ERADAO(object):
    """
    Data access object to ENA's ERA schema instances (predominantly ERAPRO)
    """

    def __init__(self, data_access_object):
        """
        Constructor
        """
        self.data_access_object = data_access_object

    def retrieve_submitted_files(self, study_id):
        """
            Returns a list of submitted files.
        :type study_id: str
        :return: list[dict]
        """

        query = "SELECT study_id, sample_id, run_id, library_layout, data_file_path, tax_id, library_strategy, library_source, data_file_role FROM v_mgp_run_file WHERE study_id = '{}' and DATA_FILE_ROLE = 'SUBMISSION_FILE' and (DATA_FILE_FORMAT='FASTA' OR DATA_FILE_FORMAT='FASTQ') and DATA_FILE_PATH not like 'err%'".format(
            study_id)
        return self.data_access_object._runQuery(query)

    def retrieve_generated_files(self, study_id):
        """
            Returns a list of generated files.
        :type study_id: str
        :return:
        """

        query = "SELECT study_id, sample_id, run_id, library_layout, data_file_path, tax_id, library_strategy, library_source, data_file_role FROM v_mgp_run_file WHERE study_id = '{}' and DATA_FILE_ROLE = 'GENERATED_FILE' and DATA_FILE_FORMAT='FASTQ' and DATA_FILE_PATH not like 'err%'".format(
            study_id)
        return self.data_access_object._runQuery(query)

    def retrieve_assembly_metadata(self, study_id):
        """
        Returns assembly metadata
        :param study_id:
        :return:
        """

        query = "SELECT study_id, project_id, sample_id, analysis_id, data_file_format, data_file_path, tax_id, bytes, md5 FROM v_mgp_assembly_file WHERE data_file_format = 'FASTA' and study_id = '{}'".format(
            study_id)
        return self.data_access_object._runQuery(query)

    def retrieve_webin_accound_id(self, study_id):
        """
        Returns ENA's submission account id
        :param study_id:
        :return:
        """

        query = "SELECT SUBMISSION_ACCOUNT_ID FROM ERA.V_MGP_ASSEMBLY_STUDY WHERE study_id = '{}'".format(
            study_id)
        return self.data_access_object._runQuery(query)


class EmgDBError(Exception):
    """EMG Data access exception"""

    def __init__(self, msg):
        self.msg = msg