__author__ = 'HD'

"""
Created on 31 Jan 2018

@author: Hubert DENISE
"""


class ENADAO(object):
    """
    Data access object to ENA's ENA schema instances (predominantly ENAPRO)
    """

    def __init__(self, data_access_object):
        """
        Constructor
        """
        self.data_access_object = data_access_object

    def retrieve_assembly_data(self, project_id):
        """
        Returns assembly data
        :param project_id:
        :return:
        """

        query = "SELECT project_acc, sample_id, assembly_id, wgs_acc, gc_id, contig_cnt from gcs_assembly where project_acc = '{}' AND CANCELLED IS NULL".format(
            project_id)
        return self.data_access_object._runQuery(query)


class EmgDBError(Exception):
    """EMG Data access exception"""

    def __init__(self, msg):
        self.msg = msg
