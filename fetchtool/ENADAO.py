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
            project_id
        )
        return self.data_access_object._runQuery(query)


class EmgDBError(Exception):
    """EMG Data access exception"""

    def __init__(self, msg):
        self.msg = msg
