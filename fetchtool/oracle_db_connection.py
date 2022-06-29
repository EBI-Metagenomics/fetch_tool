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


import cx_Oracle


class OracleDBConnection:
    """Context manager class for oracle DB connection"""

    def __init__(self, user, password, host, port, instance):
        self.user = user
        self.password = password
        self.host = host
        self.port = port
        self.instance = instance
        self.connection = None

    def __enter__(self):
        if self.connection is not None:
            raise RuntimeError("Connection already exists")
        connStr = None
        if not self.host and not self.port:
            connStr = "%s/%s@%s" % (self.user, self.password, self.instance)
        else:
            connStr = "%s/%s@%s:%s/%s" % (
                self.user,
                self.password,
                self.host,
                self.port,
                self.instance,
            )
        self.connection = cx_Oracle.connect(connStr)
        return self.connection

    def __exit__(self, ext_type, exc_value, traceback):
        self.connection.close()
        self.connection = None


if __name__ == "__main__":
    pass
