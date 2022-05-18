"""
Created on 1 Apr 2015

@author: maq
@author: Maxim Scheremetjew
"""

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
