import cx_Oracle
import logging

__author__ = "maxim"

"""
Created on 07/12/2015

@author: Maxim Scheremetjew
"""


class OracleDataAccessObject:
    """Database access object"""

    def __init__(self, dbConnection):
        """
        Constructor
        """
        self.dbConnection = dbConnection

    def close(self):
        self.dbConnection.close()

    def _runQuery(self, query, args=None):
        """Runs the query"""
        with self.dbConnection as c:
            cursor = c.cursor()
            if args:
                cursor.execute(query, args)
            else:
                cursor.execute(query)
            results = []
            result = cursor.fetchone()
            while result is not None:
                numCols = len(result)
                item = {}
                for col in range(0, numCols):
                    key = cursor.description[col][0]
                    value = result[col]
                    try:
                        isinstance(value, str)
                    except UnicodeEncodeError:
                        value = value.encode("utf-8", "replace")
                    try:
                        # handle oracle clob datatypes
                        value = result[col].read()
                    except AttributeError:
                        pass
                    item[key] = value
                    # item = {cursor.description[col][0]: result[col] for col in range(0, numCols)}
                results.append(item)
                result = cursor.fetchone()
            cursor.close()
            return results

    def _runInsertStatement(self, insert_stmt, data=None):
        logging.error("Should not be running an INSERT to an Oracle DB")
        exit(1)
        try:
            """Runs the given insert statement"""
            logging.debug("Running the following insert statement: " + insert_stmt)
            with self.dbConnection as c:
                cursor = c.cursor()
                # insert_stmt = c.escape_string(insert_stmt)
                if data:
                    cursor.execute(insert_stmt, data)
                else:
                    cursor.execute(insert_stmt)
                c.commit()
                cursor.close()
        except cx_Oracle.DatabaseError as exception:
            (error,) = exception
            logging.error("Oracle error: ", error.message)
        except UnicodeEncodeError as unicodeException:
            error = unicodeException
            logging.error("Encoding error: ", error.message)


if __name__ == "__main__":
    pass
