import sqlite3
import pandas as pd
import matplotlib.pyplot as plt


def connect_db(db):
    """
    :param db: the data base name in text
    :return: a connection with database and a cursor
    """
    conn = sqlite3.connect(db)
    c = conn.cursor()
    return conn, c


def retrieve_data(conn, sheet, sym, val, sym1=None, val1=None, sym2=None, val2=None):
    """
    :param conn: a connection pointed to database
    :param sheet: a datasheet that is going to be retrieved
    :param sym: a column name that targets the data scope
    :param val: the column value
    :param sym1: position argument
    :param val1: position argument
    :param sym2: position argument
    :param val2: position argument
    :return: a dataframe
    """

    if sym1 is None and sym2 is None:
        df = pd.read_sql_query("SELECT * "
                               "FROM {} "
                               "WHERE {}={};".format(sheet, sym, val), conn)
        return df

    elif sym1 is not None and sym2 is None:
        df = pd.read_sql_query("SELECT * "
                               "FROM {} "
                               "WHERE {}={} AND {}={}};".format(sheet, sym, val, sym1, val1), conn)
        return df

    elif sym1 is not None and sym2 is not None:
        df = pd.read_sql_query("SELECT * "
                               "FROM {} "
                               "WHERE {}={} AND {}={} AND {}={}};".format(sheet, sym, val, sym1,
                                                                          val1, sym2, val2), conn)
        return df


def plot_data(df, xx, yy):
    """
    :param df: a dataframe retrieved from database
    :param xx: x axis label
    :param yy: y axis label
    :return: plot
    """
    figure = plt.figure()
    ax = figure.add_subplot(111)
    ax.plot(df[xx], df[yy])
    plt.show()
    return


def disconnect_db(conn):
    """
    :param conn: a sqlite database connection
    :return: None
    """
    conn.commit()
    conn.close()


# Excample: DB initializes and manipulates SQLite3 databases.
class DB(object):
    def __init__(self, database='ESP.db', statements=None):
        """Initialize a new or connect to an existing database. Accept setup statements to be executed.
        """

        # the database filename
        self.database = database
        # holds incomplete statements
        if statements is None:
            statements = []
        # indicates if selected data is to be returned or printed
        self.display = False

        self.connect()

        # execute setup satements
        self.execute(statements)

        self.close()

    def connect(self):
        """Connect to the SQLite3 database."""

        self.connection = sqlite3.connect(self.database)
        self.cursor = self.connection.cursor()
        self.connected = True
        self.statement = ''

    def close(self):
        """Close the SQLite3 database."""

        self.connection.commit()
        self.connection.close()
        self.connected = False

    def incomplete(self, statement):
        """Concatenate clauses until a complete statement is made."""

        self.statement += statement
        if self.statement.count(';') > 1:
            print ('An error has occurerd: ' + 'You may only execute one statement at a time.')
            print('For the statement: %s' % self.statement)
            self.statement = ''
        if sqlite3.complete_statement(self.statement):
            # the statement is not incomplete, it's complete
            return False
        else:
            # the statement is incomplete
            return True

    def execute(self, statements):
        """Execute complete SQL statements.

        Incomplete statements are concatenated to self.statement until they
        are complete.

        Selected data is returned as a list of query results. Example:

        for result in db.execute(queries):
            for row in result:
                print row
        """

        queries = []
        close = False
        if not self.connected:
            # open a previously closed connection
            self.connect()
            # mark the connection to be closed once complete
            close = True
        if type(statements) == str:
            # all statements must be in a list
            statements = [statements]
        for statement in statements:
            if self.incomplete(statement):
                # the statement is incomplete
                continue
            # the statement is complete
            try:
                statement = self.statement.strip()
                # reset the test statement
                self.statement = ''
                self.cursor.execute(statement)
                # retrieve selected data
                data = self.cursor.fetchall()
                if statement.upper().startswith('SELECT'):
                    # append query results
                    queries.append(data)

            except sqlite3.Error as error:
                print('An error occurred:', error.args[0])
                print('For the statement:', statement)

        # only close the connection if opened in this function
        if close:
            self.close()
        # print results for all queries
        if self.display:
            for result in queries:
                if result:
                    for row in result:
                        print(row)
                else:
                    print(result)
        # return results for all queries
        else:
            return queries

    def terminal(self):
        """A simple SQLite3 terminal.

        The terminal will concatenate incomplete statements until they are
        complete.
        """

        self.connect()
        self.display = True

        print('SQLite3 terminal for %s. Press enter for commands.' % self.database)

        while True:
            statement = input('')
            if statement == '':
                user = input(
                    'Type discard, exit (commit), or press enter (commit): ')
                if not user:
                    self.connection.commit()
                elif user == 'discard':
                    self.connect()
                elif user == 'exit':
                    break
            self.execute(statement)

        self.display = False
        self.close()


if __name__ == "__main__":
    db = "ESP.db"
    conn, c = connect_db(db)
    sheet = "BHI_P100ESP"
    symbol = "RPM"
    value = 3000
    result = retrieve_data(conn, sheet, symbol, value)
    print(result)

    disconnect_db(conn)

