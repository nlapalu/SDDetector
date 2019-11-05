#!/usr/bin/env python

import logging
import sqlite3

class SqliteDB(object):

    def __init__(self, dbfile='', logLevel='DEBUG', copy=False):
        
        self.dbfile = dbfile

        if logLevel:
            self.setLoggingLevel(logLevel)

        if dbfile:
            if copy == False:
                self.conn = sqlite3.connect(dbfile)
        else:
            dbfile='/tmp/tmpSqlite.db' 
            self.conn = sqlite3.connect(dbfile)
            logging.info("SQLiteDB in {}".format(dbfile))

    def setLoggingLevel(self, logLevel):
        logging.basicConfig(level=logLevel)


    def getListOfTables(self):
        pass
