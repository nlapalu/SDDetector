#!/usr/bin/env python

import os

from SDDetector.DB.SqliteDB import SqliteDB


class GeneDB(SqliteDB):

    def __init__(self, dbfile='', logLevel='ERROR'):
        """GeneDB Constructor"""

        SqliteDB.__init__(self, ddfile=dbfile, logLevel=logLevel)
        self._createDBSchema()

    def _createDBSchema(self):
        """Create Database Schema"""

        self._createGeneTable()
        self._createTranscriptTable()
#        self._createExonTable() incomplete implem 1 exon could be shared by several mRNA
        self._createCDSTable()

    def _createGeneTable(self):
        """Create table Gene"""

        self.conn.execute('''create table gene (id text primary key not null,
                                                seqid text not null,
                                                start int not null,
                                                end int not null,
                                                strand int not null);''')
        self.conn.commit()

    def _createTranscriptTable(self):
        """Create table Transcript"""

        self.conn.execute('''create table transcript (id text primary key not null,
                                                      seqid text not null,
                                                      start int not null,
                                                      end int not null,
                                                      strand int not null,
                                                      gene_id text not null);''')
        self.conn.commit()

    def _createExonTable(self):
        """Create table Exon"""

        self.conn.execute('''create table exon (id text primary key not null,
                                                seqid text not null,
                                                start int not null,
                                                end int not null,
                                                strand int not null,
                                                transcript_ids text not null,
                                                gene_id text);''')
        self.conn.commit()
  
    
    def _createCDSTable(self):
        """Create table CDS"""

        self.conn.execute('''create table cds (id int primary key not null,
                                               cds_id text not null,
                                               seqid text not null,
                                               start int not null,
                                               end int not null,
                                               strand int not null,
                                               transcript_id text not null);''')
        self.conn.commit()

    def insertlGenes(self, lGenes):
        """Insert a list of genes"""

        self.conn.executemany('''insert into gene (id, query, sbjct,
                                 qstart, qend, sstart,send, length, identities,
                                 qstrand, sstrand) values (?,?,?,?,?,?,?,?,?,?,?)''', [(algmt.id, algmt.query,algmt.sbjct, algmt.qstart, algmt.qend, algmt.sstart, algmt.send, algmt.length, algmt.identities, algmt.qstrand, algmt.sstrand) for algmt in lAlignments ])
        self.commit()

     
