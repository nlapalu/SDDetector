#!/usr/bin/env python

import os
import logging

from SDDetector.Db.SqliteDB import SqliteDB

from SDDetector.Entities.Gene import Gene
from SDDetector.Entities.Transcript import Transcript
from SDDetector.Entities.CDS import CDS

class GeneDB(SqliteDB):

    def __init__(self, dbfile='', logLevel='ERROR'):
        """GeneDB Constructor"""

        SqliteDB.__init__(self, dbfile=dbfile, logLevel=logLevel)
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

#    def _createExonTable(self):
#        """Create table Exon"""

#        self.conn.execute('''create table exon (id text primary key not null,
#                                                seqid text not null,
#                                                start int not null,
#                                                end int not null,
#                                                strand int not null,
#                                                transcript_ids text not null,
#                                                gene_id text);''')
#        self.conn.commit()
  
    
    def _createCDSTable(self):
        """Create table CDS"""

        self.conn.execute('''create table cds (cds_id text not null,
                                               seqid text not null,
                                               start int not null,
                                               end int not null,
                                               strand int not null,
                                               transcript_id text not null);''')
        self.conn.commit()

    def insertlGenes(self, lGenes):
        """Insert a list of genes"""

        lTranscripts = []
        lCDS = []

        for gene in lGenes:
            lTranscripts.extend(gene.lTranscripts)
            for transcript in gene.lTranscripts: 
                lCDS.extend(transcript.lCDS)

        logging.info('Inserting gene features in db')

        self.conn.executemany('''insert into gene (id, seqid, start,end,strand) values (?,?,?,?,?)''', [(gene.id,gene.seqid,gene.start,gene.end,gene.strand) for gene in lGenes ])

        logging.info('Inserting mRNA features in db')

        self.conn.executemany('''insert into transcript (id, seqid, start,end,strand,gene_id) values (?,?,?,?,?,?)''', [(transcript.id,transcript.seqid,transcript.start,transcript.end,transcript.strand, transcript.gene_id) for transcript in lTranscripts ])

        logging.info('Inserting CDS features in db')

        self.conn.executemany('''insert into CDS (cds_id, seqid, start,end,strand,transcript_id) values (?,?,?,?,?,?)''', [(cds.cds_id,cds.seqid,cds.start,cds.end,cds.strand,cds.transcript_id) for cds in lCDS ])

        self.commit()


    def selectAllGenes(self):
        """Select all genes"""


        lGenes = []
        dGenes = {}
        lTranscripts = []
        dTranscripts = {}
        lCDS = []
        dCDS = {}

        cursor = self.conn.execute('''select id, seqid, start, end, strand from gene''') 
        for row in cursor:        
            dGenes[row[0]] = Gene(row[0],row[1],row[2],row[3],row[4])

        cursor = self.conn.execute('''select id, seqid, start,end,strand,gene_id from transcript''')
        for row in cursor:
            transcript = Transcript(row[0],row[1],row[2],row[3],row[4],row[5])
            dTranscripts[row[0]] = transcript

            if len(dGenes[transcript.gene_id].lTranscripts) > 0:
                dGenes[transcript.gene_id].lTranscripts.append(transcript)
            else:
                dGenes[transcript.gene_id].lTranscripts = [transcript]

        cursor = self.conn.execute('''select cds_id, seqid, start,end,strand,transcript_id from cds order by start''')
        for row in cursor:
            cds = CDS(row[0],row[1],row[2],row[3],row[4],row[5])

            if len(dTranscripts[cds.transcript_id].lCDS) > 0:
                dTranscripts[cds.transcript_id].lCDS.append(cds)
            else:
                dTranscripts[cds.transcript_id].lCDS = [cds]

        return dGenes.values()


    def getlGenesFromCoordinates(self, seqid, start, end):
        """Get genes overlapping a defined region"""

        ###########################
        ####### TODO reactoring with selectAllGenes + polymorphirsm
        ### Extract gene included in dup / not overlapping


        lGenes = []
        dGenes = {}
        lTranscripts = []
        dTranscripts = {}
        lCDS = []
        dCDS = {}

#        cursor = self.conn.execute('''select id, seqid, start, end, strand from gene where seqid = \'{}\' and start < {} and end > {}'''.format(seqid,end,start)) 
        cursor = self.conn.execute('''select id, seqid, start, end, strand from gene where seqid = \'{}\' and start > {} and end < {} order by start'''.format(seqid,start,end)) 
        for row in cursor:        
            dGenes[row[0]] = Gene(row[0],row[1],row[2],row[3],row[4])

#        cursor = self.conn.execute('''select id, seqid, start,end,strand,gene_id from transcript where seqid = \'{}\' and start < {} and end > {}'''.format(seqid,end,start))
        if dGenes:

            for i in dGenes:
                print i

            cursor = self.conn.execute('''select id, seqid, start,end,strand,gene_id from transcript where seqid = \'{}\' and start > {} and end < {} order by start'''.format(seqid,start,end)) 
            for row in cursor:
                transcript = Transcript(row[0],row[1],row[2],row[3],row[4],row[5])

                if transcript.gene_id in dGenes:
                    dTranscripts[row[0]] = transcript
                    if len(dGenes[transcript.gene_id].lTranscripts) > 0:
                        dGenes[transcript.gene_id].lTranscripts.append(transcript)
                    else:
                        dGenes[transcript.gene_id].lTranscripts = [transcript]

#        cursor = self.conn.execute('''select cds_id, seqid, start,end,strand,transcript_id from cds where seqid = \'{}\' and start < {} and end > {}'''.format(seqid,end,start)) 
        if dTranscripts:

            for i in dTranscripts:
                print i

            cursor = self.conn.execute('''select cds_id, seqid, start,end,strand,transcript_id from cds where seqid = \'{}\' and start > {} and end < {} order by start'''.format(seqid,start,end)) 
            for row in cursor:
                cds = CDS(row[0],row[1],row[2],row[3],row[4],row[5])

                if cds.transcript_id in dTranscripts:
                    if len(dTranscripts[cds.transcript_id].lCDS) > 0:
                        dTranscripts[cds.transcript_id].lCDS.append(cds)
                    else:
                        dTranscripts[cds.transcript_id].lCDS = [cds]

        return dGenes.values()

         
    def  commit(self):
        """Commit transactions"""

        self.conn.commit()


    def deleteDB(self):
        """Delete the DB"""

        self.conn.close()
        os.remove(self.dbfile)

