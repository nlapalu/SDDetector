#!/usr/bin/env python

import xml.parsers.expat
import logging
import sys

from SDDetector.Entities.Alignment import Alignment

class BlastXMLParserExpat(object):

    def __init__(self, inputBlastXMLFile=""):
        """Constructor"""


        if inputBlastXMLFile:
            self.setFile(inputBlastXMLFile)

        self.lAlgmts = []
        self.elmt = None
        self.query = None
        self.sbjct = None
        self.qstart = None
        self.qend = None
        self.sstart = None
        self.send = None
        self.qframe = None
        self.sframe = None
        self.length = None
        self.ident = None
        self.nbAlgmts = 0
        self.lRegionIndex = []
        self.dRegionAlgmts = {}

    def setFile(self,filename):
        """set input xml blast file"""

        try:
            self.inputBlastXMLFile = file(filename)
        except IOError as err:
            logging.error('Error reading input Blast XML File: {}'.format(err))
            sys.exit(1)

    def parseToSDDetector(self):
        """parsing performed to feed SDDetector structures"""

        p = xml.parsers.expat.ParserCreate()
        p.buffer_text = True #required to avoid multi-callback of CharacterDataHandler (see doc)
        p.buffer_size = 1000000000
        p.StartElementHandler = self.startElementSDDetector
        p.CharacterDataHandler = self.charDataSDDetector 
        p.EndElementHandler = self.endElementSDDetector
        p.ParseFile(self.inputBlastXMLFile)

    def startElementSDDetector(self, name, attrs):
        """start"""

        self.elmt = name

    def charDataSDDetector(self, data):
        """data"""

        if self.elmt == 'Iteration_query-def':
            self.qid = data
        elif self.elmt == 'Hit_id':
            self.sid = data
        elif self.elmt == 'Hsp_query-from':
            self.qstart = int(data)
        elif self.elmt == 'Hsp_query-to':
            self.qend = int(data)
        elif self.elmt == 'Hsp_hit-from':
            self.sstart = int(data)
        elif self.elmt == 'Hsp_hit-to':
            self.send = int(data)
        elif self.elmt == 'Hsp_query-frame':
            self.qframe = int(data)
        elif self.elmt == 'Hsp_hit-frame':
            self.sframe = int(data)
        elif self.elmt == 'Hsp_identity':
            self.ident = int(data)
        elif self.elmt == 'Hsp_align-len':
            self.length = int(data)
        elif self.elmt == 'Hsp_qseq':
            self.qseq = data
        elif self.elmt == 'Hsp_hseq':
            self.sseq = data
        else:
            pass
            
        
    def endElementSDDetector(self, name):
        """end"""

        self.elmt = None
        if name == 'Hsp' and not self.lRegionIndex:
            self.addAlignment()
        elif name == 'Hsp' and self.lRegionIndex:
            self.getAlignmentFromRegions()
    
    
    def addAlignment(self):
        """add alignment to the list"""

        self.nbAlgmts += 1

        if self.sframe == 1:
            self.lAlgmts.append(Alignment(self.qid,self.sid,self.qstart,
                                self.qend,self.sstart,self.send,self.length,
                                self.ident, self.qframe, self.sframe, id=self.nbAlgmts))
        elif self.sframe == -1:
            self.lAlgmts.append(Alignment(self.qid,self.sid,self.qstart,
                                self.qend,self.send,self.sstart,self.length,
                                self.ident, self.qframe, self.sframe, id=self.nbAlgmts))
        else:
            logging.error('Blast Parsing: Unknown strand')

    def getAlignmentFromRegions(self):
        """...."""

        index = (self.sid,self.sstart,self.send,self.qid,self.qstart,self.qend)
        if index in self.lRegionIndex:
            self.dRegionAlgmts[index] = (self.sseq,self.qseq)
            idx = self.lRegionIndex.index(index)
            del self.lRegionIndex[idx]
            
            

    def getAllAlignments(self):
        """all alignments"""

        self.parseToSDDetector()
        return self.lAlgmts

    def getAlignmentsFromTupleOfRegions(self, lRegions):

        for reg1,reg2 in lRegions:
            if reg1.strand == -1:
                self.lRegionIndex.append((reg1.seq,reg1.end,reg1.start,reg2.seq,reg2.start,reg2.end))
            else:
                self.lRegionIndex.append((reg1.seq,reg1.start,reg1.end,reg2.seq,reg2.start,reg2.end))

        self.dRegionAlgmts = { i: () for i in self.lRegionIndex}

        lRegionIndexParsed = self.lRegionIndex[:]

        self.parseToSDDetector()

        lRegionAlgmts = [ () for i in range(len(lRegions)) ]
        for i,val in enumerate(lRegionIndexParsed):
            if val in self.dRegionAlgmts:
                lRegionAlgmts[i] = self.dRegionAlgmts[val]
            else:
               logging.error('Error parsing, missing sequence alignment for {}'.format(lRegionIndexParsed[i])) 
        

        for i,a in enumerate(lRegionAlgmts):
            if len(lRegionAlgmts[i]) != 2:
                logging.error('Error missing sequence alignment for {}'.format(lRegionIndexParsed[i])) 

        return lRegionAlgmts



if __name__ == "__main__":
    parser = BlastXMLParserExpat()

    # test
    parser.setFile("/media/backup/ubucluster/save/nlapalu/TAIR10_chr1.xml")
    lalgmt = parser.getAllAlignments()
    print lalgmt[0]
    print len(lalgmt)
