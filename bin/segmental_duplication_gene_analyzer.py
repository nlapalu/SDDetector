#!/usr/bin/env python

import argparse
import logging

from SDDetector.Parser.Gff.GffDuplicationParser import GffDuplicationParser
from SDDetector.Parser.Gff.GffGeneParser import GffGeneParser
from SDDetector.Parser.Blast.BlastXMLParser import BlastXMLParser

from SDDetector.Db.GeneDB import GeneDB


class GeneAnalyzer(object):


    def __init__(self, SDFile='',BlastXMLFile='',GeneFile='', logLevel='ERROR'):
        """Constructor"""

        self.SDFile = SDFile
        self.BlastXMLFile = BlastXMLFile
        self.GeneFile = GeneFile
        self.logLevel = logLevel
        logging.basicConfig(level=self.logLevel)

        self.lDuplications = []


    def parseAttributesFromArgsCLI(self):
        """Parse arguments from command line"""

        version = 0.1
        description = "SDAnalyzer: analyzes segmental duplications in genome"
        parser = argparse.ArgumentParser(description=description)
        parser.add_argument("SDFile", help="Segmental Duplication gff3 file = output of SDDetector (filtered or not)", type=str)
        parser.add_argument("BlastXMLFile", help="Input Blast XML file", type=str)
        parser.add_argument("GeneFile", help="gene annotation file in gff3 format", type=str)
        #parser.add_argument("outputFile", help="Output File in gff3 format", type=str)

        parser.add_argument("-v", "--verbosity", type=int, choices=[1,2,3],
                            help="increase output verbosity 1=error, 2=info, 3=debug")

        args = parser.parse_args()
        self._setAttributesFromArgsCLI(args)


    def _setAttributesFromArgsCLI(self, args):
        """Set attributes from argparse"""
 
        if args.verbosity == 1:
            self.logLevel = 'ERROR'
        if args.verbosity == 2:
            self.logLevel = 'INFO'
        if args.verbosity == 3:
            self.logLevel = 'DEBUG'
        logging.getLogger().setLevel(self.logLevel)

        self.SDFile = args.SDFile
        self.BlastXMLFile = args.BlastXMLFile
        self.GeneFile = args.GeneFile


    def runAnalyze(self):
        """run analyze"""

        iGffDuplicationParser = GffDuplicationParser(self.SDFile)
        #self.lDuplications = iGffDuplicationParser.getNonRedondantDuplications()
        lRegions = []
        for dup in self.lDuplications:
            for region in dup.lRegions:
                lRegions.append(region)
        
        iBlastXMLParser = BlastXMLParser(self.BlastXMLFile)
        #lAlignmentTuples = iBlastXMLParser.getAlignmentsFromTupleOfRegions(lRegions)
        lAlignmentTuples = []

        index = 0
        for dup in self.lDuplications:
            lAlgmts = []
            for region in dup.lRegions:
                lAlgmts.append((lAlignmentTuples[index][0],lAlignmentTuples[index][1]))
                index += 1
            dup.lSeqAlgmts = lAlgmts
                
        
        #lDuplicationBoundaries = self.extractDuplicationBoundaries()
        print self.GeneFile
        iGffGeneParser = GffGeneParser(self.GeneFile)
        self.db = GeneDB(dbfile='gene.db')
        self.db.insertlGenes(iGffGeneParser.getAllGenes())
        
        #lGenesInDuplications = iGffGeneParser.getGenesIncludedInBoundaries(lDuplication)


    def extractDuplicationBoundaries(self):
        """extract duplication boundaries"""

        for dup in self.Duplication:
            pass
        

if __name__ == "__main__":

    analyzer = GeneAnalyzer()
    analyzer.parseAttributesFromArgsCLI()
    analyzer.runAnalyze()
