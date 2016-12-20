#!/usr/bin/env python

import argparse
import logging
import os

from SDDetector.version import __version__
from SDDetector.Db.AlignDB import AlignDB
from SDDetector.Utils.AlignmentChainer import AlignmentChainer
from SDDetector.Parser.Blast.BlastTabParser import BlastTabParser
from SDDetector.Parser.Blast.BlastXMLParserExpat import BlastXMLParserExpat


class Detector(object):

    def __init__(self, db='', inputFile='', inputFormat='', outputFile='',logLevel='ERROR'):
        """Constructor"""

        self.dbFile = db
        self.inputFile = inputFile
        self.inputFormat = inputFormat
        self.outputFile = outputFile
        self.db = None
        self.parser = None
        self.logLevel = logLevel
        logging.basicConfig(level=self.logLevel)

    def __del__(self):
        """Destructor"""

        if (self.db and self.dbFile != ':memory:'):
            self.db.deleteDB()

    def parseAttributesFromArgsCLI(self):
        """Parse arguments from command line"""

        program = 'SDDetector'
        version = __version__
        description = "SDDetector: detects segmental duplications in genome"

        parser = argparse.ArgumentParser(prog=program)
        parser = argparse.ArgumentParser(description=description)
        parser.add_argument('--version', action='version', version='{} {}'.format(program,version))

        parser.add_argument("inputFile", help="Input Blast file", type=str)
        parser.add_argument("inputFormat", help="Input Blast format [xml|tab]", type=str)
        parser.add_argument("outputFile", help="Output File in gff3 format", type=str)
        parser.add_argument("db", help="file to store SQLite db, ex: mydb.db, or ':memory:' \
                                        to keep db in RAM", type=str)

        parser.add_argument("-v", "--verbosity", type=int, choices=[1,2,3],
                            help="increase output verbosity 1=error, 2=info, 3=debug")
        parser.add_argument("-i", "--minIdent", type=float, default=0.9,
                            help="minimum alignment percentage identity, [default=0.9], range [0-1]")
        parser.add_argument("-g", "--maxGap", type=int, default=3000,
                            help="maximum gap size allowed to chain two fragments, in bp \
                            [default=3000]")
        parser.add_argument("-l", "--chainLength", type=int, default=5000,
                            help="minimum chain length (without gaps) to keep, in bp \
                            [default=5000]")
        parser.add_argument("-t", "--matchLength", type=int, default=0,
                            help="minimum match length to consider [default=0]")
        parser.add_argument("-p","--matchOverlap", type=int, default=0,
                            help="maximum overlap to consider between 2 alignments \
                            before removing suboptimal alignment [default=0]")
        parser.add_argument("-r", "--keepOverDup", action="store_true", default=False,
                            help="keep overlapping duplications, [default=False] only the longest is kept")
        parser.add_argument("-s", "--keepInternSimDup", action="store_true", default=False,
                            help="keep duplication with internal similarity, [default=False]")
        parser.add_argument("-a", "--exportall", action="store_true",
                            help="export the content of the db at each step: all, removing \
                            self-matches, identity threshold, suboptimal matches")
        parser.add_argument("-b", "--bed", action="store_true",
                            help="export also results in bed format")

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

        self.inputFormat = args.inputFormat.lower()

        if self.inputFormat not in ['xml','tab']:
            raise Exception('untractable blast format')
            exit(1)

        self.inputFile= args.inputFile
        self.outputFile = args.outputFile

        self.dbFile = args.db
        if self.dbFile == ':memory:':
            logging.debug("SQLite db stored in memory")
        else:
            logging.debug("SQLite db stored in {}".format(self.dbFile))

        if (args.minIdent < 0 or args.minIdent > 1):
            raise Exception('minimum identity not in range [0-1], example: 0.95')
        else:
            self.minIdentity = args.minIdent
        self.maxGap = args.maxGap
        self.chainLength = args.chainLength
        self.matchLength = args.matchLength
        self.matchOverlap = args.matchOverlap
        self.keepOverDup = args.keepOverDup
        self.keepInternSimDup = args.keepInternSimDup
        self.exportDBAllSteps = args.exportall
        self.exportBed = args.bed

    def runSDDetection(self):
        """run segmental duplication detection"""

        self.db = AlignDB(dbfile=self.dbFile)

        logging.info('Parsing alignments')
        self.parseAlignments()
        logging.info('Loading alignments into database')
        if self.matchLength > 0:
            logging.info('Minimum size to consider for alignment: {}pb'.format(self.matchLength))
        self.loadAlignmentsInDb()
        if self.exportDBAllSteps:
            logging.info('Exporting matches after loading in database in gff3 format, file: {}.loading'.format(self.outputFile))
            self.exportMatches('{}.loading'.format(self.outputFile))
        logging.info('Removing self-self matches')
        self.detectAndRemoveSelfMatchAlignments()
        if self.exportDBAllSteps:
            logging.info('Exporting matches after removing self-matches in gff3 format, file: {}.selfmatch'.format(self.outputFile))
            self.exportMatches('{}.selfmatch'.format(self.outputFile))
        logging.info('Removing Alignment below the identity threshold: {}'.format(self.minIdentity))
        self.detectAndRemoveAlignmentBelowIdentityThreshold()
        if self.exportDBAllSteps:
            logging.info('Exporting matches after removing mimimal identity in gff3 format, file: {}.minidentity'.format(self.outputFile))
            self.exportMatches('{}.minidentity'.format(self.outputFile))
        logging.info('Removing suboptimal matches with maximum overlap {}'.format(self.matchOverlap))
        self.detectAndRemoveSuboptimalAlignments()
        if self.exportDBAllSteps:
            logging.info('Exporting matches after removing suboptimal alignments in gff3 format, file: {}.suboptimal'.format(self.outputFile))
            self.exportMatches('{}.suboptimal'.format(self.outputFile))
        logging.info('Chaining alignments with parameters: maximum Gap between fargments = {} bp, minimum chain length = {} bp'.format(self.maxGap, self.chainLength))
        self.chainAlignments(maxGap=self.maxGap, chainLength=self.chainLength)
        if not self.keepOverDup:
            logging.info('Removing overlapping duplications: only the longest one is keep')
            self.removeOverlappingDuplications()
        if not self.keepInternSimDup:
            logging.info('Removing intra-sequence duplications with internal similarity')
            self.removeDuplicationWithInternalSimilarity()
        logging.info('Exporting chains in gff3 format, file: {}'.format(self.outputFile))
        self.exportChains(self.outputFile)
        if self.exportBed:
            logging.info('Exporting chains in bed format, file: {}.bed'.format(self.outputFile))
            self.exportChains('{}.bed'.format(self.outputFile), format='bed')


    def parseAlignments(self):
        """Parse Alignments"""

#        if self.inputFormat == 'xml' and Biopython_available == True:
#            #self.parser = BlastXMLParser(self.inputFile)
#        elif self.inputFormat == 'xml' and Biopython_available == False:
#            raise Exception('BioPython is not installed, xml parsing not available. \
#                             Please choose tab format, or install BioPython')
        if self.inputFormat == 'xml':
            self.parser = BlastXMLParserExpat(self.inputFile)
        elif self.inputFormat == 'tab':
            self.parser = BlastTabParser(self.inputFile)
        else:
            raise Exception('untractable blast format')
            exit(1)


    def loadAlignmentsInDb(self):
        """Load Blast Alignments in Database"""

        self.db.insertlAlignments(self.parser.getAllAlignments(), self.matchLength)


    def detectAndRemoveSelfMatchAlignments(self):
        """Detect and remove identical alignment

           Detect self-matches, intra-chromosomal
           hits with sbjct coordinates identical to
           query coordinates

        """

        lAlgmts = self.db.selectSelfSelfMatchAlgmts()
        self.db.deletelAlignments(lAlgmts)
        self.db.commit()


    def detectAndRemoveAlignmentBelowIdentityThreshold(self, threshold=0.90):
        """Detect and remove alignment below defined identity threshold

           Threshold set to 0.9 by default

        """

        lAlgmts = self.db.selectAlgmtsBelowIdentityThreshold()
        self.db.deletelAlignments(lAlgmts)
        self.db.commit()

    def detectAndRemoveSuboptimalAlignments(self):
        """Detect and remove suboptimal alignment"""

        lAlgmts = self.db.selectSuboptimalAlgmts(self.matchOverlap)

        self.db.deletelAlignments(lAlgmts)
        self.db.commit()

    def chainAlignments(self, maxGap=3000, chainLength=5000):
        """Chain Alignments

           - sort by sbjct, query ann coordiantes
           - chaine if distance < 3000 nt
           - remove if chained alignment < 5000 bp
           - ??? distance on query ???

        """

        lSbjcts = self.db.selectAllSbjcts()
        lQueries = self.db.selectAllQueries()
        lSelectedChains = []
        for sbjct in lSbjcts:
            for query in lQueries:
                logging.debug('Chaining Alignment with subject: {} and query: {}'.format(sbjct, query))
                lAlgmts = self.db.selectAlignmentsWithDefinedSbjctAndQueryOrderBySbjctCoord(sbjct,query)
                chainer = AlignmentChainer(self.db, maxGap=maxGap)
                chainer.chainAlignments(lAlgmts)

                for chain in chainer.lChains:
                    if chain.getLength() > chainLength:
                        lSelectedChains.append(chain)

        chainer2 = AlignmentChainer(self.db, maxGap=maxGap)
        self.lSortedChains = chainer2.sortListOfChains(lSelectedChains)

    def removeOverlappingDuplications(self):
        """Remove Overlapping Duplications"""

        chainer = AlignmentChainer(self.db)
        self.lSortedChains = chainer.removeOverlappingChains(self.lSortedChains)

    def removeDuplicationWithInternalSimilarity(self):
        """Remove Duplication with internal similarity"""

        chainer = AlignmentChainer(self.db)
        self.lSortedChains = chainer.removeChainsWithInternalSimilarity(self.lSortedChains)

    def exportChains(self, outputFile, format='gff3'):
        """Export Chains in gff3|bed format"""

        if format not in ['gff3','bed']:
            raise Exception('format {} is not supported for export'.format(format))

        with open(outputFile, 'w') as f:
            for id, chain in enumerate(self.lSortedChains):
                f.write(chain.convertChain(id+1, format))
        f.close()


    def exportMatches(self, outputFile, format='gff3'):
        """Export Matches in gff3|bed format"""

        if format not in ['gff3','bed']:
            raise Exception('format {} is not supported for export'.format(format))

        with open(outputFile, 'w') as f:
            for algmt in self.db.selectAllAlignments():
                f.write(algmt.convertToGff3())
        f.close()


if __name__ == "__main__":

    detector = Detector()
    detector.parseAttributesFromArgsCLI()
    detector.runSDDetection()
