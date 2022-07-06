#!/usr/bin/env python

"""
   segmental_duplication_detector
"""

import argparse
import logging
import os
import sys
import multiprocessing
import sqlite3

from SDDetector.version import __version__
from SDDetector.Db.AlignDB import AlignDB
from SDDetector.Utils.AlignmentChainer import AlignmentChainer
from SDDetector.Parser.Blast.BlastTabParser import BlastTabParser
from SDDetector.Parser.Blast.BlastXMLParserExpat import BlastXMLParserExpat




def working_process(func, args):
    result = func(*args)
    return result

def chainSbjctQueryAlgmts(db, maxGap, chainLength, sbjct, query, lAlgmts):

    lSelectedChains = []
    logging.debug('Chaining Alignment with subject: {} and query: {}'.format(sbjct, query))
    chainer = AlignmentChainer(db, maxGap=maxGap)
    chainer.chainAlignments(lAlgmts, multiproc=True)
    nb_selected_chains = 0
    for chain in chainer.lChains:
        if chain.getLength() > chainLength:
            lSelectedChains.append(chain)
            nb_selected_chains += 1
    logging.info('Selecting {} chains with subject: {} and query: {}'.format(nb_selected_chains , sbjct, query))

    return lSelectedChains



class Detector(object):


    def __init__(self, db='', inputFile='', inputFormat='', outputFile='',
                 minIdentity=0.9, maxGap=3000, chainLength=5000, matchLength=0,
                 matchOverlap=0, keepOverDup=False, keepInternSimDup=False,
                 exportDBAllSteps=False, exportBed=False, procs=1, logLevel='ERROR'):
        """Constructor"""

        self.dbFile = db
        self.inputFile = inputFile
        self.inputFormat = inputFormat
        self.outputFile = outputFile
        self.minIdentity = minIdentity
        self.maxGap = maxGap
        self.chainLength = chainLength
        self.matchLength = matchLength
        self.matchOverlap = matchOverlap
        self.keepOverDup = keepOverDup
        self.keepInternSimDup = keepInternSimDup
        self.exportDBAllSteps = exportDBAllSteps
        self.exportBed = exportBed
        self.db = None
        self.parser = None
        self.procs = procs
        self.logLevel = logLevel
        logging.basicConfig(level=self.logLevel)


    def __del__(self):
        """Destructor"""

        if (self.db and self.dbFile != ':memory:'):
            self.db.deleteDB()


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
        logging.info('Removing self-matches')
        self.detectAndRemoveSelfMatchAlignments()
        if self.exportDBAllSteps:
            logging.info('Exporting matches after removing self-matches in gff3 format, file: {}.selfmatch'.format(self.outputFile))
            self.exportMatches('{}.selfmatch'.format(self.outputFile))
        logging.info('Removing alignments below the identity threshold: {}'.format(self.minIdentity))
        self.detectAndRemoveAlignmentBelowIdentityThreshold()
        if self.exportDBAllSteps:
            logging.info('Exporting matches after removing mimimal identity in gff3 format, file: {}.minidentity'.format(self.outputFile))
            self.exportMatches('{}.minidentity'.format(self.outputFile))
        logging.info('Removing suboptimal matches with maximum overlap {}'.format(self.matchOverlap))
        self.detectAndRemoveSuboptimalAlignments()
        if self.exportDBAllSteps:
            logging.info('Exporting matches after removing suboptimal alignments in gff3 format, file: {}.suboptimal'.format(self.outputFile))
            self.exportMatches('{}.suboptimal'.format(self.outputFile))
        logging.info('Chaining alignments with parameters: maximum Gap between fragments = {} bp, minimum chain length = {} bp'.format(self.maxGap, self.chainLength))
        self.chainAlignments(maxGap=self.maxGap, chainLength=self.chainLength)
        if not self.keepOverDup:
            logging.info('Removing overlapping duplications: only the longest one is kept')
            self.removeOverlappingDuplications()
        if not self.keepInternSimDup:
            logging.info('Removing intra-sequence duplications with internal similarity')
            self.removeDuplicationWithInternalSimilarity()

        logging.info('Pairing chains')
        self.pairingChains()

        logging.info('Exporting chains in gff3 format, file: {}'.format(self.outputFile))
        self.exportChains(self.outputFile)
        if self.exportBed:
            logging.info('Exporting chains in bed format, file: {}.bed'.format(self.outputFile))
            self.exportChains('{}.bed'.format(self.outputFile), format='bed')


    def parseAlignments(self):
        """Parse Alignments"""

        if self.inputFormat == 'xml':
            self.parser = BlastXMLParserExpat(self.inputFile)
        elif self.inputFormat == 'tab':
            self.parser = BlastTabParser(self.inputFile)
        else:
            raise Exception('untractable blast format')
            exit(1)


    def loadAlignmentsInDb(self):
        """Load Blast Alignments in Database"""

        inserted_sbjcts = set()
        filtered_algmts = []
        u = self.parser.getAllAlignments()
        for algmt in sorted(u, key=lambda x: (x.sbjct ,x.sstart)):
            inserted_sbjcts.add(algmt.sbjct)
            if algmt.sbjct == algmt.query:
                if algmt.sstart < algmt.qstart:
                    filtered_algmts.append(algmt)
            elif algmt.sbjct != algmt.query:
                if algmt.query not in inserted_sbjcts:
                    filtered_algmts.append(algmt)
        self.db.insertlAlignments(filtered_algmts, self.matchLength)


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


    def _chainSbjctQueryAlgmts(self, sbjct, query):

        lSelectedChains = []
        logging.debug('Chaining Alignment with subject: {} and query: {}'.format(sbjct, query))
        lAlgmts = self.db.selectAlignmentsWithDefinedSbjctAndQueryOrderBySbjctCoord(sbjct,query)
        chainer = AlignmentChainer(self.db, maxGap=maxGap)
        chainer.chainAlignments(lAlgmts)
        nb_selected_chains = 0
        for chain in chainer.lChains:
            if chain.getLength() > chainLength:
                lSelectedChains.append(chain)
                nb_selected_chains += 1
        logging.info('Selecting {} chains with subject: {} and query: {}'.format(nb_selected_chains , sbjct, query))

        return lSelectedChains

    def cpSqliteDb(self):

        print self.db.dbfile


    def chainAlignments(self, maxGap=3000, chainLength=5000):
        """Chain Alignments"""

        lSbjcts = sorted(self.db.selectAllSbjcts())
        lQueries = sorted(self.db.selectAllQueries())

        if lSbjcts != lQueries:
            #raise 'error in list of subjects / queries of algmts'
            logging.info("Some sequences without any alignments")

        lSelectedChains = []

        if self.procs > 1:

            lDbs = []
            with open('dump.sql', 'w') as f:
                for line in self.db.conn.iterdump():
                    f.write('%s\n' % line)
            f.close()
            for  i in range(0,self.procs):
                dest = sqlite3.connect('sddetector.{}.db'.format(i))
                cur = dest.cursor()
                f = open('dump.sql','r')
                sql = f.read()
                f.close()
                cur.executescript(sql)
                db = AlignDB('sddetector.{}.db'.format(i), copy=True)
                lDbs.append(db)

            TASKS = []
            idx = 0
            for i, sbjct in enumerate(lSbjcts):
                for j, query in enumerate(lQueries[i:]):
                    lAlgmts = self.db.selectAlignmentsWithDefinedSbjctAndQueryOrderBySbjctCoord(sbjct,query)
                    TASKS.append((chainSbjctQueryAlgmts,(lDbs[idx],maxGap,chainLength,sbjct, query, lAlgmts)))
                    if idx == len(lDbs) -1 :
                        pool = multiprocessing.Pool(self.procs)
                        results = [pool.apply_async(working_process,t) for t in TASKS]

                        for i,r in enumerate(results):
                            r_chains = r.get()
                            lSelectedChains.extend(r_chains)
                        TASKS = []
                        idx = -1
                    idx += 1

            if TASKS:
                pool = multiprocessing.Pool(self.procs)
                results = [pool.apply_async(working_process,t) for t in TASKS]
                for i,r in enumerate(results):
                    r_chains = r.get()
                    lSelectedChains.extend(r_chains)

            for db in lDbs:
                os.remove(db.dbfile)
            os.remove('dump.sql')
        else:

            for i, sbjct in enumerate(lSbjcts):
                for j, query in enumerate(lQueries[i:]):
                    logging.debug('Chaining Alignment with subject: {} and query: {}'.format(sbjct, query))
                    lAlgmts = self.db.selectAlignmentsWithDefinedSbjctAndQueryOrderBySbjctCoord(sbjct,query)
                    chainer = AlignmentChainer(self.db, maxGap=maxGap)
                    chainer.chainAlignments(lAlgmts)
                    nb_selected_chains = 0
                    for chain in chainer.lChains:
                        if chain.getLength() > chainLength:
                            lSelectedChains.append(chain)
                            nb_selected_chains += 1
                    logging.info('Selecting {} chains with subject: {} and query: {}'.format(nb_selected_chains , sbjct, query))

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


    def pairingChains(self):

        chainer = AlignmentChainer(self.db)
        self.lSortedChains = chainer.pairingChains(self.lSortedChains)

    def exportChains(self, outputFile, format='gff3'):
        """Export Chains in gff3|bed format"""

        if format not in ['gff3','bed']:
            raise Exception('format {} is not supported for export'.format(format))

        with open(outputFile, 'w') as f:
            for id, chain in enumerate(self.lSortedChains):
                #f.write(chain.convertChain(id+1, format))
                f.write(chain.convertChain(chain.id, format))
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
                        help="keep duplication with internal similarity, Usefull for \
                        closed tandem duplication connected in one duplication [default=False]")
    parser.add_argument("-a", "--exportall", action="store_true", default=False,
                        help="export the content of the db at each step: all, removing \
                        self-matches, identity threshold, suboptimal matches")
    parser.add_argument("-b", "--bed", action="store_true", default=False,
                        help="export also results in bed format")
    parser.add_argument("--procs", type=int, default=1,
                        help="number of processors [default=1]")
    args = parser.parse_args()

    logLevel = 'ERROR'
    if args.verbosity == 1:
        logLevel = 'ERROR'
    if args.verbosity == 2:
        logLevel = 'INFO'
    if args.verbosity == 3:
        logLevel = 'DEBUG'
    logging.getLogger().setLevel(logLevel)

    inputFormat = args.inputFormat.lower()

    if inputFormat not in ['xml','tab']:
        raise Exception('untractable blast format')
        sys.exit(1)

    if args.db == ':memory:':
        logging.info("SQLite db stored in memory")
    else:
        logging.info("SQLite db stored in {}".format(args.db))

    if (args.minIdent < 0 or args.minIdent > 1):
        raise Exception('minimum identity not in range [0-1], example: 0.95')

    if args.procs > 1:
        args.procs = min(args.procs, multiprocessing.cpu_count()-1)
        logging.info("Multi-processs requested: {} procs will be used".format(args.procs))

    detector = Detector(args.db, args.inputFile, inputFormat, args.outputFile,
                        minIdentity=args.minIdent, maxGap=args.maxGap,
                        chainLength=args.chainLength, matchLength=args.matchLength,
                        matchOverlap=args.matchOverlap, keepOverDup=args.keepOverDup,
                        keepInternSimDup=args.keepInternSimDup,
                        exportDBAllSteps=args.exportall, exportBed=args.bed, procs=args.procs,
                        logLevel=logLevel)
    detector.runSDDetection()
