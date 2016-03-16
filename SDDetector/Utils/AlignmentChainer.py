#!/usr/bin/env python

import logging

from SDDetector.Entities.Chain import Chain

class AlignmentChainer(object):

    def __init__(self, db, maxGap=3000, logLevel='ERROR'):
        """AlignmentChainer constructor"""

        self.db = db
        self.maxGap = maxGap
        self.logLevel = logLevel
        self.dIndex = {}
        self.lChains = []

        logging.basicConfig(level=logLevel)


    def chainAlignments(self, lAlgmts):
        """Build the list of chains and keep position of alignments in chains"""

        for algmt in lAlgmts:
            if algmt.id not in self.dIndex:
                chain = Chain([algmt])
                index = len(self.lChains)
                self.lChains.append(chain)
                self.dIndex[algmt.id] = [index]

            lChainIdsCurrentAlgmt = self.dIndex[algmt.id]
	    lProximalAlgmts = self.db.selectProximalAlgmts(algmt.id, self.maxGap)
            for proxAlgmt in lProximalAlgmts:
                for chainId in lChainIdsCurrentAlgmt:
                    if self.distanceBetweenQueryAlgmts(algmt,proxAlgmt) < self.maxGap:
                        if proxAlgmt.id not in self.lChains[chainId].getIdListOfAlgmts():
                            self.lChains[chainId].lAlgmts.append(proxAlgmt)
                            if proxAlgmt.id in self.dIndex:
                                self.dIndex[proxAlgmt.id].append(chainId)
                            else:
                               self.dIndex[proxAlgmt.id] = [chainId]
        return


    def distanceBetweenQueryAlgmts(self, algmt1, algmt2):
        """Compute and return the distance between 2 alignments"""

        if algmt1.qstart < algmt2.qstart:
            return algmt2.qstart - algmt1.qend
        else:
            return algmt1.qstart - algmt2.qend


    def sortListOfChains(self, lChains):
        """Return a sorted list of chains"""
 
        ltmp = []
        for i, chain in enumerate(lChains):
            lAlgmts = chain.sortListOfAlgmts()
            for j,algmt in enumerate(lAlgmts):
                if j == 0:
                    sstartMin = algmt.sstart
                    sendMax = algmt.send 
                if algmt.sstart < sstartMin:
                    sstartMin = algmt.sstart
                if algmt.send > sendMax:
                    sendMax = algmt.send
            ltmp.append((algmt.sbjct, sstartMin, sendMax, algmt.query, i))
        ltmp.sort(key=lambda x: ('{0:0>150}'.format(x[0]).lower(), x[1], x[2], '{0:0>150}'.format(x[3])))
        return [ lChains[row[4]] for row in ltmp  ]

