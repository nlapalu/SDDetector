#!/usr/bin/env python

import os
import logging
import re

class FastaFileIndexer(object):

    def __init__(self, filename, logLevel='ERROR'):
        """FastaIndexer constructor"""

        self.logLevel = logLevel

        logging.basicConfig(level=logLevel)

        self.filename = filename
        self.lSeq = []
        self.dSeq = {}
        self.dIndexSeq = {}
        

    def read(self):
        """read"""

        mySeq = []
        currentSeq = ''
        with open(self.filename, 'r') as f:
            for line in f:
                if line == '\n':
                    next
                m = re.match('>(.*)',line)
                if m:
               
                    if currentSeq:
                        self.dSeq[currentSeq] += ''.join(mySeq)
                        mySeq = []

                    currentSeq = re.split('[\s\|]+',m.group(1))[0]
                    self.lSeq.append(currentSeq)
                    self.dSeq[currentSeq] = ''
                else:
                    mySeq.append(line.replace('\n',''))
            self.dSeq[currentSeq] += ''.join(mySeq)

        
    def getMaxSize(self):
        """Return the size of the longest sequence"""

        maxSize = 0
        for seq in self.dSeq:
            if len(self.dSeq[seq]) > maxSize:
                maxSize = len(self.dSeq[seq])
        return maxSize

    def getSeq(seq,start,stop,strand=1):
        """get a sequence, if strand -, return the reverse complement"""

        seq = self.fp[self.dSeq][start:stop]
        if strand == -1:
            pass

        return seq

    def getlSeq():
        """get a list of""" 
        pass

    def getSeqbyRegions(self):
        """get list of regions from a list of Regions"""
        pass

    def delete(self):
        """delete index"""

        logging.info('Deleting index: {}.idx'.format(self.filename)) 
        os.remove('{}.idx'.format(self.filename))
