#!/usr/bin/env python

import logging
import re

from SDDetector.Entities.Duplication import Duplication
from SDDetector.Entities.Region import Region


class GffDuplicationParser(object):

    def __init__(self, inputDuplicationGffFile=""):
        """Constructor"""
        self.inputDuplicationGffFile = inputDuplicationGffFile

    def getAllDuplications(self):
        """Return list of all Duplications, with duplicates = read from alignments"""

        lDuplications = []
        currentDupID = None
        seq1Dup = ()
        lRegions = []
 
        with open(self.inputDuplicationGffFile,  'r') as input :
            for line in input:
                if not re.match('^#',line):
                    values = line.split('\t')
                    if values[2] == 'match':
                        
                        # create and insert current Duplication before deletion
                        if currentDupID != None:
                            duplication = self._buildDuplication(seq1Dup,lRegions)

                            # create Duplication
                            lDuplications.append(duplication)
                            currentDupID = None
                            seq1Dup = ()
                            lRegions = []
 
                        # new match = new Duplication
                        m = re.match(r".*ID=([^;]*);{0,1}.*",values[8])
                        if m:
                            currentDupID = m.group(1) 
                            seq1Dup = (values[0],int(values[3]),int(values[4]))
                        else:
                            raise Exception('Cannot find tag ID for match feature in gff file: {}'.format(self.inputDuplicationGffFile))
               
                    # new match_part 
                    if values[2] == 'match_part':
                        m1 = re.match(r".*Parent=([^;]*);{0,1}.*",values[8])
                        if m1:
                            if currentDupID == m1.group(1):
                                reg1Strand = 1
                                if values[6] == '-':
                                    reg1Strand = -1
                                reg1 = Region(values[0],int(values[3]),int(values[4]),reg1Strand)
                                
                                m2 = re.match(r".*Target=\s*(\w+)\s+(\d+)\s+(\d+)\s*;.*",values[8])
                                if m2:
                                    reg2 = Region(m2.group(1),int(m2.group(2)),int(m2.group(3)),1)
                                    lRegions.append((reg1,reg2))
                                else:
                                    raise Exception('No Target tag for match_part feature in gff file: {} \nline: {}'.format(self.inputDuplicationGffFile,line))
                            else:
                                raise Exception('Error missing match feature \'{}\' for match_part features'.format(m.group(1)))
                        else:
                            raise Exception('Cannot find tag Parent for match_part feature in gff file: {}\nline: {}'.format(self.inputDuplicationGffFile,line))
                
        # add last duplication
        duplication = self._buildDuplication(seq1Dup,lRegions)
        lDuplications.append(duplication)

        input.close()
        return lDuplications

    def _buildDuplication(self,seq1Dup,lRegions):
        """build duplication"""

        # define seq2 boundaries
        seq2 = None
        seq2Start = 999999999999
        seq2End = 0
        for reg1,reg2 in lRegions:
            if reg2.start < seq2Start:
                seq2Start = reg2.start
            if reg2.end > seq2End:
                seq2End = reg2.end
                seq2 = reg2.seq
                            
        seq2Dup = (seq2,int(seq2Start),int(seq2End))

        # create and return Duplication
        return Duplication(seq1Dup[0],seq1Dup[1],seq1Dup[2],seq2Dup[0],seq2Dup[1],seq2Dup[2],lRegions)


    def getNonRedondantDuplications(self):
        """Return list of non-redondant duplication: 1 for 2 sequences"""

        lDuplications = self.getAllDuplications()
        lNonRedDuplications = []
        lindex = []
        for dup in lDuplications:
            if (dup.seq2,dup.start2,dup.end2,dup.seq1,dup.start1,dup.end1) not in lindex:
                lNonRedDuplications.append(dup)
                lindex.append((dup.seq1,dup.start1,dup.end1,dup.seq2,dup.start2,dup.end2))

        return lNonRedDuplications

if __name__ == "__main__":
    gffDuplicationParser = GffDuplicationParser()
