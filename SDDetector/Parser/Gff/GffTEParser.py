#!/usr/bin/env python

import logging
import re

from SDDetector.Entities.Feature import Feature 
#from SDDetector.Entities.Region import Region


class GffTEParser(object):

    def __init__(self, inputTEGffFile=""):
        """Constructor"""

        self.inputTEGffFile = inputTEGffFile

    def getAllTEs(self):
        """Return list of all TEs"""

        lTEs = []
 
        with open(self.inputTEGffFile,  'r') as input :
            for line in input:
                if not re.match('^#',line):
                    values = line.split('\t')
                    # new match_part 
                    if values[2] == 'match_part':
                        m = re.match(r".*ID=([^;]*);.*",values[8])
                        if m:
                            strand = 1
                            if values[6] == '-':
                                strand = -1
                            lTEs.append(Feature(m.group(1),values[0] ,int(values[3]),int(values[4]),strand,'TE'))
                        else:
                            raise Exception('Cannot find tag ID for match_part feature in gff file: {}'.format(self.inputTEGffFile))

        input.close()
        return lTEs


if __name__ == "__main__":
    gffTEParser = GffTEParser()
