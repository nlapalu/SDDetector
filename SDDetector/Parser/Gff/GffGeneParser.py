#!/usr/bin/env python

import logging
import re


class GffGeneParser(object):

    def __init__(self, inputGffFile=""):
        """Constructor"""

        self.inputGffFile = inputGffFile
        self.status = 'nonparsed'

    def getAllGenes(self):
        """Get all genes"""
        
        if self.status == 'nonparsed':
            self._parse('gene')
        else:
            return self.lGenes

    def parse()

        self.status = 'parsed'




#    def _parse(self, type='all'):
#       """Parse gff features"""

#        if  
