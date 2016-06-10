#!/usr/bin/env python

class GeneLink(object):

    def __init__(self, dup, gene1, gene2):
        """Constructor"""

        self.dup = dup
        self.gene1 = gene1
        self.gene2 = gene2
        self.shortlinkCoordinates = self.getShortlinkCoordinates()
        self.extendedlinkCoordinates ) self.getExtendedlinkCoordinates()
