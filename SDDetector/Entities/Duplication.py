#!/usr/bin/env python

class Duplication(object):

    def __init__(self, seq1, start1, end1, seq2, start2, end2,lRegions=[],lSeqAlgmts=[]):
        """Constructor"""

        self.seq1 = seq1
        self.start1 = start1
        self.end1 = end1
        self.seq2 = seq2
        self.start2 = start2
        self.end2 = end2
        self.lRegions =  lRegions # liste of tulpe of 2 objects Region
        self.lSeqAlgmts = lSeqAlgmts # liste of tuple of 2 string  = sequence aligned

    def __eq__(self, other):
        """Equality on all args"""
      
        return (self.seq1,self.start1,self.end1,self.seq2,self.start2,self.end2,self.lRegions,self.lSeqAlgmts) == (other.seq1,other.start1,other.end1,other.seq2,other.start2,other.end2,other.lRegions,other.lSeqAlgmts)


    def __repr__(self):
        """representation"""
         
        return('{}-{}-{}-{}-{}-{}-{}-{}'.format(self.seq1,self.start1,self.end1,self.seq2,self.start2,self.end2,self.lRegions,self.lSeqAlgmts))
