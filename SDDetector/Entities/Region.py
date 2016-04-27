#!/usr/bin/env python

class Region(object):

    __slots__=['seq','start','end','strand']

    def __init__(self, seq, start, end, strand):
        """Region constructor"""

        self.seq = seq
        self.start = start
        self.end = end
        self.strand = strand

    def __eq__(self, other):
        """Equality on all args"""

        return (self.seq,self.start,self.end,self.strand) == (other.seq,other.start,other.end,other.strand)
      
 
