#!/usrbin/env python

class Feature(object):

    def __init__(self, id, seqid, start, end, strand, type):
        """Feature constructor"""

        self.id = id
        self.seqid = seqid
        self.start = start
        self.end = end
        self.strand = strand
        self.type = type

    def __eq__(self, other):
        """Equality on all args"""
      
        return ((self.id,self.seqid,self.start,self.end,self.strand, self.type) == (other.id, other.seqid, other.start, other.end, other.strand, other.type))

    def __repr__(self):
        """Feature representation"""

        return 'Feature: {}-{}-{}-{}-{}-{}'.format(self.type,self.id,self.seqid,self.start,self.end,self.strand)
