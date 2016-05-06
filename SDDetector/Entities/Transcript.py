#!/usrbin/env python

class Transcript(object):


    def __init__(self, id, seqid, start, end, strand, gene_id):
        """Transcript constructor"""

        self.id = id
        self.seqid = seqid
        self.start = start
        self.end = end
        self.strand = strand
        self.gene_id = gene_id

    def __eq__(self, other):
        """Equality on all args"""
      
        return ((self.id,self.seqid,self.start,self.end,self.strand,self.gene_id) == (other.id, other.seqid, other.start, other.end, other.strand, other.gene_id))

    def __repr__(self):
        """Transcript representation"""

        return 'Transcript: {}'.format(self.id)



