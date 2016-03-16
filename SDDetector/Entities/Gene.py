#!/usrbin/env python

class Gene(object):


    def __init__(self, id, seqid, start, end, strand, lTranscripts=[]):
        """Gene constructor"""

        self.id = id
        self.seqid = seqid
        self.start = start
        self.end = end
        self.strand = strand
        self.lTranscripts = lTranscripts

    def __eq__(self, other):
        """Equality on all args"""
      
        return ((self.id,self.seqid,self.start,self.end,self.strand, self.lTranscripts) == (other.id, other.seqid, other.start, other.end, other.strand, other.lTranscripts))

    def __repr__(self):
        """Gene representation"""

        return 'Gene: {}-{}-{}-{}-{}-{}'.format(self.id,self.seqid,self.start,self.end,self.strand, self.lTranscripts)



