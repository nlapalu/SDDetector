#!/usrbin/env python

class CDS(object):


    def __init__(self, id, cds_id, seqid, start, end, strand, transcript_id):
        """Transcript constructor"""

        self.id = id
        self.cds_id = cds_id
        self.seqid = seqid
        self.start = start
        self.end = end
        self.strand = strand
        self.transcript_id = transcript_id

    def __eq__(self, other):
        """Equality on all args"""
      
        return ((self.id,self.cds_id,self.seqid,self.start,self.end,self.strand,self.transcript_id) == (other.id, other.cds_id,other.seqid, other.start, other.end, other.strand, other.transcript_id))

    def __repr__(self):
        """CDS representation"""

        return 'CDS: {}'.format(self.id)


