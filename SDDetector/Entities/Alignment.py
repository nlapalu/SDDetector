#!/usr/bin/env python

class Alignment(object):

    __slots__ = ['query','sbjct','qstart','qend','sstart','send',
                 'length', 'identities', 'qstrand', 'sstrand','id']

    def __init__(self, query, sbjct, qstart, qend, sstart, send, length, identities, qstrand, sstrand, id=''):
        """Alignment constructor"""

        self.query = query
        self.sbjct = sbjct
        self.qstart = qstart
        self.qend = qend
        self.sstart = sstart
        self.send = send
        self.length = length
        self.identities = identities
        self.qstrand = qstrand
        self.sstrand = sstrand
        self.id = id

    def __eq__(self, other):
        """Equality on all args"""
      
        return (self.query,self.sbjct,self.qstart,self.qend,self.sstart,self.send,self.length,self.identities,self.qstrand,self.sstrand, self.id) == (other.query,other.sbjct,other.qstart,other.qend,other.sstart,other.send,other.length,other.identities,other.qstrand,other.sstrand, other.id)

    def __repr__(self):
        """representation"""

        return ('{}-{}-{}-{}-{}-{}-{}-{}-{}-{}-{}'.format(self.query,self.sbjct,self.qstart,self.qend,self.sstart,self.send,self.length,self.identities,self.qstrand,self.sstrand, self.id))

    def convertToGff3(self):
        """Write alignment in gff3 format"""

        if self.sstrand == 1:
            strand = '+'
        else:
            strand = '-'

        gff3line = '{}\tSQLiteDB\tmatch\t{}\t{}\t.\t{}\t.\tID=match{};Target={} {} {}\n'.format(self.sbjct,self.sstart,self.send,strand,self.id,self.query,self.qstart,self.qend)
        return gff3line

