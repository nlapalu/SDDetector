#!/usr/bin/env python

from SDDetector.Entities.Region import Region

class GeneLink(object):

    def __init__(self, dup, gene1, gene2):
        """Constructor"""

        self.dup = dup
        self.gene1 = gene1
        self.gene2 = gene2
#        self.shortlinkCoordinates = self.getShortlinkCoordinates()
#        self.extendedlinkCoordinates = self.getExtendedlinkCoordinates()

    def getGeneShortlinkCoordinates(self):
        """Return the boundaries of alignment with no gene extra-base"""

        region1 = self.dup.getSeqAlignment(self.gene1.seqid,self.gene1.start,self.gene1.end) 
        region2 = self.dup.getSeqAlignment(self.gene2.seqid,self.gene2.start,self.gene2.end)

        return (region1,region2)


    def getGeneExtendedlinkCoordinates(self):
        """"
            Return the boundaries of the alignment with all gene bases. Overlapping with other 
            is possible if a the largest gene recover at least 2 genes
        """

        pass

    def getCDSAlignment(self):
        """ pass """

        start = None
        end = None
        strand = 1
        startReverse = self.dup.dSeqToSeq[self.gene1.seqid][self.gene1.start][1]
        endReverse = self.dup.dSeqToSeq[self.gene2.seqid][self.gene2.start][1]

        if startReverse > endReverse:
            strand == -1

        if strand == -1:
            if startReverse > self.gene2.end:
                start = self.gene1.start
            else:
                start = self.gene1.start - (self.gene2.end - startReverse)
            if endReverse > self.gene2.start:
                end = self.gene1.start + (self.gene2.start - endReverse)
            else:
                end = self.gene1.start
        
        else:
            if startReverse > self.gene2.start:
                start = self.gene1.start - (startReverse-self.gene2.start)
            else:
                start = self.gene1.start
            if endReverse > self.gene2.end:
                end = self.gene1.end
            else:
                end = self.gene1.end + (self.gene2.end - endReverse)

        print start
        print end 


        for i in range(start,end):
            algmtGene1 = self.dup.getSeqAlignment(self.gene1.seqid,start,end)
            algmtGene2 = self.dup.getSeqAlignment(self.gene2.seqid,self.dup.dSeqToSeq[self.gene1.seqid][start][1],self.dup.dSeqToSeq[self.gene1.seqid][end][1])
        #for CDS in gene1.lTranscripts[0].lCDS:
         #   pass

        #for CDS in gene2

        return (algmtGene1,algmtGene2)
        # il faut retourner un truc du genre

        # gtggtacATGTGATAAtgtgctgtg
        # ATGGTATAT.

        # la modif est du au 1 A

        # il faut partir ATG et travaill par 3 NTs
        # )
