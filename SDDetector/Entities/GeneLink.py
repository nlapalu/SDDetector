#!/usr/bin/env python

from SDDetector.Entities.Region import Region

class GeneLink(object):

    def __init__(self, dup, gene1, gene2):
        """Constructor"""

        self.dup = dup
        self.gene1 = gene1
        self.gene2 = gene2

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
        """ Return alignment between CDS """

        start = None
        end = None
        strand = 1
        startReverse = None
        endReverse = None
        startDeletion = ''
        endDeletion = ''

        startReverse = self.dup.dSeqToSeq[self.gene1.seqid][self.gene1.start][1]
        endReverse = self.dup.dSeqToSeq[self.gene1.seqid][self.gene1.end][1]

        while (startReverse == None):
            startReverse = self.dup.dSeqToSeq[self.gene1.seqid][self.gene1.start+1][1]
            startDeletion +='-'
        while (endReverse == None):
            endReverse = self.dup.dSeqToSeq[self.gene1.seqid][self.gene1.end-1][1]
            endDeletion +='-'

        if startReverse > endReverse:
            strand = -1
            
        if strand == -1:
            if startReverse > self.gene2.end:
                start = self.gene1.start
            else:
                start =  self.dup.dSeqToSeq[self.gene2.seqid][self.gene2.end][1]
            if endReverse > self.gene2.start:
                end = self.dup.dSeqToSeq[self.gene2.seqid][self.gene2.start][1]
            else:
                end = self.gene1.end
        else:
            if startReverse > self.gene2.start:
                start = self.dup.dSeqToSeq[self.gene2.seqid][self.gene2.start][1]
            else:
                start = self.gene1.start
            if endReverse > self.gene2.end:
                end = self.gene1.end
            else:
                end = self.dup.dSeqToSeq[self.gene2.seqid][self.gene2.end][1]

        algmtGene1Start = start 
        algmtGene1End = end

        algmtGene2Start = min(self.dup.dSeqToSeq[self.gene1.seqid][start][1],self.dup.dSeqToSeq[self.gene1.seqid][end][1])
        algmtGene2End = max(self.dup.dSeqToSeq[self.gene1.seqid][start][1],self.dup.dSeqToSeq[self.gene1.seqid][end][1])
        #FIX-BUG - None pos if deletion.
        #print "algmtGene2Start {}:  val: {} - {}".format(algmtGene2Start,self.dup.dSeqToSeq[self.gene1.seqid][start][1],self.dup.dSeqToSeq[self.gene1.seqid][end][1])
 
        #print "algmtGene2End {}:  val: {} - {}".format(algmtGene2End,self.dup.dSeqToSeq[self.gene1.seqid][start][1],self.dup.dSeqToSeq[self.gene1.seqid][end][1])
        

        algmtGene1, algmtGene1Strand = self.dup.getSeqAlignment(self.gene1.seqid,algmtGene1Start,algmtGene1End)
        algmtGene2, algmtGene2Strand = self.dup.getSeqAlignment(self.gene2.seqid,algmtGene2Start,algmtGene2End)

        
        if not algmtGene1 or not algmtGene2:
            return (None,None,None,None)

        if startDeletion:
            algmtGene2 = startDeletion + algmtGene2
        if endDeletion:
            algmtGene2 = algmtGene2 + endDeletion


        lPosCDSGene1 = []
        lPosCDSGene2 = []

        for cds in self.gene1.lTranscripts[0].lCDS:
            for pos in range(cds.start, cds.end+1):
               lPosCDSGene1.append(pos) 
        for cds in self.gene2.lTranscripts[0].lCDS:
            for pos in range(cds.start, cds.end+1):
               lPosCDSGene2.append(pos)

        # get CDS positions gene 1
        if algmtGene1Strand == -1:
            index = 0
            loopEnd = algmtGene1End
            while loopEnd != algmtGene1Start-1:
                if algmtGene1[index] != '-':
                    if loopEnd in lPosCDSGene1:
                        algmtGene1 = algmtGene1[:index] + algmtGene1[index].upper() + algmtGene1[index+1:]
                    else:
                        algmtGene1 = algmtGene1[:index] + algmtGene1[index].lower() + algmtGene1[index+1:]
                    loopEnd -= 1
                if algmtGene1[index] == '-':
                    next
                index += 1
              
        if  algmtGene1Strand == 1:
            index = 0
            loopEnd = algmtGene1Start
            while loopEnd != algmtGene1End+1:
                if algmtGene1[index] != '-':
                    if loopEnd in lPosCDSGene1:
                        algmtGene1 = algmtGene1[:index] + algmtGene1[index].upper() + algmtGene1[index+1:]
                    else:
                        algmtGene1 = algmtGene1[:index] + algmtGene1[index].lower() + algmtGene1[index+1:]
                    loopEnd += 1
                if algmtGene1[index] == '-':
                    next
                index += 1

        # get CDS positions gene 2
        if algmtGene2Strand == -1:
            index = 0
            loopEnd = algmtGene2End
            while loopEnd != algmtGene2Start-1:
                if algmtGene2[index] != '-':
                    if loopEnd in lPosCDSGene2:
                        algmtGene2 = algmtGene2[:index] + algmtGene2[index].upper() + algmtGene2[index+1:]
                    else:
                        algmtGene2 = algmtGene2[:index] + algmtGene2[index].lower() + algmtGene2[index+1:]
                    loopEnd -= 1
                if algmtGene2[index] == '-':
                    next
                index += 1
              
        if  algmtGene2Strand == 1:
            index = 0
            loopEnd = algmtGene2Start
            while loopEnd != algmtGene2End+1:
                if algmtGene2[index] != '-':
                    if loopEnd in lPosCDSGene2:
                        algmtGene2 = algmtGene2[:index] + algmtGene2[index].upper() + algmtGene2[index+1:]
                    else:
                        algmtGene2 = algmtGene2[:index] + algmtGene2[index].lower() + algmtGene2[index+1:]
                    loopEnd += 1
                if algmtGene2[index] == '-':
                    next
                index += 1
 
        return (algmtGene1,algmtGene2, Region(self.gene1.seqid,algmtGene1Start,algmtGene1End,algmtGene1Strand),Region(self.gene2.seqid,algmtGene2Start,algmtGene2End,algmtGene2Strand))

    def getEffect(self):
        """Return mutations between CDS alignments"""


        algmtGene1, algmtGene2, r1, r2 = self.getCDSAlignment()

        if algmtGene1 == None or algmtGene2 == None:
            return (None, None, None, None)

        if (self.gene1.strand == -1 and r1.strand == 1) or (self.gene1.strand == 1 and r1.strand == -1):
            algmtGene1 = self.reverseComplement(algmtGene1)
        if (self.gene2.strand == -1 and r2.strand == 1) or (self.gene2.strand == 1 and r2.strand == -1):
            algmtGene2 = self.reverseComplement(algmtGene2)

        lStrMutations = []
        lMutations = [' '] * len(algmtGene1)
        lCodonsPosGene1 = [' '] * len(algmtGene1)
        lCodonsPosGene2 = [' '] * len(algmtGene2)
        codon1 = []
        codon2 = []
        lCodon1 = []
        lCodon2 = []

        if self.gene1.strand == 1:
            algmt1Start, algmt1End = (r1.start, r1.end)
        else:
            algmt1Start, algmt1End = (r1.end, r1.start)
        if self.gene2.strand == 1:
            algmt2Start, algmt2End = (r2.start, r2.end)
        else:
            algmt2Start, algmt2End = (r2.end, r2.start)    

        indexAlgmt1 = algmt1Start
        indexAlgmt2 = algmt2Start


        for i,val in enumerate(algmtGene1):
            if algmtGene1[i].upper() != algmtGene2[i].upper():
                lStrMutations.append("{} : {} - {} --- ({}-{},{}-{})\n".format(i+1,algmtGene1[i],algmtGene2[i],r1.seq,indexAlgmt1,r2.seq,indexAlgmt2))
                lMutations[i] = '|'

            if algmtGene1[i] != '-':
                if self.gene1.strand == -1:
                    indexAlgmt1 -= 1
                else:
                    indexAlgmt1 += 1
            if algmtGene2[i] != '-':
                if self.gene2.strand == -1:
                    indexAlgmt2 -= 1
                else:
                    indexAlgmt2 += 1

            if algmtGene1[i].isupper():
                lCodonsPosGene1[i] = '-'
            if algmtGene2[i].isupper():
                lCodonsPosGene2[i] = '-'
            if algmtGene1[i].isupper():
                codon1.append(algmtGene1[i])
                if len(codon1) == 3:
                   lCodon1.append(self.translate(''.join(codon1)))
                   codon1 = []
            if algmtGene2[i].isupper():
                codon2.append(algmtGene2[i])
                if len(codon2) == 3:
                   lCodon2.append(self.translate(''.join(codon2)))
                   codon2 = []

        indexCodon = 1
        lCodonPosAlgmt1 = []
        for i in lCodonsPosGene1:
            if i == '-':
                indexCodon += 1
                if indexCodon == 3:
                   if lCodon1: 
                       lCodonPosAlgmt1.append(lCodon1.pop(0))
                   indexCodon = 0
                else:
                   lCodonPosAlgmt1.append(' ') 
            else:
                lCodonPosAlgmt1.append(' ')

        indexCodon = 1
        lCodonPosAlgmt2 = []
        for i in lCodonsPosGene2:
            if i == '-':
                indexCodon += 1
                if indexCodon == 3:
                   if lCodon2:
                       lCodonPosAlgmt2.append(lCodon2.pop(0))
                   indexCodon = 0
                else:
                   lCodonPosAlgmt2.append(' ') 
            else:
                lCodonPosAlgmt2.append(' ')

         
        lAlignEffect = []
        lAlignEffect.append(''.join(lCodonPosAlgmt1))
        lAlignEffect.append(''.join(lCodonsPosGene1))
        lAlignEffect.append(algmtGene1)
        lAlignEffect.append(''.join(lMutations))
        lAlignEffect.append(algmtGene2)
        lAlignEffect.append(''.join(lCodonsPosGene2))
        lAlignEffect.append(''.join(lCodonPosAlgmt2))

        return (lAlignEffect,lStrMutations,r1,r2)
 
    def reverseComplement(self,algmt):
        """reverse complement table"""

        table = { 'a':'t',
                  't':'a',
                  'c':'g',
                  'g':'c',
                  'A':'T',
                  'T':'A',
                  'C':'G',
                  'G':'C',
                  '-':'-'}
        algmtR = []
        for base in algmt[::-1]:
            algmtR.append(table[base])

        return ''.join(algmtR)


    def translate(self,codon):
        """codon translation table"""

        table = { 'TTT' : 'F',
                  'TTC' : 'F',
                  'TTA' : 'L',
                  'TTG' : 'L',
                  'TCT' : 'S',
                  'TCC' : 'S',
                  'TCA' : 'S',
                  'TCG' : 'S',
                  'TAT' : 'Y',
                  'TAC' : 'Y',
                  'TAA' : '*',
                  'TAG' : '*',
                  'TGT' : 'C',
                  'TGC' : 'C',
                  'TGA' : '*',
                  'TGG' : 'W',
                  'CTT' : 'L',
                  'CTC' : 'L',
                  'CTA' : 'L',
                  'CTG' : 'L',
                  'CCT' : 'P',
                  'CCC' : 'P',
                  'CCA' : 'P',
                  'CCG' : 'P',
                  'CAT' : 'H',
                  'CAC' : 'H',
                  'CAA' : 'Q',
                  'CAG' : 'Q',
                  'CGT' : 'R',
                  'CGC' : 'R',
                  'CGA' : 'R',
                  'CGG' : 'A',
                  'ATT' : 'I',
                  'ATC' : 'I',
                  'ATA' : 'I',
                  'ATG' : 'M',
                  'ACT' : 'T',
                  'ACC' : 'T',
                  'ACA' : 'T',
                  'ACG' : 'T',
                  'AAT' : 'N',
                  'AAC' : 'N',
                  'AAA' : 'K',
                  'AAG' : 'K',
                  'AGT' : 'S',
                  'AGC' : 'S',
                  'AGA' : 'R',
                  'AGG' : 'R',
                  'GTT' : 'V',
                  'GTC' : 'V',
                  'GTA' : 'V',
                  'GTG' : 'V',
                  'GCT' : 'A',
                  'GCC' : 'A',
                  'GCA' : 'A',
                  'GCG' : 'A',
                  'GAT' : 'D',
                  'GAC' : 'D',
                  'GAA' : 'E',
                  'GAG' : 'E',
                  'GGT' : 'G',
                  'GGC' : 'G',
                  'GGA' : 'G',
                  'GGG' : 'G'}

        return table[codon]
