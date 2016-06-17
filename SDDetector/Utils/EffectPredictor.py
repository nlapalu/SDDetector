#!/usr/bin/env python


class EffectPredictor(object):

    effects = []
    # structure
    # AA changes
    # snp 
    # indel 

    def __init__(self):
        """Constructor"""

       # = geneLink


    def run(self):
        pass


    def extractVariantFromAlignment(self,seq1,seq2):
        """def

            Sequence analyzed only on CDS part.
            If gene splitting, the CDS part is analyzed in 
            case of a symetric CDS
        """

        nbSNPs = 0
        nbInDels = 0
        seq1PreviousNt = None
        seq2PreviousNt = None
        for i in range(0, len(seq1)):
            if seq1[i] != seq2[i] and seq1[i] != '.' and seq2[i] != '.':
                nbSNPs += 1
            if seq1[i] == '.' and seq1PreviousNt != '.':
                nbInDels += 1
            if seq2[i] == '.' and seq2PreviousNt != '.':
                nbInDels += 1
            seq1PreviousNt = seq1[i]
            seq2PreviousNt = seq2[i]

        return (nbSNPs,nbInDels)




    def basicStructureAnalysis(self,gene1,gene2):
        """ exons structure

            possible return value: 
            - None : no effect, same structure nb exon/ size
            - SIZE : At least one exon as not the same size, the number of exon is the same
            - EXON : The number of exons is different, and possibly the size
        """

        lGene1CDSSizes = []
        lGene2CDSSizes = [] 
        for CDS in gene1.lTranscripts[0].lCDS:
            lGene1CDSSizes.append(CDS.end-CDS.start)
        for CDS in gene2.lTranscripts[0].lCDS:
            lGene2CDSSizes.append(CDS.end-CDS.start)
        if gene1.lTranscripts[0].strand == -1:
            lGene1CDSSizes.reverse()
        if gene2.lTranscripts[0].strand == -1:
            lGene2CDSSizes.reverse()
        
        if len(lGene1CDSSizes) == len(lGene2CDSSizes):
            if sum(lGene1CDSSizes) != sum(lGene2CDSSizes):
                return 'SIZE'
            else:
                return None
        else:
            return 'EXON'
 
