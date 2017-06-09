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

    def getCDSAlignmentOld(self):
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
        indexStart = self.gene1.start
        indexEnd = self.gene1.end
        while (startReverse == None):
            indexStart += 1
            startReverse = self.dup.dSeqToSeq[self.gene1.seqid][indexStart][1]
            startDeletion +='-'
        while (endReverse == None):
            indexEnd -= 1
            endReverse = self.dup.dSeqToSeq[self.gene1.seqid][indexEnd][1]
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
        algmtGene2 = algmtGene2Strand = None
        if algmtGene2Start == None or algmtGene2End == None:

            val1 = self.dup.dSeqToSeq[self.gene1.seqid][start][1]
            val2 = self.dup.dSeqToSeq[self.gene1.seqid][end][1]
            deltaVal1 = 0
            while val1 == None:
                deltaVal1 -=1
                val1 = self.dup.dSeqToSeq[self.gene1.seqid][start-deltaVal1][1]
                print "val1:{}".format(val1)
            deltaVal2 = 0
            while val2 == None:
                deltaVal2 -=1
                val2 = self.dup.dSeqToSeq[self.gene1.seqid][end-deltaVal2][1]
                print "val2:{}".format(val2)
            deltaGene2Start = deltaGene2End = 0
            
            if val1 < val2:
                algmtGene2Start = val1
                deltaGene2Start = deltaVal1
                algmtGene2End = val2
                deltaGene2End = deltaVal2
            else:
                algmtGene2Start = val2
                deltaGene2Start = deltaVal2
                algmtGene2End = val1
                deltaGene2End = deltaVal1
            
            algmtGene2, algmtGene2Strand = self.dup.getSeqAlignmentWithIndels(self.gene2.seqid,algmtGene2Start,algmtGene2End,deltaGene2Start, deltaGene2End)
            print len(algmtGene2), algmtGene2Start, algmtGene2End,deltaGene2Start, deltaGene2End
        else:
            algmtGene2, algmtGene2Strand = self.dup.getSeqAlignment(self.gene2.seqid,algmtGene2Start,algmtGene2End)

        algmtGene1, algmtGene1Strand = self.dup.getSeqAlignment(self.gene1.seqid,algmtGene1Start,algmtGene1End)

        
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

        # clean alignment with "-" in beginning
        #if len(algmtGene1) != len(algmtGene2):
#            if algmtGene1[0] == '-':
#                algmtGene1 = algmtGene1[1:]
#            if algmtGene2[0] == '-':
#                algmtGene2 = algmtGene2[1:]
      #      while algmtGene1[0] == '-':
      #          algmtGene1 = algmtGene1[1:]
      #      while algmtGene2[0] == '-':
      #          algmtGene2 = algmtGene2[1:]

        if len(algmtGene1) != len(algmtGene2):
                print "\n"
                print algmtGene1
                print algmtGene2
                print "\n"

        while len(algmtGene1) != len(algmtGene2):
            if algmtGene1[0] == '-':

                algmtGene1 = algmtGene1[1:]
            if algmtGene2[0] == '-':
                algmtGene2 = algmtGene2[1:]



        if len(algmtGene1) != len(algmtGene2):
            print algmtGene1
            print algmtGene2
        

 
        return (algmtGene1,algmtGene2, Region(self.gene1.seqid,algmtGene1Start,algmtGene1End,algmtGene1Strand),Region(self.gene2.seqid,algmtGene2Start,algmtGene2End,algmtGene2Strand))


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
        indexStart = self.gene1.start
        indexEnd = self.gene1.end
        while (startReverse == None):
            indexStart += 1
            startReverse = self.dup.dSeqToSeq[self.gene1.seqid][indexStart][1]
            startDeletion +='-'
        while (endReverse == None):
            indexEnd -= 1
            endReverse = self.dup.dSeqToSeq[self.gene1.seqid][indexEnd][1]
            endDeletion +='-'

        if startReverse > endReverse:
            strand = -1
            
        if strand == -1:
            if startReverse > self.gene2.end:
                start = self.gene1.start
            else:
                start =  self.dup.dSeqToSeq[self.gene2.seqid][self.gene2.end][1]
                endDeletion = ''
            if endReverse > self.gene2.start:
                end = self.dup.dSeqToSeq[self.gene2.seqid][self.gene2.start][1]
                startDeletion = ''
            else:
                end = self.gene1.end
        else:

            if startReverse > self.gene2.start:
                start = self.dup.dSeqToSeq[self.gene2.seqid][self.gene2.start][1]
                startDeletion = ''
            else:
                start = self.gene1.start
            if endReverse > self.gene2.end:
                end = self.gene1.end
            else:
                end = self.dup.dSeqToSeq[self.gene2.seqid][self.gene2.end][1]
                endDeletion = ''


        algmtGene1Start = start
        algmtGene1End = end

        if algmtGene1Start == None or algmtGene1End == None:
            logging.error('BUG implementation with this condition - please contact me with your data')
            sys.exit(0)


        indexStart = start
        indexEnd = end 
        startReverse2 = self.dup.dSeqToSeq[self.gene1.seqid][indexStart][1]
        endReverse2 = self.dup.dSeqToSeq[self.gene1.seqid][indexEnd][1]

        while (startReverse2 == None):
            indexStart += 1
            startReverse2 = self.dup.dSeqToSeq[self.gene1.seqid][indexStart][1]
        while (endReverse2 == None):
            indexEnd -= 1
            endReverse2 = self.dup.dSeqToSeq[self.gene1.seqid][indexEnd][1]


#        algmtGene2Start = min(self.dup.dSeqToSeq[self.gene1.seqid][start][1],self.dup.dSeqToSeq[self.gene1.seqid][end][1])
#        algmtGene2End = max(self.dup.dSeqToSeq[self.gene1.seqid][start][1],self.dup.dSeqToSeq[self.gene1.seqid][end][1])
        algmtGene2Start = min(startReverse2, endReverse2)
        algmtGene2End = max(startReverse2, endReverse2)

        #FIX-BUG - None pos if deletion.
        #print "algmtGene2Start {}:  val: {} - {}".format(algmtGene2Start,self.dup.dSeqToSeq[self.gene1.seqid][start][1],self.dup.dSeqToSeq[self.gene1.seqid][end][1])
 
        #print "algmtGene2End {}:  val: {} - {}".format(algmtGene2End,self.dup.dSeqToSeq[self.gene1.seqid][start][1],self.dup.dSeqToSeq[self.gene1.seqid][end][1])


#        algmtGene2 = algmtGene2Strand = None
#        if algmtGene2Start == None or algmtGene2End == None:

#            val1 = self.dup.dSeqToSeq[self.gene1.seqid][start][1]
#            val2 = self.dup.dSeqToSeq[self.gene1.seqid][end][1]
#            deltaVal1 = 0
#            while val1 == None:
#                deltaVal1 -=1
#                val1 = self.dup.dSeqToSeq[self.gene1.seqid][start-deltaVal1][1]
#                print "val1:{}".format(val1)
#            deltaVal2 = 0
#            while val2 == None:
#                deltaVal2 -=1
#                val2 = self.dup.dSeqToSeq[self.gene1.seqid][end-deltaVal2][1]
#                print "val2:{}".format(val2)
#            deltaGene2Start = deltaGene2End = 0
            
#            if val1 < val2:
#                algmtGene2Start = val1
#                deltaGene2Start = deltaVal1
#                algmtGene2End = val2
#                deltaGene2End = deltaVal2
#            else:
#                algmtGene2Start = val2
#                deltaGene2Start = deltaVal2
#                algmtGene2End = val1
#                deltaGene2End = deltaVal1
            
#            algmtGene2, algmtGene2Strand = self.dup.getSeqAlignmentWithIndels(self.gene2.seqid,algmtGene2Start,algmtGene2End,deltaGene2Start, deltaGene2End)
#            print len(algmtGene2), algmtGene2Start, algmtGene2End,deltaGene2Start, deltaGene2End
#        else:
#            algmtGene2, algmtGene2Strand = self.dup.getSeqAlignment(self.gene2.seqid,algmtGene2Start,algmtGene2End)

#        algmtGene1, algmtGene1Strand = self.dup.getSeqAlignment(self.gene1.seqid,algmtGene1Start,algmtGene1End)


        
        algmtGene1, algmtGene2, algmtGene1Strand, algmtGene2Strand  = self.dup.getSeqAlignmentNew(self.gene1.seqid,algmtGene1Start,algmtGene1End)


        
        if not algmtGene1 or not algmtGene2:
            return (None,None,None,None)

#        if startDeletion:
#            algmtGene2 = startDeletion + algmtGene2
#        if endDeletion:
#            algmtGene2 = algmtGene2 + endDeletion


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

        # clean alignment with "-" in beginning
        #if len(algmtGene1) != len(algmtGene2):
#            if algmtGene1[0] == '-':
#                algmtGene1 = algmtGene1[1:]
#            if algmtGene2[0] == '-':
#                algmtGene2 = algmtGene2[1:]
      #      while algmtGene1[0] == '-':
      #          algmtGene1 = algmtGene1[1:]
      #      while algmtGene2[0] == '-':
      #          algmtGene2 = algmtGene2[1:]

        if len(algmtGene1) != len(algmtGene2):
            logging.error("Algmts not same length")
            sys.exit(1)




 
        return (algmtGene1,algmtGene2, Region(self.gene1.seqid,algmtGene1Start,algmtGene1End,algmtGene1Strand),Region(self.gene2.seqid,algmtGene2Start,algmtGene2End,algmtGene2Strand))

    def getEffect(self):
        """Return mutations between CDS alignments"""


        algmtGene1, algmtGene2, r1, r2 = self.getCDSAlignment()

        #print algmtGene1
        #print algmtGene2
        #print r1.strand
        #print r2.strand

        if algmtGene1 == None or algmtGene2 == None:
            return (None, None, None, None)

#        if (self.gene1.strand == -1 and r1.strand == 1) or (self.gene1.strand == 1 and r1.strand == -1):
#            algmtGene1 = self.reverseComplement(algmtGene1)
#        if (self.gene2.strand == -1 and r2.strand == 1) or (self.gene2.strand == 1 and r2.strand == -1):
#            algmtGene2 = self.reverseComplement(algmtGene2)

        lStrMutations = []
        lMutations = [' '] * len(algmtGene1)
        lCodonsPosGene1 = [' '] * len(algmtGene1)
        lCodonsPosGene2 = [' '] * len(algmtGene2)
        codon1 = []
        codon2 = []
        lCodon1 = []
        lCodon2 = []

#        if self.gene1.strand == 1:
#            algmt1Start, algmt1End = (r1.start, r1.end)
#        else:
#            algmt1Start, algmt1End = (r1.end, r1.start)
#        if self.gene2.strand == 1:
#            algmt2Start, algmt2End = (r2.start, r2.end)
#        else:
#            algmt2Start, algmt2End = (r2.end, r2.start)    


        if r1.strand == 1:
            algmt1Start, algmt1End = (r1.start, r1.end)
        else:
            algmt1Start, algmt1End = (r1.end, r1.start)
        if r2.strand == 1:
            algmt2Start, algmt2End = (r2.start, r2.end)
        else:
            algmt2Start, algmt2End = (r2.end, r2.start)    

        indexAlgmt1 = algmt1Start
        indexAlgmt2 = algmt2Start


        # debug to remove
        if len(algmtGene1) != len(algmtGene2):
            print len(algmtGene1)
            print len(algmtGene2)
 

        # debug to remove


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
#            if algmtGene1[i].isupper():
#                codon1.append(algmtGene1[i])
#                if len(codon1) == 3:
#                   lCodon1.append(self.translate(''.join(codon1)))
#                   codon1 = []
#            if algmtGene2[i].isupper():
#                codon2.append(algmtGene2[i])
#                if len(codon2) == 3:
#                   lCodon2.append(self.translate(''.join(codon2)))
#                   codon2 = []


#v2
#            if algmtGene1[i].isupper():
#                codon1.append(algmtGene1[i])
#                if len(codon1) == 3:
#                   if r1.strand != self.gene1.strand:
#                       lCodon1.append(self.translate(self.reverseComplement(''.join(codon1))))
#                   else:
#                       lCodon1.append(self.translate(''.join(codon1)))
#                   codon1 = []
#            if algmtGene2[i].isupper():
#                codon2.append(algmtGene2[i])
#                if len(codon2) == 3:
#                   if self.gene2.id == "Lmb_jn3_12525_G":
#                       print codon2
#                   if r2.strand != self.gene2.strand:
#                       lCodon2.append(self.translate(self.reverseComplement(''.join(codon2))))
#                   else:
#                       lCodon2.append(self.translate(''.join(codon2)))
#                   codon2 = []

            if algmtGene1[i].isupper():
                codon1.append(algmtGene1[i])
                #if len(codon1) == 3:
                #    lCodon1.append(''.join(codon1))
                #    codon1 = []
            if algmtGene2[i].isupper():
                codon2.append(algmtGene2[i])
                #if len(codon2) == 3:
                #    lCodon2.append(''.join(codon2))
                #    codon2 = []

        if r1.strand != self.gene1.strand:
            #print [ i for i in range(len(lCodon1),0,-3)]
            #lCodon1 = [ self.translate(self.reverseComplement(lCodon1[i])) for i in range(len(lCodon1),0,-1)]
            #codon1 = list(self.reverseComplement(''.join(codon1)))

            #if self.gene1.id == "Lmb_jn3_01284_G":
            #    print codon1
            #print [ codon1[i:i-3] for i in range(len(codon1),0,-3) ]
            #lCodon1 = [ self.translate(self.reverseComplement(''.join(codon1[i:i-3]))) for i in range(len(codon1),0,-3) ]
            for i in range(len(codon1),0,-3):
                codon = ''.join(codon1[max(0,i-3):i])
                if len(codon) == 3:
                    lCodon1.append(self.translate(self.reverseComplement(codon)))
                if len(codon) > 0 and len(codon) < 3:
                    lCodon1.append('')

            lCodon1 = lCodon1[::-1] 
        else:
            #print [ i for i in range(len(lCodon1),0,-3)]
            #lCodon1 = [ self.translate(lCodon1[i]) for i in range(0,len(lCodon1))]
            #lCodon1 = [ self.translate(i) for i in lCodon1]
            #lCodon1 = [ self.translate(''.join(codon1[i:i+3])) for i in range(0,len(codon1),3) ]

            for i in range(0,len(codon1),3):
                codon = ''.join(codon1[i:i+3])
                if len(codon) == 3:
                    lCodon1.append(self.translate(codon))
                if len(codon) > 0 and len(codon) < 3:
                    lCodon1.append('')
            

        if r2.strand != self.gene2.strand:
            #lCodon2 = [ self.translate(self.reverseComplement(lCodon2[i])) for i in range(len(lCodon2),0,-1)]
            #lCodon2 = [ self.translate(self.reverseComplement(i)) for i in lCodon2[::-1] ]
#            lCodon2 = [ self.translate(self.reverseComplement(''.join(codon2[i:i-3]))) for i in range(len(codon2),0,-3) ]
            for i in range(len(codon2),0,-3):
                codon = ''.join(codon2[max(0,i-3):i])
                if len(codon) == 3:
                    lCodon2.append(self.translate(self.reverseComplement(codon)))
                if len(codon) > 0 and len(codon) < 3:
                    lCodon2.append('')

            lCodon2 = lCodon2[::-1] 
        else:
            #lCodon2 = [ self.translate(lCodon2[i]) for i in range(0,len(lCodon2))]
            #lCodon2 = [ self.translate(i) for i in lCodon2 ]
 #           lCodon2 = [ self.translate(''.join(codon2[i:i+3])) for i in range(0,len(codon2),3) ]
            


            for i in range(0,len(codon2),3):
                codon = ''.join(codon2[i:i+3])
                if len(codon) == 3:
                    lCodon2.append(self.translate(codon))
                if len(codon) > 0 and len(codon) < 3:
                    lCodon2.append('')
            

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

        # debug prot
        if self.gene1.strand != r1.strand:
            print ">{}\n{}".format(self.gene1.id,''.join(lCodonPosAlgmt1[::-1]).replace(' ',''))
        else:
            print ">{}\n{}".format(self.gene1.id,''.join(lCodonPosAlgmt1).replace(' ',''))
        if self.gene2.strand != r2.strand:
            print ">{}\n{}".format(self.gene2.id,''.join(lCodonPosAlgmt2[::-1]).replace(' ',''))
        else:
            print ">{}\n{}".format(self.gene2.id,''.join(lCodonPosAlgmt2).replace(' ',''))



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
                  'CGG' : 'R',
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
