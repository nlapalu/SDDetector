#!/usrbin/env python

class Chain(object):


    def __init__(self, lAlgmts=[], id=''):
        """Chain constructor"""

        self.lAlgmts = lAlgmts
        self.id = id
        self.length = None

    def getLength(self):
        """Compute and return chain length"""

        length = 0
        for i, algmt in enumerate(self.lAlgmts):
            if i == 0:
               sstartMin = algmt.sstart
               sendMax = algmt.send 
            if algmt.sstart < sstartMin:
                sstartMin = algmt.sstart
            if algmt.send > sendMax:
                sendMax = algmt.send

        size = sendMax-sstartMin + 1
        myLength = [0] * size

        for algmt in self.lAlgmts:
            for j in range((algmt.sstart-sstartMin),(algmt.send-sstartMin + 1)):
                myLength[j] = 1

        self.length = myLength.count(1)

        return self.length

    def getSStart(self):

        return min([ algmt.sstart for algmt in self.lAlgmts ])

    def getSEnd(self):

        return max([ algmt.send for algmt in self.lAlgmts ])

    def getQStart(self):

        return min([ algmt.qstart for algmt in self.lAlgmts ])

    def getQEnd(self):

        return max([ algmt.qend for algmt in self.lAlgmts ])

    def getAlgmtMaxLength(self):
        """return the size of the longest alignment"""

        return max([ algmt.length for algmt in self.lAlgmts ])
        

    def getNbAlgmts(self):
        """return the number of alignments"""

        return len(self.lAlgmts)


    def getNbGaps(self):
        """Compute nb Gaps"""

        nbGaps = 0
        for algmt in self.lAlgmts:
            nbGaps += algmt.length*2 - ((algmt.qend-algmt.qstart) + (algmt.send-algmt.sstart) + 2  )
        return nbGaps

    def getNbSNPs(self):
        """Compute SNP"""

        nbSNPs = 0
        for algmt in self.lAlgmts:
            nbSNPs += (algmt.length - algmt.identities)
        
        return nbSNPs - self.getNbGaps()

    def getIdListOfAlgmts(self):
        """return the ID list of the alignments"""

        return [ algmt.id for algmt in self.lAlgmts ]


    def sortListOfAlgmts(self):
        """Sort the list of Alignments by Sbjct, coordinates"""

        self.lAlgmts.sort(key=lambda algmt: algmt.sstart) 
        return self.lAlgmts

    def deleteListOfAlgmts(self,lIds):
        """Delete a list of alignments from the chain"""

        self.lAlgmts = [ algmt for algmt in self.lAlgmts if algmt.id not in lIds ]


    def convertChain(self, id, format='gff3'):
        """Convert a chain to gff3|bed lines"""

        lAlgmts = self.sortListOfAlgmts()
        lines = []
        identities = 0
        minSeqs = 0
        for i,algmt in enumerate(lAlgmts):
            if i == 0:
               sstartMin = algmt.sstart
               sendMax = algmt.send 
            if algmt.sstart < sstartMin:
                sstartMin = algmt.sstart
            if algmt.send > sendMax:
                sendMax = algmt.send
            if algmt.sstrand == 1:
                strand = '+'
            else:
                strand = '-'

            minSeq = min(algmt.qend-algmt.qstart,algmt.send-algmt.sstart) + 1
            minSeqs +=minSeq
            identities += algmt.identities
 
            if format == 'gff3':
                lines.append('{}\tSDDetector\tmatch_part\t{}\t{}\t.\t{}\t.\tID=match{};Parent=chain{};Target={} {} {};length={};identities={};identity_percentage={:.3f}\n'.format(algmt.sbjct,algmt.sstart,algmt.send,strand,algmt.id,id,algmt.query,algmt.qstart,algmt.qend,algmt.length,algmt.identities,(algmt.identities/float(minSeq))))
            if format == 'bed':
                lines.append('{}\t'.format(algmt.sbjct))
        if format == 'gff3':
            lines.insert(0, '{}\tSDDetector\tmatch\t{}\t{}\t.\t.\t.\tID=chain{};length={};nbIndels={};nbSNPs={};identity_percentage={:.3f}\n'.format(algmt.sbjct,sstartMin,sendMax,id,self.getLength(),self.getNbGaps(),self.getNbSNPs(),identities/float(minSeqs)))

        return ''.join(lines)

    def __eq__(self, other):
        """Equality on all args"""
      
        return (self.sortListOfAlgmts() == other.sortListOfAlgmts())

    def __repr__(self):
        """Chain representation"""

        return 'lAlgmts: {}'.format(','.join([ str(algmt.id) for algmt in self.lAlgmts ]))



