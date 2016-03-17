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

        size = sendMax-sstartMin
        myLength = [0] * size

        for algmt in self.lAlgmts:
            for j in range((algmt.sstart-sstartMin),(algmt.send-sstartMin)):
                myLength[j] = 1

        self.length = myLength.count(1)

        return self.length


    def getIdListOfAlgmts(self):
        """return the ID list of the alignments"""

        return [ algmt.id for algmt in self.lAlgmts ]


    def sortListOfAlgmts(self):
        """Sort the list of Alignments by Sbjct, coordinates"""

        self.lAlgmts.sort(key=lambda algmt: algmt.sstart) 
        return self.lAlgmts


    def convertChain(self, id, format='gff3'):
        """Convert a chain to gff3|bed lines"""

        lAlgmts = self.sortListOfAlgmts()
        lines = []
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

            if format == 'gff3':
                lines.append('{}\tSegDupAna\tmatch_part\t{}\t{}\t.\t{}\t.\tID=match{};Target={} {} {};length={};identities={};identity_percentage={:.2f}\n'.format(algmt.sbjct,algmt.sstart,algmt.send,strand,algmt.id,algmt.query,algmt.qstart,algmt.qend,algmt.length,algmt.identities,(algmt.identities/float(algmt.length))))
            if format == 'bed':
                lines.append('{}\t'.format(algmt.sbjct))
        if format == 'gff3':
            lines.insert(0, '{}\tSegDupAna\tmatch\t{}\t{}\t.\t.\t.\tID=chain{};length={}\n'.format(algmt.sbjct,sstartMin,sendMax,id,self.getLength()))

        return ''.join(lines)

    def __eq__(self, other):
        """Equality on all args"""
      
        return (self.sortListOfAlgmts() == other.sortListOfAlgmts())

    def __repr__(self):
        """Chain representation"""

        return 'lAlgmts: {}'.format(','.join([ str(algmt.id) for algmt in self.lAlgmts ]))



