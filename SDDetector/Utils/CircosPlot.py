#!/usr/bin/env python

import argparse
import logging
import os
import math

from SDDetector.Entities.GeneLink import GeneLink
from SDDetector.Parser.Gff.GffDuplicationParser import GffDuplicationParser
from SDDetector.Parser.Gff.GffGeneParser import GffGeneParser
from SDDetector.Parser.Gff.GffTEParser import GffTEParser
from SDDetector.Parser.Blast.BlastXMLParser import BlastXMLParser
from SDDetector.Utils.FastaFileIndexer import FastaFileIndexer
from SDDetector.Db.GeneDB import GeneDB

class CircosPlot(object):

    def __init__(self, GenomeFile='',SDFile='', GeneFile='', TEFile='', BlastXMLFile='', logLevel='ERROR'):
        """Constuctor"""

        self.GenomeFile = GenomeFile
        self.SDFile = SDFile
        self.GeneFile = GeneFile
        self.TEFile = TEFile
        self.BlastXMLFile = BlastXMLFile

    def writeDataFiles(self):
        """Write all data files"""

        self.writeSegDupDataFile(self)

    def writeConfigurationFiles(self):
        """Write all circos configuration files"""
        pass

    def draw(self):
        pass

    def writeCircosConf(self):
        """Write circos.conf file"""

        conf = 'circos.conf'

        with open(conf,'w') as f:

            unit = 100000

            # default parameters
            f.write('<image>\n<<include etc/image.conf>>\n</image>\n')
            f.write('<<include etc/colors_fonts_patterns.conf>>\n')
            f.write('<<include etc/housekeeping.conf>>\n')

            # ideogram
            f.write('<ideogram>\n') 
            f.write('<spacing>\n')
            f.write('default = 0.005r\n')
            f.write('break   = 0.3r\n')
            f.write('axis_break_at_edge = yes\n')
            f.write('axis_break         = yes\n')
            f.write('axis_break_style   = 2\n')
            f.write('<break_style 2>\n')
            f.write('stroke_color     = black\n')
            f.write('stroke_thickness = 5p\n')
            f.write('thickness        = 2r\n')
            f.write('</break>\n')
            f.write('</spacing>\n')
            f.write('radius = 0.9r\n')
            f.write('thickness = 20p\n')
            f.write('fill = yes\n')
            f.write('show_label = yes\n')
            f.write('label_fonts = default\n')
            f.write('label_radius = dims(image,radius)-60p\n')
            f.write('label_size = 30\n')
            f.write('label_parallel = yes\n')
            f.write('</ideogram>\n')

            # ticks
            f.write('show_ticks = yes\n')
            f.write('show_tick_labels = yes\n')
            f.write('<ticks>\n')
            f.write('radius           = 1r\n')
            f.write('color            = black\n')
            f.write('thickness        = 2p\n')
            f.write('multiplier       = 1e-3\n')
            f.write('format           = %d\n')
            f.write('<tick>\n')
            f.write('spacing        = 0.05u\n')
            f.write('size           = 12p\n')
            f.write('thickness      = 2p\n')
            f.write('color = black\n')
            f.write('show_label     = yes\n')
            f.write('label_size     = 14p\n')
            f.write('label_offset   = 3p\n')
            f.write('format         = %d\n')
            f.write('</tick>\n')
            f.write('<tick>\n')
            f.write('spacing = 0.01u\n')
            f.write('size = 4p\n')
            f.write('thickness = 2p\n')
            f.write('color = dgrey\n')
            f.write('show_label = no\n')
            f.write('</tick>\n')
            f.write('</ticks>\n')

            # sequences
            f.write('karyotype = genome.txt\n')
            f.write('chromosomes_units = {}\n'.format(unit))
            f.write('chromosomes = {}\n'.format(self.getSequencesToDraw()))
            f.write('chromosomes_breaks = {}\n'.format(self.getRegionsNotToDraw(unit)))
            f.write('chromosomes_display_default = no\n')

            # links
            f.write('<links>\n')
            f.write('ribbon = yes\n')
            f.write('flat = yes\n')
            # segmental duplications
            f.write('<link>\n')
            f.write('file = segdup.txt\n')
            f.write('color = vvlgrey\n')
            f.write('radius = 0.8r\n')
            f.write('bezier_radius = 0r\n')
            f.write('thickness = 2\n')
            f.write('</link>\n')
            # genes
            f.write('<link>\n')
            f.write('file = gene-link.txt\n')
            f.write('color = blue\n')
            f.write('radius = 0.8r\n')
            f.write('bezier_radius = 0r\n')
            f.write('thickness = 2\n')
            f.write('</link>\n')
            f.write('</links>\n')

            # genes
            f.write('<plots>\n')
            f.write('<plot>\n')
            f.write('type = tile\n')
            f.write('layers_overflow = collapse\n')
            f.write('file = gene.txt\n')
            f.write('r1 = 0.83r\n')
            f.write('r0 = 0.81r\n')
            f.write('layers = 2\n')
            f.write('orientation = center\n')
            f.write('margin      = 0.02u\n')
            f.write('thickness   = 12\n')
            f.write('padding     = 1\n')
            f.write('color = green\n')
            f.write('</plot>\n')

            f.write('<plot>\n')
            f.write('type             = text\n')
            f.write('color            = black\n')
            f.write('file             = gene.txt\n')
            f.write('r0 = 0.84r\n')
            f.write('r1 = 0.99r\n')
            f.write('show_links     = yes\n')
            f.write('link_dims      = 4p,4p,8p,4p,4p\n')
            f.write('link_thickness = 2p\n')
            f.write('link_color     = red\n')
            f.write('label_size   = 20p\n')
            f.write('label_font   = condensed\n')
            f.write('padding  = 0p\n')
            f.write('rpadding = 0p\n')
            f.write('<backgrounds>\n')
            f.write('<background>\n')
            f.write('color = vvlgrey\n')
            f.write('</background>\n')
            f.write('</backgrounds>\n')
            f.write('</plot>\n')


            # genes
            f.write('<plot>\n')
            f.write('type = tile\n')
            f.write('layers_overflow = collapse\n')
            f.write('file = TE.txt\n')
            f.write('r1 = 0.83r\n')
            f.write('r0 = 0.81r\n')
            f.write('layers = 2\n')
            f.write('orientation = center\n')
            f.write('margin      = 0.02u\n')
            f.write('thickness   = 12\n')
            f.write('padding     = 1\n')
            f.write('color = red\n')
            f.write('</plot>\n')
            f.write('</plots>\n') 


        f.close()

    def getSequencesToDraw(self):
        """pass"""

        dSequences = { seq[0] : '-' for seq in self.lRegionsToDraw }
        return ';'.join([x for x in self.lSequences if x in dSequences])
         

    def getRegionsToDraw(self, unit):
        """pass"""

        return ';'.join(sorted([ '{}:{:.3f}-{:.3f}'.format(i[0],(i[1]/unit)-((i[2]-i[1])/(4*unit)),(i[2]/unit)+((i[2]-i[1])/(4*unit))) for i in self.lRegionsToDraw ],key=lambda x : (x[0],x[1])))
       
    def getRegionsNotToDraw(self, unit):
        """pass"""

        dSeq = { seq[0] : [0] for seq in self.lRegionsToDraw }

        lRegionsToDraw = [self.lRegionsToDraw[0]]
        for i in self.lRegionsToDraw[1:]:
           modif = False
           #print i
           for idx,j in enumerate(lRegionsToDraw) :
               if i[0] == j[0]:
                   if (i[1] < j[1] and i[2] < j[1]) or (i[1] > j[2] and i[2] > j[2]):
                       next 
                   elif (i[1] < j[1] and i[2] < j[2] and i[2] > j[1]):
                       #j[1] = i[1]
                       lRegionsToDraw[idx] = (j[0],i[1],j[2])
                       modif = True
                       print 1
                   elif (i[1] > j[1] and i[1] < j[2] and i[2] > j[2]):
                       #j[2] = i[2]
                       lRegionsToDraw[idx] = (j[0],j[1],i[2])
                       modif = True
                       print 2
                   elif (i[1] < j[1] and i[2] > j[2]):
                       #j[1] = i[1]
                       #j[2] = i[2]
                       lRegionsToDraw[idx] = (j[0],i[1],i[2])
                       modif = True
                       print 3
                   elif (i[1] > j[1] and i[2] < j[2]):
                       modif = True
                       print 4
                       next
           if modif == False: 
               lRegionsToDraw.append(i)
           #print modif
           #print lRegionsToDraw
       
        #print "HAH" 
        #print lRegionsToDraw
        #print "tete"
        #print self.lRegionsToDraw



        for i in lRegionsToDraw:
            dSeq[i[0]].extend([(i[1]/unit)-((i[2]-i[1])/(4*unit)),(i[2]/unit)+((i[2]-i[1])/(4*unit))])

        lRegions = []
        for i in dSeq:
            dSeq[i].append(')')
            dSeq[i].sort()
            for v, w in zip(dSeq[i][::2],dSeq[i][1::2]):
                t = '-{}:{:.3f}-'.format(i,v)
                if w == ')':
                    t += '{}'.format(w)
                else:
                    t += '{:.3f}'.format(w)
                lRegions.append(t)

        return ';'.join(lRegions)
             


    def writeSeqDataFile(self):
        """Write Sequence data file=Karyotype"""

        seqdatafile = 'genome.txt'
        parser = FastaFileIndexer(self.GenomeFile)
        parser.read()
        self.lSequences = parser.lSeq
 
        with open(seqdatafile,'w') as f:
            for i,seq in enumerate(parser.lSeq):
                f.write('chr - {} {} {} {} {}\n'.format(seq,seq,0,len(parser.dSeq[seq]),self.getColorByIndex(i)))
        f.close()
        return seqdatafile


    def getColorByIndex(self,index):
        """Return a color"""

        lColors = ['180,60,45','180,130,45','180,180,45','130,180,50','65,180,50',
                   '45,180,170','45,115,180','115,45,180','195,55,135','194,55,80',
                   '245,25,10','240,230,15','83,245,20','12,245,145','12,235,242',
                   '10,100,242','128,15,242','227,15,242','242,12,104','242,240,208']
        return lColors[index % 20]
        

    def writeSegDupDataFile(self):
        """Write Segmental Duplication data file"""
 
        sddatafile = 'segdup.txt'
        self.lRegionsToDraw = []
        parser = GffDuplicationParser(self.SDFile)
        lNonRedDuplications = parser.getNonRedondantDuplications()
        with open(sddatafile,'w') as f:
            for dup in lNonRedDuplications:
                f.write('{} {} {} {} {} {}\n'.format(dup.seq1, dup.start1, dup.end1, dup.seq2, dup.start2, dup.end2))
                self.lRegionsToDraw.append((dup.seq1,float(dup.start1),float(dup.end1)))
                self.lRegionsToDraw.append((dup.seq2,float(dup.start2),float(dup.end2)))
        f.close()
        
        return sddatafile 

    def writeGeneDataFile(self):
        """Write gene data file"""

        genedatafile = 'gene.txt'
        parser =  GffGeneParser(self.GeneFile)
        lGenes = parser.getAllGenes()
        with open(genedatafile,'w') as f:
            for gene in lGenes:
                f.write('{} {} {} {}\n'.format(gene.seqid,gene.start,gene.end,gene.id))
        f.close()

        return genedatafile

    def writeGeneLinkDataFile(self):
        """Write gene links data file"""

        parserDup = GffDuplicationParser(self.SDFile)
        lDuplications = parserDup.getNonRedondantDuplications()
        lRegions = []
        for dup in lDuplications:
            for region in dup.lRegions:
                lRegions.append(region)
 
        parserBlast = BlastXMLParser(self.BlastXMLFile)
        lAlignmentTuples = parserBlast.getAlignmentsFromTupleOfRegions(lRegions)
        #print lAlignmentTuples
        #print lRegions


        index = 0
        for dup in lDuplications:
            lAlgmts = []
            for region in dup.lRegions:
#                print region
#                print index
#                print lAlignmentTuples[index][0]
#                print lAlignmentTuples[index][1]
                lAlgmts.append((lAlignmentTuples[index][0],lAlignmentTuples[index][1]))
                index += 1
            dup.lSeqAlgmts = lAlgmts
            dup.dSeqToSeq = dup.getdSeqToSeq() 

        parserGene = GffGeneParser(self.GeneFile)
        self.db = GeneDB(dbfile='gene.db')
        self.db.insertlGenes(parserGene.getAllGenes())
       
        lLinks = [] 
        for dup in lDuplications:
            (lGeneSeq1,lGeneSeq2) = self._extractGeneInDuplication(dup)
            lLinks.extend(self._buildGeneLinks(lGeneSeq1,lGeneSeq2,dup))


 
        genelinkdatafile = 'gene-link.txt'
        with open(genelinkdatafile,'w') as f:
            for link in lLinks:
                f.write('{} {} {} {} {} {}\n'.format(link.gene1.seqid, link.gene1.start, link.gene1.end, link.gene2.seqid, link.gene2.start, link.gene2.end))
        f.close()
       
        ###TODO### Remove from here
        for link in lLinks:
            # analyse CDS Share Alignment
            print 'Gene: ({},{}); sequence: ({},{}); strand: ({},{})'.format(link.gene1.id, link.gene2.id,link.gene1.seqid,link.gene2.seqid,link.gene1.strand,link.gene2.strand)
            #print 'Algmt: {}'.format(link.getCDSAlignment())
            #print 'Effect: {}'.format(link.getEffect())
            lAlignEffect, nbMutations, r1, r2 = link.getEffect()

            nbBases = len(lAlignEffect[0])
            size = 60
            indexSize = 0
            indexBase = 0
            algmtGene = ''

            if link.gene1.strand == 1:
                algmt1Start, algmt1End = (r1.start, r1.end)
            else:
                algmt1Start, algmt1End = (r1.end, r1.start)
            if link.gene2.strand == 1:
                algmt2Start, algmt2End = (r2.start, r2.end)
            else:
                algmt2Start, algmt2End = (r2.end, r2.start)
            

            start1 = algmt1Start
            start2 = algmt2Start
            end1 = 0
            end2 = 0
            while indexBase < nbBases:
                nbHyphen1 = lAlignEffect[2][indexBase:indexBase+size].count('-')
                nbHyphen2 = lAlignEffect[4][indexBase:indexBase+size].count('-')
                
                if link.gene1.strand == -1:
                    end1 = start1-size-nbHyphen1
                else:
                    end1 = start1+size-nbHyphen1
                if link.gene2.strand == -1:
                    end2 = start2-size-nbHyphen2
                else:
                    end2 = start2+size-nbHyphen2
                
                scale1 = str(start1) + ' '*(size-len(str(start1))-len(str(end1))) + str(end1)
                scale2 = str(start2) + ' '*(size-len(str(start2))-len(str(end2))) + str(end2)


                algmtGene += '{}\n{}\n{}\n{}\n{}\n{}\n{}\n{}\n{}\n\n'.format(scale1,lAlignEffect[0][indexBase:indexBase+size],lAlignEffect[1][indexBase:indexBase+size],lAlignEffect[2][indexBase:indexBase+size],lAlignEffect[3][indexBase:indexBase+size],lAlignEffect[4][indexBase:indexBase+size],lAlignEffect[5][indexBase:indexBase+size],lAlignEffect[6][indexBase:indexBase+size],scale2)
                indexBase += size
                start1 = end1+1
                start2 = end2+1

            print algmtGene


 
        return genelinkdatafile 

    def _buildGeneLinks(self,lGeneSeq1,lGeneSeq2,dup):
        """build"""

        lLinks = []
        for gene1 in lGeneSeq1:
            print gene1.id
            (seq2ID,val1) = dup.dSeqToSeq[gene1.seqid][gene1.start]
            (seq2ID,val2) = dup.dSeqToSeq[gene1.seqid][gene1.end]
            seq2Start = min(val1,val2)
            seq2End = max(val1,val2)
#            print  'start {} end {}'.format(seq2Start,seq2End)
            for gene2 in lGeneSeq2:
#                print 'gene2 {} start {} end {}'.format(gene2.seqid, gene2.start, gene2.end)
                if (gene2.start < seq2Start and gene2.end < seq2Start) or (gene2.start > seq2End and gene2.end > seq2End):
                    next
                else:
                   print 'gene1 -  gene2 : {} {}'.format(gene1.id,gene2.id)
                   lLinks.append(GeneLink(dup=dup,gene1=gene1,gene2=gene2)) 
        return lLinks        
        # todo set : + logging.debug
      

    def _extractGeneInDuplication(self,dup):
        """extract """

        lGeneSeq1 = self.db.getlGenesFromCoordinates(dup.seq1,dup.start1,dup.end1)
        lGeneSeq2 = self.db.getlGenesFromCoordinates(dup.seq2,dup.start2,dup.end2)

        return (lGeneSeq1,lGeneSeq2)


    def writeSimilarityDataFile(self):
        """Write similarity fiel"""

        similaritydatafile = 'similarity.txt'
        parserDup = GffDuplicationParser(self.SDFile)
        lDuplications = parserDup.getNonRedondantDuplications()
        lRegions = []
        for dup in lDuplications:
            for region in dup.lRegions:
                lRegions.append(region)

        parserBlast = BlastXMLParser(self.BlastXMLFile)
        lAlignmentTuples = parserBlast.getAlignmentsFromTupleOfRegions(lRegions)

        index = 0
        for dup in lDuplications:
            lAlgmts = []
            for region in dup.lRegions:
                lAlgmts.append((lAlignmentTuples[index][0],lAlignmentTuples[index][1]))
                index += 1
            dup.lSeqAlgmts = lAlgmts
            dup.dSeqToSeq = dup.getdSeqToSeq() 

        length = 100
        overlap = 0
        lStatus = []
        with open(similaritydatafile, 'w') as f:
            for dup in lDuplications:
                for i,(reg1,reg2) in enumerate(dup.lRegions):
                    for j,d in enumerate(dup.lSeqAlgmts[i][0]):
                        if dup.lSeqAlgmts[i][0][j] != dup.lSeqAlgmts[i][1][j]:
                            lStatus.append(0)
                        else:
                            lStatus.append(1)

                    lvalues = self._slidingWindow(lStatus,length,overlap)
                    for h,value in enumerate(lvalues):
                        if reg1.strand == 1:
                            f.write("{} {} {} {}\n".format(reg1.seq,reg1.start+(h*(length-overlap)),reg1.start+(h*(length-overlap))+length,value))
                        elif reg1.strand == -1:
                            f.write("{} {} {} {}\n".format(reg1.seq,reg1.end-(h*(length-overlap)),reg1.end-(h*(length-overlap))-length,value))
                    for h,value in enumerate(lvalues):
                        if reg2.strand == 1:
                            f.write("{} {} {} {}\n".format(reg2.seq,reg2.start+(h*(length-overlap)),reg2.start+(h*(length-overlap))+length,value))
                        elif reg1.strand == -1:
                            f.write("{} {} {} {}\n".format(reg2.seq,reg2.end-(h*(length-overlap)),reg2.end-(h*(length-overlap))-length,value))

                    lStatus = []

        f.close()
        return similaritydatafile
        
           

    def _slidingWindow(self,lStatus,length,overlap):
        """sliding window"""

        lvalues = []
        for i in range(0,len(lStatus),length-overlap):
            if (i+length) < (len(lStatus)-1):
                lvalues.append(sum(lStatus[i:i+length])/float(length)*100)
        return lvalues 

 
    def writeTEDataFile(self):
        """"Write TEs data file"""

        tedatafile = 'TE.txt'
        parser = GffTEParser(self.TEFile)
        lTEs = parser.getAllTEs()
        with open(tedatafile,'w') as f:
            for TE in lTEs:
                f.write('{} {} {} {}\n'.format(TE.seqid,TE.start,TE.end,TE.id))
        f.close()

        return tedatafile

