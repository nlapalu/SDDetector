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

    def __init__(self, logLevel='ERROR'):
        """Constuctor"""
 
        self.lRegionsToDraw = []
        self.lSeqNames = []
        self.lGeneLinks = []
        self.lGenes = []
        self.lTEs = []
        self.similarity = False

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
            f.write('break   = 0.04r\n')
            f.write('axis_break_at_edge = yes\n')
            f.write('axis_break         = yes\n')
            f.write('axis_break_style   = 2\n')
            lSeqToDraw = self.getSequencesToDraw().split(';') 
            f.write('<pairwise {} {}>\n'.format(lSeqToDraw[0],lSeqToDraw[-1]))
            f.write('spacing = 4r\n')
            f.write('</paiwise>\n')

            f.write('<break_style 2>\n')
            f.write('#stroke_color = black\n')
            f.write('stroke_thickness = 5p\n')
            f.write('thickness = 2r\n')
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
            if self.lDuplications:
                f.write('<link>\n')
                f.write('file = segdup.txt\n')
                f.write('color = vvlgrey\n')
                f.write('radius = 0.73r\n')
                f.write('bezier_radius = 0r\n')
                f.write('thickness = 2\n')
                f.write('</link>\n')
            # gene-links
            if self.lGeneLinks:
                f.write('<link>\n')
                f.write('file = gene-link.txt\n')
                f.write('color = blue\n')
                f.write('radius = 0.73r\n')
                f.write('bezier_radius = 0r\n')
                f.write('thickness = 2\n')
                f.write('stroke_color = black\n')
                f.write('stroke_thickness = 1\n')
                f.write('</link>\n')
            f.write('</links>\n')

            # genes
            f.write('<plots>\n')
            if self.lGenes:
                f.write('<plot>\n')
                f.write('type = tile\n')
                f.write('layers_overflow = collapse\n')
                f.write('file = gene.txt\n')
                f.write('r1 = 0.84r\n')
                f.write('r0 = 0.81r\n')
                f.write('layers = 2\n')
                f.write('orientation = center\n')
                f.write('margin = 0.02u\n')
                f.write('thickness = 12\n')
                f.write('padding = 1\n')
                f.write('color = green\n')
                f.write('stroke_color = black\n')
                f.write('stroke_thickness = 1\n')
                f.write('<backgrounds>\n')
                f.write('<background>\n')
                f.write('color = vvlgrey\n')
                f.write('</background>\n')
                f.write('</backgrounds>\n')
                f.write('</plot>\n')
                f.write('<plot>\n')
                f.write('type = text\n')
                f.write('color = black\n')
                f.write('file = gene.txt\n')
                f.write('r0 = 0.84r\n')
                f.write('r1 = 0.99r\n')
                f.write('show_links = yes\n')
                f.write('link_dims = 4p,4p,8p,4p,4p\n')
                f.write('link_thickness = 2p\n')
                f.write('link_color = red\n')
                f.write('label_size = 20p\n')
                f.write('label_font  = condensed\n')
                f.write('padding = 0p\n')
                f.write('rpadding = 0p\n')
                f.write('<backgrounds>\n')
                f.write('<background>\n')
                f.write('color = vvlgrey\n')
                f.write('</background>\n')
                f.write('</backgrounds>\n')
                f.write('</plot>\n')

            # TEs
            if self.lTEs:
                f.write('<plot>\n')
                f.write('type = tile\n')
                f.write('layers_overflow = collapse\n')
                f.write('file = TE.txt\n')
                f.write('r1 = 0.84r\n')
                f.write('r0 = 0.81r\n')
                f.write('layers = 2\n')
                f.write('orientation = center\n')
                f.write('margin      = 0.02u\n')
                f.write('thickness   = 12\n')
                f.write('padding     = 1\n')
                f.write('color = red\n')
                f.write('stroke_color = black\n')
                f.write('stroke_thickness = 1\n')
                f.write('</plot>\n')
            
            # Similarity
            if self.similarity:
                f.write('<plot>\n')
                f.write('type = line\n')
                f.write('thickness = 2\n')
                f.write('file = similarity.txt\n')
                f.write('color = black\n')
                f.write('min = 85\n')
                f.write('max = 102\n')
                f.write('r0 = 0.74r\n')
                f.write('r1 = 0.8r\n')
                f.write('<backgrounds>\n')
                f.write('<background>\n')
                f.write('color     = vvlgreen\n')
                f.write('y0        = 95\n')
                f.write('</background>\n')
                f.write('<background>\n')
                f.write('color     = vvlred\n')
                f.write('y1        = 90\n')
                f.write('</background>\n')
                f.write('</backgrounds>\n')
                f.write('<axes>\n')
                f.write('<axis>\n')
                f.write('color     = dgrey\n')
                f.write('thickness = 1\n')
                f.write('spacing   = 0.1r\n')
                f.write('</axis>\n')
                f.write('</axes>\n')
                f.write('<rules>\n')
                f.write('<rule>\n')
                f.write('condition    = var(value) > 95\n')
                f.write('color        = dgreen\n')
                f.write('#fill_color   = dgreen_a1\n')
                f.write('</rule>\n')
                f.write('<rule>\n')
                f.write('condition    = var(value) < 90\n')
                f.write('color        = dred\n')
                f.write('#fill_color   = dred_a1\n')
                f.write('</rule>\n')
                f.write('</rules>\n')
                f.write('</plot>\n')
            f.write('</plots>\n') 

        f.close()

    def getSequencesToDraw(self):
        """pass"""

        dSequences = { seq[0] : '-' for seq in self.lRegionsToDraw }
        return ';'.join([x for x in self.lSeqNames if x in dSequences])
         

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
             


    def writeSeqDataFile(self, lSeqs, SeqDataFile):
        """Write Sequence data file=Karyotype"""

        self.lSeqNames = [ seq[0] for seq in lSeqs ]
 
        with open(SeqDataFile,'w') as f:
            for i,seq in enumerate(lSeqs):
                f.write('chr - {} {} {} {} {}\n'.format(seq[0],seq[0],0,seq[1],self.getColorByIndex(i)))
        f.close()
        return SeqDataFile


    def getColorByIndex(self,index):
        """Return a color"""

        lColors = ['180,60,45','180,130,45','180,180,45','130,180,50','65,180,50',
                   '45,180,170','45,115,180','115,45,180','195,55,135','194,55,80',
                   '245,25,10','240,230,15','83,245,20','12,245,145','12,235,242',
                   '10,100,242','128,15,242','227,15,242','242,12,104','242,240,208']
        return lColors[index % 20]
        

    def writeSegDupDataFile(self, lDuplications, SDDataFile):
        """Write Segmental Duplication data file"""
 
        self.lRegionsToDraw = []
        self.lDuplications = lDuplications
        with open(SDDataFile,'w') as f:
            for dup in lDuplications:
                f.write('{} {} {} {} {} {}\n'.format(dup.seq1, dup.start1, dup.end1, dup.seq2, dup.start2, dup.end2))
                self.lRegionsToDraw.append((dup.seq1,float(dup.start1),float(dup.end1)))
                self.lRegionsToDraw.append((dup.seq2,float(dup.start2),float(dup.end2)))
        f.close()
        
        return SDDataFile 

    def writeGeneDataFile(self, lGenes, GeneDataFile):
        """Write gene data file"""

        self.lGenes = lGenes
        with open(GeneDataFile,'w') as f:
            for gene in lGenes:
                f.write('{} {} {} {}\n'.format(gene.seqid,gene.start,gene.end,gene.id))
        f.close()

        return GeneDataFile

    def writeGeneLinkDataFile(self, lGeneLinks, GeneLinkDataFile):
        """Write gene links data file"""

        self.lGeneLinks = lGeneLinks
        with open(GeneLinkDataFile,'w') as f:
            for link in lGeneLinks:
                f.write('{} {} {} {} {} {}\n'.format(link.gene1.seqid, link.gene1.start, link.gene1.end, link.gene2.seqid, link.gene2.start, link.gene2.end))
        f.close()
       
        return GeneLinkDataFile


    def writeSimilarityDataFile(self, lDuplications, SimilarityDataFile):
        """Write similarity file"""

        self.lDuplications = lDuplications
        length = 100
        overlap = 0
        lStatus = []
        with open(SimilarityDataFile, 'w') as f:
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
        self.similarity = True
        return SimilarityDataFile
        
           

    def _slidingWindow(self,lStatus,length,overlap):
        """sliding window"""

        lvalues = []
        for i in range(0,len(lStatus),length-overlap):
            if (i+length) < (len(lStatus)-1):
                lvalues.append(sum(lStatus[i:i+length])/float(length)*100)
        return lvalues 

 
    def writeTEDataFile(self, lTEs, TEDataFile):
        """"Write TEs data file"""

        self.lTEs = lTEs
        with open(TEDataFile,'w') as f:
            for TE in lTEs:
                f.write('{} {} {} {}\n'.format(TE.seqid,TE.start,TE.end,TE.id))
        f.close()

        return TEDataFile

