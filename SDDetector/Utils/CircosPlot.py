#!/usr/bin/env python

import argparse
import logging
import os

from SDDetector.Parser.Gff.GffDuplicationParser import GffDuplicationParser
from SDDetector.Parser.Gff.GffGeneParser import GffGeneParser
from SDDetector.Parser.Gff.GffTEParser import GffTEParser
from SDDetector.Utils.FastaFileIndexer import FastaFileIndexer

class CircosPlot(object):

    def __init__(self, GenomeFile='',SDFile='', GeneFile='', TEFile='', logLevel='ERROR'):
        """Constuctor"""

        self.GenomeFile = GenomeFile
        self.SDFile = SDFile
        self.GeneFile = GeneFile
        self.TEFile = TEFile

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
            f.write('<ticks>\n')
            f.write('radius           = 1r\n')
            f.write('color            = black\n')
            f.write('thickness        = 2p\n')
            f.write('multiplier       = 1e-2\n')
            f.write('format           = %d\n')
            f.write('<tick>\n')
            f.write('spacing        = 0.01u\n')
            f.write('size           = 15p\n')
            f.write('show_label     = yes\n')
            f.write('label_size     = 20p\n')
            f.write('label_offset   = 10p\n')
            f.write('format         = %d\n')
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
            f.write('label_size   = 24p\n')
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
           print i
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
           print modif
           print lRegionsToDraw
       
        print "HAH" 
        print lRegionsToDraw
        print "tete"
        print self.lRegionsToDraw



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
            for seq in parser.lSeq:
                f.write('chr - {} {} {} {} {}\n'.format(seq,seq,0,len(parser.dSeq[seq]),seq))
        f.close()
        return seqdatafile
        

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
 
        genelinkdatafile = 'gene-link.txt'
        #self.lRegionsToDraw = []
        #parser = GffDuplicationParser(self.SDFile)
        #lNonRedDuplications = parser.getNonRedondantDuplications()
        #with open(genelinkdatafile,'w') as f:
        #    for link in lLinks:
        #        f.write('{} {} {} {} {} {}\n'.format(link.seq1, link.start1, link.end1, dup.seq2, dup.start2, dup.end2))
        #        self.lRegionsToDraw.append((dup.seq1,float(dup.start1),float(dup.end1)))
        #        self.lRegionsToDraw.append((dup.seq2,float(dup.start2),float(dup.end2)))
        #f.close()
        
        return genelinkdatafile 

    def writeSimilarityFile(self):
        """Write similarity fiel"""

        similaritydatafile = 'similarity.txt'
        parser = ""
        lSimilarities = ""

 
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

