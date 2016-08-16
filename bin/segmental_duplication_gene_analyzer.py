#!/usr/bin/env python

import argparse
import logging
import os
import sys

from SDDetector.Entities.Region import Region
from SDDetector.Entities.GeneLink import GeneLink

from SDDetector.Parser.Gff.GffDuplicationParser import GffDuplicationParser
from SDDetector.Parser.Gff.GffGeneParser import GffGeneParser
from SDDetector.Parser.Gff.GffTEParser import GffTEParser
from SDDetector.Parser.Blast.BlastXMLParser import BlastXMLParser

from SDDetector.Db.GeneDB import GeneDB

from SDDetector.Utils.CircosPlot import CircosPlot
from SDDetector.Utils.FastaFileIndexer import FastaFileIndexer

try:
    from SDDetector.Parser.Blast.BlastXMLParser import BlastXMLParser
except ImportError:
    raise Exception('BioPython is not installed, xml parsing not available. \
                    Please install BioPython if you want to use this tool')
    exit(1)


class Analyzer(object):

    def __init__(self, SDFile='', BlastXMLFile='', GeneFile='', outputFile='', logLevel='ERROR'):
        """Constructor"""

        self.SDFile = SDFile
        self.BlastXMLFile = BlastXMLFile
        self.GeneFile = GeneFile
        self.outputFile = outputFile
        self.logLevel = logLevel
        logging.basicConfig(level=self.logLevel)

        self.lDuplications = []


    def __del__(self):
        """Destructor"""

        if os.path.exists('gene.db'): 
            os.remove('gene.db')


    def parseAttributesFromArgsCLI(self):
        """Parse arguments from command line"""

        program = 'SDAnalyzer'
        version = 0.1
        description = "SDAnalyzer: analyzes segmental duplications in genome"

        parser = argparse.ArgumentParser(prog=program)
        parser = argparse.ArgumentParser(description=description)
        parser.add_argument('--version', action='version', version='{} {}'.format(program,version))

        parser.add_argument("SDFile", help="Segmental Duplication gff3 file = output of SDDetector (filtered or not)", type=str)
        parser.add_argument("BlastXMLFile", help="Input Blast XML file", type=str)
        parser.add_argument("GeneFile", help="gene annotation file in gff3 format", type=str)
        parser.add_argument("outputFile", help="Output File", type=str)

        parser.add_argument("-v", "--verbosity", type=int, choices=[1,2,3],
                            help="increase output verbosity 1=error, 2=info, 3=debug")
        parser.add_argument("-t", "--TEFile", type=str, default=None, help="Transposable \
                            elements / Repeats file in gff3 format")
        parser.add_argument("-g", "--GenomeFile", type=str, default=None, help="Genome \
                            fasta file, required for circos plot")
        parser.add_argument("--circos", action="store_true", help="Write circos \
                            configuration file and associated data files")

        args = parser.parse_args()
        self._setAttributesFromArgsCLI(args)


    def _setAttributesFromArgsCLI(self, args):
        """Set attributes from argparse"""
 
        if args.verbosity == 1:
            self.logLevel = 'ERROR'
        if args.verbosity == 2:
            self.logLevel = 'INFO'
        if args.verbosity == 3:
            self.logLevel = 'DEBUG'
        logging.getLogger().setLevel(self.logLevel)

        self.SDFile = args.SDFile
        self.BlastXMLFile = args.BlastXMLFile
        self.GeneFile = args.GeneFile
        self.TEFile = args.TEFile
        self.circos = args.circos
        self.GenomeFile = args.GenomeFile
        self.outputFile = args.outputFile


    def getPolymorphismEffect(self):
        """Analyze polymorphism between genes and return list of variants"""
       
        with open(self.outputFile,'w') as f:
            logging.info('Writing polymorphism effect in {}'.format(self.outputFile))
            for link in self.lGeneLinks:
                # analyse CDS Share Alignment
                f.write('Gene: ({},{}); sequence: ({},{}); strand: ({},{})'.format(link.gene1.id, link.gene2.id,link.gene1.seqid,link.gene2.seqid,link.gene1.strand,link.gene2.strand))
                #print 'Algmt: {}'.format(link.getCDSAlignment())
                #print 'Effect: {}'.format(link.getEffect())
                lAlignEffect, lMutations, r1, r2 = link.getEffect()
                for strMutation in lMutations:
                    f.write(strMutation)

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

                f.write(algmtGene)
        f.close()




    def runAnalyze(self):
        """run analyze"""

        logging.info('Parsing Duplication gff file')
        iGffDuplicationParser = GffDuplicationParser(self.SDFile)
        self.lDuplications = iGffDuplicationParser.getNonRedondantDuplications()
        lRegions = []
        for dup in self.lDuplications:
            for region in dup.lRegions:
                lRegions.append(region)

        logging.info('Parsing Blast xml file')
        iBlastXMLParser = BlastXMLParser(self.BlastXMLFile)
        lAlignmentTuples = iBlastXMLParser.getAlignmentsFromTupleOfRegions(lRegions)
        #print lAlignmentTuples

        index = 0
        for dup in self.lDuplications:
            lAlgmts = []
            for region in dup.lRegions:
                ##TODO: bug here if blast.xml has not sdd regions
                lAlgmts.append((lAlignmentTuples[index][0],lAlignmentTuples[index][1]))
                index += 1
            dup.lSeqAlgmts = lAlgmts
            dup.dSeqToSeq = dup.getdSeqToSeq() 
                
        logging.info('Parsing Gene gff file')
        iGffGeneParser = GffGeneParser(self.GeneFile)
        self.db = GeneDB(dbfile='gene.db')
        self.lGenes = iGffGeneParser.getAllGenes()
        self.db.insertlGenes(self.lGenes)
        
        self.lGeneLinks = [] 
        for dup in self.lDuplications:
            (lGeneSeq1,lGeneSeq2) = self._extractGeneInDuplication(dup)
            self.lGeneLinks.extend(self._buildGeneLinks(lGeneSeq1,lGeneSeq2,dup))

        self.getPolymorphismEffect()

        if self.circos:

            if self.GenomeFile:
                logging.info('Indexing Genome fasta file')
                iFastaGenomeParser = FastaFileIndexer(self.GenomeFile)
                iFastaGenomeParser.read()
                ##TODO: modify parser : implement an iterator on sequence/ obj seq with size ?
                lSeqNames = iFastaGenomeParser.lSeq
                self.lSeqs = [ (seq, len(iFastaGenomeParser.dSeq[seq])) for seq in lSeqNames ]
            else:
                logging.error('Missing Genome File - required for Circos plot')
                sys.exit(1)
        
            if self.TEFile:
                logging.debug('Parsing TE gff file')
                parser = GffTEParser(self.TEFile)
                self.lTEs = parser.getAllTEs()

            logging.info('Generating circos files')
            self.writeCircosPlotFiles()



    def writeCircosPlotFiles(self):
        """write circos files"""

        cPlot = CircosPlot()

        if self.lSeqs:
            SeqDataFile = 'genome.txt'
            logging.info('Writing circos sequence data file in {}'.format(SeqDataFile))
            cPlot.writeSeqDataFile(self.lSeqs, SeqDataFile)

        if self.lDuplications:
            SDDataFile = 'segdup.txt' 
            logging.info('Writing circos SD data file in {}'.format(SDDataFile))
            cPlot.writeSegDupDataFile(self.lDuplications, SDDataFile)
            SimilarityDataFile = 'similarity.txt'
            logging.info('Writing circos similarity data file in {}'.format(SimilarityDataFile))
            cPlot.writeSimilarityDataFile(self.lDuplications, SimilarityDataFile)

        if self.lGenes:
            GeneDataFile = 'gene.txt'
            logging.info('Writing circos gene data file in {}'.format(GeneDataFile))
            cPlot.writeGeneDataFile(self.lGenes, GeneDataFile)

        if self.lGeneLinks:
            GeneLinkDataFile = 'gene-link.txt'
            logging.info('Writing circos gene-link data file in {}'.format(GeneLinkDataFile))
            cPlot.writeGeneLinkDataFile(self.lGeneLinks, GeneLinkDataFile)
       
        if self.lTEs:
            TEDataFile = 'TE.txt'
            logging.info('Writing circos TE/Repeat data file in {}'.format(TEDataFile))
            cPlot.writeTEDataFile(self.lTEs, TEDataFile)

        CircosConfFile = 'circos.conf'
        logging.info('Writing circos configuration file in {}'.format(CircosConfFile))
        cPlot.writeCircosConf()



    def _buildGeneLinks(self,lGeneSeq1,lGeneSeq2,dup):
        """build"""

        lLinks = []
        for gene1 in lGeneSeq1:
            #print gene1.id
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
                   #print 'gene1 -  gene2 : {} {}'.format(gene1.id,gene2.id)
                   lLinks.append(GeneLink(dup=dup,gene1=gene1,gene2=gene2)) 
        return lLinks        
        # todo set : + logging.debug

      
    def _extractGeneInDuplication(self,dup):
        """extract """

        lGeneSeq1 = self.db.getlGenesFromCoordinates(dup.seq1,dup.start1,dup.end1)
        lGeneSeq2 = self.db.getlGenesFromCoordinates(dup.seq2,dup.start2,dup.end2)

        return (lGeneSeq1,lGeneSeq2)

if __name__ == "__main__":

    analyzer = Analyzer()
    analyzer.parseAttributesFromArgsCLI()
    analyzer.runAnalyze()
