#!/usr/bin/env python

import os
import unittest
import filecmp

from SDDetector.Entities.Duplication import Duplication
from SDDetector.Entities.Region import Region
from SDDetector.Entities.Gene import Gene
from SDDetector.Entities.GeneLink import GeneLink
from SDDetector.Entities.Feature import Feature
from SDDetector.Utils.CircosPlot import CircosPlot

class TestCircosPlot(unittest.TestCase):

    def setUp(self):
         self.plot = CircosPlot()

    def tearDown(self):
        pass        

    def test_writeCircosConf(self):
        """Test writeCircosConf"""
       
        lSeqs = [('seq1',20000),('seq2',30000)] 
        GenomeDataFile = self.plot.writeSeqDataFile(lSeqs, 'genome.txt')
        lRegions = [(Region('seq1',100,220,1),Region('seq2',100,220,1)),
                    (Region('seq1',1200,1300,-1),Region('seq2',1300,1400,1))] 
        lAlgmts = [('ATGCATGCATGCATGCATGCATGCATGCATGCATGCAGGCATGCATGCATGCATGCATGAATGCATGCATGCATGCATGCATGCATGCATGCATGCATGAATATGTGTAGTGAGTCGTCCC',
                    'ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGAATGCATGCATGCATGCATGCATGCATGCATGCATGCATGAATGATGTACGATATAGCCCAC'),
                   ('ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGAATGCATGCATGCATGCATGCATGCATGCATGCATGCATGAA',
                    'ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGAATGCATGCATGCATGCATGCATGCATGCATGCATGCATGAA')]
        lDuplications = [Duplication('seq1',1,5000,'seq2',1,6000,lRegions,lAlgmts)]
        SDDataFile = self.plot.writeSegDupDataFile(lDuplications,'segdup.txt') 
        self.plot.writeCircosConf()
        self.assertTrue(filecmp.cmp('circos.conf','test-data/circos.conf'))
        lGenes = [Gene('GENE1','seq1',12,600,1),Gene('GENE2','seq2',100,1000,-1)]
        GeneDataFile = self.plot.writeGeneDataFile(lGenes,'gene.txt') 
        lTEs = [Feature('TE1','seq1',1000,2000,1,'TE'),Feature('TE2','seq2',4000,4500,-1,'TE')]
        TEDataFile = self.plot.writeTEDataFile(lTEs,'TE.txt')
        self.plot.writeCircosConf()
        self.assertTrue(filecmp.cmp('circos.conf','test-data/circos2.conf'))
        SimilarityDataFile = self.plot.writeSimilarityDataFile(lDuplications,'similarity.txt')
        self.plot.writeCircosConf()
        self.assertTrue(filecmp.cmp('circos.conf','test-data/circos3.conf'))
        os.remove('circos.conf')

    def test_writeSeqDataFile(self):
        """Test writeSeqDataFile"""

        lSeqs = [('seq1',4),('seq2',8)] 
        GenomeDataFile = self.plot.writeSeqDataFile(lSeqs, 'genome.txt')
        self.assertTrue(filecmp.cmp(GenomeDataFile,'test-data/genome.txt'))
        os.remove('genome.txt')

    def test_writeSegDupDataFile(self):
        """Test writeSegDupDataFile"""

        lDuplications = [Duplication('seq1',5,5000,'seq2',10,5000)]
        SDDataFile = self.plot.writeSegDupDataFile(lDuplications,'segdup.txt') 
        self.assertTrue(filecmp.cmp(SDDataFile,'test-data/segdup.txt'))
        os.remove('segdup.txt')

    def test_writeGeneDataFile(self):
        """Test writeGeneDataFile"""

        lGenes = [Gene('GENE1','seq1',12,600,1),Gene('GENE2','seq2',100,1000,-1)]
        GeneDataFile = self.plot.writeGeneDataFile(lGenes,'gene.txt') 
        self.assertTrue(filecmp.cmp(GeneDataFile,'test-data/gene.txt'))
        os.remove('gene.txt')

    def test_writeGeneLinkDataFile(self):
        """Test writeGeneLinkDataFile"""

        iDup = Duplication('seq1',5,5000,'seq2',10,5000)
        iGene1 = Gene('GENE1','seq1',10,100,1)
        iGene2 = Gene('GENE2','seq2',100,190,-1)
        lGeneLinks = [GeneLink(dup=iDup,gene1=iGene1,gene2=iGene2)]
        GeneLinkDataFile = self.plot.writeGeneLinkDataFile(lGeneLinks,'gene-link.txt')
        self.assertTrue(filecmp.cmp(GeneLinkDataFile,'test-data/gene-link.txt'))
        os.remove('gene-link.txt')

    def test_writeTEDataFile(self):
        """Test writeTEDataFile"""

        lTEs = [Feature('TE1','seq1',1000,2000,1,'TE'),Feature('TE2','seq2',4000,4500,-1,'TE')]
        TEDataFile = self.plot.writeTEDataFile(lTEs,'TE.txt')
        self.assertTrue(filecmp.cmp(TEDataFile, 'test-data/TE.txt'))
        os.remove('TE.txt')

    def test_writeSimilarityDataFile(self):
        """Test writeSimilarityDataFile"""

        lRegions = [(Region('seq1',100,220,1),Region('seq2',100,220,1)),
                    (Region('seq1',1200,1300,-1),Region('seq2',1300,1400,1))] 
        lAlgmts = [('ATGCATGCATGCATGCATGCATGCATGCATGCATGCAGGCATGCATGCATGCATGCATGAATGCATGCATGCATGCATGCATGCATGCATGCATGCATGAATATGTGTAGTGAGTCGTCCC',
                    'ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGAATGCATGCATGCATGCATGCATGCATGCATGCATGCATGAATGATGTACGATATAGCCCAC'),
                   ('ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGAATGCATGCATGCATGCATGCATGCATGCATGCATGCATGAA',
                    'ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGAATGCATGCATGCATGCATGCATGCATGCATGCATGCATGAA')]
        lDuplications = [Duplication('seq1',1,5000,'seq2',1,6000,lRegions,lAlgmts)]
        SimilarityDataFile = self.plot.writeSimilarityDataFile(lDuplications,'similarity.txt')
        self.assertTrue(filecmp.cmp(SimilarityDataFile,'test-data/similarity.txt'))
        os.remove('similarity.txt')


if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(TestCircosPlot)
    unittest.TextTestRunner(verbosity=2).run(suite)
