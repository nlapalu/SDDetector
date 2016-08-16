#!/usr/bin/env python

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
#        GenomeFile = ('test-data/test2.fasta')
#        GenomeFile = ('test-data/test.fasta')
#        SDFile = ('test-data/sdd4.gff3')
#        SDFile = ('test-data/sdd3.gff3')
#        GeneFile = ('test-data/gene2.gff3')
#        GeneFile = ('test-data/gene.gff3')
#        TEFile = ('test-data/TE2.gff3')
#        TEFile = ('test-data/TE.gff3')
#        BlastXMLFile = ('test-data/blast4.xml')
#        BlastXMLFile = ('test-data/blast2.xml')
#        self.plot = CircosPlot(GenomeFile=GenomeFile,SDFile=SDFile,GeneFile=GeneFile,TEFile=TEFile,BlastXMLFile=BlastXMLFile)
         self.plot = CircosPlot()

    def tearDown(self):
        pass

    def test_writeCircosConf(self):
        """Test writeCircosConf"""
       
#        self.plot.writeSeqDataFile() 
#        self.plot.writeSegDupDataFile()
#        self.plot.writeGeneDataFile()
#        self.plot.writeTEDataFile() 
#        self.plot.writeCircosConf()

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


    def test_writeSeqDataFile(self):
        """Test writeSeqDataFile"""

        lSeqs = [('seq1',4),('seq2',8)] 
        GenomeDataFile = self.plot.writeSeqDataFile(lSeqs, 'genome.txt')
        self.assertTrue(filecmp.cmp(GenomeDataFile,'test-data/genome.txt'))

    def test_writeSegDupDataFile(self):
        """Test writeSegDupDataFile"""

        lDuplications = [Duplication('seq1',5,5000,'seq2',10,5000)]
        SDDataFile = self.plot.writeSegDupDataFile(lDuplications,'segdup.txt') 
        self.assertTrue(filecmp.cmp(SDDataFile,'test-data/segdup.txt'))


    def test_writeGeneDataFile(self):
        """Test writeGeneDataFile"""

        lGenes = [Gene('GENE1','seq1',12,600,1),Gene('GENE2','seq2',100,1000,-1)]
        GeneDataFile = self.plot.writeGeneDataFile(lGenes,'gene.txt') 
        self.assertTrue(filecmp.cmp(GeneDataFile,'test-data/gene.txt'))

    def test_writeGeneLinkDataFile(self):
        """Test writeGeneLinkDataFile"""

        iDup = Duplication('seq1',5,5000,'seq2',10,5000)
        iGene1 = Gene('GENE1','seq1',10,100,1)
        iGene2 = Gene('GENE2','seq2',100,190,-1)
        lGeneLinks = [GeneLink(dup=iDup,gene1=iGene1,gene2=iGene2)]
        GeneLinkDataFile = self.plot.writeGeneLinkDataFile(lGeneLinks,'gene-link.txt')
        self.assertTrue(filecmp.cmp(GeneLinkDataFile,'test-data/gene-link.txt'))

    def test_writeTEDataFile(self):
        """Test writeTEDataFile"""

        lTEs = [Feature('TE1','seq1',1000,2000,1,'TE'),Feature('TE2','seq2',4000,4500,-1,'TE')]
        TEDataFile = self.plot.writeTEDataFile(lTEs,'TE.txt')
        self.assertTrue(filecmp.cmp(TEDataFile, 'test-data/TE.txt'))


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


if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(TestCircosPlot)
    unittest.TextTestRunner(verbosity=2).run(suite)
