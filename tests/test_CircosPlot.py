#!/usr/bin/env python

import unittest
import filecmp

from SDDetector.Utils.CircosPlot import CircosPlot

class TestCircosPlot(unittest.TestCase):

    def setUp(self):
        GenomeFile = ('test-data/test2.fasta')
#        SDFile = ('test-data/sdd4.gff3')
        SDFile = ('test-data/sdd3.gff3')
        GeneFile = ('test-data/gene2.gff3')
        TEFile = ('test-data/TE2.gff3')
#        BlastXMLFile = ('test-data/blast4.xml')
        BlastXMLFile = ('test-data/blast2.xml')
        self.plot = CircosPlot(GenomeFile=GenomeFile,SDFile=SDFile,GeneFile=GeneFile,TEFile=TEFile,BlastXMLFile=BlastXMLFile)

    def tearDown(self):
        pass

    def test_writeCircosConf(self):
        """Test writeCircosConf"""
       
        self.plot.writeSeqDataFile() 
        self.plot.writeSegDupDataFile()
        self.plot.writeGeneDataFile()
        self.plot.writeTEDataFile() 
        self.plot.writeCircosConf()

    def test_writeSeqDataFile(self):
        """Test writeSeqDataFile"""

        GenomeDataFile = self.plot.writeSeqDataFile()
        self.assertTrue(filecmp.cmp(GenomeDataFile,'test-data/genome.txt'))

    def test_writeSegDupDataFile(self):
        """Test writeSegDupDataFile"""

        SDDataFile = self.plot.writeSegDupDataFile() 
        self.assertTrue(filecmp.cmp(SDDataFile,'test-data/segdup.txt'))


    def test_writeGeneDataFile(self):
        """Test writeGeneDataFile"""

        GeneDataFile = self.plot.writeGeneDataFile() 
        self.assertTrue(filecmp.cmp(GeneDataFile,'test-data/gene.txt'))

    def test_writeGeneLinkDataFile(self):
        """Test writeGeneLinkDataFile"""

        GeneLinkDataFile = self.plot.writeGeneLinkDataFile()
        self.assertTrue(filecmp.cmp(GeneLinkDataFile,'test-data/gene-link.txt'))

    def test_writeTEDataFile(self):
        """Test writeTEDataFile"""

        TEDataFile = self.plot.writeTEDataFile()
        self.assertTrue(filecmp.cmp(TEDataFile, 'test-data/TE.txt'))


    def test_writeSimilarityDataFile(self):
        """Test writeSimilarityDataFile"""

        #SimilarityDataFile = self.plot.writeSimilarityDataFile()
        #self.assertTrue(filecmp.cmp(SimilarityDataFile,'test-data/similarity.txt'))
        pass


if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(TestCircosPlot)
    unittest.TextTestRunner(verbosity=2).run(suite)
