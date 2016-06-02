#!/usr/bin/env python

import unittest
import filecmp

from SDDetector.Utils.CircosPlot import CircosPlot

class TestCircosPlot(unittest.TestCase):

    def setUp(self):
        GenomeFile = ('test-data/test2.fasta')
        SDFile = ('test-data/sdd2.gff3')
        GeneFile = ('test-data/gene.gff3')
        self.plot = CircosPlot(GenomeFile=GenomeFile,SDFile=SDFile,GeneFile=GeneFile)

    def tearDown(self):
        pass

    def test_writeCircosConf(self):
        """Test writeCircosConf"""
       
        self.plot.writeSeqDataFile() 
        self.plot.writeSegDupDataFile()
        self.plot.writeGeneDataFile() 
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

if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(TestCircosPlot)
    unittest.TextTestRunner(verbosity=2).run(suite)
