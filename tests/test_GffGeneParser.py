#!/usr/bin/env python


import unittest

from SDDetector.Entities.Gene import Gene
from SDDetector.Entities.Transcript import Transcript
from SDDetector.Entities.CDS import CDS
from SDDetector.Parser.Gff.GffGeneParser import GffGeneParser

class TestGffGeneParser(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_getAllGenes(self):
        """Test getAllGenes method"""

        iGffGeneParser = GffGeneParser("test-data/gene.gff3")
        lGenes = [Gene('G00001','Chr1',23988,24919,-1,[Transcript('G00001.1','Chr1',23988,24919,-1,'G00001',[CDS('G00001.1_cds_1','Chr1',23988,24083, -1, 'G00001.1'),CDS('G00001.1_cds_1','Chr1',24274,24427,-1,'G00001.1'),CDS('G00001.1_cds_1','Chr1',24489,24919,-1,'G00001.1')])])]

        self.assertEqual(iGffGeneParser.getAllGenes()[0],lGenes[0])

if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(TestGffGeneParser)
    unittest.TextTestRunner(verbosity=2).run(suite)
