#!/usr/bin/env python

import unittest

from SDDetector.Utils.EffectPredictor import EffectPredictor
from SDDetector.Entities.Gene import Gene
from SDDetector.Entities.Transcript import Transcript
from SDDetector.Entities.CDS import CDS

class TestEffectPredictor(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_basicStructureAnalysis(self):
        """ test """

        EP = EffectPredictor()
        g1 = Gene('G00001','Chr1',23988,24919,-1,[Transcript('G00001.1','Chr1',23988,24919,-1,'G00001',[CDS('G00001.1_cds_1','Chr1',23988,24083, -1, 'G00001.1'),CDS('G00001.1_cds_1','Chr1',24274,24427,-1,'G00001.1'),CDS('G00001.1_cds_1','Chr1',24489,24919,-1,'G00001.1')])])
        g2 = Gene('G00002','Chr1',23988,24919,-1,[Transcript('G00002.1','Chr1',23988,24919,-1,'G00002',[CDS('G00002.1_cds_1','Chr1',23988,24083, -1, 'G00002.1'),CDS('G00002.1_cds_1','Chr1',24274,24427,-1,'G00002.1'),CDS('G00002.1_cds_1','Chr1',24489,24919,-1,'G00002.1')])])
        g3 = Gene('G00002','Chr1',23988,24920,-1,[Transcript('G00002.1','Chr1',23988,24920,-1,'G00002',[CDS('G00002.1_cds_1','Chr1',23988,24083, -1, 'G00002.1'),CDS('G00002.1_cds_1','Chr1',24274,24427,-1,'G00002.1'),CDS('G00002.1_cds_1','Chr1',24489,24920,-1,'G00002.1')])])
        g4 = Gene('G00002','Chr1',23988,24427,-1,[Transcript('G00002.1','Chr1',23988,24427,-1,'G00002',[CDS('G00002.1_cds_1','Chr1',23988,24083, -1, 'G00002.1'),CDS('G00002.1_cds_1','Chr1',24274,24427,-1,'G00002.1')])])
        self.assertFalse(EP.basicStructureAnalysis(g1,g2))
        self.assertEquals(EP.basicStructureAnalysis(g1,g3),'SIZE')
        self.assertEquals(EP.basicStructureAnalysis(g1,g4),'EXON')


    def test_extractVariantFromAlignment(self):
        """ test """

        EP = EffectPredictor()

        seq1 = 'ATGGCAGCAGCAAGCTAA'
        seq2 = 'ATGGCAGCAGCAATGCAA'

        self.assertEquals((3,0),EP.extractVariantFromAlignment(seq1,seq2))
        
        seq3 = 'ATG..GCAGC.AGCAAAAAGCTAA'
        seq4 = 'ATGGGGCAGCTAGC...AATGCAA'

        self.assertEquals((3,3),EP.extractVariantFromAlignment(seq3,seq4))


if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(TestEffectPredictor)
    unittest.TextTestRunner(verbosity=2).run(suite)
