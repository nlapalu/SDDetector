#!/usr/bin/env python

import unittest

from SDDetector.Entities.Duplication import Duplication
from SDDetector.Entities.Region import Region
from SDDetector.Entities.GeneLink import GeneLink
from SDDetector.Entities.Gene import Gene
from SDDetector.Entities.Transcript import Transcript
from SDDetector.Entities.CDS import CDS

class TestGeneLink(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_getCDSAlignment(self):
        """test getCDSAlignment"""

        gene1 = Gene('G00001','Chr1',1,27,1,[Transcript('G00001.1','Chr1',1,27,1,'G00001',[CDS('G00001.1_cds_1','Chr1',1,6,1, 'G00001.1'),CDS('G00001.1_cds_1','Chr1',13,21,1,'G00001.1')])])
        gene2 = Gene('G00002','Chr5',1,27,1,[Transcript('G00002.1','Chr5',1,27,1,'G00002',[CDS('G00002.1_cds_1','Chr5',1,6,1, 'G00002.1'),CDS('G00002.1_cds_1','Chr5',13,27,1,'G00002.1')])])

        gl = GeneLink(Duplication('Chr1',1,58,'Chr5',1,60,[(Region('Chr1',1,58,1),Region('Chr5',1,60,1))],[('ATGTATTCTATCTCATGTTAATGCTAATACTAGTCATGATCAGATACGATGATGAT--TA','ATGTATTCTATCTCATGTTACTGCTAATACTAGTCATGATCAGATACGATGATGATCATA')]),gene1,gene2)

        self.assertEquals(('ATGTATtctatcTCATGTTAAtgctaa','ATGTATtctatcTCATGTTACTGCTAA',Region('Chr1',1,27,1),Region('Chr5',1,27,1)),gl.getCDSAlignment())

        gene1 = Gene('G00001','Chr1',1,27,1,[Transcript('G00001.1','Chr1',1,27,1,'G00001',[CDS('G00001.1_cds_1','Chr1',1,6,1, 'G00001.1'),CDS('G00001.1_cds_1','Chr1',13,21,1,'G00001.1')])])
        gene2 = Gene('G00002','Chr5',27,60,-1,[Transcript('G00002.1','Chr5',27,60,-1,'G00002',[CDS('G00002.1_cds_1','Chr5',27,39,-1, 'G00002.1'),CDS('G00002.1_cds_1','Chr5',48,60,-1,'G00002.1')])])

        gl = GeneLink(Duplication('Chr1',1,58,'Chr5',1,60,[(Region('Chr1',1,58,-1),Region('Chr5',1,60,1))],[('ATGTATTCTATCTCATGTTAATGCTAATACTAGTCATGATCAGATACGATGATGAG--TA','ATGTATTCTATCTCATGTTACTGCTAATACTAGTCATGATCAGATACGATGATGATCATA')]),gene1,gene2)

        self.assertEquals(('atactagtcatGATCAGATAcgatgaTGAG--TA','ATACTAGTCATGAtcagatacGATGATGATCATA',Region('Chr1',1,32,-1),Region('Chr5',27,60,1)),gl.getCDSAlignment())


if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(TestGeneLink)
    unittest.TextTestRunner(verbosity=2).run(suite)
