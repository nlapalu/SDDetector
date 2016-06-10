#!/usr/bin/env python

import unittest

from SDDetector.Entities.Feature import Feature
from SDDetector.Parser.Gff.GffTEParser import GffTEParser

class TestGffTEParser(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_getAllTEs(self):
        """Test getAllTEs method"""

        iGffTEParser = GffTEParser("test-data/TE.gff3")
        lTEs = [Feature('TE','mp1029-1_seq2_RLX-chim_CollD-B-R14-Map20','seq2',17438,17470,-1),Feature('TE','mp1029-2_seq2_RLX-chim_CollD-B-R14-Map20','seq2',17438,17783,-1)]

        self.assertEqual(iGffTEParser.getAllTEs(),lTEs)


if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(TestGffTEParser)
    unittest.TextTestRunner(verbosity=2).run(suite)
