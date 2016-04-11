#!/usr/bin/env python

import sqlite3
import unittest

from SDDetector.Entities.Alignment import Alignment
from SDDetector.Parser.Blast.BlastTabParser import BlastTabParser

class TestBlastTabParser(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_getAllAlignments(self):
        """Test getAllAlignments method"""

        iBlastTabParser = BlastTabParser("test-data/blast.tab")
        lAlignments = [Alignment('seq1','seq1',1,167,1,167,167,167,1,1,id=1),
                       Alignment('seq1','seq2',1,164,237,400,164,164,1,-1,id=2),
                       Alignment('seq1','seq2',50,113,54,118,66,63,1,1,id=3),
                       Alignment('seq2','seq2',1,400,1,400,400,400,1,1,id=4),
                       Alignment('seq2','seq2',288,351,54,118,66,63,1,-1,id=5),
                       Alignment('seq2','seq2',54,118,288,351,66,63,1,-1,id=6),
                       Alignment('seq2','seq1',237,400,1,164,164,164,1,-1,id=7),
                       Alignment('seq2','seq1',54,118,50,113,66,63,1,1,id=8)]

        self.assertEqual(iBlastTabParser.getAllAlignments(),lAlignments)


if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(TestBlastTabParser)
    unittest.TextTestRunner(verbosity=2).run(suite)
