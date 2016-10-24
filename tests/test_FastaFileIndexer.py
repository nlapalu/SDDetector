#!/usr/bin/env python

import os
import unittest

from SDDetector.Utils.FastaFileIndexer import FastaFileIndexer

class TestFastaFileIndexer(unittest.TestCase):

    def setUp(self):
        self.index = FastaFileIndexer(filename='test-data/seq.fasta', logLevel='DEBUG')

    def tearDown(self):
        pass

    def test_read(self):
        self.index.read()
        self.assertEquals(['seq1','seq2'],self.index.lSeq)
        self.assertEquals(167,len(self.index.dSeq['seq1']))
        self.assertEquals(400,len(self.index.dSeq['seq2']))

 
    def test_write(self):
        pass
#        self.index.read()
#        self.index.write()
#        self.index.delete()


if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(TestFastaFileIndexer)
    unittest.TextTestRunner(verbosity=2).run(suite)
