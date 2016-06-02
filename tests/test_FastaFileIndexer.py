#!/usr/bin/env python

import unittest

from SDDetector.Utils.FastaFileIndexer import FastaFileIndexer

class TestFastaFileIndexer(unittest.TestCase):

    def setUp(self):
        self.index = FastaFileIndexer(filename='test-data/test.fasta', logLevel='DEBUG')

    def tearDown(self):
        self.index.delete()

    def test_read(self):
        #self.index.read()
        pass
 
    def test_write(self):
        self.index.read()
        self.index.write()


if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(TestFastaFileIndexer)
    unittest.TextTestRunner(verbosity=2).run(suite)
