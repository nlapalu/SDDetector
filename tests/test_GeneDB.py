#!/usr/bin/env python

import unittest

from SDDetector.Db.GeneDB import GeneDB

from SDDetector.Entities.Gene import Gene
from SDDetector.Entities.Transcript import Transcript
from SDDetector.Entities.CDS import CDS

class TestGeneDB(unittest.TestCase):

    def setUp(self):
        self.db = GeneDB(dbfile='test.db', logLevel='DEBUG')

    def tearDown(self):
        self.db.deleteDB()

    def test_selectAllGenes(self):
        """Test selectAllGenes"""

        gene1 = Gene('G00001','Chr1',23988,24919,-1,[Transcript('G00001.1','Chr1',23988,24919,-1,'G00001',[CDS('G00001.1_cds_1','Chr1',23988,24083, -1, 'G00001.1'),CDS('G00001.1_cds_1','Chr1',24274,24427,-1,'G00001.1'),CDS('G00001.1_cds_1','Chr1',24489,24919,-1,'G00001.1')])])

        lGenes = [gene1]
        # self.db.deleteAllGenes()  TODO
        self.db.insertlGenes(lGenes)

        self.assertEquals([gene1],self.db.selectAllGenes())

    def test_getlGenesFromCoordinates(self):
        """Test getlGenesFromCoordinates"""

        gene1 = Gene('G00001','Chr1',23988,24919,-1,[Transcript('G00001.1','Chr1',23988,24919,-1,'G00001',[CDS('G00001.1_cds_1','Chr1',23988,24083, -1, 'G00001.1'),CDS('G00001.1_cds_1','Chr1',24274,24427,-1,'G00001.1'),CDS('G00001.1_cds_1','Chr1',24489,24919,-1,'G00001.1')])])
        gene2 = Gene('G00002','Chr1',239880,249190,-1,[Transcript('G00002.1','Chr1',239880,249190,-1,'G00002',[CDS('G00002.1_cds_1','Chr1',239880,240830, -1, 'G00002.1'),CDS('G00002.1_cds_1','Chr1',242740,244270,-1,'G00002.1'),CDS('G00002.1_cds_1','Chr1',244890,249190,-1,'G00002.1')])])

        lGenes = [gene1,gene2]
        # self.db.deleteAllGenes()  TODO
        self.db.insertlGenes(lGenes)

        self.assertEquals([gene2],self.db.getlGenesFromCoordinates('Chr1',230000,250000))


if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(TestGeneDB)
    unittest.TextTestRunner(verbosity=2).run(suite)
