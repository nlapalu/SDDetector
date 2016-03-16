#!/usr/bin/env python

import unittest

from SDDetector.Db.AlignDB import AlignDB
from SDDetector.Entities.Alignment import Alignment
from SDDetector.Utils.AlignmentChainer import AlignmentChainer
from SDDetector.Entities.Chain import Chain

class TestAlignmentChainer(unittest.TestCase):

    def setUp(self):
        self.db = AlignDB(dbfile='test.db', logLevel='DEBUG')

    def tearDown(self):
        self.db.deleteDB()

    def test_distanceBetweenQueryAlgmts(self):
        """Test distanceBetweenQueryAlgmts"""

        algmt1 = Alignment('q1','s1',20,30,10,20,10,10,1,1,id=1)
        algmt2 = Alignment('q1','s1',100,110,50,60,10,10,1,1,id=2)
        algmt3 = Alignment('q1','s1',100,110,10,20,10,10,1,1,id=3)
        algmt4 = Alignment('q1','s1',20,30,50,60,10,10,1,1,id=4)
        algmt5 = Alignment('q1','s1',10,30,10,20,10,10,1,1,id=5)
        algmt6 = Alignment('q1','s1',20,50,50,60,10,10,1,1,id=6)
        algmt7 = Alignment('q1','s1',10,50,10,20,10,10,1,1,id=7)
        algmt8 = Alignment('q1','s1',20,30,50,60,10,10,1,1,id=8)
        algmt9 = Alignment('q1','s1',20,50,10,20,10,10,1,1,id=9)
        algmt10 = Alignment('q1','s1',10,30,50,60,10,10,1,1,id=10)
        algmt11 = Alignment('q1','s1',20,30,10,20,10,10,1,1,id=11)
        algmt12 = Alignment('q1','s1',10,50,50,60,10,10,1,1,id=12)

        iAlgmtChainer = AlignmentChainer(self.db)
        self.assertEquals(70,iAlgmtChainer.distanceBetweenQueryAlgmts(algmt1,algmt2), "1 before 2")
        self.assertEquals(70,iAlgmtChainer.distanceBetweenQueryAlgmts(algmt3,algmt4), "2 before 1")
        self.assertEquals(-10,iAlgmtChainer.distanceBetweenQueryAlgmts(algmt5,algmt6), "1 start , 2 overlap")
        self.assertEquals(-30,iAlgmtChainer.distanceBetweenQueryAlgmts(algmt7,algmt8), "1 start before, 2 nested")
        self.assertEquals(-10,iAlgmtChainer.distanceBetweenQueryAlgmts(algmt9,algmt10), "2 starrt before, 1 overlap")
        self.assertEquals(-30,iAlgmtChainer.distanceBetweenQueryAlgmts(algmt11,algmt12), "2 start before, 1 nested")
        

    def test_chainAlignments_Case1(self):
        """Test chainAlignment case1: simple chain"""
         
        al1 = Alignment('q1','s1',1000,1200,2000,2200,200,200,1,1, id=1) 
        al2 = Alignment('q1','s1',1500,2000,3000,3500,500,500,1,1, id=2)
        lAlignments = [al1, al2]
        self.db.deleteAllAlignments()
        self.db.insertlAlignments(lAlignments)
        iAlgmtChainer = AlignmentChainer(self.db)
        iAlgmtChainer.chainAlignments(lAlignments)
        self.assertEquals([Chain([al1,al2])],iAlgmtChainer.lChains)
        self.assertEquals({1:[0], 2:[0]},iAlgmtChainer.dIndex)

    def test_chainAlignments_Case2(self):
        """Test chainAlignment case2: variable constraints"""

        al1 = Alignment('q1','s1',1000,1200,2000,2200,200,200,1,1, id=1)
        al2 = Alignment('q1','s1',1500,2000,3000,3500,500,500,1,1, id=2)
        al3 = Alignment('q1','s1',10000,13000,150000,153000,3000,3000,1,-1,id=3)
        al4 = Alignment('q1','s1',123000,124000,160000,161000,1000,1000,1,1,id=4)
        lAlignments = [al1, al2, al3, al4]
        self.db.deleteAllAlignments()
        self.db.insertlAlignments(lAlignments)
        iAlgmtChainer = AlignmentChainer(self.db)
        iAlgmtChainer.chainAlignments(lAlignments)
        self.assertEquals([Chain([al1,al2]),Chain([al3]),Chain([al4])],iAlgmtChainer.lChains)
        self.assertEquals({1:[0], 2:[0], 3:[1], 4:[2]},iAlgmtChainer.dIndex)


    def test_chainAlignments_Case3(self):
        """Test chainAlignment case3: complex chain with interleaved coordinates"""

        al1 = Alignment('q1','s1',1000,1200,2000,2200,200,200,1,1, id=1)
        al2 = Alignment('q1','s1',1500,2000,3000,3500,500,500,1,1, id=2)
        al3 = Alignment('q1','s1',10000,13000,150000,153000,3000,3000,1,-1,id=3)
        al4 = Alignment('q1','s1',11000,14000,154000,157000,3000,3000,1,-1,id=4)
        al5 = Alignment('q1','s1',123000,124000,160000,161000,1000,1000,1,1,id=5)
        lAlignments = [al1, al2, al3, al4, al5]
        self.db.deleteAllAlignments()
        self.db.insertlAlignments(lAlignments)
        iAlgmtChainer = AlignmentChainer(self.db)
        iAlgmtChainer.chainAlignments(lAlignments)
        self.assertEquals([Chain([al1,al2]),Chain([al3,al4]),Chain([al5])],iAlgmtChainer.lChains)
        self.assertEquals({1:[0], 2:[0], 3:[1], 4:[1], 5:[2]},iAlgmtChainer.dIndex)


if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(TestAlignmentChainer)
    unittest.TextTestRunner(verbosity=2).run(suite)
