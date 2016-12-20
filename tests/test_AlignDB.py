#!/usr/bin/env python

import unittest

from SDDetector.Db.AlignDB import AlignDB
from SDDetector.Entities.Alignment import Alignment

class TestAlignDB(unittest.TestCase):

    def setUp(self):
        self.db = AlignDB(dbfile='test.db', logLevel='DEBUG')

    def tearDown(self):
        self.db.deleteDB()

    def test_selectSelfSelfMatchAlgmts(self):
        """Test selectSelfSelfMatchAlgmts

              1            200
        s1    |==============|
              1            200
        s1    |==============|

        """

        lAlignments = [Alignment('s1','s1',1,200,1,200,200,200,1,1, id=1),
                       Alignment('q1','s1',750,800,125,175,50,50,1,1, id=2),
                       Alignment('q1','s1',3025,3100,160,235,75,75,1,1, id=3)]
        self.db.deleteAllAlignments()
        self.db.insertlAlignments(lAlignments)
        self.assertEquals([1],self.db.selectSelfSelfMatchAlgmts())

    def test_selectAlgmtsBelowIdentityThreshold(self):
        """Test selectAlgmtsBelowIdentityThreshold"""

        lAlignments = [Alignment('q1','s1',4004934,4006803,15489,13621,1879,1562,1,-1, id=1),
                       Alignment('q1','s1',1655049,1656926,15497,13621,1878,1874,1,-1, id=2)]
        self.db.deleteAllAlignments()
        self.db.insertlAlignments(lAlignments)
        self.assertEquals([1],self.db.selectAlgmtsBelowIdentityThreshold())

    def test_selectProximalAlgmts(self):
        """Test selectProximalAlgmts"""

        al1 = Alignment('q1','s1',400,600,300,500,200,200,1,1, id=1)
        al2 = Alignment('q1','s1',800,900,1000,1100,100,100,1,-1, id=2)
        al3 = Alignment('q1','s1',1000,1100,18000,19000,100,100,1,-1, id=3)
        al4 = Alignment('q1','s1',10000,13000,150000,153000,3000,3000,1,-1,id=4)
        al5 = Alignment('q1','s1',11000,14000,154000,157000,3000,3000,1,-1,id=5)
        al6 = Alignment('q1','s1',800,900,1000,1100,100,100,1,1, id=6)

        lAlignments = [al1, al2, al3, al4, al5, al6]
        self.db.deleteAllAlignments()
        self.db.insertlAlignments(lAlignments)

        self.assertEquals([al6],self.db.selectProximalAlgmts(1))
        self.assertEquals([],self.db.selectProximalAlgmts(1,maxGap=200))
        self.assertEquals([],self.db.selectProximalAlgmts(2))
        self.assertEquals([],self.db.selectProximalAlgmts(3))
        self.assertEquals([al3],self.db.selectProximalAlgmts(2,maxGap=18000))
        self.assertEquals([al5],self.db.selectProximalAlgmts(4))

    def test_selectSuboptimalAlgmts_case1(self):
        """Test selectSuboptimalAlgmts"""

        al1 = Alignment('q1','s1',1040,1120,40,120,80,80,1,1, id=1)
        al2 = Alignment('q1','s1',1000,1070,10,70,60,60,1,-1, id=2)
        al3 = Alignment('q1','s1',1060,1100,60,100,40,40,1,-1, id=3)
        al4 = Alignment('q1','s1',1110,1160,110,160,50,50,1,-1,id=4)
        al5 = Alignment('q1','s1',2000,2040,110,150,40,40,1,-1,id=5)
        lAlignments = [al1,al2,al3,al4,al5]
        self.db.deleteAllAlignments()
        self.db.insertlAlignments(lAlignments)

        self.assertEquals([2,3,4],self.db.selectSuboptimalAlgmts(0))

    def test_selectSuboptimalAlgmtsOverlap(self):
        """Test selectSuboptimalAlgmts with overlap"""

        al2 = Alignment('q1','s1',1000,1070,10,70,60,60,1,-1, id=2)
        al3 = Alignment('q1','s1',1060,1100,60,100,40,40,1,-1, id=3)
        lAlignments = [al2,al3]
        self.db.deleteAllAlignments()
        self.db.insertlAlignments(lAlignments)
        self.assertEquals([3],self.db.selectSuboptimalAlgmts(0))
        self.assertEquals([],self.db.selectSuboptimalAlgmts(11))
        self.assertEquals([3],self.db.selectSuboptimalAlgmts(10))


    def test_selectSuboptimalAlgmts(self):
        """Fix bug"""

        al1 = Alignment('Chr1','Chr1',8798183,8824158,8812562,8838548,25998,25998,1,1, id=1)
        al2 = Alignment('Chr1','Chr1',8798183, 8809780, 8826940, 8838548,11620,11620,1,1,id=2)
        lAlignments = [al1,al2]
        self.db.deleteAllAlignments()
        self.db.insertlAlignments(lAlignments)

        self.assertEquals([2],self.db.selectSuboptimalAlgmts(0))


#8798183	8824158	.	+	.	ID=match21;Parent=chain6;Target=Chr1 8812562 8838548;length=25998;identities=25903;identity_percentage=0.997
#8798183	8809780	.	+	.	ID=match281;Parent=chain7;Target=Chr1 8826940 8838548;length=11620;identities=11529;identity_percentage=0.994

    def test_selectQueryOnlySuboptimalAlgmts(self):

        al1 = Alignment('Chr1','Chr1',12796459,12799562,12806596,12809714,3123,3123,1,1,1)
        al2 = Alignment('Chr1','Chr1',12796459,12799028,12810507,12813088,2592,2592,1,1,2)
        lAlignments = [al1, al2]
        self.db.deleteAllAlignments()
        self.db.insertlAlignments(lAlignments)

        self.assertEquals([2],self.db.selectQueryOnlySuboptimalAlgmts([1,2]))

    def test_selectQueryOnlySuboptimalAlgmtsCase2(self):

        al1 = Alignment('Chr1','Chr1',8839058,8839291,8810281,8810514,237,237,1,1,1)
        al2 = Alignment('Chr1','Chr1',8839058,8839291,8824659,8824892,237,237,1,1,2)
        lAlignments = [al1, al2]
        self.db.deleteAllAlignments()
        self.db.insertlAlignments(lAlignments)

<<<<<<< HEAD
        self.assertEquals([1,2],sorted(self.db.selectQueryOnlySuboptimalAlgmts([1,2])))
=======
<<<<<<< HEAD
        self.assertEquals(sorted([2,1]),self.db.selectQueryOnlySuboptimalAlgmts([1,2]))
=======
        self.assertEquals([1,2],sorted(self.db.selectQueryOnlySuboptimalAlgmts([1,2])))
>>>>>>> eebe332c5cce0d1e6ee1a022995181b7e151de68
>>>>>>> 4e32a8fa5a027a1ab3853714e538aa98bd7adca3


    def test_selectAlignmentsWithDefinedIdOrderBySbjctCoord(self):
        """Test selectAlignmentsWithDefinedIdOrderBySbjctCoord"""

        al1 = Alignment('q1','s1',400,600,300,500,200,200,1,1, id=1)
        al2 = Alignment('q1','s1',800,900,1000,1100,100,100,1,-1, id=2)
        al3 = Alignment('q1','s1',1000,1100,18000,19000,100,100,1,-1, id=3)
        al4 = Alignment('q1','s1',10000,13000,150000,153000,3000,3000,1,-1,id=4)
        al5 = Alignment('q1','s1',11000,14000,154000,157000,3000,3000,1,-1,id=5)

        lAlignments = [al1,al2,al3,al4,al5]
        self.db.deleteAllAlignments()
        self.db.insertlAlignments(lAlignments)

        self.assertEquals([al1,al3,al4,al5],self.db.selectAlignmentsWithDefinedIdOrderBySbjctCoord([5,4,3,1]))

if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(TestAlignDB)
    unittest.TextTestRunner(verbosity=2).run(suite)
