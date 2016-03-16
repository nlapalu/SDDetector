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
        self.assertEquals(69,iAlgmtChainer.distanceBetweenQueryAlgmts(algmt1,algmt2), "1 before 2")
        self.assertEquals(69,iAlgmtChainer.distanceBetweenQueryAlgmts(algmt3,algmt4), "2 before 1")
        self.assertEquals(-11,iAlgmtChainer.distanceBetweenQueryAlgmts(algmt5,algmt6), "1 start , 2 overlap")
        self.assertEquals(-31,iAlgmtChainer.distanceBetweenQueryAlgmts(algmt7,algmt8), "1 start before, 2 nested")
        self.assertEquals(-11,iAlgmtChainer.distanceBetweenQueryAlgmts(algmt9,algmt10), "2 starrt before, 1 overlap")
        self.assertEquals(-31,iAlgmtChainer.distanceBetweenQueryAlgmts(algmt11,algmt12), "2 start before, 1 nested")


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
        self.assertEquals([Chain([al1,al2]),Chain([al3]),Chain([al4]),Chain([al5])],iAlgmtChainer.lChains)
        self.assertEquals({1:[0], 2:[0], 3:[1], 4:[2], 5:[3]},iAlgmtChainer.dIndex)

    def test_removeOverlappingChains(self):
        """Test removeOverlappingChains"""

        al1 = Alignment('q1','s1',1000,1201,2000,2201,202,202,1,1, id=1)
        al2 = Alignment('q1','s1',1500,2001,3000,3501,502,502,1,1, id=2)
        al3 = Alignment('q1','s1',1000,1200,2100,2300,201,201,1,1,id=3)
        al4 = Alignment('q1','s1',1500,2000,3100,3600,501,501,1,1,id=4)
        al5 = Alignment('q1','s1',123000,124000,160000,161000,1001,1001,1,1,id=5)
        al1b = Alignment('s1','q1',2000,2201,1000,1201,202,202,1,1, id=16)
        al2b = Alignment('s1','q1',3000,3501,1500,2001,502,502,1,1, id=11)
        al3b = Alignment('s1','q1',2100,2300,1000,1200,201,201,1,1,id=31)
        al4b = Alignment('s1','q1',3100,3600,1500,2000,501,501,1,1,id=41)
        lAlignments = [al1, al2, al3, al4, al5,al1b, al2b, al3b, al4b]
        self.db.deleteAllAlignments()
        self.db.insertlAlignments(lAlignments)
        iAlgmtChainer = AlignmentChainer(self.db)
        iAlgmtChainer.chainAlignments(lAlignments)
        #for i in iAlgmtChainer.lChains:
        #    print i.convertChain(1)

    def test_removeOverlappingChains_case2(self):


        al1 = Alignment('Chr1','Chr1',15088003,15088582,15253139,15253732,594,594,1,1,id=1)
        al2 = Alignment('Chr1','Chr1',15094505,15096966,15254239,15256725,2499,2499,1,1,id=2)
        al3 = Alignment('Chr1','Chr1',15090989,15092375,15257818,15259247,1432,1432,1,1,id=3)
        al4 = Alignment('Chr1','Chr1',15098601,15102401,15259265,15263155,3904,3904,1,1,id=4)
        al5 = Alignment('Chr1','Chr1',15088003,15090006,15257754,15259751,2021,2021,1,1,id=5)
        al6 = Alignment('Chr1','Chr1',15090989,15092375,15257818,15259247,1432,1432,1,1,id=6)
        al7 = Alignment('Chr1','Chr1',15094655,15097013,15258796,15261191,2410,2410,1,1,id=7)
        al8 = Alignment('Chr1','Chr1',15098601,15102401,15259265,15263155,3904,3904,1,1,id=8)
        lAlignments = [al1, al2, al3, al4, al5, al6, al7, al8]
        self.db.deleteAllAlignments()
        self.db.insertlAlignments(lAlignments)
        iAlgmtChainer = AlignmentChainer(self.db)
        iAlgmtChainer.chainAlignments(lAlignments)
#        for i in iAlgmtChainer.lChains:
#            print i.convertChain(1)


#Chr1	SDDetector	match_part	15253139	15253732	.	+	.	ID=match90878;Parent=chain37;Target=Chr1 15088003 15088582;length=594;identities=556;identity_percentage=0.959
#Chr1	SDDetector	match_part	15254239	15256725	.	+	.	ID=match21649;Parent=chain37;Target=Chr1 15094505 15096966;length=2499;identities=2321;identity_percentage=0.943
#Chr1	SDDetector	match_part	15257818	15259247	.	+	.	ID=match50016;Parent=chain37;Target=Chr1 15090989 15092375;length=1432;identities=1293;identity_percentage=0.932
#Chr1	SDDetector	match_part	15259265	15263155	.	+	.	ID=match11790;Parent=chain37;Target=Chr1 15098601 15102401;length=3904;identities=3502;identity_percentage=0.921
#Chr1	SDDetector	match	15257754	15263155	.	.	.	ID=chain38;length=5402;nbIndels=268;nbSNPs=658;identity_percentage=0.926
#Chr1	SDDetector	match_part	15257754	15259751	.	+	.	ID=match29498;Parent=chain38;Target=Chr1 15088003 15090006;length=2021;identities=1882;identity_percentage=0.942
#Chr1	SDDetector	match_part	15257818	15259247	.	+	.	ID=match50016;Parent=chain38;Target=Chr1 15090989 15092375;length=1432;identities=1293;identity_percentage=0.932
#Chr1	SDDetector	match_part	15258796	15261191	.	+	.	ID=match27930;Parent=chain38;Target=Chr1 15094655 15097013;length=2410;identities=2164;identity_percentage=0.917
#Chr1	SDDetector	match_part	15259265	15263155	.	+	.	ID=match11790;Parent=chain38;Target=Chr1 15098601 15102401;length=3904;identities=3502;identity_percentage=0.921




#Chr1	SDDetector	match	15228832	15313131	.	.	.	ID=chain32;length=83144;nbIndels=171;nbSNPs=558;identity_percentage=0.993
#Chr1	SDDetector	match_part	15228832	15232676	.	+	.	ID=match9745;Parent=chain32;Target=Chr1 15322960 15326783;length=3867;identities=3590;identity_percentage=0.939
#Chr1	SDDetector	match_part	15229778	15248735	.	+	.	ID=match36;Parent=chain32;Target=Chr1 15326899 15345854;length=18958;identities=18956;identity_percentage=1.000
#Chr1	SDDetector	match_part	15239918	15240542	.	+	.	ID=match84887;Parent=chain32;Target=Chr1 15346528 15347152;length=625;identities=592;identity_percentage=0.947
#Chr1	SDDetector	match_part	15239966	15240542	.	+	.	ID=match105642;Parent=chain32;Target=Chr1 15321661 15322236;length=578;identities=520;identity_percentage=0.903
#Chr1	SDDetector	match_part	15249811	15298024	.	+	.	ID=match12;Parent=chain32;Target=Chr1 15347912 15396127;length=48223;identities=48170;identity_percentage=0.999
#Chr1	SDDetector	match_part	15253631	15253944	.	+	.	ID=match135366;Parent=chain32;Target=Chr1 15321254 15321567;length=315;identities=287;identity_percentage=0.914
#Chr1	SDDetector	match_part	15258216	15259751	.	+	.	ID=match51879;Parent=chain32;Target=Chr1 15318940 15320435;length=1548;identities=1364;identity_percentage=0.912
#Chr1	SDDetector	match_part	15259354	15259837	.	+	.	ID=match100030;Parent=chain32;Target=Chr1 15322413 15322898;length=487;identities=470;identity_percentage=0.971
#Chr1	SDDetector	match_part	15259359	15259751	.	+	.	ID=match130178;Parent=chain32;Target=Chr1 15320789 15321180;length=396;identities=355;identity_percentage=0.906
#Chr1	SDDetector	match_part	15298106	15312962	.	+	.	ID=match52;Parent=chain32;Target=Chr1 15396188 15411044;length=14857;identities=14857;identity_percentage=1.000
#Chr1	SDDetector	match_part	15310409	15313131	.	+	.	ID=match12094;Parent=chain32;Target=Chr1 15411105 15413833;length=2730;identities=2694;identity_percentage=0.989
#Chr1	SDDetector	match	15229778	15313131	.	.	.	ID=chain33;length=82198;nbIndels=146;nbSNPs=592;identity_percentage=0.993
#Chr1	SDDetector	match_part	15229778	15248735	.	+	.	ID=match36;Parent=chain33;Target=Chr1 15326899 15345854;length=18958;identities=18956;identity_percentage=1.000
#Chr1	SDDetector	match_part	15235255	15239076	.	+	.	ID=match9861;Parent=chain33;Target=Chr1 15322960 15326783;length=3843;identities=3557;identity_percentage=0.931
#Chr1	SDDetector	match_part	15239918	15240542	.	+	.	ID=match84887;Parent=chain33;Target=Chr1 15346528 15347152;length=625;identities=592;identity_percentage=0.947
#Chr1	SDDetector	match_part	15239966	15240542	.	+	.	ID=match105642;Parent=chain33;Target=Chr1 15321661 15322236;length=578;identities=520;identity_percentage=0.903
#Chr1	SDDetector	match_part	15249811	15298024	.	+	.	ID=match12;Parent=chain33;Target=Chr1 15347912 15396127;length=48223;identities=48170;identity_percentage=0.999
#Chr1	SDDetector	match_part	15253631	15253944	.	+	.	ID=match135366;Parent=chain33;Target=Chr1 15321254 15321567;length=315;identities=287;identity_percentage=0.914
#Chr1	SDDetector	match_part	15258216	15259751	.	+	.	ID=match51879;Parent=chain33;Target=Chr1 15318940 15320435;length=1548;identities=1364;identity_percentage=0.912
#Chr1	SDDetector	match_part	15259354	15259837	.	+	.	ID=match100030;Parent=chain33;Target=Chr1 15322413 15322898;length=487;identities=470;identity_percentage=0.971
#Chr1	SDDetector	match_part	15259359	15259751	.	+	.	ID=match130178;Parent=chain33;Target=Chr1 15320789 15321180;length=396;identities=355;identity_percentage=0.906
#Chr1	SDDetector	match_part	15298106	15312962	.	+	.	ID=match52;Parent=chain33;Target=Chr1 15396188 15411044;length=14857;identities=14857;identity_percentage=1.000
#Chr1	SDDetector	match_part	15310409	15313131	.	+	.	ID=match12094;Parent=chain33;Target=Chr1 15411105 15413833;length=2730;identities=2694;identity_percentage=0.989









    def test_removeInternalAlignments(self):
        """..."""

        al1 = Alignment('Chr1','Chr1',12806596,12809714,12796459,12799562,3123,3123,1,1,id=1)
        al2 = Alignment('Chr1','Chr1',12810507,12813088,12796459,12799028,2592,2592,1,1,id=2)
        lAlignments = [al1, al2]
        self.db.deleteAllAlignments()
        self.db.insertlAlignments(lAlignments)
        iAlgmtChainer = AlignmentChainer(self.db)
        iAlgmtChainer.chainAlignments(lAlignments)



if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(TestAlignmentChainer)
    unittest.TextTestRunner(verbosity=2).run(suite)
