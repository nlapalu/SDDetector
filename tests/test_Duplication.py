#!/usr/bin/env python

import unittest

from SDDetector.Entities.Duplication import Duplication
from SDDetector.Entities.Region import Region

class TestDuplication(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_dSeqToSeq(self):
        """Duplication with alignment"""

        dSeqToSeq = {'seq1': {1: ('seq1', 1), 2: ('seq1', 2),
                              3: ('seq1', 3), 4: ('seq1', 4),
                              5: ('seq1', 5)}}

        dup = Duplication('seq1', 1, 5, 'seq1', 1, 5, [(Region('seq1', 1, 5, 1),
                          Region('seq1', 1, 5, 1))], [('ATGAT', 'ATGAT')]) 

        self.maxDiff = None
        self.assertEqual(dSeqToSeq,dup.dSeqToSeq)

        dSeqToSeq = {'seq2': {54: ('seq2', 302), 55: ('seq2', 301),
                              56: ('seq2', 300), 57: ('seq2', None),
                              58: ('seq2', 299), 59: ('seq2', 298),
                              60: ('seq2', 297), 61: ('seq2', 295),
                              62: ('seq2', 294), 63: ('seq2', 293),
                              64: ('seq2', None), 65: ('seq2', 292),
                              66: ('seq2', 291), 67: ('seq2', 290),
                              68: ('seq2', 289), 69: ('seq2', 288),
                              288: ('seq2', 69), 289: ('seq2', 68),
                              290: ('seq2', 67), 291: ('seq2', 66),
                              292: ('seq2', 65), 293: ('seq2', 63),
                              294: ('seq2', 62), 295: ('seq2', 61),
                              296: ('seq2', None), 297: ('seq2', 60),
                              298: ('seq2', 59), 299: ('seq2', 58),
                              300: ('seq2', 56), 301: ('seq2', 55),
                              302: ('seq2', 54)}}

        dup =  Duplication('seq2', 288, 302, 'seq2', 54, 69, [(Region('seq2',288, 302, -1),
                           Region('seq2', 54, 69, 1))],[('TGA-AGTCGCT-GTTTT', 'TGATAGT-GCTGGTTTT')])
 
        self.maxDiff = None
        self.assertEqual(dSeqToSeq,dup.dSeqToSeq)

        dSeqToSeq = {'seq1': {1: ('seq2', 241), 2: ('seq2', 240),
                              3: ('seq2', 239), 4: ('seq2', 238),
                              5: ('seq2', 237)},
                     'seq2': {237: ('seq1', 5), 238: ('seq1', 4),
                              239: ('seq1', 3), 240: ('seq1', 2),
                              241: ('seq1', 1)}}

        dup = Duplication('seq1', 1, 5, 'seq2', 237, 241, [(Region('seq1', 1, 5, -1),
                          Region('seq2', 237, 241, 1))],[('TCCTA', 'TCCTA')])

        self.maxDiff = None
        self.assertEqual(dSeqToSeq,dup.dSeqToSeq)

        dSeqToSeq = {'Chr1': {1: ('Chr5', 60), 2: ('Chr5', 59),
                              3: ('Chr5', 56), 4: ('Chr5', 55),
                              5: ('Chr5', 54), 6: ('Chr5', 53),
                              7: ('Chr5', 52), 8: ('Chr5', 51),
                              9: ('Chr5', 50), 10: ('Chr5', 49),
                              11: ('Chr5', 48), 12: ('Chr5', 47),
                              13: ('Chr5', 46), 14: ('Chr5', 45),
                              15: ('Chr5', 44), 16: ('Chr5', 43),
                              17: ('Chr5', 42), 18: ('Chr5', 41),
                              19: ('Chr5', 40), 20: ('Chr5', 39),
                              21: ('Chr5', 38), 22: ('Chr5', 37),
                              23: ('Chr5', 36), 24: ('Chr5', 35),
                              25: ('Chr5', 34), 26: ('Chr5', 33),
                              27: ('Chr5', 32), 28: ('Chr5', 31),
                              29: ('Chr5', 30), 30: ('Chr5', 29),
                              31: ('Chr5', 28), 32: ('Chr5', 27),
                              33: ('Chr5', 26), 34: ('Chr5', 25),
                              35: ('Chr5', 24), 36: ('Chr5', 23),
                              37: ('Chr5', 22), 38: ('Chr5', 21),
                              39: ('Chr5', 20), 40: ('Chr5', 19),
                              41: ('Chr5', 18), 42: ('Chr5', 17),
                              43: ('Chr5', 16), 44: ('Chr5', 15),
                              45: ('Chr5', 14), 46: ('Chr5', 13),
                              47: ('Chr5', 12), 48: ('Chr5', 11),
                              49: ('Chr5', 10), 50: ('Chr5', 9),
                              51: ('Chr5', 8), 52: ('Chr5', 7),
                              53: ('Chr5', 6), 54: ('Chr5', 5),
                              55: ('Chr5', 4), 56: ('Chr5', 3),
                              57: ('Chr5', 2), 58: ('Chr5', 1)},
                    'Chr5': {1: ('Chr1', 58), 2: ('Chr1', 57),
                              3: ('Chr1', 56), 4: ('Chr1', 55),
                              5: ('Chr1', 54), 6: ('Chr1', 53),
                              7: ('Chr1', 52), 8: ('Chr1', 51),
                              9: ('Chr1', 50), 10: ('Chr1', 49),
                              11: ('Chr1', 48), 12: ('Chr1', 47),
                              13: ('Chr1', 46), 14: ('Chr1', 45),
                              15: ('Chr1', 44), 16: ('Chr1', 43),
                              17: ('Chr1', 42), 18: ('Chr1', 41),
                              19: ('Chr1', 40), 20: ('Chr1', 39),
                              21: ('Chr1', 38), 22: ('Chr1', 37),
                              23: ('Chr1', 36), 24: ('Chr1', 35),
                              25: ('Chr1', 34), 26: ('Chr1', 33),
                              27: ('Chr1', 32), 28: ('Chr1', 31),
                              29: ('Chr1', 30), 30: ('Chr1', 29),
                              31: ('Chr1', 28), 32: ('Chr1', 27),
                              33: ('Chr1', 26), 34: ('Chr1', 25),
                              35: ('Chr1', 24), 36: ('Chr1', 23),
                              37: ('Chr1', 22), 38: ('Chr1', 21),
                              39: ('Chr1', 20), 40: ('Chr1', 19),
                              41: ('Chr1', 18), 42: ('Chr1', 17),
                              43: ('Chr1', 16), 44: ('Chr1', 15),
                              45: ('Chr1', 14), 46: ('Chr1', 13),
                              47: ('Chr1', 12), 48: ('Chr1', 11),
                              49: ('Chr1', 10), 50: ('Chr1', 9),
                              51: ('Chr1', 8), 52: ('Chr1', 7),
                              53: ('Chr1', 6), 54: ('Chr1', 5),
                              55: ('Chr1', 4), 56: ('Chr1', 3),
                              57: ('Chr1', None), 58: ('Chr1', None),
                              59: ('Chr1', 2), 60: ('Chr1', 1)}}

        dup = Duplication('Chr1' ,1 ,58 ,'Chr5' ,1 ,60 , [(Region('Chr1' ,1 ,58 ,-1),Region('Chr5',1,60,1))],[('ATGTATTCTATCTCATGTTAATGCTAATACTAGTCATGATCAGATACGATGATGAG--TA','ATGTATTCTATCTCATGTTACTGCTAATACTAGTCATGATCAGATACGATGATGATCATA')])
 

        self.maxDiff = None
        self.assertEqual(dSeqToSeq,dup.dSeqToSeq)

    def test_getSeqAlignment(self):
        """test"""

        dup = Duplication('seq1',1,164,'seq2',237,400,[(Region('seq1',1,164,1),Region('seq2',237,400,1))],[('TCCTAGCCGATCTAATCGGAGACTATCTAGGCGACATATAAAATTTTGACGAGCTATAATTATCCTCTAGGCTAATCGAGACCACCGCTATCCTACAGCCAAAACAGCGACTTCATCATCGAGCCGCACTATCCTCTACGGCGAGCGCTCTCTACGAGCATCAT','TCCTAGCCGATCTAATCGGAGACTATCTAGGCGACATATAAAATTTTGACGAGCTATAATTATCCTCTAGGCTAATCGAGACCACCGCTATCCTACAGCCAAAACAGCGACTTCATCATCGAGCCGCACTATCCTCTACGGCGAGCGCTCTCTACGAGCATCAT')])
        
        self.assertEquals('TCCT',dup.getSeqAlignment('seq1',1,4)[0])
 
        dup = Duplication('seq1',1,164,'seq2',237,400,[(Region('seq1',1,164,-1),Region('seq2',237,400,1))],[('TCCTAGCCGATCTAATCGGAGACTATCTAGGCGACATATAAAATTTTGACGAGCTATAATTATCCTCTAGGCTAATCGAGACCACCGCTATCCTACAGCCAAAACAGCGACTTCATCATCGAGCCGCACTATCCTCTACGGCGAGCGCTCTCTACGAGCATCAT','TCCTAGCCGATCTAATCGGAGACTATCTAGGCGACATATAAAATTTTGACGAGCTATAATTATCCTCTAGGCTAATCGAGACCACCGCTATCCTACAGCCAAAACAGCGACTTCATCATCGAGCCGCACTATCCTCTACGGCGAGCGCTCTCTACGAGCATCAT')])
        
        self.assertEquals('TCAT',dup.getSeqAlignment('seq1',1,4)[0]) 

        dup = Duplication('seq4',1,40,'seq5',20,30,[(Region('seq4',1,10,-1),Region('seq5',20,25,1)),(Region('seq4',35,40,-1),Region('seq5',25,30,1))],[('ATATATATAT','AT----ATAT'),('ATGT-TT','AT-TTTG')])

        self.assertEquals('ATGT',dup.getSeqAlignment('seq4',37,40)[0]) 
        self.assertEquals('T----ATAT',dup.getSeqAlignment('seq5',21,25)[0]) 


if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(TestDuplication)
    unittest.TextTestRunner(verbosity=2).run(suite)
