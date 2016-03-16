#!/usr/bin/env python

import unittest

from SDDetector.Entities.Duplication import Duplication
from SDDetector.Entities.Region import Region
from SDDetector.Parser.Gff.GffDuplicationParser import GffDuplicationParser

class TestGffDuplicationParser(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_getAllDuplications(self):
        """Test getAllDuplication method"""

        iGffDuplicationParser = GffDuplicationParser("test-data/sdd.gff3")
        lDuplications = [Duplication('seq1',6010863,6029759,'seq8',4391356,4410272,[(Region('seq1',6010863,6029759,1),Region('seq8',4391356,4410272,1))]),
                         Duplication('seq2',26727,32020,'seq11',521201,524615,[(Region('seq2',26727,29266,-1),Region('seq11',522092,524615,1)),(Region('seq2',31119,32020,-1),Region('seq11',521201,522101,1))]),
                         Duplication('seq8',4391356,4410272,'seq1',6010863,6029759,[(Region('seq8',4391356,4410272,1),Region('seq1',6010863,6029759,1))]),
                         Duplication('seq11',26582,33594,'seq11',584193,591205,[(Region('seq11',26582,33594,-1),Region('seq11',584193,591205,1))]),
                         Duplication('seq11',38277,40563,'seq11',554466,556516,[(Region('seq11',38277,38402,1),Region('seq11',554466,554591,1)),(Region('seq11',38511,40563,-1),Region('seq11',554467,556516,1))]),
                         Duplication('seq11',521201,524615,'seq2',26727,32020,[(Region('seq11',521201,522101,-1),Region('seq2',31119,32020,1)),(Region('seq11',522092,524615,-1),Region('seq2',26727,29266,1))]),
                         Duplication('seq11',554466,556516,'seq11',38277,40563,[(Region('seq11',554466,554591,1),Region('seq11',38277,38402,1)),(Region('seq11',554467,556516,-1),Region('seq11',38511,40563,1))]),
                         Duplication('seq11',584193,591205,'seq11',26582,33594,[(Region('seq11',584193,591205,-1),Region('seq11',26582,33594,1))])]


        self.assertEqual(iGffDuplicationParser.getAllDuplications(),lDuplications)

    def test_getNonRedondantDuplications(self):
        """Test getNonRedondantDuplications method"""

        iGffDuplicationParser = GffDuplicationParser("test-data/sdd.gff3")
        lDuplications = [Duplication('seq1',6010863,6029759,'seq8',4391356,4410272,[(Region('seq1',6010863,6029759,1),Region('seq8',4391356,4410272,1))]),
                         Duplication('seq2',26727,32020,'seq11',521201,524615,[(Region('seq2',26727,29266,-1),Region('seq11',522092,524615,1)),(Region('seq2',31119,32020,-1),Region('seq11',521201,522101,1))]),
                         Duplication('seq11',26582,33594,'seq11',584193,591205,[(Region('seq11',26582,33594,-1),Region('seq11',584193,591205,1))]),
                         Duplication('seq11',38277,40563,'seq11',554466,556516,[(Region('seq11',38277,38402,1),Region('seq11',554466,554591,1)),(Region('seq11',38511,40563,-1),Region('seq11',554467,556516,1))])]

        self.assertEqual(iGffDuplicationParser.getNonRedondantDuplications(),lDuplications)


if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(TestGffDuplicationParser)
    unittest.TextTestRunner(verbosity=2).run(suite)
