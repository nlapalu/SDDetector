import unittest

from distutils.cmd import Command
from distutils.core import setup
from unittest import TextTestRunner, TestLoader, TestSuite

from tests.test_AlignDB import TestAlignDB
from tests.test_AlignmentChainer import TestAlignmentChainer
from tests.test_BlastTabParser import TestBlastTabParser

try:
    from tests.test_BlastXMLParser import TestBlastXMLParser
    Biopython_available = True
except ImportError:
    Biopython_available = False

class TestSuite(Command):
    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):

        lTestCases = [TestAlignDB,TestAlignmentChainer, TestBlastTabParser]
        if Biopython_available:
            lTestCases.append(TestBlastXMLParser)
        for case in lTestCases:
            suite = unittest.TestLoader().loadTestsFromTestCase(case)
            t = TextTestRunner(verbosity = 2)
            t.run(suite)

execfile('SDDetector/version.py'). 

setup(name='SDDetector',
      version = __version__,
      description='Segmental Duplication Detection',
      url='http://github.com/nlapalu/sddetector',
      author='Nicolas Lapalu',
      author_email='nicolas.lapalu@versailles.inra.fr',
      license='GPL V3',
      scripts=['bin/segmental_duplication_detector.py'],
      packages=['SDDetector.Db','SDDetector.Utils','SDDetector.Parser',
                'SDDetector.Parser.Blast', 'SDDetector.Entities','SDDetector', 'tests'],
      data_files=[('test-data', ['test-data/blast.xml','test-data/blast.tab']),
                  ('.',['README.md'])],
      cmdclass={
          'test': TestSuite
      }
      )
