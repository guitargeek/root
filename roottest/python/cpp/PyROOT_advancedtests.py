# File: roottest/python/cpp/PyROOT_advancedtests.py
# Author: Wim Lavrijsen (LBNL, WLavrijsen@lbl.gov)
# Created: 06/04/05
# Last: 04/27/16

"""C++ advanced language interface unit tests for PyROOT package.

NOTE: most of the original test cases in this file were removed because they
are covered upstream in the cppyy test suite (bindings/pyroot/cppyy/cppyy/test/),
which descends from the same original PyROOT tests: see test_advancedcpp.py,
test_templates.py, test_lowlevel.py, test_operators.py, test_fragile.py,
test_regression.py and test_stltypes.py. What remains is only the global
variables test, whose upstream equivalent
(test_advancedcpp.py::test21_access_to_global_variables) is currently
marked xfail.
"""

import sys, os, unittest
sys.path.append(os.path.dirname(os.path.dirname(__file__)))

import ROOT
from ROOT import gROOT
from common import *

__all__ = [
   'Cpp08GlobalVariables',
]

gROOT.LoadMacro( "AdvancedCpp.C+" )


### Check access to global variables =========================================
# NOTE: kept although covered by cppyy's
# test_advancedcpp.py::test21_access_to_global_variables, because that
# upstream test is currently marked xfail. Delete once the xfail is resolved.
class Cpp08GlobalVariables( MyTestCase ):
   def test1DoubleArray( self ):
      """Verify access to array of doubles"""

      self.assertEqual( ROOT.myGlobalDouble, 12. )
      self.assertRaises( IndexError, ROOT.myGlobalArray.__getitem__, 500 )

   def test2WriteGlobalInstances( self ):
      """Verify writability of global instances"""

      NS_PR_Lumi = ROOT.NS_PR_Lumi

      def verify( func, name, val ):
         self.assertEqual( func(),                    val )
         self.assertEqual( getattr( ROOT, name ),     val )

      verify( ROOT.PR_GetLumi1, "PR_Lumi1", "::1 C++ global lumi" )

      ROOT.PR_Lumi1 = "::1 python global lumi"
      verify( ROOT.PR_GetLumi1, "PR_Lumi1", "::1 python global lumi" )

      ROOT.PR_Lumi2 = "::2 python global lumi"
      verify( ROOT.PR_GetLumi2, "PR_Lumi2", "::2 python global lumi" )

      def verify( func, name, val ):
         self.assertEqual( func(),                      val )
         self.assertEqual( getattr( NS_PR_Lumi, name ), val )

      verify( NS_PR_Lumi.PR_GetLumi1, "PR_Lumi1", "NS::1 C++ global lumi" )

      NS_PR_Lumi.PR_Lumi1 = "NS::1 python global lumi"
      verify( NS_PR_Lumi.PR_GetLumi1, "PR_Lumi1", "NS::1 python global lumi" )

      NS_PR_Lumi.PR_Lumi2 = "NS::2 python global lumi"
      verify( NS_PR_Lumi.PR_GetLumi2, "PR_Lumi2", "NS::2 python global lumi" )


## actual test run
if __name__ == '__main__':
   from MyTextTestRunner import MyTextTestRunner

   loader = unittest.TestLoader()
   testSuite = loader.loadTestsFromModule( sys.modules[ __name__ ] )

   runner = MyTextTestRunner( verbosity = 2 )
   result = not runner.run( testSuite ).wasSuccessful()

   sys.exit( result )
