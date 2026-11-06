# File: roottest/python/regression/PyROOT_regressiontests.py
# Author: Wim Lavrijsen (LBNL, WLavrijsen@lbl.gov)
# Created: 01/02/07
# Last: 04/26/16

"""Regression tests, lacking a better place, for PyROOT package.

NOTE: several of the original test cases in this file were removed because
they are covered upstream in the cppyy test suite
(bindings/pyroot/cppyy/cppyy/test/): see test_advancedcpp.py,
test_templates.py, test_datatypes.py, test_fragile.py, test_streams.py,
test_regression.py and test_crossinheritance.py. What remains are
ROOT-specific tests and tests without upstream equivalents.
"""

import platform
import sys, os, unittest
sys.path.append(os.path.dirname(os.path.dirname(__file__)))

if not os.path.exists('Amir.py'):
    os.chdir(os.path.dirname(__file__))

try:
   import commands
   WEXITSTATUS = os.WEXITSTATUS
except ImportError:
   import subprocess as commands
   def WEXITSTATUS(arg): return arg

original_preload = os.environ.get('LD_PRELOAD', None)

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = False
from ROOT import gROOT
from ROOT import TClass, TFile
from ROOT import TVector3, TGraph, TMatrixD
import cppyy

cleaned_preload = os.environ.get('LD_PRELOAD', None)

from common import *

__all__ = [
   'Regression01TwiceImportStar',
   'Regression03OldCrashers',
   'Regression04Threading',
   'Regression09TVector3Pythonize',
   'Regression11GlobalsLookup',
   'Regression12WriteTGraph',
   'Regression17MatrixD',
]


### "from ROOT import *" done in import-*-ed module ==========================
from Amir import *


class Regression01TwiceImportStar( MyTestCase ):
   def test1FromROOTImportStarInModule( self ):
      """Test handling of twice 'from ROOT import*'"""

      x = TestTChain()        # TestTChain defined in Amir.py


# NOTE: the former Regression02PyException was migrated to the cppjit test
# suite: test_regression.py::test63_python_exception_from_cpp_function.


### Several tests that used to cause crashes =================================
# NOTE: the former tests 1, 3, 4 and 5 of this class were migrated to the
# cppjit test suite: test_advancedcpp.py::test34_temporary_with_custom_operator_new,
# test_pythonization.py::test11_namespaced_class_with_iterators,
# test_operators.py::test18_dereference_operator_returning_self and
# test_fragile.py::test35_meta_class_attribute_lookup.
class Regression03OldCrashers( MyTestCase ):
   def test2UsageOfTQClassInstance( self ):
      """Calls on a TQClass instance"""

      self.assertEqual( TClass.GetClass("TQClass").GetName(), "TQClass" )


### Test the condition under which to (not) start the GUI thread =============
class Regression04Threading( MyTestCase ):

   hasThread = gROOT.IsBatch() and 5 or 6   # can't test if no display ...
   noThread  = 5

   def test1SpecialCasegROOT( self ):
      """Test the special role that gROOT plays vis-a-vis threading"""

      # Restore original LD_PRELOAD when running under AddressSanitizer
      if original_preload is not None:
         os.environ['LD_PRELOAD'] = original_preload

      cmd = sys.executable + "  -c 'import sys, ROOT; ROOT.gROOT; %s "\
            "sys.exit( 5 + int(\"thread\" in ROOT.__dict__) )'"
      if self.hasThread == self.noThread:
         cmd += " - -b"

      stat, out = commands.getstatusoutput( cmd % "" )
      self.assertEqual( WEXITSTATUS(stat), self.noThread )

      stat, out = commands.getstatusoutput( cmd % "ROOT.gROOT.SetBatch( 1 );" )
      self.assertEqual( WEXITSTATUS(stat), self.noThread )

      stat, out = commands.getstatusoutput( cmd % "ROOT.gROOT.SetBatch( 0 );" )
      self.assertEqual( WEXITSTATUS(stat), self.noThread )

      stat, out = commands.getstatusoutput(
         cmd % "ROOT.gROOT.ProcessLine( \"cout << 42 << endl;\" ); " )
      self.assertEqual( WEXITSTATUS(stat), self.hasThread )

      stat, out = commands.getstatusoutput( cmd % "ROOT.gDebug;" )
      self.assertEqual( WEXITSTATUS(stat), self.hasThread )

      # Restore the cleaned LD_PRELOAD for other tests
      if cleaned_preload is not None:
         os.environ['LD_PRELOAD'] = cleaned_preload

   def test2ImportStyles( self ):
      """Test different import styles vis-a-vis threading"""

      # Restore original LD_PRELOAD when running under AddressSanitizer
      if original_preload is not None:
         os.environ['LD_PRELOAD'] = original_preload

      cmd = sys.executable + " -c 'import sys; %s ;"\
            "import ROOT; sys.exit( 5 + int(\"thread\" in ROOT.__dict__) )'"
      if self.hasThread == self.noThread:
         cmd += " - -b"

      stat, out = commands.getstatusoutput( cmd % "from ROOT import gROOT" )
      self.assertEqual( WEXITSTATUS(stat), self.noThread )

      stat, out = commands.getstatusoutput( cmd % "from ROOT import gDebug" )
      self.assertEqual( WEXITSTATUS(stat), self.hasThread )

      # Restore the cleaned LD_PRELOAD for other tests
      if cleaned_preload is not None:
         os.environ['LD_PRELOAD'] = cleaned_preload

   def test3SettingOfBatchMode( self ):
      """Test various ways of preventing GUI thread startup"""

      # Restore original LD_PRELOAD when running under AddressSanitizer
      if original_preload is not None:
         os.environ['LD_PRELOAD'] = original_preload

      cmd = sys.executable + " -c '%s import ROOT, sys; sys.exit( 5+int(\"thread\" in ROOT.__dict__ ) )'"
      if self.hasThread == self.noThread:
         cmd += " - -b"

      stat, out = commands.getstatusoutput(
         cmd % 'import ROOT; ROOT.PyConfig.StartGuiThread = 0;' )
      self.assertEqual( WEXITSTATUS(stat), self.noThread )

      stat, out = commands.getstatusoutput(
         cmd % 'from ROOT import PyConfig; PyConfig.StartGuiThread = 0; from ROOT import gDebug;' )
      self.assertEqual( WEXITSTATUS(stat), self.noThread )

      stat, out = commands.getstatusoutput(
         cmd % 'from ROOT import PyConfig; PyConfig.StartGuiThread = 1; from ROOT import gDebug;' )
      self.assertEqual( WEXITSTATUS(stat), self.hasThread )

      # Restore the cleaned LD_PRELOAD for other tests
      if cleaned_preload is not None:
         os.environ['LD_PRELOAD'] = cleaned_preload


### test pythonization and operators of TVector3 ===========================================
class Regression09TVector3Pythonize( MyTestCase ):
   def test1TVector3( self ):
      """Verify TVector3 pythonization"""

      v = TVector3( 1., 2., 3.)
      self.assertEqual( list(v), [1., 2., 3. ] )

      w = 2*v
      self.assertEqual( w.x(), 2*v.x() )
      self.assertEqual( w.y(), 2*v.y() )
      self.assertEqual( w.z(), 2*v.z() )

   def test2TVector3(self):
      """Verify that using one operator* overload does not mask the others"""
      # ROOT-10278
      v = TVector3(1., 2., 3.)
      v*2
      self.assertEqual(v*v, 14.0)


# NOTE: the former Regression10CoralAttributeListIterators was migrated to
# the cppjit test suite: test_stltypes.py::test05_iteration_with_base_iterator.


### entities in the ROOT:: namespace =========================================
class Regression11GlobalsLookup( MyTestCase ):
   def test2GlobalFromROOTNamespace( self ):
      """Entities in 'ROOT::' need no explicit 'ROOT.'"""

      import ROOT
      m = ROOT.Math


### importing cout should not result in printed errors =======================
class Regression12WriteTGraph( MyTestCase ):
   def test1WriteTGraph( self ):
      """Write a TGraph object and read it back correctly"""

      gr = TGraph()
      ff = TFile( "test.root", "RECREATE" )
      ff.WriteObject( gr, "grname", "" )
      # In new PyROOT, use a nicer way to get objects in files,
      # the TDirectory::Get() pythonisation:
      ff.Get("grname")
      os.remove( "test.root" )


# NOTE: the former Regression14TPyException was migrated to the cppjit test
# suite: test_regression.py::test64_pyexception_construction.


### matrix access has to go through non-const lookup =========================
class Regression17MatrixD( MyTestCase ):
   def test1MatrixElementAssignment( self ):
      """Matrix lookup has to be non-const to allow assigment"""

      m = TMatrixD( 5, 5 )
      self.assertTrue( not 'const' in type(m[0]).__name__ )

    # test assignment
      m[1][2] = 3.
      self.assertEqual( m[1][2], 3. )

      m[1, 2] = 4.
      self.assertEqual( m[1][2], 4. )


### Tests for TGL classes ================
try:
   from ROOT import TGLLine3, TGLVertex3, TGLVector3
except ImportError:
   print("GL classes not found, skipping GL tests")
else:
   class Regression19TGL(MyTestCase):
      def test1TGLVertex3OperatorPlus(self):
         """Try invoking TGLVertex3::operator+ twice"""
         # ROOT-10166
         scatteringPoint = TGLVertex3(2., 3., 0.2)
         glvec3 = TGLVector3(1,2,3)

         vertexEnd = scatteringPoint + glvec3
         vertexEnd = scatteringPoint + glvec3

      def test2TGLLine3Constructor(self):
         """Check that the right constructor of TGLLine3 is called"""
         # ROOT-10102
         trackAfterScattering = TGLLine3(TGLVertex3(2., 3., 0.2), TGLVector3(0., 0., -20.))

         self.assertEqual(trackAfterScattering.Vector().X(), .0)
         self.assertEqual(trackAfterScattering.Vector().Y(), .0)
         self.assertEqual(trackAfterScattering.Vector().Z(), -20.0)


### Getting and setting configuration options of gEnv ================
class Regression20gEnv(MyTestCase):
   def test1GetSetValue(self):
      """Set a value with gEnv and retrieve it afterwards"""
      # ROOT-10155
      from ROOT import gEnv

      optname = "SomeOption"
      defval = -1
      self.assertEqual(gEnv.GetValue(optname, defval), defval)
      newval = 0
      gEnv.SetValue(optname, newval)
      self.assertEqual(gEnv.GetValue(optname, defval), newval)

### Tests related to cleanup of proxied objects ================
class Regression22ObjectCleanup(MyTestCase):
   def test1GetListOfGraphs(self):
      """List returned by GetListOfGraphs should not have kMustCleanup set to true"""
      # ROOT-9040
      mg = ROOT.TMultiGraph()
      tg = ROOT.TGraph()
      # The TMultiGraph will take the ownership of the added TGraphs
      ROOT.SetOwnership(tg, False)
      mg.Add(tg)

      l = mg.GetListOfGraphs()
      self.assertEqual(l.TestBit(ROOT.kMustCleanup), False)

      c = ROOT.TCanvas()
      mg.Draw()


class Regression23TFractionFitter(MyTestCase):
   def test1TFractionFitterDestruction(self):
      """Test proper destruction of TFractionFitter object"""
      # ROOT-9414
      h1 = ROOT.TH1F("h1","h1",1,0,1)
      h2 = ROOT.TH1F("h2","h2",1,0,1)
      h3 = ROOT.TH1F("h3","h3",1,0,1)

      h1.Fill(0.5)
      h2.Fill(0.5)
      h3.Fill(0.5)
      h3.Fill(0.5)

      mc = ROOT.TObjArray(2)
      mc.Add(h1)
      mc.Add(h2)

      ff = ROOT.TFractionFitter(h3, mc)
      ff.Fit()


class Regression24CppPythonInheritance(MyTestCase):
   # NOTE: the former tests 01-07, 09 and 10 of this class were removed: they
   # are covered upstream in the cppyy test suite by test_crossinheritance.py
   # (test24_non_copyable, test14_protected_access, test13_virtual_dtors_and_del,
   # test35_deletion, test02_constructor, test19/test20 multiple inheritance,
   # test22_multiple_inheritance_with_defaults, test30_access_and_overload,
   # test17_deep_hierarchy, test28_cross_deep).

   def test08ConstructorAllDefaultPars(self):
       """Invocation of a constructor that has default values for all its parameters"""
       # 6578
       class pMainFrame(ROOT.TGMainFrame):
           def __init__(self, parent, width, height ):
               ROOT.TGMainFrame.__init__(self, parent, width, height)

       window = pMainFrame(ROOT.gClient.GetRoot(), 200, 200)

   # NOTE: the former test11MultiInheritancePyCpp was migrated to the cppjit
   # test suite: test_crossinheritance.py::test40_mixed_python_cpp_bases.


# NOTE: the former Regression25MapGetItemToCall, Regression26OverloadedOperator
# and Regression27ImplicitSmartPtrOverload were migrated to the cppjit test
# suite: test_operators.py::test19_base_index_operator_with_derived_call_operator,
# test_advancedcpp.py::test35_using_of_base_operator_call and
# test_cpp11features.py::test22_smart_ptr_overload_and_rvalue_ref.

class Regression28EnumResolution(MyTestCase):
   def test01EnumResolution(self):
      """Test that enums are correctly resolved when using an alias within a class/struct"""
      # 9246
      
      cppyy.cppdef('''
         enum Enum { a,b,c };
         ''')

      ea = cppyy.gbl.Enum.a
      eb = cppyy.gbl.Enum.b
      self.assertEqual(ea, 0)
      self.assertEqual(eb, 1)

      cppyy.cppdef('''
              struct Track { using Type = Enum; };
              ''')
      trackea = cppyy.gbl.Track.Type.a
      trackeb = cppyy.gbl.Track.Type.b

      self.assertEqual(trackea, 0)
      self.assertEqual(trackeb, 1)

class Regression29TemplateFunctionOverloadSFINAE(MyTestCase):
   def test01TemplateFunctionsOverloadWithSFINAE(self):
      """Test template functions overload with SFINAE """
      # 15085
      
      cppyy.cppdef('''
#include <type_traits> 
#include <string> 
template <unsigned int  N>
class A
{
public:
//
template <unsigned int K,
typename std::enable_if<(K<N),int>::type = 0 >
std::string fun_1 () { return "ququ" ;} 
//
template <unsigned int K,
typename std::enable_if<(N<=K),int>::type = 0 >
int         fun_1 () { return K ;} 
};
// the same with functions
template <unsigned int K,unsigned int N, 
typename std::enable_if<(K<N),int>::type = 0 >
inline std::string fun_2 ( A<N>& ) { return "ququ" ;} 
//
template <unsigned int K,unsigned int N,
typename std::enable_if<(N<=K),int>::type = 0 >
inline int         fun_2 ( A<N>& ) { return K ;} 
         ''')


      output = ['ququ', 'ququ', 'ququ', 'ququ', 4, 5, 6, 7, 8]
      N = 4 
      a = cppyy.gbl.A [ N ]()

      for k in range ( 2 * N + 1 ) : 
         self.assertEqual(a.fun_1[k](), output[k])

      for k in range( 2 * N + 1  ) :
         func = cppyy.gbl.fun_2[k,N]
         self.assertEqual(func(a), output[k])

class Regression30TemplatedBinaryOperatorImplementations(MyTestCase):
   def test01TemplatedBinaryOperatorImplementations(self):
      """Test templated binary operator implementations like operator+"""
      # 15378

      cppyy.cppdef("""

template<class T>
struct Foo {
    Foo(T v) : val{v} {}
    T val;
};

template <class T>
Foo<T> operator+(Foo<T> const &a, Foo<T> const &b)
{
   return Foo<T>{a.val + b.val};
}

      """)

      klass = cppyy.gbl.Foo["int"]

      a = klass(4)
      b = klass(7)

      c = a + b
      self.assertEqual(c.val, 11)

class Regression31AccessingProtectedDataMember(MyTestCase):
   def test01AccessingProtectedDataMember(self):
      """Test accessibility of protected data members"""
      # 16852

      cppyy.cppdef("""
class MyBaseClass {
public:
   virtual ~MyBaseClass() = default;
protected:
   int protectedFunc() { return 5; }
   int _protectedData = 4;
};
""")

      class MyDerivedClass(cppyy.gbl.MyBaseClass):
         pass

      my_obj = MyDerivedClass()

      self.assertEqual(my_obj.protectedFunc(), 5)
      self.assertEqual(my_obj._protectedData, 4)  


class Regression32StdFunctionForTemplates(MyTestCase):
   def test01StdFunctionForTemplates(self):
      """Test cppyy.gbl.std.function(...) for template types like std::vector<int>, ROOT::RVec<int>"""
      # 19358

      cppyy.gbl.std.function("bool(int)")
      cppyy.gbl.std.function("bool()")
      cppyy.gbl.std.function("bool(std::vector<int>)")
      cppyy.gbl.std.function("bool(ROOT::RVec<int>)")

class Regression33NameSpaceDroppedBinding(MyTestCase):
   def test01NameSpaceDroppedBinding(self):
      """Test that the namespace is not dropped in std::function cppyy binding"""
      # 10680

      cppyy.cppdef("""
#include <functional>

namespace foo {
class Foo {
  public:
    Foo() : bar(0) {};
    int bar;
};
}

void apply(std::function<void(foo::Foo*)> f) {
    foo::Foo x;
    f(&x);
}
   """)
      cppyy.gbl.apply(lambda x: print(x.bar))

## actual test run
if __name__ == '__main__':
   from MyTextTestRunner import MyTextTestRunner

   loader = unittest.TestLoader()
   testSuite = loader.loadTestsFromModule( sys.modules[ __name__ ] )

   runner = MyTextTestRunner( verbosity = 2 )
   result = not runner.run( testSuite ).wasSuccessful()

   sys.exit( result )
