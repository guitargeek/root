# Author: Enric Tejedor CERN  06/2018

################################################################################
# Copyright (C) 1995-2018, Rene Brun and Fons Rademakers.                      #
# All rights reserved.                                                         #
#                                                                              #
# For the licensing terms see $ROOTSYS/LICENSE.                                #
# For the list of contributors see $ROOTSYS/README/CREDITS.                    #
################################################################################

r'''
/**
\class TTree
\brief \parblock \endparblock
\htmlonly
<div class="pyrootbox">
\endhtmlonly
## PyROOT

The TTree class has several additions for its use from Python, which are also
available in its subclasses e.g. TChain and TNtuple.

First, TTree instances are iterable in Python. Therefore, assuming `t` is
a TTree instance, we can do:
\code{.py}
for entry in t:
    x = entry.branch_name
    ...
\endcode

At each iteration, a new entry of the tree will be read. In the code above,
`entry` allows to access the branch values for the current entry. This can be
done with the syntax `entry.branch_name` or, if the branch name is incompatible
with Python naming rules, with e.g. "getattr(entry, '1_branch_name')".

<em>Please note</em> that iterating in Python can be slow, so only iterate over
a tree as described above if performance is not an issue or when dealing with
a small dataset. To read and process the entries of a tree in a much faster
way, please use ROOT::RDataFrame.

Second, a couple of TTree methods have been modified to facilitate their use
from Python: TTree::Branch and TTree::SetBranchAddress.

Regarding TTree::Branch, the following example shows how we can create
different types of branches of a TTree. Note that `Branch` will just link
the new branch with a given Python object, so it is still necessary to fill
such object with the desired content before calling TTree::Fill.
\code{.py}
from array import array
import numpy as np
import ROOT
from ROOT import addressof

# Basic type branch (float) - use array of length 1
n = array('f', [ 1.5 ])
t.Branch('floatb', n, 'floatb/F')

# Array branch - use array of length N
N = 10
a = array('d', N*[ 0. ])
t.Branch('arrayb', a, 'arrayb[' + str(N) + ']/D')

# Array branch - use NumPy array of length N
npa = np.array(N*[ 0. ])
t.Branch('nparrayb', npa, 'nparrayb[' + str(N) + ']/D')

# std::vector branch
v = ROOT.std.vector('double')(N*[ 0. ])
t.Branch('vectorb0', v)

# Class branch / struct in single branch
cb = ROOT.MyClass()
t.Branch('classb', cb)

# Struct as leaflist
# Assuming:
# struct MyStruct {
#   int myint;
#   float myfloat;
# };
ms = ROOT.MyStruct()
t.Branch('structll', ms, 'myint/I:myfloat/F')

# Store struct members individually
ms = ROOT.MyStruct()
# Use `addressof` to get the address of the struct members
t.Branch('myintb', addressof(ms, 'myint'), 'myint/I')
t.Branch('myfloatb', addressof(ms, 'myfloat'), 'myfloat/F')
\endcode

Concerning TTree::SetBranchAddress, below is an example of prepare
the reading of different types of branches of a TTree. Note that
`SetBranchAddress` will just link a given branch with a certain
Python object; after that, in order to read the content of such
branch for a given TTree entry `x`, TTree::GetEntry(x) must be
invoked.
\code{.py}
from array import array
import numpy as np
import ROOT

# Basic type branch (float) - use array of length 1
n = array('f', [ 0. ])
t.SetBranchAddress('floatb', n)

# Array branch - use array of length N
N = 10
a = array('d', N*[ 0. ])
t.SetBranchAddress('arrayb', a)

# Array branch - use NumPy array of length N
npa = np.array(N*[ 0. ])
t.SetBranchAddress('nparrayb', a)

# std::vector branch
v = ROOT.std.vector('double')()
t.SetBranchAddress('vectorb', v)

# Class branch
cb = ROOT.MyClass()
t.SetBranchAddress('classb', cb)

# Struct branch (both single-branch and leaf list)
ms = ROOT.MyStruct()
ds.SetBranchAddress('structb', ms)
\endcode
\htmlonly
</div>
\endhtmlonly
*/
'''

from libROOTPythonizations import AddBranchAttrSyntax, BranchPyz, CreateBufferFromAddress
from . import pythonization

"""
\brief Add pythonization for TTree::SetBranchAddress.
\param[in] self Always null, since this is a module function.
\param[in] args Pointer to a Python tuple object containing the arguments
received from Python.

Modify the behaviour of SetBranchAddress so that proxy references can be passed
as arguments from the Python side, more precisely in cases where the C++
implementation of the method expects the address of a pointer.

For example:
```
v = ROOT.std.vector('int')()
t.SetBranchAddress("my_vector_branch", v)
```
"""
def SetBranchAddressPyz(tree, name, address):

    import cppyy

    branch = tree.GetBranch(name)
    if not branch:
        raise RuntimeError("TTree::SetBranchAddress must be called with a valid branch name")

    is_leaf_list = branch.ClassName() == "TBranch"

    print(address)
    print(type(address))
    print("def", str(cppyy.addressof(address)))
    print("tru", str(cppyy.addressof(address, byref=True)))
    print("fls", str(cppyy.addressof(address, byref=False)))
    try:
        return tree.SetBranchAddress(name, cppyy.addressof(address, byref=False))
    except TypeError:
        return None

    """
    void *buf = 0;
    if (CPPInstance_Check(address)) {
       ((CPPInstance *)address)->GetDatamemberCache(); // force creation of cache

       if (((CPPInstance *)address)->fFlags & CPPInstance::kIsReference || is_leaf_list)
          buf = (void *)((CPPInstance *)address)->GetObject();
       else
          buf = (void *)&(((CPPInstance *)address)->GetObjectRaw());
    } else
       Utility::GetBuffer(address, '*', 1, buf, false);

    if buf:
       auto res = tree->SetBranchAddress(CPyCppyy_PyText_AsString(name), buf);
       return PyInt_FromLong(res);
    """

    return None

# TTree iterator
def _TTree__iter__(self):
    i = 0
    bytes_read = self.GetEntry(i)
    while 0 < bytes_read:
        yield self
        i += 1
        bytes_read = self.GetEntry(i)

    if bytes_read == -1:
        raise RuntimeError("TTree I/O error")

def _SetBranchAddress(self, *args):
    res = None
    # Modify the behaviour if args is (const char*, void*)
    if len(args) == 2:
        res = SetBranchAddressPyz(self, *args)

    if res is None:
        # Fall back to the original implementation for the rest of overloads
        res = self._OriginalSetBranchAddress(*args)

    return res

def _Branch(self, *args):
    # Modify the behaviour if args is one of:
    # ( const char*, void*, const char*, Int_t = 32000 )
    # ( const char*, const char*, T**, Int_t = 32000, Int_t = 99 )
    # ( const char*, T**, Int_t = 32000, Int_t = 99 )
    res = BranchPyz(self, *args)

    if res is None:
        # Fall back to the original implementation for the rest of overloads
        res = self._OriginalBranch(*args)

    return res

@pythonization('TTree')
def pythonize_ttree(klass, name):
    # Parameters:
    # klass: class to be pythonized
    # name: string containing the name of the class

    # Pythonizations that are common to TTree and its subclasses.
    # To avoid duplicating the same logic in the pythonizors of
    # the subclasses, inject the pythonizations for all the target
    # classes here.

    # Pythonic iterator
    klass.__iter__ = _TTree__iter__

    # tree.branch syntax
    AddBranchAttrSyntax(klass)

    # SetBranchAddress
    klass._OriginalSetBranchAddress = klass.SetBranchAddress
    klass.SetBranchAddress = _SetBranchAddress

    # Branch
    klass._OriginalBranch = klass.Branch
    klass.Branch = _Branch

@pythonization('TChain')
def pythonize_tchain(klass):
    # Parameters:
    # klass: class to be pythonized

    # TChain needs to be explicitly pythonized because it redefines
    # SetBranchAddress in C++. As a consequence, TChain does not
    # inherit TTree's pythonization for SetBranchAddress, which
    # needs to be injected to TChain too. This is not the case for
    # other classes like TNtuple, which will inherit all the
    # pythonizations added here for TTree.

    # SetBranchAddress
    klass._OriginalSetBranchAddress = klass.SetBranchAddress
    klass.SetBranchAddress = _SetBranchAddress
