# Author: Enric Tejedor CERN  06/2018

################################################################################
# Copyright (C) 1995-2018, Rene Brun and Fons Rademakers.                      #
# All rights reserved.                                                         #
#                                                                              #
# For the licensing terms see $ROOTSYS/LICENSE.                                #
# For the list of contributors see $ROOTSYS/README/CREDITS.                    #
################################################################################

r"""
\pythondoc TTree

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

Two methods of TTree have been pythonized to facilitate their: TTree::Branch and
TTree::SetBranchAddress.

### Pythonization of TTree::Branch

The following example shows how we can create different types of branches of a TTree.
`Branch` links the new branch with a given Python object. It is therefore possible to
fill such object with the desired content before calling TTree::Fill.

\code{.py}
from array import array
import numpy as np
import ROOT

# We create the file and the tree
with ROOT.TFile("outfile.root", "RECREATE") as ofile:
    t = ROOT.TTree("mytree", "mytree")

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
    cb = ROOT.TH1D("myHisto", "myHisto", 64, -4, 4)
    # This could have been any class known to ROOT, also custom
    #cb = ROOT.MyCustomClass()
    t.Branch('classb', cb)

    # Struct as leaflist. This is interpreted on the fly,
    # but could be known to ROOT by other means, such as
    # header inclusion or dictionary load.
    ROOT.gInterpreter.Declare('''
    struct MyStruct {
    int myint;
    float myfloat;
    };
    ''')
    ms = ROOT.MyStruct()
    t.Branch('structll', ms, 'myint/I:myfloat/F')

    # Store struct members individually
    ms = ROOT.MyStruct()
    # Use the `addressof` function in the ROOT module
    # to get the address of the struct members
    t.Branch('myintb', ROOT.addressof(ms, 'myint'), 'myint/I')
    t.Branch('myfloatb', ROOT.addressof(ms, 'myfloat'), 'myfloat/F')

    # Let's write one entry in our tree
    t.Fill()
    # Finally flush the content of the tree to the file
    t.Write()
\endcode

### Pythonization of TTree::SetBranchAddress

This section is to be considered for advanced users. Simple event
loops reading tree entries in Python can be performed as shown above.

Below an example is shown of reading different types tree branches.
Note that `SetBranchAddress` will just link a given branch with a
certain Python object; after that, in order to read the content of such
branch for a given TTree entry `x`, TTree::GetEntry(x) must be
invoked.

\code{.py}
from array import array
import numpy as np
import ROOT

with ROOT.TFile('outfile.root') as infile:

    t = infile['mytree']

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
    cb = ROOT.TH1D()
    # Any other class known to ROOT would have worked
    #cb = ROOT.MyClass()
    t.SetBranchAddress('classb', cb)

    # Struct as leaflist. This is interpreted on the fly,
    # but could be known to ROOT by other means, such as
    # header inclusion or dictionary load.
    ROOT.gInterpreter.Declare('''
    struct MyStruct {
    int myint;
    float myfloat;
    };
    ''')
    ms = ROOT.MyStruct()
    t.SetBranchAddress('structll', ms)

    t.GetEntry(0)
\endcode

\endpythondoc
"""

from . import pythonization
from ._memory_utils import (
    _constructor_releasing_ownership,
    _SetDirectory_SetOwnership,
    _should_give_up_ownership,
)
from ._rvec import _get_cpp_type_from_numpy_type


# Lazily declared C++ helpers for accessing raw TTree/TBranch/TLeaf pointers
# that cppyy would otherwise try to convert into Python types (e.g. char* to
# Python string). They are registered with Cling on first use.
_ttree_helpers_declared = False


def _declare_ttree_helpers():
    global _ttree_helpers_declared
    if _ttree_helpers_declared:
        return
    import cppyy
    cppyy.cppdef(r"""
namespace _pyroot_ttree {
   inline std::intptr_t branch_address(TBranch* b) {
      return (std::intptr_t)b->GetAddress();
   }
   inline std::intptr_t leaf_value_pointer(TLeaf* l) {
      return (std::intptr_t)l->GetValuePointer();
   }
   inline std::intptr_t deref_voidp(std::intptr_t addr) {
      return (std::intptr_t)(*(void**)addr);
   }
   inline std::intptr_t branch_element_object(TBranchElement* be) {
      return (std::intptr_t)be->GetObject();
   }
   inline Long_t branch_element_offset(TBranchElement* be) {
      return ((TStreamerElement*)be->GetInfo()->GetElements()->At(be->GetID()))->GetOffset();
   }
}
""")
    _ttree_helpers_declared = True


def _get_python_buffer_address(pyobject):
    """Return the address of a Python buffer-like object (array.array, numpy
    array, bytearray, ...) as an integer, or None if no address can be obtained.
    Text-like objects (str, bytes) are explicitly rejected.
    """
    import ctypes

    # Exclude text-like objects (policy decision, matches the old C++ helper)
    if isinstance(pyobject, (bytes, str)):
        return None

    # numpy arrays and array-interface objects
    iface = getattr(pyobject, "__array_interface__", None)
    if iface is not None:
        return iface["data"][0]

    # array.array objects
    if hasattr(pyobject, "buffer_info"):
        return pyobject.buffer_info()[0]

    # bytearray
    if isinstance(pyobject, bytearray):
        return ctypes.addressof((ctypes.c_char * len(pyobject)).from_buffer(pyobject))

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


def _pythonize_branch_addr(branch, addr_orig):
    """Helper for the SetBranchAddress pythonization, extracting the relevant
    address from a Python object if possible.
    """
    import ctypes

    import ROOT

    is_leaf_list = branch.IsA() is ROOT.TBranch.Class()

    if is_leaf_list:
        # If the branch is a leaf list, SetBranchAddress expects the
        # address of the object that has the corresponding data members.
        return ctypes.c_void_p(ROOT._cppyy.addressof(instance=addr_orig, byref=False))

    # Otherwise, SetBranchAddress is expecting a pointer to the address of
    # the object, and the pointer needs to stay alive. Therefore, we create
    # a container for the pointer and cache it in the original cppyy proxy.
    addr_view = ROOT.array["std::intptr_t", 1]([ROOT._cppyy.addressof(instance=addr_orig, byref=False)])

    if not hasattr(addr_orig, "_set_branch_cached_pointers"):
        addr_orig._set_branch_cached_pointers = []
    addr_orig._set_branch_cached_pointers.append(addr_view)

    # Finally, we have to return the address of the container
    return ctypes.c_void_p(ROOT._cppyy.addressof(instance=addr_view, byref=False))


def _get_cpp_type_from_array_typecode(typecode):
    # Complete list from https://docs.python.org/3/library/array.html
    c_type_names = {
        "b": "signed char",
        "B": "unsigned char",
        "u": "wchar_t",
        "h": "signed short",
        "H": "unsigned short",
        "i": "signed int",
        "I": "unsigned int",
        "l": "signed long",
        "L": "unsigned long",
        "q": "signed long long",
        "Q": "unsigned long long",
        "f": "float",
        "d": "double",
    }
    return c_type_names[typecode]


def _determine_data_type(addr):
    """Figure out data_type in case addr is a numpy.ndarray or array.array."""

    # For NumPy arrays
    if hasattr(addr, "__array_interface__"):
        return _get_cpp_type_from_numpy_type(addr.__array_interface__["typestr"][1:])

    # For the builtin array library
    if hasattr(addr, "buffer_info"):
        return _get_cpp_type_from_array_typecode(addr.typecode)

    return None


def _SetBranchAddress(self, bname, addr, *args, **kwargs):
    """
    Pythonization for TTree::SetBranchAddress.

    Modify the behaviour of SetBranchAddress so that proxy references can be passed
    as arguments from the Python side, more precisely in cases where the C++
    implementation of the method expects the address of a pointer.

    For example:
    ```
    v = ROOT.std.vector('int')()
    t.SetBranchAddress("my_vector_branch", v)
    ```
    """
    import cppyy

    import ROOT

    branch = self.GetBranch(bname)

    # Pythonization for cppyy proxies (of type CPPInstance)
    if isinstance(addr, ROOT._cppyy.types.Instance):
        addr = _pythonize_branch_addr(branch, addr)

    # Figure out data_type in case addr is a numpy.ndarray or array.array
    data_type = _determine_data_type(addr)

    if data_type is None:
        return self._OriginalSetBranchAddress(bname, addr, *args, **kwargs)

    # In the case the data_type is available, we would like to call the
    # template overload of SetBranchAddress instantiatied for that type.
    # However, there are two such overloads candidates:
    #
    #   template <class T> int TTree::SetBranchAddress(const char *bname, T **add, ...);
    #   template <class T> int TTree::SetBranchAddress(const char *bname, T *add, ...);
    #
    # The cppyy bindings can't make a meaningful selection here as Python is
    # lacking pointer semantics, so it considers both overloads as valid
    # choices. In the past, we just happened to be lucky that it tried the T *
    # overload first, which is the one we need. But as cppyy becomes more
    # strict about overload resolution ambiguity errors, this won't work
    # anymore. That's why we re-implement what happens in the template overload
    # on the Python side.

    cl = ROOT.TClass.GetClass[data_type]()
    tp = ROOT.kOther_t
    if not cl:
        tp = ROOT.TDataType.GetType(cppyy.typeid(getattr(ROOT, data_type)))

    # Extract the TBranch ptr argument if available
    tbranch_ptr = ROOT.nullptr
    if len(args) > 0:
        tbranch_ptr = args[0]
    elif "ptr" in kwargs:
        tbranch_ptr = kwargs["ptr"]

    return self._OriginalSetBranchAddress(bname, addr, ptr=tbranch_ptr, realClass=cl, datatype=tp, isptr=False)


def _try_branch_leaflist_overload(tree, args):
    """Try matching ``Branch(name, address, leaflist, bufsize=32000)``.

    Returns the created TBranch on success, or None on signature mismatch.
    """
    import ctypes
    import ROOT

    if len(args) < 3 or len(args) > 4:
        return None
    name, address, leaflist = args[0], args[1], args[2]
    if not isinstance(name, str) or not isinstance(leaflist, str):
        return None
    bufsize = None
    if len(args) == 4:
        bufsize = args[3]
        if not isinstance(bufsize, int):
            return None

    if isinstance(address, ROOT._cppyy.types.Instance):
        buf_addr = ROOT._cppyy.addressof(instance=address, byref=False)
    else:
        buf_addr = _get_python_buffer_address(address)
    if not buf_addr:
        return None

    buf = ctypes.c_void_p(buf_addr)
    if bufsize is None:
        return tree._OriginalBranch(name, buf, leaflist)
    return tree._OriginalBranch(name, buf, leaflist, bufsize)


def _try_branch_ptr_to_ptr_overloads(tree, args):
    """Try matching one of the two ``Branch`` ptr-to-ptr overloads:

    - ``Branch(name, classname, address, bufsize=32000, splitlevel=99)``
    - ``Branch(name, address, bufsize=32000, splitlevel=99)``

    Returns the created TBranch on success, or None on signature mismatch.
    """
    import ctypes
    import ROOT

    if not args or not isinstance(args[0], str):
        return None
    name = args[0]

    classname = None
    if len(args) >= 2 and isinstance(args[1], str):
        # ( name, classname, address, [bufsize, splitlevel] )
        if len(args) < 3 or len(args) > 5:
            return None
        classname = args[1]
        address = args[2]
        extra = args[3:]
    else:
        # ( name, address, [bufsize, splitlevel] )
        if len(args) < 2 or len(args) > 4:
            return None
        address = args[1]
        extra = args[2:]

    if any(not isinstance(x, int) for x in extra):
        return None

    # The Branch(name, classname, T**, ...) overloads need a stable T** that
    # lives long enough. For a CPPInstance, we allocate a 1-element array
    # holding the object's address and pass &array[0]. The array is cached on
    # the proxy to keep it alive while the tree holds the pointer.
    if isinstance(address, ROOT._cppyy.types.Instance):
        if classname is None:
            classname = type(address).__cpp_name__
        obj_addr = ROOT._cppyy.addressof(instance=address, byref=False)
        addr_holder = ROOT.array["std::intptr_t", 1]([obj_addr])
        if not hasattr(address, "_branch_cached_pointers"):
            address._branch_cached_pointers = []
        address._branch_cached_pointers.append(addr_holder)
        buf_addr = ROOT._cppyy.addressof(instance=addr_holder, byref=False)
    else:
        if classname is None:
            return None
        buf_addr = _get_python_buffer_address(address)

    if not buf_addr:
        return None

    buf = ctypes.c_void_p(buf_addr)
    return tree._OriginalBranch(name, classname, buf, *extra)


def _Branch(self, *args):
    """Pythonization for TTree::Branch.

    Modify the behaviour of Branch so that proxy references can be passed as
    arguments from the Python side, more precisely in cases where the C++
    implementation of the method expects the address of a pointer. The
    following signatures are handled:

    - ( const char*, void*, const char*, Int_t = 32000 )
    - ( const char*, const char*, T**, Int_t = 32000, Int_t = 99 )
    - ( const char*, T**, Int_t = 32000, Int_t = 99 )
    """
    res = _try_branch_leaflist_overload(self, args)
    if res is not None:
        return res

    res = _try_branch_ptr_to_ptr_overloads(self, args)
    if res is not None:
        return res

    return self._OriginalBranch(*args)


def _search_for_branch(tree, name):
    branch = tree.GetBranch(name)
    if not branch:
        # Sub-branches often carry a trailing '.'
        branch = tree.GetBranch(name + ".")
    return branch


def _search_for_leaf(tree, name, branch):
    leaf = tree.GetLeaf(name)
    if branch and not leaf:
        leaf = branch.GetLeaf(name)
        if not leaf:
            leaves = branch.GetListOfLeaves()
            if leaves.GetSize() and leaves.First() == leaves.Last():
                leaf = leaves.At(0)
    return leaf


def _get_multi_dims(title):
    """Parse a TLeaf title of the form ``name[dim1][dim2]...`` into a list of
    static dimensions.
    """
    dims = []
    i = 0
    while True:
        lb = title.find("[", i)
        if lb < 0:
            break
        rb = title.find("]", lb)
        if rb < 0:
            break
        sub = title[lb + 1:rb]
        if sub:
            dims.append(int(sub))
        i = rb + 1
    return dims


def _resolve_branch(tree, name, branch):
    """Return ``(address, type_name)`` for the object held by *branch*, where
    *address* is an integer pointer suitable for ``bind_object`` and
    *type_name* is the C++ class name without trailing '*'. If the branch
    cannot be represented as a wrapped object, return ``(None, "")``.
    """
    import ROOT

    _declare_ttree_helpers()
    helpers = ROOT._cppyy.gbl._pyroot_ttree

    TBranchElement = ROOT.TBranchElement
    TBranchObject = ROOT.TBranchObject

    # Partial split-object branch: a TBranchElement that points into a member
    # of its parent object.
    if branch.InheritsFrom(TBranchElement.Class()):
        be = ROOT._cppyy.bind_object(
            ROOT._cppyy.addressof(instance=branch, byref=False), TBranchElement
        )
        cur_cls = be.GetCurrentClass()
        tgt_cls = be.GetTargetClass()
        if cur_cls and cur_cls != tgt_cls and be.GetID() >= 0:
            offset = helpers.branch_element_offset(be)
            obj_addr = helpers.branch_element_object(be)
            return obj_addr + int(offset), cur_cls.GetName()

    # Full object branch (TBranchElement or TBranchObject)
    branch_cls = branch.IsA()
    if branch_cls is TBranchElement.Class() or branch_cls is TBranchObject.Class():
        baddr = helpers.branch_address(branch)
        if baddr:
            # *(void**)branch.GetAddress() - dereference one level
            return helpers.deref_voidp(baddr), branch.GetClassName()

        # Try leaf, otherwise mark as failure (with a typed null object)
        leaves = branch.GetListOfLeaves()
        if (not tree.GetLeaf(name)) and (
            not (leaves.GetSize() and leaves.First() == leaves.Last())
        ):
            return 0, branch.GetClassName()

    return None, ""


def _wrap_leaf(leaf):
    """Wrap a TLeaf into a Python value or LowLevelView, mirroring the
    behaviour of the previous C++ ``WrapLeaf`` helper.
    """
    import ROOT

    _declare_ttree_helpers()
    helpers = ROOT._cppyy.gbl._pyroot_ttree

    type_name = leaf.GetTypeName()
    n_static = leaf.GetLenStatic()
    leaf_count = leaf.GetLeafCount()

    if n_static > 1 or leaf_count:
        is_static = n_static > 1
        title = leaf.GetTitle()
        # Multi-dimensional static array case: dims come from the title.
        if title.count("[") >= 2:
            dims = _get_multi_dims(title)
        else:
            dims = [leaf.GetNdata()]

        # The address is either the branch address or the leaf's value pointer.
        address = 0
        branch = leaf.GetBranch()
        if branch:
            address = helpers.branch_address(branch)
        if not address:
            address = helpers.leaf_value_pointer(leaf)

        ext_type = type_name + ("[]" if is_static else "*")
        return ROOT._cppyy._backend._create_value_from_memory(ext_type, address, dims)

    if helpers.leaf_value_pointer(leaf):
        # Value type. For TLeafElement and TLeafObject, the "value pointer" is
        # actually a pointer to the object pointer, so dereference once.
        TLeafElement = ROOT.TLeafElement
        TLeafObject = ROOT.TLeafObject
        leaf_cls = leaf.IsA()

        if leaf_cls is TLeafElement.Class() or leaf_cls is TLeafObject.Class():
            address = helpers.deref_voidp(helpers.leaf_value_pointer(leaf))
        else:
            address = helpers.leaf_value_pointer(leaf)
        return ROOT._cppyy._backend._create_value_from_memory(type_name, address)

    return None


def _TTree__getattr__(self, key):
    """
    Allow branches to be accessed as attributes of a tree.

    Allow access to branches/leaves as if they were Python data attributes of
    the tree (e.g. mytree.branch).

    Parameters:
    self (TTree): The instance of the TTree object from which the attribute is being retrieved.
    key (str): The name of the branch to retrieve from the TTree object.
    """
    import ROOT

    tree = self

    # Resolve alias
    aliased = tree.GetAlias(key)
    name = aliased if aliased else key

    # Search for the branch first (typical for object branches).
    branch = _search_for_branch(tree, name)
    if branch:
        addr, type_name = _resolve_branch(tree, name, branch)
        if type_name:
            return ROOT._cppyy.bind_object(addr, type_name)

    # Otherwise try a leaf.
    leaf = _search_for_leaf(tree, name, branch)
    if leaf:
        wrapped = _wrap_leaf(leaf)
        if wrapped is not None:
            return wrapped

    raise AttributeError(
        "'{}' object has no attribute '{}'".format(tree.IsA().GetName(), name)
    )


def _TTree_CloneTree(self, *args, **kwargs):
    """
    Forward the arguments to the C++ function and give up ownership if the
    TTree is attached to a TFile, which is the owner in that case.
    """
    import ROOT

    out_tree = self._CloneTree(*args, **kwargs)
    if _should_give_up_ownership(out_tree):
        ROOT.SetOwnership(out_tree, False)

    return out_tree


@pythonization("TTree")
def pythonize_ttree(klass, name):
    # Parameters:
    # klass: class to be pythonized
    # name: string containing the name of the class

    # Functions that need to drop the ownership if the current directory is a TFile

    klass._cpp_constructor = klass.__init__
    klass.__init__ = _constructor_releasing_ownership

    klass._CloneTree = klass.CloneTree
    klass.CloneTree = _TTree_CloneTree

    # Pythonizations that are common to TTree and its subclasses.
    # To avoid duplicating the same logic in the pythonizors of
    # the subclasses, inject the pythonizations for all the target
    # classes here.

    # Pythonic iterator
    klass.__iter__ = _TTree__iter__

    # tree.branch syntax
    klass.__getattr__ = _TTree__getattr__

    # SetBranchAddress
    klass._OriginalSetBranchAddress = klass.SetBranchAddress
    klass.SetBranchAddress = _SetBranchAddress

    # Branch
    klass._OriginalBranch = klass.Branch
    klass.Branch = _Branch

    klass._Original_SetDirectory = klass.SetDirectory
    klass.SetDirectory = _SetDirectory_SetOwnership


@pythonization("TChain")
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


@pythonization("TNtuple")
def pythonize_tntuple(klass):

    # The constructor needs to be explicitly pythonized for derived classes.
    klass._cpp_constructor = klass.__init__
    klass.__init__ = _constructor_releasing_ownership
