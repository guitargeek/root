"""Compatibility shim: the cppyy package lives inside ROOT as ROOT._cppyy.

Importing 'cppyy' (or any 'cppyy.*' submodule) yields the very same module
objects as 'ROOT._cppyy' (or 'ROOT._cppyy.*'), so both import styles share a
single module universe: state is shared, and each submodule is executed only
once, no matter under which of the two names it is first imported.
"""

import importlib
import importlib.util
import sys


class _AliasLoader:
    """Loader that resolves an aliased 'cppyy.*' name to its ROOT._cppyy module."""

    def create_module(self, spec):
        return importlib.import_module("ROOT._cppyy" + spec.name[len("cppyy") :])

    def exec_module(self, module):
        # The module was already executed under its canonical name.
        pass


class _AliasFinder:
    """Meta-path finder that redirects 'cppyy.*' imports to 'ROOT._cppyy.*'."""

    def find_spec(self, fullname, path=None, target=None):
        if fullname.startswith("cppyy."):
            return importlib.util.spec_from_loader(fullname, _AliasLoader())
        return None


sys.meta_path.insert(0, _AliasFinder())

# Replace this shim with the actual implementation package, so that
# sys.modules['cppyy'] is the same module object as sys.modules['ROOT._cppyy'].
sys.modules[__name__] = importlib.import_module("ROOT._cppyy")
