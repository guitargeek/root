from ._utils import _kwargs_to_roocmdargs

RooFit = None

def _RooFit():
    """Does lazy import of RooFit namespace and then returns it.

    We can't do `from cppyy.gbl import RooFit` directly, because then we would
    see the RooFit banner every time pyROOT is initialized.

    Code inspired by:
    https://wiki.python.org/moin/PythonSpeed/PerformanceTips#Import_Statement_Overhead
    """
    global RooFit
    if RooFit is None:
        from cppyy.gbl import RooFit
    return RooFit

def LineColor(*args, **kwargs):
    return _RooFit()._LineColor(*args, **kwargs)

def FitOptions(*args,**kwargs):
    args, kwargs = _kwargs_to_roocmdargs(*args, **kwargs)
    return _RooFit()._FitOptions(*args, **kwargs)
