from cppyy.gbl import RooFit

from ._utils import _kwargs_to_roocmdargs

def LineColor(*args, **kwargs):
    return RooFit._LineColor(*args, **kwargs)

def FitOptions(*args,**kwargs):
    args, kwargs = _kwargs_to_roocmdargs(*args, **kwargs)
    return RooFit._FitOptions(*args, **kwargs)
