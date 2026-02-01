# Import the internal implementation
import ROOT._cppyy as _impl

# Re-export everything
globals().update(_impl.__dict__)

# Make introspection nicer
__all__ = getattr(_impl, "__all__", tuple(_impl.__dict__.keys()))
__doc__ = _impl.__doc__
__file__ = _impl.__file__
__path__ = getattr(_impl, "__path__", None)
__package__ = __name__

# Ensure sys.modules consistency
sys.modules[__name__] = _impl
