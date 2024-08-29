from .. import pythonization

from ._utils import _kwargs_to_roocmdargs


def _SPlot_init(self, *args, **kwargs):
    r"""The SPlot constructor is pythonized with the command argument pythonization.
    The keywords must correspond to the CmdArgs of the constructor.
    """

    # Pad args with default parameter values, so we can add the extra fit
    # keyword arguments at the very end.
    args = args + ([], True, False, "")[len(args) - 5 :]  # last number is first idx of default arg

    args, kwargs = _kwargs_to_roocmdargs(*args, **kwargs)
    self._init(*args, **kwargs)


def _SPlot_AddSWeight(self, *args, **kwargs):
    r"""The SPlot::AddSWeight function is pythonized with the command argument
    pythonization.
    """

    # Pad args with default parameter values, so we can add the extra fit
    # keyword arguments at the very end.
    args = args + ([], True)[len(args) - 2 :]  # last number is first idx of default arg

    args, kwargs = _kwargs_to_roocmdargs(*args, **kwargs)
    self._AddSWeight(*args, **kwargs)


@pythonization("SPlot", ns="RooStats")
def pythonize_splot(klass):
    klass._init = klass.__init__
    klass.__init__ = _SPlot_init

    klass._AddSWeight = klass.AddSWeight
    klass.AddSWeight = _SPlot_AddSWeight
