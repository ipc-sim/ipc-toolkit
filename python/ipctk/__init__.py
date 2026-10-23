# The bindings live in the compiled ipctk._ipctk extension. We re-export all of
# it here, so ipctk.BarrierPotential and friends work exactly as they did when
# ipctk was a single extension module. Making ipctk a package lets us ship
# optional pure-Python submodules next to the bindings, such as ipctk.torch,
# which users import explicitly.

from . import _ipctk
from ._ipctk import *  # noqa: F401,F403
from ._ipctk import __doc__, __version__  # noqa: F401


def __getattr__(name):
    # The extension resolves deprecated names in its own module-level
    # __getattr__ (with a DeprecationWarning). A star import never copies
    # those, so we forward any name we don't define to the extension.
    return getattr(_ipctk, name)
