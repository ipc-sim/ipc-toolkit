"""IPC Toolkit: reusable functions to integrate Incremental Potential Contact
(IPC) into a simulation.

The bindings live in the compiled ``ipctk._ipctk`` extension. We re-export them
here so ``import ipctk`` gives the same flat namespace as always, while the
package can also ship pure-Python modules (e.g., ``ipctk.demo.rigid``).
"""

import sys as _sys

from . import _ipctk
from ._ipctk import *  # noqa: F401,F403
from ._ipctk import __doc__  # noqa: F401


def _package_version():
    """The installed package's version, e.g., "2.0.0a0".

    The compiled module's version comes from CMake's project(VERSION), which
    only holds numbers, so it drops pre-release tags like "a0". We prefer the
    version pip installed from pyproject.toml. A different ipctk may be
    installed in the environment than the one we're running (e.g., an in-tree
    build on PYTHONPATH), so we only trust the metadata if that distribution
    installed this exact compiled module. Otherwise, we fall back to the
    compiled version.
    """
    import importlib.metadata
    import pathlib

    try:
        dist = importlib.metadata.distribution("ipctk")
    except importlib.metadata.PackageNotFoundError:
        return _ipctk.__version__
    this_module = pathlib.Path(_ipctk.__file__).resolve()
    for file in dist.files or []:
        if file.name.startswith("_ipctk") and file.parent.name == "ipctk":
            if pathlib.Path(dist.locate_file(file)).resolve() == this_module:
                return dist.version
    return _ipctk.__version__


__version__ = _package_version()
del _package_version

# The extension resolves deprecated names in a module-level __getattr__. Star
# imports skip underscore names, so we forward it explicitly.
from ._ipctk import __getattr__  # noqa: F401

# The extension registers its C++ submodules (filib, ogc, cuda, ...) in
# sys.modules under their public names, including a C++ "ipctk.demo". We drop
# that one so the import system finds our ipctk.demo package instead, which
# re-exports the C++ names alongside Python demos like ipctk.demo.rigid. We
# also drop the star-imported attribute, because `from . import demo` would
# otherwise return it without importing the package.
_sys.modules.pop("ipctk.demo", None)
globals().pop("demo", None)
from . import demo  # noqa: E402
