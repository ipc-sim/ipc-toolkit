"""Demo body-dynamics simulator built from IPC Toolkit components.

This package re-exports the C++ demo bindings (``Simulator``,
``KinematicDriver``, ...) and holds runnable demos, e.g.::

    python -m ipctk.demo.rigid [scene.glb [settings.py]]

The demos need extra packages: ``pip install ipctk[demo]``.
"""

from .. import _ipctk

# The C++ demo bindings are optional (IPC_TOOLKIT_BUILD_DEMO). Without them this
# package is empty, and the demos report the missing build option.
_bindings = getattr(_ipctk, "demo", None)
if _bindings is not None:
    globals().update(
        {k: v for k, v in vars(_bindings).items() if not k.startswith("_")}
    )
del _ipctk
