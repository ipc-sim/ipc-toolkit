"""Interactive rigid body simulation demo (in-tree launcher).

The demo itself now ships with the package as ipctk.demo.rigid. We keep this
launcher so running it from a source checkout still works without installing:

    python rigid.py [scene.glb [settings.py]]

Once ipctk is installed, use ``python -m ipctk.demo.rigid`` instead.
"""

from find_ipctk import ipctk  # noqa: F401  (puts a built ipctk on sys.path)
from ipctk.demo.rigid import main

if __name__ == "__main__":
    main()
