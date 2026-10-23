import sys
import pathlib
import importlib.util
from importlib.machinery import PathFinder


def _is_built(spec):
    """Whether spec points at an importable ipctk, not just its Python sources.

    The pure-Python half of the package lives in python/ipctk/, and that
    directory's parent is on sys.path whenever this file is imported. Until the
    extension is built into it (e.g., by an in-place build), importing it fails
    on the missing _ipctk, so we look past it to a build directory instead.
    """
    if spec is None:
        return False
    if spec.submodule_search_locations is None:
        return True  # ipctk from before it became a package
    return any(
        any(pathlib.Path(location).glob("_ipctk.*"))
        for location in spec.submodule_search_locations
    )


# First, try to find the system-installed version
if _is_built(importlib.util.find_spec("ipctk")):
    print("Using system-installed ipctk module")
else:
    repo_root = pathlib.Path(__file__).parents[1]
    possible_paths = [
        pathlib.Path(".").resolve(),
        pathlib.Path("python").resolve(),
        repo_root / "build" / "python",
        repo_root / "build" / "python" / "python",
        repo_root / "build" / "release" / "python",
        repo_root / "build" / "debug" / "python",
    ]

    for path in possible_paths:
        if _is_built(PathFinder.find_spec("ipctk", [str(path)])):
            # We insert at the front so this build wins over the unbuilt
            # sources in python/ipctk/.
            sys.path.insert(0, str(path))
            print(f"Using found ipctk module at {path}")
            break
    else:
        raise ImportError("Could not find the ipctk module")

# We use a regular import (instead of loading the module by hand) so ipctk is
# registered in sys.modules, which `import ipctk.torch` needs.
import ipctk  # noqa: E402
