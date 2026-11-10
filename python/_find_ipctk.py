import importlib
import importlib.util
import pathlib
import sys


def _is_built(spec):
    """Whether ``spec`` is an ipctk package with its compiled _ipctk inside.

    The source tree's python/ipctk/ is the same package without the
    extension, so we must not stop at the first ``ipctk`` on sys.path.
    """
    if spec is None or not spec.submodule_search_locations:
        return False
    return any(
        any(pathlib.Path(location).glob("_ipctk*"))
        for location in spec.submodule_search_locations
    )


# First, try to find the system-installed version
if _is_built(importlib.util.find_spec("ipctk")):
    print("Using system-installed ipctk module")
else:
    repo_root = pathlib.Path(__file__).parents[1]
    possible_paths = [
        repo_root / "build" / "python" / "python",
        repo_root / "build" / "release" / "python",
        repo_root / "build" / "debug" / "python",
    ]

    for path in possible_paths:
        if path.exists():
            # We insert at the front so this wins over a source-only ipctk/
            # that may be on sys.path (e.g., when running from python/).
            sys.path.insert(0, str(path))
            if _is_built(importlib.util.find_spec("ipctk")):
                print(f"Using found ipctk module at {path}")
                break
            sys.path.pop(0)  # Remove the path if ipctk is not found
    else:
        raise ImportError("Could not find the ipctk module")

# We import through the normal machinery (not exec_module) so ipctk lands in
# sys.modules and its relative imports (ipctk._ipctk, ipctk.demo) resolve.
ipctk = importlib.import_module("ipctk")
