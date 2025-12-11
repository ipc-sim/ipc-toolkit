import sys
import pathlib
import importlib.util

# First, try to find the system-installed version
spec = importlib.util.find_spec("ipctk")
if spec is not None:
    ipctk = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(ipctk)
    print("Using system-installed ipctk module")
else:
    repo_root = pathlib.Path(__file__).parents[1]
    possible_paths = [
        pathlib.Path(".").resolve(),
        pathlib.Path("python").resolve(),
        repo_root / "build" / "python",
        repo_root / "build" / "release" / "python",
        repo_root / "build" / "debug" / "python",
    ]

    for path in possible_paths:
        if path.exists():
            sys.path.append(str(path))
            spec = importlib.util.find_spec("ipctk")
            if spec is not None:
                print(f"Using found ipctk module at {path}")
                ipctk = importlib.util.module_from_spec(spec)
                spec.loader.exec_module(ipctk)
                break
            else:
                sys.path.pop()  # Remove the path if ipctk is not found
    else:
        raise ImportError("Could not find the ipctk module")
