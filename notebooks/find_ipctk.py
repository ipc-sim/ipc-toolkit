try:
    import ipctk  # Try to import the built module
except ImportError:
    import sys
    import pathlib
    repo_root = pathlib.Path(__file__).parents[1]
    possible_paths = [
        pathlib.Path("python").resolve(),
        repo_root / "build" / "python",
        repo_root / "build" / "release" / "python",
        repo_root / "build" / "debug" / "python",
    ]
    for path in possible_paths:
        if path.exists() and len(list(path.glob("ipctk.*"))) > 0:
            sys.path.append(str(path))
            break
    else:
        raise ImportError("Could not find the ipctk module")
    print(f"Using found ipctk module at {path}")
    import ipctk  # Try again
