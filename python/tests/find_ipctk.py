try:
    import ipctk  # Try to import the built module
except ImportError:
    import sys
    import pathlib
    sys.path.append(
        str(pathlib.Path(__file__).parents[2] / "build" / "release" / "python"))
    import ipctk  # Try again
