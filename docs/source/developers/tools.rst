Tools for Developers
=====================

Using Pre-Commit Hooks
----------------------

Use the ``.pre-commit-config.yaml`` file to apply clang-format on before commits.

Steps:
1. ``pip install pre-commit``
2. ``pre-commit install``
3. ``git add . && git commit -m "My commit message"``

Hopefully, that will avoid upsetting the GitHub checks :slightly_smiling_face:

Clang-Tidy
----------

To run locally:

1. `Install clang-tidy on mac <https://stackoverflow.com/questions/53111082/how-to-install-clang-tidy-on-macos/78243685#78243685>`_
2. ``cd build; cmake -DCMAKE_EXPORT_COMPILE_COMMANDS=ON ..``: this exports a JSON describing how all the files will be compiled
3. ``cd ..; run-clang-tidy -quiet -p build $(find src -name '*.cpp')``

You can also tack on a ``-j <N>`` to run-clang-tidy to make it parallel.
