C++
===

.. role:: cpp(code)
   :language: c++
.. role:: cmake(code)
   :language: cmake

.. image:: https://github.com/ipc-sim/ipc-toolkit/actions/workflows/continuous.yml/badge.svg
   :target: https://github.com/ipc-sim/ipc-toolkit/actions/workflows/continuous.yml
   :alt: Build
.. image:: https://github.com/ipc-sim/ipc-toolkit/actions/workflows/docs.yml/badge.svg
   :target: https://ipctk.xyz/
   :alt: Docs
.. image:: https://codecov.io/github/ipc-sim/ipc-toolkit/graph/badge.svg?token=9BR6GPKRY8
   :target: https://codecov.io/github/ipc-sim/ipc-toolkit
.. image:: https://img.shields.io/github/license/ipc-sim/ipc-toolkit.svg?color=blue
   :target: https://github.com/ipc-sim/ipc-toolkit/blob/main/LICENSE
   :alt: License

.. include:: ../../README.md
    :parser: myst_parser.sphinx_
    :start-after: <!--- BEGIN C++ README 1 --->
    :end-before: <!--- BEGIN C++ README 2 --->

.. tip::
   If your :cmake:`IPC_TOOLKIT_GIT_TAG` is a tag (e.g. ``v1.2.0``), then you can use the :cmake:`FetchContent_Declare` argument :cmake:`GIT_SHALLOW TRUE` to download only a single commit. Otherwise, you should use the default :cmake:`GIT_SHALLOW FALSE`.

.. include:: ../../README.md
    :parser: myst_parser.sphinx_
    :start-after: <!--- BEGIN C++ README 2 --->
    :end-before: <!--- END C++ README --->