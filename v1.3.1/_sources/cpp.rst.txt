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
   If your :cmake:`IPC_TOOLKIT_GIT_TAG` is a tag (e.g. ``v1.3.1``), then you can use the :cmake:`FetchContent_Declare` argument :cmake:`GIT_SHALLOW TRUE` to download only a single commit. Otherwise, you should use the default :cmake:`GIT_SHALLOW FALSE`.

.. include:: ../../README.md
    :parser: myst_parser.sphinx_
    :start-after: <!--- BEGIN C++ README 2 --->
    :end-before: <!--- FILIB DEPENDENCY NOTE --->

.. warning::
    ``filib`` is licensed under `LGPL-2.1 <https://github.com/zfergus/filib/blob/main/LICENSE>`_ and as such it is required to be dynamically linked. Doing so automatically is a challenge, so by default we use static linkage. Enabling dynaic linkage requires copying the ``.so``/``.dylib``/``.dll`` file to the binary directory or system path. To enable this, set the CMake option :cmake:`FILIB_BUILD_SHARED_LIBS` to :cmake:`ON` and add this CMake code to copy the shared libaray object to the binary directory:

    .. code-block:: cmake

        # Copy shared lib to the output directory
        add_custom_command(
            TARGET ${MY_EXE_TARGET} POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_RUNTIME_DLLS:${MY_EXE_TARGET}> $<TARGET_FILE_DIR:${MY_EXE_TARGET}>
            COMMAND_EXPAND_LISTS
        )

    where ``${MY_EXE_TARGET}`` is the name of your executable target. If you know a better way to handle this, please `let us know <https://github.com/ipc-sim/ipc-toolkit/discussions>`_!

    If you would rather avoid LGPL code entirely, you can disable filib by setting :cmake:`IPC_TOOLKIT_WITH_FILIB` to :cmake:`OFF`. With this option disabled, CMake will not download or use any of filib's code.

.. include:: ../../README.md
    :parser: myst_parser.sphinx_
    :start-after: <!--- FILIB DEPENDENCY NOTE --->
    :end-before: <!--- END C++ README --->