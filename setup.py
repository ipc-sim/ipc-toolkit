import setuptools
import os
import re
import sys
import sysconfig
import platform
import subprocess
import pathlib

from distutils.version import LooseVersion
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext


class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=""):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = pathlib.Path(sourcedir).resolve()


class CMakeBuild(build_ext):
    def run(self):
        if platform.system() == "Darwin":
            self.build_temp = self.build_temp.replace("build", "build.nosync")

        try:
            out = subprocess.check_output(["cmake", "--version"])
        except OSError:
            raise RuntimeError(
                "CMake must be installed to build the following extensions: " +
                ", ".join(e.name for e in self.extensions))

        if platform.system() == "Windows":
            cmake_version = LooseVersion(re.search(r"version\s*([\d.]+)",
                                                   out.decode()).group(1))
            if cmake_version < "3.2.0":
                raise RuntimeError("CMake >= 3.2.0 is required on Windows")

        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        extdir = pathlib.Path(self.get_ext_fullpath(
            ext.name)).resolve().parent
        cmake_args = ["-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=" + str(extdir),
                      "-DPYTHON_EXECUTABLE=" + sys.executable,
                      "-DIPC_TOOLKIT_BUILD_UNIT_TESTS=OFF",
                      "-DIPC_TOOLKIT_WITH_PYTHON=ON"]

        cfg = "Debug" if self.debug else "Release"
        build_args = ["--config", cfg]
        cmake_args += ["-DCMAKE_BUILD_TYPE=" + cfg]

        if platform.system() == "Windows":
            cmake_args += [
                "-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{}={}".format(cfg.upper(), extdir)]
            if os.environ.get("CMAKE_GENERATOR") != "NMake Makefiles":
                if sys.maxsize > 2**32:
                    cmake_args += ["-A", "x64"]
            # build_args += ["--", "/m"]
        else:
            build_args += ["--", "-j4"]

        env = os.environ.copy()
        env['CXXFLAGS'] = '{} -DVERSION_INFO=\\"{}\\"'.format(
            env.get('CXXFLAGS', ''),
            self.distribution.get_version())
        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)
        subprocess.check_call(['cmake', ext.sourcedir] + cmake_args,
                              cwd=self.build_temp, env=env)
        subprocess.check_call(['cmake', '--build', '.'] + build_args,
                              cwd=self.build_temp)
        print()  # Add an empty line for cleaner output


try:
    with open("README.md", "r") as fh:
        long_description = fh.read()
except FileNotFoundError:
    long_description = "IPC Toolkit"

setuptools.setup(
    name="ipc-toolkit",
    version="0.0.1",
    author="Zachary Ferguson",
    author_email="zfergus@nyu.edu",
    description="A set of reusable functions to integrate IPC into an existing simulation.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/ipc-sim/ipc-toolkit",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: C++",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 2.7",
        "Topic :: Games/Entertainment :: Simulation",
        "License :: OSI Approved :: MIT License",
    ],
    ext_modules=[CMakeExtension('ipctk')],
    cmdclass=dict(build_ext=CMakeBuild),
    zip_safe=False,
    install_requires=[
        'numpy',
        'scipy'
    ]
)
