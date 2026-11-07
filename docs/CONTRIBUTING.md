# Contributing

Contributions are welcome, and they are greatly appreciated! Every
little bit helps, and credit will always be given.

## Types of Contributions

### Report Bugs

Report bugs at <https://github.com/ipc-sim/ipc-toolkit/issues>.

If you are reporting a bug, please include:

* Your operating system name and version.
* Your compiler name and version.
* Any details about your local setup that might be helpful in troubleshooting.
* Detailed steps to reproduce the bug.

### Fix Bugs

Look through the GitHub issues for bugs. Anything tagged with "bug"
and "help wanted" is open to whoever wants to implement it.

### Implement Features

Look through the GitHub issues for features. Anything tagged with "enhancement"
and "help wanted" is open to whoever wants to implement it. Those that are
tagged with "first-timers-only" is suitable for those getting started in open-source software.

### Write Documentation

IPC Toolkit could always use more documentation, whether as part of the
official docs, in code comments, or even on the web in blog posts, articles,
and such.

### Submit Feedback

The best way to send feedback is to file an issue at <https://github.com/ipc-sim/ipc-toolkit/issues>.

If you are proposing a feature:

* Explain in detail how it would work.
* Keep the scope as narrow as possible, to make it easier to implement.
* Remember that this is a volunteer-driven project, and that contributions
  are welcome :)

## Continuous Integration

### CUDA Jobs

The CUDA workflow builds the CUDA code on Linux and Windows. On Linux it also
runs the `[gpu]` and `[cuda]` unit tests, but GitHub-hosted runners have no GPU,
so they run on a simulated NVIDIA T4
([PantheonSim](https://github.com/pantheongpu/setup-pantheonsim)). A green CUDA
job therefore means the kernels compile and give correct results, with these
limits:

* It does not measure performance.
* It skips the CCD and STQ tests on the cloth and puffer ball meshes, which take
  several minutes each on a simulated GPU.

If your change touches CUDA code, please also run the full GPU test suite on
real hardware.

## Contributors

[![Top contributors to the IPC Toolkit.](https://contrib.rocks/image?repo=ipc-sim/ipc-toolkit)](https://github.com/ipc-sim/ipc-toolkit/graphs/contributors?all=1)
