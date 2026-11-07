"""Check a CUDA wheel on a machine with no GPU or NVIDIA driver.

CI runs this against the wheels built by pypi.yml. The runners have no GPU, so
the wheel must still import, report that CUDA is unavailable, and raise (not
crash) when its GPU classes are used anyway.
"""
from find_ipctk import ipctk
from test_cuda import test_cuda_is_available_predicts_gpu_broad_phase

assert hasattr(ipctk.cuda, "LBVH"), "this wheel was built without CUDA"
assert not ipctk.cuda.is_available(), "expected no CUDA device on this machine"
test_cuda_is_available_predicts_gpu_broad_phase()
print("CUDA wheel imports without a GPU and raises on GPU use")
