#include "eigen_abi.hpp"

namespace ipc {

// We define this in a .cpp file, instead of in the header, so that it reports
// the flags the library itself was compiled with. It calls this file's own copy
// of the internal-linkage eigen_abi(), not the consumer's.
EigenAbi library_eigen_abi() { return eigen_abi(); }

} // namespace ipc
