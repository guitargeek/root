//--------------------------------------------------------------------*- C++ -*-
// cling-compat: LLJIT's generic IR platform, fixed for nested initialization.
//
// This file is dual-licensed: you can choose to license it under the University
// of Illinois Open Source License or the GNU Lesser General Public License. See
// LICENSE.TXT for details.
//------------------------------------------------------------------------------

#ifndef CLING_COMPAT_JIT_PLATFORM_H
#define CLING_COMPAT_JIT_PLATFORM_H

#include "llvm/ExecutionEngine/Orc/Core.h"

namespace llvm {
namespace orc {
class LLJIT;
}
} // namespace llvm

namespace cling {
  ///\brief Like llvm::orc::setUpGenericLLVMIRPlatform(), but initialize()
  /// may be called while it runs (see JITPlatform.cpp).
  llvm::Expected<llvm::orc::JITDylibSP> setUpJITPlatform(llvm::orc::LLJIT& J);
} // namespace cling

#endif // CLING_COMPAT_JIT_PLATFORM_H
