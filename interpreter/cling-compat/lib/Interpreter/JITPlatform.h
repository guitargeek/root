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

#include <functional>

namespace llvm {
namespace orc {
class LLJIT;
}
} // namespace llvm

namespace cling {
  ///\brief Called around running the static initializers of JIT-ed code,
  /// i.e. user code (but not around compiling it): ROOT releases its
  /// interpreter lock while user code runs.
  struct JITPlatformHooks {
    std::function<void*()> EnterUserCode;
    std::function<void(void*)> ReturnFromUserCode;
  };

  ///\brief Like llvm::orc::setUpGenericLLVMIRPlatform(), but initialize()
  /// may be called while it runs (see JITPlatform.cpp).
  llvm::Expected<llvm::orc::JITDylibSP>
  setUpJITPlatform(llvm::orc::LLJIT& J, JITPlatformHooks Hooks);
} // namespace cling

#endif // CLING_COMPAT_JIT_PLATFORM_H
