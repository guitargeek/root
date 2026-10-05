//------------------------------------------------------------------------------
// cling-compat: port of cling's NullDerefProtectionTransformer.
//
// This file is dual-licensed: you can choose to license it under the University
// of Illinois Open Source License or the GNU Lesser General Public License. See
// LICENSE.TXT for details.
//------------------------------------------------------------------------------

#ifndef CLING_COMPAT_NULL_DEREF_PROTECTION_H
#define CLING_COMPAT_NULL_DEREF_PROTECTION_H

#include "llvm/ADT/DenseMap.h"

namespace clang {
  class Decl;
  class DirectoryEntry;
}

namespace cling {
  class Interpreter;

  /// With `--ptrcheck`, guards pointer dereferences in user code with a call
  /// to cling_runtime_internal_throwIfInvalidPointer(), which throws an
  /// InvalidDerefException (reported by ROOT) instead of crashing.
  class NullDerefProtectionTransformer {
    Interpreter& m_Interp;

    /// Whether to visit a Decl coming from a file in a given directory.
    llvm::DenseMap<const clang::DirectoryEntry*, bool> m_ShouldVisitDir;

    /// Whether the declaration should be visited and possibly transformed.
    bool shouldTransform(const clang::Decl* D);

  public:
    NullDerefProtectionTransformer(Interpreter& I);

    void Transform(clang::Decl* D);
  };
} // namespace cling

#endif // CLING_COMPAT_NULL_DEREF_PROTECTION_H
