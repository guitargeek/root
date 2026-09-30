//--------------------------------------------------------------------*- C++ -*-
// cling-compat: stand-in for cling's DeclUnloader.
//
// This file is dual-licensed: you can choose to license it under the University
// of Illinois Open Source License or the GNU Lesser General Public License. See
// LICENSE.TXT for details.
//------------------------------------------------------------------------------

#ifndef CLING_COMPAT_DECL_UNLOADER_STUB_H
#define CLING_COMPAT_DECL_UNLOADER_STUB_H

#include "clang/AST/DeclBase.h"

namespace clang {
  class Sema;
}

namespace cling {
  /// \brief cling removed failed declarations (e.g. template instantiations
  /// that did not compile) from the AST again. clang-repl has no such
  /// facility for individual declarations, so they are only marked invalid.
  /// FIXME: GAP - selective AST unloading (see GAPS.md).
  inline bool UnloadDecl(clang::Sema* /*S*/, clang::Decl* D) {
    if (D)
      D->setInvalidDecl();
    return true;
  }
} // namespace cling

#endif // CLING_COMPAT_DECL_UNLOADER_STUB_H
