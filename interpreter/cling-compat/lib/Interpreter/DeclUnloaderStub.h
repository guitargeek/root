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
#include "clang/AST/DeclTemplate.h"

namespace clang {
  class Sema;
}

namespace cling {
  namespace internal {
    template <class EntryType>
    void removeSpecializationImpl(llvm::FoldingSetVector<EntryType>& Specs,
                                  const EntryType* Entry) {
      // FoldingSetVector cannot remove a node: rebuild it without Entry.
      llvm::FoldingSetVector<EntryType> Keep;
      for (auto& S : Specs) {
        if (&S != Entry) {
          S.SetNextInBucket(nullptr); // Avoid the assertion in InsertNode().
          Keep.InsertNode(&S);
        }
      }
      std::swap(Specs, Keep);
    }

    /// Access to the specializations of a template (protected), as in
    /// cling's DeclUnloader.
    class ClassTemplateDeclExt : public clang::ClassTemplateDecl {
    public:
      static void
      removeSpecialization(clang::ClassTemplateDecl* Self,
                           const clang::ClassTemplateSpecializationDecl* Spec) {
        removeSpecializationImpl(
            static_cast<ClassTemplateDeclExt*>(Self)->getCommonPtr()
                ->Specializations,
            Spec);
      }
    };

    class VarTemplateDeclExt : public clang::VarTemplateDecl {
    public:
      static void
      removeSpecialization(clang::VarTemplateDecl* Self,
                           const clang::VarTemplateSpecializationDecl* Spec) {
        removeSpecializationImpl(
            static_cast<VarTemplateDeclExt*>(Self)->getCommonPtr()
                ->Specializations,
            Spec);
      }
    };

    class FunctionTemplateDeclExt : public clang::FunctionTemplateDecl {
    public:
      static void
      removeSpecialization(clang::FunctionTemplateDecl* Self,
                           const clang::FunctionTemplateSpecializationInfo* Spec) {
        removeSpecializationImpl(
            static_cast<FunctionTemplateDeclExt*>(Self)->getCommonPtr()
                ->Specializations,
            Spec);
      }
    };
  } // namespace internal

  /// \brief Remove an implicit instantiation of a class, variable or function
  /// template from its template's specializations, so that the next use
  /// instantiates it again. The declaration itself stays in the AST. Partial
  /// and explicit specializations are left alone.
  /// \returns whether \p D was removed.
  inline bool forgetInstantiation(clang::Decl* D) {
    using namespace clang;
    if (auto* CTSD = llvm::dyn_cast_or_null<ClassTemplateSpecializationDecl>(D)) {
      if (llvm::isa<ClassTemplatePartialSpecializationDecl>(CTSD) ||
          isTemplateExplicitInstantiationOrSpecialization(
              CTSD->getSpecializationKind()))
        return false;
      CTSD = llvm::cast<ClassTemplateSpecializationDecl>(CTSD->getCanonicalDecl());
      internal::ClassTemplateDeclExt::removeSpecialization(
          CTSD->getSpecializedTemplate(), CTSD);
      return true;
    }
    if (auto* VTSD = llvm::dyn_cast_or_null<VarTemplateSpecializationDecl>(D)) {
      if (llvm::isa<VarTemplatePartialSpecializationDecl>(VTSD) ||
          isTemplateExplicitInstantiationOrSpecialization(
              VTSD->getSpecializationKind()))
        return false;
      VTSD = llvm::cast<VarTemplateSpecializationDecl>(VTSD->getCanonicalDecl());
      internal::VarTemplateDeclExt::removeSpecialization(
          VTSD->getSpecializedTemplate(), VTSD);
      return true;
    }
    if (auto* FD = llvm::dyn_cast_or_null<FunctionDecl>(D)) {
      FD = FD->getCanonicalDecl();
      FunctionTemplateSpecializationInfo* Info =
          FD->getTemplateSpecializationInfo();
      if (!Info || isTemplateExplicitInstantiationOrSpecialization(
                       Info->getTemplateSpecializationKind()))
        return false;
      internal::FunctionTemplateDeclExt::removeSpecialization(
          Info->getTemplate(), Info);
      return true;
    }
    return false;
  }

  /// \brief Forget the implicit instantiation of a class template
  /// specialization that failed (e.g. a member of incomplete type), like
  /// cling's DeclUnloader: the next use of the specialization instantiates it
  /// again, maybe successfully because the missing definitions exist by then.
  /// The declaration itself stays in the AST.
  /// \returns whether \p D was such a specialization.
  inline bool forgetFailedInstantiation(clang::Decl* D) {
    if (!D || !D->isInvalidDecl() ||
        !llvm::isa<clang::ClassTemplateSpecializationDecl>(D))
      return false;
    return forgetInstantiation(D);
  }

  /// \brief cling removed failed declarations (e.g. template instantiations
  /// that did not compile) from the AST again. clang-repl has no such
  /// facility for individual declarations: failed instantiations of class
  /// templates are forgotten (see forgetFailedInstantiation()), anything
  /// else is only marked invalid.
  /// FIXME: GAP - selective AST unloading (see GAPS.md).
  inline bool UnloadDecl(clang::Sema* /*S*/, clang::Decl* D) {
    if (D) {
      D->setInvalidDecl();
      forgetFailedInstantiation(D);
    }
    return true;
  }
} // namespace cling

#endif // CLING_COMPAT_DECL_UNLOADER_STUB_H
