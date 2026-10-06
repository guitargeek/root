//--------------------------------------------------------------------*- C++ -*-
// cling-compat: records the declarations of each incremental input.
//
// This file is dual-licensed: you can choose to license it under the University
// of Illinois Open Source License or the GNU Lesser General Public License. See
// LICENSE.TXT for details.
//------------------------------------------------------------------------------

#ifndef CLING_COMPAT_DECL_COLLECTOR_H
#define CLING_COMPAT_DECL_COLLECTOR_H

#include "clang/AST/ASTConsumer.h"
#include "clang/AST/ASTMutationListener.h"

#include "llvm/ADT/DenseMap.h"
#include "llvm/ADT/SmallVector.h"

#include <memory>

namespace clang {
  class Preprocessor;
}

namespace cling {
  class EvaluateTSynthesizer;
  class NullDerefProtectionTransformer;
  class Interpreter;
  class Transaction;

  ///\brief Passive ASTConsumer that appends everything clang hands to the
  /// consumers into the innermost open cling::Transaction.
  ///
  /// Unlike cling's DeclCollector it does not forward to CodeGen (clang-repl
  /// multiplexes CodeGen itself, after this consumer). The only AST
  /// transformation is cling's dynamic scopes (EvaluateTSynthesizer).
  ///
  /// As an ASTMutationListener, it also records the template specializations
  /// created while a transaction is open, so that rolling the transaction back
  /// can forget them (see endTransaction()). That includes function templates:
  /// a function specialization made in the same transaction refers to the
  /// class specializations forgotten with it (e.g. in its parameter types),
  /// which the next use would no longer match.
  class DeclCollector : public clang::ASTConsumer,
                        public clang::ASTMutationListener {
    class PPAdapter;

    const Interpreter& m_Interp;
    std::unique_ptr<EvaluateTSynthesizer> m_DynamicLookup;
    std::unique_ptr<NullDerefProtectionTransformer> m_NullDerefProtection;

    ///\brief Template specializations created while each open transaction
    /// was the innermost one.
    llvm::DenseMap<const Transaction*, llvm::SmallVector<clang::Decl*, 4>>
        m_CreatedSpecializations;
    void addedSpecialization(clang::Decl* D);

    Transaction* getTransaction() const;
    bool comesFromASTReader(clang::DeclGroupRef DGR) const;
    void append(clang::DeclGroupRef DGR, int CCI) const;

  public:
    DeclCollector(const Interpreter& I);
    ~DeclCollector() override;

    ///\brief Register the macro tracking with the preprocessor.
    void Setup(clang::Preprocessor& PP);

    bool HandleTopLevelDecl(clang::DeclGroupRef DGR) override;
    void HandleInterestingDecl(clang::DeclGroupRef DGR) override;
    void HandleTagDeclDefinition(clang::TagDecl* TD) override;
    void HandleVTable(clang::CXXRecordDecl* RD) override;
    void CompleteTentativeDefinition(clang::VarDecl* VD) override;
    void HandleCXXImplicitFunctionInstantiation(clang::FunctionDecl* D) override;
    void HandleCXXStaticMemberVarInstantiation(clang::VarDecl* D) override;

    clang::ASTMutationListener* GetASTMutationListener() override {
      return this;
    }
    void AddedCXXTemplateSpecialization(
        const clang::ClassTemplateDecl* TD,
        const clang::ClassTemplateSpecializationDecl* D) override;
    void AddedCXXTemplateSpecialization(
        const clang::VarTemplateDecl* TD,
        const clang::VarTemplateSpecializationDecl* D) override;
    void AddedCXXTemplateSpecialization(const clang::FunctionTemplateDecl* TD,
                                        const clang::FunctionDecl* D) override;

    ///\brief \p T ends. If it is rolled back (\p RolledBack), forget the
    /// template instantiations it created, like cling's TransactionUnloader
    /// did: they may have been made while their arguments were incomplete
    /// (e.g. a type trait evaluated to false), and every later use would find
    /// them again.
    void endTransaction(const Transaction& T, bool RolledBack);
  };
} // namespace cling

#endif // CLING_COMPAT_DECL_COLLECTOR_H
