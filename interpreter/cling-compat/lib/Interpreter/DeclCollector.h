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
  class DeclCollector : public clang::ASTConsumer {
    class PPAdapter;

    const Interpreter& m_Interp;
    std::unique_ptr<EvaluateTSynthesizer> m_DynamicLookup;
    std::unique_ptr<NullDerefProtectionTransformer> m_NullDerefProtection;

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
  };
} // namespace cling

#endif // CLING_COMPAT_DECL_COLLECTOR_H
