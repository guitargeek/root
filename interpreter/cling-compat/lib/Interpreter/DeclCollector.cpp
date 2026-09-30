//------------------------------------------------------------------------------
// cling-compat: records the declarations of each incremental input.
//
// This file is dual-licensed: you can choose to license it under the University
// of Illinois Open Source License or the GNU Lesser General Public License. See
// LICENSE.TXT for details.
//------------------------------------------------------------------------------

#include "DeclCollector.h"

#include "cling/Interpreter/Interpreter.h"
#include "cling/Interpreter/Transaction.h"

#include "clang/AST/Decl.h"
#include "clang/AST/DeclCXX.h"
#include "clang/Lex/PPCallbacks.h"
#include "clang/Lex/Preprocessor.h"
#include "clang/Lex/Token.h"

using namespace clang;

namespace cling {

  ///\brief Records macro (un)definitions into the current transaction.
  class DeclCollector::PPAdapter : public clang::PPCallbacks {
    const DeclCollector* m_Parent;

    void MacroDirective(const clang::Token& MacroNameTok,
                        const clang::MacroDirective* MD) {
      if (Transaction* T = m_Parent->getTransaction())
        T->append(Transaction::MacroDirectiveInfo(
            MacroNameTok.getIdentifierInfo(), MD));
    }

  public:
    PPAdapter(const DeclCollector* P) : m_Parent(P) {}

    void MacroDefined(const clang::Token& MacroNameTok,
                      const clang::MacroDirective* MD) final {
      MacroDirective(MacroNameTok, MD);
    }

    void MacroUndefined(const clang::Token& MacroNameTok,
                        const clang::MacroDefinition& /*MD*/,
                        const clang::MacroDirective* Undef) final {
      if (Undef)
        MacroDirective(MacroNameTok, Undef);
    }
  };

  DeclCollector::~DeclCollector() = default;

  void DeclCollector::Setup(clang::Preprocessor& PP) {
    PP.addPPCallbacks(std::make_unique<PPAdapter>(this));
  }

  Transaction* DeclCollector::getTransaction() const {
    if (m_Interp.m_OpenTransactions.empty())
      return nullptr;
    return m_Interp.m_OpenTransactions.back();
  }

  void DeclCollector::append(DeclGroupRef DGR, int CCI) const {
    if (DGR.isNull())
      return;
    // Decls that clang produces while no transaction is open (e.g. while
    // clang::Interpreter parses its runtime) are not tracked.
    if (Transaction* T = getTransaction())
      T->append(Transaction::DelayCallInfo(
          DGR, static_cast<Transaction::ConsumerCallInfo>(CCI)));
  }

  bool DeclCollector::HandleTopLevelDecl(DeclGroupRef DGR) {
    append(DGR, Transaction::kCCIHandleTopLevelDecl);
    return true;
  }

  void DeclCollector::HandleInterestingDecl(DeclGroupRef DGR) {
    append(DGR, Transaction::kCCIHandleInterestingDecl);
  }

  void DeclCollector::HandleTagDeclDefinition(TagDecl* TD) {
    append(DeclGroupRef(TD), Transaction::kCCIHandleTagDeclDefinition);
  }

  void DeclCollector::HandleVTable(CXXRecordDecl* RD) {
    append(DeclGroupRef(RD), Transaction::kCCIHandleVTable);
  }

  void DeclCollector::CompleteTentativeDefinition(VarDecl* VD) {
    append(DeclGroupRef(VD), Transaction::kCCICompleteTentativeDefinition);
  }

  void DeclCollector::HandleCXXImplicitFunctionInstantiation(FunctionDecl* D) {
    append(DeclGroupRef(D),
           Transaction::kCCIHandleCXXImplicitFunctionInstantiation);
  }

  void DeclCollector::HandleCXXStaticMemberVarInstantiation(VarDecl* D) {
    append(DeclGroupRef(D),
           Transaction::kCCIHandleCXXStaticMemberVarInstantiation);
  }

} // namespace cling
