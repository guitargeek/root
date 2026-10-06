//------------------------------------------------------------------------------
// cling-compat: records the declarations of each incremental input.
//
// This file is dual-licensed: you can choose to license it under the University
// of Illinois Open Source License or the GNU Lesser General Public License. See
// LICENSE.TXT for details.
//------------------------------------------------------------------------------

#include "DeclCollector.h"

#include "DeclUnloaderStub.h"
#include "DynamicLookup.h"
#include "NullDerefProtection.h"

#include "cling/Interpreter/Interpreter.h"
#include "cling/Interpreter/InvocationOptions.h"
#include "cling/Interpreter/Transaction.h"
#include "cling/Utils/AST.h"

#include "clang/AST/Decl.h"
#include "clang/AST/ASTContext.h"
#include "clang/AST/DeclCXX.h"
#include "clang/AST/DeclTemplate.h"
#include "clang/AST/ExprCXX.h"
#include "clang/AST/Stmt.h"
#include "clang/Lex/PPCallbacks.h"
#include "clang/Lex/Preprocessor.h"
#include "clang/Lex/Token.h"
#include "clang/Sema/Sema.h"

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

  DeclCollector::DeclCollector(const Interpreter& I) : m_Interp(I) {}

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

  /// Unnamed macros that `return` are compiled as wrapper functions with a
  /// deduced return type (see MetaProcessor::readInputFromFile()). They may
  /// return on some paths only, as cling's void wrappers could: make falling
  /// off the end return a value-initialized result instead of being
  /// undefined.
  static void addFinalReturn(Decl* D, const Interpreter& Interp) {
    auto* FD = dyn_cast<FunctionDecl>(D);
    if (!FD || !utils::Analyze::IsWrapper(FD) ||
        Interp.getDiagnostics().hasErrorOccurred())
      return;
    auto* Body = dyn_cast_or_null<CompoundStmt>(FD->getBody());
    QualType RetTy = FD->getReturnType();
    if (!Body || !RetTy->isScalarType())
      return;
    if (!Body->body_empty() && isa<ReturnStmt>(Body->body_back()))
      return;
    ASTContext& Ctx = FD->getASTContext();
    SourceLocation Loc = Body->getRBracLoc();
    Expr* Zero = new (Ctx) CXXScalarValueInitExpr(
        RetTy, Ctx.getTrivialTypeSourceInfo(RetTy, Loc), Loc);
    llvm::SmallVector<Stmt*, 16> Stmts(Body->body());
    Stmts.push_back(ReturnStmt::Create(Ctx, Loc, Zero, nullptr));
    FD->setBody(CompoundStmt::Create(Ctx, Stmts, FPOptionsOverride(),
                                     Body->getLBracLoc(), Loc));
  }

  bool DeclCollector::HandleTopLevelDecl(DeclGroupRef DGR) {
    // Dynamic scopes: replace what uses names that are resolved at runtime
    // before CodeGen sees it. TClingCallbacks marks the functions and
    // top-level statements that need it.
    if (m_Interp.isDynamicLookupEnabled() &&
        !m_Interp.getDiagnostics().hasErrorOccurred()) {
      for (Decl* D : DGR) {
        if (!D->hasAttrs())
          continue;
        if (!m_DynamicLookup)
          m_DynamicLookup =
              std::make_unique<EvaluateTSynthesizer>(&m_Interp.getSema());
        if (!m_DynamicLookup->Transform(D))
          return false;
      }
    }
    // `--ptrcheck`: guard pointer dereferences in user code (after dynamic
    // scopes, like cling's transformers).
    // (No transaction is open while clang::Interpreter parses its runtime,
    // before the interpreter is set up.)
    Transaction* T = getTransaction();
    if (m_Interp.getOptions().PtrCheck && T &&
        T->getCompilationOpts().CheckPointerValidity &&
        !m_Interp.getDiagnostics().hasErrorOccurred()) {
      if (!m_NullDerefProtection)
        m_NullDerefProtection =
            std::make_unique<NullDerefProtectionTransformer>(
                const_cast<Interpreter&>(m_Interp));
      for (Decl* D : DGR)
        m_NullDerefProtection->Transform(D);
    }
    for (Decl* D : DGR)
      addFinalReturn(D, m_Interp);
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

  void DeclCollector::addedSpecialization(Decl* D) {
    if (const Transaction* T = getTransaction())
      m_CreatedSpecializations[T].push_back(D);
  }

  void DeclCollector::AddedCXXTemplateSpecialization(
      const ClassTemplateDecl* /*TD*/,
      const ClassTemplateSpecializationDecl* D) {
    if (!isa<ClassTemplatePartialSpecializationDecl>(D))
      addedSpecialization(const_cast<ClassTemplateSpecializationDecl*>(D));
  }

  void DeclCollector::AddedCXXTemplateSpecialization(
      const VarTemplateDecl* /*TD*/, const VarTemplateSpecializationDecl* D) {
    if (!isa<VarTemplatePartialSpecializationDecl>(D))
      addedSpecialization(const_cast<VarTemplateSpecializationDecl*>(D));
  }

  void DeclCollector::AddedCXXTemplateSpecialization(
      const FunctionTemplateDecl* /*TD*/, const FunctionDecl* D) {
    addedSpecialization(const_cast<FunctionDecl*>(D));
  }

  void DeclCollector::endTransaction(const Transaction& T, bool RolledBack) {
    auto Itr = m_CreatedSpecializations.find(&T);
    if (Itr == m_CreatedSpecializations.end())
      return;
    llvm::SmallVector<Decl*, 4> Created = std::move(Itr->second);
    m_CreatedSpecializations.erase(Itr);
    if (!RolledBack)
      return;
    // Newest first, like an unloader.
    for (Decl* D : llvm::reverse(Created))
      forgetInstantiation(D);
  }

} // namespace cling
