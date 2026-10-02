//------------------------------------------------------------------------------
// cling-compat: cling::Interpreter implemented on top of clang::Interpreter.
//
// Incremental parsing, code generation, JIT execution, undo and value
// capturing are done by clang-repl. This file only adapts cling's interface
// (transactions, callbacks, values, file and library loading) that ROOT uses.
//
// Parts are adapted from cling's Interpreter.cpp (Lukasz Janyst, Vassil
// Vassilev, Axel Naumann and others).
// This file is dual-licensed: you can choose to license it under the University
// of Illinois Open Source License or the GNU Lesser General Public License. See
// LICENSE.TXT for details.
//------------------------------------------------------------------------------

#include "cling/Interpreter/Interpreter.h"

#include "ClingUtils.h"
#include "DeclCollector.h"
#include "EnterUserCodeRAII.h"
#include "FilteringDiagConsumer.h"
#include "ForwardDeclPrinter.h"
#include "JITPlatform.h"
#include "MultiplexInterpreterCallbacks.h"

#include "cling/Interpreter/AutoloadCallback.h"
#include "cling/Interpreter/CIFactory.h"
#include "cling/Interpreter/ClangInternalState.h"
#include "cling/Interpreter/CompilationOptions.h"
#include "cling/Interpreter/DynamicLibraryManager.h"
#include "cling/Interpreter/Exception.h"
#include "cling/Interpreter/InterpreterCallbacks.h"
#include "cling/Interpreter/LookupHelper.h"
#include "cling/Interpreter/Transaction.h"
#include "cling/Interpreter/Value.h"
#include "cling/Utils/AST.h"
#include "cling/Utils/Diagnostics.h"
#include "cling/Utils/Output.h"
#include "cling/Utils/Paths.h"
#include "cling/Utils/Platform.h"
#include "cling/Utils/SourceNormalization.h"

// Compatibility.h must come first: it includes CppInterOpInterpreter.h
// once the types that header needs are declared.
#include "Compatibility.h"
#include "CppInterOpInterpreter.h"

#include "clang/AST/ASTContext.h"
#include "clang/AST/Decl.h"
#include "clang/AST/DeclCXX.h"
#include "clang/AST/GlobalDecl.h"
#include "clang/Basic/Diagnostic.h"
#include "clang/Basic/DiagnosticParse.h"
#include "clang/Basic/DiagnosticSema.h"
#include "clang/Basic/SourceManager.h"
#include "clang/CodeGen/ModuleBuilder.h"
#include "clang/Frontend/ASTConsumers.h"
#include "clang/Frontend/CompilerInstance.h"
#include "clang/Frontend/MultiplexConsumer.h"
#include "clang/Frontend/Utils.h"
#include "clang/Interpreter/CodeCompletion.h"
#include "clang/Interpreter/Interpreter.h"
#include "clang/Interpreter/PartialTranslationUnit.h"
#include "clang/Interpreter/Value.h"
#include "clang/Lex/HeaderSearch.h"
#include "clang/Lex/MacroInfo.h"
#include "clang/Lex/Preprocessor.h"
#include "clang/Parse/Parser.h"
#include "clang/Sema/Sema.h"
#include "clang/Serialization/ASTReader.h"
#include "clang/Serialization/ASTWriter.h"
#include "clang/Serialization/PCHContainerOperations.h"

#include "clang/Interpreter/IncrementalExecutor.h"

#include "llvm/ExecutionEngine/JITLink/JITLinkMemoryManager.h"
#include "llvm/ExecutionEngine/Orc/Debugging/DebuggerSupport.h"
#include "llvm/ExecutionEngine/Orc/ExecutionUtils.h"
#include "llvm/ExecutionEngine/Orc/JITTargetMachineBuilder.h"
#include "llvm/ExecutionEngine/Orc/LLJIT.h"
#include "llvm/ExecutionEngine/Orc/ObjectLinkingLayer.h"
#include "llvm/ADT/STLExtras.h"
#include "llvm/ADT/SetVector.h"
#include "llvm/ADT/SmallPtrSet.h"
#include "llvm/Passes/PassBuilder.h"
#include "llvm/Support/DynamicLibrary.h"
#include "llvm/Support/Error.h"
#include "llvm/Support/Path.h"
#include "llvm/Support/Process.h"
#include "llvm/Support/SaveAndRestore.h"
#include "llvm/Support/raw_ostream.h"

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <utility>
#include <vector>

using namespace clang;

// clang-repl's value capturing runtime (clang/lib/Interpreter). JIT-ed code
// calls these by name; clang-repl expects the host process to export them.
extern "C" {
void* __clang_Interpreter_SetValueWithAlloc(void* This, void* OutVal,
                                            void* OpaqueType);
void __clang_Interpreter_SetValueNoAlloc(void* This, void* OutVal,
                                         void* OpaqueType, ...);
int __cxa_atexit(void (*Func)(void*), void* Arg, void* DSO);
}
struct __clang_Interpreter_NewTag;
void* operator new(size_t, void*, __clang_Interpreter_NewTag) noexcept;

// The header keeps clang-repl's types opaque (see m_ClangInterpOpaque).
#define m_ClangInterp (static_cast<clang::Interpreter*>(m_ClangInterpOpaque))

namespace {
  /// Interpreters with a JIT, for deinitializeJITsAtExit().
  std::vector<cling::Interpreter*>& interpretersWithJIT() {
    // Never destroyed: it is used by an atexit handler.
    static auto* Interps = new std::vector<cling::Interpreter*>();
    return *Interps;
  }

  /// Runs the JIT's deinitializers (static destructors and atexit functions
  /// of JIT-ed code). clang-repl does this in its destructor, and compiles
  /// the JIT's platform support code (__lljit_run_atexits) for it if that
  /// was not needed before. At process exit that is too late: CodeGen uses
  /// function-local statics (e.g. SDNode::getValueTypeList()) that are
  /// destroyed in reverse order of their first use, and ROOT destroys the
  /// interpreter from TROOT's destructor. Instructions selected after their
  /// destruction refer to freed memory (random "Cannot select", "scalarize",
  /// "scalar-to-vector" errors and crashes at exit). Hence the platform code
  /// is compiled when the JIT is created and the deinitializers run from an
  /// atexit handler registered after that. Deinitializing again (in
  /// clang-repl's destructor) only runs what was registered since.
  void deinitializeJIT(cling::Interpreter& Interp) {
    llvm::orc::LLJIT* JIT = Interp.getExecutionEngine();
    if (!JIT)
      return;
    // The lookups must not reach client generators such as ROOT's library
    // autoloading.
    if (llvm::orc::JITDylib* GenJD = JIT->getExecutionSession()
            .getJITDylibByName("<cling-compat generators>"))
      JIT->getMainJITDylib().removeFromLinkOrder(*GenJD);
    if (llvm::Error Err = JIT->deinitialize(JIT->getMainJITDylib()))
      llvm::logAllUnhandledErrors(std::move(Err), llvm::errs(),
                                  "cling-compat: deinitializing the JIT: ");
  }

  /// Registered with atexit once the first JIT has compiled code; handlers
  /// run in reverse order of registration.
  void deinitializeJITsAtExit() {
    for (cling::Interpreter* Interp : interpretersWithJIT()) {
      Interp->runAtExitFuncs();
      deinitializeJIT(*Interp);
    }
  }

  /// The __cxa_atexit of JIT-ed code. With the process' own, destructors of
  /// JIT-ed static objects would only run at the very end of exit(), after
  /// ROOT's cleanups deleted what they refer to (e.g. a unique_ptr<TCanvas>
  /// whose canvas ROOT already deleted). Like cling, record them with the
  /// interpreter: runAndRemoveStaticDestructors() runs them first.
  int jitCxaAtExit(void (*Func)(void*), void* Arg, void* DSO) {
    if (interpretersWithJIT().empty())
      return ::__cxa_atexit(Func, Arg, DSO);
    interpretersWithJIT().front()->AddAtExitFunc(Func, Arg);
    return 0;
  }

  /// From cling's IncrementalJIT: a JITLinkMemoryManager that never frees
  /// its allocations. ROOT keeps pointers to JIT-ed code, e.g. TClingCallFunc's
  /// wrappers, which survive unloading the inputs after them (root-project/root
  /// issue 10898).
  class NeverFreeingJITLinkMemoryManager
      : public llvm::jitlink::InProcessMemoryManager {
  public:
    using InProcessMemoryManager::InProcessMemoryManager;

    void deallocate(std::vector<FinalizedAlloc> Allocs,
                    OnDeallocatedFunction OnDeallocated) override {
      // Release the allocations (resets their addresses, as ~FinalizedAlloc
      // asserts) without freeing them.
      for (auto& Alloc : Allocs)
        Alloc.release();
      OnDeallocated(llvm::Error::success());
    }
  };

  /// The input whose code runs on this thread (see getCurrentTransaction()).
  /// Per thread: code runs without ROOT's interpreter lock.
  cling::Transaction*& executingTransaction() {
    static thread_local cling::Transaction* T = nullptr;
    return T;
  }

  std::deque<clang::Value>* keptAlive(void* P) {
    return static_cast<std::deque<clang::Value>*>(P);
  }

  /// Print an llvm::Error the way cling reported failures, unless it only
  /// says that parsing failed (clang already printed the diagnostics).
  void reportError(llvm::Error Err, const char* Context) {
    llvm::handleAllErrors(std::move(Err), [&](const llvm::ErrorInfoBase& EIB) {
      const std::string Msg = EIB.message();
      if (llvm::StringRef(Msg).starts_with("Parsing failed"))
        return;
      cling::errs() << "cling-compat: " << Context << ": " << Msg << "\n";
    });
  }

  /// Remove trailing whitespace and semicolons: clang-repl captures the value
  /// of the last expression statement only if its semicolon is missing.
  std::string stripTrailingSemicolons(llvm::StringRef Input) {
    llvm::StringRef S = Input.rtrim();
    while (S.ends_with(";"))
      S = S.drop_back().rtrim();
    return S.str();
  }

  /// cling compiled most inputs inside a wrapper function: it inserted
  /// "void __cling_Un1Qu3N(void* vpClingValue) {\n " at \p WrapPoint. Lay
  /// out \p Source the same way, so that diagnostics keep cling's lines and
  /// columns (and reference outputs of tests match).
  std::string withClingWrapperLayout(std::string Source, size_t WrapPoint) {
    if (WrapPoint < Source.size())
      Source.insert(WrapPoint, "\n ");
    return Source;
  }

  /// Whether \p S is a single statement ending in ';' that may be an
  /// expression statement: no other ';' outside of brackets and literals,
  /// no trailing '}' and no leading statement or declaration keyword.
  bool mayBeSingleExpressionStatement(llvm::StringRef S) {
    if (!S.ends_with(";"))
      return false;
    llvm::StringRef Body = S.drop_back().rtrim();
    if (Body.empty() || Body.ends_with("}") || Body.starts_with("#") ||
        Body.starts_with("."))
      return false;
    int Depth = 0;
    for (size_t I = 0, N = Body.size(); I < N; ++I) {
      char C = Body[I];
      if (C == '"' || C == '\'') {
        for (++I; I < N && Body[I] != C; ++I)
          if (Body[I] == '\\')
            ++I;
      } else if (C == '(' || C == '[' || C == '{') {
        ++Depth;
      } else if (C == ')' || C == ']' || C == '}') {
        --Depth;
      } else if (C == ';' && Depth == 0) {
        return false;
      }
    }
    llvm::StringRef First = Body.take_while(
        [](char C) { return isalnum(C) || C == '_'; });
    static const char* const Keywords[] = {
        "if", "for", "while", "do", "switch", "return", "using", "namespace",
        "typedef", "template", "class", "struct", "union", "enum", "extern",
        "static", "goto", "break", "continue", "try", "throw", "case",
        "default", "delete", "static_assert"};
    for (const char* K : Keywords)
      if (First == K)
        return false;
    return true;
  }

  /// Forwards everything to its consumers, but HandleTranslationUnit() only
  /// when it is not called by clang-repl at the end of an incremental input:
  /// a module file must be written once, when rootcling asks for it.
  class ExplicitEndOfTUConsumer : public MultiplexConsumer {
    const bool& m_InIncrementalParse;

  public:
    ExplicitEndOfTUConsumer(std::vector<std::unique_ptr<ASTConsumer>> C,
                            const bool& InIncrementalParse)
        : MultiplexConsumer(std::move(C)),
          m_InIncrementalParse(InIncrementalParse) {}

    void HandleTranslationUnit(ASTContext& Ctx) override {
      if (!m_InIncrementalParse)
        MultiplexConsumer::HandleTranslationUnit(Ctx);
    }
  };
} // unnamed namespace

namespace cling {

  //--------------------------------------------------------------------------
  // PushTransactionRAII
  //--------------------------------------------------------------------------

  Interpreter::PushTransactionRAII::PushTransactionRAII(const Interpreter* i)
    : m_Interpreter(i) {
    CompilationOptions CO = m_Interpreter->makeDefaultCompilationOpts();
    CO.ResultEvaluation = 0;
    CO.DynamicScoping = 0;
    m_Transaction = m_Interpreter->beginTransaction(CO);
  }

  Interpreter::PushTransactionRAII::~PushTransactionRAII() {
    pop();
  }

  void Interpreter::PushTransactionRAII::pop() const {
    if (!m_Transaction)
      return;
    Transaction* T = m_Transaction;
    const_cast<PushTransactionRAII*>(this)->m_Transaction = nullptr;
    m_Interpreter->endTransaction(T,
                            T->getIssuedDiags() != Transaction::kErrors);
  }

  //--------------------------------------------------------------------------
  // Transactions
  //--------------------------------------------------------------------------

  Transaction*
  Interpreter::beginTransaction(const CompilationOptions& CO) const {
    auto* T = new Transaction(CO, getSema());
    if (!m_OpenTransactions.empty()) {
      // Remember the nesting for printing / debugging; nested transactions
      // are nevertheless committed on their own.
      T->setParent(m_OpenTransactions.back());
    }
    m_OpenTransactions.push_back(T);
    return T;
  }

  void Interpreter::endTransaction(Transaction* T, bool Commit) const {
    assert(!m_OpenTransactions.empty() && m_OpenTransactions.back() == T &&
           "Ending a transaction that is not the innermost one");
    // Like cling's commitTransaction: instantiate what was used while the
    // transaction was open. Inside Parse() clang-repl does this at the end of
    // the input; outside, pending instantiations would otherwise wait forever
    // (clang-repl's per-input instantiation scope stashes them away).
    if (Commit && !m_InIncrementalParse && m_OpenTransactions.size() == 1 &&
        !getSema().getDiagnostics().hasErrorOccurred())
      getSema().PerformPendingInstantiations();
    m_OpenTransactions.pop_back();
    T->setParent(nullptr);

    auto* This = const_cast<Interpreter*>(this);
    InterpreterCallbacks* CB = This->getCallbacks();
    if (!Commit) {
      T->setState(Transaction::kRolledBack);
      if (CB)
        CB->TransactionRollback(*T);
      delete T;
      return;
    }

    T->setState(Transaction::kCompleted);
    if (T->empty() && !T->getPTU()) {
      delete T;
      return;
    }

    T->setState(Transaction::kCommitted);
    if (!This->m_Transactions.empty())
      This->m_Transactions.back()->setNext(T);
    This->m_Transactions.emplace_back(T);
    if (CB)
      CB->TransactionCommitted(*T);
  }

  const Transaction* Interpreter::getFirstTransaction() const {
    return m_Transactions.empty() ? nullptr : m_Transactions.front().get();
  }

  const Transaction* Interpreter::getLastTransaction() const {
    return m_Transactions.empty() ? nullptr : m_Transactions.back().get();
  }

  const Transaction* Interpreter::getLastWrapperTransaction() const {
    // clang-repl does not wrap statements into functions.
    return nullptr;
  }

  const Transaction* Interpreter::getCurrentTransaction() const {
    if (!m_OpenTransactions.empty())
      return m_OpenTransactions.back();
    return executingTransaction();
  }

  const Transaction* Interpreter::getLatestTransaction() const {
    if (const Transaction* T = getCurrentTransaction())
      return T;
    return getLastTransaction();
  }

  //--------------------------------------------------------------------------
  // Construction
  //--------------------------------------------------------------------------

  Interpreter::Interpreter(int argc, const char* const* argv,
                           const char* llvmdir,
                           const ModuleFileExtensions& moduleExtensions,
                           void* extraLibHandle, bool noRuntime)
    : m_Opts(argc, argv), m_KeptAliveValues(new std::deque<clang::Value>()) {

    if (m_Opts.Help || m_Opts.ShowVersion)
      return;

    std::unique_ptr<CompilerInstance> CI =
        CIFactory::createCIForRepl(m_Opts, llvmdir, moduleExtensions);
    if (!CI)
      return;

    m_SyntaxOnly =
        CI->getFrontendOpts().ProgramAction == frontend::ParseSyntaxOnly;

    // Like cling's BackendPasses: at O0 and O1, skip the inliner passes. In
    // the default pipeline, the function simplification passes run inside
    // the inliner, so they are skipped too. JIT'd code mostly runs once; with
    // the full O1 pipeline, RDataFrame tutorials ran ~15% longer than with
    // cling.
    if (CI->getCodeGenOpts().OptimizationLevel <= 1)
      CI->getCodeGenOpts().PassBuilderCallbacks.push_back(
          [](llvm::PassBuilder& PB) {
            PB.getPassInstrumentationCallbacks()
                ->registerShouldRunOptionalPassCallback(
                    [](llvm::StringRef P, llvm::Any) {
                      return !llvm::is_contained(
                          {"ModuleInlinerWrapperPass",
                           "InlineAdvisorAnalysisPrinterPass",
                           "PartialInlinerPass", "buildInlinerPipeline",
                           "ModuleInlinerPass", "InlinerPass",
                           "InlineAdvisorAnalysis",
                           "PartiallyInlineLibCallsPass",
                           "RelLookupTableConverterPass",
                           "InlineCostAnnotationPrinterPass",
                           "InlineSizeEstimatorAnalysisPrinterPass",
                           "InlineSizeEstimatorAnalysis"},
                          P);
                    });
          });

    // Load any requested plugins. FIXME: clang-repl has no plugin support
    // (cling used to load e.g. clad here).

    // clang-repl calls this once the Preprocessor and ASTContext exist.
    auto ConsumerFactory = [this](CompilerInstance& CI)
        -> std::unique_ptr<ASTConsumer> {
      std::vector<std::unique_ptr<ASTConsumer>> Consumers;
      auto Collector = std::make_unique<DeclCollector>(*this);
      m_DeclCollector = Collector.get();
      Consumers.push_back(std::move(Collector));
      if (auto Writer = createModuleWriter(CI))
        Consumers.push_back(std::move(Writer));
      return std::make_unique<MultiplexConsumer>(std::move(Consumers));
    };

    // Like cling's IncrementalJIT: the JIT's TargetMachine must match what
    // CodeGen was told (triple, CPU, features, optimization level).
    // clang-repl's default uses the host's CPU and features instead, which
    // disagrees with the IR's function attributes.
    auto IEB = std::make_unique<clang::IncrementalExecutorBuilder>();
    if (!m_SyntaxOnly) {
      llvm::orc::JITTargetMachineBuilder JTMB(CI->getTarget().getTriple());
      JTMB.setCPU(CI->getTargetOpts().CPU);
      JTMB.addFeatures(CI->getTargetOpts().Features);
      JTMB.getOptions().MCOptions.ABIName = CI->getTarget().getABI().str();
      // LLJIT picks the small code model for JITLink if none is set. All
      // inputs go into one JITDylib, where weak symbols (like
      // DW.ref.__gxx_personality_v0, referenced from .eh_frame) are merged
      // and can end up more than 2 GB away from their users: "relocation
      // ... out of range of Delta32 fixup" and crashes when unwinding.
      // cling had the large code model (the x86-64 JIT default).
      if (CI->getTarget().getTriple().getArch() == llvm::Triple::x86_64)
        JTMB.setCodeModel(llvm::CodeModel::Large);
      switch (CI->getCodeGenOpts().OptimizationLevel) {
      case 0: JTMB.setCodeGenOptLevel(llvm::CodeGenOptLevel::None); break;
      case 1: JTMB.setCodeGenOptLevel(llvm::CodeGenOptLevel::Less); break;
      case 3: JTMB.setCodeGenOptLevel(llvm::CodeGenOptLevel::Aggressive); break;
      default: JTMB.setCodeGenOptLevel(llvm::CodeGenOptLevel::Default); break;
      }
      auto JB = std::make_unique<llvm::orc::LLJITBuilder>();
      JB->setJITTargetMachineBuilder(std::move(JTMB));
      JB->setPrePlatformSetup([](llvm::orc::LLJIT& J) {
        llvm::consumeError(llvm::orc::enableDebuggerSupport(J));
        return llvm::Error::success();
      });
      // LLJIT's generic platform, but with re-entrant initialization: ROOT
      // loads libraries while linking JIT-ed code, and their initializers
      // run more code (see JITPlatform.cpp).
      // ROOT releases its interpreter lock while user code runs, i.e. the
      // initializers (clang-repl runs an input's statements as
      // initializers), but not while it is compiled: another thread would
      // parse and generate code meanwhile, in the same LLVMContext.
      JB->setPlatformSetUp([this](llvm::orc::LLJIT& J) {
        JITPlatformHooks Hooks;
        Hooks.EnterUserCode = [this]() -> void* {
          if (InterpreterCallbacks* CB = getCallbacks())
            return CB->EnteringUserCode();
          return nullptr;
        };
        Hooks.ReturnFromUserCode = [this](void* StateInfo) {
          if (InterpreterCallbacks* CB = getCallbacks())
            CB->ReturnedFromUserCode(StateInfo);
        };
        return setUpJITPlatform(J, std::move(Hooks));
      });
      JB->setObjectLinkingLayerCreator(
          [](llvm::orc::ExecutionSession& ES)
              -> llvm::Expected<std::unique_ptr<llvm::orc::ObjectLayer>> {
            unsigned PageSize = llvm::cantFail(llvm::sys::Process::getPageSize());
            return std::make_unique<llvm::orc::ObjectLinkingLayer>(
                ES, std::make_unique<NeverFreeingJITLinkMemoryManager>(PageSize));
          });
      IEB->JITBuilder = std::move(JB);
    }

    llvm::Expected<std::unique_ptr<clang::Interpreter>> InterpOrErr = [&] {
      // The runtime clang-repl parses at startup is an incremental input too.
      llvm::SaveAndRestore<bool> InParse(m_InIncrementalParse, true);
      return clang::Interpreter::create(std::move(CI), std::move(IEB),
                                        ConsumerFactory);
    }();
    if (!InterpOrErr) {
      reportError(InterpOrErr.takeError(), "cannot create clang::Interpreter");
      m_DeclCollector = nullptr;
      return;
    }
    m_ClangInterpOpaque = InterpOrErr->get();
    m_Repl = std::make_unique<CppInternal::Interpreter>(std::move(*InterpOrErr));
    if (llvm::orc::LLJIT* JIT = getExecutionEngine()) {
      // See deinitializeJIT().
      if (auto Sym = JIT->lookup("__lljit_run_atexits"); !Sym)
        reportError(Sym.takeError(), "cannot compile the JIT's platform code");
      static const bool AtExitRegistered =
          std::atexit(&deinitializeJITsAtExit) == 0;
      (void)AtExitRegistered;
      interpretersWithJIT().push_back(this);
    }

    CompilerInstance& TheCI = *getCI();
    // Like cling: filter diagnostics of incremental inputs. The filter takes
    // over the DiagnosticsEngine's client (and its ownership).
    m_DiagConsumer.reset(new FilteringDiagConsumer(TheCI.getDiagnostics(),
                                                   /*Own=*/false));
    m_DeclCollector->Setup(TheCI.getPreprocessor());
    CIFactory::finalizeForRepl(TheCI, m_Opts);

    // Honor set of `-mllvm` options.
    if (!m_Opts.CompilerOpts.LLVMArgs.empty()) {
      std::vector<const char*> Args(m_Opts.CompilerOpts.LLVMArgs.size() + 1);
      Args[0] = "cling-compat";
      for (size_t I = 0, N = m_Opts.CompilerOpts.LLVMArgs.size(); I < N; ++I)
        Args[I + 1] = m_Opts.CompilerOpts.LLVMArgs[I].c_str();
      llvm::cl::ParseCommandLineOptions(Args.size(), Args.data());
    }

    m_LookupHelper.reset(new LookupHelper(&m_ClangInterp->getParser(), this));

    if (!m_SyntaxOnly) {
      m_DyLibManager.reset(new DynamicLibraryManager());
      m_DyLibManager->initializeDyld([](llvm::StringRef) { return false; });
      for (const std::string& P : m_Opts.LibSearchPath)
        m_DyLibManager->addSearchPath(P);

      if (llvm::orc::LLJIT* JIT = getExecutionEngine()) {
        // clang-repl's runtime lives in (and is hidden inside) libCling;
        // tell the JIT where it is.
        llvm::orc::MangleAndInterner Mangle(JIT->getExecutionSession(),
                                            JIT->getDataLayout());
        llvm::orc::SymbolMap Runtime;
        Runtime[Mangle("__cxa_atexit")] = {
            llvm::orc::ExecutorAddr::fromPtr(&jitCxaAtExit),
            llvm::JITSymbolFlags::Exported};
        Runtime[Mangle("__clang_Interpreter_SetValueWithAlloc")] = {
            llvm::orc::ExecutorAddr::fromPtr(
                &__clang_Interpreter_SetValueWithAlloc),
            llvm::JITSymbolFlags::Exported};
        Runtime[Mangle("__clang_Interpreter_SetValueNoAlloc")] = {
            llvm::orc::ExecutorAddr::fromPtr(
                &__clang_Interpreter_SetValueNoAlloc),
            llvm::JITSymbolFlags::Exported};
        using NewTagNew = void* (*)(size_t, void*, __clang_Interpreter_NewTag);
        Runtime[Mangle("_ZnwmPv26__clang_Interpreter_NewTag")] = {
            llvm::orc::ExecutorAddr::fromPtr(
                reinterpret_cast<void*>(static_cast<NewTagNew>(&operator new))),
            llvm::JITSymbolFlags::Exported};
        llvm::cantFail(JIT->getMainJITDylib().define(
            llvm::orc::absoluteSymbols(std::move(Runtime))));
      }

      // Resolve symbols from the library that hosts us (libCling is loaded
      // with RTLD_LOCAL, so its symbols are not visible to dlsym(RTLD_DEFAULT)).
      if (extraLibHandle) {
        if (llvm::orc::LLJIT* JIT = getExecutionEngine()) {
          char Prefix = JIT->getDataLayout().getGlobalPrefix();
          JIT->getMainJITDylib().addGenerator(
              std::make_unique<llvm::orc::DynamicLibrarySearchGenerator>(
                  llvm::sys::DynamicLibrary(extraLibHandle), Prefix));
        }
      }
    }

    if (getSema().getLangOpts().Modules) {
      // Explicitly create the ASTReader now; creating it later would replace
      // the ExternalSemaSources installed by the callbacks.
      if (!TheCI.getASTReader())
        TheCI.createASTReader();
      // cling's default callbacks: autoload forward-declared entities.
      setCallbacks(std::make_unique<AutoloadCallback>(this, false));
    }

    Initialize(noRuntime, m_SyntaxOnly);
  }

  void Interpreter::Initialize(bool NoRuntime, bool SyntaxOnly) {
    // Unlike cling, no wrapper and value printing runtime is needed; only
    // make cling::runtime::gCling & co. available to interpreted code.
    if (NoRuntime || !getSema().getLangOpts().CPlusPlus)
      return;

    std::string Code;
    llvm::raw_string_ostream Strm(Code);
    Strm << "#include <new>\n"
            "#include <cling/Interpreter/RuntimeUniverse.h>\n";
    if (!SyntaxOnly) {
      Strm << "namespace cling { class Interpreter; namespace runtime { "
              "Interpreter* gCling=(Interpreter*)" << (void*)this << ";\n"
              "RuntimeOptions* gClingOpts=(RuntimeOptions*)"
           << (void*)&m_RuntimeOptions << ";}}\n";
    }
    Strm.flush();
    if (declare(Code) != kSuccess)
      cling::errs() << "cling-compat: failed to initialize the runtime\n";
  }

  std::unique_ptr<ASTConsumer>
  Interpreter::createModuleWriter(CompilerInstance& CI) {
    // Adapted from cling's CIFactory (itself from GenerateModuleAction).
    const CompilerOptions& COpts = m_Opts.CompilerOpts;
    if (!COpts.CxxModules || COpts.ModuleName.empty())
      return nullptr;

    llvm::SmallString<256> Output;
    llvm::sys::path::append(Output, COpts.CachePath,
                            COpts.ModuleName + ".pcm");
    std::unique_ptr<llvm::raw_pwrite_stream> OS =
        CI.createOutputFile(Output, /*Binary=*/true,
                            /*RemoveFileOnSignal=*/false,
                            /*useTemporary=*/true,
                            /*CreateMissingDirectories=*/true);
    if (!OS)
      return nullptr;

    auto PCHBuff = std::make_shared<PCHBuffer>();
    std::vector<std::unique_ptr<ASTConsumer>> Consumers;
    Consumers.push_back(std::make_unique<PCHGenerator>(
        CI.getPreprocessor(), CI.getModuleCache(), Output, /*Sysroot=*/"",
        PCHBuff, CI.getCodeGenOpts(), CI.getFrontendOpts().ModuleFileExtensions,
        /*AllowASTWithErrors=*/false,
        /*IncludeTimestamps=*/+CI.getFrontendOpts().BuildingImplicitModule));
    Consumers.push_back(CI.getPCHContainerWriter().CreatePCHContainerGenerator(
        CI, "", Output.str().str(), std::move(OS), PCHBuff));

    // Tell clang which module we are building so that it does not build it
    // on demand when a header of this module is included.
    CI.getLangOpts().CurrentModule = COpts.ModuleName;
    CI.getLangOpts().setCompilingModule(LangOptions::CMK_ModuleMap);
    SourceManager& SM = CI.getSourceManager();
    SM.pushModuleBuildStack(COpts.ModuleName, FullSourceLoc(SourceLocation(), SM));

    return std::make_unique<ExplicitEndOfTUConsumer>(std::move(Consumers),
                                                     m_InIncrementalParse);
  }

  Interpreter::~Interpreter() {
    llvm::erase(interpretersWithJIT(), this);
    // The LookupHelper's parse buffers and the transactions refer to the AST;
    // destroy them before the clang::Interpreter.
    m_LookupHelper.reset();
    m_Callbacks.reset();
    m_DiagConsumer.reset(); // restores the original DiagnosticConsumer
    m_Transactions.clear();
    delete keptAlive(m_KeptAliveValues);
    m_KeptAliveValues = nullptr;
    // clang-repl deinitializes the JIT again when it is destroyed.
    deinitializeJIT(*this);
    m_Repl.reset();
  }

  bool Interpreter::isValid() const { return m_ClangInterp; }

  const char* Interpreter::getVersion() { return "cling-compat (clang-repl)"; }

  bool Interpreter::isInSyntaxOnlyMode() const { return m_SyntaxOnly; }

  CompilerInstance* Interpreter::getCI() const {
    return m_ClangInterp ? m_ClangInterp->getCompilerInstance() : nullptr;
  }

  CompilerInstance* Interpreter::getCIOrNull() const { return getCI(); }

  Sema& Interpreter::getSema() const { return getCI()->getSema(); }

  DiagnosticsEngine& Interpreter::getDiagnostics() const {
    return getCI()->getDiagnostics();
  }

  const Parser& Interpreter::getParser() const {
    return m_ClangInterp->getParser();
  }

  Parser& Interpreter::getParser() { return m_ClangInterp->getParser(); }

  SourceLocation Interpreter::getNextAvailableLoc() const {
    // Every call must return a new location, like cling's
    // getNextAvailableUniqueSourceLoc(): clang caches e.g. the last module
    // import by its location (CompilerInstance::loadModule), and the
    // LookupHelper parses each snippet "included" from such a location.
    static constexpr unsigned BufferSize = 1 << 16;
    SourceManager& SM = getCI()->getSourceManager();
    if (!m_UniqueLocBufferStart || m_UniqueLocOffset >= BufferSize - 1) {
      static const std::string Spaces(BufferSize, ' ');
      FileID FID = SM.createFileID(
          llvm::MemoryBuffer::getMemBuffer(Spaces, "<cling-compat locations>",
                                           /*RequiresNullTerminator=*/false),
          SrcMgr::C_User);
      m_UniqueLocBufferStart = SM.getLocForStartOfFile(FID).getRawEncoding();
      m_UniqueLocOffset = 0;
    }
    return SourceLocation::getFromRawEncoding(m_UniqueLocBufferStart)
        .getLocWithOffset(m_UniqueLocOffset++);
  }

  SourceLocation Interpreter::getSourceLocation(bool /*skipWrapper*/) const {
    const Transaction* T = getLatestTransaction();
    if (!T)
      return SourceLocation();
    return T->getSourceStart(getCI()->getSourceManager());
  }

  void Interpreter::replaceDiagnosticConsumer(DiagnosticConsumer* Consumer,
                                              bool Own) {
    static_cast<FilteringDiagConsumer&>(*m_DiagConsumer)
        .setTargetConsumer(Consumer, Own);
  }

  bool Interpreter::hasReplacedDiagnosticConsumer() const {
    return static_cast<FilteringDiagConsumer&>(*m_DiagConsumer)
               .getTargetConsumer() != nullptr;
  }

  CompilationOptions Interpreter::makeDefaultCompilationOpts() const {
    CompilationOptions CO;
    CO.DeclarationExtraction = 0;
    CO.EnableShadowing = 0;
    CO.ValuePrinting = CompilationOptions::VPDisabled;
    CO.CodeGeneration = !isInSyntaxOnlyMode();
    CO.DynamicScoping = isDynamicLookupEnabled();
    CO.Debug = isPrintingDebug();
    CO.IgnorePromptDiags = 0;
    CO.CheckPointerValidity = !isRawInputEnabled();
    CO.OptLevel = getDefaultOptLevel();
    return CO;
  }

  //--------------------------------------------------------------------------
  // Include paths, macros, state dumping
  //--------------------------------------------------------------------------

  void Interpreter::AddIncludePaths(llvm::StringRef PathStr, const char* Delm) {
    CompilerInstance* CI = getCI();
    HeaderSearchOptions& HOpts = CI->getHeaderSearchOpts();

    // Save the current number of entries
    size_t Idx = HOpts.UserEntries.size();
    utils::AddIncludePaths(PathStr, HOpts, Delm);

    Preprocessor& PP = CI->getPreprocessor();
    SourceManager& SM = PP.getSourceManager();
    FileManager& FM = SM.getFileManager();
    HeaderSearch& HSearch = PP.getHeaderSearchInfo();
    const bool isFramework = false;

    // Add all the new entries into Preprocessor
    for (const size_t N = HOpts.UserEntries.size(); Idx < N; ++Idx) {
      const HeaderSearchOptions::Entry& E = HOpts.UserEntries[Idx];
      if (auto DE = FM.getOptionalDirectoryRef(E.Path))
        HSearch.AddSearchPath(DirectoryLookup(*DE, SrcMgr::C_User, isFramework),
                              E.Group == frontend::Angled);
    }
  }

  void Interpreter::AddIncludePath(llvm::StringRef PathsStr) {
    return AddIncludePaths(PathsStr, nullptr);
  }

  void Interpreter::GetIncludePaths(llvm::SmallVectorImpl<std::string>& incpaths,
                                    bool withSystem, bool withFlags) {
    utils::CopyIncludePaths(getCI()->getHeaderSearchOpts(), incpaths,
                            withSystem, withFlags);
  }

  void Interpreter::DumpIncludePath(llvm::raw_ostream* S) const {
    utils::DumpIncludePaths(getCI()->getHeaderSearchOpts(),
                            S ? *S : cling::outs(),
                            true /*withSystem*/, true /*withFlags*/);
  }

  void Interpreter::DumpDynamicLibraryInfo(llvm::raw_ostream* S) const {
    if (auto DLM = getDynamicLibraryManager())
      DLM->dump(S);
  }

  void Interpreter::dump(llvm::StringRef what, llvm::StringRef filter) {
    llvm::raw_ostream &where = cling::log();
    PushTransactionRAII RAII(this);
    if (what == "asttree") {
      std::unique_ptr<clang::ASTConsumer> printer =
        clang::CreateASTDumper(nullptr /*Dump to stdout.*/,
                               filter, true  /*DumpDecls*/,
                               false /*Deserialize*/,
                               false /*DumpLookups*/,
                               false /*DumpDeclTypes*/,
                               ADOF_Default /*DumpFormat*/);
      printer->HandleTranslationUnit(getSema().getASTContext());
    } else if (what == "ast")
      getSema().getASTContext().PrintStats();
    else if (what == "decl")
      ClangInternalState::printLookupTables(where, getSema().getASTContext());
    else if (what == "undo") {
      for (const auto& T : m_Transactions)
        T->printStructureBrief();
    }
  }

  void Interpreter::storeInterpreterState(const std::string& /*name*/) const {
    // FIXME: cling-compat: .storeState needs cling's ClangInternalState
    // bookkeeping of llvm::Modules, which clang-repl owns.
    cling::errs() << "cling-compat: .storeState is not supported\n";
  }

  void Interpreter::compareInterpreterState(const std::string& /*name*/) const {
    cling::errs() << "cling-compat: .compareState is not supported\n";
  }

  void Interpreter::printIncludedFiles(llvm::raw_ostream& Out) const {
    ClangInternalState::printIncludedFiles(Out, getCI()->getSourceManager());
  }

  const MacroInfo* Interpreter::getMacro(llvm::StringRef Macro) const {
    clang::Preprocessor& PP = getCI()->getPreprocessor();
    if (IdentifierInfo* II = PP.getIdentifierInfo(Macro)) {
      if (II->isOutOfDate())
        PP.getExternalSource()->updateOutOfDateIdentifier(*II);
      MacroDefinition MDef = PP.getMacroDefinition(II);
      return MDef.getMacroInfo();
    }
    return nullptr;
  }

  std::string Interpreter::getMacroValue(llvm::StringRef Macro,
                                         const char* Trim) const {
    std::string Value;
    if (const MacroInfo* MI = getMacro(Macro)) {
      for (const clang::Token& Tok : MI->tokens()) {
        llvm::SmallString<64> Buffer;
        Macro = getCI()->getPreprocessor().getSpelling(Tok, Buffer);
        if (!Value.empty())
          Value += " ";
        Value += Trim ? Macro.trim(Trim).str() : Macro.str();
      }
    }
    return Value;
  }

  void Interpreter::createUniqueName(std::string& Out) {
    llvm::raw_string_ostream Strm(Out);
    Strm << utils::Synthesize::UniquePrefix << m_UniqueCounter++;
  }

  bool Interpreter::isUniqueName(llvm::StringRef name) {
    return name.starts_with(utils::Synthesize::UniquePrefix);
  }

  //--------------------------------------------------------------------------
  // Compilation and execution
  //--------------------------------------------------------------------------

  void Interpreter::convertValue(void* FromOpaque, Value& To) {
    clang::Value& From = *static_cast<clang::Value*>(FromOpaque);
    if (!From.isValid()) {
      To = Value();
      return;
    }
    QualType QT = From.getType();
    if (From.isVoid() || QT->isVoidType()) {
      To = Value(getCI()->getASTContext().VoidTy, *this);
      return;
    }

    switch (From.getKind()) {
#define X(type, name)                                                          \
    case clang::Value::K_##name: {                                             \
      To = Value(QT, *this);                                                   \
      type Val = From.get##name();                                             \
      std::memcpy(To.getPtrAddress(), &Val, sizeof(Val));                      \
      return;                                                                  \
    }
    REPL_BUILTIN_TYPES
#undef X
    case clang::Value::K_PtrOrObj: {
      QualType Canon = QT.getCanonicalType();
      if (getenv("CLING_COMPAT_DEBUG"))
        llvm::errs() << "cling-compat: convertValue " << QT.getAsString()
                     << " ptr " << From.getPtr() << " manual "
                     << From.isManuallyAlloc() << "\n";
      if (Canon->isRecordType() || Canon->isConstantArrayType() ||
          Canon->isMemberPointerType()) {
        // The object lives in storage owned by the clang::Value. cling::Value
        // cannot own it, so keep the clang::Value alive.
        // FIXME: this leaks until the interpreter dies.
        keptAlive(m_KeptAliveValues)->push_back(std::move(From));
        To = Value::MakeNonOwning(QT, *this, keptAlive(m_KeptAliveValues)->back().getPtr());
        return;
      }
      To = Value(QT, *this);
      *To.getPtrAddress() = From.getPtr();
      return;
    }
    default:
      To = Value();
      return;
    }
  }

  Interpreter::CompilationResult
  Interpreter::ParseAndExecute(const std::string& input,
                               const CompilationOptions& CO, Value* CV,
                               bool print, Transaction** OutT) {
    if (OutT)
      *OutT = nullptr;
    if (CV)
      *CV = Value();
    if (!m_ClangInterp)
      return kFailure;

    Transaction* T = beginTransaction(CO);
    llvm::Expected<clang::PartialTranslationUnit&> PTUOrErr = [&] {
      llvm::SaveAndRestore<bool> InParse(m_InIncrementalParse, true);
      FilteringDiagConsumer::RAAI IgnorePromptDiags(*m_DiagConsumer,
                                                    CO.IgnorePromptDiags);
      // See process(). Only for this input, not for nested ones.
      bool PrevCapture = m_ClangInterp->capturesValueOfExprStmts();
      m_ClangInterp->setCaptureValueOfExprStmts(
          std::exchange(m_CaptureValueOfExprStmt, false));
      auto Result = m_ClangInterp->Parse(input);
      m_ClangInterp->setCaptureValueOfExprStmts(PrevCapture);
      return Result;
    }();
    if (!PTUOrErr) {
      reportError(PTUOrErr.takeError(), "failed to parse");
      endTransaction(T, /*Commit=*/false);
      return kFailure;
    }
    clang::PartialTranslationUnit& PTU = *PTUOrErr;
    T->setPTU(&PTU);
    const bool HasDecls = !T->empty();
    // Commit before running the static initializers, as cling did: ROOT
    // needs to know about new classes and globals while they run.
    endTransaction(T, /*Commit=*/true);
    if (OutT && HasDecls)
      *OutT = T;

    if (isInSyntaxOnlyMode() || !PTU.TheModule)
      return kSuccess;

    {
      // The JIT platform releases ROOT's lock while the code runs (see
      // setUpJITPlatform()). Like cling, the input stays the current
      // transaction while it runs:
      // e.g. TCling::AutoParse() ties the payloads it parses to it, so that
      // undoing the input allows parsing them again.
      llvm::SaveAndRestore<Transaction*> Executing(executingTransaction(), T);
      if (llvm::Error Err = m_ClangInterp->Execute(PTU)) {
        reportError(std::move(Err), "failed to execute");
        return kFailure;
      }
    }

    clang::Value LastValue;
    m_ClangInterp->takeLastValue(LastValue);
    if (!LastValue.isValid() || (!CV && !print))
      return kSuccess;

    Value Result;
    convertValue(&LastValue, Result);
    if (print && Result.isValid() && !Result.isVoid())
      Result.dump();
    if (CV)
      *CV = std::move(Result);
    return kSuccess;
  }

  Interpreter::CompilationResult
  Interpreter::process(const std::string& input, Value* V /* = 0 */,
                       Transaction** T /* = 0 */,
                       bool disableValuePrinting /* = false*/) {
    CompilationOptions CO = makeDefaultCompilationOpts();
    if (isRawInputEnabled())
      return ParseAndExecute(input, CO, nullptr, /*print=*/false, T);

    CO.ValuePrinting = disableValuePrinting ? CompilationOptions::VPDisabled
                                            : CompilationOptions::VPAuto;
    CO.ResultEvaluation = (bool)V;
    // Like cling's EvaluateInternal: no "unused expression" noise.
    CO.IgnorePromptDiags = 1;

    // cling accepts a declaration without trailing ';' at the prompt and
    // prints the declared variable (`auto h = new TH1F(...)`); clang-repl
    // only allows that for expression statements. Try as-is (keeping the
    // diagnostics), and on "expected ';' after top level declarator" retry
    // with the semicolon, then print the variable.
    // FIXME: GAP - clang-repl could support this directly.
    llvm::StringRef Trimmed = llvm::StringRef(input).trim();
    // TCling prefixes prompt input with `#line 1 "ROOT_prompt_N"`.
    while (Trimmed.starts_with("#line"))
      Trimmed = Trimmed.split('\n').second.trim();

    // cling also returned the value of a last expression ending in ';'
    // (without printing it): TCling::ProcessLine("&f;") returns the address
    // of f, e.g. for RooFit's code generation. clang-repl captures a value
    // only if the ';' is missing; ask it to capture it anyway. (Dropping the
    // ';' made declarations like `std::map<int,std::string> p;` fail to
    // parse once, and the retry lacked the template instantiations made by
    // the failed attempt.)
    llvm::SaveAndRestore<bool> CaptureValue(m_CaptureValueOfExprStmt, false);
    if (V && mayBeSingleExpressionStatement(Trimmed)) {
      m_CaptureValueOfExprStmt = true;
      disableValuePrinting = true;
      CO.ValuePrinting = CompilationOptions::VPDisabled;
    }
    // Where cling would have started its wrapper function (see
    // withClingWrapperLayout()).
    std::string Input = input;
    size_t WrapPoint = utils::getWrapPoint(Input, getSema().getLangOpts());
    Input = withClingWrapperLayout(std::move(Input), WrapPoint);
    // cling wrapped statements into a function, where `return` is valid (and
    // ROOT uses `root -e return`). clang-repl has no enclosing function.
    // FIXME: GAP - `return` in top-level statements (e.g. unnamed macros).
    if (Trimmed == "return" || Trimmed == "return;")
      return kSuccess;
    if (Trimmed.empty() || Trimmed.ends_with(";") || Trimmed.ends_with("}") ||
        Trimmed.starts_with("#"))
      return ParseAndExecute(Input, CO, V, !disableValuePrinting, T);

    /// Forwards all diagnostics but the missing ';' one, which it records.
    class MissingSemiDetector : public DiagnosticConsumer {
      DiagnosticsEngine& m_Diags;
      DiagnosticConsumer* m_Prev;
      bool m_OwnsPrev;
    public:
      bool Seen = false;
      MissingSemiDetector(DiagnosticsEngine& Diags)
          : m_Diags(Diags), m_Prev(Diags.getClient()),
            m_OwnsPrev(Diags.ownsClient()) {
        if (m_OwnsPrev)
          Diags.takeClient().release();
        Diags.setClient(this, /*ShouldOwnClient=*/false);
      }
      ~MissingSemiDetector() override {
        m_Diags.setClient(m_Prev, m_OwnsPrev);
      }
      void HandleDiagnostic(DiagnosticsEngine::Level Level,
                            const Diagnostic& Info) override {
        if (Info.getID() == diag::err_invalid_token_after_toplevel_declarator) {
          Seen = true;
          DiagnosticConsumer::HandleDiagnostic(Level, Info); // count it
          return;
        }
        DiagnosticConsumer::HandleDiagnostic(Level, Info);
        m_Prev->HandleDiagnostic(Level, Info);
      }
      void clear() override {
        DiagnosticConsumer::clear();
        m_Prev->clear();
      }
      void BeginSourceFile(const LangOptions& LO, const Preprocessor* PP) override {
        m_Prev->BeginSourceFile(LO, PP);
      }
      void EndSourceFile() override { m_Prev->EndSourceFile(); }
      void finish() override { m_Prev->finish(); }
    };

    CompilationResult Res;
    bool MissingSemi = false;
    {
      MissingSemiDetector Detector(getDiagnostics());
      Res = ParseAndExecute(Input, CO, V, !disableValuePrinting, T);
      MissingSemi = Res == kFailure && Detector.Seen;
    }
    if (!MissingSemi)
      return Res;

    Transaction* DeclT = nullptr;
    CompilationOptions DeclCO = CO;
    DeclCO.ValuePrinting = CompilationOptions::VPDisabled;
    Res = ParseAndExecute(Input + ";", DeclCO, nullptr, /*print=*/false,
                          &DeclT);
    if (T)
      *T = DeclT;
    if (Res != kSuccess || !DeclT)
      return Res;

    const VarDecl* LastVar = nullptr;
    for (auto I = DeclT->decls_begin(), E = DeclT->decls_end(); I != E; ++I)
      if (I->m_Call == Transaction::kCCIHandleTopLevelDecl)
        for (Decl* D : I->m_DGR)
          if (auto* VD = dyn_cast<VarDecl>(D))
            LastVar = VD;
    if (!LastVar || (disableValuePrinting && !V))
      return Res;
    return ParseAndExecute(LastVar->getQualifiedNameAsString(), CO, V,
                           !disableValuePrinting, nullptr);
  }

  Interpreter::CompilationResult
  Interpreter::parse(const std::string& input, Transaction** T) const {
    // FIXME: clang-repl always generates code; the PTU is simply not executed.
    CompilationOptions CO = makeDefaultCompilationOpts();
    CO.CodeGeneration = 0;
    Transaction* Tr = beginTransaction(CO);
    llvm::Expected<clang::PartialTranslationUnit&> PTUOrErr = [&] {
      llvm::SaveAndRestore<bool> InParse(m_InIncrementalParse, true);
      return m_ClangInterp->Parse(input);
    }();
    if (!PTUOrErr) {
      reportError(PTUOrErr.takeError(), "failed to parse");
      endTransaction(Tr, false);
      return kFailure;
    }
    Tr->setPTU(&*PTUOrErr);
    bool HasDecls = !Tr->empty();
    endTransaction(Tr, true);
    if (T)
      *T = HasDecls ? Tr : nullptr;
    return kSuccess;
  }

  Interpreter::CompilationResult
  Interpreter::declare(const std::string& input, Transaction** T/*=0 */) {
    CompilationOptions CO = makeDefaultCompilationOpts();
    CO.CheckPointerValidity = 0;
    return ParseAndExecute(input, CO, nullptr, /*print=*/false, T);
  }

  Interpreter::CompilationResult
  Interpreter::evaluate(const std::string& input, Value& V) {
    CompilationOptions CO = makeDefaultCompilationOpts();
    CO.ResultEvaluation = 1;
    CO.CheckPointerValidity = 0;
    CO.IgnorePromptDiags = 1;
    return ParseAndExecute(
        withClingWrapperLayout(stripTrailingSemicolons(input), 0), CO, &V,
        /*print=*/false);
  }

  Interpreter::CompilationResult
  Interpreter::echo(const std::string& input, Value* V /* = 0 */) {
    CompilationOptions CO = makeDefaultCompilationOpts();
    CO.ValuePrinting = CompilationOptions::VPEnabled;
    CO.ResultEvaluation = (bool)V;
    CO.IgnorePromptDiags = 1;
    return ParseAndExecute(
        withClingWrapperLayout(stripTrailingSemicolons(input), 0), CO, V,
        /*print=*/true);
  }

  Interpreter::CompilationResult
  Interpreter::execute(const std::string& input) {
    CompilationOptions CO = makeDefaultCompilationOpts();
    CO.DynamicScoping = 0;
    return ParseAndExecute(withClingWrapperLayout(input, 0), CO, nullptr,
                           /*print=*/false);
  }

  Interpreter::CompilationResult
  Interpreter::parseForModule(const std::string& input) {
    CompilationOptions CO = makeDefaultCompilationOpts();
    CO.CodeGenerationForModule = 1;

    DiagnosticsEngine& Diag = getDiagnostics();
    Diag.setSeverity(clang::diag::warn_field_is_uninit,
                     clang::diag::Severity::Ignored, SourceLocation());
    CompilationResult Result = ParseAndExecute(input, CO, nullptr, false);
    Diag.setSeverity(clang::diag::warn_field_is_uninit,
                     clang::diag::Severity::Warning, SourceLocation());
    return Result;
  }

  Interpreter::CompilationResult Interpreter::emitAllDecls(Transaction* T) {
    // Declarations created by hand (e.g. ROOT's special objects) must reach
    // CodeGen. Hand them to clang-repl's consumers; the code lands in the
    // module of the input being parsed (or the next one).
    std::unique_ptr<Transaction> Owner(T);
    if (isInSyntaxOnlyMode())
      return kSuccess;
    ASTConsumer& Consumer = getSema().getASTConsumer();
    for (auto I = T->decls_begin(), E = T->decls_end(); I != E; ++I)
      if (I->m_Call == Transaction::kCCIHandleTopLevelDecl)
        Consumer.HandleTopLevelDecl(I->m_DGR);
    return kSuccess;
  }

  Interpreter::CompilationResult
  Interpreter::codeComplete(const std::string& line, size_t& cursor,
                            std::vector<std::string>& completions) const {
    // FIXME: cling created a second CompilerInstance through CIFactory for
    // clang::ReplCodeCompleter. Not yet ported.
    (void)line; (void)cursor; (void)completions;
    return kFailure;
  }

  bool Interpreter::loadModule(const std::string& moduleName,
                               bool complain /*= true*/) {
    assert(getCI()->getLangOpts().Modules
           && "Function only relevant when C++ modules are turned on!");

    Preprocessor& PP = getCI()->getPreprocessor();
    HeaderSearch &HS = PP.getHeaderSearchInfo();

    if (Module *M = HS.lookupModule(moduleName, SourceLocation(),
                                    /*AllowSearch*/true,
                                    /*AllowExtraSearch*/ true))
      return loadModule(M, complain);

    if (complain)
      llvm::errs() << "Module " << moduleName << " not found.\n";

    return false;
  }

  bool Interpreter::loadModule(clang::Module* M, bool complain /* = true*/) {
    assert(getCI()->getLangOpts().Modules
           && "Function only relevant when C++ modules are turned on!");
    assert(M && "Module missing");
    if (getSema().isModuleVisible(M))
      return true;

    Preprocessor& PP = getCI()->getPreprocessor();
    IdentifierInfo* II = PP.getIdentifierInfo(M->Name);
    SourceLocation ValidLoc = getNextAvailableLoc();

    // Don't push the ImportDecl into a ClassDecl or whatever - who knows which
    // context triggered this import!
    DeclContext *OldDC = getSema().CurContext;
    getSema().CurContext = getSema().getASTContext().getTranslationUnitDecl();

    auto Path = IdentifierLoc(ValidLoc, II);
    Module* Mod = getSema().getModuleLoader().loadModule(
        ValidLoc, Path, Module::AllVisible, /*IsInclusionDirective=*/false);
    bool success = Mod && !getSema()
                               .ActOnModuleImport(ValidLoc, /*ExportLoc*/ {},
                                                  ValidLoc, Mod, Path)
                               .isInvalid();
    getSema().CurContext = OldDC;

    if (getenv("CLING_COMPAT_DEBUG"))
      llvm::errs() << "cling-compat: loadModule(" << M->Name << "): loaded "
                   << (Mod != nullptr) << ", success " << success
                   << ", import loc valid " << ValidLoc.isValid()
                   << ", visible " << getSema().isModuleVisible(M)
                   << ", same module " << (Mod == M) << ", Mod visible "
                   << (Mod && getSema().isModuleVisible(Mod)) << ", visID "
                   << M->getVisibilityID() << "/"
                   << (Mod ? Mod->getVisibilityID() : 0) << "\n";

    if (success) {
      // Also make the module visible in the preprocessor to export its macros.
      PP.makeModuleVisible(M, ValidLoc);
      return success;
    }

    if (complain) {
      if (M->IsSystem)
        llvm::errs() << "Failed to load module " << M->Name << "\n";
      else
        llvm::outs() << "Failed to load module " << M->Name << "\n";
    }

    return false;
  }

  //--------------------------------------------------------------------------
  // Files and libraries
  //--------------------------------------------------------------------------

  std::string Interpreter::lookupFileOrLibrary(llvm::StringRef file) {
    std::string canonicalFile = DynamicLibraryManager::normalizePath(file);
    if (canonicalFile.empty())
      canonicalFile = file.str();

    bool isAngled = false;
    ConstSearchDirIterator FromDir = nullptr;
    const FileEntry* FromFile = nullptr;
    ConstSearchDirIterator* CurDir = nullptr;
    Preprocessor& PP = getCI()->getPreprocessor();
    SourceLocation fileNameLoc;
    auto FE = PP.LookupFile(fileNameLoc, canonicalFile, isAngled, FromDir,
                            FromFile, CurDir, /*SearchPath*/nullptr,
                            /*RelativePath*/ nullptr, /*suggestedModule*/nullptr,
                            /*IsMapped*/ nullptr, /*IsFrameworkFound*/ nullptr,
                            /*SkipCache*/ false, /*OpenFile*/ false,
                            /*CacheFail*/ false);
    if (FE)
      return FE->getName().str();
    if (auto* DLM = getDynamicLibraryManager())
      return DLM->lookupLibrary(canonicalFile);
    return "";
  }

  Interpreter::CompilationResult
  Interpreter::loadLibrary(const std::string& filename, bool lookup) {
    DynamicLibraryManager* DLM = getDynamicLibraryManager();
    if (!DLM)
      return kFailure;
    std::string canonicalLib;
    if (lookup)
      canonicalLib = DLM->lookupLibrary(filename);

    const std::string &library = lookup ? canonicalLib : filename;
    if (!library.empty()) {
      switch (DLM->loadLibrary(library, /*permanent*/false, /*resolved*/true)) {
      case DynamicLibraryManager::kLoadLibSuccess: // Intentional fall through
      case DynamicLibraryManager::kLoadLibAlreadyLoaded:
        return kSuccess;
      case DynamicLibraryManager::kLoadLibNotFound:
        assert(0 && "Cannot find library with existing canonical name!");
        return kFailure;
      default:
        // Not a source file (canonical name is non-empty) but can't load.
        return kFailure;
      }
    }
    return kMoreInputExpected;
  }

  Interpreter::CompilationResult
  Interpreter::loadHeader(const std::string& filename,
                          Transaction** T /*= 0*/) {
    std::string code;
    code += "#include \"" + filename + "\"";

    CompilationOptions CO = makeDefaultCompilationOpts();
    CO.CheckPointerValidity = 1;
    return ParseAndExecute(code, CO, nullptr, /*print=*/false, T);
  }

  Interpreter::CompilationResult
  Interpreter::loadFile(const std::string& filename,
                        bool allowSharedLib /*=true*/,
                        Transaction** T /*= 0*/) {
    if (allowSharedLib) {
      CompilationResult result = loadLibrary(filename, true);
      if (result != kMoreInputExpected)
        return result;
    }
    return loadHeader(filename, T);
  }

  /// Like cling's TransactionUnloader: undefine the macros that \p T defined
  /// and forget the contents of the files it read. Loading an unloaded file
  /// again then reads it again (instead of skipping it because its include
  /// guard is still defined), and sees its changes.
  static void unloadFromPreprocessor(Transaction& T, Sema& S) {
    Preprocessor& PP = S.getPreprocessor();
    SourceManager& SM = S.getSourceManager();
    llvm::SmallSetVector<FileID, 4> Files;
    auto Collect = [&](SourceLocation Loc) {
      if (Loc.isInvalid())
        return;
      FileID FID = SM.getFileID(SM.getSpellingLoc(Loc));
      if (FID.isInvalid() || FID == SM.getMainFileID())
        return;
      // Not the inputs themselves (included from the main file).
      SourceLocation IncludeLoc = SM.getIncludeLoc(FID);
      if (IncludeLoc.isInvalid() ||
          SM.getFileID(IncludeLoc) == SM.getMainFileID())
        return;
      Files.insert(FID);
    };
    for (auto MI = T.rmacros_begin(), ME = T.rmacros_end(); MI != ME; ++MI) {
      Collect(MI->m_MD->getLocation());
      // An #undef has no definition to remove.
      if (MI->m_MD->getMacroInfo())
        PP.removeMacro(MI->m_II, const_cast<MacroDirective*>(MI->m_MD));
    }
    for (auto I = T.decls_begin(), E = T.decls_end(); I != E; ++I)
      for (Decl* D : I->m_DGR)
        Collect(D->getBeginLoc());
    for (FileID FID : Files)
      SM.invalidateCache(FID);
  }

  /// Like cling's DeclUnloader::MaybeRemoveDeclFromModule(): make CodeGen
  /// forget the mangled names of the functions and variables of \p T.
  /// Otherwise CodeGen keeps their deferred definitions, and emits them when
  /// the name is used again, e.g. by the same template declared again after
  /// the unloaded file is loaded once more ("definition with same mangled
  /// name as another definition").
  static void forgetCodeGenDecls(const Transaction& T, CodeGenerator& CG) {
    llvm::SmallPtrSet<const Decl*, 32> Seen;
    auto Forget = [&](GlobalDecl GD) {
      std::string Name;
      utils::Analyze::maybeMangleDeclName(GD, Name);
      if (!Name.empty())
        CG.forgetDecl(Name);
    };
    std::function<void(Decl*)> Visit = [&](Decl* D) {
      if (!D || D->isInvalidDecl() || !Seen.insert(D).second)
        return;
      if (auto* FD = dyn_cast<FunctionDecl>(D)) {
        if (FD->isDependentContext())
          return;
        if (auto* CD = dyn_cast<CXXConstructorDecl>(FD)) {
          Forget(GlobalDecl(CD, Ctor_Complete));
          Forget(GlobalDecl(CD, Ctor_Base));
        } else if (auto* DD = dyn_cast<CXXDestructorDecl>(FD)) {
          Forget(GlobalDecl(DD, Dtor_Deleting));
          Forget(GlobalDecl(DD, Dtor_Complete));
          Forget(GlobalDecl(DD, Dtor_Base));
        } else {
          Forget(GlobalDecl(FD));
        }
      } else if (auto* VD = dyn_cast<VarDecl>(D)) {
        if (VD->hasGlobalStorage() && !VD->isStaticLocal() &&
            !VD->getDeclContext()->isDependentContext() &&
            !VD->getType()->isDependentType())
          Forget(GlobalDecl(VD));
      } else if (isa<NamespaceDecl>(D) || isa<LinkageSpecDecl>(D) ||
                 isa<CXXRecordDecl>(D)) {
        auto* DC = cast<DeclContext>(D);
        if (!DC->isDependentContext())
          for (Decl* Member : DC->decls())
            Visit(Member);
      }
    };
    for (auto I = T.decls_begin(), E = T.decls_end(); I != E; ++I)
      for (Decl* D : I->m_DGR)
        Visit(D);
  }

  void Interpreter::unload(Transaction& T) {
    // clang-repl can only undo the most recent inputs, in order. That is all
    // ROOT needs: unloading a file (`.U`, TSystem::Unload(), ACLiC) reverts
    // every transaction since the file was loaded, newest first.
    // FIXME: cling could unload any transaction (TransactionUnloader).
    if (getLastTransaction() != &T) {
      cling::errs() << "cling-compat: can only unload the last input\n";
      return;
    }

    for (unsigned i = 0; i < kNumTransactions; ++i)
      if (m_CachedTrns[i] == &T)
        m_CachedTrns[i] = nullptr;
    if (executingTransaction() == &T)
      executingTransaction() = nullptr;

    // Transactions without an input of their own (PushTransactionRAII, e.g.
    // template instantiations made by a lookup) have nothing that clang-repl
    // could undo: their declarations stay, like the instantiations made by
    // an input that is undone. Just drop them. Unloading must always pop the
    // transaction; callers revert until a given transaction is the last one.
    const bool HasInput = T.getPTU();
    if (HasInput)
      if (InterpreterCallbacks* callbacks = getCallbacks())
        callbacks->TransactionUnloaded(T);

    if (getOptions().ErrorOut) {
      // Tag the transaction as "won't need to be committed" (ROOT-10798).
      T.setState(Transaction::kRolledBack);
      return;
    }

    unloadFromPreprocessor(T, getSema());
    if (HasInput)
      if (CodeGenerator* CG = m_ClangInterp->getCodeGen())
        forgetCodeGenDecls(T, *CG);
    if (HasInput) {
      // Destroy the input's globals while their code is still there.
      runAtExitFuncs(T);
      if (InterpreterCallbacks* callbacks = getCallbacks())
        callbacks->TransactionRollback(T);
      if (llvm::Error Err = m_ClangInterp->Undo(1)) {
        reportError(std::move(Err), "failed to unload");
        T.setState(Transaction::kRolledBackWithErrors);
      } else {
        T.setState(Transaction::kRolledBack);
      }
    } else {
      T.setState(Transaction::kRolledBack);
    }
    m_Transactions.pop_back();
    if (!m_Transactions.empty())
      m_Transactions.back()->setNext(nullptr);
  }

  void Interpreter::unload(unsigned numberOfTransactions) {
    // `.undo N` reverts the last N inputs, and the transactions without an
    // input of their own (e.g. made by lookups) in between: cling had
    // nested those into the current transaction.
    for (unsigned i = 0; i < numberOfTransactions;) {
      Transaction* T = const_cast<Transaction*>(getLastTransaction());
      if (!T || T == getFirstTransaction()) {
        cling::errs() << "cling: Can't unload first transaction!  Unloaded "
                      << i << " of " << numberOfTransactions << "\n";
        return;
      }
      const bool IsInput = T->getPTU();
      unload(*T);
      if (IsInput)
        ++i;
    }
  }

  void Interpreter::runAndRemoveStaticDestructors() {
    // TROOT::EndOfProcessCleanups() calls this (via TCling::ResetGlobals())
    // before it deletes ROOT's lists of canvases, files, etc.: JIT-ed globals
    // that own such objects must be destroyed first. The atexit handler of
    // deinitializeJIT() then finds nothing left to run.
    runAtExitFuncs();
    deinitializeJIT(*this);
  }

  void Interpreter::runAndRemoveStaticDestructors(unsigned) {}

  void Interpreter::enableDynamicLookup(bool value /*=true*/) {
    if (!m_DynamicLookupDeclared && value && !isInSyntaxOnlyMode()) {
      // No dynlookup for the dynlookup header!
      m_DynamicLookupEnabled = false;
      if (declare("#include \"cling/Interpreter/DynamicLookupRuntimeUniverse.h\"")
          != kSuccess)
        return;
      m_DynamicLookupDeclared = true;
    }

    // Enable it *after* parsing the headers.
    m_DynamicLookupEnabled = value;
  }

  void
  Interpreter::addGenerator(std::unique_ptr<llvm::orc::DefinitionGenerator> G) {
    llvm::orc::LLJIT* JIT = getExecutionEngine();
    if (!JIT)
      return;
    // ORC asks a JITDylib's generators before it moves on to the next
    // JITDylib in the link order. Put client generators (e.g. ROOT's
    // library autoloading, which can scan all libraries) into their own
    // JITDylib at the end of the main JITDylib's link order, so that they only
    // see symbols that neither the JIT, its platform nor the process define.
    llvm::orc::JITDylib* GenJD =
        JIT->getExecutionSession().getJITDylibByName("<cling-compat generators>");
    if (!GenJD) {
      GenJD = &JIT->getExecutionSession().createBareJITDylib(
          "<cling-compat generators>");
      JIT->getMainJITDylib().addToLinkOrder(*GenJD);
    }
    GenJD->addGenerator(std::move(G));
  }

  Value Interpreter::Evaluate(const char* expr, DeclContext* /*DC*/,
                              bool ValuePrinterReq) {
    Sema& TheSema = getCI()->getSema();
    Sema::ContextRAII pushDC(TheSema,
                             TheSema.getASTContext().getTranslationUnitDecl());

    Value Result;
    if (InterpreterCallbacks* CB = getCallbacks())
      CB->SetIsRuntime(true);
    if (ValuePrinterReq)
      echo(expr, &Result);
    else
      evaluate(expr, Result);
    if (InterpreterCallbacks* CB = getCallbacks())
      CB->SetIsRuntime(false);

    return Result;
  }

  void Interpreter::setCallbacks(std::unique_ptr<InterpreterCallbacks> C) {
    // We need it to enable LookupObject callback.
    if (!m_Callbacks) {
      m_Callbacks.reset(new MultiplexInterpreterCallbacks(this));
      if (m_DyLibManager)
        m_DyLibManager->setCallbacks(m_Callbacks.get());
    }

    static_cast<MultiplexInterpreterCallbacks*>(m_Callbacks.get())
      ->addCallback(std::move(C));
  }

  llvm::orc::LLJIT* Interpreter::getExecutionEngine() {
    if (!m_Repl || isInSyntaxOnlyMode())
      return nullptr;
    return m_Repl->getExecutionEngine();
  }

  const DynamicLibraryManager* Interpreter::getDynamicLibraryManager() const {
    return m_DyLibManager.get();
  }

  DynamicLibraryManager* Interpreter::getDynamicLibraryManager() {
    return m_DyLibManager.get();
  }

  //--------------------------------------------------------------------------
  // Symbols and functions
  //--------------------------------------------------------------------------

  void* Interpreter::compileFunction(llvm::StringRef name, llvm::StringRef code,
                                     bool ifUnique, bool withAccessControl) {
    if (isInSyntaxOnlyMode())
      return nullptr;

    if (ifUnique) {
      if (void* Addr = (void*)getAddressOfGlobal(name))
        return Addr;
    }

    // See cling's DeclareCFunction for why this diagnostic is silenced.
    DiagnosticsEngine& Diag = getDiagnostics();
    Diag.setSeverity(clang::diag::ext_nested_name_member_ref_lookup_ambiguous,
                     clang::diag::Severity::Ignored, SourceLocation());

    LangOptions& LO = const_cast<LangOptions&>(getCI()->getLangOpts());
    bool savedAccessControl = LO.AccessControl;
    LO.AccessControl = withAccessControl;
    CompilationResult CR = declare(code.str());
    LO.AccessControl = savedAccessControl;

    Diag.setSeverity(clang::diag::ext_nested_name_member_ref_lookup_ambiguous,
                     clang::diag::Severity::Warning, SourceLocation());

    if (CR != kSuccess)
      return nullptr;

    return getAddressOfGlobal(name);
  }

  void* Interpreter::compileDtorCallFor(const clang::RecordDecl* RD) {
    void* &addr = m_DtorWrappers[RD];
    if (addr)
      return addr;

    if (const CXXRecordDecl *CXX = dyn_cast<CXXRecordDecl>(RD)) {
      // Don't generate a stub for a destructor that does nothing
      if (CXX->hasIrrelevantDestructor())
        return nullptr;
    }

    std::string funcname;
    llvm::raw_string_ostream fn(funcname);
    fn << "__cling_Destruct_" << RD;
    fn.flush();

    std::string code;
    llvm::raw_string_ostream c(code);
    c << "extern \"C\" void " << funcname << "(void* obj){(("
      << utils::TypeName::GetFullyQualifiedName(
             RD->getASTContext().getCanonicalTagType(RD), RD->getASTContext())
      << "*)obj)->~" << RD->getNameAsString() << "();}";
    c.flush();

    // ifUniq = false: we know it's unique, no need to check.
    addr = compileFunction(funcname, code, false /*ifUniq*/,
                           false /*withAccessControl*/);
    return addr;
  }

  void* Interpreter::getAddressOfGlobal(const GlobalDecl& GD,
                                        bool* fromJIT /*=0*/) const {
    std::string mangledName;
    utils::Analyze::maybeMangleDeclName(GD, mangledName);
    return getAddressOfGlobal(mangledName, fromJIT);
  }

  void* Interpreter::getAddressOfGlobal(llvm::StringRef SymName,
                                        bool* fromJIT /*=0*/) const {
    if (fromJIT)
      *fromJIT = false;
    if (isInSyntaxOnlyMode() || !m_Repl)
      return nullptr;

    // Symbols of the process and of loaded libraries first, then the JIT.
    // FIXME: cling asked the JIT first, but only for symbols it had emitted.
    if (void* Addr = llvm::sys::DynamicLibrary::SearchForAddressOfSymbol(
            SymName.str()))
      return Addr;

    auto AddrOrErr = m_ClangInterp->getSymbolAddressFromLinkerName(SymName);
    if (!AddrOrErr) {
      llvm::consumeError(AddrOrErr.takeError()); // okay to be missing
      // The definition might only have been handed to CodeGen so far (e.g.
      // cling's value printer builds its wrappers directly in the AST).
      if (!const_cast<Interpreter*>(this)->flushPendingCodeGen())
        return nullptr;
      AddrOrErr = m_ClangInterp->getSymbolAddressFromLinkerName(SymName);
      if (!AddrOrErr) {
        llvm::consumeError(AddrOrErr.takeError());
        return nullptr;
      }
    }
    if (fromJIT)
      *fromJIT = true;
    return AddrOrErr->toPtr<void*>();
  }

  bool Interpreter::flushPendingCodeGen() {
    // clang-repl hands the llvm::Module to the JIT at the end of Parse(); code
    // that CodeGen generated for declarations passed to it outside of Parse()
    // waits for the next input. Parse (and run) an empty one.
    // FIXME: clang-repl could offer a way to emit and execute pending CodeGen.
    if (isInSyntaxOnlyMode() || m_InIncrementalParse || !m_ClangInterp)
      return false;
    if (m_OpenTransactions.empty())
      getSema().PerformPendingInstantiations();
    llvm::Expected<clang::PartialTranslationUnit&> PTU = [&] {
      llvm::SaveAndRestore<bool> InParse(m_InIncrementalParse, true);
      return m_ClangInterp->Parse("");
    }();
    if (!PTU) {
      llvm::consumeError(PTU.takeError());
      return false;
    }
    if (!PTU->TheModule ||
        (PTU->TheModule->empty() && PTU->TheModule->global_empty()))
      return false;
    if (llvm::Error Err = m_ClangInterp->Execute(*PTU)) {
      reportError(std::move(Err), "failed to execute");
      return false;
    }
    return true;
  }

  void Interpreter::AddAtExitFunc(void (*Func) (void*), void* Arg) {
    std::lock_guard<std::mutex> Lock(m_AtExitFuncsMutex);
    m_AtExitFuncs.push_back({Func, Arg, executingTransaction()});
  }

  void Interpreter::runAtExitFuncs(const Transaction& T) {
    // Newest first; a destructor may register more.
    while (true) {
      AtExitFunc F;
      {
        std::lock_guard<std::mutex> Lock(m_AtExitFuncsMutex);
        auto It = std::find_if(m_AtExitFuncs.rbegin(), m_AtExitFuncs.rend(),
                               [&](const AtExitFunc& A) { return A.m_T == &T; });
        if (It == m_AtExitFuncs.rend())
          return;
        F = *It;
        m_AtExitFuncs.erase(std::next(It).base());
      }
      (*F.m_Func)(F.m_Arg);
    }
  }

  void Interpreter::runAtExitFuncs() {
    // Run in reverse order of registration, like atexit().
    std::vector<AtExitFunc> Funcs;
    {
      std::lock_guard<std::mutex> Lock(m_AtExitFuncsMutex);
      Funcs = std::move(m_AtExitFuncs);
      m_AtExitFuncs.clear();
    }
    for (auto I = Funcs.rbegin(), E = Funcs.rend(); I != E; ++I)
      (*I->m_Func)(I->m_Arg);
    // The static destructors of JIT-ed code (registered through __cxa_atexit
    // with the JIT's platform support) run in deinitializeJITsAtExit().
  }

  //--------------------------------------------------------------------------
  // Value printing helpers
  //--------------------------------------------------------------------------

  namespace valuePrinterInternal {
    void declarePrintValue(Interpreter &Interp);
  }

  std::string Interpreter::toString(const char* type, void* obj) {
    LockCompilationDuringUserCodeExecutionRAII LCDUCER(*this);
    cling::valuePrinterInternal::declarePrintValue(*this);
    std::string buf, ret;
    llvm::raw_string_ostream ss(buf);
    ss << "*((std::string*)" << &ret << ") = cling::printValue((" << type << "*)"
       << obj << ");";
    CompilationResult result = process(ss.str().c_str());
    if (result != cling::Interpreter::kSuccess)
      llvm::errs() << "Error in Interpreter::toString: the input " << ss.str()
                   << " cannot be evaluated";

    return ret;
  }

  //--------------------------------------------------------------------------
  // Forward declarations (rootmap / rootcling)
  //--------------------------------------------------------------------------

  void Interpreter::GenerateAutoLoadingMap(llvm::StringRef inFile,
                                           llvm::StringRef outFile,
                                           bool enableMacros,
                                           bool enableLogs) {
    const char *const dummy="cling_fwd_declarator";
    // Create an interpreter without any runtime, producing the fwd decls.
    std::string llvmdir
      = getCI()->getHeaderSearchOpts().ResourceDir + "/../../../";
    cling::Interpreter fwdGen(1, &dummy, llvmdir.c_str(),
                              /*moduleExtensions*/ {},
                              /*ExtraLibHandle*/ nullptr, /*noRuntime=*/true);

    // Copy the same header search options to the new instance.
    Preprocessor& fwdGenPP = fwdGen.getCI()->getPreprocessor();
    HeaderSearchOptions headerOpts = getCI()->getHeaderSearchOpts();
    clang::ApplyHeaderSearchOptions(fwdGenPP.getHeaderSearchInfo(), headerOpts,
                                    fwdGenPP.getLangOpts(),
                                    fwdGenPP.getTargetInfo().getTriple());

    std::string includeFile = std::string("#include \"") + inFile.str() + "\"";
    Transaction* T = nullptr;
    if (fwdGen.parse(includeFile, &T) != kSuccess || !T)
      return;

    std::error_code EC;
    llvm::raw_fd_ostream out(outFile.data(), EC,
                             llvm::sys::fs::OpenFlags::OF_None);
    llvm::raw_fd_ostream log((outFile + ".skipped").str().c_str(),
                             EC, llvm::sys::fs::OpenFlags::OF_None);
    log << "Generated for :" << inFile << "\n";
    forwardDeclare(*T, fwdGenPP, fwdGen.getCI()->getSema().getASTContext(),
                   out, enableMacros, enableLogs ? &log : nullptr);
  }

  void Interpreter::forwardDeclare(Transaction& T, Preprocessor& P,
                                   clang::ASTContext& Ctx,
                                   llvm::raw_ostream& out,
                                   bool enableMacros /*=false*/,
                                   llvm::raw_ostream* logs /*=0*/,
                                   IgnoreFilesFunc_t ignoreFiles) const {
    llvm::raw_null_ostream null;
    if (!logs)
      logs = &null;

    ForwardDeclPrinter visitor(out, *logs, P, Ctx, T, 0, false, ignoreFiles);
    visitor.printStats();

    T.setState(Transaction::kCommitted);
  }

} // end namespace cling
