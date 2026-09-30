//--------------------------------------------------------------------*- C++ -*-
// cling-compat: the subset of cling's Interpreter API that ROOT uses,
// implemented on top of clang::Interpreter (clang-repl).
//
// Original cling interface by Lukasz Janyst, Vassil Vassilev, Axel Naumann.
// This file is dual-licensed: you can choose to license it under the University
// of Illinois Open Source License or the GNU Lesser General Public License. See
// LICENSE.TXT for details.
//------------------------------------------------------------------------------

#ifndef CLING_INTERPRETER_H
#define CLING_INTERPRETER_H

/// Defined when the cling API is provided by cling-compat on clang-repl.
#define CLING_COMPAT_CLANG_REPL 1

#ifndef LLVM_PATH
#define LLVM_PATH nullptr
#endif

#include "cling/Interpreter/InvocationOptions.h"
#include "cling/Interpreter/RuntimeOptions.h"

#include "llvm/ADT/StringRef.h"

#include <cstdlib>
#include <deque>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

namespace llvm {
  class raw_ostream;
  template <typename T> class SmallVectorImpl;
  namespace orc {
    class DefinitionGenerator;
    class LLJIT;
  }
}

namespace clang {
  class ASTConsumer;
  class ASTContext;
  class CompilerInstance;
  class Decl;
  class DeclContext;
  class DiagnosticConsumer;
  class DiagnosticsEngine;
  class FunctionDecl;
  class GlobalDecl;
  class MacroInfo;
  class Module;
  class ModuleFileExtension;
  class NamedDecl;
  class Parser;
  class Preprocessor;
  class PresumedLoc;
  class QualType;
  class RecordDecl;
  class Sema;
  class SourceLocation;
  struct PartialTranslationUnit;
}

namespace CppInternal {
  class Interpreter;
}

namespace cling {
  namespace runtime {
    namespace internal {
      class LifetimeHandler;
    }
  }
  class ClangInternalState;
  class CompilationOptions;
  class DeclCollector;
  class DynamicLibraryManager;
  class InterpreterCallbacks;
  class LookupHelper;
  class Transaction;
  class Value;

  ///\brief Class that implements the interpreter-like behavior. It manages the
  /// incremental compilation, now delegated to clang::Interpreter.
  ///
  class Interpreter {
  public:
     // IgnoreFilesFunc_t takes a const reference to avoid having to
     // include the actual definition of PresumedLoc.
     using IgnoreFilesFunc_t = bool (*)(const clang::PresumedLoc&);

    ///\brief Pushes a new transaction, which will collect the decls that came
    /// within the scope of the RAII object. Calls commit transaction at
    /// destruction.
    class PushTransactionRAII {
    private:
      Transaction* m_Transaction;
      const Interpreter* m_Interpreter;
    public:
      PushTransactionRAII(const Interpreter* i);
      ~PushTransactionRAII();
      void pop() const;
    };

    ///\brief Describes the return result of the different routines that do the
    /// incremental compilation.
    ///
    enum CompilationResult {
      kSuccess,
      kFailure,
      kMoreInputExpected
    };

    ///\brief Describes the result of running a function.
    ///
    enum ExecutionResult {
      kExeSuccess,
      kExeNoCodeGen,
      kExeFirstError,
      kExeFunctionNotCompiled = kExeFirstError,
      kExeUnresolvedSymbols,
      kExeCompilationError,
      kExeUnkownFunction,
      kExeNoModule,
      kNumExeResults
    };

    ///\brief Flags that provide additional information about the context in
    /// which parsing occurs. These flags can be bitwise-OR'ed together.
    enum InputFlags {
      ///\brief Input comes from an external file
      kInputFromFile = 0x01,
      ///\brief If `kInputFromFile == 1`, whether `Interpreter::process()` is
      /// called once per line
      kIFFLineByLine = 0x02
    };

    ///\brief A RAII object that temporarily sets `Interpreter::m_InputFlags`
    /// and restores it on destruction.
    struct InputFlagsRAII {
      Interpreter &m_Interp;
      unsigned m_OldFlags;
      InputFlagsRAII(Interpreter &I, unsigned F)
        : m_Interp(I), m_OldFlags(I.m_InputFlags) { I.m_InputFlags = F; }
      ~InputFlagsRAII() { m_Interp.m_InputFlags = m_OldFlags; }
    };

    using ModuleFileExtensions =
        std::vector<std::shared_ptr<clang::ModuleFileExtension>>;

  private:
    ///\brief Interpreter invocation options.
    InvocationOptions m_Opts;

    ///\brief CppInterOp's wrapper around the clang::Interpreter. It owns the
    /// clang::Interpreter, which is shared with CppInterOp (and thus with
    /// PyROOT's cppjit) through Cpp::UseExternalInterpreter().
    std::unique_ptr<CppInternal::Interpreter> m_Repl;

    ///\brief The clang::Interpreter, owned by m_Repl. Kept opaque (like the
    /// other clang-repl types) so that `using namespace clang; using
    /// namespace cling;` in clients stays unambiguous.
    void* m_ClangInterpOpaque = nullptr;

    ///\brief Collects the decls of each input into the current Transaction.
    /// Owned by the clang::Interpreter's consumer chain.
    DeclCollector* m_DeclCollector = nullptr;

    ///\brief Cling's reflection information query.
    std::unique_ptr<LookupHelper> m_LookupHelper;

    ///\brief The dynamic library manager (a copy of cling's).
    std::unique_ptr<DynamicLibraryManager> m_DyLibManager;

    ///\brief Interpreter callbacks.
    std::unique_ptr<InterpreterCallbacks> m_Callbacks;

    ///\brief All transactions that were committed, in order.
    std::vector<std::unique_ptr<Transaction>> m_Transactions;

    ///\brief The currently open transactions (innermost last).
    mutable std::vector<Transaction*> m_OpenTransactions;

    ///\brief clang::Values of record type that back cling::Values handed out
    /// by process() and evaluate(). cling::Value cannot own a clang::Value, so
    /// these are kept alive until the interpreter dies.
    void* m_KeptAliveValues = nullptr; // std::deque<clang::Value>*

    ///\brief Cache of compiled destructors wrappers.
    std::unordered_map<const clang::RecordDecl*, void*> m_DtorWrappers;

    ///\brief Counter used when we need unique names.
    mutable unsigned long long m_UniqueCounter = 0;

    ///\brief Source of unique source locations (see getNextAvailableLoc()):
    /// the raw encoding of the start of a virtual buffer, and the next offset.
    mutable unsigned m_UniqueLocBufferStart = 0;
    mutable unsigned m_UniqueLocOffset = 0;

    bool m_PrintDebug = false;
    bool m_DynamicLookupEnabled = false;
    bool m_RawInputEnabled = false;
    bool m_SyntaxOnly = false;
    unsigned m_InputFlags{};
    cling::runtime::RuntimeOptions m_RuntimeOptions;
    int m_OptLevel = 2;

    ///\brief cling's filter of incremental-input diagnostics; owns the
    /// original DiagnosticConsumer and forwards to a replacement, if any.
    std::unique_ptr<clang::DiagnosticConsumer> m_DiagConsumer;

    ///\brief Whether clang-repl is currently parsing an input. clang-repl
    /// ends the "translation unit" after every input; the module writer only
    /// honors explicit end-of-TU calls (rootcling), made outside of parsing.
    mutable bool m_InIncrementalParse = false;

    ///\brief Emit and run what CodeGen generated outside of an input.
    ///\returns whether anything was emitted.
    bool flushPendingCodeGen();

    ///\brief With -fmodule-name, the consumers writing the module file.
    std::unique_ptr<clang::ASTConsumer>
    createModuleWriter(clang::CompilerInstance& CI);

    ///\brief Registered atexit functions.
    std::vector<std::pair<void (*)(void*), void*>> m_AtExitFuncs;

    enum {
      kStdStringTransaction = 0, // Transaction known to contain std::string
      kNumTransactions
    };
    mutable const Transaction* m_CachedTrns[kNumTransactions] = {};

    ///\brief Parse and (unless in syntax-only mode) execute \p input as one
    /// incremental input, collecting its decls into a new Transaction.
    ///
    ///\param [out] CV - if non-null, receives the value of the last
    ///             expression statement (when its semicolon is missing).
    ///\param [in] print - whether to print the value of that expression.
    ///
    CompilationResult ParseAndExecute(const std::string& input,
                                      const CompilationOptions& CO,
                                      Value* CV, bool print,
                                      Transaction** T = nullptr);

    ///\brief Convert the value captured by clang-repl to a cling::Value.
    void convertValue(void* From /*clang::Value**/, Value& To);

    ///\brief Begin/end a Transaction that collects decls.
    Transaction* beginTransaction(const CompilationOptions& CO) const;
    void endTransaction(Transaction* T, bool Commit) const;

    void Initialize(bool NoRuntime, bool SyntaxOnly);

  public:
    ///\brief Constructor for Interpreter.
    ///
    ///\param[in] argc - no. of args.
    ///\param[in] argv - arguments passed when driver is invoked.
    ///\param[in] llvmdir - the resource directory
    ///\param[in] extraLibHandle - resolve symbols also from this dylib
    ///\param[in] noRuntime - flag to control the presence of runtime universe
    ///
    Interpreter(int argc, const char* const* argv, const char* llvmdir = LLVM_PATH,
                const ModuleFileExtensions& moduleExtensions = {},
                void *extraLibHandle = nullptr, bool noRuntime = false);

    virtual ~Interpreter();

    ///\brief Whether the Interpreter is setup and ready to be used.
    bool isValid() const;

    const InvocationOptions& getOptions() const { return m_Opts; }
    InvocationOptions& getOptions() { return m_Opts; }

    const cling::runtime::RuntimeOptions& getRuntimeOptions() const { return m_RuntimeOptions; }
    cling::runtime::RuntimeOptions& getRuntimeOptions() { return m_RuntimeOptions; }

    LookupHelper& getLookupHelper() const { return *m_LookupHelper; }

    const clang::Parser& getParser() const;
    clang::Parser& getParser();


    ///\brief CppInterOp's handle of the interpreter (to share it).
    CppInternal::Interpreter* getCppInterOpInterpreter() const {
      return m_Repl.get();
    }

    ///\brief Returns the current or last Transactions source location.
    clang::SourceLocation getSourceLocation(bool skipWrapper = true) const;

    ///\brief Returns the next available valid free source location.
    clang::SourceLocation getNextAvailableLoc() const;

    ///\brief true if -fsyntax-only flag passed.
    bool isInSyntaxOnlyMode() const;

    static const char* getVersion();

    void createUniqueName(std::string& out);
    bool isUniqueName(llvm::StringRef name);

    void AddIncludePaths(llvm::StringRef PathsStr, const char* Delim = ":");
    void AddIncludePath(llvm::StringRef PathsStr);
    void GetIncludePaths(llvm::SmallVectorImpl<std::string>& incpaths,
                         bool withSystem, bool withFlags);

    ///\brief Call printValue( "(type*)" + obj ) and return string
    std::string toString(const char* type, void* obj);

    void DumpIncludePath(llvm::raw_ostream* S = nullptr) const;
    void DumpDynamicLibraryInfo(llvm::raw_ostream* S = nullptr) const;
    void dump(llvm::StringRef what, llvm::StringRef filter);
    void storeInterpreterState(const std::string& name) const;
    void compareInterpreterState(const std::string& name) const;
    void printIncludedFiles (llvm::raw_ostream& out) const;

    ///\brief Compiles the given input (declarations and/or statements) and
    /// runs it. Prints the value of a trailing expression without semicolon
    /// unless \p disableValuePrinting.
    CompilationResult process(const std::string& input, Value* V = nullptr,
                              Transaction** T = nullptr,
                              bool disableValuePrinting = false);

    ///\brief Parses input line, which doesn't contain statements. No code
    /// generation is done.
    CompilationResult parse(const std::string& input,
                            Transaction** T = nullptr) const;

    bool loadModule(const std::string& moduleName, bool complain = true);
    bool loadModule(clang::Module* M, bool complain = true);

    CompilationResult parseForModule(const std::string& input);

    CompilationResult codeComplete(const std::string& line, size_t& cursor,
                                   std::vector<std::string>& completions) const;

    CompilationResult declare(const std::string& input, Transaction** T = nullptr);
    CompilationResult evaluate(const std::string& input, Value& V);
    CompilationResult echo(const std::string& input, Value* V = nullptr);
    CompilationResult execute(const std::string& input);
    CompilationResult emitAllDecls(Transaction* T);

    std::string lookupFileOrLibrary(llvm::StringRef file);
    CompilationResult loadLibrary(const std::string& filename,
                                  bool lookup = true);
    CompilationResult loadHeader(const std::string& filename,
                                 Transaction** T = nullptr);
    CompilationResult loadFile(const std::string& filename,
                               bool allowSharedLib = true,
                               Transaction** T = nullptr);

    ///\brief Unloads (forgets) a transaction from AST and JITed symbols.
    /// Only possible for the most recent inputs (clang::Interpreter::Undo).
    void unload(Transaction& T);
    void unload(unsigned numberOfTransactions);

    void runAndRemoveStaticDestructors();
    void runAndRemoveStaticDestructors(unsigned numberOfTransactions);

    bool isPrintingDebug() const { return m_PrintDebug; }
    void enablePrintDebug(bool print = true) { m_PrintDebug = print; }

    void enableDynamicLookup(bool value = true);
    bool isDynamicLookupEnabled() const { return m_DynamicLookupEnabled; }

    bool isRawInputEnabled() const { return m_RawInputEnabled; }
    void enableRawInput(bool raw = true) { m_RawInputEnabled = raw; }

    unsigned getInputFlags() const { return m_InputFlags; }
    void setInputFlags(unsigned value) { m_InputFlags = value; }

    int getDefaultOptLevel() const { return m_OptLevel; }
    void setDefaultOptLevel(int optLevel) { m_OptLevel = optLevel; }

    clang::CompilerInstance* getCI() const;
    clang::CompilerInstance* getCIOrNull() const;
    clang::Sema& getSema() const;
    clang::DiagnosticsEngine& getDiagnostics() const;

    void replaceDiagnosticConsumer(clang::DiagnosticConsumer* Consumer, bool Own = false);
    bool hasReplacedDiagnosticConsumer() const;

    CompilationOptions makeDefaultCompilationOpts() const;

    void addGenerator(std::unique_ptr<llvm::orc::DefinitionGenerator> G);

    ///\brief Evaluates given expression within given declaration context.
    Value Evaluate(const char* expr, clang::DeclContext* DC,
                   bool ValuePrinterReq = false);

    void setCallbacks(std::unique_ptr<InterpreterCallbacks> C);
    const InterpreterCallbacks* getCallbacks() const {return m_Callbacks.get();}
    InterpreterCallbacks* getCallbacks() { return m_Callbacks.get(); }

    llvm::orc::LLJIT* getExecutionEngine();

    const DynamicLibraryManager* getDynamicLibraryManager() const;
    DynamicLibraryManager* getDynamicLibraryManager();

    const Transaction* getFirstTransaction() const;
    const Transaction* getLastTransaction() const;
    const Transaction* getLastWrapperTransaction() const;
    const Transaction* getCurrentTransaction() const;
    const Transaction* getLatestTransaction() const;

    const Transaction*& getStdStringTransaction() const {
      return m_CachedTrns[kStdStringTransaction];
    }

    void* compileFunction(llvm::StringRef name, llvm::StringRef code,
                          bool ifUniq = true, bool withAccessControl = true);

    void* compileDtorCallFor(const clang::RecordDecl* RD);

    void* getAddressOfGlobal(const clang::GlobalDecl& D,
                             bool* fromJIT = nullptr) const;
    void* getAddressOfGlobal(llvm::StringRef SymName, bool* fromJIT = nullptr) const;

    const clang::MacroInfo* getMacro(llvm::StringRef Name) const;
    std::string getMacroValue(llvm::StringRef Name,
                              const char* Strip = "\"") const;

    void AddAtExitFunc(void (*Func) (void*), void* Arg);
    void runAtExitFuncs();

    void GenerateAutoLoadingMap(llvm::StringRef inFile, llvm::StringRef outFile,
                                bool enableMacros = false, bool enableLogs = true);

    void forwardDeclare(Transaction& T, clang::Preprocessor& P,
                        clang::ASTContext& Ctx,
                        llvm::raw_ostream& out,
                        bool enableMacros = false,
                        llvm::raw_ostream* logs = nullptr,
                        IgnoreFilesFunc_t ignoreFiles =
                          [](const clang::PresumedLoc&) { return false;}) const;

    friend class runtime::internal::LifetimeHandler;
    friend class DeclCollector;
  };
} // namespace cling

#endif // CLING_INTERPRETER_H
