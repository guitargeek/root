//===--- Interpreter.h - Incremental Compilation and Execution---*- C++ -*-===//
//
// Part of the LLVM Project, under the Apache License v2.0 with LLVM Exceptions.
// See https://llvm.org/LICENSE.txt for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//===----------------------------------------------------------------------===//
//
// This file defines the component which performs incremental code
// compilation and execution.
//
//===----------------------------------------------------------------------===//

#ifndef LLVM_CLANG_INTERPRETER_INTERPRETER_H
#define LLVM_CLANG_INTERPRETER_INTERPRETER_H

#include "clang/AST/GlobalDecl.h"
#include "clang/Interpreter/IncrementalExecutor.h"
#include "clang/Interpreter/PartialTranslationUnit.h"
#include "clang/Interpreter/Value.h"

#include "llvm/ADT/DenseMap.h"
#include "llvm/ExecutionEngine/JITSymbol.h"
#include "llvm/ExecutionEngine/Orc/ExecutorProcessControl.h"
#include "llvm/ExecutionEngine/Orc/Shared/ExecutorAddress.h"
#include "llvm/Support/Error.h"
#include <cstdint>
#include <functional>
#include <map>
#include <memory>
#include <mutex>
#include <thread>
#include <vector>

namespace llvm {
namespace orc {
class ThreadSafeContext;
} // namespace orc
} // namespace llvm

namespace clang {

namespace driver {
class Compilation;
} // namespace driver

class CodeGenerator;
class CompilerInstance;
class CXXRecordDecl;
class Decl;
class IncrementalParser;
class IncrementalCUDADeviceParser;
class Parser;

/// Create a pre-configured \c CompilerInstance for incremental processing.
class IncrementalCompilerBuilder {
  using DriverCompilationFn = llvm::Error(const driver::Compilation &);

public:
  IncrementalCompilerBuilder() {}

  void SetCompilerArgs(const std::vector<const char *> &Args) {
    UserArgs = Args;
  }

  void SetTargetTriple(std::string TT) { TargetTriple = TT; }

  // General C++
  llvm::Expected<std::unique_ptr<CompilerInstance>> CreateCpp();

  // Offload options
  void SetOffloadArch(llvm::StringRef Arch) { OffloadArch = Arch; };

  // CUDA specific
  void SetCudaSDK(llvm::StringRef path) { CudaSDKPath = path; };

  // Hand over the compilation.
  void SetDriverCompilationCallback(std::function<DriverCompilationFn> C) {
    CompilationCB = C;
  }

  llvm::Expected<std::unique_ptr<CompilerInstance>> CreateCudaHost();
  llvm::Expected<std::unique_ptr<CompilerInstance>> CreateCudaDevice();

private:
  llvm::Expected<std::unique_ptr<CompilerInstance>>
  create(std::string TT, std::vector<const char *> &ClangArgv);

  llvm::Expected<std::unique_ptr<CompilerInstance>> createCuda(bool device);

  std::vector<const char *> UserArgs;
  std::optional<std::string> TargetTriple;

  llvm::StringRef OffloadArch;
  llvm::StringRef CudaSDKPath;

  std::optional<std::function<DriverCompilationFn>> CompilationCB;
};

class IncrementalAction;
class InProcessPrintingASTConsumer;

/// Provides top-level interfaces for incremental compilation and execution.
class Interpreter {
  friend class Value;
  friend InProcessPrintingASTConsumer;

  std::unique_ptr<llvm::orc::ThreadSafeContext> TSCtx;
  /// Long-lived, incremental parsing action.
  std::unique_ptr<IncrementalAction> Act;
  std::unique_ptr<IncrementalParser> IncrParser;
  std::unique_ptr<IncrementalExecutor> IncrExecutor;

  // An optional parser for CUDA offloading
  std::unique_ptr<IncrementalCUDADeviceParser> DeviceParser;

  // An optional action for CUDA offloading
  std::unique_ptr<IncrementalAction> DeviceAct;

  /// List containing information about each incrementally parsed piece of code.
  std::list<PartialTranslationUnit> PTUs;

  unsigned InitPTUSize = 0;

  // This member holds the last result of the value printing. It's a class
  // member because we might want to access it after more inputs. If no value
  // printing happens, it's in an invalid state.
  Value LastValue;

  /// Inputs may run concurrently on several threads; each thread gets its own
  /// result. The generated code passes &LastValue, see getLastValueSlot().
  std::mutex LastValuesMutex;
  std::map<std::thread::id, Value> LastValues;

  /// See setCodeGenFilter().
  std::function<bool(const Decl *)> CodeGenFilter;

  /// See setCaptureValueOfExprStmts().
  bool CaptureValueOfExprStmts = false;

  /// Compiler instance performing the incremental compilation.
  std::unique_ptr<CompilerInstance> CI;

  /// An optional compiler instance for CUDA offloading
  std::unique_ptr<CompilerInstance> DeviceCI;

protected:
  // Derived classes can use an extended interface of the Interpreter.
  Interpreter(std::unique_ptr<CompilerInstance> Instance, llvm::Error &Err,
              std::unique_ptr<IncrementalExecutorBuilder> IEB = nullptr,
              std::unique_ptr<clang::ASTConsumer> Consumer = nullptr,
              std::function<std::unique_ptr<ASTConsumer>(CompilerInstance &)>
                  ConsumerFactory = nullptr);

  // Create the internal IncrementalExecutor, or re-create it after calling
  // ResetExecutor().
  llvm::Error CreateExecutor();

  // Delete the internal IncrementalExecutor. This causes a hard shutdown of the
  // JIT engine. In particular, it doesn't run cleanup or destructors.
  void ResetExecutor() { IncrExecutor.reset(); }

public:
  virtual ~Interpreter();
  static llvm::Expected<std::unique_ptr<Interpreter>>
  create(std::unique_ptr<CompilerInstance> CI,
         std::unique_ptr<IncrementalExecutorBuilder> IEB = nullptr);
  /// Creates an additional ASTConsumer once the CompilerInstance has its
  /// Preprocessor and ASTContext (i.e. when the incremental action creates
  /// its consumers).
  using ASTConsumerFactory =
      std::function<std::unique_ptr<ASTConsumer>(CompilerInstance &)>;

  /// Like create(CI, IEB), but also passes every top-level declaration (and
  /// the other ASTConsumer callbacks) to the consumer made by \p Factory,
  /// e.g. to track what each incremental input declared or to write a
  /// module file. Value printing keeps working.
  static llvm::Expected<std::unique_ptr<Interpreter>>
  create(std::unique_ptr<CompilerInstance> CI,
         std::unique_ptr<IncrementalExecutorBuilder> IEB,
         ASTConsumerFactory Factory);
  static llvm::Expected<std::unique_ptr<Interpreter>>
  createWithCUDA(std::unique_ptr<CompilerInstance> CI,
                 std::unique_ptr<CompilerInstance> DCI);

  const ASTContext &getASTContext() const;
  ASTContext &getASTContext();
  const CompilerInstance *getCompilerInstance() const;
  CompilerInstance *getCompilerInstance();

  /// The parser used for the incremental inputs. Clients that parse
  /// additional snippets on the side (e.g. to look up a type from its name)
  /// must use this parser, as a second Parser cannot share the Preprocessor.
  Parser &getParser();
  llvm::Expected<IncrementalExecutor &> getExecutionEngine();

  llvm::Expected<PartialTranslationUnit &> Parse(llvm::StringRef Code);
  llvm::Error Execute(PartialTranslationUnit &T);
  llvm::Error ParseAndExecute(llvm::StringRef Code, Value *V = nullptr);

  /// Moves the value of the last executed expression statement that was
  /// missing its semicolon into \p V (leaving it invalid if there is none).
  /// For clients that call Parse() and Execute() separately.
  void takeLastValue(Value &V);

  /// The value of the last input executed on this thread, where the value
  /// printing runtime stores the result of an input.
  Value &getLastValueSlot();

  /// Also capture the value of top-level expression statements that end in
  /// a semicolon (for takeLastValue(); such values are never printed). For
  /// clients that need the value of an input like `&f;`: dropping the
  /// semicolon would turn a declaration into a parse error.
  void setCaptureValueOfExprStmts(bool Capture) {
    CaptureValueOfExprStmts = Capture;
  }
  bool capturesValueOfExprStmts() const { return CaptureValueOfExprStmts; }

  /// Undo N previous incremental inputs.
  llvm::Error Undo(unsigned N = 1);

  /// The code generator of the inputs (null in syntax-only mode). Clients
  /// that undo inputs use it to make CodeGen forget their declarations.
  CodeGenerator *getCodeGen() const;

  /// Link a dynamic library
  llvm::Error LoadDynamicLibrary(const char *name);

  /// \returns the \c ExecutorAddr of a \c GlobalDecl. This interface uses
  /// the CodeGenModule's internal mangling cache to avoid recomputing the
  /// mangled name.
  llvm::Expected<llvm::orc::ExecutorAddr> getSymbolAddress(GlobalDecl GD) const;

  /// \returns the \c ExecutorAddr of a given name as written in the IR.
  llvm::Expected<llvm::orc::ExecutorAddr>
  getSymbolAddress(llvm::StringRef IRName) const;

  /// \returns the \c ExecutorAddr of a given name as written in the object
  /// file.
  llvm::Expected<llvm::orc::ExecutorAddr>
  getSymbolAddressFromLinkerName(llvm::StringRef LinkerName) const;

  const IncrementalExecutorBuilder &getIncrementalExecutorBuilder() const {
    return *IncrExecutorBuilder;
  }

  /// The inputs parsed so far, oldest first. Undo() removes them from the
  /// end.
  const std::list<PartialTranslationUnit> &getPTUs() const { return PTUs; }

  /// Decides which declarations are passed to CodeGen, including the members
  /// of namespaces and linkage specifications; e.g. to skip definitions that
  /// a loaded library provides. Without a filter, all declarations are.
  void setCodeGenFilter(std::function<bool(const Decl *)> Filter) {
    CodeGenFilter = std::move(Filter);
  }
  const std::function<bool(const Decl *)> &getCodeGenFilter() const {
    return CodeGenFilter;
  }

private:
  size_t getEffectivePTUSize() const;
  void markUserCodeStart();

  // A cache for the compiled destructors used to for de-allocation of managed
  // clang::Values.
  mutable llvm::DenseMap<CXXRecordDecl *, llvm::orc::ExecutorAddr> Dtors;

  std::array<Expr *, 4> ValuePrintingInfo = {0};

  std::unique_ptr<IncrementalExecutorBuilder> IncrExecutorBuilder;

  /// @}
  /// @name Value and pretty printing support
  /// @{

  std::string ValueDataToString(const Value &V) const;
  std::string ValueTypeToString(const Value &V) const;

  llvm::Expected<Expr *> convertExprToValue(Expr *E);

  // When we deallocate clang::Value we need to run the destructor of the type.
  // This function forces emission of the needed dtor.
  llvm::Expected<llvm::orc::ExecutorAddr>
  CompileDtorCall(CXXRecordDecl *CXXRD) const;
};
} // namespace clang

#endif // LLVM_CLANG_INTERPRETER_INTERPRETER_H
