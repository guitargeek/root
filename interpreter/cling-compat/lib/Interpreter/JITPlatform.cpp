//------------------------------------------------------------------------------
// cling-compat: LLJIT's generic IR platform, fixed for nested initialization.
//
// A copy of GenericLLVMIRPlatformSupport and setUpGenericLLVMIRPlatform() from
// llvm/lib/ExecutionEngine/Orc/LLJIT.cpp (LLVM 22.1.8; Apache-2.0 WITH
// LLVM-exception), which ROOT may take from the system. The only change
// (marked "cling-compat") makes initialize() re-entrant; it is proposed
// upstream as clang-repl/pending/06-lljit-nested-initialize.patch. Remove
// this file once the LLVM that ROOT uses has the fix.
//------------------------------------------------------------------------------

#include "JITPlatform.h"

#include "llvm/ADT/ScopeExit.h"
#include "llvm/Analysis/TargetLibraryInfo.h"
#include "llvm/ExecutionEngine/Orc/EHFrameRegistrationPlugin.h"
#include "llvm/ExecutionEngine/Orc/ExecutionUtils.h"
#include "llvm/ExecutionEngine/Orc/LLJIT.h"
#include "llvm/ExecutionEngine/Orc/ObjectLinkingLayer.h"
#include "llvm/ExecutionEngine/Orc/UnwindInfoRegistrationPlugin.h"
#include "llvm/IR/GlobalVariable.h"
#include "llvm/IR/IRBuilder.h"
#include "llvm/IR/Module.h"
#include "llvm/Support/Debug.h"
#include "llvm/Support/FormatVariadic.h"
#include "llvm/Transforms/Utils/ModuleUtils.h"

#include <map>
#include <thread>

#define DEBUG_TYPE "orc"

using namespace llvm;
using namespace llvm::orc;

namespace {

/// Adds helper function decls and wrapper functions that call the helper with
/// some additional prefix arguments.
///
/// E.g. For wrapper "foo" with type i8(i8, i64), helper "bar", and prefix
/// args i32 4 and i16 12345, this function will add:
///
/// declare i8 @bar(i32, i16, i8, i64)
///
/// define i8 @foo(i8, i64) {
/// entry:
///   %2 = call i8 @bar(i32 4, i16 12345, i8 %0, i64 %1)
///   ret i8 %2
/// }
///
Function *addHelperAndWrapper(Module &M, StringRef WrapperName,
                              FunctionType *WrapperFnType,
                              GlobalValue::VisibilityTypes WrapperVisibility,
                              StringRef HelperName,
                              ArrayRef<Value *> HelperPrefixArgs) {
  std::vector<Type *> HelperArgTypes;
  for (auto *Arg : HelperPrefixArgs)
    HelperArgTypes.push_back(Arg->getType());
  llvm::append_range(HelperArgTypes, WrapperFnType->params());
  auto *HelperFnType =
      FunctionType::get(WrapperFnType->getReturnType(), HelperArgTypes, false);
  auto *HelperFn = Function::Create(HelperFnType, GlobalValue::ExternalLinkage,
                                    HelperName, M);

  auto *WrapperFn = Function::Create(
      WrapperFnType, GlobalValue::ExternalLinkage, WrapperName, M);
  WrapperFn->setVisibility(WrapperVisibility);

  auto *EntryBlock = BasicBlock::Create(M.getContext(), "entry", WrapperFn);
  IRBuilder<> IB(EntryBlock);

  std::vector<Value *> HelperArgs;
  llvm::append_range(HelperArgs, HelperPrefixArgs);
  for (auto &Arg : WrapperFn->args())
    HelperArgs.push_back(&Arg);
  auto *HelperResult = IB.CreateCall(HelperFn, HelperArgs);
  if (HelperFn->getReturnType()->isVoidTy())
    IB.CreateRetVoid();
  else
    IB.CreateRet(HelperResult);

  return WrapperFn;
}

class GenericLLVMIRPlatformSupport;

/// orc::Platform component of Generic LLVM IR Platform support.
/// Just forwards calls to the GenericLLVMIRPlatformSupport class below.
class GenericLLVMIRPlatform : public Platform {
public:
  GenericLLVMIRPlatform(GenericLLVMIRPlatformSupport &S) : S(S) {}
  Error setupJITDylib(JITDylib &JD) override;
  Error teardownJITDylib(JITDylib &JD) override;
  Error notifyAdding(ResourceTracker &RT,
                     const MaterializationUnit &MU) override;
  Error notifyRemoving(ResourceTracker &RT) override {
    // Noop -- Nothing to do (yet).
    return Error::success();
  }

private:
  GenericLLVMIRPlatformSupport &S;
};

/// This transform parses llvm.global_ctors to produce a single initialization
/// function for the module, records the function, then deletes
/// llvm.global_ctors.
class GlobalCtorDtorScraper {
public:
  GlobalCtorDtorScraper(GenericLLVMIRPlatformSupport &PS,
                        StringRef InitFunctionPrefix,
                        StringRef DeInitFunctionPrefix)
      : PS(PS), InitFunctionPrefix(InitFunctionPrefix),
        DeInitFunctionPrefix(DeInitFunctionPrefix) {}
  Expected<ThreadSafeModule> operator()(ThreadSafeModule TSM,
                                        MaterializationResponsibility &R);

private:
  GenericLLVMIRPlatformSupport &PS;
  StringRef InitFunctionPrefix;
  StringRef DeInitFunctionPrefix;
};

/// Generic IR Platform Support
///
/// Scrapes llvm.global_ctors and llvm.global_dtors and replaces them with
/// specially named 'init' and 'deinit'. Injects definitions / interposes for
/// some runtime API, including __cxa_atexit, dlopen, and dlclose.
class GenericLLVMIRPlatformSupport : public LLJIT::PlatformSupport {
public:
  GenericLLVMIRPlatformSupport(LLJIT &J, JITDylib &PlatformJD,
                               cling::JITPlatformHooks Hooks)
      : J(J), InitFunctionPrefix(J.mangle("__orc_init_func.")),
        DeInitFunctionPrefix(J.mangle("__orc_deinit_func.")),
        Hooks(std::move(Hooks)) {

    getExecutionSession().setPlatform(
        std::make_unique<GenericLLVMIRPlatform>(*this));

    setInitTransform(J, GlobalCtorDtorScraper(*this, InitFunctionPrefix,
                                              DeInitFunctionPrefix));

    SymbolMap StdInterposes;

    StdInterposes[J.mangleAndIntern("__lljit.platform_support_instance")] = {
        ExecutorAddr::fromPtr(this), JITSymbolFlags::Exported};
    StdInterposes[J.mangleAndIntern("__lljit.cxa_atexit_helper")] = {
        ExecutorAddr::fromPtr(registerCxaAtExitHelper), JITSymbolFlags()};

    cantFail(PlatformJD.define(absoluteSymbols(std::move(StdInterposes))));
    cantFail(setupJITDylib(PlatformJD));
    cantFail(J.addIRModule(PlatformJD, createPlatformRuntimeModule()));
  }

  ExecutionSession &getExecutionSession() { return J.getExecutionSession(); }

  /// Adds a module that defines the __dso_handle global.
  Error setupJITDylib(JITDylib &JD) {

    // Add per-jitdylib standard interposes.
    SymbolMap PerJDInterposes;
    PerJDInterposes[J.mangleAndIntern("__lljit.run_atexits_helper")] = {
        ExecutorAddr::fromPtr(runAtExitsHelper), JITSymbolFlags()};
    PerJDInterposes[J.mangleAndIntern("__lljit.atexit_helper")] = {
        ExecutorAddr::fromPtr(registerAtExitHelper), JITSymbolFlags()};
    cantFail(JD.define(absoluteSymbols(std::move(PerJDInterposes))));

    auto Ctx = std::make_unique<LLVMContext>();
    auto M = std::make_unique<Module>("__standard_lib", *Ctx);
    M->setDataLayout(J.getDataLayout());

    auto *Int64Ty = Type::getInt64Ty(*Ctx);
    auto *DSOHandle = new GlobalVariable(
        *M, Int64Ty, true, GlobalValue::ExternalLinkage,
        ConstantInt::get(Int64Ty, reinterpret_cast<uintptr_t>(&JD)),
        "__dso_handle");
    DSOHandle->setVisibility(GlobalValue::DefaultVisibility);
    DSOHandle->setInitializer(
        ConstantInt::get(Int64Ty, ExecutorAddr::fromPtr(&JD).getValue()));

    auto *GenericIRPlatformSupportTy =
        StructType::create(*Ctx, "lljit.GenericLLJITIRPlatformSupport");

    auto *PlatformInstanceDecl = new GlobalVariable(
        *M, GenericIRPlatformSupportTy, true, GlobalValue::ExternalLinkage,
        nullptr, "__lljit.platform_support_instance");

    auto *VoidTy = Type::getVoidTy(*Ctx);
    addHelperAndWrapper(
        *M, "__lljit_run_atexits", FunctionType::get(VoidTy, {}, false),
        GlobalValue::HiddenVisibility, "__lljit.run_atexits_helper",
        {PlatformInstanceDecl, DSOHandle});

    auto *IntTy = Type::getIntNTy(*Ctx, sizeof(int) * CHAR_BIT);
    auto *AtExitCallbackPtrTy = PointerType::getUnqual(*Ctx);
    auto *AtExit = addHelperAndWrapper(
        *M, "atexit", FunctionType::get(IntTy, {AtExitCallbackPtrTy}, false),
        GlobalValue::HiddenVisibility, "__lljit.atexit_helper",
        {PlatformInstanceDecl, DSOHandle});
    Attribute::AttrKind AtExitExtAttr =
        TargetLibraryInfo::getExtAttrForI32Return(J.getTargetTriple());
    if (AtExitExtAttr != Attribute::None)
      AtExit->addRetAttr(AtExitExtAttr);

    return J.addIRModule(JD, ThreadSafeModule(std::move(M), std::move(Ctx)));
  }

  Error notifyAdding(ResourceTracker &RT, const MaterializationUnit &MU) {
    auto &JD = RT.getJITDylib();
    if (auto &InitSym = MU.getInitializerSymbol())
      InitSymbols[&JD].add(InitSym, SymbolLookupFlags::WeaklyReferencedSymbol);
    else {
      // If there's no identified init symbol attached, but there is a symbol
      // with the GenericIRPlatform::InitFunctionPrefix, then treat that as
      // an init function. Add the symbol to both the InitSymbols map (which
      // will trigger a lookup to materialize the module) and the InitFunctions
      // map (which holds the names of the symbols to execute).
      for (auto &KV : MU.getSymbols())
        if ((*KV.first).starts_with(InitFunctionPrefix)) {
          InitSymbols[&JD].add(KV.first,
                               SymbolLookupFlags::WeaklyReferencedSymbol);
          InitFunctions[&JD].add(KV.first);
        } else if ((*KV.first).starts_with(DeInitFunctionPrefix)) {
          DeInitFunctions[&JD].add(KV.first);
        }
    }
    return Error::success();
  }

  Error initialize(JITDylib &JD) override {
    LLVM_DEBUG({
      dbgs() << "GenericLLVMIRPlatformSupport getting initializers to run\n";
    });
    // cling-compat: materializing the inputs may run code that calls
    // initialize() again, e.g. ROOT's autoloading generator, which loads a
    // library whose initializers declare and run more code (ROOT-10426).
    // Each call gets a frame that records the init functions registered
    // while it materializes, so that a nested call runs only those: the ones
    // of the outer calls are still being materialized and would never become
    // ready (deadlock).
    // The frames are per thread: initializers run without ROOT's lock (see
    // below), and another thread may initialize its inputs meanwhile.
    getExecutionSession().runSessionLocked(
        [&]() { InitFrames[std::this_thread::get_id()].emplace_back(); });
    scope_exit PopFrame([&]() {
      getExecutionSession().runSessionLocked([&]() {
        auto Itr = InitFrames.find(std::this_thread::get_id());
        Itr->second.pop_back();
        if (Itr->second.empty())
          InitFrames.erase(Itr);
      });
    });
    if (auto Initializers = getInitializers(JD)) {
      LLVM_DEBUG(
          { dbgs() << "GenericLLVMIRPlatformSupport running initializers\n"; });
      // cling-compat: the initializers are user code (clang-repl runs the
      // statements of an input as initializers); compiling them is not.
      void *UserCodeState =
          Hooks.EnterUserCode ? Hooks.EnterUserCode() : nullptr;
      scope_exit ReturnFromUserCode([&]() {
        if (Hooks.ReturnFromUserCode)
          Hooks.ReturnFromUserCode(UserCodeState);
      });
      for (auto InitFnAddr : *Initializers) {
        LLVM_DEBUG({
          dbgs() << "  Running init " << formatv("{0:x16}", InitFnAddr)
                 << "...\n";
        });
        auto *InitFn = InitFnAddr.toPtr<void (*)()>();
        InitFn();
      }
    } else
      return Initializers.takeError();
    return Error::success();
  }

  Error deinitialize(JITDylib &JD) override {
    LLVM_DEBUG({
      dbgs() << "GenericLLVMIRPlatformSupport getting deinitializers to run\n";
    });
    if (auto Deinitializers = getDeinitializers(JD)) {
      LLVM_DEBUG({
        dbgs() << "GenericLLVMIRPlatformSupport running deinitializers\n";
      });
      for (auto DeinitFnAddr : *Deinitializers) {
        LLVM_DEBUG({
          dbgs() << "  Running deinit " << formatv("{0:x16}", DeinitFnAddr)
                 << "...\n";
        });
        auto *DeinitFn = DeinitFnAddr.toPtr<void (*)()>();
        DeinitFn();
      }
    } else
      return Deinitializers.takeError();

    return Error::success();
  }

  void registerInitFunc(JITDylib &JD, SymbolStringPtr InitName) {
    getExecutionSession().runSessionLocked([&]() {
      InitFunctions[&JD].add(InitName);
      auto Itr = InitFrames.find(std::this_thread::get_id());
      if (Itr != InitFrames.end())
        Itr->second.back()[&JD].add(InitName);
    });
  }

  void registerDeInitFunc(JITDylib &JD, SymbolStringPtr DeInitName) {
    getExecutionSession().runSessionLocked(
        [&]() { DeInitFunctions[&JD].add(DeInitName); });
  }

private:
  Expected<std::vector<ExecutorAddr>> getInitializers(JITDylib &JD) {
    if (auto Err = issueInitLookups(JD))
      return std::move(Err);

    DenseMap<JITDylib *, SymbolLookupSet> LookupSymbols;
    std::vector<JITDylibSP> DFSLinkOrder;

    if (auto Err = getExecutionSession().runSessionLocked([&]() -> Error {
          if (auto DFSLinkOrderOrErr = JD.getDFSLinkOrder())
            DFSLinkOrder = std::move(*DFSLinkOrderOrErr);
          else
            return DFSLinkOrderOrErr.takeError();

          // cling-compat: a nested call (see initialize()) only runs the
          // init functions registered while it materialized its inputs.
          auto &Frames = InitFrames[std::this_thread::get_id()];
          if (Frames.size() > 1) {
            auto &Frame = Frames.back();
            for (auto &NextJD : DFSLinkOrder) {
              auto FItr = Frame.find(NextJD.get());
              if (FItr == Frame.end())
                continue;
              auto &Remaining = InitFunctions[NextJD.get()];
              for (auto &KV : FItr->second) {
                SymbolStringPtr Name = KV.first;
                Remaining.remove_if(
                    [&](const SymbolStringPtr &N, SymbolLookupFlags) {
                      return N == Name;
                    });
              }
              LookupSymbols[NextJD.get()] = std::move(FItr->second);
              Frame.erase(FItr);
            }
            return Error::success();
          }

          for (auto &NextJD : DFSLinkOrder) {
            auto IFItr = InitFunctions.find(NextJD.get());
            if (IFItr != InitFunctions.end()) {
              LookupSymbols[NextJD.get()] = std::move(IFItr->second);
              InitFunctions.erase(IFItr);
            }
          }
          return Error::success();
        }))
      return std::move(Err);

    LLVM_DEBUG({
      dbgs() << "JITDylib init order is [ ";
      for (auto &JD : llvm::reverse(DFSLinkOrder))
        dbgs() << "\"" << JD->getName() << "\" ";
      dbgs() << "]\n";
      dbgs() << "Looking up init functions:\n";
      for (auto &KV : LookupSymbols)
        dbgs() << "  \"" << KV.first->getName() << "\": " << KV.second << "\n";
    });

    auto &ES = getExecutionSession();
    auto LookupResult = Platform::lookupInitSymbols(ES, LookupSymbols);

    if (!LookupResult)
      return LookupResult.takeError();

    std::vector<ExecutorAddr> Initializers;
    while (!DFSLinkOrder.empty()) {
      auto &NextJD = *DFSLinkOrder.back();
      DFSLinkOrder.pop_back();
      auto InitsItr = LookupResult->find(&NextJD);
      if (InitsItr == LookupResult->end())
        continue;
      for (auto &KV : InitsItr->second)
        Initializers.push_back(KV.second.getAddress());
    }

    return Initializers;
  }

  Expected<std::vector<ExecutorAddr>> getDeinitializers(JITDylib &JD) {
    auto &ES = getExecutionSession();

    auto LLJITRunAtExits = J.mangleAndIntern("__lljit_run_atexits");

    DenseMap<JITDylib *, SymbolLookupSet> LookupSymbols;
    std::vector<JITDylibSP> DFSLinkOrder;

    if (auto Err = ES.runSessionLocked([&]() -> Error {
          if (auto DFSLinkOrderOrErr = JD.getDFSLinkOrder())
            DFSLinkOrder = std::move(*DFSLinkOrderOrErr);
          else
            return DFSLinkOrderOrErr.takeError();

          for (auto &NextJD : DFSLinkOrder) {
            auto &JDLookupSymbols = LookupSymbols[NextJD.get()];
            auto DIFItr = DeInitFunctions.find(NextJD.get());
            if (DIFItr != DeInitFunctions.end()) {
              LookupSymbols[NextJD.get()] = std::move(DIFItr->second);
              DeInitFunctions.erase(DIFItr);
            }
            JDLookupSymbols.add(LLJITRunAtExits,
                                SymbolLookupFlags::WeaklyReferencedSymbol);
          }
          return Error::success();
        }))
      return std::move(Err);

    LLVM_DEBUG({
      dbgs() << "JITDylib deinit order is [ ";
      for (auto &JD : DFSLinkOrder)
        dbgs() << "\"" << JD->getName() << "\" ";
      dbgs() << "]\n";
      dbgs() << "Looking up deinit functions:\n";
      for (auto &KV : LookupSymbols)
        dbgs() << "  \"" << KV.first->getName() << "\": " << KV.second << "\n";
    });

    auto LookupResult = Platform::lookupInitSymbols(ES, LookupSymbols);

    if (!LookupResult)
      return LookupResult.takeError();

    std::vector<ExecutorAddr> DeInitializers;
    for (auto &NextJD : DFSLinkOrder) {
      auto DeInitsItr = LookupResult->find(NextJD.get());
      assert(DeInitsItr != LookupResult->end() &&
             "Every JD should have at least __lljit_run_atexits");

      auto RunAtExitsItr = DeInitsItr->second.find(LLJITRunAtExits);
      if (RunAtExitsItr != DeInitsItr->second.end())
        DeInitializers.push_back(RunAtExitsItr->second.getAddress());

      for (auto &KV : DeInitsItr->second)
        if (KV.first != LLJITRunAtExits)
          DeInitializers.push_back(KV.second.getAddress());
    }

    return DeInitializers;
  }

  /// Issue lookups for all init symbols required to initialize JD (and any
  /// JITDylibs that it depends on).
  Error issueInitLookups(JITDylib &JD) {
    DenseMap<JITDylib *, SymbolLookupSet> RequiredInitSymbols;
    std::vector<JITDylibSP> DFSLinkOrder;

    if (auto Err = getExecutionSession().runSessionLocked([&]() -> Error {
          if (auto DFSLinkOrderOrErr = JD.getDFSLinkOrder())
            DFSLinkOrder = std::move(*DFSLinkOrderOrErr);
          else
            return DFSLinkOrderOrErr.takeError();

          for (auto &NextJD : DFSLinkOrder) {
            auto ISItr = InitSymbols.find(NextJD.get());
            if (ISItr != InitSymbols.end()) {
              RequiredInitSymbols[NextJD.get()] = std::move(ISItr->second);
              InitSymbols.erase(ISItr);
            }
          }
          return Error::success();
        }))
      return Err;

    return Platform::lookupInitSymbols(getExecutionSession(),
                                       RequiredInitSymbols)
        .takeError();
  }

  static void registerCxaAtExitHelper(void *Self, void (*F)(void *), void *Ctx,
                                      void *DSOHandle) {
    LLVM_DEBUG({
      dbgs() << "Registering cxa atexit function " << (void *)F << " for JD "
             << (*static_cast<JITDylib **>(DSOHandle))->getName() << "\n";
    });
    static_cast<GenericLLVMIRPlatformSupport *>(Self)->AtExitMgr.registerAtExit(
        F, Ctx, DSOHandle);
  }

  static void registerAtExitHelper(void *Self, void *DSOHandle, void (*F)()) {
    LLVM_DEBUG({
      dbgs() << "Registering atexit function " << (void *)F << " for JD "
             << (*static_cast<JITDylib **>(DSOHandle))->getName() << "\n";
    });
    static_cast<GenericLLVMIRPlatformSupport *>(Self)->AtExitMgr.registerAtExit(
        reinterpret_cast<void (*)(void *)>(F), nullptr, DSOHandle);
  }

  static void runAtExitsHelper(void *Self, void *DSOHandle) {
    LLVM_DEBUG({
      dbgs() << "Running atexit functions for JD "
             << (*static_cast<JITDylib **>(DSOHandle))->getName() << "\n";
    });
    static_cast<GenericLLVMIRPlatformSupport *>(Self)->AtExitMgr.runAtExits(
        DSOHandle);
  }

  // Constructs an LLVM IR module containing platform runtime globals,
  // functions, and interposes.
  ThreadSafeModule createPlatformRuntimeModule() {
    auto Ctx = std::make_unique<LLVMContext>();
    auto M = std::make_unique<Module>("__standard_lib", *Ctx);
    M->setDataLayout(J.getDataLayout());

    auto *GenericIRPlatformSupportTy =
        StructType::create(*Ctx, "lljit.GenericLLJITIRPlatformSupport");

    auto *PlatformInstanceDecl = new GlobalVariable(
        *M, GenericIRPlatformSupportTy, true, GlobalValue::ExternalLinkage,
        nullptr, "__lljit.platform_support_instance");

    auto *IntTy = Type::getIntNTy(*Ctx, sizeof(int) * CHAR_BIT);
    auto *BytePtrTy = PointerType::getUnqual(*Ctx);
    auto *CxaAtExitCallbackPtrTy = PointerType::getUnqual(*Ctx);

    auto *CxaAtExit = addHelperAndWrapper(
        *M, "__cxa_atexit",
        FunctionType::get(IntTy, {CxaAtExitCallbackPtrTy, BytePtrTy, BytePtrTy},
                          false),
        GlobalValue::DefaultVisibility, "__lljit.cxa_atexit_helper",
        {PlatformInstanceDecl});
    Attribute::AttrKind CxaAtExitExtAttr =
        TargetLibraryInfo::getExtAttrForI32Return(J.getTargetTriple());
    if (CxaAtExitExtAttr != Attribute::None)
      CxaAtExit->addRetAttr(CxaAtExitExtAttr);

    return ThreadSafeModule(std::move(M), std::move(Ctx));
  }

  LLJIT &J;
  std::string InitFunctionPrefix;
  std::string DeInitFunctionPrefix;
  DenseMap<JITDylib *, SymbolLookupSet> InitSymbols;
  DenseMap<JITDylib *, SymbolLookupSet> InitFunctions;
  DenseMap<JITDylib *, SymbolLookupSet> DeInitFunctions;
  /// cling-compat: per thread, one per active initialize() call, innermost
  /// last: the init functions registered while it runs.
  std::map<std::thread::id, std::vector<DenseMap<JITDylib *, SymbolLookupSet>>>
      InitFrames;
  cling::JITPlatformHooks Hooks;
  ItaniumCXAAtExitSupport AtExitMgr;
};

Error GenericLLVMIRPlatform::setupJITDylib(JITDylib &JD) {
  return S.setupJITDylib(JD);
}

Error GenericLLVMIRPlatform::teardownJITDylib(JITDylib &JD) {
  return Error::success();
}

Error GenericLLVMIRPlatform::notifyAdding(ResourceTracker &RT,
                                          const MaterializationUnit &MU) {
  return S.notifyAdding(RT, MU);
}

Expected<ThreadSafeModule>
GlobalCtorDtorScraper::operator()(ThreadSafeModule TSM,
                                  MaterializationResponsibility &R) {
  auto Err = TSM.withModuleDo([&](Module &M) -> Error {
    auto &Ctx = M.getContext();
    auto *GlobalCtors = M.getNamedGlobal("llvm.global_ctors");
    auto *GlobalDtors = M.getNamedGlobal("llvm.global_dtors");

    auto RegisterCOrDtors = [&](GlobalVariable *GlobalCOrDtors,
                                bool isCtor) -> Error {
      // If there's no llvm.global_c/dtor or it's just a decl then skip.
      if (!GlobalCOrDtors || GlobalCOrDtors->isDeclaration())
        return Error::success();
      std::string InitOrDeInitFunctionName;
      if (isCtor)
        raw_string_ostream(InitOrDeInitFunctionName)
            << InitFunctionPrefix << M.getModuleIdentifier();
      else
        raw_string_ostream(InitOrDeInitFunctionName)
            << DeInitFunctionPrefix << M.getModuleIdentifier();

      MangleAndInterner Mangle(PS.getExecutionSession(), M.getDataLayout());
      auto InternedInitOrDeInitName = Mangle(InitOrDeInitFunctionName);
      if (auto Err = R.defineMaterializing(
              {{InternedInitOrDeInitName, JITSymbolFlags::Callable}}))
        return Err;

      auto *InitOrDeInitFunc = Function::Create(
          FunctionType::get(Type::getVoidTy(Ctx), {}, false),
          GlobalValue::ExternalLinkage, InitOrDeInitFunctionName, &M);
      InitOrDeInitFunc->setVisibility(GlobalValue::HiddenVisibility);
      std::vector<std::pair<Function *, unsigned>> InitsOrDeInits;
      auto COrDtors = isCtor ? getConstructors(M) : getDestructors(M);

      for (auto E : COrDtors)
        InitsOrDeInits.push_back(std::make_pair(E.Func, E.Priority));
      llvm::stable_sort(InitsOrDeInits, llvm::less_second());

      auto *InitOrDeInitFuncEntryBlock =
          BasicBlock::Create(Ctx, "entry", InitOrDeInitFunc);
      IRBuilder<> IB(InitOrDeInitFuncEntryBlock);
      for (auto &KV : InitsOrDeInits)
        IB.CreateCall(KV.first);
      IB.CreateRetVoid();

      if (isCtor)
        PS.registerInitFunc(R.getTargetJITDylib(), InternedInitOrDeInitName);
      else
        PS.registerDeInitFunc(R.getTargetJITDylib(), InternedInitOrDeInitName);

      GlobalCOrDtors->eraseFromParent();
      return Error::success();
    };

    if (auto Err = RegisterCOrDtors(GlobalCtors, true))
      return Err;
    if (auto Err = RegisterCOrDtors(GlobalDtors, false))
      return Err;

    return Error::success();
  });

  if (Err)
    return std::move(Err);

  return std::move(TSM);
}

/// Inactive Platform Support
///
/// Explicitly disables platform support. JITDylibs are not scanned for special
/// init/deinit symbols. No runtime API interposes are injected.

} // end anonymous namespace

Expected<JITDylibSP> cling::setUpJITPlatform(LLJIT &J,
                                             JITPlatformHooks Hooks) {
  LLVM_DEBUG(
      { dbgs() << "Setting up GenericLLVMIRPlatform support for LLJIT\n"; });
  auto ProcessSymbolsJD = J.getProcessSymbolsJITDylib();
  if (!ProcessSymbolsJD)
    return make_error<StringError>(
        "Native platforms require a process symbols JITDylib",
        inconvertibleErrorCode());

  auto &PlatformJD = J.getExecutionSession().createBareJITDylib("<Platform>");
  PlatformJD.addToLinkOrder(*ProcessSymbolsJD);

  if (auto *OLL = dyn_cast<ObjectLinkingLayer>(&J.getObjLinkingLayer())) {

    bool UseEHFrames = true;

    // Enable compact-unwind support if possible.
    if (J.getTargetTriple().isOSDarwin() ||
        J.getTargetTriple().isOSBinFormatMachO()) {

      // Check if the bootstrap map says that we should force eh-frames:
      // Older libunwinds require this as they don't have a dynamic
      // registration API for compact-unwind.
      std::optional<bool> ForceEHFrames;
      if (auto Err = J.getExecutionSession().getBootstrapMapValue<bool, bool>(
              "darwin-use-ehframes-only", ForceEHFrames))
        return Err;
      if (ForceEHFrames.has_value())
        UseEHFrames = *ForceEHFrames;
      else
        UseEHFrames = false;

      // If UseEHFrames hasn't been set then we're good to use compact-unwind.
      if (!UseEHFrames) {
        if (auto UIRP =
                UnwindInfoRegistrationPlugin::Create(J.getExecutionSession())) {
          OLL->addPlugin(std::move(*UIRP));
          LLVM_DEBUG(dbgs() << "Enabled compact-unwind support.\n");
        } else
          return UIRP.takeError();
      }
    }

    // Otherwise fall back to standard unwind registration.
    if (UseEHFrames) {
      auto &ES = J.getExecutionSession();
      if (auto EHFP = EHFrameRegistrationPlugin::Create(ES)) {
        OLL->addPlugin(std::move(*EHFP));
        LLVM_DEBUG(dbgs() << "Enabled eh-frame support.\n");
      } else
        return EHFP.takeError();
    }
  }

  J.setPlatformSupport(
      std::make_unique<GenericLLVMIRPlatformSupport>(J, PlatformJD,
                                                     std::move(Hooks)));

  return &PlatformJD;
}
