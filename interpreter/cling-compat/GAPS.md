# cling-compat: gaps between cling and clang-repl

cling-compat implements the part of cling's API that ROOT uses on top of
`clang::Interpreter` (clang-repl). This file lists what that took: patches to
clang and LLVM, corners that were cut, and features that do not work yet.
Each entry is a candidate for an upstream patch or a ROOT-side decision.

Build with `-Dclang_repl=ON` (and `-Dclad=OFF`, `-Dcuda=OFF` for now).

## Status

Full ctest run (2026-10-02, Linux x86-64, system LLVM 22.1.8 with the bundled
clang): 3582 of 3637 tests pass. With cling, 3607 of 3638 pass.

The largest groups of remaining failures (from the triage of 2026-10-01; the
crashes of RooFit `rf104`, TMVA SOFIE and `runlongArgTest` are fixed since):

| Tests | Cause |
| --- | --- |
| 10 | Crashes: null dereference protection is missing (3), single causes (RooFit `rf104`, TMVA SOFIE, `runlongArgTest`, ...) |
| 7 | Wrong exit code |
| 4 | JupyROOT: the kernel does not respond |
| 3 | Redefinition at the prompt |
| 2 | Template autoloading (`TmpltNoSpec`) |

## What the layer contains

| Part | Origin | Why it is still needed |
| --- | --- | --- |
| `Interpreter.cpp`, `DeclCollector.cpp` | new | Adapts cling's API (transactions, callbacks, values, unloading) to clang-repl |
| `JITPlatform.cpp` | LLVM, one fix | LLJIT's generic IR platform with re-entrant initialization (see below) |
| `DynamicLookup.cpp` | cling, ported | Dynamic scopes: statements with names resolved at runtime |
| `CIFactory.cpp` | cling, cut down | Compiler arguments, host include discovery, modulemap overlays for libc/std |
| `FilteringDiagConsumer.h` | cling (IncrementalParser.cpp) | Owns the original diagnostic client (rootcling relies on that), filters prompt noise, supports `replaceDiagnosticConsumer` |
| `LookupHelper.cpp`, `Utils/AST.cpp` | cling, unchanged | Name lookup and type normalization used by TCling and rootcling |
| `MetaProcessor/*` | cling, small changes | The `.L`, `.x`, `.q`, ... prompt commands; unnamed macros |
| `DynamicLibraryManager*.cpp` | cling, unchanged | Library search and symbol-to-library resolution for autoloading |
| `Value.cpp`, `ValuePrinter.cpp` | cling, small change | `cling::Value` and the `printValue` machinery ROOT's types overload |
| `ForwardDeclPrinter.cpp` | cling, small change | Forward declarations for rootmap files (rootcling) |
| `Transaction.cpp` | cling, small change | The per-input record of declarations that TCling consumes |

Not carried over from cling: IncrementalParser, IncrementalExecutor,
IncrementalJIT, BackendPasses, most AST transformers (wrapping, value
extraction and printing synthesis, `auto` synthesis, null dereference checks,
definition shadowing, decl extraction), DeclUnloader, TransactionUnloader,
CUDA device compiler, UserInterface.

## Clang patches (interpreter/llvm-project/clang)

All are candidates for upstream clang-repl.

1. `Interpreter::create(CI, IEB, ASTConsumerFactory)`: lets a client add an
   ASTConsumer created once the Preprocessor exists (rootcling needs a
   `PCHGenerator` there). Also fixes that passing a consumer to the
   protected constructor silently disabled value printing.
2. `Interpreter::getParser()`: LookupHelper parses type and function names on
   the side with the incremental parser; a second `Parser` cannot share the
   `Preprocessor`.
3. `Interpreter::takeLastValue(Value&)`: lets clients call `Parse()` and
   `Execute()` separately (ROOT must see new declarations before static
   initializers run) and still get the value of a trailing expression.
4. Value printing: a `void` expression without `;` (`f()`) was replaced by
   the value-capturing call and never evaluated. Now `(f(), capture)`
   (llvm/llvm-project#219800).
5. `ASTWriter`: writing a PCH or module from an incremental
   CompilerInstance. All TranslationUnitDecls (one per input) are the
   predefined TU, the decls of all of them are written, and the TU lookup
   table comes from the primary context.
6. `ASTReader`: declarations read from AST files belong to the primary
   TranslationUnitDecl, not to whichever input's TU is current.
7. Typeinfo and vtables (`_ZT*`) are made weak before a module goes to the
   JIT: a class whose key function is defined in a later input, or that a
   dictionary library defines too, gave duplicate definitions (cling had
   the same hack).
8. A `TopLevelStmtDecl` never gets a null statement (an invalid statement
   becomes a null statement); several clients crashed on it. Same code in
   LLVM main.
9. For clients that resolve names of top-level statements themselves:
   - the first name of a statement is looked up again in the statement's
     scope (the declaration/statement disambiguation had marked `f` in
     `f();` as undeclared at file scope);
   - type-dependent expressions are left to the client instead of being
     value-printed;
   - `Interpreter::setCaptureValueOfExprStmts()` captures (not prints) the
     value of an expression statement that ends in `;` (`ProcessLine("&f;")`).
10. Include stacks end at clang-repl's main file in incremental mode: no
    "In file included from <<< inputs >>>:1:" before each diagnostic.
11. `Interpreter::Undo()` makes the declarations of the undone input
    unreachable: they were found again through the scope chain, rebuilt
    lookup tables, enumerators of unscoped enums, and members of reopened
    namespaces. `DeclContext::removeDecl()` tolerates a missing lookup
    entry. `Interpreter::getCodeGen()` lets clients make CodeGen forget
    undone definitions. Like upstream's `ASTDeclUnmerger`, the undone
    redeclarations are withdrawn from their redeclaration chains (and a class
    loses an undone definition), and the declaration they had replaced in the
    lookup tables and the scope chain is found again (`gSystem->Unload()` of a
    library whose payload forward-declared a class defined earlier). The
    tables are updated without the external source, which in ROOT autoparses.
12. Local entities of top-level statements (e.g. the closure type of a
    lambda) are mangled as locals of a function `__stmt__N`; with
    optimization, the second statement with a lambda crashed. Backport from
    LLVM main.
13. The value of an input is kept per thread (`getLastValueSlot()`): ROOT
    runs JIT-ed code of several threads at once, and their `Calc()` results
    overwrote each other.
14. `clang::Value`'s move constructor no longer releases the storage it took
    over (heap corruption with values of class type). Backport from LLVM
    main.

## LLVM patch (carried as a copy)

ROOT may take LLVM from the system, so this one cannot live in
`interpreter/llvm-project/llvm`. `JITPlatform.cpp` is a copy of LLJIT's
`GenericLLVMIRPlatformSupport` with the fix; it is proposed upstream as
`clang-repl/pending/06-lljit-nested-initialize.patch` (outside the
repository). Remove the copy once ROOT's LLVM has it.

- **Re-entrant `LLJIT::initialize()`.** Linking an input can make ROOT's
  autoloading generator load a library whose initializers declare and run
  more code (dictionary registration, `TROOT::ProcessLine()`; ROOT-10426).
  The nested `initialize()` waited for the init functions of the outer input,
  which were still being materialized: a deadlock. Each call now records the
  init functions registered while it materializes, and a nested call runs
  only those.

The copy also has a hook to release ROOT's interpreter lock while the
initializers run (they are the user's statements), but not while they are
compiled: another thread would generate code in the same `LLVMContext`.

## Things ROOT had to adapt to

- **One TranslationUnitDecl per input.** Code that iterates "the" TU must
  iterate all of them: `TMetaUtils::GetAllTranslationUnits()`, used by the
  dictionary scanners, TCling's module index and its replay of existing
  declarations. More places will need it (TClingClassInfo and
  TClingTypedefInfo iterators, `.class`/`.typedef` listings).
- **Top-level statements are not in functions.** TClingCallbacks treats a
  `TopLevelStmtDecl`, and (clang-repl only) the initializer of a global
  declared at the prompt, like cling's wrapper functions: ROOT special
  objects and dynamic scopes work there.

## Workarounds in cling-compat (each one a gap)

- **Locals in `evaluate()`.** cling compiled `Calc()` inputs into a function,
  so `Calc("R__READ_LOCKGUARD(m); f(10)")` declared a local guard. In
  clang-repl, top-level statements have no local scope: the guard would be a
  global that never unlocks. Inputs with several statements are put into
  `[]() -> decltype(auto) { ...; return last; }()`.

- **JIT setup.** clang-repl builds its JIT for the host CPU; CodeGen uses the
  CompilerInstance's target, and the mismatch crashed instruction selection.
  cling-compat configures the JIT like cling: same target, large code model
  on x86-64 (weak symbols such as `DW.ref.__gxx_personality_v0` ended up out
  of range of 32-bit relocations), and a memory manager that never frees
  code (ROOT keeps pointers to JIT-ed wrappers; root-project/root#10898).
- **Runtime symbols.** clang-repl expects the host process to export its
  value-capturing runtime (`__clang_Interpreter_SetValue*`, a placement
  `operator new`). In libCling they are hidden; cling-compat defines them
  in the JIT. The JIT also gets its own `__cxa_atexit`, which records the
  destructors per input (see unloading).
- **Shutdown order.** clang-repl runs the JIT's deinitializers when it is
  destroyed, which is after LLVM's CodeGen statics are gone: compiling the
  platform's code then failed at random. cling-compat compiles it up front,
  and runs the static destructors of JIT-ed code before ROOT's end-of-process
  cleanups (like cling).
- **Unique source locations.** clang caches the last module import by its
  location (`CompilerInstance::loadModule`); synthesized imports must each
  get a new location. cling-compat hands out offsets in a virtual buffer.
- **Emitting code outside of an input.** Declarations handed to CodeGen
  outside `Parse()` (cling's value printer, ROOT's special objects) only
  reach the JIT with the next input. cling-compat parses an empty input to
  flush when it cannot find a symbol. Pending instantiations queued outside
  `Parse()` are stashed away by clang-repl's per-input instantiation scope;
  cling-compat performs them when a transaction ends.
- **Generator order.** ORC asks a JITDylib's generators before the next
  JITDylib; ROOT's library autoloading generator (which can scan every
  library on the system) must go into a JITDylib at the end of the link
  order, and be detached before clang-repl deinitializes the JIT.
- **ExternalSemaSource.** With C++ modules, Sema's external source is already
  a `MultiplexExternalSemaSource`; ROOT's callbacks (autoloading, special
  objects, dynamic scopes) are added to it.
- **Declarations without `;`** (`auto h = new TH1F(...)`) and a bare
  `return` at the prompt: clang-repl rejects them; cling-compat retries
  with `;` and prints the declared variable.
- **Unnamed macros.** Without `return`, their statements stay top-level code
  (their variables are globals, as with cling's declaration extraction).
  With `return`, they become a wrapper function with a deduced return type;
  `-Wreturn-type` is off for it and falling off the end returns a
  value-initialized result (cling's wrappers were `void`).
- **Diagnostics layout.** cling compiled most inputs inside a wrapper
  function, so their diagnostics showed the input on the next line,
  indented by one space; roottest's reference files compare that.
  cling-compat inserts `"\n "` where cling inserted its wrapper.
- **Diagnostics ownership.** rootcling replaces the DiagnosticConsumer and
  expects the engine not to own the original one (cling's filter owned it).
- **Forward declarations** for rootmap files: the include walk stops at
  clang-repl's inputs (it reached the main file, `__main_argc_argv`).
- **Transactions are not nested.** An input is the current transaction
  while its code runs (TCling::AutoParse() relies on it); transactions made
  by lookups outside an input have no input of their own and are dropped
  when unloading.

## Features ported from cling

- **Dynamic scopes** (`DynamicLookup.cpp`, cling's EvaluateTSynthesizer):
  statements with names that only exist at runtime are compiled and run at
  runtime, in functions (unnamed macros with `return`), top-level statements
  and initializers of globals. Fixed on the way (cling has these bugs):
  `std::cout << unknown` printed as `std::((*(ostream*)0x...))`; loop
  conditions; `return unknown();` in a function with a deduced return type.
  Not supported: a global of class type or `auto` initialized from an unknown
  name (cling's LifetimeHandler only works for locals).
- **ROOT special objects** (`ntuple->Draw()` for an object in a file), at the
  prompt, in macros and in initializers of globals.
- **Unloading** (`.U`, `.undo`, `TSystem::Unload()`, ACLiC reloading). As
  with cling, unloading a file reverts all inputs after it, newest first,
  which is what clang-repl's `Undo()` can do. Like cling's
  TransactionUnloader: macros are undefined (include guards), cached file
  contents dropped, CodeGen forgets the definitions, and the destructors of
  the input's globals run. `.undo N` counts inputs.

## Missing features

- **Unloading a transaction that is not the last one** (cling's
  TransactionUnloader could). ROOT does not need it so far. `UnloadDecl`
  (used by LookupHelper to remove failed template instantiations) only marks
  the declaration invalid.
- **Null dereference protection** (cling's NullDerefProtectionTransformer):
  dereferencing a null pointer in interpreted code crashes instead of
  raising an error.
- **Redefinition at the prompt** (`AllowRedefinition`, cling's
  DefinitionShadower): not implemented.
- **`return` in top-level code** other than in unnamed macros and a bare
  `return` at the prompt.
- **Value ownership.** A `cling::Value` of class type returned by `process()`
  or `evaluate()` refers to storage owned by a `clang::Value` that is kept
  alive until the interpreter is destroyed (a leak per such value).
- **Values with several threads.** clang-repl captures the value of an
  input in one member of `clang::Interpreter`; inputs that run concurrently
  (without ROOT's lock) can overwrite each other's value.
- **`parse()`** generates code anyway (clang-repl always runs CodeGen); the
  module is just not executed.
- **Code completion** (`codeComplete`) is not ported (`roottest-root-rint-TabCom`).
- **`.storeState` / `.compareState`** are not supported.
- **Plugins** (clad) are not loaded: `-Dclad=OFF`.
- **CUDA**: not ported: `-Dcuda=OFF`.
- **Declarations made through CppInterOp** (PyROOT's `cppdef`) bypass the
  transactions, so TCling is not told about them.
