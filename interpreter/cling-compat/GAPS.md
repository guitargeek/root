# cling-compat: gaps between cling and clang-repl

cling-compat implements the part of cling's API that ROOT uses on top of
`clang::Interpreter` (clang-repl). This file lists what that took: patches to
clang and LLVM, corners that were cut, and features that do not work yet.
Each entry is a candidate for an upstream patch or a ROOT-side decision.

Build with `-Dclang_repl=ON` (and `-Dclad=OFF`, `-Dcuda=OFF` for now).

## Status

Full ctest runs (Linux x86-64, system LLVM 22.1.8 with the bundled clang,
both builds on master ad94d01ad4f in the same environment): 3609 of 3639
tests pass with clang-repl (2026-10-06), 3632 of 3639 with cling
(2026-10-05). 20 tests fail only with clang-repl (the run had 4 more
failures, timing and network, that pass when rerun):

| Tests | Cause |
| --- | --- |
| 5 | Redefinition at the prompt (cling's shadowing) |
| 3 | `MakeProject` |
| 1 | Implicit declaration of `x = expr` at the prompt |
| 11 | Single causes (code completion, function bodies, crashes in `foam_demopers` and `multidimfit`, `selectorreader-Classes`, ...) |

## What the layer contains

| Part | Origin | Why it is still needed |
| --- | --- | --- |
| `Interpreter.cpp`, `DeclCollector.cpp` | new | Adapts cling's API (transactions, callbacks, values, unloading) to clang-repl |
| `JITPlatform.cpp` | LLVM, one fix | LLJIT's generic IR platform with re-entrant initialization (see below) |
| `DynamicLookup.cpp` | cling, ported | Dynamic scopes: statements with names resolved at runtime |
| `NullDerefProtection.cpp` | cling, small change | `--ptrcheck`: errors instead of crashes on null dereferences |
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
extraction and printing synthesis, `auto` synthesis, definition shadowing,
decl extraction), DeclUnloader, TransactionUnloader,
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
15. Value capture keeps the expression's cleanups (`ExprWithCleanups`): the
    temporaries of a captured statement (e.g. `h(S())`) were destroyed at
    the end of the function running the top-level statements, and CodeGen
    no longer emitted their destructors (JIT: symbol not found). Upstream
    has the same bug for printed values (main asserts in
    `moveLazyEmissionStates`); ROOT captures every statement.
16. Nested inputs: ROOT parses more code while an input is parsed (e.g. a
    dictionary payload when a template instantiation needs a class
    complete). The nested input took CodeGen's module, with the code of the
    outer input so far, and ran it: it referred to definitions the outer
    input had not emitted yet. Now a nested input neither finishes nor
    takes the module; its code runs with the outer input's, as cling ran
    nested transactions with their parent.
17. `Interpreter::getPTUs()`: not every input has a cling transaction
    (CppInterOp compiles its wrappers with clang-repl directly). To unload a
    transaction, cling-compat undoes the inputs that came after it as well,
    and needs to know how many there are.
18. `Interpreter::setCodeGenFilter()`: decides which declarations reach
    CodeGen. cling-compat uses it like cling's DeclCollector: the
    non-inline functions and variables of a dictionary payload are not
    generated again, as the library has them (the JIT-ed copy shadowed the
    library's, e.g. after ACLiC compiled a macro that was interpreted
    before).

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
- **Rolled-back transactions.** clang-repl cannot remove the declarations of
  a transaction that failed without an input of its own (e.g. LookupHelper
  instantiating `pair<string, Inner<int>>` while `Inner<int>` is incomplete).
  Like cling's TransactionUnloader, cling-compat removes the class,
  variable and function template specializations created in it from their
  templates (it records them as an ASTMutationListener), and failed class
  template instantiations; the next use instantiates them again. Function
  specializations must go too: they refer to the class specializations made
  with them (e.g. `__is_complete_or_unbounded(__type_identity<T>)`), and the
  next use, with a new class specialization, would no longer match them. Like cling, it also
  clears clang's cache of constraint satisfaction
  (`UnsubstitutedConstraintSatisfactionCache`) on rollback and unloading: it
  kept `is_default_constructible_v<Inner<int>>` false after `Inner<int>`
  was defined. The declarations themselves stay in the AST.
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
- **Null dereference protection** (`NullDerefProtection.cpp`, cling's
  NullDerefProtectionTransformer): with `root --ptrcheck`, dereferencing an
  invalid pointer in user code raises an error instead of crashing. Unlike
  cling's inputs, clang-repl's are memory buffers without a file, and a
  top-level statement whose value is captured has no begin location; both
  count as user code.
- **Unresolved symbols** (cling's `diagnoseUnresolvedSymbols()`): an input
  that uses a symbol nobody defines reports it with cling's message ("You
  are probably missing the definition of ...") and is rolled back (unloaded)
  like with cling, if it is still the last input. The JIT platform forgets
  the initializers of inputs that failed to link; they would fail every
  later input again. Unlike cling, a function compiled while the symbol was
  missing stays broken in the JIT after the symbol gets defined.
- **ROOT special objects** (`ntuple->Draw()` for an object in a file), at the
  prompt, in macros and in initializers of globals.
- **Unloading** (`.U`, `.undo`, `TSystem::Unload()`, ACLiC reloading). As
  with cling, unloading a file reverts all inputs after it, newest first,
  which is what clang-repl's `Undo()` can do. Like cling's
  TransactionUnloader: macros are undefined (include guards), cached file
  contents dropped, CodeGen forgets the definitions, and the destructors of
  the input's globals run. Inputs without a transaction after the unloaded
  one (CppInterOp's wrappers, `flushPendingCodeGen()`) are undone with it.
  Implicit instantiations of templates the input did not declare stay in
  the AST, so CodeGen keeps them to emit them again on their next use.
  `.undo N` counts inputs.

## Missing features

- **Unloading a transaction that is not the last one** (cling's
  TransactionUnloader could). ROOT does not need it so far. `UnloadDecl`
  (used by LookupHelper to remove failed template instantiations) only marks
  the declaration invalid, except for class template instantiations (see
  "Rolled-back transactions").
- **Redefinition at the prompt** (`AllowRedefinition`, cling's
  DefinitionShadower): not implemented.
- **`return` in top-level code** other than in unnamed macros and a bare
  `return` at the prompt.
- **Value ownership.** A `cling::Value` of class type returned by `process()`
  or `evaluate()` refers to storage owned by a `clang::Value` that is kept
  alive until the interpreter is destroyed (a leak per such value).
- **`parse()`** generates code anyway (clang-repl always runs CodeGen); the
  module is just not executed.
- **Code completion** (`codeComplete`) is not ported (`roottest-root-rint-TabCom`).
- **`.storeState` / `.compareState`** are not supported.
- **Plugins** (clad) are not loaded: `-Dclad=OFF`.
- **CUDA**: not ported: `-Dcuda=OFF`.
- **Declarations made through CppInterOp** (PyROOT's `cppdef`) bypass the
  transactions, so TCling is not told about them.
