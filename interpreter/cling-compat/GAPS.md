# cling-compat: gaps between cling and clang-repl

cling-compat implements the part of cling's API that ROOT uses on top of
`clang::Interpreter` (clang-repl). This file lists what that took: patches to
clang, corners that were cut, and features that do not work yet. Each entry
is a candidate for an upstream patch or a ROOT-side decision.

Build with `-Dclang_repl=ON` (and `-Dclad=OFF`, `-Dcuda=OFF` for now).

## What the layer contains

| Part | Origin | Why it is still needed |
| --- | --- | --- |
| `Interpreter.cpp`, `DeclCollector.cpp` | new | Adapts cling's API (transactions, callbacks, values) to clang-repl |
| `CIFactory.cpp` | cling, cut down | Compiler arguments, host include discovery, modulemap overlays for libc/std |
| `FilteringDiagConsumer.h` | cling (IncrementalParser.cpp) | Owns the original diagnostic client (rootcling relies on that), filters prompt noise, supports `replaceDiagnosticConsumer` |
| `LookupHelper.cpp`, `Utils/AST.cpp` | cling, unchanged | Name lookup and type normalization used by TCling and rootcling |
| `MetaProcessor/*` | cling, unchanged | The `.L`, `.x`, `.q`, ... prompt commands |
| `DynamicLibraryManager*.cpp` | cling, unchanged | Library search and symbol-to-library resolution for autoloading |
| `Value.cpp`, `ValuePrinter.cpp` | cling, small change | `cling::Value` and the `printValue` machinery ROOT's types overload |
| `ForwardDeclPrinter.cpp` | cling, unchanged | Forward declarations for rootmap files (rootcling) |
| `Transaction.cpp` | cling, small change | The per-input record of declarations that TCling consumes |

Not carried over from cling: IncrementalParser, IncrementalExecutor,
IncrementalJIT, BackendPasses, all AST transformers (wrapping, value
extraction and printing synthesis, dynamic lookup, `auto` synthesis, null
dereference checks, definition shadowing, decl extraction), DeclUnloader,
TransactionUnloader, CUDA device compiler, UserInterface.

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
   the value-capturing call and never evaluated. Now `(f(), capture)`.
5. `ASTWriter`: writing a PCH or module from an incremental
   CompilerInstance. All TranslationUnitDecls (one per input) are the
   predefined TU, the decls of all of them are written, and the TU lookup
   table comes from the primary context.
6. `ASTReader`: declarations read from AST files belong to the primary
   TranslationUnitDecl, not to whichever input's TU is current.

## Things ROOT had to adapt to

- **One TranslationUnitDecl per input.** Code that iterates "the" TU must
  iterate all of them: `TMetaUtils::GetAllTranslationUnits()`, used by the
  dictionary scanners, TCling's module index and its replay of existing
  declarations. More places will need it (TClingClassInfo and
  TClingTypedefInfo iterators, `.class`/`.typedef` listings).

## Workarounds in cling-compat (each one a gap)

- **JIT target.** clang-repl builds its JIT for the host CPU; CodeGen uses
  the CompilerInstance's target. The mismatch crashed instruction selection.
  cling-compat passes an `IncrementalExecutorBuilder` configured like cling.
- **Runtime symbols.** clang-repl expects the host process to export its
  value-capturing runtime (`__clang_Interpreter_SetValue*`, a placement
  `operator new`). In libCling they are hidden; cling-compat defines them
  in the JIT.
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
- **Declarations without `;`** (`auto h = new TH1F(...)`) and a bare
  `return` at the prompt: clang-repl rejects them; cling-compat retries
  with `;` and prints the declared variable.
- **Diagnostics ownership.** rootcling replaces the DiagnosticConsumer and
  expects the engine not to own the original one (cling's filter owned it).

## Missing features

- **Selective unloading.** `unload(T)` only works for the most recent input
  (`clang::Interpreter::Undo`). `UnloadDecl` (used by LookupHelper to remove
  failed template instantiations) only marks the declaration invalid.
- **`return` in top-level code** other than a bare `return` (e.g. in unnamed
  macros): cling wrapped inputs into functions.
- **Dynamic scopes.** `enableDynamicLookup()` is accepted but does nothing.
- **ROOT "special objects"** (`h->Draw()` for a histogram `h` in a file):
  TClingCallbacks only resolves them inside function bodies; clang-repl's
  top-level statements are not in one.
- **Redefinition at the prompt** (`AllowRedefinition`): not implemented.
- **Value ownership.** A `cling::Value` of class type returned by `process()`
  or `evaluate()` refers to storage owned by a `clang::Value` that is kept
  alive until the interpreter is destroyed (a leak per such value).
- **Static destructors / atexit.** `runAtExitFuncs()` only runs functions
  registered through `AddAtExitFunc`; destructors of JIT-ed globals run when
  clang-repl's LLJIT is torn down.
- **`parse()`** generates code anyway (clang-repl always runs CodeGen); the
  module is just not executed.
- **Code completion** (`codeComplete`) is not ported.
- **`.storeState` / `.compareState`** are not supported.
- **Plugins** (clad) are not loaded: `-Dclad=OFF`.
- **CUDA**: not ported: `-Dcuda=OFF`.
- **Declarations made through CppInterOp** (PyROOT's `cppdef`) bypass the
  transactions, so TCling is not told about them.
