//------------------------------------------------------------------------------
// CLING - the C++ LLVM-based InterpreterG :)
// author:  Vassil Vassilev <vvasilev@cern.ch>
//
// This file is dual-licensed: you can choose to license it under the University
// of Illinois Open Source License or the GNU Lesser General Public License. See
// LICENSE.TXT for details.
//------------------------------------------------------------------------------

#include "cling/Interpreter/InterpreterCallbacks.h"

#include "cling/Interpreter/Interpreter.h"

#include "clang/AST/ASTContext.h"
#include "clang/Basic/ASTSourceDescriptor.h"
#include "clang/Frontend/CompilerInstance.h"
#include "clang/Lex/Preprocessor.h"
#include "clang/Sema/CodeCompleteConsumer.h"
#include "clang/Sema/MultiplexExternalSemaSource.h"
#include "clang/Sema/Sema.h"
#include "clang/Serialization/ASTReader.h"

#include <optional>

using namespace clang;

namespace cling {

  ///\brief Translates 'interesting' for the interpreter into interpreter
  /// callback.
  ///
  class InterpreterPPCallbacks : public PPCallbacks {
  private:
    cling::InterpreterCallbacks* m_Callbacks;
  public:
    InterpreterPPCallbacks(InterpreterCallbacks* C) : m_Callbacks(C) { }
    ~InterpreterPPCallbacks() override {}

    void InclusionDirective(clang::SourceLocation HashLoc,
                            const clang::Token &IncludeTok,
                            llvm::StringRef FileName,
                            bool IsAngled,
                            clang::CharSourceRange FilenameRange,
                            clang::OptionalFileEntryRef File,
                            llvm::StringRef SearchPath,
                            llvm::StringRef RelativePath,
                            const clang::Module *Imported,
                            bool ModuleImported,
                            SrcMgr::CharacteristicKind FileType) override {
      m_Callbacks->InclusionDirective(HashLoc, IncludeTok, FileName,
                                      IsAngled, FilenameRange, File,
                                      SearchPath, RelativePath, Imported,
                                      ModuleImported, FileType);
    }

    void EnteredSubmodule(clang::Module* M,
                          clang::SourceLocation ImportLoc,
                          bool ForPragma) override {
      m_Callbacks->EnteredSubmodule(M, ImportLoc, ForPragma);
    }

    bool FileNotFound(llvm::StringRef FileName) override {
      if (m_Callbacks)
        return m_Callbacks->FileNotFound(FileName);

      // Returning true would mean that the preprocessor should silently skip
      // this file.
      return false;
    }
  };

  /// \brief wraps an ExternalASTSource in an ExternalSemaSource. No functional
  /// difference between the original source and this wrapper intended.
  class ExternalASTSourceWrapper : public ExternalSemaSource {
    ExternalASTSource* m_Source;

  public:
    ExternalASTSourceWrapper(ExternalASTSource* Source) : m_Source(Source) {
      assert(m_Source && "Can't wrap nullptr ExternalASTSource");
    }

    Decl* GetExternalDecl(GlobalDeclID ID) override {
      return m_Source->GetExternalDecl(ID);
    }

    Selector GetExternalSelector(uint32_t ID) override {
      return m_Source->GetExternalSelector(ID);
    }

    uint32_t GetNumExternalSelectors() override {
      return m_Source->GetNumExternalSelectors();
    }

    Stmt* GetExternalDeclStmt(uint64_t Offset) override {
      return m_Source->GetExternalDeclStmt(Offset);
    }

    CXXCtorInitializer**
    GetExternalCXXCtorInitializers(uint64_t Offset) override {
      return m_Source->GetExternalCXXCtorInitializers(Offset);
    }

    CXXBaseSpecifier* GetExternalCXXBaseSpecifiers(uint64_t Offset) override {
      return m_Source->GetExternalCXXBaseSpecifiers(Offset);
    }

    void updateOutOfDateIdentifier(const IdentifierInfo& II) override {
      m_Source->updateOutOfDateIdentifier(II);
    }

    bool
    FindExternalVisibleDeclsByName(const DeclContext* DC, DeclarationName Name,
                                   const DeclContext* OriginalDC) override {
      return m_Source->FindExternalVisibleDeclsByName(DC, Name, OriginalDC);
    }

    bool LoadExternalSpecializations(const Decl* D, bool OnlyPartial) override {
      return m_Source->LoadExternalSpecializations(D, OnlyPartial);
    }

    bool LoadExternalSpecializations(
        const Decl* D, ArrayRef<TemplateArgument> TemplateArgs) override {
      return m_Source->LoadExternalSpecializations(D, TemplateArgs);
    }

    void completeVisibleDeclsMap(const DeclContext* DC) override {
      m_Source->completeVisibleDeclsMap(DC);
    }

    Module* getModule(unsigned ID) override { return m_Source->getModule(ID); }

    std::optional<ASTSourceDescriptor>
    getSourceDescriptor(unsigned ID) override {
      return m_Source->getSourceDescriptor(ID);
    }

    ExtKind hasExternalDefinitions(const Decl* D) override {
      return m_Source->hasExternalDefinitions(D);
    }

    void
    FindExternalLexicalDecls(const DeclContext* DC,
                             llvm::function_ref<bool(Decl::Kind)> IsKindWeWant,
                             SmallVectorImpl<Decl*>& Result) override {
      m_Source->FindExternalLexicalDecls(DC, IsKindWeWant, Result);
    }

    void FindFileRegionDecls(FileID File, unsigned Offset, unsigned Length,
                             SmallVectorImpl<Decl*>& Decls) override {
      m_Source->FindFileRegionDecls(File, Offset, Length, Decls);
    }

    void CompleteRedeclChain(const Decl* D) override {
      m_Source->CompleteRedeclChain(D);
    }

    void CompleteType(TagDecl* Tag) override { m_Source->CompleteType(Tag); }

    void CompleteType(ObjCInterfaceDecl* Class) override {
      m_Source->CompleteType(Class);
    }

    void ReadComments() override { m_Source->ReadComments(); }

    void StartedDeserializing() override { m_Source->StartedDeserializing(); }

    void FinishedDeserializing() override { m_Source->FinishedDeserializing(); }

    void StartTranslationUnit(ASTConsumer* Consumer) override {
      m_Source->StartTranslationUnit(Consumer);
    }

    void PrintStats() override { m_Source->PrintStats(); }

    bool layoutRecordType(
        const RecordDecl* Record, uint64_t& Size, uint64_t& Alignment,
        llvm::DenseMap<const FieldDecl*, uint64_t>& FieldOffsets,
        llvm::DenseMap<const CXXRecordDecl*, CharUnits>& BaseOffsets,
        llvm::DenseMap<const CXXRecordDecl*, CharUnits>& VirtualBaseOffsets)
        override {
      return m_Source->layoutRecordType(Record, Size, Alignment, FieldOffsets,
                                        BaseOffsets, VirtualBaseOffsets);
    }
  };

  ///\brief Translates 'interesting' for the interpreter ExternalSemaSource
  /// events into interpreter callbacks.
  ///
  class InterpreterExternalSemaSource : public clang::ExternalSemaSource {
  protected:
    ///\brief The interpreter callback which are subscribed for the events.
    ///
    /// Usually the callbacks is the owner of the class and the interpreter owns
    /// the callbacks so they can't be out of sync. Eg we notifying the wrong
    /// callback class.
    ///
    InterpreterCallbacks* m_Callbacks; // we don't own it.

    Sema* m_Sema; // we don't own // FIXME: once we remove ForgetSema delete.

  public:
    InterpreterExternalSemaSource(InterpreterCallbacks* C)
      : m_Callbacks(C), m_Sema(nullptr) {}

    ~InterpreterExternalSemaSource() override {
      // FIXME: Another gross hack due to the missing multiplexing AST external
      // source see Interpreter::setCallbacks.
      if (m_Sema) {
        ASTContext& C = m_Sema->getASTContext();
        // tell the owning ptr to not delete it, the callbacks will delete it.
        if (C.ExternalSource.get() == this)
          C.ExternalSource.resetWithoutRelease();
      }
    }

    void InitializeSema(Sema& S) override {
      m_Sema = &S;
    }

    void ForgetSema() override {
      m_Sema = nullptr;
    }

    InterpreterCallbacks* getCallbacks() const { return m_Callbacks; }

    /// \brief Provides last resort lookup for failed unqualified lookups.
    ///
    /// This gets translated into InterpreterCallback's call.
    ///
    ///\param[out] R The recovered symbol.
    ///\param[in] S The scope in which the lookup failed.
    ///
    ///\returns true if a suitable declaration is found.
    ///
    bool LookupUnqualified(clang::LookupResult& R, clang::Scope* S) override {
      if (m_Callbacks) {
        return m_Callbacks->LookupObject(R, S);
      }

      return false;
    }

    bool FindExternalVisibleDeclsByName(const clang::DeclContext* DC,
                                        clang::DeclarationName Name,
                                        const DeclContext *OriginalDC) override {
      if (m_Callbacks)
        return m_Callbacks->LookupObject(DC, Name);

      return false;
    }

    // Silence warning virtual function was hidden.
    using ExternalASTSource::CompleteType;
    void CompleteType(TagDecl* Tag) override {
      if (m_Callbacks)
        m_Callbacks->LookupObject(Tag);
    }

    void UpdateWithNewDeclsFwd(const DeclContext *DC, DeclarationName Name,
                               llvm::ArrayRef<NamedDecl*> Decls) {
      SetExternalVisibleDeclsForName(DC, Name, Decls);
    }
  };

  InterpreterCallbacks::InterpreterCallbacks(Interpreter* interp,
                             bool enableExternalSemaSourceCallbacks/* = false*/,
                                             bool enablePPCallbacks/* = false*/)
    : m_Interpreter(interp), m_ExternalSemaSource(nullptr), m_PPCallbacks(nullptr),
      m_IsRuntime(false) {
    Sema& SemaRef = interp->getSema();
    ASTReader* Reader = m_Interpreter->getCI()->getASTReader().get();
    ExternalSemaSource* externalSemaSrc = SemaRef.getExternalSource();
    // With C++ modules, clang can have made Sema's source a multiplexer
    // (e.g. of the module ASTReader) before we get here; add ours to it.
    if (enableExternalSemaSourceCallbacks)
      if (!externalSemaSrc || externalSemaSrc == Reader ||
          isa<MultiplexExternalSemaSource>(externalSemaSrc)) {
        // If the ExternalSemaSource is the PCH reader we still need to insert
        // our listener.
        m_ExternalSemaSource = new InterpreterExternalSemaSource(this);
        m_ExternalSemaSource->InitializeSema(SemaRef);
        m_Interpreter->getSema().addExternalSource(m_ExternalSemaSource);

        // Overwrite the ExternalASTSource.
        clang::ASTContext& Ctx = SemaRef.getASTContext();
        auto ExistingSource = Ctx.getExternalSource();
        // If we already have source, we need to create a multiplexer with the
        // existing source.
        if (ExistingSource) {
          // Make sure the context is not deleting the existing source.
          // FIXME: We should delete this, but looking at the other TODO's in
          // the destructor we can't easily free our callbacks from here...
          Ctx.ExternalSource.resetWithoutRelease();

          // Wrap the existing source in a wrapper so that it becomes an
          // external sema source. This way we can use the existing multiplexer
          // for this.
          auto wrapper = new ExternalASTSourceWrapper(ExistingSource);

          // Wrap both the existing source and our source. We give our own
          // source preference to the existing one.
          IntrusiveRefCntPtr<ExternalASTSource> S;
          S = new MultiplexExternalSemaSource(m_ExternalSemaSource, wrapper);

          Ctx.setExternalSource(S);
        } else {
          // We don't have an existing source, so just set our own source.
          Ctx.setExternalSource(m_ExternalSemaSource);
        }
    }

    if (enablePPCallbacks) {
      Preprocessor& PP = m_Interpreter->getCI()->getPreprocessor();
      m_PPCallbacks = new InterpreterPPCallbacks(this);
      PP.addPPCallbacks(std::unique_ptr<InterpreterPPCallbacks>(m_PPCallbacks));
    }
  }

  // pin the vtable here
  InterpreterCallbacks::~InterpreterCallbacks() {
    // FIXME: we have to remove the external source at destruction time. Needs
    // further tweaks of the patch in clang. This will be done later once the
    // patch is in clang's mainline.
  }

  void InterpreterCallbacks::SetIsRuntime(bool val) {
    m_IsRuntime = val;
  }

  ExternalSemaSource*
  InterpreterCallbacks::getInterpreterExternalSemaSource() const {
    return m_ExternalSemaSource;
  }

  bool InterpreterCallbacks::FileNotFound(llvm::StringRef) {
    // Default implementation is no op.
    return false;
  }

  bool InterpreterCallbacks::LookupObject(LookupResult&, Scope*) {
    // Default implementation is no op.
    return false;
  }

  bool InterpreterCallbacks::LookupObject(const DeclContext*, DeclarationName) {
    // Default implementation is no op.
    return false;
  }

  bool InterpreterCallbacks::LookupObject(TagDecl*) {
    // Default implementation is no op.
    return false;
  }

  void InterpreterCallbacks::UpdateWithNewDecls(const DeclContext *DC,
                                                DeclarationName Name,
                                             llvm::ArrayRef<NamedDecl*> Decls) {
    if (m_ExternalSemaSource)
      m_ExternalSemaSource->UpdateWithNewDeclsFwd(DC, Name, Decls);
  }
} // end namespace cling
