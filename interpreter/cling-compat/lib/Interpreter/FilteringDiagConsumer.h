//--------------------------------------------------------------------*- C++ -*-
// cling-compat: cling's diagnostic filter for incremental inputs, moved here
// from cling's IncrementalParser.cpp.
//
// This file is dual-licensed: you can choose to license it under the University
// of Illinois Open Source License or the GNU Lesser General Public License. See
// LICENSE.TXT for details.
//------------------------------------------------------------------------------

#ifndef CLING_COMPAT_FILTERING_DIAG_CONSUMER_H
#define CLING_COMPAT_FILTERING_DIAG_CONSUMER_H

#include "cling/Utils/Diagnostics.h"

#include "clang/Basic/Diagnostic.h"
#include "clang/Basic/DiagnosticSema.h"

#include "llvm/ADT/PointerIntPair.h"

#include <stack>

namespace cling {
  using namespace clang;

  /// \brief Overrides the current DiagnosticConsumer to supress many warnings
  /// issued as a result of incremental compilation (see `HandleDiagnostic()`).
  ///
  /// Diagnostics passing the filter are, by default, forwarded to the previous
  /// DiagnosticConsumer instance.  A different consumer may be specified via
  /// `setTargetConsumer()`.
  /// In that case, given that internal state might be updated as part of
  /// `{Begin,End}SourceFile` (e.g. in TextDiagnosticPrinter), calls to such
  /// functions will be forwarded to both, the user-specified and the original
  /// consumer; however `HandleDiagnostic()` calls shall only be seen by the
  /// former.
  ///
  /// On destruction, the original (i.e. overridden consumer) is restored.
  ///
  class FilteringDiagConsumer : public cling::utils::DiagnosticsOverride {
    std::stack<bool> m_IgnorePromptDiags;
    llvm::PointerIntPair<DiagnosticConsumer*, 1, bool /*Own*/> m_Target{};

    void SyncDiagCountWithTarget() {
      NumWarnings = m_PrevClient.getNumWarnings();
      NumErrors = m_PrevClient.getNumErrors();
    }

    void BeginSourceFile(const LangOptions &LangOpts,
                         const Preprocessor *PP=nullptr) override {
      m_PrevClient.BeginSourceFile(LangOpts, PP);
      if (auto C = m_Target.getPointer())
        C->BeginSourceFile(LangOpts, PP);
    }

    void EndSourceFile() override {
      m_PrevClient.EndSourceFile();
      if (auto C = m_Target.getPointer())
        C->EndSourceFile();
      SyncDiagCountWithTarget();
    }

    void finish() override {
      m_PrevClient.finish();
      if (auto C = m_Target.getPointer())
        C->finish();
      SyncDiagCountWithTarget();
    }

    void clear() override {
      m_PrevClient.clear();
      if (auto C = m_Target.getPointer())
        C->clear();
      SyncDiagCountWithTarget();
    }

    bool IncludeInDiagnosticCounts() const override {
      return m_PrevClient.IncludeInDiagnosticCounts();
    }

    void HandleDiagnostic(DiagnosticsEngine::Level DiagLevel,
                          const Diagnostic &Info) override {
      if (Info.getID() == diag::warn_falloff_nonvoid) {
        DiagLevel = DiagnosticsEngine::Error;
      }
      if (Ignoring()) {
        if (Info.getID() == diag::warn_unused_expr
            || Info.getID() == diag::warn_unused_result
            || Info.getID() == diag::warn_unused_call
            || Info.getID() == diag::warn_unused_comparison)
          return; // ignore!
        // Unlike cling, clang-repl does not wrap inputs into void functions:
        // ext_return_has_expr ('return x;' at the prompt) is a real error
        // now and must not be suppressed (clang-repl would go on and emit
        // the invalid statement).
      }

      // In principle, for simplicity, we preserve the old behavior of
      // delivering diagnostics to just one consumer (that is why we don't emit
      // to both), but we allow the "sink" to be changed.
      // Note, however, that consumers might update their internal state in
      // calls to, e.g. `BeginSourceFile()` or `EndSourceFile()` (actually,
      // `TextDiagnosticPrinter` is an example of this), so in order to be able
      // to restore the original consumer, we need to keep forwarding these
      // calls also to `m_PrevClient` (see above).
      if (auto C = m_Target.getPointer())
        C->HandleDiagnostic(DiagLevel, Info);
      else
        m_PrevClient.HandleDiagnostic(DiagLevel, Info);
      SyncDiagCountWithTarget();
    }

    bool Ignoring() const {
      return !m_IgnorePromptDiags.empty() && m_IgnorePromptDiags.top();
    }

  public:
    FilteringDiagConsumer(DiagnosticsEngine& Diags, bool Own)
      : DiagnosticsOverride(Diags, Own) {}
    ~FilteringDiagConsumer() override { setTargetConsumer(nullptr); }

    /// \brief Sets the DiagnosticConsumer that sees `HandleDiagnostic()` calls.
    /// \param[in] Consumer - The target DiagnosticConsumer, or `nullptr` to
    ///    revert to original client.
    /// \param[in] Own - Whether we own the pointee
    ///
    void setTargetConsumer(DiagnosticConsumer* Consumer, bool Own = false) {
      if (m_Target.getInt())
        if (auto C = m_Target.getPointer())
          delete C;

      m_Target.setPointer(Consumer);
      m_Target.setInt(Own);
    }

    DiagnosticConsumer* getTargetConsumer() const
    { return m_Target.getPointer(); }

    struct RAAI {
      FilteringDiagConsumer& m_Client;
      RAAI(DiagnosticConsumer& F, bool Ignore) :
       m_Client(static_cast<FilteringDiagConsumer&>(F)) {
        m_Client.m_IgnorePromptDiags.push(Ignore);
      }
      ~RAAI() { m_Client.m_IgnorePromptDiags.pop(); }
    };
  };
} // namespace cling

#endif // CLING_COMPAT_FILTERING_DIAG_CONSUMER_H
