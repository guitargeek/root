// @(#)root/utils/src:$Id$
// Author: Philippe Canal November 2011 ; originated from Zdenek Culik   16/04/2010 and Velislava Spasova.

/*************************************************************************
 * Copyright (C) 1995-2011, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/rootcint.            *
 *************************************************************************/

#include "Scanner.h"
#include "clang/AST/ASTConsumer.h"
#include "clang/AST/DeclLookups.h"
#include "clang/Basic/SourceLocation.h"
#include "clang/Basic/SourceManager.h"
#include "llvm/ADT/SmallPtrSet.h"
#include "llvm/ADT/SmallSet.h"
#include "llvm/ADT/SmallVector.h"
#include "clang/Sema/Sema.h"
#include "clang/Frontend/CompilerInstance.h"

#include "cling/Interpreter/Interpreter.h"
#include "llvm/Support/Path.h"

#include "TClassEdit.h"

#include <algorithm>
#include <cctype>
#include <iostream>
#include <sstream> // class ostringstream
#include <type_traits>
#include "llvm/ADT/StringExtras.h"

#include "SelectionRules.h"

namespace {

   class RPredicateIsSameNamespace
   {
   private:
      clang::NamespaceDecl *fTarget;
   public:
      RPredicateIsSameNamespace(clang::NamespaceDecl *target) : fTarget(target) {}

      bool operator()(const RScanner::AnnotatedNamespaceDecl& element)
      {
         return (fTarget == element);
      }
   };

template<class T>
inline static bool IsElementPresent(const std::vector<T> &v, const T &el){
   return std::find(v.begin(),v.end(),el) != v.end();
}

template<class T>
inline static bool IsElementPresent(const std::vector<const T*> &v, T *el){
   return std::find(v.begin(),v.end(),el) != v.end();
}

}

using namespace ROOT;
using namespace clang;

const char* RScanner::fgClangDeclKey = "ClangDecl"; // property key used for connection with Clang objects
const char* RScanner::fgClangFuncKey = "ClangFunc"; // property key for demangled names

int RScanner::fgAnonymousClassCounter = 0;
int RScanner::fgBadClassCounter = 0;
int RScanner::fgAnonymousEnumCounter  = 0;

std::map <clang::Decl*, std::string> RScanner::fgAnonymousClassMap;
std::map <clang::Decl*, std::string> RScanner::fgAnonymousEnumMap;

////////////////////////////////////////////////////////////////////////////////
/// Regular constructor setting up the scanner to search for entities
/// matching the 'rules'.

RScanner::RScanner (SelectionRules &rules,
                    EScanType stype,
                    const cling::Interpreter &interpret,
                    ROOT::TMetaUtils::TNormalizedCtxt &normCtxt,
                    unsigned int verbose /* = 0 */) :
  fVerboseLevel(verbose),
  fSourceManager(nullptr),
  fInterpreter(interpret),
  fRecordDeclCallback(nullptr),
  fNormCtxt(normCtxt),
  fSelectionRules(rules),
  fScanType(stype),
  fFirstPass(true)
{
   // Build the cache for all selection rules
   fSelectionRules.FillCache();
   CollectMemberRulePrefixes();

   for (int i = 0; i <= fgDeclLast; i ++)
      fDeclTable [i] = false;

   for (int i = 0; i <= fgTypeLast; i ++)
      fTypeTable [i] = false;

   fLastDecl = nullptr;
}

////////////////////////////////////////////////////////////////////////////////

RScanner::~RScanner ()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Whether we can actually visit this declaration, i.e. if it is reachable
/// via name lookup.
///
/// RScanner shouldn't touch decls for which this method returns false as we
/// call Sema methods on those declarations. Those will fail in strange way as
/// they assume those decls are already visible.
///
/// The main problem this is supposed to prevent is when we use C++ modules and
/// have hidden declarations in our AST. Usually they can't be found as they are
/// hidden from name lookup until their module is actually imported, but as the
/// RecursiveASTVisitor is not supposed to be restricted by lookup limitations,
/// it still reaches those hidden declarations.
bool RScanner::shouldVisitDecl(clang::NamedDecl *D)
{
   if (auto M = D->getOwningModule()) {
      return fInterpreter.getSema().isModuleVisible(M);
   }
   return true;
}

////////////////////////////////////////////////////////////////////////////////

inline void* ToDeclProp(clang::Decl* item)
{
   /* conversion and type check used by AddProperty */
   return item;
}

////////////////////////////////////////////////////////////////////////////////

inline size_t APIntToSize(const llvm::APInt& num)
{
   return *num.getRawData();
}

////////////////////////////////////////////////////////////////////////////////

inline long APIntToLong(const llvm::APInt& num)
{
   return *num.getRawData();
}

////////////////////////////////////////////////////////////////////////////////

inline std::string APIntToStr(const llvm::APInt& num)
{
   return llvm::toString(num, /*radix*/10, /*signed*/true);
}

////////////////////////////////////////////////////////////////////////////////

inline std::string IntToStr(int num)
{
   std::string txt = "";
   txt += num;
   return txt;
}

////////////////////////////////////////////////////////////////////////////////

inline std::string IntToStd(int num)
{
   std::ostringstream stream;
   stream << num;
   return stream.str();
}

////////////////////////////////////////////////////////////////////////////////

inline std::string Message(const std::string &msg, const std::string &location)
{
   std::string loc = location;

   if (loc == "")
      return msg;
   else
      return loc + " " + msg;
}

////////////////////////////////////////////////////////////////////////////////

void RScanner::ShowInfo(const std::string &msg, const std::string &location) const
{
   const std::string message = Message(msg, location);
   std::cout << message << std::endl;
}

////////////////////////////////////////////////////////////////////////////////

void RScanner::ShowWarning(const std::string &msg, const std::string &location) const
{
   const std::string message = Message(msg, location);
   std::cout << message << std::endl;
}

////////////////////////////////////////////////////////////////////////////////

void RScanner::ShowError(const std::string &msg, const std::string &location) const
{
   const std::string message = Message(msg, location);
   std::cout << message << std::endl;
}

////////////////////////////////////////////////////////////////////////////////

void RScanner::ShowTemplateInfo(const std::string &msg, const std::string &location) const
{
   std::string loc = location;
   if (loc == "")
      loc = GetLocation (fLastDecl);
   ShowWarning(msg, loc);
}

////////////////////////////////////////////////////////////////////////////////

std::string RScanner::GetSrcLocation(clang::SourceLocation L) const
{
   std::string location = "";
   llvm::raw_string_ostream stream(location);
   L.print(stream, *fSourceManager);
   return stream.str();
}

////////////////////////////////////////////////////////////////////////////////

std::string RScanner::GetLocation(clang::Decl* D) const
{
   if (D == nullptr)
   {
      return "";
   }
   else
   {
      std::string location = "";
      llvm::raw_string_ostream stream(location);
      D->getLocation().print(stream, *fSourceManager);
      return stream.str();
   }
}

////////////////////////////////////////////////////////////////////////////////

std::string RScanner::GetName(clang::Decl* D) const
{
   std::string name = "";
   // std::string kind = D->getDeclKindName();

   if (clang::NamedDecl* ND = dyn_cast <clang::NamedDecl> (D)) {
      name = ND->getQualifiedNameAsString();
   }

   return name;
}

////////////////////////////////////////////////////////////////////////////////

inline std::string AddSpace(const std::string &txt)
{
   if (txt == "")
      return "";
   else
      return txt + " ";
}

////////////////////////////////////////////////////////////////////////////////

void RScanner::DeclInfo(clang::Decl* D) const
{
   std::string location = GetLocation(D);
   std::string kind = D->getDeclKindName();
   std::string name = GetName(D);
   ShowInfo("Scan: " + kind + " declaration " + name, location);
}

////////////////////////////////////////////////////////////////////////////////
/// unknown - this kind of declaration was not known to programmer

void RScanner::UnknownDecl(clang::Decl* D, const std::string &txt) const
{
   std::string location = GetLocation(D);
   std::string kind = D->getDeclKindName();
   std::string name = GetName(D);
   ShowWarning("Unknown " + AddSpace(txt) + kind + " declaration " + name, location);
}

////////////////////////////////////////////////////////////////////////////////
/// unexpected - this kind of declaration is unexpected (in concrete place)

void RScanner::UnexpectedDecl(clang::Decl* D, const std::string &txt) const
{
   std::string location = GetLocation(D);
   std::string kind = D->getDeclKindName();
   std::string name = GetName(D);
   ShowWarning("Unexpected " + kind + " declaration " + name, location);
}

////////////////////////////////////////////////////////////////////////////////
/// unsupported - this kind of declaration is probably not used (in current version of C++)

void RScanner::UnsupportedDecl(clang::Decl* D, const std::string &txt) const
{
   std::string location = GetLocation(D);
   std::string kind = D->getDeclKindName();
   std::string name = GetName(D);
   ShowWarning("Unsupported " + AddSpace(txt) + kind + " declaration " + name, location);
}

////////////////////////////////////////////////////////////////////////////////
/// unimportant - this kind of declaration is not stored into reflex

void RScanner::UnimportantDecl(clang::Decl* D, const std::string &txt) const
{
}

////////////////////////////////////////////////////////////////////////////////
/// information about item, that should be implemented

void RScanner::UnimplementedDecl(clang::Decl* D, const std::string &txt)
{
   clang::Decl::Kind k = D->getKind();

   bool show = true;
   if (k <= fgDeclLast) {
      if (fDeclTable [k])
         show = false; // already displayed
      else
         fDeclTable [k] = true;
   }

   if (show)
   {
      std::string location = GetLocation(D);
      std::string kind = D->getDeclKindName();
      std::string name = GetName(D);
      std::string msg = "Unimplemented ";
      if (txt == "") {
         msg +=  "declaration";
      } else {
         msg += txt;
      }
      msg += ": ";
      msg += kind;
      msg += " ";
      msg += name;
      ShowWarning(msg,location);
   }
}

////////////////////////////////////////////////////////////////////////////////

void RScanner::UnknownType(clang::QualType qual_type) const
{
   std::string location = GetLocation(fLastDecl);
   std::string kind = qual_type.getTypePtr()->getTypeClassName();
   ShowWarning("Unknown " + kind + " type " + qual_type.getAsString(), location);
}

////////////////////////////////////////////////////////////////////////////////

void RScanner::UnsupportedType(clang::QualType qual_type) const
{
   std::string location = GetLocation(fLastDecl);
   std::string kind = qual_type.getTypePtr()->getTypeClassName();
   ShowWarning("Unsupported " + kind + " type " + qual_type.getAsString(), location);
}

////////////////////////////////////////////////////////////////////////////////

std::string RScanner::GetEnumName(clang::EnumDecl* D) const
{
   std::string enum_name = D->getQualifiedNameAsString();

   if (! D->getDeclName ()) {
      if (fgAnonymousEnumMap.find (D) != fgAnonymousEnumMap.end())
      {
         // already encountered anonymous enumeration type
         enum_name = fgAnonymousEnumMap [D];
      }
      else
      {
         fgAnonymousEnumCounter ++;
         enum_name = "_ANONYMOUS_ENUM_" + IntToStd(fgAnonymousEnumCounter) + "_";  // !?
         fgAnonymousEnumMap [D] = enum_name;
         // ShowInfo ("anonymous enum " + enum_name, GetLocation (D));
      }
   }

   return enum_name;
}

////////////////////////////////////////////////////////////////////////////////

std::string RScanner::ExprToStr(clang::Expr* expr) const
{
   clang::LangOptions lang_opts;
   clang::PrintingPolicy print_opts(lang_opts); // !?

   std::string text = "";
   llvm::raw_string_ostream stream(text);

   expr->printPretty(stream, nullptr, print_opts);

   return stream.str();
}

////////////////////////////////////////////////////////////////////////////////

std::string RScanner::ConvTemplateName(clang::TemplateName& N) const
{
   clang::LangOptions lang_opts;
   clang::PrintingPolicy print_opts(lang_opts);  // !?

   std::string text = "";
   llvm::raw_string_ostream stream(text);

   N.print(stream, print_opts);

   return stream.str();
}

////////////////////////////////////////////////////////////////////////////////

std::string RScanner::FuncParameters(clang::FunctionDecl* D) const
{
   std::string result = "";

   for (clang::FunctionDecl::param_iterator I = D->param_begin(), E = D->param_end(); I != E; ++I) {
      clang::ParmVarDecl* P = *I;

      if (result != "")
         result += ";";  // semicolon, not comma, important

      std::string type = P->getType().getAsString();
      std::string name = P->getNameAsString();

      result += type + " " + name;

      // NO if (P->hasDefaultArg ()) // check hasUnparsedDefaultArg () and hasUninstantiatedDefaultArg ()
      if (P->getInit()) {
         std::string init_value = ExprToStr(P->getDefaultArg());
         result += "=" + init_value;
      }
   }

   return result;
}

////////////////////////////////////////////////////////////////////////////////

std::string RScanner::FuncParameterList(clang::FunctionDecl* D) const
{
   std::string result = "";

   for (clang::FunctionDecl::param_iterator I = D->param_begin(), E = D->param_end(); I != E; ++I) {
      clang::ParmVarDecl* P = *I;

      if (result != "")
         result += ",";

      std::string type = P->getType().getAsString();
      result += type;
   }

   return "(" + result + ")";
}

////////////////////////////////////////////////////////////////////////////////
/// This method visits a namespace node

bool RScanner::VisitNamespaceDecl(clang::NamespaceDecl* N)
{
   // We don't need to visit this while creating the big PCM
   if (fScanType == EScanType::kOnePCM)
      return true;

   if (!shouldVisitDecl(N))
      return true;

   // in case it is implicit we don't create a builder
   // [Note: Can N be nullptr?, is so 'ShouldVisitDecl' should test or we should test sooner]
   if((N && N->isImplicit()) || !N){
      return true;
   }

   bool ret = true;

   const ClassSelectionRule *selected = fSelectionRules.IsDeclSelected(N);
   if (selected) {

      clang::DeclContext* primary_ctxt = N->getPrimaryContext();
      clang::NamespaceDecl* primary = llvm::dyn_cast<clang::NamespaceDecl>(primary_ctxt);

      RPredicateIsSameNamespace pred(primary);
      if ( find_if(fSelectedNamespaces.begin(),fSelectedNamespaces.end(),pred) == fSelectedNamespaces.end() ) {
         // The namespace is not already registered.

         if (fVerboseLevel > 0) {
            std::string qual_name;
            GetDeclQualName(N,qual_name);
            //      std::cout<<"\tSelected namespace -> " << qual_name << " ptr " << (void*)N <<   " decl ctxt " << (void*)N->getPrimaryContext() << " classname " <<primary->getNameAsString() << "\n";
            std::cout<<"\tSelected namespace -> " << qual_name << "\n";
         }
         fSelectedNamespaces.push_back(AnnotatedNamespaceDecl(primary,selected->GetIndex(),selected->RequestOnlyTClass()));
      }
      ret = true;
   }

   return ret;
}

////////////////////////////////////////////////////////////////////////////////

bool RScanner::VisitRecordDecl(clang::RecordDecl* D)
{
   if (!shouldVisitDecl(D))
      return true;

   // This method visits a class node
   return TreatRecordDeclOrTypedefNameDecl(D);


}

////////////////////////////////////////////////////////////////////////////////

int RScanner::AddAnnotatedRecordDecl(const ClassSelectionRule* selected,
                                      const clang::Type* req_type,
                                      const clang::RecordDecl* recordDecl,
                                      const std::string& attr_name,
                                      const clang::TypedefNameDecl* typedefNameDecl,
                                      unsigned int indexOffset)
{

   bool has_attr_name = selected->HasAttributeName();

   if (recordDecl->isUnion() &&
       0 != ROOT::TMetaUtils::GetClassVersion(recordDecl,fInterpreter)) {
      std::string normName;
      TMetaUtils::GetNormalizedName(normName,
                                    recordDecl->getASTContext().getCanonicalTagType(recordDecl),
                                    fInterpreter,
                                    fNormCtxt);
      ROOT::TMetaUtils::Error(nullptr,"Union %s has been selected for I/O. This is not supported. Interactive usage of unions is supported, as all C++ entities, without the need of dictionaries.\n",normName.c_str());
      return 1;
   }

   // clang-format off
   if (has_attr_name) {
      fSelectedClasses.emplace_back(selected->GetIndex() + indexOffset,
                                    req_type,
                                    recordDecl,
                                    attr_name.c_str(),
                                    selected->RequestStreamerInfo(),
                                    selected->RequestNoStreamer(),
                                    selected->RequestNoInputOperator(),
                                    selected->RequestOnlyTClass(),
                                    selected->RequestedVersionNumber(),
                                    selected->RequestedRNTupleSerializationMode(),
                                    selected->RequestedRNTupleSoARecord(),
                                    fInterpreter,
                                    fNormCtxt);
   } else {
      fSelectedClasses.emplace_back(selected->GetIndex() + indexOffset,
                                    recordDecl,
                                    selected->RequestStreamerInfo(),
                                    selected->RequestNoStreamer(),
                                    selected->RequestNoInputOperator(),
                                    selected->RequestOnlyTClass(),
                                    selected->RequestedVersionNumber(),
                                    selected->RequestedRNTupleSerializationMode(),
                                    selected->RequestedRNTupleSoARecord(),
                                    fInterpreter,
                                    fNormCtxt);
   }
   // clang-format on

   if (fVerboseLevel > 0) {
      std::string qual_name;
      GetDeclQualName(recordDecl,qual_name);
      std::string normName;
      TMetaUtils::GetNormalizedName(normName,
                                    recordDecl->getASTContext().getCanonicalTagType(recordDecl),
                                    fInterpreter,
                                    fNormCtxt);
      std::string typedef_qual_name;
      std::string typedefMsg;
      if (typedefNameDecl){
         GetDeclQualName(typedefNameDecl,typedef_qual_name);
         typedefMsg = "(through typedef/alias " + typedef_qual_name + ") ";
      }

      std::cout << "Selected class "
      << typedefMsg
      << "-> "
      << qual_name
      << " for ROOT: "
      << normName
      << "\n";
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////

bool RScanner::TreatRecordDeclOrTypedefNameDecl(clang::TypeDecl* typeDecl)
{
   // For every class is created a new class buider irrespectful of weather the
   // class is internal for another class declaration or not.
   // RecordDecls and TypedefDecls (or RecordDecls!) are treated.
   // We follow two different codepaths if the typeDecl is a RecordDecl or
   // a TypedefDecl. If typeDecl is a TypedefDecl, recordDecl becomes the
   // underlying RecordDecl.
   // This is done to leverage the selections rule matching in SelectionRules
   // which works basically with names.
   // At the end of the method, if the typedef name is matched, an AnnotatedRecordDecl
   // with the underlying RecordDecl is fed to the machinery.

   const clang::RecordDecl* recordDecl = clang::dyn_cast<clang::RecordDecl>(typeDecl);
   const clang::TypedefNameDecl* typedefNameDecl = clang::dyn_cast<clang::TypedefNameDecl>(typeDecl);

   // If typeDecl is not a RecordDecl, try to fetch the RecordDecl behind the TypedefDecl
   if (!recordDecl && typedefNameDecl) {
      recordDecl = ROOT::TMetaUtils::GetUnderlyingRecordDecl(typedefNameDecl->getUnderlyingType());
   }

   // If at this point recordDecl is still NULL, we have a problem
   if (!recordDecl) {
      ROOT::TMetaUtils::Warning("RScanner::TreatRecordDeclOrTypeNameDecl",
       "Could not cast typeDecl either to RecordDecl or could not get RecordDecl underneath typedef.\n");
      return true;
   }

   // Do not select unnamed records.
   if (!recordDecl->getIdentifier())
      return true;

   // Do not select dependent types.
   if (recordDecl->isDependentType())
      return true;

   if (fScanType == EScanType::kOnePCM && ROOT::TMetaUtils::IsStdClass(*recordDecl))
      return true;


   // At this point, recordDecl must be a RecordDecl pointer.

   if (fRecordDeclCallback) {
      // Pass on any declaration.   This is usually used to record dependency.
      // Since rootcint see C++ compliant header files, we can assume that
      // if a forward declaration or declaration has been inserted, the
      // classes for which we are creating a dictionary will be using
      // them either directly or indirectly.   Any false positive can be
      // resolved by removing the spurrious dependency in the (user) header
      // files.
      fRecordDeclCallback(recordDecl);
   }

   // in case it is implicit or a forward declaration, we are not interested.
   if(recordDecl->isImplicit() || !recordDecl->isCompleteDefinition()) {
      return true;
   }

   // Never select the class templates themselves.
   const clang::CXXRecordDecl *cxxdecl = llvm::dyn_cast<clang::CXXRecordDecl>(recordDecl);
   if (cxxdecl && cxxdecl->getDescribedClassTemplate ()) {
      return true;
   }

   const ClassSelectionRule *selectedFromTypedef = typedefNameDecl ? fSelectionRules.IsDeclSelected(typedefNameDecl) : nullptr;

   const ClassSelectionRule *selectedFromRecDecl = fSelectionRules.IsDeclSelected(recordDecl, false /* exclude typedef rules*/);

   const ClassSelectionRule *selected = typedefNameDecl ? selectedFromTypedef : selectedFromRecDecl;

   if (! selected) return true; // early exit. Nothing more to be done.

   // Selected through typedef but excluded with concrete classname
   bool excludedFromRecDecl = false;
   if ( selectedFromRecDecl )
      excludedFromRecDecl = selectedFromRecDecl->GetSelected() == BaseSelectionRule::kNo;

   if (selected->GetSelected() != BaseSelectionRule::kYes || excludedFromRecDecl)
      return true;

   // Save the typedef
   if (selectedFromTypedef) {
      if (!IsElementPresent(fSelectedTypedefs, typedefNameDecl))
         fSelectedTypedefs.push_back(typedefNameDecl);
      // Don't generate a dictionary for the class underlying a typedef found
      // for a file name match (eg. "defined_in")
      if (!selectedFromRecDecl && selectedFromTypedef->HasAttributeFileName())
         return true;
   }

   if (selected->IsFromTypedef()) {
      if (typedefNameDecl && !IsElementPresent(fSelectedTypedefs, typedefNameDecl))
         fSelectedTypedefs.push_back(typedefNameDecl);
      return true;
   }

   if (typedefNameDecl)
      ROOT::TMetaUtils::Info("RScanner::TreatRecordDeclOrTypedefNameDecl",
                              "Typedef is selected %s.\n", typedefNameDecl->getNameAsString().c_str());

   // For the case kNo, we could (but don't) remove the node from the pcm
   // For the case kDontCare, the rule is just a place holder and we are actually trying to exclude some of its children
   // (this is used only in the selection xml case).

   // Reject the selection of std::pair on the ground that it is trivial
   // and can easily be recreated from the AST information.
   if (recordDecl->getName() == "pair") {
      const clang::NamespaceDecl *nsDecl = llvm::dyn_cast<clang::NamespaceDecl>(recordDecl->getDeclContext());
      if (!nsDecl){
         ROOT::TMetaUtils::Error("RScanner::TreatRecordDeclOrTypedefNameDecl",
                                 "Cannot convert context of RecordDecl called pair into a namespace.\n");
         return true;
      }
      const clang::NamespaceDecl *nsCanonical = nsDecl->getCanonicalDecl();
      if (nsCanonical && nsCanonical == fInterpreter.getCI()->getSema().getStdNamespace()) {
         if (selected->HasAttributeFileName() || selected->HasAttributeFilePattern()) {
            return true;
         }
      }
   }

   // Insert in the selected classes if not already there
   // We need this check since the same class can be selected through its name or typedef
   bool rcrdDeclNotAlreadySelected = fselectedRecordDecls.insert((RecordDecl*)recordDecl->getCanonicalDecl()).second;
   if (!fFirstPass && !rcrdDeclNotAlreadySelected) {
      // Diagnose conflicting selection rules:
      auto declSelRuleMapIt = fDeclSelRuleMap.find(recordDecl->getCanonicalDecl());
      if (declSelRuleMapIt != fDeclSelRuleMap.end() &&
            declSelRuleMapIt->second != selected) {
         std::string normName;
         TMetaUtils::GetNormalizedName(normName,
                                       recordDecl->getASTContext().getCanonicalTagType(recordDecl),
                                       fInterpreter,
                                       fNormCtxt);

         auto previouslyMatchingRule = (const ClassSelectionRule*)declSelRuleMapIt->second;
         int previouslineno = previouslyMatchingRule->GetLineNumber();

         std::string cleanFileName =  llvm::sys::path::filename(selected->GetSelFileName()).str();
         auto lineno = selected->GetLineNumber();
         auto rulesAreCompatible = SelectionRulesUtils::areEqual<ClassSelectionRule>(selected, previouslyMatchingRule, true /*moduloNameOrPattern*/);
         if (!rulesAreCompatible){
            std::stringstream message;
            if (lineno > 1) message << "Selection file " << cleanFileName << ", lines "
                                    << lineno << " and " << previouslineno << ". ";
            message << "Attempt to select a class "<< normName << " with two rules which have incompatible attributes. "
                  << "The attributes such as transiency might not be correctly propagated to the typesystem of ROOT.\n";
            selected->Print(message);
            message << "Conflicting rule already matched:\n";
            previouslyMatchingRule->Print(message);
            ROOT::TMetaUtils::Warning(nullptr,"%s\n", message.str().c_str());
         }
      }
   }

   fDeclSelRuleMap[recordDecl->getCanonicalDecl()] = selected;

   if (!rcrdDeclNotAlreadySelected || fFirstPass)
      return true;

   // Before adding the decl to the selected ones, check its access.
   // We do not yet support I/O of private or protected classes.
   // See ROOT-7450.
   // Additionally, private declarations lead to uncompilable code, so just ignore (ROOT-9112).
   if (recordDecl->getAccess() == AS_private || recordDecl->getAccess() == AS_protected) {
      // Don't warn about types selected by "everything in that file".
      auto isFileSelection = selected->HasAttributeFileName() &&
                           selected->HasAttributePattern() &&
                           "*" == selected->GetAttributePattern();
      if (!isFileSelection) {
         std::string normName;
         TMetaUtils::GetNormalizedName(normName,
                                       recordDecl->getASTContext().getCanonicalTagType(recordDecl),
                                       fInterpreter,
                                       fNormCtxt);
         auto msg = "Class or struct %s was selected but its dictionary cannot be generated: "
                  "this is a private or protected class and this is not supported. No direct "
                  "I/O operation of %s instances will be possible.\n";
         ROOT::TMetaUtils::Warning(nullptr,msg,normName.c_str(),normName.c_str());
      }
      return true;
   }

   // Replace on the fly the type if the type for IO is different for example
   // in presence of unique_ptr<T> or collections thereof.
   // The following lines are very delicate: we need to preserve the special
   // ROOT opaque typedefs.
   auto req_type = selected->GetRequestedType();
   clang::QualType thisType(req_type, 0);
   std::string attr_name = selected->GetAttributeName().c_str();

   auto sc = AddAnnotatedRecordDecl(selected, req_type, recordDecl, attr_name, typedefNameDecl);
   if (sc != 0) {
      return false;
   }

   if (auto CTSD = llvm::dyn_cast<clang::ClassTemplateSpecializationDecl>(recordDecl)) {
      fDelayedAnnotatedRecordDecls.emplace_back(DelayedAnnotatedRecordDeclInfo{selected, CTSD, typedefNameDecl});
   }

   return true;
}

void RScanner::AddDelayedAnnotatedRecordDecls()
{
   for (auto &&info: fDelayedAnnotatedRecordDecls) {
      const clang::Type *thisType = info.fSelected->GetRequestedType();
      if (!thisType)
         thisType = info.fDecl->getASTContext().getCanonicalTagType(info.fDecl)->getTypePtr();
      const clang::CXXRecordDecl *recordDecl = info.fDecl;
      auto nameTypeForIO = ROOT::TMetaUtils::GetNameTypeForIO(clang::QualType(thisType, 0), fInterpreter, fNormCtxt);
      auto typeForIO = nameTypeForIO.second;
      // It could be that we have in hands a type which is not a class, e.g.
      // in presence of unique_ptr<T> we got a T with T=double.
      if (typeForIO.getTypePtr() == thisType)
         continue;
      if (auto recordDeclForIO = typeForIO->getAsCXXRecordDecl()) {
         const auto canRecordDeclForIO = recordDeclForIO->getCanonicalDecl();
         if (!fselectedRecordDecls.insert(canRecordDeclForIO).second)
            continue;
         recordDecl = canRecordDeclForIO;
         fDeclSelRuleMap[recordDecl] = info.fSelected;
         thisType = typeForIO.getTypePtr();
      }

      AddAnnotatedRecordDecl(info.fSelected, thisType, recordDecl,
                             nameTypeForIO.first, info.fTypedefNameDecl, 1000);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Visitor for every TypedefNameDecl, i.e. aliases and typedefs
/// We check three conditions before trying to match the name:
/// 1) If we are creating a big PCM
/// 2) If the underlying decl is a RecordDecl
/// 3) If the typedef is eventually contained in the std namespace

bool RScanner::VisitTypedefNameDecl(clang::TypedefNameDecl* D)
{
   if (fScanType == EScanType::kOnePCM)
      return true;

   if (!shouldVisitDecl(D))
      return true;

   const clang::DeclContext *ctx = D->getDeclContext();

   bool isInStd=false;
   if (ctx) {
      const clang::NamedDecl *parent = llvm::dyn_cast<clang::NamedDecl> (ctx);
      isInStd = parent && 0 == parent->getQualifiedNameAsString().compare(0,5,"std::");
      }

   if (ROOT::TMetaUtils::GetUnderlyingRecordDecl(D->getUnderlyingType()) &&
       !isInStd){
      TreatRecordDeclOrTypedefNameDecl(D);
   }

    return true;
}

////////////////////////////////////////////////////////////////////////////////

bool RScanner::VisitEnumDecl(clang::EnumDecl* D)
{
   if (fScanType == EScanType::kOnePCM)
      return true;

   if (!shouldVisitDecl(D))
      return true;

   if(fSelectionRules.IsDeclSelected(D) &&
      !IsElementPresent(fSelectedEnums, D)){ // Removal of duplicates.
      fSelectedEnums.push_back(D);
   }

   return true;
}

////////////////////////////////////////////////////////////////////////////////

bool RScanner::VisitVarDecl(clang::VarDecl* D)
{
   if (!D->hasGlobalStorage() ||
       fScanType == EScanType::kOnePCM)
      return true;

   if (!shouldVisitDecl(D))
      return true;

   if(fSelectionRules.IsDeclSelected(D)){
      fSelectedVariables.push_back(D);
   }

   return true;
}

////////////////////////////////////////////////////////////////////////////////
/// Nothing to be done here

bool RScanner::VisitFieldDecl(clang::FieldDecl* D)
{
   return true;

//    bool ret = true;
//
//    if(fSelectionRules.IsDeclSelected(D)){
//
//       // if (fVerboseLevel > 0) {
// //      std::string qual_name;
// //      GetDeclQualName(D,qual_name);
// //      std::cout<<"\tSelected field -> " << qual_name << "\n";
//       // }
//    }
//    else {
//    }
//
//    return ret;
}

////////////////////////////////////////////////////////////////////////////////

bool RScanner::VisitFunctionDecl(clang::FunctionDecl* D)
{
   if (fScanType == EScanType::kOnePCM)
      return true;

   if (!shouldVisitDecl(D))
      return true;

   if(clang::FunctionDecl::TemplatedKind::TK_FunctionTemplate == D->getTemplatedKind())
      return true;

   if(fSelectionRules.IsDeclSelected(D)){
      fSelectedFunctions.push_back(D);
   }

   return true;
}

////////////////////////////////////////////////////////////////////////////////
/// Whether it is enough for the AST visitors to traverse only the members of
/// DC that are already loaded (see ROOT::Internal::GetLoadedDecls()), instead
/// of all of them. This is the case for classes, enums, functions, etc. that
/// were deserialized from a C++ module, e.g. the one of a dictionary we depend
/// on. Iterating over all members of such a DeclContext with decls() makes the
/// ASTReader deserialize every one of them, which dominates the memory usage
/// and run time of rootcling.
///
/// Namespaces, linkage specs and export decls are always traversed fully. So
/// is everything if modules are disabled (e.g. with a PCH), since then the
/// ASTReader doesn't add the decls it loads to the lookup tables.

bool ROOT::Internal::CanTraverseOnlyLoadedDecls(const clang::DeclContext *DC)
{
   auto *D = llvm::dyn_cast_or_null<clang::Decl>(DC);
   return D && D->isFromASTFile() && D->getASTContext().getLangOpts().Modules &&
          !llvm::isa<clang::NamespaceDecl, clang::LinkageSpecDecl, clang::ExportDecl>(D);
}

////////////////////////////////////////////////////////////////////////////////
/// Collect the decls lexically contained in DC that are already loaded,
/// without deserializing the other members of DC. These are:
///
///  - the lexical decls of DC that are loaded already,
///  - the named decls in the lookup table of DC: when the ASTReader
///    deserializes the first declaration of a named member for any reason
///    (e.g. name lookup, or because it is needed by a template
///    instantiation), it adds it to the lookup table of its DeclContext, so
///    that it can be merged with declarations of the same entity from other
///    modules (see clang::ASTDeclReader::FindExistingResult's destructor,
///    which calls DeclContext::makeDeclVisibleInContextImpl()). A name lookup
///    might replace it by another declaration of the same entity,
///  - the first declarations of the templates and the definitions of the
///    classes found above, if they are in DC,
///  - the explicit specializations of the member templates found above, which
///    are not in the lookup table, and which RecursiveASTVisitor would only
///    visit as lexical decls of DC.
///
/// Getting the specializations loads the lazy specializations of these member
/// templates, like RecursiveASTVisitor does when it visits their
/// instantiations anyway.
///
/// The decls are sorted by source location, so they are visited in a
/// deterministic order that does not depend on the lookup table, which is a
/// hash map, and in the same order as in a full traversal.

void ROOT::Internal::GetLoadedDecls(clang::DeclContext *DC, llvm::SmallVectorImpl<clang::Decl *> &decls)
{
   llvm::SmallPtrSet<clang::Decl *, 32> seen;
   auto add = [&](clang::Decl *Child) {
      if (seen.insert(Child).second)
         decls.push_back(Child);
   };
   auto addIfInDC = [&](clang::Decl *Child) {
      if (Child && Child->getLexicalDeclContext() == DC)
         add(Child);
   };
   // RecursiveASTVisitor traverses the instantiations of a template from its
   // first declaration, and the RScanner selects a class at its definition.
   // Add these if they are in DC, so that they are visited at their position
   // in DC, like in a full traversal. This matters if a class is selected by
   // several rules, e.g. by its name and by a typedef: the first one wins.
   auto addWithFirstDeclOrDefinition = [&](clang::Decl *Child) {
      addIfInDC(Child);
      if (auto *TD = llvm::dyn_cast<clang::RedeclarableTemplateDecl>(Child))
         addIfInDC(TD->getCanonicalDecl());
      else if (auto *RD = llvm::dyn_cast<clang::CXXRecordDecl>(Child))
         addIfInDC(RD->getDefinition());
   };

   // Collect the loaded decls first: getting the first declaration or the
   // definition of a decl might deserialize declarations, which the ASTReader
   // adds to the lookup table of DC that we iterate over.
   llvm::SmallVector<clang::Decl *, 32> loaded;
   for (clang::Decl *Child : DC->noload_decls())
      loaded.push_back(Child);
   for (clang::DeclContextLookupResult lookupResult : DC->noload_lookups(/*PreserveInternalState=*/true)) {
      for (clang::NamedDecl *ND : lookupResult)
         loaded.push_back(ND);
   }
   for (clang::Decl *Child : loaded)
      addWithFirstDeclOrDefinition(Child);

   auto addExplicitSpecializations = [&](auto *TD) {
      for (auto *Spec : TD->specializations()) {
         using SpecDecl = std::remove_pointer_t<decltype(Spec)>;
         for (auto *RD : Spec->redecls()) {
            const clang::TemplateSpecializationKind kind = llvm::cast<SpecDecl>(RD)->getSpecializationKind();
            if (kind != clang::TSK_ImplicitInstantiation && kind != clang::TSK_Undeclared &&
                RD->getLexicalDeclContext() == DC)
               add(RD);
         }
      }
   };
   for (std::size_t i = 0, n = decls.size(); i < n; ++i) {
      if (auto *CTD = llvm::dyn_cast<clang::ClassTemplateDecl>(decls[i]))
         addExplicitSpecializations(CTD);
      else if (auto *VTD = llvm::dyn_cast<clang::VarTemplateDecl>(decls[i]))
         addExplicitSpecializations(VTD);
   }

   // Sort by source location. Use the location of the name, not the beginning
   // of the decl: the parser adds a class to DC when it sees its name, and a
   // typedef etc. after its declarator, so the class in
   // `typedef struct A {} B;` comes before the typedef in the list of lexical
   // decls of DC. Decls without a valid location go first. Decls at the same
   // location are sorted by their ID in the AST file, except for the ones that
   // were not deserialized, which go first and keep their order in the list of
   // lexical decls.
   // A member class of a class template specialization that was declared before
   // its definition gets the location of the definition in the template when
   // it is instantiated. Use the location of the declaration it was instantiated
   // from, which is where a full traversal visits it.
   auto getLoc = [](clang::Decl *D) {
      if (auto *RD = llvm::dyn_cast<clang::CXXRecordDecl>(D)) {
         if (clang::CXXRecordDecl *pattern = RD->getInstantiatedFromMemberClass())
            return pattern->getLocation();
      }
      return D->getLocation();
   };
   const clang::SourceManager &SM = DC->getParentASTContext().getSourceManager();
   std::stable_sort(decls.begin(), decls.end(), [&SM, &getLoc](clang::Decl *A, clang::Decl *B) {
      const clang::SourceLocation locA = getLoc(A);
      const clang::SourceLocation locB = getLoc(B);
      if (locA.isInvalid() != locB.isInvalid())
         return locA.isInvalid();
      if (locA.isValid() && locA != locB) {
         if (SM.isBeforeInTranslationUnit(locA, locB))
            return true;
         if (SM.isBeforeInTranslationUnit(locB, locA))
            return false;
      }
      if (A->isFromASTFile() != B->isFromASTFile())
         return !A->isFromASTFile();
      return A->isFromASTFile() && A->getGlobalID() < B->getGlobalID();
   });
}

////////////////////////////////////////////////////////////////////////////////
/// The decls that the RScanner visits in a class etc. from an AST file are
/// limited to the ones that are already loaded (see
/// ROOT::Internal::GetLoadedDecls()). This is enough for the class selection
/// rules that select by name: SelectionRules::SearchNames() has looked up
/// these classes before the scan, which loaded them. Other selection rules
/// might select members that are not loaded, so collect the beginnings of the
/// qualified names they can match, up to the first wildcard or template
/// argument list. These are:
///
///  - class selection rules with a pattern, or with a name that was not found
///    by SearchNames(),
///  - enum and variable selection rules (e.g. for static data members).
///
/// Rules on the file name can select any decl, so with these all members are
/// traversed. Function selection rules don't select member functions, and the
/// only functions lexically declared in classes are friend functions, which
/// are not considered here. Class and enum rules that don't select anything
/// are not relevant either, but variable rules are, because
/// SelectionRules::AreAllSelectionRulesUsed() warns about the ones that don't
/// match any variable.

void RScanner::CollectMemberRulePrefixes()
{
   fMemberRulePrefixes.clear();
   fTraverseAllMembers = fSelectionRules.GetHasFileNameRule();

   auto addPrefix = [this](const std::string &nameOrPattern) {
      std::string prefix;
      for (char c : nameOrPattern) {
         if (c == '*' || c == '\\' || c == '<')
            break;
         if (!std::isspace(static_cast<unsigned char>(c)))
            prefix += c;
      }
      if (prefix.compare(0, 2, "::") == 0)
         prefix.erase(0, 2);
      if (prefix.empty())
         fTraverseAllMembers = true;
      else
         fMemberRulePrefixes.push_back(prefix);
   };

   auto processRule = [&](const BaseSelectionRule &rule, bool nameIsLoaded, bool needsMatchIfNotSelecting) {
      if (rule.GetSelected() == BaseSelectionRule::kNo && !needsMatchIfNotSelecting)
         return;
      if (rule.HasAttributeFileName() || rule.HasAttributeFilePattern())
         fTraverseAllMembers = true;
      else if (rule.HasAttributePattern())
         addPrefix(rule.GetAttributePattern());
      else if (rule.HasAttributeName() && !nameIsLoaded)
         addPrefix(rule.GetAttributeName());
   };

   for (const ClassSelectionRule &rule : fSelectionRules.GetClassSelectionRules()) {
      const clang::CXXRecordDecl *decl = rule.GetCXXRecordDecl();
      processRule(rule, decl && decl != reinterpret_cast<const clang::CXXRecordDecl *>(-1),
                  /*needsMatchIfNotSelecting=*/false);
   }
   for (const BaseSelectionRule &rule : fSelectionRules.GetEnumSelectionRules())
      processRule(rule, false, /*needsMatchIfNotSelecting=*/false);
   for (const BaseSelectionRule &rule : fSelectionRules.GetVariableSelectionRules())
      processRule(rule, false, /*needsMatchIfNotSelecting=*/true);
}

////////////////////////////////////////////////////////////////////////////////
/// Whether one of the selection rules might select a member of D that is not
/// loaded yet. See CollectMemberRulePrefixes().

bool RScanner::MayHaveSelectedUnloadedMembers(const clang::Decl &D) const
{
   if (fTraverseAllMembers)
      return true;
   if (fMemberRulePrefixes.empty())
      return false;

   // Build the qualified name of D up to the first template or unnamed scope,
   // both with and without inline namespaces, since the selection rules might
   // be matched against either.
   llvm::SmallVector<const clang::Decl *, 8> scopes;
   for (const clang::Decl *S = &D; S && !llvm::isa<clang::TranslationUnitDecl>(S);
        S = llvm::dyn_cast<clang::Decl>(S->getDeclContext()))
      scopes.push_back(S);

   std::string name;
   std::string nameNoInline;
   bool complete = true;
   for (auto it = scopes.rbegin(); it != scopes.rend(); ++it) {
      const clang::Decl *S = *it;
      if (llvm::isa<clang::LinkageSpecDecl, clang::ExportDecl>(S))
         continue;
      std::string component;
      bool isInline = false;
      if (auto *NS = llvm::dyn_cast<clang::NamespaceDecl>(S)) {
         component = NS->isAnonymousNamespace() ? "(anonymous namespace)" : NS->getNameAsString();
         isInline = NS->isInline();
      } else if (auto *RD = llvm::dyn_cast<clang::CXXRecordDecl>(S)) {
         component = RD->getNameAsString();
         if (component.empty() || llvm::isa<clang::ClassTemplateSpecializationDecl>(RD) ||
             RD->getDescribedClassTemplate())
            complete = false;
      } else if (auto *ND = llvm::dyn_cast<clang::NamedDecl>(S)) {
         component = ND->getNameAsString();
         if (!llvm::isa<clang::TagDecl>(ND) || component.empty())
            complete = false;
      } else {
         complete = false;
      }
      name += component;
      if (!isInline)
         nameNoInline += component;
      if (!complete)
         break;
      name += "::";
      if (!isInline)
         nameNoInline += "::";
   }

   // If the name is complete, all members have names starting with name, and
   // otherwise with an unknown continuation of it.
   auto mayMatch = [](const std::string &prefix, const std::string &scopeName) {
      const std::size_t n = std::min(prefix.size(), scopeName.size());
      return prefix.compare(0, n, scopeName, 0, n) == 0;
   };
   for (const std::string &prefix : fMemberRulePrefixes) {
      if (mayMatch(prefix, name) || mayMatch(prefix, nameNoInline))
         return true;
   }
   return false;
}

////////////////////////////////////////////////////////////////////////////////
/// Whether D might not be traversed because it is lexically in a class etc.
/// from a C++ module of which only the loaded members are traversed.

bool RScanner::IsInPartlyTraversedDeclContext(const clang::Decl &D) const
{
   for (const clang::DeclContext *DC = D.getLexicalDeclContext(); DC; DC = DC->getLexicalParent()) {
      if (ROOT::Internal::CanTraverseOnlyLoadedDecls(DC) &&
          !MayHaveSelectedUnloadedMembers(*llvm::cast<clang::Decl>(DC)))
         return true;
   }
   return false;
}

////////////////////////////////////////////////////////////////////////////////
/// RecursiveASTVisitor only traverses the instantiations of a template when it
/// visits the canonical, i.e. first, declaration of the template. This might
/// not be visited if it is in a class etc. from a C++ module of which only the
/// loaded members are traversed, e.g. for a member template that is defined
/// outside of the class, or for a template that is first declared as a friend.
/// So when visiting another declaration of such a template, traverse the
/// instantiations of its canonical declaration, unless this has been done.

template <class TemplateDecl>
bool RScanner::TraverseInstantiationsFromCanonicalDecl(TemplateDecl *D)
{
   TemplateDecl *canonical = D->getCanonicalDecl();
   if (!canonical->isFromASTFile())
      return true;
   if (D == canonical) {
      // RecursiveASTVisitor has traversed its instantiations.
      fTemplatesWithTraversedInstantiations.insert(canonical);
      return true;
   }
   if (!IsInPartlyTraversedDeclContext(*canonical) || !fTemplatesWithTraversedInstantiations.insert(canonical).second)
      return true;
   return TraverseTemplateInstantiations(canonical);
}

bool RScanner::TraverseClassTemplateDecl(clang::ClassTemplateDecl *D)
{
   return RecursiveASTVisitor<RScanner>::TraverseClassTemplateDecl(D) && TraverseInstantiationsFromCanonicalDecl(D);
}

bool RScanner::TraverseVarTemplateDecl(clang::VarTemplateDecl *D)
{
   return RecursiveASTVisitor<RScanner>::TraverseVarTemplateDecl(D) && TraverseInstantiationsFromCanonicalDecl(D);
}

bool RScanner::TraverseFunctionTemplateDecl(clang::FunctionTemplateDecl *D)
{
   return RecursiveASTVisitor<RScanner>::TraverseFunctionTemplateDecl(D) && TraverseInstantiationsFromCanonicalDecl(D);
}

////////////////////////////////////////////////////////////////////////////////

bool RScanner::TraverseDeclContextHelper(DeclContext *DC)
{
   bool ret = true;

   if (!DC)
      return true;

   clang::Decl* D = dyn_cast<clang::Decl>(DC);
   // skip implicit decls
   if (D && D->isImplicit()){
      return true;
   }

   // For classes etc. from C++ modules, only traverse the members that are
   // already loaded, unless the selection rules might select others.
   if (ROOT::Internal::CanTraverseOnlyLoadedDecls(DC) && !MayHaveSelectedUnloadedMembers(*D)) {
      llvm::SmallVector<clang::Decl *, 32> loadedDecls;
      ROOT::Internal::GetLoadedDecls(DC, loadedDecls);
      for (clang::Decl *Child : loadedDecls) {
         if (!TraverseDecl(Child))
            return false;
      }
      return true;
   }

   if (fScanType == EScanType::kOnePCM){
      const clang::NamespaceDecl *parent = llvm::dyn_cast<clang::NamespaceDecl> (DC);
      if (parent && 0 == parent->getQualifiedNameAsString().compare(0,5,"std::"))
         return true;
      }

   for (DeclContext::decl_iterator Child = DC->decls_begin(), ChildEnd = DC->decls_end();
        ret && (Child != ChildEnd); ++Child) {
      ret=TraverseDecl(*Child);
   }

   return ret;

}

////////////////////////////////////////////////////////////////////////////////

bool RScanner::GetDeclName(clang::Decl* D, std::string& name) const
{
   clang::NamedDecl* N = dyn_cast<clang::NamedDecl> (D);

   if (N) {
      name = N->getNameAsString();
      return true;
   }
   else {
      name = "UNNAMED";
      return false;
   }
}

////////////////////////////////////////////////////////////////////////////////

bool RScanner::GetDeclQualName(const clang::Decl* D, std::string& qual_name)
{
   auto N = dyn_cast<const clang::NamedDecl> (D);

   if (N) {
      llvm::raw_string_ostream stream(qual_name);
      N->getNameForDiagnostic(stream,D->getASTContext().getPrintingPolicy(),true); // qual_name = N->getQualifiedNameAsString();
      return true;
   }
   else {
      return false;
   }
}

////////////////////////////////////////////////////////////////////////////////

bool RScanner::GetFunctionPrototype(clang::Decl* D, std::string& prototype) const {
   if (!D) {
      return false;
   }

   clang::FunctionDecl* F = dyn_cast<clang::FunctionDecl> (D);

   if (F) {

      prototype = "";
      for (clang::FunctionDecl::param_iterator I = F->param_begin(), E = F->param_end(); I != E; ++I) {
         clang::ParmVarDecl* P = *I;

         if (prototype != "")
            prototype += ",";

         //std::string type = P->getType().getAsString();
         std::string type = P->getType().getAsString();
         if (type.at(type.length()-1) == '*') {
            type.at(type.length()-2) = '*';
            type.erase(type.length()-1);
         }
         prototype += type;
      }

      prototype = "(" + prototype + ")";
      return true;
   }
   else {
      ShowWarning("can't convert Decl to FunctionDecl","");
      return false;
   }
}

////////////////////////////////////////////////////////////////////////////////

void RScanner::Scan(const clang::ASTContext &C)
{
   fSourceManager = &C.getSourceManager();

//    if (fVerboseLevel >= 3) fSelectionRules.PrintSelectionRules();

   if (fVerboseLevel > 0 && fSelectionRules.GetHasFileNameRule())  {
      std::cout<<"File name detected"<<std::endl;
   }

   cling::Interpreter::PushTransactionRAII RAII(const_cast<cling::Interpreter *>(&fInterpreter));

   if (fScanType == EScanType::kTwoPasses)
      TraverseDecl(C.getTranslationUnitDecl());

   fFirstPass=false;
   fselectedRecordDecls.clear();
   fSelectedEnums.clear();
   fSelectedTypedefs.clear();
   fSelectedVariables.clear();
   fSelectedFunctions.clear();
   fTemplatesWithTraversedInstantiations.clear();
   TraverseDecl(C.getTranslationUnitDecl());

   // The RecursiveASTVisitor uses range-based for; we must not modify the AST
   // during iteration / visitation. Instead, buffer the lookups that could
   // potentially create new template specializations, and handle them here:
   AddDelayedAnnotatedRecordDecls();
}


////////////////////////////////////////////////////////////////////////////////
/// Set the callback to the RecordDecl and return the previous one.

RScanner::DeclCallback RScanner::SetRecordDeclCallback(RScanner::DeclCallback callback)
{
   DeclCallback old = fRecordDeclCallback;
   fRecordDeclCallback = callback;
   return old;
}
