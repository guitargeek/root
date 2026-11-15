#ifdef __CLING__

#pragma link C++ class MembersOfModuleClassRedecl+;
// Instantiations of templates whose first declaration is in a class from the
// other C++ module
#pragma link C++ class MembersOfModuleClass::OuterWithRedecls::DefinedOutside<float>+;
#pragma link C++ class MembersOfModuleClass::OuterWithRedecls::Redeclared<int>+;
#pragma link C++ class MembersOfModuleClass::FriendTmplFirstDeclaredAsFriend<int>+;
// Nested classes in classes from the other C++ module that are declared before
// they are defined
#pragma link C++ class MembersOfModuleClass::TmplWithForwardDecl<int>::Member+;
#pragma link C++ class MembersOfModuleClass::DerivedWithUsing::Inner+;
// Classes selected both by name and by a typedef. The rule for the name wins,
// so they get the automatic streamer.
#pragma link C++ class MembersOfModuleClass::OuterWithTypedef::Tmpl<int>+;
#pragma link C++ class MembersOfModuleClass::OuterWithTypedef::TmplInt;
#pragma link C++ class MembersOfModuleClass::OuterWithTypedefStruct::Inner+;
#pragma link C++ class MembersOfModuleClass::OuterWithTypedefStruct::Alias;
// Excludes static data members of a class from the other C++ module. This
// dictionary is generated with -failOnWarnings, so it fails if the rule is
// unused.
#pragma link off global MembersOfModuleClass::OuterWithStatic::s*;

#endif
