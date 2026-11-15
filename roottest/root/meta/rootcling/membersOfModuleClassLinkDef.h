#ifdef __CLING__

#pragma link C++ class MembersOfModuleClass::Outer+;
#pragma link C++ class MembersOfModuleClass::OuterWithSpecialization+;
#pragma link C++ class MembersOfModuleClass::OuterWithRedecls+;
#pragma link C++ class MembersOfModuleClass::UsesTmplWithForwardDecl+;
#pragma link C++ class MembersOfModuleClass::DerivedWithUsing+;
#pragma link C++ class MembersOfModuleClass::UsesOuterWithTypedef+;

#endif
