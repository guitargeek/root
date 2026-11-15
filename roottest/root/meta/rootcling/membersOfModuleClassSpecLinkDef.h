#ifdef __CLING__

#pragma link C++ class MembersOfModuleClassSpec+;
// An explicit specialization of a member class template, declared in the
// class from the other C++ module
#pragma link C++ class MembersOfModuleClass::OuterWithSpecialization::Tmpl<char>+;

#endif
