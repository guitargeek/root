ROOTTEST_ADD_TESTDIRS()

# Issue #5174
ROOTTEST_ADD_TEST(noClassSelected
                  COMMAND ${ROOT_rootcling_CMD} -f root_5174.Dict.cc ${CMAKE_CURRENT_SOURCE_DIR}/root_5174.h
                  PASSRC 1
                  ERRREF noClassSelected.eref)

# issue #13543
ROOTTEST_GENERATE_DICTIONARY(xmlusedByRootclingDict
                             Classy.h
                             LINKDEF xmlusedByRootcling_sel.xml
                             FIXTURES_SETUP root-meta-roocling-xmlusedByRootclingDict-fixture)

ROOTTEST_ADD_TEST(xmlusedByRootcling
                  MACRO xmlusedByRootclingDict.C
                  FIXTURES_REQUIRED root-meta-roocling-xmlusedByRootclingDict-fixture)

# Nested class templates of classes from other C++ modules
ROOTTEST_GENERATE_DICTIONARY(nestedTemplateOfModuleClassDict
                             nestedTemplateOfModuleClass.h
                             LINKDEF nestedTemplateOfModuleClassLinkDef.h
                             FIXTURES_SETUP root-meta-rootcling-nestedTemplateOfModuleClassDict-fixture)

ROOTTEST_ADD_TEST(nestedTemplateOfModuleClass
                  MACRO nestedTemplateOfModuleClass.C
                  FIXTURES_REQUIRED root-meta-rootcling-nestedTemplateOfModuleClassDict-fixture)

# Selection rules that select members of classes from another C++ module
ROOTTEST_GENERATE_DICTIONARY(membersOfModuleClassDict
                             membersOfModuleClass.h
                             LINKDEF membersOfModuleClassLinkDef.h
                             FIXTURES_SETUP root-meta-rootcling-membersOfModuleClassDict-fixture)

ROOTTEST_GENERATE_DICTIONARY(membersOfModuleClassPatternDict
                             membersOfModuleClassPattern.h
                             LINKDEF membersOfModuleClassPattern_selection.xml
                             FIXTURES_REQUIRED root-meta-rootcling-membersOfModuleClassDict-fixture
                             FIXTURES_SETUP root-meta-rootcling-membersOfModuleClassPatternDict-fixture)

ROOTTEST_ADD_TEST(membersOfModuleClassPattern
                  MACRO membersOfModuleClass.C
                  MACROARG "\"membersOfModuleClassPatternDict\", \"MembersOfModuleClass::Outer::Inner1,MembersOfModuleClass::Outer::Inner2\""
                  FIXTURES_REQUIRED root-meta-rootcling-membersOfModuleClassPatternDict-fixture)

ROOTTEST_GENERATE_DICTIONARY(membersOfModuleClassRedeclDict
                             membersOfModuleClassRedecl.h
                             LINKDEF membersOfModuleClassRedeclLinkDef.h
                             OPTIONS -failOnWarnings
                             FIXTURES_REQUIRED root-meta-rootcling-membersOfModuleClassDict-fixture
                             FIXTURES_SETUP root-meta-rootcling-membersOfModuleClassRedeclDict-fixture)

ROOTTEST_ADD_TEST(membersOfModuleClassRedecl
                  MACRO membersOfModuleClass.C
                  MACROARG "\"membersOfModuleClassRedeclDict\", \"MembersOfModuleClass::OuterWithRedecls::DefinedOutside<float>,MembersOfModuleClass::OuterWithRedecls::Redeclared<int>,MembersOfModuleClass::FriendTmplFirstDeclaredAsFriend<int>,MembersOfModuleClass::TmplWithForwardDecl<int>::Member,MembersOfModuleClass::BaseWithForwardDecl::Inner\", \"MembersOfModuleClass::OuterWithTypedef::Tmpl<int>,MembersOfModuleClass::OuterWithTypedefStruct::Inner\""
                  FIXTURES_REQUIRED root-meta-rootcling-membersOfModuleClassRedeclDict-fixture)

# GCC doesn't accept explicit specializations in class scope.
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  ROOTTEST_GENERATE_DICTIONARY(membersOfModuleClassSpecDict
                               membersOfModuleClassSpec.h
                               LINKDEF membersOfModuleClassSpecLinkDef.h
                               FIXTURES_REQUIRED root-meta-rootcling-membersOfModuleClassDict-fixture
                               FIXTURES_SETUP root-meta-rootcling-membersOfModuleClassSpecDict-fixture)

  ROOTTEST_ADD_TEST(membersOfModuleClassSpec
                    MACRO membersOfModuleClass.C
                    MACROARG "\"membersOfModuleClassSpecDict\", \"MembersOfModuleClass::OuterWithSpecialization::Tmpl<char>\""
                    FIXTURES_REQUIRED root-meta-rootcling-membersOfModuleClassSpecDict-fixture)
endif()

# ROOT-9335
ROOTTEST_ADD_TEST(classDelCtor
                  MACRO classDelCtor.C+)

# ROOT-7765
ROOTTEST_GENERATE_DICTIONARY(includeInLinkdefDict
                             LINKDEF classInLinkdef.h
                             FIXTURES_SETUP root-meta-roocling-includeInLinkdefDict-fixture)

ROOTTEST_ADD_TEST(includeInLinkdef
                  MACRO includeInLinkdef.C
                  COPY_TO_BUILDDIR headerForLdef.h
                  FIXTURES_REQUIRED root-meta-roocling-includeInLinkdefDict-fixture)

# ROOT-8492
ROOTTEST_ADD_TEST(selectUnion
                  COMMAND ${ROOT_rootcling_CMD} -f NestedClassOff.Dict.cc ${CMAKE_CURRENT_SOURCE_DIR}/NestedClassOff.h ${CMAKE_CURRENT_SOURCE_DIR}/NestedClassOff.LinkDef.h
                  OUTREF selectUnion.ref
                  PASSRC 1)

# ROOT-10798
ROOTTEST_ADD_TEST(ROOT10798
                  COMMAND ${ROOT_rootcling_CMD} -f ROOT10798Dict.cxx ${CMAKE_CURRENT_SOURCE_DIR}/ROOT10798LinkDef.h)

# Issue #18833
ROOTTEST_GENERATE_DICTIONARY(streamerInfoStdFunctionDict
                             streamerInfoStdFunction.h
                             LINKDEF streamerInfoStdFunction.xml
                             FIXTURES_SETUP root-meta-roocling-streamerInfoStdFunctionDict-fixture)

ROOTTEST_ADD_TEST(streamerInfoStdFunction
                  MACRO streamerInfoStdFunction.C
                  FIXTURES_REQUIRED root-meta-roocling-streamerInfoStdFunctionDict-fixture)

# Issue #18982
ROOTTEST_ADD_TEST(selectTemplateInvalidArg
                  COMMAND ${ROOT_rootcling_CMD} -f selectTemplateInvalidArg.Dict.cc ${CMAKE_CURRENT_SOURCE_DIR}/selectTemplateInvalidArg.h ${CMAKE_CURRENT_SOURCE_DIR}/selectTemplateInvalidArg.LinkDef.h
                  ERRREF selectTemplateInvalidArg.ref)

# Test reproducing a failure seen in CMS builds with LLVM 20
# ROOT get confused normalizing variadic template arguments
ROOTTEST_GENERATE_DICTIONARY(VariadicTemplateConflict variadicTemplate.h LINKDEF variadicTemplate.xml)
