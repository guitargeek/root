# The header of the dictionary is passed "as included", relative to an include
# directory, like in projects that use ROOT. The modulemap that is generated in
# the build directory must refer to the header such that clang finds it,
# otherwise rootcling fails (see #13689).
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/inc)
ROOT_GENERATE_DICTIONARY(G__IncludeDirModule sub/IncludeDirModule.h MODULE IncludeDirModule LINKDEF LinkDef.h)
ROOTTEST_LINKER_LIBRARY(IncludeDirModule G__IncludeDirModule.cxx TEST LIBRARIES ROOT::Core)

add_test(NAME roottest-root-meta-rootcling-modules-include-dir-build
         COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR}
                                  ${build_config}
                                  --target IncludeDirModule${fast}
                                  -- ${always-make})
set_property(TEST roottest-root-meta-rootcling-modules-include-dir-build
             PROPERTY FIXTURES_SETUP root-meta-rootcling-modules-include-dir-build)
if(GeneratorNeedsBuildSerialization)
  if(GeneratorNeedsResourecLock)
    set_property(TEST roottest-root-meta-rootcling-modules-include-dir-build APPEND PROPERTY RESOURCE_LOCK CMAKE_BUILD)
  endif()
  set_property(TEST roottest-root-meta-rootcling-modules-include-dir-build APPEND PROPERTY FIXTURES_REQUIRED CMAKE_BUILD_ALL)
endif()

ROOTTEST_ADD_TEST(cxxmodules-include-dir
                  MACRO execIncludeDirModule.C
                  ROOTEXE_OPTS -e "(void)gSystem->Load(\"libIncludeDirModule\")"
                  PASSREGEX "fI = 42"
                  FIXTURES_REQUIRED root-meta-rootcling-modules-include-dir-build)
