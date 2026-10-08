# The modules dep1 and dep2 define N::S differently. Clang only reports this
# when it reads both modules while generating the module of the dictionary
# that uses them, and rootcling must then fail (see #12788).
foreach(dep dep1 dep2)
  ROOTTEST_ADD_TEST(cxxmodules-conflicting-dependencies-${dep}
                    COMMAND ${ROOT_rootcling_CMD} -f ${CMAKE_CURRENT_BINARY_DIR}/${dep}.cxx -cxxmodule
                            -s ${CMAKE_CURRENT_BINARY_DIR}/lib${dep}.so ${dep}.h
                    PRECMD ${CMAKE_COMMAND} -E remove ${CMAKE_CURRENT_BINARY_DIR}/${dep}.pcm
                    WORKING_DIR ${CMAKE_CURRENT_SOURCE_DIR}
                    FIXTURES_SETUP root-meta-rootcling-modules-conflicting-${dep}-fixture)
endforeach()

ROOTTEST_ADD_TEST(cxxmodules-conflicting-dependencies-error
                  COMMAND ${ROOT_rootcling_CMD} -f ${CMAKE_CURRENT_BINARY_DIR}/user.cxx -cxxmodule
                          -s ${CMAKE_CURRENT_BINARY_DIR}/libuser.so -m dep1 -m dep2 user.h
                  PRECMD ${CMAKE_COMMAND} -E remove ${CMAKE_CURRENT_BINARY_DIR}/user.pcm
                  WORKING_DIR ${CMAKE_CURRENT_SOURCE_DIR}
                  PASSRC 1
                  FIXTURES_REQUIRED root-meta-rootcling-modules-conflicting-dep1-fixture
                                    root-meta-rootcling-modules-conflicting-dep2-fixture)
