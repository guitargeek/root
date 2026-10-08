# A header of a submodule doesn't exist. rootcling must report it once, for the
# submodule that declares it, even though clang also attaches it to the
# top-level module. The header of a submodule whose requirements are not met
# is not needed by clang, so it must not be reported (see #13689).
ROOTTEST_ADD_TEST(cxxmodules-submodule-headers-error
                  COMMAND ${ROOT_rootcling_CMD} -f ${CMAKE_CURRENT_BINARY_DIR}/out.cxx -cxxmodule
                          -s ${CMAKE_CURRENT_BINARY_DIR}/libout.so present.h
                  PRECMD ${CMAKE_COMMAND} -E remove ${CMAKE_CURRENT_BINARY_DIR}/out.pcm
                  WORKING_DIR ${CMAKE_CURRENT_SOURCE_DIR}
                  PASSREGEX "Couldn't find header 'missing.h' of module 'out.missing'"
                  FAILREGEX "of module 'out' |unavailable.h")
