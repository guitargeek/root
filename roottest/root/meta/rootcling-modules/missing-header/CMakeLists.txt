# The module map refers to a header that doesn't exist, so the module can't be
# written. rootcling must report the missing header (see #13689).
ROOTTEST_ADD_TEST(cxxmodules-missing-header-error
                  COMMAND ${ROOT_rootcling_CMD} -f ${CMAKE_CURRENT_BINARY_DIR}/out.cxx -cxxmodule present.h
                  PRECMD ${CMAKE_COMMAND} -E remove ${CMAKE_CURRENT_BINARY_DIR}/out.pcm
                  WORKING_DIR ${CMAKE_CURRENT_SOURCE_DIR}
                  PASSREGEX "Couldn't find header 'missing.h' of module 'out'")
