# The module map has no link declaration for the dictionary library, so ROOT
# could never autoload it. rootcling must warn about it (see #12495).
ROOTTEST_ADD_TEST(cxxmodules-missing-link-warning
                  COMMAND ${ROOT_rootcling_CMD} -f ${CMAKE_CURRENT_BINARY_DIR}/out.cxx -cxxmodule
                          -s ${CMAKE_CURRENT_BINARY_DIR}/libout.so nolink.h
                  PRECMD ${CMAKE_COMMAND} -E remove ${CMAKE_CURRENT_BINARY_DIR}/out.pcm ${CMAKE_CURRENT_BINARY_DIR}/libout_rdict.pcm
                  WORKING_DIR ${CMAKE_CURRENT_SOURCE_DIR}
                  PASSREGEX "The module 'out' has no link declaration.*e.g. link \"libout.so\"")
