# The module map refers to a header relative to an include path instead of
# relative to the directory of the module map, so clang doesn't find it.
# rootcling must explain why and fail (see #13689).
ROOTTEST_ADD_TEST(cxxmodules-relative-header-error
                  COMMAND ${ROOT_rootcling_CMD} -f ${CMAKE_CURRENT_BINARY_DIR}/out.cxx -cxxmodule
                          -s ${CMAKE_CURRENT_BINARY_DIR}/libout.so -Iinc sub/relhdr.h
                  PRECMD ${CMAKE_COMMAND} -E remove ${CMAKE_CURRENT_BINARY_DIR}/out.pcm ${CMAKE_CURRENT_BINARY_DIR}/libout_rdict.pcm
                  WORKING_DIR ${CMAKE_CURRENT_SOURCE_DIR}
                  PASSREGEX "Couldn't find header 'sub/relhdr.h' of module 'out'.*relative to the directory of the modulemap")
