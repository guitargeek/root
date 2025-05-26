# Copyright (C) 1995-2025, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

# Dictionary list
set(CPPYY_TEST_DICTS
advancedcpp
advancedcpp2
conversions
crossinheritance
datatypes
doc_helper
example01
fragile
overloads
pythonizables
std_streams
templates
)

foreach(DICT ${CPPYY_TEST_DICTS})
    set(CMAKE_ROOTTEST_NOROOTMAP OFF)
    # only disable rootmap for all but example01 and stltypes which require them
    if(NOT DICT STREQUAL "example01" AND NOT DICT STREQUAL "stltypes")
        set(CMAKE_ROOTTEST_NOROOTMAP ON)
    endif()

    set(dictname ${DICT}Dict)

    if(DICT STREQUAL "overloads" OR DICT STREQUAL "example01" )
        set(build_module FALSE)
        # Generate dictionary WITHOUT PCM/module
        ROOT_GENERATE_DICTIONARY(
          G__${dictname}
          ${CMAKE_CURRENT_SOURCE_DIR}/${DICT}.h
          LINKDEF ${CMAKE_CURRENT_SOURCE_DIR}/${DICT}.xml
        )
    else()
        set(build_module TRUE)
        # Generate dictionary WITH PCM/module
        ROOT_GENERATE_DICTIONARY(
          G__${dictname}
          ${CMAKE_CURRENT_SOURCE_DIR}/${DICT}.h
          LINKDEF ${CMAKE_CURRENT_SOURCE_DIR}/${DICT}.xml
          MODULE ${dictname}
        )
    endif()

    # Create necessary shared libraries for the tests
    add_library(${DICT}Dict SHARED ${CMAKE_CURRENT_SOURCE_DIR}/${DICT}.cxx G__${DICT}Dict.cxx)
    target_include_directories(${DICT}Dict PUBLIC ${Python3_INCLUDE_DIRS} ${CMAKE_BINARY_DIR}/include)
    target_link_libraries(${DICT}Dict PUBLIC ROOT::Core)
    target_compile_options(${DICT}Dict PRIVATE "-w")
    if(MSVC)
      target_link_libraries(${DICT}Dict PUBLIC Python3::Python)
      set_target_properties(${DICT}Dict PROPERTIES
          PREFIX "lib"
          IMPORT_PREFIX "lib"
      )
      #set_target_properties(${DICT}Dict PROPERTIES PREFIX "")
      set_target_properties(${DICT}Dict PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS TRUE)
    elseif(APPLE)
      target_link_libraries(${DICT}Dict PUBLIC -Wl,-bind_at_load -Wl,-w -Wl,-undefined -Wl,dynamic_lookup)
    else()
      target_link_libraries(${DICT}Dict PUBLIC -Wl,--unresolved-symbols=ignore-all)
    endif()

    # Like in ROOTTEST_GENERATE_DICTIONARY
    if(MSVC AND NOT CMAKE_GENERATOR MATCHES Ninja)
      add_custom_command(TARGET ${dictname} POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_BINARY_DIR}/$<CONFIG>/lib${dictname}.dll
                                         ${CMAKE_CURRENT_BINARY_DIR}/lib${dictname}.dll
        COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_BINARY_DIR}/$<CONFIG>/lib${dictname}.lib
                                         ${CMAKE_CURRENT_BINARY_DIR}/lib${dictname}.lib)
      if(build_module)
        add_custom_command(TARGET ${dictname} POST_BUILD
          COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_BINARY_DIR}/lib${dictname}_rdict.pcm
                                           ${CMAKE_CURRENT_BINARY_DIR}/$<CONFIG>/lib${dictname}_rdict.pcm)
      endif()
    endif()
endforeach()
unset(CMAKE_ROOTTEST_NOROOTMAP)

# Tests list
set(CPPYY_TESTS
test_aclassloader
test_advancedcpp
test_api
test_boost
test_conversions
test_crossinheritance
test_datatypes
test_doc_features
test_eigen
test_fragile
test_lowlevel
test_numba
test_overloads
test_pythonify
test_pythonization
test_regression
test_streams
test_templates
test_concurrent)

set(enable_test_cpp11features FALSE)
if(NOT CMAKE_CXX_STANDARD EQUAL 23)
    # test_cpp11features crashes when run with C++23
    if(APPLE)
        execute_process(
            COMMAND sw_vers -productVersion
            OUTPUT_VARIABLE MACOS_VERSION
            OUTPUT_STRIP_TRAILING_WHITESPACE
        )
        if(MACOS_VERSION VERSION_LESS "26.0")
            # Parts of this test suite fail with std::make_unique 
            set(enable_test_cpp11features TRUE)
        endif()
    else()
        set(enable_test_cpp11features TRUE)
    endif()
endif()

set(common_args GENERIC ENVIRONMENT ${ld_library_path}=${CMAKE_CURRENT_BINARY_DIR} PYTHON_DEPS pytest)

# Configure the tests
if(NOT MSVC)
  foreach(TEST ${CPPYY_TESTS})
    ROOT_ADD_PYUNITTEST(${TEST} ${TEST}.py ${common_args})
  endforeach()
endif()

ROOT_ADD_PYUNITTEST(test_operators test_operators.py
    ${common_args}
    COPY_TO_BUILDDIR operators.cxx operators.h
    PRECMD root -b -q -l -e .L\ operators.cxx+
)

if(enable_test_cpp11features)
  ROOT_ADD_PYUNITTEST(test_cpp11features test_cpp11features.py
      ${common_args}
      COPY_TO_BUILDDIR cpp11features.cxx cpp11features.h
      PRECMD root -b -q -l -e .L\ cpp11features.cxx+
  )
endif()

if(NOT APPLE)
  # test_stltypes completely fails on OSX
  # because it is unable to load the necessary shared libraries
  ROOT_ADD_PYUNITTEST(test_stltypes test_stltypes.py
      ${common_args}
      COPY_TO_BUILDDIR stltypes.cxx stltypes.h
      PRECMD root -b -q -l -e .L\ stltypes.cxx+
  )
endif()


# The leakcheck test is disabled due to its sporadic nature, especially fragile on VMs
# ROOT_ADD_PYUNITTEST(test_leakcheck test_leakcheck.py ${common_args})
