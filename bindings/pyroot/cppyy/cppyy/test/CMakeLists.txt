# Copyright (C) 1995-2025, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

# Dictionary list
set(CPPYY_TEST_DICTS
advancedcpp
advancedcpp2
pythonizables
)

foreach(DICT ${CPPYY_TEST_DICTS})
    set(CMAKE_ROOTTEST_NOROOTMAP ON)

    set(dictname ${DICT}Dict)

    # Generate dictionary WITH PCM/module
    ROOT_GENERATE_DICTIONARY(
      G__${dictname}
      ${CMAKE_CURRENT_SOURCE_DIR}/${DICT}.h
      LINKDEF ${CMAKE_CURRENT_SOURCE_DIR}/${DICT}.xml
      MODULE ${dictname}
    )

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

endforeach()
unset(CMAKE_ROOTTEST_NOROOTMAP)

# Tests list
set(CPPYY_TESTS
test_advancedcpp # requires two dictionaries
test_api
test_boost
test_eigen
test_numba
test_pythonization
test_regression
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

if(MSVC)
  set(common_args GENERIC PYTHON_DEPS pytest)
else()
  set(common_args GENERIC ENVIRONMENT ${ld_library_path}=${CMAKE_CURRENT_BINARY_DIR} PYTHON_DEPS pytest)
endif()

# Configure the tests
foreach(TEST ${CPPYY_TESTS})
  ROOT_ADD_PYUNITTEST(${TEST} ${TEST}.py ${common_args})
endforeach()

function(ROOT_ADD_PYUNITTEST_ACLIC TESTNAME)
  set(options)
  set(oneValueArgs)
  set(multiValueArgs DICTNAME)
  cmake_parse_arguments(ARG
    "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN}
  )

  # Common test args
  set(_common_args ${common_args})

  set(_precmd $<TARGET_FILE:root.exe> -b -q -l -e ".L ${ARG_DICTNAME}.cxx+")

  ROOT_ADD_PYUNITTEST(${TESTNAME} ${TESTNAME}.py
    ${_common_args}
    COPY_TO_BUILDDIR ${ARG_DICTNAME}.cxx ${ARG_DICTNAME}.h
    PRECMD ${_precmd}
  )
endfunction()

# doesn't pass on Windows
ROOT_ADD_PYUNITTEST_ACLIC(
  test_doc_features
  DICTNAME doc_helper
)

ROOT_ADD_PYUNITTEST_ACLIC(
  test_operators
  DICTNAME operators
)

ROOT_ADD_PYUNITTEST_ACLIC(
  test_streams
  DICTNAME std_streams
)

ROOT_ADD_PYUNITTEST_ACLIC(
  test_crossinheritance
  DICTNAME crossinheritance
)

ROOT_ADD_PYUNITTEST_ACLIC(
  test_conversions
  DICTNAME conversions
)

ROOT_ADD_PYUNITTEST_ACLIC(
  test_fragile
  DICTNAME fragile
)

ROOT_ADD_PYUNITTEST(test_pythonify test_pythonify.py
  ${common_args}
  COPY_TO_BUILDDIR example01.cxx example01.h
  PRECMD $<TARGET_FILE:root.exe> -b -q -l -e ".L example01.cxx+"
  FIXTURES_SETUP example01_dict_fixture
)

ROOT_ADD_PYUNITTEST(test_lowlevel test_lowlevel.py
  ${common_args}
  COPY_TO_BUILDDIR datatypes.cxx datatypes.h
  PRECMD $<TARGET_FILE:root.exe> -b -q -l -e ".L datatypes.cxx+"
  FIXTURES_SETUP datatypes_dict_fixture
)

if(NOT MSVC OR win_broken_tests)
  ROOT_ADD_PYUNITTEST(test_datatypes test_datatypes.py ${common_args} FIXTURES_REQUIRED datatypes_dict_fixture)

  ROOT_ADD_PYUNITTEST_ACLIC(
    test_overloads
    DICTNAME overloads
  )
endif()

ROOT_ADD_PYUNITTEST(test_aclassloader test_aclassloader.py ${common_args} FIXTURES_REQUIRED example01_dict_fixture)

# doesn't pass on Windows
ROOT_ADD_PYUNITTEST_ACLIC(
  test_templates
  DICTNAME templates
)

# doesn't pass on Windows
if(enable_test_cpp11features)
  ROOT_ADD_PYUNITTEST_ACLIC(
    test_cpp11features
    DICTNAME cpp11features
  )
endif()

if(NOT APPLE)
  # test_stltypes completely fails on OSX
  # because it is unable to load the necessary shared libraries
  ROOT_ADD_PYUNITTEST_ACLIC(
    test_stltypes
    DICTNAME stltypes
  )
endif()

# The leakcheck test is disabled due to its sporadic nature, especially fragile on VMs
# ROOT_ADD_PYUNITTEST(test_leakcheck test_leakcheck.py ${common_args})
