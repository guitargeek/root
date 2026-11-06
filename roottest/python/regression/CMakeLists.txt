if(pyroot)
if(NOT MSVC OR win_broken_tests)
  ROOTTEST_ADD_TEST(regression
                    MACRO PyROOT_regressiontests.py OPTS -b
                    LABELS longtest)
endif()

  ROOTTEST_ADD_TEST(gh_21536
                    MACRO gh_21536.py)

  ROOTTEST_GENERATE_REFLEX_DICTIONARY(gh_24176
                                      gh_24176.h
                                      SELECTION gh_24176_selection.xml
                                      FIXTURES_SETUP python-regression-gh_24176-dict-fixture)

  ROOTTEST_ADD_TEST(gh_24176_enum
                    MACRO gh_24176.py
                    OPTS enum
                    FIXTURES_REQUIRED python-regression-gh_24176-dict-fixture)

  ROOTTEST_ADD_TEST(gh_24176_variable
                    MACRO gh_24176.py
                    OPTS variable
                    FIXTURES_REQUIRED python-regression-gh_24176-dict-fixture)
endif()
