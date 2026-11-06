if(pyroot)
if(NOT MSVC OR win_broken_tests)
  ROOTTEST_ADD_TEST(regression
                    MACRO PyROOT_regressiontests.py OPTS -b
                    LABELS longtest)
endif()

  ROOTTEST_ADD_TEST(gh_21536
                    MACRO gh_21536.py)
endif()
