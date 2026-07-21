if(pyroot)
  ROOTTEST_ADD_TEST(basic
                    MACRO PyROOT_basictests.py
                    COPY_TO_BUILDDIR SimpleClass.C)
endif()
