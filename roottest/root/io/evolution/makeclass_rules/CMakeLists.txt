ROOTTEST_GENERATE_DICTIONARY(libRoottestIoEvolutionMakeClassRuleV1
                             v1/MakeClassRule.h
                             LINKDEF v1/MakeClassRuleLinkDef.h
                             NO_CXXMODULE NO_ROOTMAP
                             FIXTURES_SETUP root-io-evolution-makeclass_rules-lib1)

ROOTTEST_GENERATE_DICTIONARY(libRoottestIoEvolutionMakeClassRuleV2
                             v2/MakeClassRule.h
                             LINKDEF v2/MakeClassRuleLinkDef.h
                             NO_CXXMODULE NO_ROOTMAP
                             FIXTURES_SETUP root-io-evolution-makeclass_rules-lib2)

ROOTTEST_ADD_TEST(write
                  MACRO execWriteMakeClassRule.C
                  OUTREF execWriteMakeClassRule.ref
                  LABELS io
                  FIXTURES_REQUIRED root-io-evolution-makeclass_rules-lib1
                  FIXTURES_SETUP root-io-evolution-makeclass_rules-write)

ROOTTEST_ADD_TEST(read
                  MACRO execReadMakeClassRule.C
                  OUTREF execReadMakeClassRule.ref
                  LABELS io
                  FIXTURES_REQUIRED root-io-evolution-makeclass_rules-lib2 root-io-evolution-makeclass_rules-write)
