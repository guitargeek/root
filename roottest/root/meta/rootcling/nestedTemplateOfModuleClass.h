#ifndef ROOTTEST_NESTEDTEMPLATEOFMODULECLASS_H
#define ROOTTEST_NESTEDTEMPLATEOFMODULECLASS_H

// ROOT::Detail::TCollectionProxyInfo comes from the C++ module of Core, which
// this dictionary depends on. Its nested class templates get instantiated for
// this dictionary only. rootcling doesn't look at all the members of classes
// from other modules, but it must still find these instantiations.
#include "TCollectionProxyInfo.h"

struct NestedTemplateArg {
   int fX = 0;
};

#endif
