#ifndef MAKECLASSRULE_H
#define MAKECLASSRULE_H

#include <Rtypes.h>

// Version 1 of the class, used to write the file.
struct MakeClassRule {
   Int_t fDropped = 0;
   Int_t fRun = 0;
   Int_t fSrc = 0;

   ClassDef(MakeClassRule, 1)
};

#endif
