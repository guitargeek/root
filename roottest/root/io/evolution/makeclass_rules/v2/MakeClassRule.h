#ifndef MAKECLASSRULE_H
#define MAKECLASSRULE_H

#include <Rtypes.h>

// Version 2 of the class, used to read the file. Compared to version 1,
// fDropped was removed, so the in-memory offsets of fRun and fSrc differ from
// their offsets in the onfile object.
struct MakeClassRule {
   Int_t fRun = 0;
   Int_t fSrc = 0;
   Int_t fTwice = 0;

   ClassDef(MakeClassRule, 2)
};

#endif
