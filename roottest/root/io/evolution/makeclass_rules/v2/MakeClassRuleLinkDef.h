#ifdef __CLING__

#pragma link C++ class MakeClassRule+;

// Rule targeting a member that is also read from the file.
#pragma read sourceClass="MakeClassRule" version="[1]" targetClass="MakeClassRule" \
   source="Int_t fRun" target="fRun" code="{ fRun = onfile.fRun + 1000; }"

// Rule using a member as source that is otherwise read unchanged.
#pragma read sourceClass="MakeClassRule" version="[1]" targetClass="MakeClassRule" \
   source="Int_t fSrc" target="fTwice" code="{ fTwice = 2 * onfile.fSrc; }"

#endif
