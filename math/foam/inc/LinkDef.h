#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class TFoamIntegrand+;
#pragma link C++ class TFoamMaxwt+;
#pragma link C++ class TFoamVect+;
#pragma link C++ class TFoamCell+;
#pragma link C++ class TFoam+;
#pragma link C++ class TFoamSampler+;
// Note that just like in the rule for version 2 and beyond, SetCells() must be
// called in the end for each cell.
#pragma read sourceClass="TFoam" targetClass="TFoam" version="[1]" \
  source="Int_t fNCells; TFoamCell **fCells; TRefArray *fCellsAct" target="fNCells,fCells,fCellsAct"\
  include="TRefArray.h" \
  code="{fNCells = onfile.fNCells; \
         fCells = onfile.fCells; \
         onfile.fCells = nullptr; \
         fCellsAct.clear(); \
         for (Int_t i=0; i < onfile.fCellsAct->GetEntries(); ++i) { \
            const TObject* cellp = onfile.fCellsAct->At(i); \
            for (Int_t j=0; j < fNCells; ++j) { \
               if (cellp == fCells[j]) { \
                 fCellsAct.push_back(j); \
                 break; \
               } \
            } \
         } \
         if(fCells) { \
           for (Int_t i=0; i < fNCells; ++i) { \
              fCells[i]->SetCells(fCells); \
           } \
         } \
  }";
// To avoid duplicating the cells array, the TFoamCell member that points to
// the full cell array is filled in this custom rule. Otherwise, IO would not
// know that they are supposed to point to the same array.
#pragma read sourceClass="TFoam" targetClass="TFoam" version="[2-]" \
  source="Int_t fNCells; TFoamCell **fCells" target="fNCells,fCells"\
  code="{fNCells = onfile.fNCells; \
         fCells = onfile.fCells; \
         onfile.fCells = nullptr; \
         if(fCells) { \
           for (Int_t i=0; i < fNCells; ++i) { \
              fCells[i]->SetCells(fCells); \
           } \
         } \
  }";
#endif
