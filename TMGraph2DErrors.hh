#ifndef TMGraph2DErrors_hh
#define TMGraph2DErrors_hh

#include <TGraph2DErrors.h>

class TMGraph2DErrors : public TGraph2DErrors
{
public:
  
  TMGraph2DErrors();
  TMGraph2DErrors(Int_t n, const Int_t *p);
  TMGraph2DErrors(Int_t n, const Float_t *p);
  TMGraph2DErrors(Int_t n, const Double_t *p);

  ClassDef(TMGraph2DErrors,0)
};

#endif


