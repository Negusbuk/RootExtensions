#ifndef TMGraphErrors_hh
#define TMGraphErrors_hh

#include <TGraphErrors.h>

class TMGraphErrors : public TGraphErrors
{
  public:

  TMGraphErrors();
  TMGraphErrors(Int_t n, const Int_t *p);
  TMGraphErrors(Int_t n, const Float_t *p);
  TMGraphErrors(Int_t n, const Double_t *p);
  
  ClassDef(TMGraphErrors,0)
};

#endif


