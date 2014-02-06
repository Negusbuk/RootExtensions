#ifndef TMGraph2D_hh
#define TMGraph2D_hh

#include <TGraph2D.h>

class TMGraph2D : public TGraph2D
{
public:
  
  TMGraph2D();
  TMGraph2D(Int_t n, const Int_t *p);
  TMGraph2D(Int_t n, const Float_t *p);
  TMGraph2D(Int_t n, const Double_t *p);

  ClassDef(TMGraph2D,0)
};

#endif


