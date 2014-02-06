
#include "TMGraphErrors.hh"

ClassImp(TMGraphErrors)

TMGraphErrors::TMGraphErrors()
  :TGraphErrors()
{

}

TMGraphErrors::TMGraphErrors(Int_t n, const Int_t *p)
  :TGraphErrors(n)
{
  if (!p) return;
  for (Int_t i=0;i<n;i++) {
    fX[i]  = p[i*4];
    fY[i]  = p[i*4+1];
    fEX[i] = p[i*4+2];
    fEY[i] = p[i*4+3];
  }
}

TMGraphErrors::TMGraphErrors(Int_t n, const Float_t *p)
  :TGraphErrors(n)
{
  if (!p) return;
  for (Int_t i=0;i<n;i++) {
    fX[i]  = p[i*4];
    fY[i]  = p[i*4+1];
    fEX[i] = p[i*4+2];
    fEY[i] = p[i*4+3];
  }
}

TMGraphErrors::TMGraphErrors(Int_t n, const Double_t *p)
  :TGraphErrors(n)
{
  if (!p) return;
  for (Int_t i=0;i<n;i++) {
    fX[i]  = p[i*4];
    fY[i]  = p[i*4+1];
    fEX[i] = p[i*4+2];
    fEY[i] = p[i*4+3];
  }
}
