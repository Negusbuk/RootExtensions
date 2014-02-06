
#include "TMGraph2DErrors.hh"

ClassImp(TMGraph2DErrors)

TMGraph2DErrors::TMGraph2DErrors()
  :TGraph2DErrors()
{
  
}

TMGraph2DErrors::TMGraph2DErrors(Int_t n, const Int_t *p)
  :TGraph2DErrors(n)
{
  if (!p) return;
  for (Int_t i=0;i<n;i++) {
    SetPoint(i,p[i*6],p[i*6+1],p[i*6+2]);
    SetPointError(i,p[i*6+3],p[i*6+4],p[i*6+5]);
  }
}

TMGraph2DErrors::TMGraph2DErrors(Int_t n, const Float_t *p)
  :TGraph2DErrors(n)
{
  if (!p) return;
  for (Int_t i=0;i<n;i++) {
    SetPoint(i,p[i*6],p[i*6+1],p[i*6+2]);
    SetPointError(i,p[i*6+3],p[i*6+4],p[i*6+5]);
  }
}

TMGraph2DErrors::TMGraph2DErrors(Int_t n, const Double_t *p)
  :TGraph2DErrors(n)
{
  if (!p) return;
  for (Int_t i=0;i<n;i++) {
    SetPoint(i,p[i*6],p[i*6+1],p[i*6+2]);
    SetPointError(i,p[i*6+3],p[i*6+4],p[i*6+5]);
  }
}
