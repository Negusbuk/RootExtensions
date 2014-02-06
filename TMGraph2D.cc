
#include "TMGraph2D.hh"

ClassImp(TMGraph2D)

TMGraph2D::TMGraph2D()
  :TGraph2D()
{

}

TMGraph2D::TMGraph2D(Int_t n, const Int_t *p)
  :TGraph2D(n)
{
  if (!p) return;
  for (Int_t i=0;i<n;i++) {
    SetPoint(i,p[i*3],p[i*3+1],p[i*3+2]);
  }
}

TMGraph2D::TMGraph2D(Int_t n, const Float_t *p)
  :TGraph2D(n)
{
  if (!p) return;
  for (Int_t i=0;i<n;i++) {
    SetPoint(i,p[i*3],p[i*3+1],p[i*3+2]);
  }
}

TMGraph2D::TMGraph2D(Int_t n, const Double_t *p)
  :TGraph2D(n)
{
  if (!p) return;
  for (Int_t i=0;i<n;i++) {
    SetPoint(i,p[i*3],p[i*3+1],p[i*3+2]);
  }
}
