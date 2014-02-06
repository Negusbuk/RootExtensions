
#include "TMGraph.hh"

ClassImp(TMGraph)

TMGraph::TMGraph()
  :TGraph()
{

}

TMGraph::TMGraph(Int_t n, const Int_t *p)
  :TGraph(n)
{
  if (!p) return;
  for (Int_t i=0;i<n;i++) {
    fX[i] = p[i*2];
    fY[i] = p[i*2+1];
  }
}

TMGraph::TMGraph(Int_t n, const Float_t *p)
  :TGraph(n)
{
  if (!p) return;
  for (Int_t i=0;i<n;i++) {
    fX[i] = p[i*2];
    fY[i] = p[i*2+1];
  }
}

TMGraph::TMGraph(Int_t n, const Double_t *p)
  :TGraph(n)
{
  if (!p) return;
  for (Int_t i=0;i<n;i++) {
    fX[i] = p[i*2];
    fY[i] = p[i*2+1];
  }
}
