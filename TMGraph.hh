#ifndef TMGraph_hh
#define TMGraph_hh

#include <TGraph.h>

class TMGraph : public TGraph
{
public:
  
  TMGraph();
  TMGraph(Int_t n, const Int_t *p);
  TMGraph(Int_t n, const Float_t *p);
  TMGraph(Int_t n, const Double_t *p);

  ClassDef(TMGraph,0)
};

#endif


