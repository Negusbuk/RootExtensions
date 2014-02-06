#ifndef TMPalette_hh
#define TMPalette_hh

#include <TNamed.h>

class TMPalette : public TNamed
{
public:
  
  TMPalette(const char * name);
  virtual ~TMPalette();
  
  void                         CreateGradientColorTable(UInt_t Number, 
							Double_t* Length, 
							Double_t* Red, 
							Double_t* Green,
							Double_t* Blue,
							UInt_t NColors);

  virtual void                 cd();
  static void                  SetPalette(const char * name);
  
protected:

  UInt_t                      fNColors;
  Int_t                      *fColors;

  ClassDef(TMPalette,0);
};

#endif


