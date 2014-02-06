#include <math.h>

#include <iostream>

#include <TROOT.h>
#include <TClass.h>
#include <TStyle.h>
#include <TColor.h>

#include "TMPalette.hh"

ClassImp(TMPalette);

TMPalette::TMPalette(const char * name)
  :TNamed(name, name)
{
  fNColors = 0;
  fColors = 0;

  TObject * obj = gROOT->GetListOfSpecials()->First();
  while (obj) {
    if (strcmp(obj->GetName(), name)==0) {
      gROOT->GetListOfSpecials()->Remove(obj);
      delete obj;
      break;
    }
    obj = gROOT->GetListOfSpecials()->After(obj);
  }

  gROOT->GetListOfSpecials()->Add(this);
}

TMPalette::~TMPalette()
{
  if (fColors) delete [] fColors;
}
  
void TMPalette::CreateGradientColorTable(UInt_t Number, 
					 Double_t* Length, 
					 Double_t* Red, 
					 Double_t* Green,
					 Double_t* Blue,
					 UInt_t NColors)
{
  // Linear gradient color table:
  // Red, Green and Blue are several RGB colors with values from 0.0 .. 1.0.
  // Their number is "Intervals".
  // Length is the length of the color interval between the RGB-colors:
  // Imaging the whole gradient goes from 0.0 for the first RGB color to 1.0
  // for the last RGB color, then each "Length"-entry in between stands for
  // the length of the intervall between the according RGB colors.
  //
  // This definition is similar to the povray-definition of gradient
  // color tables.
  //
  // In order to create a color table do the following:
  // Define the RGB Colors:
  // > UInt_t Number = 5;
  // > Double_t Red[5]   = { 0.00, 0.09, 0.18, 0.09, 0.00 };
  // > Double_t Green[5] = { 0.01, 0.02, 0.39, 0.68, 0.97 };
  // > Double_t Blue[5]  = { 0.17, 0.39, 0.62, 0.79, 0.97 };
  // Define the length of the (color)-interval between this points
  // > Double_t Stops[5] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  // i.e. the color interval between Color 2 and Color 3 is
  // 0.79 - 0.62 => 17 % of the total palette area between these colors
  //
  //  Original code by Andreas Zoglauer <zog@mpe.mpg.de>

  UInt_t g, c;
  
  if (fColors) delete [] fColors;

  fNColors = 0;
  fColors = new Int_t[NColors+1];

  UInt_t nColorsGradient;
  TColor * color;
  Int_t highestIndex = 0;

  // Check if all RGB values are between 0.0 and 1.0 and
  // Length goes from 0.0 to 1.0 in increasing order.
  for (c = 0; c < Number; c++) {
    if (Red[c] < 0 || Red[c] > 1.0 ||
	Green[c] < 0 || Green[c] > 1.0 ||
	Blue[c] < 0 || Blue[c] > 1.0 ||
	Length[c] < 0 || Length[c] > 1.0) {
      //Error("CreateGradientColorTable",
      //      "All RGB colors and interval lengths have to be between 0.0 and 1.0");
      
      delete [] fColors;
      fColors = 0;
      
      return;
    }
    if (c >= 1) {
      if (Length[c-1] > Length[c]) {
	//Error("CreateGradientColorTable",
	//      "The interval lengths have to be in increasing order");
	
	delete [] fColors;
	fColors = 0;

	return;
      }
    }
  }
  
  TColor::InitializeColors();

  // Search for the highest color index not used in ROOT:
  // We do not want to overwrite some colors...
  TSeqCollection * colorTable = gROOT->GetListOfColors();
  if ((color = (TColor *) colorTable->Last()) != 0) {
    if (color->GetNumber() > highestIndex) {
      highestIndex = color->GetNumber();
    }
    while ((color = (TColor *) (colorTable->Before(color))) != 0) {
      if (color->GetNumber() > highestIndex) {
	highestIndex = color->GetNumber();
      }
    }
  }
  highestIndex++;

  // Now create the colors and add them to the default palette:
  
  // For each defined gradient...
  for (g = 1; g < Number; g++) {
    // create the colors...
    nColorsGradient = (Int_t) (floor(NColors*Length[g]) - floor(NColors*Length[g-1]));
    for (c = 0; c < nColorsGradient; c++) {
      color = new TColor(highestIndex,
			 Red[g-1] + c * (Red[g] - Red[g-1])/ nColorsGradient,
			 Green[g-1] + c * (Green[g] - Green[g-1])/ nColorsGradient,
			 Blue[g-1] + c * (Blue[g] - Blue[g-1])/ nColorsGradient,
			 Form("%s_%d", GetName(), fNColors));
      fColors[fNColors] = highestIndex;
      fNColors++;
      highestIndex++;
    }
  }

  fNColors--;

  TColor::SetPalette(fNColors, fColors);
}
 
void TMPalette::cd()
{
  gStyle->SetNumberContours(fNColors);
  gStyle->SetPalette(fNColors, fColors);
}
 
void TMPalette::SetPalette(const char * name)
{
  TObject * obj = gROOT->GetListOfSpecials()->First();
  while (obj) {
    if (strcmp(obj->GetName(), name)==0 && 
	obj->IsA()->InheritsFrom(TMPalette::Class())) {
      ((TMPalette*)obj)->cd();
      break;
    }
    obj = gROOT->GetListOfSpecials()->After(obj);
  }
}
