#ifndef TArcArrow_hh
#define TArcArrow_hh

#include <TEllipse.h>
#include <TObject.h>
#include <TAttLine.h>
#include <TAttFill.h>

class TArcArrow : public TEllipse
{
public:
  
  TArcArrow();
  TArcArrow(Double_t x1, Double_t y1,
	    Double_t r1, Double_t r2=0,
	    Double_t phimin=0, Double_t phimax=360,
	    Double_t theta=0, Float_t arrowsize=0.05,
	    Option_t *option=">");
  TArcArrow(const TArcArrow &arc);
  virtual ~TArcArrow();
  void Copy(TObject &obj) const;

  Double_t       GetX1() const {return fX1;}
  Double_t       GetY1() const {return fY1;}
  Double_t       GetR1() const {return fR1;}
  Double_t       GetR2() const {return fR2;}
  Double_t       GetPhimin() const {return fPhimin;}
  Double_t       GetPhimax() const {return fPhimax;}
  Double_t       GetTheta() const  {return fTheta;}
  Float_t        GetAngle() const {return fAngle;}
  Float_t        GetArrowSize() const {return fArrowSize;}
  Option_t      *GetOption() const { return fOption.Data();}
  
  virtual void   SetX1(Double_t x1) {fX1=x1;} // *MENU*
  virtual void   SetY1(Double_t y1) {fY1=y1;} // *MENU*
  virtual void   SetR1(Double_t r1) {fR1=r1;} // *MENU*
  virtual void   SetR2(Double_t r2) {fR2=r2;} // *MENU*
  virtual void   SetPhimin(Double_t phi=0)   {fPhimin=phi;} // *MENU*
  virtual void   SetPhimax(Double_t phi=360) {fPhimax=phi;} // *MENU*
  virtual void   SetTheta(Double_t theta=0) {fTheta=theta;} // *MENU*
  virtual void   SetAngle(Float_t angle=60) {fAngle=angle;} // *MENU*
  virtual void   SetArrowSize(Float_t arrowsize=0.05) {fArrowSize=arrowsize;} // *MENU*
  virtual void   SetOption(Option_t *option=">"){ fOption = option;} // *MENU*
  
  virtual Int_t  DistancetoPrimitive(Int_t px, Int_t py);
  virtual void   ExecuteEvent(Int_t event, Int_t px, Int_t py);

  virtual void   Draw(Option_t *option="");
  virtual void   DrawArcArrow(Double_t x1, Double_t y1,
			      Double_t r1,Double_t r2,
			      Double_t phimin, Double_t phimax,
			      Double_t theta, Float_t arrowsize=0.0,
			      Option_t *option="");
  
  virtual void   Paint(Option_t *option="");
  virtual void   PaintArcArrow(Double_t x1, Double_t y1, 
			       Double_t r1, Double_t r2,
			       Double_t phimin, Double_t phimax,
			       Double_t theta, Float_t arrowsize=0.0,
			       Option_t *option=">");
  
protected:
  
  Double_t    fX1;           //X coordinate of centre
  Double_t    fY1;           //Y coordinate of centre
  Double_t    fR1;           //first radius
  Double_t    fR2;           //second radius
  Double_t    fPhimin;       //Minimum angle (degrees)
  Double_t    fPhimax;       //Maximum angle (degrees)
  Double_t    fTheta;        //Rotation angle (degrees)
  Float_t     fAngle;        //Arrow opening angle (degrees)
  Float_t     fArrowSize;    //Arrow Size
  TString     fOption;       //Arrow shapes
  
  static Float_t      fgDefaultAngle;        //default Arrow opening angle (degrees)
  static Float_t      fgDefaultArrowSize;    //default Arrow Size
  static TString      fgDefaultOption;       //default Arrow shapes

  ClassDef(TArcArrow, 1);
};

#endif
