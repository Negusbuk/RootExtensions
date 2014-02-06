#include <iostream>

#include <Riostream.h>
#include <TROOT.h>
#include <TVirtualPad.h>
#include <TMath.h>
#include <TClass.h>

#include "TArcArrow.hh"

const Double_t kPI = 3.14159265358979323846;

Float_t TArcArrow::fgDefaultAngle      = 60;
Float_t TArcArrow::fgDefaultArrowSize  = 0.05;
TString TArcArrow::fgDefaultOption     = ">";

ClassImp(TArcArrow);

TArcArrow::TArcArrow()
  :TEllipse()
{
  fX1 = 0;
  fY1 = 0;
  fR1 = 1;
  fR2 = 1;
  fPhimin = 0;
  fPhimax = 360;
  fTheta  = 0;
  
  fAngle = fgDefaultAngle;
  fgDefaultArrowSize = fgDefaultArrowSize;
  fOption = fgDefaultOption;
}

//______________________________________________________________________________
TArcArrow::TArcArrow(Double_t x1, Double_t y1,
		     Double_t r1, Double_t r2,
		     Double_t phimin, Double_t phimax,
		     Double_t theta, Float_t arrowsize,
		     Option_t *option)
  :TEllipse()
{
  // Arc  normal constructor.
  //
  //  x1,y1  : coordinates of centre of arc
  //  r1     : arc radius
  //  phimin : min and max angle in degrees (default is 0-->360)
  //  phimax :
  //
  //  When a circle sector only is drawn, the lines connecting the center
  //  of the arc to the edges are drawn by default. One can specify
  //  the drawing option "only" to not draw these lines.

  fX1     = x1;
  fY1     = y1;
  fR1     = r1;
  fR2     = r2;
  fPhimin = phimin;
  fPhimax = phimax;
  fTheta  = theta;
  if (r2 <= 0) fR2 = fR1;
  
  fAngle       = fgDefaultAngle;
  fArrowSize   = arrowsize;
  fOption      = option;
  SetFillColor(GetLineColor());
  SetFillStyle(1001);
}

//______________________________________________________________________________
TArcArrow::TArcArrow(const TArcArrow &arc)
  :TEllipse(arc)
{
  ((TArcArrow&)arc).Copy(*this);
}


//______________________________________________________________________________
TArcArrow::~TArcArrow()
{
  // Arc default destructor.
}


//______________________________________________________________________________
void TArcArrow::Copy(TObject &obj) const
{
  TEllipse::Copy(obj);
  ((TArcArrow&)obj).fX1 = fX1;
  ((TArcArrow&)obj).fY1 = fY1;
  ((TArcArrow&)obj).fR1 = fR1;
  ((TArcArrow&)obj).fR2 = fR2;
  ((TArcArrow&)obj).fPhimin = fPhimin;
  ((TArcArrow&)obj).fPhimax = fPhimax;
  ((TArcArrow&)obj).fTheta  = fTheta;
  ((TArcArrow&)obj).fAngle  = fAngle;
  ((TArcArrow&)obj).fArrowSize  = fArrowSize;
  ((TArcArrow&)obj).fOption  = fOption;
}

Int_t TArcArrow::DistancetoPrimitive(Int_t px, Int_t py)
{
   // Compute distance from point px,py to an ellipse.
   //
   //  Compute the closest distance of approach from point px,py to this ellipse.
   //  The distance is computed in pixels units.

   Double_t x = gPad->PadtoX(gPad->AbsPixeltoX(px));
   Double_t y = gPad->PadtoY(gPad->AbsPixeltoY(py));

   Double_t dxnr = x - fX1;
   Double_t dynr = y - fY1;

   Double_t ct = TMath::Cos(kPI*GetTheta()/180.0);
   Double_t st = TMath::Sin(kPI*GetTheta()/180.0);

   Double_t dx =  dxnr*ct + dynr*st;
   Double_t dy = -dxnr*st + dynr*ct;

   Double_t r1 = fR1;
   Double_t r2 = fR2;

   if (dx == 0 || r1 == 0 || r2 == 0) return 9999;
   Double_t distp = TMath::Sqrt(dx*dx + dy*dy);

   Double_t tana = dy/dx;
   tana *= tana;
   Double_t distr = TMath::Sqrt((1+tana)/(1.0/(r1*r1) + tana/(r2*r2)));
   Int_t dist = 9999;
   if (GetFillColor() && GetFillStyle()) {
      if (distr > distp) dist = 0;
   } else {
      if (TMath::Abs(distr-distp)/(r1+r2) < 0.01) dist = 0;
   }
   return dist;
}

void TArcArrow::ExecuteEvent(Int_t event, Int_t px, Int_t py)
{
  // Execute action corresponding to one event.
  //
  //  This member function is called when a line is clicked with the locator
  //
  //  If Left button clicked on one of the line end points, this point
  //     follows the cursor until button is released.
  //
  //  if Middle button clicked, the line is moved parallel to itself
  //     until the button is released.
  //
  //  NOTE that support for log scale is not implemented

  Int_t kMaxDiff = 10;
  const Int_t kMinSize = 25;
  const Int_t np = 40;
  static Int_t x[np+2], y[np+2];
  static Int_t px1,py1,npe,r1,r2,sav1,sav2;
  static Int_t pxold, pyold;
  static Int_t sig,impair;
  Int_t i, dpx, dpy;
  Double_t angle,dx,dy,dphi,ct,st,fTy,fBy,fLx,fRx;
  static Bool_t pTop, pL, pR, pBot, pINSIDE;
  static Int_t pTx,pTy,pLx,pLy,pRx,pRy,pBx,pBy;

  if (!gPad->IsEditable()) return;

  switch (event) {

  case kButton1Down:
    gVirtualX->SetLineColor(-1);
    TAttLine::Modify();
    dphi = (fPhimax-fPhimin)*kPI/(180*np);
    ct   = TMath::Cos(kPI*fTheta/180);
    st   = TMath::Sin(kPI*fTheta/180);
    for (i=0;i<np;i++) {
      angle = fPhimin*kPI/180 + Double_t(i)*dphi;
      dx    = fR1*TMath::Cos(angle);
      dy    = fR2*TMath::Sin(angle);
      x[i]  = gPad->XtoAbsPixel(fX1 + dx*ct - dy*st);
      y[i]  = gPad->YtoAbsPixel(fY1 + dx*st + dy*ct);
    }
    if (fPhimax-fPhimin >= 360 ) {
      x[np] = x[0];
      y[np] = y[0];
      npe = np;
    } else {
      x[np]   = gPad->XtoAbsPixel(fX1);
      y[np]   = gPad->YtoAbsPixel(fY1);
      x[np+1] = x[0];
      y[np+1] = y[0];
      npe = np + 1;
    }
    impair = 0;
    px1 = gPad->XtoAbsPixel(fX1);
    py1 = gPad->YtoAbsPixel(fY1);
    pTx = pBx = px1;
    pLy = pRy = py1;
    pTy = gPad->YtoAbsPixel(fR2+fY1);
    pBy = gPad->YtoAbsPixel(-fR2+fY1);
    pLx = gPad->XtoAbsPixel(-fR1+fX1);
    pRx = gPad->XtoAbsPixel(fR1+fX1);
    r2 = (pBy-pTy)/2;
    r1 = (pRx-pLx)/2;
    gVirtualX->DrawLine(pRx+4, py1+4, pRx-4, py1+4);
    gVirtualX->DrawLine(pRx-4, py1+4, pRx-4, py1-4);
    gVirtualX->DrawLine(pRx-4, py1-4, pRx+4, py1-4);
    gVirtualX->DrawLine(pRx+4, py1-4, pRx+4, py1+4);
    gVirtualX->DrawLine(pLx+4, py1+4, pLx-4, py1+4);
    gVirtualX->DrawLine(pLx-4, py1+4, pLx-4, py1-4);
    gVirtualX->DrawLine(pLx-4, py1-4, pLx+4, py1-4);
    gVirtualX->DrawLine(pLx+4, py1-4, pLx+4, py1+4);
    gVirtualX->DrawLine(px1+4, pBy+4, px1-4, pBy+4);
    gVirtualX->DrawLine(px1-4, pBy+4, px1-4, pBy-4);
    gVirtualX->DrawLine(px1-4, pBy-4, px1+4, pBy-4);
    gVirtualX->DrawLine(px1+4, pBy-4, px1+4, pBy+4);
    gVirtualX->DrawLine(px1+4, pTy+4, px1-4, pTy+4);
    gVirtualX->DrawLine(px1-4, pTy+4, px1-4, pTy-4);
    gVirtualX->DrawLine(px1-4, pTy-4, px1+4, pTy-4);
    gVirtualX->DrawLine(px1+4, pTy-4, px1+4, pTy+4);
    // No break !!!

  case kMouseMotion:
    px1 = gPad->XtoAbsPixel(fX1);
    py1 = gPad->YtoAbsPixel(fY1);
    pTx = pBx = px1;
    pLy = pRy = py1;
    pTy = gPad->YtoAbsPixel(fR2+fY1);
    pBy = gPad->YtoAbsPixel(-fR2+fY1);
    pLx = gPad->XtoAbsPixel(-fR1+fX1);
    pRx = gPad->XtoAbsPixel(fR1+fX1);
    pTop = pL = pR = pBot = pINSIDE = kFALSE;
    if ((TMath::Abs(px - pTx) < kMaxDiff) &&
	(TMath::Abs(py - pTy) < kMaxDiff)) {             // top edge
      pTop = kTRUE;
      gPad->SetCursor(kTopSide);
    }
    else
      if ((TMath::Abs(px - pBx) < kMaxDiff) &&
          (TMath::Abs(py - pBy) < kMaxDiff)) {             // bottom edge
	pBot = kTRUE;
	gPad->SetCursor(kBottomSide);
      }
      else
	if ((TMath::Abs(py - pLy) < kMaxDiff) &&
	    (TMath::Abs(px - pLx) < kMaxDiff)) {             // left edge
	  pL = kTRUE;
	  gPad->SetCursor(kLeftSide);
	}
	else
	  if ((TMath::Abs(py - pRy) < kMaxDiff) &&
	      (TMath::Abs(px - pRx) < kMaxDiff)) {             // right edge
	    pR = kTRUE;
	    gPad->SetCursor(kRightSide);
	  }
	  else {pINSIDE= kTRUE; gPad->SetCursor(kMove); }
    pxold = px;  pyold = py;

    break;

  case kButton1Motion:
    gVirtualX->DrawLine(pRx+4, py1+4, pRx-4, py1+4);
    gVirtualX->DrawLine(pRx-4, py1+4, pRx-4, py1-4);
    gVirtualX->DrawLine(pRx-4, py1-4, pRx+4, py1-4);
    gVirtualX->DrawLine(pRx+4, py1-4, pRx+4, py1+4);
    gVirtualX->DrawLine(pLx+4, py1+4, pLx-4, py1+4);
    gVirtualX->DrawLine(pLx-4, py1+4, pLx-4, py1-4);
    gVirtualX->DrawLine(pLx-4, py1-4, pLx+4, py1-4);
    gVirtualX->DrawLine(pLx+4, py1-4, pLx+4, py1+4);
    gVirtualX->DrawLine(px1+4, pBy+4, px1-4, pBy+4);
    gVirtualX->DrawLine(px1-4, pBy+4, px1-4, pBy-4);
    gVirtualX->DrawLine(px1-4, pBy-4, px1+4, pBy-4);
    gVirtualX->DrawLine(px1+4, pBy-4, px1+4, pBy+4);
    gVirtualX->DrawLine(px1+4, pTy+4, px1-4, pTy+4);
    gVirtualX->DrawLine(px1-4, pTy+4, px1-4, pTy-4);
    gVirtualX->DrawLine(px1-4, pTy-4, px1+4, pTy-4);
    gVirtualX->DrawLine(px1+4, pTy-4, px1+4, pTy+4);
    for (i=0;i<npe;i++) gVirtualX->DrawLine(x[i], y[i], x[i+1], y[i+1]);
    if (pTop) {
      sav1 = py1;
      sav2 = r2;
      py1 += (py - pyold)/2;
      r2 -= (py - pyold)/2;
      if (TMath::Abs(pyold-py)%2==1) impair++;
      if (py-pyold>0) sig=+1;
      else sig=-1;
      if (impair==2) { impair = 0; py1 += sig; r2 -= sig;}
      if (py1 > pBy-kMinSize) {py1 = sav1; r2 = sav2; py = pyold;}
    }
    if (pBot) {
      sav1 = py1;
      sav2 = r2;
      py1 += (py - pyold)/2;
      r2 += (py - pyold)/2;
      if (TMath::Abs(pyold-py)%2==1) impair++;
      if (py-pyold>0) sig=+1;
      else sig=-1;
      if (impair==2) { impair = 0; py1 += sig; r2 += sig;}
      if (py1 < pTy+kMinSize) {py1 = sav1; r2 = sav2; py = pyold;}
    }
    if (pL) {
      sav1 = px1;
      sav2 = r1;
      px1 += (px - pxold)/2;
      r1 -= (px - pxold)/2;
      if (TMath::Abs(pxold-px)%2==1) impair++;
      if (px-pxold>0) sig=+1;
      else sig=-1;
      if (impair==2) { impair = 0; px1 += sig; r1 -= sig;}
      if (px1 > pRx-kMinSize) {px1 = sav1; r1 = sav2; px = pxold;}
    }
    if (pR) {
      sav1 = px1;
      sav2 = r1;
      px1 += (px - pxold)/2;
      r1 += (px - pxold)/2;
      if (TMath::Abs(pxold-px)%2==1) impair++;
      if (px-pxold>0) sig=+1;
      else sig=-1;
      if (impair==2) { impair = 0; px1 += sig; r1 += sig;}
      if (px1 < pLx+kMinSize) {px1 = sav1; r1 = sav2; px = pxold;}
    }
    if (pTop || pBot || pL || pR) {
      gVirtualX->SetLineColor(-1);
      TAttLine::Modify();
      dphi = (fPhimax-fPhimin)*kPI/(180*np);
      ct   = TMath::Cos(kPI*fTheta/180);
      st   = TMath::Sin(kPI*fTheta/180);
      for (i=0;i<np;i++) {
	angle = fPhimin*kPI/180 + Double_t(i)*dphi;
	dx    = r1*TMath::Cos(angle);
	dy    = r2*TMath::Sin(angle);
	x[i]  = px1 + Int_t(dx*ct - dy*st);
	y[i]  = py1 + Int_t(dx*st + dy*ct);
      }
      if (fPhimax-fPhimin >= 360 ) {
	x[np] = x[0];
	y[np] = y[0];
	npe = np;
      } else {
	x[np]   = px1;
	y[np]   = py1;
	x[np+1] = x[0];
	y[np+1] = y[0];
	npe = np + 1;
      }
      for (i=0;i<npe;i++) gVirtualX->DrawLine(x[i], y[i], x[i+1], y[i+1]);
    }
    if (pINSIDE) {
      dpx  = px-pxold;  dpy = py-pyold;
      px1 += dpx; py1 += dpy;
      for (i=0;i<=npe;i++) { x[i] += dpx; y[i] += dpy;}
      for (i=0;i<npe;i++) gVirtualX->DrawLine(x[i], y[i], x[i+1], y[i+1]);
    }
    pTx = pBx = px1;
    pRx = px1+r1;
    pLx = px1-r1;
    pRy = pLy = py1;
    pTy = py1-r2;
    pBy = py1+r2;
    gVirtualX->DrawLine(pRx+4, py1+4, pRx-4, py1+4);
    gVirtualX->DrawLine(pRx-4, py1+4, pRx-4, py1-4);
    gVirtualX->DrawLine(pRx-4, py1-4, pRx+4, py1-4);
    gVirtualX->DrawLine(pRx+4, py1-4, pRx+4, py1+4);
    gVirtualX->DrawLine(pLx+4, py1+4, pLx-4, py1+4);
    gVirtualX->DrawLine(pLx-4, py1+4, pLx-4, py1-4);
    gVirtualX->DrawLine(pLx-4, py1-4, pLx+4, py1-4);
    gVirtualX->DrawLine(pLx+4, py1-4, pLx+4, py1+4);
    gVirtualX->DrawLine(px1+4, pBy+4, px1-4, pBy+4);
    gVirtualX->DrawLine(px1-4, pBy+4, px1-4, pBy-4);
    gVirtualX->DrawLine(px1-4, pBy-4, px1+4, pBy-4);
    gVirtualX->DrawLine(px1+4, pBy-4, px1+4, pBy+4);
    gVirtualX->DrawLine(px1+4, pTy+4, px1-4, pTy+4);
    gVirtualX->DrawLine(px1-4, pTy+4, px1-4, pTy-4);
    gVirtualX->DrawLine(px1-4, pTy-4, px1+4, pTy-4);
    gVirtualX->DrawLine(px1+4, pTy-4, px1+4, pTy+4);
    pxold = px;
    pyold = py;
    break;

  case kButton1Up:
    if (gROOT->IsEscaped()) {
      gROOT->SetEscape(kFALSE);
      break;
    }

    fX1 = gPad->AbsPixeltoX(px1);
    fY1 = gPad->AbsPixeltoY(py1);
    fBy = gPad->AbsPixeltoY(py1+r2);
    fTy = gPad->AbsPixeltoY(py1-r2);
    fLx = gPad->AbsPixeltoX(px1+r1);
    fRx = gPad->AbsPixeltoX(px1-r1);
    fR1 = TMath::Abs(fRx-fLx)/2;
    fR2 = TMath::Abs(fTy-fBy)/2;
    gPad->Modified(kTRUE);
    gVirtualX->SetLineColor(-1);
  }
}

void TArcArrow::Draw(Option_t *option)
{
   AppendPad(option);
}

void TArcArrow::DrawArcArrow(Double_t x1, Double_t y1,
			     Double_t r1, Double_t r2,
			     Double_t phimin, Double_t phimax,
			     Double_t theta, Float_t arrowsize,
			     Option_t *option)
{
  TArcArrow *newarcarrow = new TArcArrow(x1, y1, r1, r2, 
					 phimin, phimax,
					 theta, arrowsize,
					 option);
  TAttLine::Copy(*newarcarrow);
  TAttFill::Copy(*newarcarrow);
  newarcarrow->SetBit(kCanDelete);
  newarcarrow->AppendPad(option);
}

void TArcArrow::Paint(Option_t *option)
{
  Option_t *opt;
  if (option && strlen(option)) opt = option;
  else                          opt = (char*)GetOption();
  
  PaintArcArrow(fX1,fY1,fR1,fR2,fPhimin,fPhimax,fTheta, fArrowSize, opt);
}

void TArcArrow::PaintArcArrow(Double_t x1, Double_t y1, Double_t r1, Double_t r2,
			      Double_t phimin, Double_t phimax, Double_t theta,
			      Float_t arrowsize, Option_t *option)
{
  // Draw this ellipse with new coordinates.
  
  const Int_t np = 200;
  static Double_t x[np+3], y[np+3];
  TAttLine::Modify();  //Change line attributes only if necessary
  TAttFill::Modify();  //Change fill attributes only if necessary

  Double_t phi1 = TMath::Min(phimin,phimax);
  Double_t phi2 = TMath::Max(phimin,phimax);

  //set number of points approximatively proportional to the ellipse circumference
  Double_t circ = kPI*(r1+r2)*(phi2-phi1)/360;
  Int_t n = (Int_t)(np*circ/((gPad->GetX2()-gPad->GetX1())+(gPad->GetY2()-gPad->GetY1())));
  if (n < 8) n= 8;
  if (n > np) n = np;
  Double_t angle,dx,dy;
  Double_t dphi = (phi2-phi1)*kPI/(180*n);
  Double_t ct   = TMath::Cos(kPI*theta/180);
  Double_t st   = TMath::Sin(kPI*theta/180);
  Double_t length = 0;
  for (Int_t i=0;i<=n+2;i++) {
    angle = phi1*kPI/180 + Double_t(i-1)*dphi;
    dx    = r1*TMath::Cos(angle);
    dy    = r2*TMath::Sin(angle);
    x[i]  = gPad->XtoPad(x1 + dx*ct - dy*st);
    y[i]  = gPad->YtoPad(y1 + dx*st + dy*ct);
    
    if (i>0) length += TMath::Sqrt((x[i]-x[i-1])*(x[i]-x[i-1])+(y[i]-y[i-1])*(y[i]-y[i-1]));
  }
  TString opt = option;
  opt.ToLower();

  if (phi2-phi1 >= 360 ) {
    if (GetLineStyle()) gPad->PaintPolyLine(n+2, &x[1], &y[1]);
  } else {
    if (GetLineStyle()) gPad->PaintPolyLine(n+2, &x[1], &y[1]);
  }
  
  Double_t x1n, y1n, x2n, y2n;
  x1n = x[1];
  x2n = x[n+2];
  y1n = y[1];
  y2n = y[n+2];

  Double_t rSize  = 0.7*arrowsize*length;
  Double_t dSize  = rSize*TMath::Sin(TMath::Pi()*fAngle/180.);
  Double_t cosT1   = 1;
  Double_t sinT1   = 0;
  Double_t cosT2   = 1;
  Double_t sinT2   = 0;
  Double_t A1 = 0, A2 = 0;

  if (length>0) {
    A1 = TMath::ATan2(y[0]-y[2], x[0]-x[2]);
    A2 = TMath::ATan2(y[n+1]-y[n-1], x[n+1]-x[n-1]);
  
    cosT1 = TMath::Cos(A1);
    sinT1 = TMath::Sin(A1);
    cosT2 = TMath::Cos(A2);
    sinT2 = TMath::Sin(A2);
  }

  Double_t x1ar[4], y1ar[4];
  Double_t x2ar[4], y2ar[4];
  
  if (opt.BeginsWith("|-")) {

    x1ar[0] = x1n-sinT1*dSize;
    y1ar[0] = y1n+cosT1*dSize;
    x1ar[1] = x1n+sinT1*dSize;
    y1ar[1] = y1n-cosT1*dSize;
    
    gPad->PaintLine(x1ar[0],y1ar[0],x1ar[1],y1ar[1]);
    opt(0) = ' ';
  }
  
  if (opt.EndsWith("-|")) {
    
    x1ar[0] = x2n-sinT2*dSize;
    y1ar[0] = y2n+cosT2*dSize;
    x1ar[1] = x2n+sinT2*dSize;
    y1ar[1] = y2n-cosT2*dSize;
    
    gPad->PaintLine(x1ar[0],y1ar[0],x1ar[1],y1ar[1]);
    opt(opt.Length()-1) = ' ';
  }
  
  Double_t x1h = x1n;
  Double_t y1h = y1n;
  Double_t x2h = x2n;
  Double_t y2h = y2n;
  if (opt.Contains("->-") || opt.Contains("-|>-")) {
    x2h = x[n/2+2];
    y2h = y[n/2+2];

    A2 = TMath::ATan2(y[n/2+1+2]-y[n/2-1+2], x[n/2+1+2]-x[n/2-1+2]);
    cosT2 = TMath::Cos(A2);
    sinT2 = TMath::Sin(A2);
  }
  if (opt.Contains("-<-") || opt.Contains("-<|-")) {
    x2h = x[n/2-5];
    y2h = y[n/2-5];
    
    A1 = TMath::ATan2(y[n/2+1-5]-y[n/2-1-5], x[n/2+1-5]-x[n/2-1-5]);
    cosT1 = TMath::Cos(A1);
    sinT1 = TMath::Sin(A1);
  }
  
  if (opt.Contains(">")) {
    x2ar[0] = x2h - rSize*cosT2 - sinT2*dSize;
    y2ar[0] = y2h - rSize*sinT2 + cosT2*dSize;
    x2ar[1] = x2h;
    y2ar[1] = y2h;
    x2ar[2] = x2h - rSize*cosT2 + sinT2*dSize;
    y2ar[2] = y2h - rSize*sinT2 - cosT2*dSize;
    x2ar[3] = x2ar[0];
    y2ar[3] = y2ar[0];

    if (opt.Contains("|>")) {
      if (GetFillColor()) {
	gPad->PaintFillArea(3,x2ar,y2ar);
	gPad->PaintPolyLine(4,x2ar,y2ar);
      } else {
	gPad->PaintPolyLine(4,x2ar,y2ar);
      }
    } else {
      gPad->PaintPolyLine(3,x2ar,y2ar);
    }
  }

  if (opt.Contains("<")) {
    x1ar[0] = x1h - rSize*cosT1 - sinT1*dSize;
    y1ar[0] = y1h - rSize*sinT1 + cosT1*dSize;
    x1ar[1] = x1h;
    y1ar[1] = y1h;
    x1ar[2] = x1h - rSize*cosT1 + sinT1*dSize;
    y1ar[2] = y1h - rSize*sinT1 - cosT1*dSize;
    x1ar[3] = x1ar[0];
    y1ar[3] = y1ar[0];
    
    if (opt.Contains("<|")) {
      if (GetFillColor()) {
	gPad->PaintFillArea(3,x1ar,y1ar);
	gPad->PaintPolyLine(4,x1ar,y1ar);
      } else {
	gPad->PaintPolyLine(4,x1ar,y1ar);
      }
    } else {
      gPad->PaintPolyLine(3, x1ar,y1ar);
    }
  }
}
