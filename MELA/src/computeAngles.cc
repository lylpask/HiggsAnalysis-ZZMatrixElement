/** \file
 *
 *  MELA - cf. http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/sbologne/MELAproject/
 *
 *  $Date: 2013/01/31 18:10:04 $
 *  $Revision: 1.7 $
 */

#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <utility>
#include <algorithm>
#include <ZZMatrixElement/MELA/interface/computeAngles.h>
#include "TMath.h"
#include "TLorentzRotation.h"


using namespace std;


namespace mela{
  bool forbidMassiveLeptons = false;
}


/***** Leptons *****/

void mela::applyLeptonMassCorrection(bool flag){ mela::forbidMassiveLeptons = flag; }

void mela::constrainedRemoveLeptonMass(TLorentzVector& p1, TLorentzVector& p2){
  TLorentzVector nullFourVector(0, 0, 0, 0);
  if (p1==nullFourVector || p2==nullFourVector) return;

  TLorentzVector pZ_old = p1 + p2;
  TLorentzVector p1_new = p1;
  TLorentzVector p2_new = p2;

  double energy = p1.T();
  double mom = p1.P();
  double mom_new = (energy + mom) / 2.;
  if (mom != 0){
    double x_new = p1.X()*mom_new / mom;
    double y_new = p1.Y()*mom_new / mom;
    double z_new = p1.Z()*mom_new / mom;
    p1_new.SetXYZT(x_new, y_new, z_new, mom_new);
  }
  energy = p2.T();
  mom = p2.P();
  mom_new = (energy + mom) / 2.;
  if (mom != 0){
    double x_new = p2.X()*mom_new / mom;
    double y_new = p2.Y()*mom_new / mom;
    double z_new = p2.Z()*mom_new / mom;
    p2_new.SetXYZT(x_new, y_new, z_new, mom_new);
  }
  TLorentzVector pZ_new = p1_new + p2_new;
  double delta_mZ=pZ_old.M()-pZ_new.M();

  p1_new.Boost(-(pZ_new.BoostVector()));
  mom = p1_new.T();
  mom_new = mom + delta_mZ / 2.;
  if (mom!=0) p1_new *= mom_new/mom;
  p1_new.Boost((pZ_old.BoostVector()));

  p2_new.Boost(-(pZ_new.BoostVector()));
  mom = p2_new.T();
  mom_new = mom + delta_mZ / 2.;
  if (mom!=0) p2_new *= mom_new/mom;
  p2_new.Boost((pZ_old.BoostVector()));

  p1=p1_new; p2=p2_new;
}

void mela::computeAngles(
  TLorentzVector p4M11, int Z1_lept1Id,
  TLorentzVector p4M12, int Z1_lept2Id,
  TLorentzVector p4M21, int Z2_lept1Id,
  TLorentzVector p4M22, int Z2_lept2Id,
  float& costhetastar,
  float& costheta1,
  float& costheta2,
  float& Phi,
  float& Phi1
  ){
  /*
  cout << "dec begin" << endl;
  cout << "Entry:\n";
  cout << p4M11.X() << "\t" << p4M11.Y() << "\t" << p4M11.Z() << "\t" << p4M11.T() << endl;
  cout << p4M12.X() << "\t" << p4M12.Y() << "\t" << p4M12.Z() << "\t" << p4M12.T() << endl;
  cout << p4M21.X() << "\t" << p4M21.Y() << "\t" << p4M21.Z() << "\t" << p4M21.T() << endl;
  cout << p4M22.X() << "\t" << p4M22.Y() << "\t" << p4M22.Z() << "\t" << p4M22.T() << endl;
  */
  if (mela::forbidMassiveLeptons){
    if (!(fabs(Z1_lept1Id)==23 || fabs(Z1_lept1Id)==24 || fabs(Z1_lept1Id)==21 || fabs(Z1_lept1Id)==22 || fabs(Z1_lept1Id)==25 || fabs(Z1_lept2Id)==23 || fabs(Z1_lept2Id)==24 || fabs(Z1_lept2Id)==21 || fabs(Z1_lept2Id)==22 || fabs(Z1_lept2Id)==25)) mela::constrainedRemoveLeptonMass(p4M11, p4M12);
    if (!(fabs(Z2_lept1Id)==23 || fabs(Z2_lept1Id)==24 || fabs(Z2_lept1Id)==21 || fabs(Z2_lept1Id)==22 || fabs(Z2_lept1Id)==25 || fabs(Z2_lept2Id)==23 || fabs(Z2_lept2Id)==24 || fabs(Z2_lept2Id)==21 || fabs(Z2_lept2Id)==22 || fabs(Z2_lept2Id)==25)) mela::constrainedRemoveLeptonMass(p4M21, p4M22);
  }
  /*
  cout << "a.b.:\n";
  cout << p4M11.X() << "\t" << p4M11.Y() << "\t" << p4M11.Z() << "\t" << p4M11.T() << endl;
  cout << p4M12.X() << "\t" << p4M12.Y() << "\t" << p4M12.Z() << "\t" << p4M12.T() << endl;
  cout << p4M21.X() << "\t" << p4M21.Y() << "\t" << p4M21.Z() << "\t" << p4M21.T() << endl;
  cout << p4M22.X() << "\t" << p4M22.Y() << "\t" << p4M22.Z() << "\t" << p4M22.T() << endl;
  */

  //build Z 4-vectors
  TLorentzVector p4Z1 = p4M11 + p4M12;
  TLorentzVector p4Z2 = p4M21 + p4M22;

  // Sort Z1 leptons so that:
  if (
    (Z1_lept1Id*Z1_lept2Id<0 && Z1_lept1Id<0) // for OS pairs: lep1 must be the negative one
    ||
    ((Z1_lept1Id*Z1_lept2Id>0 || (Z1_lept1Id==0 && Z1_lept2Id==0)) && p4M11.Phi()<=p4M12.Phi()) //for SS pairs: use random deterministic convention
    ) swap(p4M11, p4M12);
  // Same for Z2 leptons
  if (
    (Z2_lept1Id*Z2_lept2Id<0 && Z2_lept1Id<0)
    ||
    ((Z2_lept1Id*Z2_lept2Id>0 || (Z2_lept1Id==0 && Z2_lept2Id==0)) && p4M21.Phi()<=p4M22.Phi())
    ) swap(p4M21, p4M22);

  // BEGIN THE CALCULATION

  // build H 4-vectors
  TLorentzVector p4H = p4Z1 + p4Z2;

  // -----------------------------------

  //// costhetastar
  TVector3 boostX = -(p4H.BoostVector());
  TLorentzVector thep4Z1inXFrame(p4Z1);
  TLorentzVector thep4Z2inXFrame(p4Z2);
  thep4Z1inXFrame.Boost(boostX);
  thep4Z2inXFrame.Boost(boostX);
  TVector3 theZ1X_p3 = TVector3(thep4Z1inXFrame.X(), thep4Z1inXFrame.Y(), thep4Z1inXFrame.Z());
  TVector3 theZ2X_p3 = TVector3(thep4Z2inXFrame.X(), thep4Z2inXFrame.Y(), thep4Z2inXFrame.Z());
  costhetastar = theZ1X_p3.CosTheta();

  TVector3 boostV1(0, 0, 0);
  TVector3 boostV2(0, 0, 0);
  //// --------------------------- costheta1
  if (!(fabs(Z1_lept1Id)==21 || fabs(Z1_lept1Id)==22 || fabs(Z1_lept2Id)==21 || fabs(Z1_lept2Id)==22)){
    boostV1 = -(p4Z1.BoostVector());
    if (boostV1.Mag()>=1.) {
      cout << "Warning: Mela::computeAngles: Z1 boost with beta=1, scaling down" << endl;
      boostV1*=0.9999/boostV1.Mag();
    }
    TLorentzVector p4M11_BV1(p4M11);
    TLorentzVector p4M12_BV1(p4M12);
    TLorentzVector p4M21_BV1(p4M21);
    TLorentzVector p4M22_BV1(p4M22);
    p4M11_BV1.Boost(boostV1);
    p4M12_BV1.Boost(boostV1);
    p4M21_BV1.Boost(boostV1);
    p4M22_BV1.Boost(boostV1);

    TLorentzVector p4V2_BV1 = p4M21_BV1 + p4M22_BV1;
    //// costheta1
    costheta1 = -p4V2_BV1.Vect().Unit().Dot(p4M11_BV1.Vect().Unit());
  }
  else costheta1 = 0;

  //// --------------------------- costheta2
  if (!(fabs(Z2_lept1Id)==21 || fabs(Z2_lept1Id)==22 || fabs(Z2_lept2Id)==21 || fabs(Z2_lept2Id)==22)){
    boostV2 = -(p4Z2.BoostVector());
    if (boostV2.Mag()>=1.) {
      cout << "Warning: Mela::computeAngles: Z2 boost with beta=1, scaling down" << endl;
      boostV2*=0.9999/boostV2.Mag();
    }
    TLorentzVector p4M11_BV2(p4M11);
    TLorentzVector p4M12_BV2(p4M12);
    TLorentzVector p4M21_BV2(p4M21);
    TLorentzVector p4M22_BV2(p4M22);
    p4M11_BV2.Boost(boostV2);
    p4M12_BV2.Boost(boostV2);
    p4M21_BV2.Boost(boostV2);
    p4M22_BV2.Boost(boostV2);

    TLorentzVector p4V1_BV2 = p4M11_BV2 + p4M12_BV2;
    //// costheta2
    costheta2 = -p4V1_BV2.Vect().Unit().Dot(p4M21_BV2.Vect().Unit());
  }
  else costheta2 = 0;

  //// --------------------------- Phi and Phi1 (old phistar1 - azimuthal production angle)
  TLorentzVector p4M11_BX(p4M11);
  TLorentzVector p4M12_BX(p4M12);
  TLorentzVector p4M21_BX(p4M21);
  TLorentzVector p4M22_BX(p4M22);

  p4M11_BX.Boost(boostX);
  p4M12_BX.Boost(boostX);
  p4M21_BX.Boost(boostX);
  p4M22_BX.Boost(boostX);
  TLorentzVector p4V1_BX = p4M11_BX + p4M12_BX;

  TVector3 beamAxis(0, 0, 1);
  TVector3 p3V1_BX = p4V1_BX.Vect().Unit();
  TVector3 normal1_BX = (p4M11_BX.Vect().Cross(p4M12_BX.Vect())).Unit();
  TVector3 normal2_BX = (p4M21_BX.Vect().Cross(p4M22_BX.Vect())).Unit();
  TVector3 normalSC_BX = (beamAxis.Cross(p3V1_BX)).Unit();


  //// Phi
  float tmpSgnPhi = p3V1_BX.Dot(normal1_BX.Cross(normal2_BX));
  float sgnPhi = 0;
  if (fabs(tmpSgnPhi)>0.) sgnPhi = tmpSgnPhi/fabs(tmpSgnPhi);
  float dot_BX12 = normal1_BX.Dot(normal2_BX);
  if (fabs(dot_BX12)>=1.) dot_BX12 *= 1./fabs(dot_BX12);
  Phi = sgnPhi * acos(-1.*dot_BX12);


  //// Phi1
  float tmpSgnPhi1 = p3V1_BX.Dot(normal1_BX.Cross(normalSC_BX));
  float sgnPhi1 = 0;
  if (fabs(tmpSgnPhi1)>0.) sgnPhi1 = tmpSgnPhi1/fabs(tmpSgnPhi1);
  float dot_BX1SC = normal1_BX.Dot(normalSC_BX);
  if (fabs(dot_BX1SC)>=1.) dot_BX1SC *= 1./fabs(dot_BX1SC);
  Phi1 = sgnPhi1 * acos(dot_BX1SC);

  if (isnan(costhetastar) || isnan(costheta1) || isnan(costheta2) || isnan(Phi) || isnan(Phi1)){
    cout << "WARNING: NaN in computeAngles: "
      << costhetastar << " "
      << costheta1  << " "
      << costheta2  << " "
      << Phi  << " "
      << Phi1  << " " << endl;
    cout << "   boostV1: " <<boostV1.Pt() << " " << boostV1.Eta() << " " << boostV1.Phi() << " " << boostV1.Mag() << endl;
    cout << "   boostV2: " <<boostV2.Pt() << " " << boostV2.Eta() << " " << boostV2.Phi() << " " << boostV2.Mag() << endl;
  }
}

void mela::computeAnglesCS(
  TLorentzVector p4M11, int Z1_lept1Id,
  TLorentzVector p4M12, int Z1_lept2Id,
  TLorentzVector p4M21, int Z2_lept1Id,
  TLorentzVector p4M22, int Z2_lept2Id,
  float pbeam,
  float& costhetastar,
  float& costheta1,
  float& costheta2,
  float& Phi,
  float& Phi1
  ){
  if (mela::forbidMassiveLeptons){
    if (!(fabs(Z1_lept1Id)==23 || fabs(Z1_lept1Id)==24 || fabs(Z1_lept1Id)==21 || fabs(Z1_lept1Id)==22 || fabs(Z1_lept1Id)==25 || fabs(Z1_lept2Id)==23 || fabs(Z1_lept2Id)==24 || fabs(Z1_lept2Id)==21 || fabs(Z1_lept2Id)==22 || fabs(Z1_lept2Id)==25)) mela::constrainedRemoveLeptonMass(p4M11, p4M12);
    if (!(fabs(Z2_lept1Id)==23 || fabs(Z2_lept1Id)==24 || fabs(Z2_lept1Id)==21 || fabs(Z2_lept1Id)==22 || fabs(Z2_lept1Id)==25 || fabs(Z2_lept2Id)==23 || fabs(Z2_lept2Id)==24 || fabs(Z2_lept2Id)==21 || fabs(Z2_lept2Id)==22 || fabs(Z2_lept2Id)==25)) mela::constrainedRemoveLeptonMass(p4M21, p4M22);
  }

  TVector3 LabXaxis(1.0, 0.0, 0.0);
  TVector3 LabYaxis(0.0, 1.0, 0.0);
  TVector3 LabZaxis(0.0, 0.0, 1.0);

  float Mprot = 0.938;
  float Ebeam = sqrt(pbeam*pbeam + Mprot*Mprot);

  TLorentzVector targ(0., 0., -pbeam, Ebeam);
  TLorentzVector beam(0., 0., pbeam, Ebeam);

  //build Z 4-vectors
  TLorentzVector p4Z1 = p4M11 + p4M12;
  TLorentzVector p4Z2 = p4M21 + p4M22;

  // Sort Z1 leptons so that:
  if (
    (Z1_lept1Id*Z1_lept2Id<0 && Z1_lept1Id<0) // for OS pairs: lep1 must be the negative one
    ||
    ((Z1_lept1Id*Z1_lept2Id>0 || (Z1_lept1Id==0 && Z1_lept2Id==0)) && p4M11.Phi()<=p4M12.Phi()) //for SS pairs: use random deterministic convention
    ) swap(p4M11, p4M12);
  // Same for Z2 leptons
  if (
    (Z2_lept1Id*Z2_lept2Id<0 && Z2_lept1Id<0)
    ||
    ((Z2_lept1Id*Z2_lept2Id>0 || (Z2_lept1Id==0 && Z2_lept2Id==0)) && p4M21.Phi()<=p4M22.Phi())
    ) swap(p4M21, p4M22);


  //build H 4-vectors
  TLorentzVector p4H = p4Z1 + p4Z2;
  TVector3 boostX = -(p4H.BoostVector());

  /////////////////////////////
  // Collin-Sopper calculation:
  // in the CS frame, the z-axis is along the bisectrice of one beam and the opposite of the other beam,
  // after their boost in X
  ///////////////////////////////
  // Rotation for the CS Frame

  TRotation rotationCS;

  TLorentzVector beaminX(beam);
  TLorentzVector targinX(targ);
  targinX.Boost(boostX);
  beaminX.Boost(boostX);

  //Bisectrice: sum of unit vectors (remember: you need to invert one beam vector)
  TVector3 beam_targ_bisecinX((beaminX.Vect().Unit() - targinX.Vect().Unit()).Unit());

  // Define a rotationCS Matrix, with Z along the bisectric, 
  TVector3 newZaxisCS(beam_targ_bisecinX.Unit());
  TVector3 newYaxisCS(beaminX.Vect().Unit().Cross(newZaxisCS).Unit());
  TVector3 newXaxisCS(newYaxisCS.Unit().Cross(newZaxisCS).Unit());
  rotationCS.RotateAxes(newXaxisCS, newYaxisCS, newZaxisCS);
  rotationCS.Invert();

  //// costhetastar
  TLorentzVector thep4Z1inXFrame_rotCS(p4Z1);
  TLorentzVector thep4Z2inXFrame_rotCS(p4Z2);
  thep4Z1inXFrame_rotCS.Transform(rotationCS);
  thep4Z2inXFrame_rotCS.Transform(rotationCS);
  thep4Z1inXFrame_rotCS.Boost(boostX);
  thep4Z2inXFrame_rotCS.Boost(boostX);
  TVector3 theZ1XrotCS_p3 = TVector3(thep4Z1inXFrame_rotCS.X(), thep4Z1inXFrame_rotCS.Y(), thep4Z1inXFrame_rotCS.Z());
  costhetastar = theZ1XrotCS_p3.CosTheta();

  //// --------------------------- Phi and Phi1 (old phistar1 - azimuthal production angle)
  //    TVector3 boostX = -(thep4H.BoostVector());
  TLorentzVector p4M11_BX_rotCS(p4M11);
  TLorentzVector p4M12_BX_rotCS(p4M12);
  TLorentzVector p4M21_BX_rotCS(p4M21);
  TLorentzVector p4M22_BX_rotCS(p4M22);
  p4M11_BX_rotCS.Transform(rotationCS);
  p4M12_BX_rotCS.Transform(rotationCS);
  p4M21_BX_rotCS.Transform(rotationCS);
  p4M22_BX_rotCS.Transform(rotationCS);
  p4M11_BX_rotCS.Boost(boostX);
  p4M12_BX_rotCS.Boost(boostX);
  p4M21_BX_rotCS.Boost(boostX);
  p4M22_BX_rotCS.Boost(boostX);

  TLorentzVector p4Z1_BX_rotCS = p4M11_BX_rotCS + p4M12_BX_rotCS;
  TVector3 p3V1_BX_rotCS = (p4Z1_BX_rotCS.Vect()).Unit();
  TVector3 beamAxis(0, 0, 1);
  TVector3 normal1_BX_rotCS = (p4M11_BX_rotCS.Vect().Cross(p4M12_BX_rotCS.Vect())).Unit();
  TVector3 normal2_BX_rotCS = (p4M21_BX_rotCS.Vect().Cross(p4M22_BX_rotCS.Vect())).Unit();
  TVector3 normalSC_BX_rotCS = (beamAxis.Cross(p3V1_BX_rotCS)).Unit();

  //// Phi
  float tmpSgnPhi = p3V1_BX_rotCS.Vect().Dot(normal1_BX_rotCS.Cross(normal2_BX_rotCS));
  float sgnPhi = 0;
  if (fabs(tmpSgnPhi)>0.) sgnPhi = tmpSgnPhi/fabs(tmpSgnPhi);
  float dot_BX12 = normal1_BX_rotCS.Dot(normal2_BX_rotCS);
  if (fabs(dot_BX12)>=1.) dot_BX12 *= 1./fabs(dot_BX12);
  Phi = sgnPhi * acos(-1.*dot_BX12);

  //// Phi1
  float tmpSgnPhi1 = p3V1_BX_rotCS.Vect().Dot(normal1_BX_rotCS.Cross(normalSC_BX_rotCS));
  float sgnPhi1 = 0;
  if (fabs(tmpSgnPhi1)>0.) sgnPhi1 = tmpSgnPhi1/fabs(tmpSgnPhi1);
  float dot_BX1SC = normal1_BX_rotCS.Dot(normalSC_BX_rotCS);
  if (fabs(dot_BX1SC)>=1.) dot_BX1SC *= 1./fabs(dot_BX1SC);
  Phi1 = sgnPhi1 * acos(dot_BX1SC);

  //// --------------------------- costheta1
  if (!(fabs(Z1_lept1Id)==21 || fabs(Z1_lept1Id)==22 || fabs(Z1_lept2Id)==21 || fabs(Z1_lept2Id)==22)){
    //define Z1 rotation 
    TRotation rotationZ1;
    TVector3 newZaxisZ1(thep4Z1inXFrame_rotCS.Vect().Unit());
    TVector3 newXaxisZ1(newYaxisCS.Cross(newZaxisZ1).Unit());
    TVector3 newYaxisZ1(newZaxisZ1.Cross(newXaxisZ1).Unit());
    rotationZ1.RotateAxes(newXaxisZ1, newYaxisZ1, newZaxisZ1);
    rotationZ1.Invert();

    TLorentzVector thep4Z1inXFrame_rotCS_rotZ1(thep4Z1inXFrame_rotCS);
    thep4Z1inXFrame_rotCS_rotZ1.Transform(rotationZ1);
    TVector3 boostZ1inX_rotCS_rotZ1= -(thep4Z1inXFrame_rotCS_rotZ1.BoostVector());

    TLorentzVector p4M11_BX_rotCS_rotZ1(p4M11_BX_rotCS);
    TLorentzVector p4M12_BX_rotCS_rotZ1(p4M12_BX_rotCS);
    TLorentzVector p4M21_BX_rotCS_rotZ1(p4M21_BX_rotCS);
    TLorentzVector p4M22_BX_rotCS_rotZ1(p4M22_BX_rotCS);
    p4M11_BX_rotCS_rotZ1.Transform(rotationZ1);
    p4M12_BX_rotCS_rotZ1.Transform(rotationZ1);
    p4M21_BX_rotCS_rotZ1.Transform(rotationZ1);
    p4M22_BX_rotCS_rotZ1.Transform(rotationZ1);
    p4M11_BX_rotCS_rotZ1.Boost(boostZ1inX_rotCS_rotZ1);
    p4M12_BX_rotCS_rotZ1.Boost(boostZ1inX_rotCS_rotZ1);
    p4M21_BX_rotCS_rotZ1.Boost(boostZ1inX_rotCS_rotZ1);
    p4M22_BX_rotCS_rotZ1.Boost(boostZ1inX_rotCS_rotZ1);

    TLorentzVector p4V2_BX_rotCS_rotZ1 = p4M21_BX_rotCS_rotZ1 + p4M22_BX_rotCS_rotZ1;
    //// costheta1
    costheta1 = -p4V2_BX_rotCS_rotZ1.Vect().Dot(p4M11_BX_rotCS_rotZ1.Vect())/p4V2_BX_rotCS_rotZ1.Vect().Mag()/p4M11_BX_rotCS_rotZ1.Vect().Mag();
  }
  else costheta1=0;

  //// --------------------------- costheta2
  if (!(fabs(Z2_lept1Id)==21 || fabs(Z2_lept1Id)==22 || fabs(Z2_lept2Id)==21 || fabs(Z2_lept2Id)==22)){
    //define Z2 rotation 
    TRotation rotationZ2;
    TVector3 newZaxisZ2(thep4Z2inXFrame_rotCS.Vect().Unit());
    TVector3 newXaxisZ2(newYaxisCS.Cross(newZaxisZ2).Unit());
    TVector3 newYaxisZ2(newZaxisZ2.Cross(newXaxisZ2).Unit());
    rotationZ2.RotateAxes(newXaxisZ2, newYaxisZ2, newZaxisZ2);
    rotationZ2.Invert();

    TLorentzVector thep4Z2inXFrame_rotCS_rotZ2(thep4Z2inXFrame_rotCS);
    thep4Z2inXFrame_rotCS_rotZ2.Transform(rotationZ2);
    TVector3 boostZ2inX_rotCS_rotZ2= -(thep4Z2inXFrame_rotCS_rotZ2.BoostVector());

    TLorentzVector p4M11_BX_rotCS_rotZ2(p4M11_BX_rotCS);
    TLorentzVector p4M12_BX_rotCS_rotZ2(p4M12_BX_rotCS);
    TLorentzVector p4M21_BX_rotCS_rotZ2(p4M21_BX_rotCS);
    TLorentzVector p4M22_BX_rotCS_rotZ2(p4M22_BX_rotCS);
    p4M11_BX_rotCS_rotZ2.Transform(rotationZ2);
    p4M12_BX_rotCS_rotZ2.Transform(rotationZ2);
    p4M21_BX_rotCS_rotZ2.Transform(rotationZ2);
    p4M22_BX_rotCS_rotZ2.Transform(rotationZ2);
    p4M11_BX_rotCS_rotZ2.Boost(boostZ2inX_rotCS_rotZ2);
    p4M12_BX_rotCS_rotZ2.Boost(boostZ2inX_rotCS_rotZ2);
    p4M21_BX_rotCS_rotZ2.Boost(boostZ2inX_rotCS_rotZ2);
    p4M22_BX_rotCS_rotZ2.Boost(boostZ2inX_rotCS_rotZ2);


    TLorentzVector p4V1_BX_rotCS_rotZ2= p4M11_BX_rotCS_rotZ2 + p4M12_BX_rotCS_rotZ2;
    //// costheta2
    costheta2 = -p4V1_BX_rotCS_rotZ2.Vect().Dot(p4M21_BX_rotCS_rotZ2.Vect())/p4V1_BX_rotCS_rotZ2.Vect().Mag()/p4M21_BX_rotCS_rotZ2.Vect().Mag();
  }
  else costheta2=0;

  if (isnan(costhetastar) || isnan(costheta1) || isnan(costheta2) || isnan(Phi) || isnan(Phi1)){
    cout << "WARNING: NaN in computeAngles: "
      << costhetastar << " "
      << costheta1  << " "
      << costheta2  << " "
      << Phi  << " "
      << Phi1  << " " << endl;
  }
}

/***** Jets *****/

void mela::computeVBFangles(
  float& costhetastar,
  float& costheta1,
  float& costheta2,
  float& Phi,
  float& Phi1,
  float& Q2V1,
  float& Q2V2,
  TLorentzVector p4M11, int Z1_lept1Id,
  TLorentzVector p4M12, int Z1_lept2Id,
  TLorentzVector p4M21, int Z2_lept1Id,
  TLorentzVector p4M22, int Z2_lept2Id,
  TLorentzVector jet1, int jet1Id,
  TLorentzVector jet2, int jet2Id,
  TLorentzVector* injet1, int injet1Id, // Gen. partons in lab frame
  TLorentzVector* injet2, int injet2Id
  ){

  if (mela::forbidMassiveLeptons){
    if (!(fabs(Z1_lept1Id)==23 || fabs(Z1_lept1Id)==24 || fabs(Z1_lept1Id)==21 || fabs(Z1_lept1Id)==25 || fabs(Z1_lept2Id)==23 || fabs(Z1_lept2Id)==24 || fabs(Z1_lept2Id)==21 || fabs(Z1_lept2Id)==25)) mela::constrainedRemoveLeptonMass(p4M11, p4M12);
    if (!(fabs(Z2_lept1Id)==23 || fabs(Z2_lept1Id)==24 || fabs(Z2_lept1Id)==21 || fabs(Z2_lept1Id)==25 || fabs(Z2_lept2Id)==23 || fabs(Z2_lept2Id)==24 || fabs(Z2_lept2Id)==21 || fabs(Z2_lept2Id)==25)) mela::constrainedRemoveLeptonMass(p4M21, p4M22);
  }

  TLorentzVector jet1massless, jet2massless;
  mela::computeJetMassless(jet1, jet1massless);
  mela::computeJetMassless(jet2, jet2massless);

  //build Z 4-vectors
  TLorentzVector p4Z1 = p4M11 + p4M12;
  TLorentzVector p4Z2 = p4M21 + p4M22;
  TLorentzVector pH = p4Z1+p4Z2;
  //jet1 is defined as going forwards (bigger pz), jet2 going backwards (smaller pz)
  if (jet1massless.Z() < jet2massless.Z()) { swap(jet1massless, jet2massless); swap(jet1Id, jet2Id); }

  //Find the incoming partons - first boost so the pT(HJJ) = 0, then boost away the pz.
  //This preserves the z direction.  Then assume that the partons come in this z direction.
  //This is exactly correct at LO (since pT=0 anyway).
  //Then associate the one going forwards with jet1 and the one going backwards with jet2
  TLorentzRotation movingframe;
  TLorentzVector pHJJ = pH+jet1massless+jet2massless;
  TLorentzVector pHJJ_perp(pHJJ.X(), pHJJ.Y(), 0, pHJJ.T());
  movingframe.Boost(-pHJJ_perp.BoostVector());
  pHJJ.Boost(-pHJJ_perp.BoostVector());
  movingframe.Boost(-pHJJ.BoostVector());
  pHJJ.Boost(-pHJJ.BoostVector());   //make sure to boost HJJ AFTER boosting movingframe

  TLorentzVector P1(0, 0, pHJJ.T()/2, pHJJ.T()/2);
  TLorentzVector P2(0, 0, -pHJJ.T()/2, pHJJ.T()/2);
  // Transform incoming partons back to original frame
  P1.Transform(movingframe.Inverse());
  P2.Transform(movingframe.Inverse());
  //movingframe, HJJ, and HJJ_T will not be used anymore
  if (injet1!=0 && injet2!=0){ // Handle gen. partons if they are available
    if (fabs((*injet1+*injet2).P()-pHJJ.P())<pHJJ.P()*1e-4){
      P1=*injet1;
      P2=*injet2;
      if (P1.Z() < P2.Z()){
        swap(P1, P2);
        swap(injet1Id, injet2Id);
      }
      // In the case of gen. partons, check if the intermediates are a Z or a W.
      int diff1Id = jet1Id-injet1Id;
      int diff2Id = jet2Id-injet2Id;
      if (
        !(
        (diff1Id==0 && diff2Id==0 && !(injet1Id==21 || injet2Id==21)) // Two Z bosons
        ||
        ((fabs(diff1Id)==1 || fabs(diff1Id)==3 || fabs(diff1Id)==5) && (fabs(diff2Id)==1 || fabs(diff2Id)==3 || fabs(diff2Id)==5)) // Two W bosons, do not check W+ vs W-
        )
        ){
        int diff12Id = jet1Id-injet2Id;
        int diff21Id = jet2Id-injet1Id;
        if (
          ((diff12Id==0 || diff21Id==0) && !(injet1Id==21 || injet2Id==21)) // At least one Z boson
          ||
          ((fabs(diff12Id)==1 || fabs(diff12Id)==3 || fabs(diff12Id)==5) || (fabs(diff21Id)==1 || fabs(diff21Id)==3 || fabs(diff21Id)==5)) // At least one W boson
          ){
          swap(P1, P2);
          swap(injet1Id, injet2Id);
        }
      }
    }
  }

  TLorentzRotation ZZframe;
  ZZframe.Boost(-pH.BoostVector());
  P1.Transform(ZZframe);
  P2.Transform(ZZframe);
  p4Z1.Transform(ZZframe);
  p4Z2.Transform(ZZframe);
  jet1massless.Transform(ZZframe);
  jet2massless.Transform(ZZframe);

  TLorentzVector V1 = P1-jet1massless; // V1 = (-p12) - p11 = -Z1
  TLorentzVector V2 = P2-jet2massless; // V2 = (-p22) - p21 = -Z2
  Q2V1 = -V1.M2();
  Q2V2 = -V2.M2();

  costhetastar = -V1.Vect().Unit().Dot(p4Z2.Vect().Unit());
  costheta1 = -V1.Vect().Unit().Dot(jet1massless.Vect().Unit());
  costheta2 = -V2.Vect().Unit().Dot(jet2massless.Vect().Unit());

  TVector3 normvec1 = P1.Vect().Cross(jet1massless.Vect()).Unit(); // p11 x p12 = (-p12) x p11
  TVector3 normvec2 = P2.Vect().Cross(jet2massless.Vect()).Unit(); // p21 x p22 = (-p22) x p21
  TVector3 normvec3 = V1.Vect().Cross(p4Z2.Vect()).Unit(); // == z x Z1

  double cosPhi = normvec1.Dot(normvec2);
  double sgnPhi = normvec1.Cross(normvec2).Dot(-V1.Vect());
  double cosPhi1 = normvec1.Dot(normvec3);
  double sgnPhi1 = normvec1.Cross(normvec3).Dot(-V1.Vect());
  if (fabs(cosPhi)>1) cosPhi *= 1./fabs(cosPhi);
  if (fabs(cosPhi1)>1) cosPhi1 *= 1./fabs(cosPhi1);
  Phi = TMath::Sign(acos(-cosPhi), sgnPhi);            //TMath::Sign(a,b) = |a|*(b/|b|)
  Phi1 = TMath::Sign(acos(cosPhi1), sgnPhi1);
}

void mela::computeVHangles(
  float& costhetastar,
  float& costheta1,
  float& costheta2,
  float& Phi,
  float& Phi1,
  TLorentzVector p4M11, int Z1_lept1Id,
  TLorentzVector p4M12, int Z1_lept2Id,
  TLorentzVector p4M21, int Z2_lept1Id,
  TLorentzVector p4M22, int Z2_lept2Id,
  TLorentzVector jet1, int jet1Id,
  TLorentzVector jet2, int jet2Id,
  TLorentzVector* injet1, int injet1Id, // Gen. partons in lab frame
  TLorentzVector* injet2, int injet2Id
  ){
  TLorentzVector nullFourVector(0, 0, 0, 0);

  //cout << p4M11.X() << '\t' << p4M11.Y() << '\t' << p4M11.Z() << '\t' << p4M11.T() << endl;
  //cout << p4M12.X() << '\t' << p4M12.Y() << '\t' << p4M12.Z() << '\t' << p4M12.T() << endl;
  //cout << p4M21.X() << '\t' << p4M21.Y() << '\t' << p4M21.Z() << '\t' << p4M21.T() << endl;
  //cout << p4M22.X() << '\t' << p4M22.Y() << '\t' << p4M22.Z() << '\t' << p4M22.T() << endl;
  if (mela::forbidMassiveLeptons){
    if (!(fabs(Z1_lept1Id)==23 || fabs(Z1_lept1Id)==24 || fabs(Z1_lept1Id)==21 || fabs(Z1_lept1Id)==25 || fabs(Z1_lept2Id)==23 || fabs(Z1_lept2Id)==24 || fabs(Z1_lept2Id)==21 || fabs(Z1_lept2Id)==25)) mela::constrainedRemoveLeptonMass(p4M11, p4M12);
    if (!(fabs(Z2_lept1Id)==23 || fabs(Z2_lept1Id)==24 || fabs(Z2_lept1Id)==21 || fabs(Z2_lept1Id)==25 || fabs(Z2_lept2Id)==23 || fabs(Z2_lept2Id)==24 || fabs(Z2_lept2Id)==21 || fabs(Z2_lept2Id)==25)) mela::constrainedRemoveLeptonMass(p4M21, p4M22);
  }
  //cout << endl;
  //cout << p4M11.X() << '\t' << p4M11.Y() << '\t' << p4M11.Z() << '\t' << p4M11.T() << endl;
  //cout << p4M12.X() << '\t' << p4M12.Y() << '\t' << p4M12.Z() << '\t' << p4M12.T() << endl;
  //cout << p4M21.X() << '\t' << p4M21.Y() << '\t' << p4M21.Z() << '\t' << p4M21.T() << endl;
  //cout << p4M22.X() << '\t' << p4M22.Y() << '\t' << p4M22.Z() << '\t' << p4M22.T() << endl;

  // Build Z 4-vectors
  TLorentzVector p4Z1 = p4M11 + p4M12;
  TLorentzVector p4Z2 = p4M21 + p4M22;
  TLorentzVector pH = p4Z1 + p4Z2;

  TLorentzVector jet1massless, jet2massless;
  mela::computeJetMassless(jet1, jet1massless);
  mela::computeJetMassless(jet2, jet2massless);

  // Apply convention for outgoing particles
  if (
    (jet1Id*jet2Id<0 && jet1Id<0) // for OS pairs: jet1 must be the particle
    ||
    (jet1Id*jet2Id>0 && jet1massless.Phi()<=jet2massless.Phi()) // for SS pairs: use random deterministic convention
    ){
    swap(jet1massless, jet2massless);
    swap(jet1Id, jet2Id);
  }

  //Find the incoming partons - first boost so the pT(HJJ) = 0, then boost away the pz.
  //This preserves the z direction.  Then assume that the partons come in this z direction.
  //This is exactly correct at LO (since pT=0 anyway).
  //Then associate the one going forwards with jet1 and the one going backwards with jet2
  TLorentzRotation movingframe;
  TLorentzVector pHJJ = pH+jet1massless+jet2massless;
  TLorentzVector pHJJ_perp(pHJJ.X(), pHJJ.Y(), 0, pHJJ.T());
  movingframe.Boost(-pHJJ_perp.BoostVector());
  pHJJ.Boost(-pHJJ_perp.BoostVector());
  movingframe.Boost(-pHJJ.BoostVector());
  pHJJ.Boost(-pHJJ.BoostVector());   //make sure to boost HJJ AFTER boosting movingframe

  TLorentzVector P1(0, 0, -pHJJ.T()/2, pHJJ.T()/2);
  TLorentzVector P2(0, 0, pHJJ.T()/2, pHJJ.T()/2);
  // Transform incoming partons back to the original frame
  P1.Transform(movingframe.Inverse());
  P2.Transform(movingframe.Inverse());
  //movingframe, HJJ, and HJJ_T will not be used anymore
  if (injet1!=0 && injet2!=0){ // Handle gen. partons if they are available
    if (fabs((*injet1+*injet2).P()-pHJJ.P())<=pHJJ.P()*1e-4){
      P1=*injet1;
      P2=*injet2;
      // Apply convention for incoming (!) particles
      if (
        (injet1Id*injet2Id<0 && injet1Id>0) // for OS pairs: parton 2 must be the particle
        ||
        (injet1Id*injet2Id>0 && P1.Z()>=P2.Z()) //for SS pairs: use random deterministic convention
        ){
        swap(P1, P2);
        swap(injet1Id, injet2Id);
      }
    }
  }

  // Rotate every vector such that Z1 - Z2 axis is the "beam axis" analogue of decay
  TLorentzRotation ZZframe;
  TVector3 beamAxis(0, 0, 1);
  if (p4Z1==nullFourVector || p4Z2==nullFourVector){
    TVector3 pNewAxis = (p4Z2-p4Z1).Vect().Unit(); // Let Z2 be in the z direction so that once the direction of H is reversed, Z1 is in the z direction
    if (pNewAxis != nullFourVector.Vect()){
      TVector3 pNewAxisPerp = pNewAxis.Cross(beamAxis);
      ZZframe.Rotate(acos(pNewAxis.Dot(beamAxis)), pNewAxisPerp);

      P1.Transform(ZZframe);
      P2.Transform(ZZframe);
      jet1massless = -jet1massless;
      jet2massless = -jet2massless;
      jet1massless.Transform(ZZframe);
      jet2massless.Transform(ZZframe);
      jet1massless = -jet1massless;
      jet2massless = -jet2massless;
    }
  }
  else{
    ZZframe.Boost(-pH.BoostVector());
    p4Z1.Boost(-pH.BoostVector());
    p4Z2.Boost(-pH.BoostVector());
    TVector3 pNewAxis = (p4Z2-p4Z1).Vect().Unit(); // Let Z2 be in the z direction so that once the direction of H is reversed, Z1 is in the z direction
    TVector3 pNewAxisPerp = pNewAxis.Cross(beamAxis);
    ZZframe.Rotate(acos(pNewAxis.Dot(beamAxis)), pNewAxisPerp);
    P1.Transform(ZZframe);
    P2.Transform(ZZframe);
    jet1massless = -jet1massless;
    jet2massless = -jet2massless;
    jet1massless.Transform(ZZframe);
    jet2massless.Transform(ZZframe);
    jet1massless = -jet1massless;
    jet2massless = -jet2massless;
  }

  mela::computeAngles(
    -P1, 23, // Id is 23 to avoid an attempt to remove quark mass
    -P2, 0, // Id is 0 to avoid swapping
    jet1massless, 23,
    jet2massless, 0,
    costhetastar,
    costheta1,
    costheta2,
    Phi,
    Phi1
    );
}




void mela::computeJetMassless(TLorentzVector massiveJet, TLorentzVector& masslessJet){
  double energy, p3sq, ratio;
  energy = massiveJet.T();
  p3sq = massiveJet.P();
  ratio = (p3sq>0 ? (energy / p3sq) : 1);
  masslessJet.SetPxPyPzE(massiveJet.Px()*ratio, massiveJet.Py()*ratio, massiveJet.Pz()*ratio, energy);
}

void mela::computeFakeJet(TLorentzVector realJet, TLorentzVector others, TLorentzVector& fakeJet){
  TLorentzVector masslessRealJet(0, 0, 0, 0);
  mela::computeJetMassless(realJet, masslessRealJet);

  fakeJet = others + masslessRealJet;
  fakeJet.SetVect(-fakeJet.Vect());
  fakeJet.SetE(fakeJet.P());
}

