#include <ZZMatrixElement/MELA/interface/MelaPConstant.h>
#include <iostream>
#include <vector>
#include <cmath>
#include "TROOT.h"
#include "TFile.h"
#include "TString.h"


using namespace std;


MelaPConstant::MelaPConstant(
  TVar::MatrixElement me_,
  TVar::Production prod_,
  TVar::Process proc_,
  const char* path,
  const char* spname
  ) :
  processME(me_),
  processProd(prod_),
  processProc(proc_),
  fcnLow(0),
  fcnHigh(0),
  fcnMid(0)
{
  GetFcnFromFile(path, spname);
}

MelaPConstant::~MelaPConstant(){}

void MelaPConstant::GetFcnFromFile(const char* path, const char* spname){
  TFile* fin = TFile::Open(path, "read");
  gROOT->cd();
  TString spname_core = spname;
  spname_core.Append("_Smooth");
  spname_core.Prepend("sp_tg_");
  TString spname_low = spname_core; spname_low.Prepend("lowFcn_"); fcnLow = (TF1*)fin->Get(spname_low);
  TString spname_high = spname_core; spname_high.Prepend("highFcn_"); fcnHigh = (TF1*)fin->Get(spname_high);
  TString spname_mid = spname_core; fcnMid = (TSpline3*)fin->Get(spname_mid);
  fin->Close();
}

double MelaPConstant::Eval(MelaIO* RcdME)const{
  double result=1;
  if (RcdME->melaCand->genStatus==-1) return result;

  double candMass = RcdME->melaCand->m();
  if (candMass<=0.) return 0.;
  else if (fcnLow!=0 && candMass<fcnLow->GetXmax()) result = fcnLow->Eval(candMass);
  else if (fcnHigh!=0 && candMass>fcnHigh->GetXmin()) result = fcnHigh->Eval(candMass);
  else if (fcnMid!=0) result = fcnMid->Eval(candMass);
  else return result;

  if (
    (processME==TVar::JHUGen || processME==TVar::MCFM)
    &&
    processProd==TVar::ZZGG
    &&
    processProc==TVar::HSMHiggs
    ){
    result = pow(10., result);

    double alphasVal, propagator, mh, gah;
    alphasVal = RcdME->getAlphaSatMZ();
    RcdME->getHiggsMassWidth(mh, gah, 0);
    propagator = 1./(pow(pow(candMass, 2)-pow(mh, 2), 2) + pow(mh*gah, 2));
    result *= pow(alphasVal, 2);
    result *= propagator;

    double aL1, aR1, aL2, aR2;
    RcdME->getVDaughterCouplings(aL1, aR1, 0);
    RcdME->getVDaughterCouplings(aL2, aR2, 1);
    if (fabs(aL1)>0. || fabs(aR1)>0.) result *= pow(aL1, 2)+pow(aR1, 2);
    if (fabs(aL2)>0. || fabs(aR2)>0.) result *= pow(aL2, 2)+pow(aR2, 2);
  }
  else if (
    processME==TVar::MCFM
    &&
    processProd==TVar::ZZGG
    &&
    processProc==TVar::bkgZZ
    ){
    result = pow(10., result);

    double alphasVal;
    alphasVal = RcdME->getAlphaSatMZ();
    result *= pow(alphasVal, 2);

    double aL1, aR1, aL2, aR2;
    RcdME->getVDaughterCouplings(aL1, aR1, 0);
    RcdME->getVDaughterCouplings(aL2, aR2, 1);
    if (fabs(aL1)>0. || fabs(aR1)>0.) result *= pow(aL1, 2)+pow(aR1, 2);
    if (fabs(aL2)>0. || fabs(aR2)>0.) result *= pow(aL2, 2)+pow(aR2, 2);
  }
  else if (
    processME==TVar::MCFM
    &&
    processProd==TVar::ZZQQB
    &&
    processProc==TVar::bkgZZ
    ){
    result = pow(10., result);

    double aL1, aR1, aL2, aR2;
    RcdME->getVDaughterCouplings(aL1, aR1, 0);
    RcdME->getVDaughterCouplings(aL2, aR2, 1);
    if (fabs(aL1)>0. || fabs(aR1)>0.) result *= pow(aL1, 2)+pow(aR1, 2);
    if (fabs(aL2)>0. || fabs(aR2)>0.) result *= pow(aL2, 2)+pow(aR2, 2);
  }
  else if (
    processME==TVar::JHUGen
    &&
    processProd==TVar::JQCD
    &&
    processProc==TVar::HSMHiggs
    ){
    result = pow(10., result);

    double alphasVal;
    alphasVal = RcdME->getAlphaSatMZ();
    result *= pow(alphasVal, 3);
  }
  else if (
    processME==TVar::JHUGen
    &&
    processProd==TVar::JJQCD
    &&
    processProc==TVar::HSMHiggs
    ){
    result = pow(10., result);

    double alphasVal;
    alphasVal = RcdME->getAlphaSatMZ();
    result *= pow(alphasVal, 4);
  }

  return result;
}



