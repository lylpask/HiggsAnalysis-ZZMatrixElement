#ifndef _TEVTPROB_HH_
#define _TEVTPROB_HH_
//-----------------------------------------------------------------------------
// Description: Class TEvtProb: EvtProb base class
// ------------
//
//      Event Probability Density Calculation
//
// Feb 21 2011
// Sergo Jindariani
// Yanyan Gao
//-----------------------------------------------------------------------------
#include <sstream>
#include <cstdio>
#include <string>

#include <iostream>
#include <iomanip>
#include <ostream>
#include <fstream>

#include "TObject.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TString.h"
#include "assert.h"
#include "TROOT.h"
// ME related
#include "TMCFM.hh"
#include "TCouplings.hh"
#include "TVar.hh"
#include "TUtil.hh"
#include "ZZMatrixElement/MELA/interface/HiggsCSandWidth_MELA.h"


//----------------------------------------
// Class TEvtProb
//----------------------------------------
class TEvtProb : public TObject {
  
public:
  //---------------------------------------------------------------------------
  // Constructors and Destructor
  //---------------------------------------------------------------------------
  TEvtProb() {};
  TEvtProb(const char* path, double ebeam, const char* pathtoPDFSet, int PDFMember=0);
  ~TEvtProb();
  
  //----------------------
  // Function
  //----------------------
  void SetProcess(TVar::Process tmp) { _process = tmp; }
  void SetMatrixElement(TVar::MatrixElement tmp){ _matrixElement = tmp; }
  void SetProduction(TVar::Production tmp){ _production = tmp; }
  void SetLeptonInterf(TVar::LeptonInterference tmp){ _leptonInterf = tmp; }
  void AllowSeparateWWCouplings(bool doAllow=false){ SetJHUGenDistinguishWWCouplings(doAllow); selfDSpinZeroCoupl.allow_WWZZSeparation(doAllow); }
  void ResetMCFM_EWKParameters(double ext_Gf, double ext_aemmz, double ext_mW, double ext_mZ, double ext_xW, int ext_ewscheme=3);
  void Set_LHAgrid(const char* path, int pdfmember=0);

  double XsecCalc(
    TVar::Process proc, TVar::Production production,
    const hzz4l_event_type &hzz4l_event,
    TVar::VerbosityLevel verbosity
    );

  double XsecCalc_VVXVV(
    TVar::Process proc, TVar::Production production,
    const hzz4l_event_type &hzz4l_event,
    TVar::VerbosityLevel verbosity
    );

  double XsecCalcXJJ(
    TVar::Process proc, TVar::Production production, 
    TLorentzVector p4[3],
    TVar::VerbosityLevel verbosity
    );

  double XsecCalcXJ(
    TVar::Process proc, TVar::Production production,
    TLorentzVector p4[2],
    TVar::VerbosityLevel verbosity
    );

  double XsecCalc_VX(
    TVar::Process proc, TVar::Production production,
    vh_event_type &vh_event,
    TVar::VerbosityLevel verbosity
    );

  double XsecCalc_TTX(
    TVar::Process proc, TVar::Production production,
    tth_event_type &tth_event,
    int topDecay, int topProcess,
    TVar::VerbosityLevel verbosity
    );

  // this appears to be some kind of 
  // way of setting MCFM parameters through
  // an interface defined in TMCFM.hh
  void SetHiggsMass(double mass, float wHiggs=-1);

  void SetRenFacScaleMode(TVar::EventScaleScheme renormalizationSch, TVar::EventScaleScheme factorizationSch, double ren_sf, double fac_sf);
  void ResetRenFacScaleMode(){ SetRenFacScaleMode(TVar::DefaultScaleScheme, TVar::DefaultScaleScheme, 0.5, 0.5); };

  void ResetCouplings(){ selfDSpinZeroCoupl.reset(); selfDSpinOneCoupl.reset(); selfDSpinTwoCoupl.reset(); AllowSeparateWWCouplings(false); };
  void ResetIORecord(){ RcdME.reset(); };


  // Get-functions
  SpinZeroCouplings* GetSelfDSpinZeroCouplings(){ return selfDSpinZeroCoupl.getRef(); };
  SpinOneCouplings* GetSelfDSpinOneCouplings(){ return selfDSpinOneCoupl.getRef(); };
  SpinTwoCouplings* GetSelfDSpinTwoCouplings(){ return selfDSpinTwoCoupl.getRef(); };
  MelaIO* GetIORecord(){ return RcdME.getRef(); };

private:
  //--------------------
  // Variables
  //--------------------
  TVar::Process _process;
  TVar::MatrixElement _matrixElement;
  TVar::Production _production;
  TVar::LeptonInterference _leptonInterf;
  double _hmass;
  double _hwidth;
  double EBEAM;
  event_scales_type event_scales;

  SpinZeroCouplings selfDSpinZeroCoupl;
  SpinOneCouplings selfDSpinOneCoupl;
  SpinTwoCouplings selfDSpinTwoCoupl;
  MelaIO RcdME;

  HiggsCSandWidth_MELA *myCSW_;

  ClassDef(TEvtProb,0);
};

#endif

