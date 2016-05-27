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
#include <vector>
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
#include <ZZMatrixElement/MELA/interface/HiggsCSandWidth_MELA.h>


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
  // Functions
  //----------------------
  void SetProcess(TVar::Process tmp) { process = tmp; }
  void SetMatrixElement(TVar::MatrixElement tmp){ matrixElement = tmp; }
  void SetProduction(TVar::Production tmp){ production = tmp; }
  void SetVerbosity(TVar::VerbosityLevel tmp){ verbosity = tmp; }
  void SetLeptonInterf(TVar::LeptonInterference tmp){ leptonInterf = tmp; }

  void SetCurrentCandidate(unsigned int icand);
  void SetCurrentCandidate(MELACandidate* cand);

  void AllowSeparateWWCouplings(bool doAllow=false){ SetJHUGenDistinguishWWCouplings(doAllow); selfDSpinZeroCoupl.allow_WWZZSeparation(doAllow); }
  void ResetMCFM_EWKParameters(double ext_Gf, double ext_aemmz, double ext_mW, double ext_mZ, double ext_xW, int ext_ewscheme=3);
  void Set_LHAgrid(const char* path, int pdfmember=0);

  void SetInputEvent(
    std::vector<std::pair<int, TLorentzVector>>* pDaughters,
    std::vector<std::pair<int, TLorentzVector>>* pAssociated=0,
    std::vector<std::pair<int, TLorentzVector>>* pMothers=0,
    bool isGen=false
    );
  void AppendTopCandidate(std::vector<std::pair<int, TLorentzVector>>* TopDaughters);

  double XsecCalc_XVV(
    TVar::Process process_, TVar::Production production_,
    TVar::VerbosityLevel verbosity_
    );

  double XsecCalc_VVXVV(
    TVar::Process process_, TVar::Production production_,
    TVar::VerbosityLevel verbosity_
    );

  double XsecCalcXJJ(
    TVar::Process process_, TVar::Production production_,
    TVar::VerbosityLevel verbosity_
    );

  double XsecCalcXJ(
    TVar::Process process_, TVar::Production production_,
    TVar::VerbosityLevel verbosity_
    );

  double XsecCalc_VX(
    TVar::Process process_, TVar::Production production_,
    TVar::VerbosityLevel verbosity_
    );

  double XsecCalc_TTX(
    TVar::Process process_, TVar::Production production_,
    TVar::VerbosityLevel verbosity_,
    int topProcess, int topDecay
    );

  // this appears to be some kind of
  // way of setting MCFM parameters through
  // an interface defined in TMCFM.hh
  void SetHiggsMass(double mass, double wHiggs=-1., int whichResonance=-1);

  void SetRenFacScaleMode(TVar::EventScaleScheme renormalizationSch, TVar::EventScaleScheme factorizationSch, double ren_sf, double fac_sf);
  void ResetRenFacScaleMode(){ SetRenFacScaleMode(TVar::DefaultScaleScheme, TVar::DefaultScaleScheme, 0.5, 0.5); };

  void ResetCouplings(){ selfDSpinZeroCoupl.reset(); selfDSpinOneCoupl.reset(); selfDSpinTwoCoupl.reset(); AllowSeparateWWCouplings(false); };
  void ResetIORecord(){ RcdME.reset(); };
  void ResetInputEvent();


  // Get-functions
  SpinZeroCouplings* GetSelfDSpinZeroCouplings(){ return selfDSpinZeroCoupl.getRef(); };
  SpinOneCouplings* GetSelfDSpinOneCouplings(){ return selfDSpinOneCoupl.getRef(); };
  SpinTwoCouplings* GetSelfDSpinTwoCouplings(){ return selfDSpinTwoCoupl.getRef(); };
  MelaIO* GetIORecord(){ return RcdME.getRef(); };
  MELACandidate* GetCurrentCandidate();
  int GetCurrentCandidateIndex();

private:
  //--------------------
  // Variables
  //--------------------
  TVar::Process process;
  TVar::MatrixElement matrixElement;
  TVar::Production production;
  TVar::VerbosityLevel verbosity;
  TVar::LeptonInterference leptonInterf;
  double _hmass;
  double _hwidth;
  double _h2mass;
  double _h2width;
  double EBEAM;
  HiggsCSandWidth_MELA* myCSW_;
  event_scales_type event_scales;

  SpinZeroCouplings selfDSpinZeroCoupl;
  SpinOneCouplings selfDSpinOneCoupl;
  SpinTwoCouplings selfDSpinTwoCoupl;
  MelaIO RcdME;

  MELACandidate* melaCand; // Only a pointer to the top-level (input) candList object
  std::vector<MELAParticle*> particleList; // Container of intermediate objects, for bookkeeping to delete later
  std::vector<MELACandidate*> candList; // Container of candidate objects, for bookkeeping to delete later
  std::vector<MELATopCandidate*> topCandList; // Container of candidate objects, for bookkeeping to delete later

  // Convert std::vectors to MELAPArticle* and MELACandidate* objects, stored in particleList and candList, respectively
  MELACandidate* ConvertVectorFormat(
    std::vector<std::pair<int, TLorentzVector>>* pDaughters,
    std::vector<std::pair<int, TLorentzVector>>* pAssociated=0,
    std::vector<std::pair<int, TLorentzVector>>* pMothers=0,
    bool isGen
    );
  MELATopCandidate* ConvertTopCandidate(std::vector<std::pair<int, TLorentzVector>>* TopDaughters);
  // Check if at least one input candidate is present
  bool CheckInputPresent();
  void SetRcdCandPtr();


  ClassDef(TEvtProb,0);
};

#endif

