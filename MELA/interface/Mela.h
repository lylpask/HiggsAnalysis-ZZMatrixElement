/*
************* HEADER: CMS MELA interface to MCFM/JHUGen-MELA *************
Please see the ../src/Mela.cc file for the instructions.
*/

#ifndef MELA_Mela_h
#define MELA_Mela_h

#include "TLorentzVector.h"
#include <vector>
#include <TRandom3.h>


class TFile; 
class TH1F; 
class TH2F;
class TH3F;
class RooRealVar;
class RooAbsPdf;
class RooArgSet;
class ScalarPdfFactory_ggH;
class VectorPdfFactory;
class TensorPdfFactory;
class RooqqZZ_JHU_ZgammaZZ_fast;
class newZZMatrixElement;
class TGraph;
class SuperMELA;

#include <ZZMatrixElement/MELA/interface/TVar.hh>
#include <ZZMatrixElement/MELA/interface/TEvtProb.hh>
#include <ZZMatrixElement/MELA/interface/ScalarPdfFactory_ggH.h>
#include <ZZMatrixElement/MELA/interface/VectorPdfFactory.h>
#include <ZZMatrixElement/MELA/interface/TensorPdfFactory.h>
#include <ZZMatrixElement/MELA/interface/RooqqZZ_JHU_ZgammaZZ_fast.h>

class Mela{

public:

  // Mela(){};
  Mela(int LHCsqrts=13, float mh=125); // higgs mass for supermela
  ~Mela();

  void setProcess(TVar::Process myModel, TVar::MatrixElement myME, TVar::Production myProduction);
  void setMelaHiggsMass(double myHiggsMass, int index=0);
  void setMelaHiggsWidth(double myHiggsWidth=-1, int index=0);
  void setMelaLeptonInterference(TVar::LeptonInterference myLepInterf=TVar::DefaultLeptonInterf);
  void setRemoveLeptonMasses(bool MasslessLeptonSwitch=true);
  void resetMCFM_EWKParameters(double ext_Gf, double ext_aemmz, double ext_mW, double ext_mZ, double ext_xW, int ext_ewscheme=3);
  void setRenFacScaleMode(TVar::EventScaleScheme renormalizationSch, TVar::EventScaleScheme factorizationSch, double ren_sf, double fac_sf);


  MelaIO* getIORecord(); // Full parton-by-parton ME record
  MELACandidate* getCurrentCandidate();
  int getCurrentCandidateIndex();
  std::vector<MELATopCandidate*>* getTopCandidateCollection();


  float getConstant(bool useOldggZZConstants=false);
  float getConstant_m4l(bool useOldggZZConstants=false);
  void get_PAux(float& prob); // SuperProb

  void computeP_selfDspin0(
    double selfDHvvcoupl_input[nSupportedHiggses][SIZE_HVV][2],
    float& prob,
    bool useConstant=false
    );
  void computeP_selfDspin0(
    float& prob,
    bool useConstant=false
    );

  void computeP_selfDspin1(
    double selfDZvvcoupl_input[SIZE_ZVV][2],
    float& prob,
    bool useConstant=false
    );
  void computeP_selfDspin1(
    float& prob,
    bool useConstant=false
    );

  void computeP_selfDspin2(
    double selfDGggcoupl_input[SIZE_GGG][2],
    double selfDGvvcoupl_input[SIZE_GVV][2],
    float& prob,
    bool useConstant=false
    );
  void computeP_selfDspin2(
    float& prob,
    bool useConstant=false
    );

  void computeP(
    double selfDHvvcoupl_freenorm_input[SIZE_HVV_FREENORM],
    float& prob,
    bool useConstant=true
    );
  void computeP(
    float& prob,
    bool useConstant=true
    );

  void computeD_CP(
    TVar::MatrixElement myME,
    TVar::Process myType,
    float& prob
    );


  //****VVH Spin-0****//
  void computeProdDecP(
    float& prob,
    bool useConstant=true
    );
  void computeProdDecP(
    float& prob,
    bool useConstant=true
    );

  //****HJ/HJJ/VBF Spin-0****//
  void computeProdP(
    double selfDHggcoupl_input[SIZE_HGG][2],
    double selfDHvvcoupl_input[nSupportedHiggses][SIZE_HVV][2],
    double selfDHwwcoupl_input[nSupportedHiggses][SIZE_HVV][2],
    float& prob,
    bool useConstant=true
    );
  void computeProdP(
    float& prob,
    bool useConstant=true
    );

  //****VH Spin-0****//
  void computeProdP(
    double selfDHvvcoupl_input[nSupportedHiggses][SIZE_HVV][2],
    float& prob,
    bool useConstant=true
    );
  void computeProdP(
    float& prob,
    bool includeHiggsDecay=false,
    bool useConstant=true
    );

  //***ttH Spin-0****//
  void computeProdP(
    float& prob,
    int topDecay=0,
    int topProcess=2,
    bool useConstant=true
    );


  void computePM4l(
    float mZZ,
    TVar::LeptonFlavor flavor,
    TVar::SuperMelaSyst syst,
    float& prob
    ); //SuperMela

  void computePM4l(
    TVar::SuperMelaSyst syst,
    float& prob
    ); //SuperMela


  // Access ZZMEs Calculate4Momentum
  std::vector<TLorentzVector> calculate4Momentum(double Mx, double M1, double M2, double theta, double theta1, double theta2, double Phi1, double Phi);

  // Calculation weight to correct for lepton interference
  void computeWeight(
    // return variables:
    float& w
    );

  void computeWeight(
    double selfDHvvcoupl_freenorm[SIZE_HVV_FREENORM],
    // return variables:
    float& w
    );

  RooAbsPdf* pdf;
  ScalarPdfFactory_ggH* ggSpin0Model;
  VectorPdfFactory* spin1Model;
  TensorPdfFactory* spin2Model;
  RooqqZZ_JHU_ZgammaZZ_fast* qqZZmodel;
  SuperMELA* super;
  TRandom3 *myR; // random number for resolution systematics


  RooRealVar* mzz_rrv;
  RooRealVar* z1mass_rrv;
  RooRealVar* z2mass_rrv;
  RooRealVar* costhetastar_rrv;
  RooRealVar* costheta1_rrv;
  RooRealVar* costheta2_rrv;
  RooRealVar* phi_rrv;
  RooRealVar* phi1_rrv;
  RooRealVar* Y_rrv;
  RooRealVar* upFrac_rrv;

  TGraph* vaScale_4e;
  TGraph* vaScale_4mu;
  TGraph* vaScale_2e2mu;
  TGraph* DggZZ_scalefactor;
  void setCTotalBkgGraphs(TFile* fcontainer, TGraph* tgC[]);
  void constructDggr(
    float mzz, float bkg_VAMCFM_noscale,
    float ggzz_VAMCFM_noscale, float ggHZZ_prob_pure_noscale,
    float ggHZZ_prob_int_noscale, float& myDggr
    );

  TGraph* tgtotalbkg[3];
  void computeD_gg(
    float mZZ, float mZ1, float mZ2, // input kinematics
    float costhetastar,
    float costheta1,
    float costheta2,
    float phi,
    float phi1,
    int flavor,
    TVar::MatrixElement myME,
    TVar::Process myType,
    float& prob
    );

  // Self-define arrays are now members of MELA.
  // There are a lot of them!
  //****Spin-0****//
  double selfDHvvcoupl_freenorm[SIZE_HVV_FREENORM];
  double selfDHqqcoupl[SIZE_HQQ][2];
  double selfDHggcoupl[SIZE_HGG][2];
  // The first dimension (of size [nSupportedHiggses=2]) supports a second resonance present in MCFM
  double selfDHzzcoupl[nSupportedHiggses][SIZE_HVV][2];
  double selfDHwwcoupl[nSupportedHiggses][SIZE_HVV][2];
  double selfDHzzLambda_qsq[nSupportedHiggses][4][3];
  double selfDHwwLambda_qsq[nSupportedHiggses][4][3];
  int selfDHzzCLambda_qsq[nSupportedHiggses][3];
  int selfDHwwCLambda_qsq[nSupportedHiggses][3];
  bool differentiate_HWW_HZZ;
  //****Spin-1****//
  double selfDZqqcoupl[SIZE_ZQQ][2];
  double selfDZvvcoupl[SIZE_ZVV][2];
  //****Spin-2****//
  double selfDGqqcoupl[SIZE_GQQ][2];
  double selfDGggcoupl[SIZE_GGG][2];
  double selfDGvvcoupl[SIZE_GVV][2];
  // That is a lot of them!

protected:

  // 
  // data memmbers 
  // 
  bool usePowhegTemplate_;
  int LHCsqrts;
  TVar::Process myModel_;
  TVar::MatrixElement myME_;
  TVar::Production myProduction_;
  newZZMatrixElement* ZZME;

  float auxiliaryProb;

  MELACandidate* melaCand; // Pointer to persistent TEvtProb object

  //
  // Functions
  //
  void configureAnalyticalPDFs();
  void reset_SelfDCouplings();
  void reset_PAux(); // SuperProb reset
  void reset_CandRef();
};

#endif

