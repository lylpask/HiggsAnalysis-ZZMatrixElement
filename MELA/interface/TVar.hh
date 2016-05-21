#ifndef EvtProb_VAR
#define EvtProb_VAR

#include <cstring>
#include <string>
#include <vector>
#include <pair>
#include "TLorentzVector.h"
#include "TCouplings.hh"
#include "MelaIO.h"
#include "TH2F.h"
#include "TH1F.h"

//#define EBEAM 4000.00
#define fbGeV2 0.389379E12
#define smallnumber 1e-15
#define sixteen_2Pi_to_8 3.88650230418250561e+07
#define eight_2Pi_to_5 7.83410393050320417e+04
#define four_2Pi_to_2 39.478417604357432


class TVar{
public:

  enum{
    kNoAssociated=1,
    kUseAssociated_Leptons=2, // l or nu
    kUseAssociated_Photons=3,
    kUseAssociated_Jets=5
  };
  enum VerbosityLevel {
    ERROR = 0,
    INFO = 1,
    DEBUG = 2
  };
  enum MatrixElement{
    MCFM = 0,
    JHUGen = 1,
    ANALYTICAL = 2
  };
  enum Production{
    ZZGG = 0,
    ZZQQB = 1,
    ZZQQB_STU = 2,
    ZZQQB_S = 3,
    ZZQQB_TU = 4,
    ZZINDEPENDENT= 5,
    JJGG = 6, // SBF
    JJVBF = 7, // WBF
    JH = 8, // H + 1 jet
    ZH = 9, // ZH
    WH = 10, // W(+/-)H
    ttH = 11, // ttH
    bbH = 12 // bbH
    //
  };
  enum LeptonInterference{
    DefaultLeptonInterf = 0,
    InterfOn = 1,
    InterfOff=2
  };

  enum Process{

    HSMHiggs,    //0+, call this for MCFM |H|**2.
    H0minus,    //0-
    H0hplus,    //0h+
    H0_g1prime2,   //g1=0, g1prime2=-12046.01
    H0_Zgs,
    H0_gsgs,
    H0_Zgs_PS,
    H0_gsgs_PS,
    H0_Zgsg1prime2,

    D_g1g4, // D_CP
    D_g1g4_pi_2, // D_CP_T
    D_g1g2,   // D_int
    D_g1g2_pi_2, // D_int_T
    D_g1g1prime2, // D_int_lambda1
    D_zzzg,
    D_zzgg,
    D_zzzg_PS,
    D_zzgg_PS,
    D_zzzg_g1prime2,
    D_zzzg_g1prime2_pi_2,

    H1minus,    //1-
    H1plus,    //1+

    H2_g1, //2m+, Zg, gg
    H2_g2, // 2h2+
    H2_g3, // 2h3+
    H2_g4, //2h+
    H2_g5, //2b+
    H2_g1g5, //2m+
    H2_g6, // 2h6+
    H2_g7, // 2h7+
    H2_g8, //2h-
    H2_g9, // 2h9-
    H2_g10, // 2h10-

    bkgZZ,    //qq->ZZ
    bkgZZ_SMHiggs,    //ggZZ+SMHiggs, ggZZ always calculated by MCFM, ME stands for SMHiggs ME, JHUGen: MCFM ggZZ + JHUGen SMHiggs, MCFM: MCFM (ggZZ+ SMHiggs)
    /**** For width ***/
    D_gg10,

    /***** Self Defined******/
    SelfDefine_spin0,
    SelfDefine_spin1,
    SelfDefine_spin2,

    Null
  };
  enum SuperMelaSyst{
    SMSyst_None      = 0, // nominal value
    SMSyst_ScaleUp   = 1, //Scale Uncertaintie
    SMSyst_ScaleDown = 2,
    SMSyst_ResUp     = 3, // Resolution Uncertainty
    SMSyst_ResDown   = 4
  };
  enum EventScaleScheme{
    DefaultScaleScheme,
    Fixed_mH,
    Fixed_mW,
    Fixed_mZ,
    Fixed_mWPlusmH,
    Fixed_mZPlusmH,
    Fixed_TwomtPlusmH,
    Fixed_mtPlusmH,
    Dynamic_qH,
    Dynamic_qJJH,
    Dynamic_qJJ_qH,
    Dynamic_qJ_qJ_qH,
    Dynamic_HT
  };

  //---------------------------------
  // Function
  //---------------------------------
  static TString ProcessName(int temp){
    if (temp==TVar::HSMHiggs) return TString("HSMHiggs");
    else if (temp==TVar::H0minus) return TString("H0minus");
    else if (temp==TVar::H0hplus) return TString("H0hplus");

    else if (temp==TVar::H1minus) return TString("H1minus");
    else if (temp==TVar::H1plus) return TString("H1plus");

    else if (temp==TVar::H2_g8) return TString("H2_g8");
    else if (temp==TVar::H2_g4) return TString("H2_g4");
    else if (temp==TVar::H2_g5) return TString("H2_g5");
    else if (temp==TVar::H2_g1g5) return TString("H2_g1g5");
    else if (temp == TVar::H2_g2) return TString("H2_g2");
    else if (temp == TVar::H2_g3) return TString("H2_g3");
    else if (temp == TVar::H2_g6) return TString("H2_g6");
    else if (temp == TVar::H2_g7) return TString("H2_g7");
    else if (temp == TVar::H2_g9) return TString("H2_g9");
    else if (temp == TVar::H2_g10) return TString("H2_g10");


    else if (temp==TVar::bkgZZ) return TString("bkgZZ");
    else if (temp==TVar::bkgZZ_SMHiggs) return TString("bkgZZ_SMHiggs");

    else if (temp==TVar::H0_g1prime2) return TString("H0_g1prime2");

    else if (temp==TVar::D_g1g4) return TString("D_g1g4");
    else if (temp==TVar::D_g1g4_pi_2) return TString("D_g1g4_pi_2");
    else if (temp==TVar::D_g1g2) return TString("D_g1g2");
    else if (temp==TVar::D_g1g2_pi_2) return TString("D_g1g2_pi_2");
    else if (temp==TVar::D_g1g1prime2) return TString("D_g1g1prime2");

    else if (temp==TVar::SelfDefine_spin0) return TString("SelfDefine_spin0");
    else if (temp==TVar::SelfDefine_spin1) return TString("SelfDefine_spin1");
    else if (temp==TVar::SelfDefine_spin2) return TString("SelfDefine_spin2");

    else if (temp==TVar::D_gg10) return TString("D_gg10");

    else if (temp==TVar::H0_Zgs) return TString("H0_Zgs");
    else if (temp==TVar::H0_gsgs) return TString("H0_gsgs");
    else if (temp==TVar::D_zzzg) return TString("D_zzzg");
    else if (temp==TVar::D_zzgg) return TString("D_zzgg");

    else if (temp==TVar::H0_Zgs_PS) return TString("H0_Zgs_PS");
    else if (temp==TVar::H0_gsgs_PS) return TString("H0_gsgs_PS");
    else if (temp==TVar::D_zzzg_PS) return TString("D_zzzg_PS");
    else if (temp==TVar::D_zzgg_PS) return TString("D_zzgg_PS");

    else if (temp==TVar::H0_Zgsg1prime2) return TString("H0_Zgsg1prime2");
    else if (temp==TVar::D_zzzg_g1prime2) return TString("D_zzzg_g1prime2");

    else return TString("UnKnown");
  };

  inline virtual ~TVar(){};
  ClassDef(TVar, 0)
};


struct branch_particle {
  int   PdgCode   ;
  int   Charge    ;
  double Px       ;
  double Py       ;
  double Pz       ;
  double E        ;
  double Eta      ;
  double Phi      ;
};
static const TString branch_format_particle =
  "PdgCode/I:"
  "Charge/I:"
  "Px/D:"
  "Py/D:"
  "Pz/D:"
  "E/D:"
  "Eta/D:"
  "Phi/D";

struct simple_event_record{
  std::vector<int> intermediateVid; // Origin of daughters, not associated particles
  std::pair<std::vector<int>, std::vector<TLorentzVector>> pDaughters;
  std::pair<std::vector<int>, std::vector<TLorentzVector>> pAssociated;
  std::pair<std::vector<int>, std::vector<TLorentzVector>> pMothers;
};

/*
struct hzz4l_event_type{
  int PdgCode[4];
  TLorentzVector p[4];
  int extraParticle_PdgCode[2];
  TLorentzVector extraParticle_p[2];
};
struct vh_event_type{ // ME is 2 -> 3
  int PdgCode[3];
  TLorentzVector p[3]; // H, V-daughter 1, V-daughter 2
  TLorentzVector pHdecay[4]; // H-daughter 1, H-daughter 2; optional: H-daughter 3, H-daughter 4
  int PdgCode_Hdecay[4];
};
struct tth_event_type{ // ME is 2 -> 8
  int PdgCode_tdecay[2][3];
  TLorentzVector p[7]; // H, tbar-daughters 1-3, t-daughters 1-3
};
struct mcfm_event_type{
  int PdgCode[8];
  TLorentzVector p[8];
  double pswt;
};
*/

struct event_scales_type{
  TVar::EventScaleScheme renomalizationScheme;
  TVar::EventScaleScheme factorizationScheme;
  double ren_scale_factor;
  double fac_scale_factor;
};


struct event_type{
  TLorentzVector p1,p2,ep,em,nu,nb;
  double PSWeight;
};
struct anomcoup{
  double delg1_z, delg1_g, lambda_g, lambda_z, delk_g, delk_z_,tevscale;
};
struct EffHist{
  TH2F* els_eff_mc;
  TH2F* mus_eff_mc;
};


#endif
