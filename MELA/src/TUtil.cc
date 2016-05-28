#include "ZZMatrixElement/MELA/interface/TUtil.hh"
#include <iostream>
#include <cstdio>
#include <cmath>
#include <TMath.h>


using namespace std;


void TUtil::SetEwkCouplingParameters(){
  ewscheme_.ewscheme = 3; // Switch ewscheme to full control, default is 1

/*
// Settings used in Run I MC, un-synchronized to JHUGen and Run II generation
  ewinput_.Gf_inp = 1.16639E-05;
  ewinput_.wmass_inp = 80.385;
  ewinput_.zmass_inp = 91.1876;
//  ewinput_.aemmz_inp = 7.81751e-3; // Not used in the compiled new default MCFM ewcheme=1
//  ewinput_.xw_inp = 0.23116864; // Not used in the compiled new default MCFM ewcheme=1
  ewinput_.aemmz_inp = 7.562468901984759e-3;
  ewinput_.xw_inp = 0.2228972225239183;
*/

// SETTINGS TO MATCH JHUGen ME/generator:
  ewinput_.Gf_inp=1.16639E-05;
  ewinput_.aemmz_inp=1./128.;
  ewinput_.wmass_inp=80.399;
  ewinput_.zmass_inp=91.1876;
  ewinput_.xw_inp=0.23119;

/*
// INPUT SETTINGS in default MCFM generator:
  ewscheme_.ewscheme = 1;
  ewinput_.Gf_inp=1.16639E-05;
  ewinput_.wmass_inp=80.398;
  ewinput_.zmass_inp=91.1876;
  ewinput_.aemmz_inp=7.7585538055706e-3;
  ewinput_.xw_inp=0.2312;
*/
/*
// ACTUAL SETTINGS in default MCFM generator, gives the same values as above but is more explicit:
  ewscheme_.ewscheme = 3;
  ewinput_.Gf_inp=1.16639E-05;
  ewinput_.wmass_inp=80.398;
  ewinput_.zmass_inp=91.1876;
  ewinput_.aemmz_inp=7.55638390728736e-3;
  ewinput_.xw_inp=0.22264585341299625;
*/
}
double TUtil::InterpretScaleScheme(TVar::Production production, TVar::MatrixElement matrixElement, TVar::EventScaleScheme scheme, TLorentzVector p[mxpart]){
  double Q=0;
  TLorentzVector nullFourVector(0, 0, 0, 0);
  if (scheme == TVar::Fixed_mH) Q = masses_mcfm_.hmass;
  else if (scheme == TVar::Fixed_mW) Q = masses_mcfm_.wmass;
  else if (scheme == TVar::Fixed_mZ) Q = masses_mcfm_.zmass;
  else if (scheme == TVar::Fixed_mWPlusmH) Q = (masses_mcfm_.wmass+masses_mcfm_.hmass);
  else if (scheme == TVar::Fixed_mZPlusmH) Q = (masses_mcfm_.zmass+masses_mcfm_.hmass);
  else if (scheme == TVar::Fixed_TwomtPlusmH) Q = (2.*masses_mcfm_.mt+masses_mcfm_.hmass);
  else if (scheme == TVar::Fixed_mtPlusmH) Q = masses_mcfm_.mt+masses_mcfm_.hmass;
  else if (scheme == TVar::Dynamic_qH){
    TLorentzVector pTotal = p[2]+p[3]+p[4]+p[5];
    Q = fabs(pTotal.M());
  }
  else if (scheme == TVar::Dynamic_qJJH){
    TLorentzVector pTotal = p[2]+p[3]+p[4]+p[5]+p[6]+p[7];
    Q = fabs(pTotal.M());
  }
  else if (scheme == TVar::Dynamic_qJJ_qH){
    TLorentzVector qH = p[2]+p[3]+p[4]+p[5];
    TLorentzVector qJJ = p[6]+p[7];
    Q = fabs(qH.M()+qJJ.M());
  }
  else if (scheme == TVar::Dynamic_qJ_qJ_qH){
    TLorentzVector qH = p[2]+p[3]+p[4]+p[5];
    Q = fabs(qH.M()+p[6].M()+p[7].M());
  }
  else if (scheme == TVar::Dynamic_HT){
    for (int c=2; c<mxpart; c++) Q += p[c].Pt(); // Scalar sum of all pTs
  }
  else if (scheme == TVar::DefaultScaleScheme){
    if (matrixElement==TVar::JHUGen){
      if (
        production==TVar::JJGG
        || production==TVar::JJVBF
        || production==TVar::JH
        || production==TVar::ZZGG
        || production==TVar::ZZQQB
        || production==TVar::ZZQQB_STU
        || production==TVar::ZZQQB_S
        || production==TVar::ZZQQB_TU
        || production==TVar::ZZINDEPENDENT
        || production==TVar::WH
        || production==TVar::ZH
        ){
        TLorentzVector pTotal = p[2]+p[3]+p[4]+p[5];
        Q = fabs(pTotal.M());
      }
      else if (production==TVar::ttH || production==TVar::bbH) Q = (2.*masses_mcfm_.mt+masses_mcfm_.hmass);
    }
    else if (matrixElement==TVar::MCFM){
      if (
        production==TVar::JJGG
        || production==TVar::JJVBF
        || production==TVar::JH
        || production==TVar::ZZGG
        || production==TVar::ZZQQB
        || production==TVar::ZZQQB_STU
        || production==TVar::ZZQQB_S
        || production==TVar::ZZQQB_TU
        || production==TVar::ZZINDEPENDENT
        || production==TVar::ZH
        || production==TVar::WH
        ){
        TLorentzVector pTotal = p[2]+p[3]+p[4]+p[5];
        Q = fabs(pTotal.M());
      }
      else if (
        production==TVar::ttH
        || production==TVar::bbH
        ){
        TLorentzVector pTotal = p[2]+p[3]+p[4]+p[5]+p[6]+p[7];
        Q = fabs(pTotal.M());
      }
    }
  }

  if (Q<=0){
    cout << "Scaling fails, defaulting to dynamic scheme m3456 " << endl;
    TLorentzVector pTotal = p[2]+p[3]+p[4]+p[5];
    Q = fabs(pTotal.M());
  }
  return Q;
}
void TUtil::SetAlphaS(double Q_ren, double Q_fac, double multiplier_ren, double multiplier_fac, int mynloop, int mynflav, string mypartons){
  bool hasReset=false;
  if (multiplier_ren<=0 || multiplier_fac<=0){
    cerr << "Invalid scale multipliers" << endl;
    return;
  }
  if (Q_ren<=1 || Q_fac<=1 || mynloop<=0 || mypartons.compare("Default")==0){
    if (Q_ren<0) cout << "Invalid QCD scale for alpha_s, setting to mH/2..." << endl;
    if (Q_fac<0) cout << "Invalid factorization scale, setting to mH/2..." << endl;
    Q_ren = (masses_mcfm_.hmass)*0.5;
    Q_fac = Q_ren;
    mynloop = 1;
    hasReset=true;
  }
  if (!hasReset){
    Q_ren *= multiplier_ren;
    Q_fac *= multiplier_fac;
  }

  /***** JHUGen Alpha_S *****/
  // To be implemented

  /***** MCFM Alpha_S *****/
  bool nflav_is_same = (nflav_.nflav == mynflav);
  scale_.scale = Q_ren;
  scale_.musq = Q_ren*Q_ren;
  facscale_.facscale = Q_fac;
  nlooprun_.nlooprun = mynloop;

  if (mypartons.compare("Default")!=0 && mypartons.compare("cteq6_l")!=0 && mypartons.compare("cteq6l1")!=0){
    cout << "Only default=cteq6l1 or cteq6_l are supported. Modify mela.cc symlinks, put the pdf table into data/Pdfdata and retry. Setting mypartons to Default..." << endl;
    mypartons = "Default";
  }
  // From pdfwrapper_linux.f:
  if (mypartons.compare("cteq6_l")==0) couple_.amz = 0.118;
  else if (mypartons.compare("cteq6l1")==0 || mypartons.compare("Default")==0) couple_.amz = 0.130;
  else couple_.amz = 0.118; // Add pdf as appropriate

  // For proper pdfwrapper_linux.f execution (alpha_s computation does not use pdf but examines the pdf name to initialize amz.)
  if (!nflav_is_same){
    nflav_.nflav = mynflav;

    if (mypartons.compare("Default")!=0) sprintf(pdlabel_.pdlabel, "%s", mypartons.c_str());
    else sprintf(pdlabel_.pdlabel, "%s", "cteq6l1"); // Default pdf is cteq6l1
    coupling2_();
  }
  else qcdcouple_.as = alphas_(&(scale_.scale), &(couple_.amz), &(nlooprun_.nlooprun));

  qcdcouple_.gsq = 4.0*TMath::Pi()*qcdcouple_.as;
  qcdcouple_.ason2pi = qcdcouple_.as/(2.0*TMath::Pi());
  qcdcouple_.ason4pi = qcdcouple_.as/(4.0*TMath::Pi());

  // TEST RUNNING SCALE PER EVENT:
  /*
  if(verbosity >= TVar::DEBUG){
  cout << "My pdf is: " << pdlabel_.pdlabel << endl;
  cout << "My Q_ren: " << Q_ren << " | My alpha_s: " << qcdcouple_.as << " at order " << nlooprun_.nlooprun << " with a(m_Z): " << couple_.amz << '\t'
  << "Nflav: " << nflav_.nflav << endl;
  */
}
bool TUtil::MCFM_chooser(TVar::Process process, TVar::Production production, TVar::LeptonInterference leptonInterf, MELACandidate* cand){
  MELAPArticle* V1 = cand->getSorteV(0);
  MELAPArticle* V2 = cand->getSorteV(1);
  if (V1==0 || V2==0) return false;

  unsigned int ndau = V1->getNDaughters() + V2->getNDaughters();
  unsigned int najets = cand->getNAssociatedJets();
  unsigned int naneutrinos = cand->getNAssociatedNeutrinos();
  unsigned int naleps = cand->getNAssociatedLeptons()-naneutrinos;
  unsigned int naphotons = cand->getNAssociatedPhotons();
  unsigned int naferms = naleps+naneutrinos+najets;
  unsigned int naparts = naferms+naphotons;
  bool definiteInterf=(
    ndau>3
    &&
    V1->getDaughter(0)->id==V1->getDaughter(0)->id
    &&
    V1->getDaughter(1)->id==V2->getDaughter(1)->id
    &&
    !isAnUnknownJet(V1->getDaughter(0)->id) && !isAnUnknownJet(V1->getDaughter(1)->id)
    );

  // VV->4f
  if ( // Check for support in qqZZ+0J
    ndau>=4
    &&
    process==TVar::bkgZZ
    &&
    (production == TVar::ZZQQB || production == TVar::ZZQQB_STU || production == TVar::ZZQQB_S || production == TVar::ZZQQB_TU || production == TVar::ZZINDEPENDENT)
    ){
    //81 '  f(p1)+f(p2) --> Z^0(-->mu^-(p3)+mu^+(p4)) + Z^0(-->e^-(p5)+e^+(p6))'
    //86 '  f(p1)+f(p2) --> Z^0(-->e^-(p5)+e^+(p6))+Z^0(-->mu^-(p3)+mu^+(p4)) (NO GAMMA*)'
    //90 '  f(p1)+f(p2) --> Z^0(-->e^-(p3)+e^+(p4)) + Z^0(-->e^-(p5)+e^+(p6))' 'L'

    // these settings are identical to use the chooser_() function
    npart_.npart=4;
    nqcdjets_.nqcdjets=0;

    nwz_.nwz=0;
    bveg1_mcfm_.ndim=10;
    breit_.n2=1;
    breit_.n3=1;

    breit_.mass2=masses_mcfm_.zmass;
    breit_.width2=masses_mcfm_.zwidth;
    breit_.mass3=masses_mcfm_.zmass;
    breit_.width3=masses_mcfm_.zwidth;

    zcouple_.q1=-1.0; // Pretty important for MCFM 6.8+, does not matter for earlier versions
    zcouple_.l1=zcouple_.le;
    zcouple_.r1=zcouple_.re;

    zcouple_.q2=-1.0; // Pretty important for MCFM 6.8+, does not matter for earlier versions
    zcouple_.l2=zcouple_.le;
    zcouple_.r2=zcouple_.re;

    vsymfact_.vsymfact=1.0;
    interference_.interference=false;
    if (definiteInterf && (leptonInterf==TVar::DefaultLeptonInterf || leptonInterf==TVar::InterfOn)){
      //90 '  f(p1)+f(p2) --> Z^0(-->e^-(p3)+e^+(p4)) + Z^0(-->e^-(p5)+e^+(p6))' 'L'
      vsymfact_.vsymfact=0.125; // MELA FACTOR (0.25 in MCFM 6.8)  --->   Result of just removing if(bw34_56) statements in FORTRAN code and not doing anything else
      //                vsymfact_.vsymfact=0.25; // MELA FACTOR (Same 0.25 in MCFM 6.7)
      interference_.interference=true;
    }

  }
  else if ( // Check for support in ggH+0J
    ndau>=4
    &&
    production == TVar::ZZGG
    &&
    (process==TVar::bkgZZ || process==TVar::HSMHiggs || process == TVar::bkgZZ_SMHiggs)
    ){
    // 114 '  f(p1)+f(p2) --> H(--> Z^0(mu^-(p3)+mu^+(p4)) + Z^0(e^-(p5)+e^+(p6))' 'N'
    /*
    nprocs:
    c--- 128 '  f(p1)+f(p2) --> H(--> Z^0(e^-(p3)+e^+(p4)) + Z^0(mu^-(p5)+mu^+(p6)) [top, bottom loops, exact]' 'L'
    c--- 129 '  f(p1)+f(p2) --> H(--> Z^0(e^-(p3)+e^+(p4)) + Z^0(mu^-(p5)+mu^+(p6)) [only H, gg->ZZ intf.]' 'L' -> NOT IMPLEMENTED
    c--- 130 '  f(p1)+f(p2) --> H(--> Z^0(e^-(p3)+e^+(p4)) + Z^0(mu^-(p5)+mu^+(p6)) [H squared and H, gg->ZZ intf.]' 'L'
    c--- 131 '  f(p1)+f(p2) --> Z^0(e^-(p3)+e^+(p4)) + Z^0(mu^-(p5)+mu^+(p6) [gg only, (H + gg->ZZ) squared]' 'L'
    c--- 132 '  f(p1)+f(p2) --> Z^0(e^-(p3)+e^+(p4)) + Z^0(mu^-(p5)+mu^+(p6) [(gg->ZZ) squared]' 'L'
    */

    npart_.npart=4;
    nqcdjets_.nqcdjets=0;

    bveg1_mcfm_.ndim=10;

    breit_.n2=1;
    breit_.n3=1;

    breit_.mass2 =masses_mcfm_.zmass;
    breit_.width2=masses_mcfm_.zwidth;
    breit_.mass3 =masses_mcfm_.zmass;
    breit_.width3=masses_mcfm_.zwidth;

    zcouple_.q1=-1.0;
    zcouple_.l1=zcouple_.le;
    zcouple_.r1=zcouple_.re;

    zcouple_.q2=-1.0;
    zcouple_.l2=zcouple_.le;
    zcouple_.r2=zcouple_.re;

    vsymfact_.vsymfact=1.0;
    interference_.interference=false;
    if (definiteInterf && (leptonInterf==TVar::DefaultLeptonInterf || leptonInterf==TVar::InterfOn)){
      vsymfact_.vsymfact=0.5;
      interference_.interference=true;
    }

  }
  // JJ + VV->4f
  else if ( // Check for support in qq'H+2J
    najets>=2 // Only jets since MCFM ME does not support leptons or neutrinos as associated particles
    &&
    ndau>=4
    &&
    production == TVar::JJVBF
    &&
    (process==TVar::bkgZZ || process==TVar::HSMHiggs || process == TVar::bkgZZ_SMHiggs)
    ){
    // 220 '  f(p1)+f(p2) --> Z(e-(p3),e^+(p4))Z(mu-(p5),mu+(p6)))+f(p7)+f(p8) [weak]' 'L'

    if (process==TVar::bkgZZ){
      wwbits_.Hbit[0]=0;
      wwbits_.Hbit[1]=0;
      wwbits_.Bbit[0]=1;
      wwbits_.Bbit[1]=0;
    }
    else if (process==TVar::bkgZZ_SMHiggs){
      wwbits_.Hbit[0]=1;
      wwbits_.Hbit[1]=0;
      wwbits_.Bbit[0]=1;
      wwbits_.Bbit[1]=0;
    }
    else{
      wwbits_.Hbit[0]=1;
      wwbits_.Hbit[1]=0;
      wwbits_.Bbit[0]=0;
      wwbits_.Bbit[1]=0;
    }

    npart_.npart=6;
    nwz_.nwz=2;

    zcouple_.q1=-1.0;
    zcouple_.l1=zcouple_.le;
    zcouple_.r1=zcouple_.re;

    zcouple_.q2=-1.0;
    zcouple_.l2=zcouple_.le;
    zcouple_.r2=zcouple_.re;

    bveg1_mcfm_.ndim=16;
    nqcdjets_.nqcdjets=2;

    breit_.n2=1;
    breit_.n3=1;

    breit_.mass2 =masses_mcfm_.zmass;
    breit_.width2=masses_mcfm_.zwidth;
    breit_.mass3 =masses_mcfm_.zmass;
    breit_.width3=masses_mcfm_.zwidth;

    vsymfact_.vsymfact=1.0;
    interference_.interference=false;
    if (definiteInterf && (leptonInterf==TVar::DefaultLeptonInterf || leptonInterf==TVar::InterfOn)){
      vsymfact_.vsymfact=1.0;
      interference_.interference=true;
    }

  }
  else{
    cerr <<"TUtil::MCFM_chooser: Can't identify Process: " << process << endl;
    cerr <<"TUtil::MCFM_chooser: ndau: " << ndau << '\t';
    cerr <<"TUtil::MCFM_chooser: naparts: " << naparts << '\t';
    cerr <<"TUtil::MCFM_chooser: nainvalid: " << nainvalid << '\t';
    cerr <<"TUtil::MCFM_chooser: najets: " << najets << '\t';
    cerr <<"TUtil::MCFM_chooser: naneutrinos: " << naneutrinos << '\t';
    cerr <<"TUtil::MCFM_chooser: naleps: " << naleps << '\t';
    cerr <<"TUtil::MCFM_chooser: naferms: " << naferms << '\t';
    cerr <<"TUtil::MCFM_chooser: naphotons: " << naphotons << '\t';
    cerr << endl;
    return false;
  }
  return true;
}

void TUtil::InitJHUGenMELA(const char* pathtoPDFSet, int PDFMember){
  char path_pdf_c[200];
  sprintf(path_pdf_c, "%s", pathtoPDFSet);
  int pathpdfLength = strlen(path_pdf_c);
  __modjhugen_MOD_initfirsttime(path_pdf_c, &pathpdfLength, &PDFMember);
}
void TUtil::SetJHUGenHiggsMassWidth(double MReso, double GaReso){
  const double GeV = 1./100.;
  MReso *= GeV; // GeV units in JHUGen
  GaReso *= GeV; // GeV units in JHUGen
  __modjhugenmela_MOD_sethiggsmasswidth(&MReso, &GaReso);
}
void TUtil::SetJHUGenDistinguishWWCouplings(bool doAllow){
  int iAllow = (doAllow ? 1 : 0);
  __modjhugenmela_MOD_setdistinguishwwcouplingsflag(&iAllow);
}
void TUtil::SetMCFMSpinZeroVVCouplings(bool useBSM, SpinZeroCouplings* Hcouplings, bool forceZZ){
  if (!useBSM){
    spinzerohiggs_anomcoupl_.AllowAnomalousCouplings = false;
    spinzerohiggs_anomcoupl_.distinguish_HWWcouplings = false;

    /***** REGULAR RESONANCE *****/
    //
    spinzerohiggs_anomcoupl_.cz_q1sq = 0;
    spinzerohiggs_anomcoupl_.Lambda_z11 = 100;
    spinzerohiggs_anomcoupl_.Lambda_z21 = 100;
    spinzerohiggs_anomcoupl_.Lambda_z31 = 100;
    spinzerohiggs_anomcoupl_.Lambda_z41 = 100;
    spinzerohiggs_anomcoupl_.cz_q2sq = 0;
    spinzerohiggs_anomcoupl_.Lambda_z12 = 100;
    spinzerohiggs_anomcoupl_.Lambda_z22 = 100;
    spinzerohiggs_anomcoupl_.Lambda_z32 = 100;
    spinzerohiggs_anomcoupl_.Lambda_z42 = 100;
    spinzerohiggs_anomcoupl_.cz_q12sq = 0;
    spinzerohiggs_anomcoupl_.Lambda_z10 = 100;
    spinzerohiggs_anomcoupl_.Lambda_z20 = 100;
    spinzerohiggs_anomcoupl_.Lambda_z30 = 100;
    spinzerohiggs_anomcoupl_.Lambda_z40 = 100;
    //
    spinzerohiggs_anomcoupl_.cw_q1sq = 0;
    spinzerohiggs_anomcoupl_.Lambda_w11 = 100;
    spinzerohiggs_anomcoupl_.Lambda_w21 = 100;
    spinzerohiggs_anomcoupl_.Lambda_w31 = 100;
    spinzerohiggs_anomcoupl_.Lambda_w41 = 100;
    spinzerohiggs_anomcoupl_.cw_q2sq = 0;
    spinzerohiggs_anomcoupl_.Lambda_w12 = 100;
    spinzerohiggs_anomcoupl_.Lambda_w22 = 100;
    spinzerohiggs_anomcoupl_.Lambda_w32 = 100;
    spinzerohiggs_anomcoupl_.Lambda_w42 = 100;
    spinzerohiggs_anomcoupl_.cw_q12sq = 0;
    spinzerohiggs_anomcoupl_.Lambda_w10 = 100;
    spinzerohiggs_anomcoupl_.Lambda_w20 = 100;
    spinzerohiggs_anomcoupl_.Lambda_w30 = 100;
    spinzerohiggs_anomcoupl_.Lambda_w40 = 100;
    //
    spinzerohiggs_anomcoupl_.ghz1[0] = 1; spinzerohiggs_anomcoupl_.ghz1[1] = 0;
    spinzerohiggs_anomcoupl_.ghw1[0] = 1; spinzerohiggs_anomcoupl_.ghw1[1] = 0;
    for (int im=0; im<2; im++){
      spinzerohiggs_anomcoupl_.ghz2[im] =  0;
      spinzerohiggs_anomcoupl_.ghz3[im] =  0;
      spinzerohiggs_anomcoupl_.ghz4[im] =  0;
      spinzerohiggs_anomcoupl_.ghz1_prime[im] = 0;
      spinzerohiggs_anomcoupl_.ghz1_prime2[im] = 0;
      spinzerohiggs_anomcoupl_.ghz1_prime3[im] = 0;
      spinzerohiggs_anomcoupl_.ghz1_prime4[im] = 0;
      spinzerohiggs_anomcoupl_.ghz1_prime5[im] = 0;
      spinzerohiggs_anomcoupl_.ghz1_prime6[im] = 0;
      spinzerohiggs_anomcoupl_.ghz1_prime7[im] = 0;
      spinzerohiggs_anomcoupl_.ghz2_prime[im] = 0;
      spinzerohiggs_anomcoupl_.ghz2_prime2[im] = 0;
      spinzerohiggs_anomcoupl_.ghz2_prime3[im] = 0;
      spinzerohiggs_anomcoupl_.ghz2_prime4[im] = 0;
      spinzerohiggs_anomcoupl_.ghz2_prime5[im] = 0;
      spinzerohiggs_anomcoupl_.ghz2_prime6[im] = 0;
      spinzerohiggs_anomcoupl_.ghz2_prime7[im] = 0;
      spinzerohiggs_anomcoupl_.ghz3_prime[im] = 0;
      spinzerohiggs_anomcoupl_.ghz3_prime2[im] = 0;
      spinzerohiggs_anomcoupl_.ghz3_prime3[im] = 0;
      spinzerohiggs_anomcoupl_.ghz3_prime4[im] = 0;
      spinzerohiggs_anomcoupl_.ghz3_prime5[im] = 0;
      spinzerohiggs_anomcoupl_.ghz3_prime6[im] = 0;
      spinzerohiggs_anomcoupl_.ghz3_prime7[im] = 0;
      spinzerohiggs_anomcoupl_.ghz4_prime[im] = 0;
      spinzerohiggs_anomcoupl_.ghz4_prime2[im] = 0;
      spinzerohiggs_anomcoupl_.ghz4_prime3[im] = 0;
      spinzerohiggs_anomcoupl_.ghz4_prime4[im] = 0;
      spinzerohiggs_anomcoupl_.ghz4_prime5[im] = 0;
      spinzerohiggs_anomcoupl_.ghz4_prime6[im] = 0;
      spinzerohiggs_anomcoupl_.ghz4_prime7[im] = 0;
      //
      spinzerohiggs_anomcoupl_.ghzgs1_prime2[im] = 0;
      spinzerohiggs_anomcoupl_.ghzgs2[im] = 0;
      spinzerohiggs_anomcoupl_.ghzgs3[im] = 0;
      spinzerohiggs_anomcoupl_.ghzgs4[im] = 0;
      spinzerohiggs_anomcoupl_.ghgsgs2[im] = 0;
      spinzerohiggs_anomcoupl_.ghgsgs3[im] = 0;
      spinzerohiggs_anomcoupl_.ghgsgs4[im] = 0;
      //
      spinzerohiggs_anomcoupl_.ghw2[im] =  0;
      spinzerohiggs_anomcoupl_.ghw3[im] =  0;
      spinzerohiggs_anomcoupl_.ghw4[im] =  0;
      spinzerohiggs_anomcoupl_.ghw1_prime[im] = 0;
      spinzerohiggs_anomcoupl_.ghw1_prime2[im] = 0;
      spinzerohiggs_anomcoupl_.ghw1_prime3[im] = 0;
      spinzerohiggs_anomcoupl_.ghw1_prime4[im] = 0;
      spinzerohiggs_anomcoupl_.ghw1_prime5[im] = 0;
      spinzerohiggs_anomcoupl_.ghw1_prime6[im] = 0;
      spinzerohiggs_anomcoupl_.ghw1_prime7[im] = 0;
      spinzerohiggs_anomcoupl_.ghw2_prime[im] = 0;
      spinzerohiggs_anomcoupl_.ghw2_prime2[im] = 0;
      spinzerohiggs_anomcoupl_.ghw2_prime3[im] = 0;
      spinzerohiggs_anomcoupl_.ghw2_prime4[im] = 0;
      spinzerohiggs_anomcoupl_.ghw2_prime5[im] = 0;
      spinzerohiggs_anomcoupl_.ghw2_prime6[im] = 0;
      spinzerohiggs_anomcoupl_.ghw2_prime7[im] = 0;
      spinzerohiggs_anomcoupl_.ghw3_prime[im] = 0;
      spinzerohiggs_anomcoupl_.ghw3_prime2[im] = 0;
      spinzerohiggs_anomcoupl_.ghw3_prime3[im] = 0;
      spinzerohiggs_anomcoupl_.ghw3_prime4[im] = 0;
      spinzerohiggs_anomcoupl_.ghw3_prime5[im] = 0;
      spinzerohiggs_anomcoupl_.ghw3_prime6[im] = 0;
      spinzerohiggs_anomcoupl_.ghw3_prime7[im] = 0;
      spinzerohiggs_anomcoupl_.ghw4_prime[im] = 0;
      spinzerohiggs_anomcoupl_.ghw4_prime2[im] = 0;
      spinzerohiggs_anomcoupl_.ghw4_prime3[im] = 0;
      spinzerohiggs_anomcoupl_.ghw4_prime4[im] = 0;
      spinzerohiggs_anomcoupl_.ghw4_prime5[im] = 0;
      spinzerohiggs_anomcoupl_.ghw4_prime6[im] = 0;
      spinzerohiggs_anomcoupl_.ghw4_prime7[im] = 0;
    }
    /***** END REGULAR RESONANCE *****/
    //
    /***** SECOND RESONANCE *****/
    //
    spinzerohiggs_anomcoupl_.c2z_q1sq = 0;
    spinzerohiggs_anomcoupl_.Lambda2_z11 = 100;
    spinzerohiggs_anomcoupl_.Lambda2_z21 = 100;
    spinzerohiggs_anomcoupl_.Lambda2_z31 = 100;
    spinzerohiggs_anomcoupl_.Lambda2_z41 = 100;
    spinzerohiggs_anomcoupl_.c2z_q2sq = 0;
    spinzerohiggs_anomcoupl_.Lambda2_z12 = 100;
    spinzerohiggs_anomcoupl_.Lambda2_z22 = 100;
    spinzerohiggs_anomcoupl_.Lambda2_z32 = 100;
    spinzerohiggs_anomcoupl_.Lambda2_z42 = 100;
    spinzerohiggs_anomcoupl_.c2z_q12sq = 0;
    spinzerohiggs_anomcoupl_.Lambda2_z10 = 100;
    spinzerohiggs_anomcoupl_.Lambda2_z20 = 100;
    spinzerohiggs_anomcoupl_.Lambda2_z30 = 100;
    spinzerohiggs_anomcoupl_.Lambda2_z40 = 100;
    //
    spinzerohiggs_anomcoupl_.c2w_q1sq = 0;
    spinzerohiggs_anomcoupl_.Lambda2_w11 = 100;
    spinzerohiggs_anomcoupl_.Lambda2_w21 = 100;
    spinzerohiggs_anomcoupl_.Lambda2_w31 = 100;
    spinzerohiggs_anomcoupl_.Lambda2_w41 = 100;
    spinzerohiggs_anomcoupl_.c2w_q2sq = 0;
    spinzerohiggs_anomcoupl_.Lambda2_w12 = 100;
    spinzerohiggs_anomcoupl_.Lambda2_w22 = 100;
    spinzerohiggs_anomcoupl_.Lambda2_w32 = 100;
    spinzerohiggs_anomcoupl_.Lambda2_w42 = 100;
    spinzerohiggs_anomcoupl_.c2w_q12sq = 0;
    spinzerohiggs_anomcoupl_.Lambda2_w10 = 100;
    spinzerohiggs_anomcoupl_.Lambda2_w20 = 100;
    spinzerohiggs_anomcoupl_.Lambda2_w30 = 100;
    spinzerohiggs_anomcoupl_.Lambda2_w40 = 100;
    //
    spinzerohiggs_anomcoupl_.gh2z1[0] = 0; spinzerohiggs_anomcoupl_.gh2z1[1] = 0;
    spinzerohiggs_anomcoupl_.gh2w1[0] = 0; spinzerohiggs_anomcoupl_.gh2w1[1] = 0;
    for (int im=0; im<2; im++){
      spinzerohiggs_anomcoupl_.gh2z2[im] =  0;
      spinzerohiggs_anomcoupl_.gh2z3[im] =  0;
      spinzerohiggs_anomcoupl_.gh2z4[im] =  0;
      spinzerohiggs_anomcoupl_.gh2z1_prime[im] = 0;
      spinzerohiggs_anomcoupl_.gh2z1_prime2[im] = 0;
      spinzerohiggs_anomcoupl_.gh2z1_prime3[im] = 0;
      spinzerohiggs_anomcoupl_.gh2z1_prime4[im] = 0;
      spinzerohiggs_anomcoupl_.gh2z1_prime5[im] = 0;
      spinzerohiggs_anomcoupl_.gh2z1_prime6[im] = 0;
      spinzerohiggs_anomcoupl_.gh2z1_prime7[im] = 0;
      spinzerohiggs_anomcoupl_.gh2z2_prime[im] = 0;
      spinzerohiggs_anomcoupl_.gh2z2_prime2[im] = 0;
      spinzerohiggs_anomcoupl_.gh2z2_prime3[im] = 0;
      spinzerohiggs_anomcoupl_.gh2z2_prime4[im] = 0;
      spinzerohiggs_anomcoupl_.gh2z2_prime5[im] = 0;
      spinzerohiggs_anomcoupl_.gh2z2_prime6[im] = 0;
      spinzerohiggs_anomcoupl_.gh2z2_prime7[im] = 0;
      spinzerohiggs_anomcoupl_.gh2z3_prime[im] = 0;
      spinzerohiggs_anomcoupl_.gh2z3_prime2[im] = 0;
      spinzerohiggs_anomcoupl_.gh2z3_prime3[im] = 0;
      spinzerohiggs_anomcoupl_.gh2z3_prime4[im] = 0;
      spinzerohiggs_anomcoupl_.gh2z3_prime5[im] = 0;
      spinzerohiggs_anomcoupl_.gh2z3_prime6[im] = 0;
      spinzerohiggs_anomcoupl_.gh2z3_prime7[im] = 0;
      spinzerohiggs_anomcoupl_.gh2z4_prime[im] = 0;
      spinzerohiggs_anomcoupl_.gh2z4_prime2[im] = 0;
      spinzerohiggs_anomcoupl_.gh2z4_prime3[im] = 0;
      spinzerohiggs_anomcoupl_.gh2z4_prime4[im] = 0;
      spinzerohiggs_anomcoupl_.gh2z4_prime5[im] = 0;
      spinzerohiggs_anomcoupl_.gh2z4_prime6[im] = 0;
      spinzerohiggs_anomcoupl_.gh2z4_prime7[im] = 0;
      //
      spinzerohiggs_anomcoupl_.gh2zgs1_prime2[im] = 0;
      spinzerohiggs_anomcoupl_.gh2zgs2[im] = 0;
      spinzerohiggs_anomcoupl_.gh2zgs3[im] = 0;
      spinzerohiggs_anomcoupl_.gh2zgs4[im] = 0;
      spinzerohiggs_anomcoupl_.gh2gsgs2[im] = 0;
      spinzerohiggs_anomcoupl_.gh2gsgs3[im] = 0;
      spinzerohiggs_anomcoupl_.gh2gsgs4[im] = 0;
      //
      spinzerohiggs_anomcoupl_.gh2w2[im] =  0;
      spinzerohiggs_anomcoupl_.gh2w3[im] =  0;
      spinzerohiggs_anomcoupl_.gh2w4[im] =  0;
      spinzerohiggs_anomcoupl_.gh2w1_prime[im] = 0;
      spinzerohiggs_anomcoupl_.gh2w1_prime2[im] = 0;
      spinzerohiggs_anomcoupl_.gh2w1_prime3[im] = 0;
      spinzerohiggs_anomcoupl_.gh2w1_prime4[im] = 0;
      spinzerohiggs_anomcoupl_.gh2w1_prime5[im] = 0;
      spinzerohiggs_anomcoupl_.gh2w1_prime6[im] = 0;
      spinzerohiggs_anomcoupl_.gh2w1_prime7[im] = 0;
      spinzerohiggs_anomcoupl_.gh2w2_prime[im] = 0;
      spinzerohiggs_anomcoupl_.gh2w2_prime2[im] = 0;
      spinzerohiggs_anomcoupl_.gh2w2_prime3[im] = 0;
      spinzerohiggs_anomcoupl_.gh2w2_prime4[im] = 0;
      spinzerohiggs_anomcoupl_.gh2w2_prime5[im] = 0;
      spinzerohiggs_anomcoupl_.gh2w2_prime6[im] = 0;
      spinzerohiggs_anomcoupl_.gh2w2_prime7[im] = 0;
      spinzerohiggs_anomcoupl_.gh2w3_prime[im] = 0;
      spinzerohiggs_anomcoupl_.gh2w3_prime2[im] = 0;
      spinzerohiggs_anomcoupl_.gh2w3_prime3[im] = 0;
      spinzerohiggs_anomcoupl_.gh2w3_prime4[im] = 0;
      spinzerohiggs_anomcoupl_.gh2w3_prime5[im] = 0;
      spinzerohiggs_anomcoupl_.gh2w3_prime6[im] = 0;
      spinzerohiggs_anomcoupl_.gh2w3_prime7[im] = 0;
      spinzerohiggs_anomcoupl_.gh2w4_prime[im] = 0;
      spinzerohiggs_anomcoupl_.gh2w4_prime2[im] = 0;
      spinzerohiggs_anomcoupl_.gh2w4_prime3[im] = 0;
      spinzerohiggs_anomcoupl_.gh2w4_prime4[im] = 0;
      spinzerohiggs_anomcoupl_.gh2w4_prime5[im] = 0;
      spinzerohiggs_anomcoupl_.gh2w4_prime6[im] = 0;
      spinzerohiggs_anomcoupl_.gh2w4_prime7[im] = 0;
    }
    /***** END SECOND RESONANCE *****/
  }
  else{
    spinzerohiggs_anomcoupl_.AllowAnomalousCouplings = true;
    spinzerohiggs_anomcoupl_.distinguish_HWWcouplings = (Hcouplings->separateWWZZcouplings && !forceZZ);

    /***** REGULAR RESONANCE *****/
    //
    spinzerohiggs_anomcoupl_.cz_q1sq = (Hcoupling->HzzCLambda_qsq)[0];
    spinzerohiggs_anomcoupl_.Lambda_z11 = (Hcoupling->HzzLambda_qsq)[0][0];
    spinzerohiggs_anomcoupl_.Lambda_z21 = (Hcoupling->HzzLambda_qsq)[1][0];
    spinzerohiggs_anomcoupl_.Lambda_z31 = (Hcoupling->HzzLambda_qsq)[2][0];
    spinzerohiggs_anomcoupl_.Lambda_z41 = (Hcoupling->HzzLambda_qsq)[3][0];
    spinzerohiggs_anomcoupl_.cz_q2sq = (Hcoupling->HzzCLambda_qsq)[1];
    spinzerohiggs_anomcoupl_.Lambda_z12 = (Hcoupling->HzzLambda_qsq)[0][1];
    spinzerohiggs_anomcoupl_.Lambda_z22 = (Hcoupling->HzzLambda_qsq)[1][1];
    spinzerohiggs_anomcoupl_.Lambda_z32 = (Hcoupling->HzzLambda_qsq)[2][1];
    spinzerohiggs_anomcoupl_.Lambda_z42 = (Hcoupling->HzzLambda_qsq)[3][1];
    spinzerohiggs_anomcoupl_.cz_q12sq = (Hcoupling->HzzCLambda_qsq)[2];
    spinzerohiggs_anomcoupl_.Lambda_z10 = (Hcoupling->HzzLambda_qsq)[0][2];
    spinzerohiggs_anomcoupl_.Lambda_z20 = (Hcoupling->HzzLambda_qsq)[1][2];
    spinzerohiggs_anomcoupl_.Lambda_z30 = (Hcoupling->HzzLambda_qsq)[2][2];
    spinzerohiggs_anomcoupl_.Lambda_z40 = (Hcoupling->HzzLambda_qsq)[3][2];
    //
    for (int im=0; im<2; im++){
      spinzerohiggs_anomcoupl_.ghz1[im] = (Hcoupling->Hzzcoupl)[0][im];
      spinzerohiggs_anomcoupl_.ghz2[im] = (Hcoupling->Hzzcoupl)[1][im];
      spinzerohiggs_anomcoupl_.ghz3[im] = (Hcoupling->Hzzcoupl)[2][im];
      spinzerohiggs_anomcoupl_.ghz4[im] = (Hcoupling->Hzzcoupl)[3][im];
      spinzerohiggs_anomcoupl_.ghz1_prime[im] = (Hcoupling->Hzzcoupl)[10][im];
      spinzerohiggs_anomcoupl_.ghz1_prime2[im] = (Hcoupling->Hzzcoupl)[11][im];
      spinzerohiggs_anomcoupl_.ghz1_prime3[im] = (Hcoupling->Hzzcoupl)[12][im];
      spinzerohiggs_anomcoupl_.ghz1_prime4[im] = (Hcoupling->Hzzcoupl)[13][im];
      spinzerohiggs_anomcoupl_.ghz1_prime5[im] = (Hcoupling->Hzzcoupl)[14][im];
      spinzerohiggs_anomcoupl_.ghz2_prime[im] = (Hcoupling->Hzzcoupl)[15][im];
      spinzerohiggs_anomcoupl_.ghz2_prime2[im] = (Hcoupling->Hzzcoupl)[16][im];
      spinzerohiggs_anomcoupl_.ghz2_prime3[im] = (Hcoupling->Hzzcoupl)[17][im];
      spinzerohiggs_anomcoupl_.ghz2_prime4[im] = (Hcoupling->Hzzcoupl)[18][im];
      spinzerohiggs_anomcoupl_.ghz2_prime5[im] = (Hcoupling->Hzzcoupl)[19][im];
      spinzerohiggs_anomcoupl_.ghz3_prime[im] = (Hcoupling->Hzzcoupl)[20][im];
      spinzerohiggs_anomcoupl_.ghz3_prime2[im] = (Hcoupling->Hzzcoupl)[21][im];
      spinzerohiggs_anomcoupl_.ghz3_prime3[im] = (Hcoupling->Hzzcoupl)[22][im];
      spinzerohiggs_anomcoupl_.ghz3_prime4[im] = (Hcoupling->Hzzcoupl)[23][im];
      spinzerohiggs_anomcoupl_.ghz3_prime5[im] = (Hcoupling->Hzzcoupl)[24][im];
      spinzerohiggs_anomcoupl_.ghz4_prime[im] = (Hcoupling->Hzzcoupl)[25][im];
      spinzerohiggs_anomcoupl_.ghz4_prime2[im] = (Hcoupling->Hzzcoupl)[26][im];
      spinzerohiggs_anomcoupl_.ghz4_prime3[im] = (Hcoupling->Hzzcoupl)[27][im];
      spinzerohiggs_anomcoupl_.ghz4_prime4[im] = (Hcoupling->Hzzcoupl)[28][im];
      spinzerohiggs_anomcoupl_.ghz4_prime5[im] = (Hcoupling->Hzzcoupl)[29][im];
      spinzerohiggs_anomcoupl_.ghz1_prime6[im] = (Hcoupling->Hzzcoupl)[31][im];
      spinzerohiggs_anomcoupl_.ghz1_prime7[im] = (Hcoupling->Hzzcoupl)[32][im];
      spinzerohiggs_anomcoupl_.ghz2_prime6[im] = (Hcoupling->Hzzcoupl)[33][im];
      spinzerohiggs_anomcoupl_.ghz2_prime7[im] = (Hcoupling->Hzzcoupl)[34][im];
      spinzerohiggs_anomcoupl_.ghz3_prime6[im] = (Hcoupling->Hzzcoupl)[35][im];
      spinzerohiggs_anomcoupl_.ghz3_prime7[im] = (Hcoupling->Hzzcoupl)[36][im];
      spinzerohiggs_anomcoupl_.ghz4_prime6[im] = (Hcoupling->Hzzcoupl)[37][im];
      spinzerohiggs_anomcoupl_.ghz4_prime7[im] = (Hcoupling->Hzzcoupl)[38][im];
      //
      spinzerohiggs_anomcoupl_.ghzgs1_prime2[im] = (Hcoupling->Hzzcoupl)[30][im];
      spinzerohiggs_anomcoupl_.ghzgs2[im] = (Hcoupling->Hzzcoupl)[4][im];
      spinzerohiggs_anomcoupl_.ghzgs3[im] = (Hcoupling->Hzzcoupl)[5][im];
      spinzerohiggs_anomcoupl_.ghzgs4[im] = (Hcoupling->Hzzcoupl)[6][im];
      spinzerohiggs_anomcoupl_.ghgsgs2[im] = (Hcoupling->Hzzcoupl)[7][im];
      spinzerohiggs_anomcoupl_.ghgsgs3[im] = (Hcoupling->Hzzcoupl)[8][im];
      spinzerohiggs_anomcoupl_.ghgsgs4[im] = (Hcoupling->Hzzcoupl)[9][im];
    }
    //
    if (spinzerohiggs_anomcoupl_.distinguish_HWWcouplings){
      //
      spinzerohiggs_anomcoupl_.cw_q1sq = (Hcoupling->HwwCLambda_qsq)[0];
      spinzerohiggs_anomcoupl_.Lambda_w11 = (Hcoupling->HwwLambda_qsq)[0][0];
      spinzerohiggs_anomcoupl_.Lambda_w21 = (Hcoupling->HwwLambda_qsq)[1][0];
      spinzerohiggs_anomcoupl_.Lambda_w31 = (Hcoupling->HwwLambda_qsq)[2][0];
      spinzerohiggs_anomcoupl_.Lambda_w41 = (Hcoupling->HwwLambda_qsq)[3][0];
      spinzerohiggs_anomcoupl_.cw_q2sq = (Hcoupling->HwwCLambda_qsq)[1];
      spinzerohiggs_anomcoupl_.Lambda_w12 = (Hcoupling->HwwLambda_qsq)[0][1];
      spinzerohiggs_anomcoupl_.Lambda_w22 = (Hcoupling->HwwLambda_qsq)[1][1];
      spinzerohiggs_anomcoupl_.Lambda_w32 = (Hcoupling->HwwLambda_qsq)[2][1];
      spinzerohiggs_anomcoupl_.Lambda_w42 = (Hcoupling->HwwLambda_qsq)[3][1];
      spinzerohiggs_anomcoupl_.cw_q12sq = (Hcoupling->HwwCLambda_qsq)[2];
      spinzerohiggs_anomcoupl_.Lambda_w10 = (Hcoupling->HwwLambda_qsq)[0][2];
      spinzerohiggs_anomcoupl_.Lambda_w20 = (Hcoupling->HwwLambda_qsq)[1][2];
      spinzerohiggs_anomcoupl_.Lambda_w30 = (Hcoupling->HwwLambda_qsq)[2][2];
      spinzerohiggs_anomcoupl_.Lambda_w40 = (Hcoupling->HwwLambda_qsq)[3][2];
      //
      for (int im=0; im<2; im++){
        spinzerohiggs_anomcoupl_.ghw1[im] = (Hcouplings->Hwwcoupl)[0][im];
        spinzerohiggs_anomcoupl_.ghw2[im] = (Hcouplings->Hwwcoupl)[1][im];
        spinzerohiggs_anomcoupl_.ghw3[im] = (Hcouplings->Hwwcoupl)[2][im];
        spinzerohiggs_anomcoupl_.ghw4[im] = (Hcouplings->Hwwcoupl)[3][im];
        spinzerohiggs_anomcoupl_.ghw1_prime[im] = (Hcouplings->Hwwcoupl)[10][im];
        spinzerohiggs_anomcoupl_.ghw1_prime2[im] = (Hcouplings->Hwwcoupl)[11][im];
        spinzerohiggs_anomcoupl_.ghw1_prime3[im] = (Hcouplings->Hwwcoupl)[12][im];
        spinzerohiggs_anomcoupl_.ghw1_prime4[im] = (Hcouplings->Hwwcoupl)[13][im];
        spinzerohiggs_anomcoupl_.ghw1_prime5[im] = (Hcouplings->Hwwcoupl)[14][im];
        spinzerohiggs_anomcoupl_.ghw2_prime[im] = (Hcouplings->Hwwcoupl)[15][im];
        spinzerohiggs_anomcoupl_.ghw2_prime2[im] = (Hcouplings->Hwwcoupl)[16][im];
        spinzerohiggs_anomcoupl_.ghw2_prime3[im] = (Hcouplings->Hwwcoupl)[17][im];
        spinzerohiggs_anomcoupl_.ghw2_prime4[im] = (Hcouplings->Hwwcoupl)[18][im];
        spinzerohiggs_anomcoupl_.ghw2_prime5[im] = (Hcouplings->Hwwcoupl)[19][im];
        spinzerohiggs_anomcoupl_.ghw3_prime[im] = (Hcouplings->Hwwcoupl)[20][im];
        spinzerohiggs_anomcoupl_.ghw3_prime2[im] = (Hcouplings->Hwwcoupl)[21][im];
        spinzerohiggs_anomcoupl_.ghw3_prime3[im] = (Hcouplings->Hwwcoupl)[22][im];
        spinzerohiggs_anomcoupl_.ghw3_prime4[im] = (Hcouplings->Hwwcoupl)[23][im];
        spinzerohiggs_anomcoupl_.ghw3_prime5[im] = (Hcouplings->Hwwcoupl)[24][im];
        spinzerohiggs_anomcoupl_.ghw4_prime[im] = (Hcouplings->Hwwcoupl)[25][im];
        spinzerohiggs_anomcoupl_.ghw4_prime2[im] = (Hcouplings->Hwwcoupl)[26][im];
        spinzerohiggs_anomcoupl_.ghw4_prime3[im] = (Hcouplings->Hwwcoupl)[27][im];
        spinzerohiggs_anomcoupl_.ghw4_prime4[im] = (Hcouplings->Hwwcoupl)[28][im];
        spinzerohiggs_anomcoupl_.ghw4_prime5[im] = (Hcouplings->Hwwcoupl)[29][im];
        spinzerohiggs_anomcoupl_.ghw1_prime6[im] = (Hcouplings->Hwwcoupl)[31][im];
        spinzerohiggs_anomcoupl_.ghw1_prime7[im] = (Hcouplings->Hwwcoupl)[32][im];
        spinzerohiggs_anomcoupl_.ghw2_prime6[im] = (Hcouplings->Hwwcoupl)[33][im];
        spinzerohiggs_anomcoupl_.ghw2_prime7[im] = (Hcouplings->Hwwcoupl)[34][im];
        spinzerohiggs_anomcoupl_.ghw3_prime6[im] = (Hcouplings->Hwwcoupl)[35][im];
        spinzerohiggs_anomcoupl_.ghw3_prime7[im] = (Hcouplings->Hwwcoupl)[36][im];
        spinzerohiggs_anomcoupl_.ghw4_prime6[im] = (Hcouplings->Hwwcoupl)[37][im];
        spinzerohiggs_anomcoupl_.ghw4_prime7[im] = (Hcouplings->Hwwcoupl)[38][im];
      }
    }
    else{
      //
      spinzerohiggs_anomcoupl_.cw_q1sq = (Hcoupling->HzzCLambda_qsq)[0];
      spinzerohiggs_anomcoupl_.Lambda_w11 = (Hcoupling->HzzLambda_qsq)[0][0];
      spinzerohiggs_anomcoupl_.Lambda_w21 = (Hcoupling->HzzLambda_qsq)[1][0];
      spinzerohiggs_anomcoupl_.Lambda_w31 = (Hcoupling->HzzLambda_qsq)[2][0];
      spinzerohiggs_anomcoupl_.Lambda_w41 = (Hcoupling->HzzLambda_qsq)[3][0];
      spinzerohiggs_anomcoupl_.cw_q2sq = (Hcoupling->HzzCLambda_qsq)[1];
      spinzerohiggs_anomcoupl_.Lambda_w12 = (Hcoupling->HzzLambda_qsq)[0][1];
      spinzerohiggs_anomcoupl_.Lambda_w22 = (Hcoupling->HzzLambda_qsq)[1][1];
      spinzerohiggs_anomcoupl_.Lambda_w32 = (Hcoupling->HzzLambda_qsq)[2][1];
      spinzerohiggs_anomcoupl_.Lambda_w42 = (Hcoupling->HzzLambda_qsq)[3][1];
      spinzerohiggs_anomcoupl_.cw_q12sq = (Hcoupling->HzzCLambda_qsq)[2];
      spinzerohiggs_anomcoupl_.Lambda_w10 = (Hcoupling->HzzLambda_qsq)[0][2];
      spinzerohiggs_anomcoupl_.Lambda_w20 = (Hcoupling->HzzLambda_qsq)[1][2];
      spinzerohiggs_anomcoupl_.Lambda_w30 = (Hcoupling->HzzLambda_qsq)[2][2];
      spinzerohiggs_anomcoupl_.Lambda_w40 = (Hcoupling->HzzLambda_qsq)[3][2];
      //
      for (int im=0; im<2; im++){
        spinzerohiggs_anomcoupl_.ghw1[im] = (Hcouplings->Hzzcoupl)[0][im];
        spinzerohiggs_anomcoupl_.ghw2[im] = (Hcouplings->Hzzcoupl)[1][im];
        spinzerohiggs_anomcoupl_.ghw3[im] = (Hcouplings->Hzzcoupl)[2][im];
        spinzerohiggs_anomcoupl_.ghw4[im] = (Hcouplings->Hzzcoupl)[3][im];
        spinzerohiggs_anomcoupl_.ghw1_prime[im] = (Hcouplings->Hzzcoupl)[10][im];
        spinzerohiggs_anomcoupl_.ghw1_prime2[im] = (Hcouplings->Hzzcoupl)[11][im];
        spinzerohiggs_anomcoupl_.ghw1_prime3[im] = (Hcouplings->Hzzcoupl)[12][im];
        spinzerohiggs_anomcoupl_.ghw1_prime4[im] = (Hcouplings->Hzzcoupl)[13][im];
        spinzerohiggs_anomcoupl_.ghw1_prime5[im] = (Hcouplings->Hzzcoupl)[14][im];
        spinzerohiggs_anomcoupl_.ghw2_prime[im] = (Hcouplings->Hzzcoupl)[15][im];
        spinzerohiggs_anomcoupl_.ghw2_prime2[im] = (Hcouplings->Hzzcoupl)[16][im];
        spinzerohiggs_anomcoupl_.ghw2_prime3[im] = (Hcouplings->Hzzcoupl)[17][im];
        spinzerohiggs_anomcoupl_.ghw2_prime4[im] = (Hcouplings->Hzzcoupl)[18][im];
        spinzerohiggs_anomcoupl_.ghw2_prime5[im] = (Hcouplings->Hzzcoupl)[19][im];
        spinzerohiggs_anomcoupl_.ghw3_prime[im] = (Hcouplings->Hzzcoupl)[20][im];
        spinzerohiggs_anomcoupl_.ghw3_prime2[im] = (Hcouplings->Hzzcoupl)[21][im];
        spinzerohiggs_anomcoupl_.ghw3_prime3[im] = (Hcouplings->Hzzcoupl)[22][im];
        spinzerohiggs_anomcoupl_.ghw3_prime4[im] = (Hcouplings->Hzzcoupl)[23][im];
        spinzerohiggs_anomcoupl_.ghw3_prime5[im] = (Hcouplings->Hzzcoupl)[24][im];
        spinzerohiggs_anomcoupl_.ghw4_prime[im] = (Hcouplings->Hzzcoupl)[25][im];
        spinzerohiggs_anomcoupl_.ghw4_prime2[im] = (Hcouplings->Hzzcoupl)[26][im];
        spinzerohiggs_anomcoupl_.ghw4_prime3[im] = (Hcouplings->Hzzcoupl)[27][im];
        spinzerohiggs_anomcoupl_.ghw4_prime4[im] = (Hcouplings->Hzzcoupl)[28][im];
        spinzerohiggs_anomcoupl_.ghw4_prime5[im] = (Hcouplings->Hzzcoupl)[29][im];
        spinzerohiggs_anomcoupl_.ghw1_prime6[im] = (Hcouplings->Hzzcoupl)[31][im];
        spinzerohiggs_anomcoupl_.ghw1_prime7[im] = (Hcouplings->Hzzcoupl)[32][im];
        spinzerohiggs_anomcoupl_.ghw2_prime6[im] = (Hcouplings->Hzzcoupl)[33][im];
        spinzerohiggs_anomcoupl_.ghw2_prime7[im] = (Hcouplings->Hzzcoupl)[34][im];
        spinzerohiggs_anomcoupl_.ghw3_prime6[im] = (Hcouplings->Hzzcoupl)[35][im];
        spinzerohiggs_anomcoupl_.ghw3_prime7[im] = (Hcouplings->Hzzcoupl)[36][im];
        spinzerohiggs_anomcoupl_.ghw4_prime6[im] = (Hcouplings->Hzzcoupl)[37][im];
        spinzerohiggs_anomcoupl_.ghw4_prime7[im] = (Hcouplings->Hzzcoupl)[38][im];
      }
    }
    /***** END REGULAR RESONANCE *****/
    //
    /***** SECOND RESONANCE *****/
    //
    spinzerohiggs_anomcoupl_.c2z_q1sq = (Hcoupling->H2zzCLambda_qsq)[0];
    spinzerohiggs_anomcoupl_.Lambda2_z11 = (Hcoupling->H2zzLambda_qsq)[0][0];
    spinzerohiggs_anomcoupl_.Lambda2_z21 = (Hcoupling->H2zzLambda_qsq)[1][0];
    spinzerohiggs_anomcoupl_.Lambda2_z31 = (Hcoupling->H2zzLambda_qsq)[2][0];
    spinzerohiggs_anomcoupl_.Lambda2_z41 = (Hcoupling->H2zzLambda_qsq)[3][0];
    spinzerohiggs_anomcoupl_.c2z_q2sq = (Hcoupling->H2zzCLambda_qsq)[1];
    spinzerohiggs_anomcoupl_.Lambda2_z12 = (Hcoupling->H2zzLambda_qsq)[0][1];
    spinzerohiggs_anomcoupl_.Lambda2_z22 = (Hcoupling->H2zzLambda_qsq)[1][1];
    spinzerohiggs_anomcoupl_.Lambda2_z32 = (Hcoupling->H2zzLambda_qsq)[2][1];
    spinzerohiggs_anomcoupl_.Lambda2_z42 = (Hcoupling->H2zzLambda_qsq)[3][1];
    spinzerohiggs_anomcoupl_.c2z_q12sq = (Hcoupling->H2zzCLambda_qsq)[2];
    spinzerohiggs_anomcoupl_.Lambda2_z10 = (Hcoupling->H2zzLambda_qsq)[0][2];
    spinzerohiggs_anomcoupl_.Lambda2_z20 = (Hcoupling->H2zzLambda_qsq)[1][2];
    spinzerohiggs_anomcoupl_.Lambda2_z30 = (Hcoupling->H2zzLambda_qsq)[2][2];
    spinzerohiggs_anomcoupl_.Lambda2_z40 = (Hcoupling->H2zzLambda_qsq)[3][2];
    //
    for (int im=0; im<2; im++){
      spinzerohiggs_anomcoupl_.gh2z1[im] = (Hcoupling->H2zzcoupl)[0][im];
      spinzerohiggs_anomcoupl_.gh2z2[im] = (Hcoupling->H2zzcoupl)[1][im];
      spinzerohiggs_anomcoupl_.gh2z3[im] = (Hcoupling->H2zzcoupl)[2][im];
      spinzerohiggs_anomcoupl_.gh2z4[im] = (Hcoupling->H2zzcoupl)[3][im];
      spinzerohiggs_anomcoupl_.gh2z1_prime[im] = (Hcoupling->H2zzcoupl)[10][im];
      spinzerohiggs_anomcoupl_.gh2z1_prime2[im] = (Hcoupling->H2zzcoupl)[11][im];
      spinzerohiggs_anomcoupl_.gh2z1_prime3[im] = (Hcoupling->H2zzcoupl)[12][im];
      spinzerohiggs_anomcoupl_.gh2z1_prime4[im] = (Hcoupling->H2zzcoupl)[13][im];
      spinzerohiggs_anomcoupl_.gh2z1_prime5[im] = (Hcoupling->H2zzcoupl)[14][im];
      spinzerohiggs_anomcoupl_.gh2z2_prime[im] = (Hcoupling->H2zzcoupl)[15][im];
      spinzerohiggs_anomcoupl_.gh2z2_prime2[im] = (Hcoupling->H2zzcoupl)[16][im];
      spinzerohiggs_anomcoupl_.gh2z2_prime3[im] = (Hcoupling->H2zzcoupl)[17][im];
      spinzerohiggs_anomcoupl_.gh2z2_prime4[im] = (Hcoupling->H2zzcoupl)[18][im];
      spinzerohiggs_anomcoupl_.gh2z2_prime5[im] = (Hcoupling->H2zzcoupl)[19][im];
      spinzerohiggs_anomcoupl_.gh2z3_prime[im] = (Hcoupling->H2zzcoupl)[20][im];
      spinzerohiggs_anomcoupl_.gh2z3_prime2[im] = (Hcoupling->H2zzcoupl)[21][im];
      spinzerohiggs_anomcoupl_.gh2z3_prime3[im] = (Hcoupling->H2zzcoupl)[22][im];
      spinzerohiggs_anomcoupl_.gh2z3_prime4[im] = (Hcoupling->H2zzcoupl)[23][im];
      spinzerohiggs_anomcoupl_.gh2z3_prime5[im] = (Hcoupling->H2zzcoupl)[24][im];
      spinzerohiggs_anomcoupl_.gh2z4_prime[im] = (Hcoupling->H2zzcoupl)[25][im];
      spinzerohiggs_anomcoupl_.gh2z4_prime2[im] = (Hcoupling->H2zzcoupl)[26][im];
      spinzerohiggs_anomcoupl_.gh2z4_prime3[im] = (Hcoupling->H2zzcoupl)[27][im];
      spinzerohiggs_anomcoupl_.gh2z4_prime4[im] = (Hcoupling->H2zzcoupl)[28][im];
      spinzerohiggs_anomcoupl_.gh2z4_prime5[im] = (Hcoupling->H2zzcoupl)[29][im];
      spinzerohiggs_anomcoupl_.gh2z1_prime6[im] = (Hcoupling->H2zzcoupl)[31][im];
      spinzerohiggs_anomcoupl_.gh2z1_prime7[im] = (Hcoupling->H2zzcoupl)[32][im];
      spinzerohiggs_anomcoupl_.gh2z2_prime6[im] = (Hcoupling->H2zzcoupl)[33][im];
      spinzerohiggs_anomcoupl_.gh2z2_prime7[im] = (Hcoupling->H2zzcoupl)[34][im];
      spinzerohiggs_anomcoupl_.gh2z3_prime6[im] = (Hcoupling->H2zzcoupl)[35][im];
      spinzerohiggs_anomcoupl_.gh2z3_prime7[im] = (Hcoupling->H2zzcoupl)[36][im];
      spinzerohiggs_anomcoupl_.gh2z4_prime6[im] = (Hcoupling->H2zzcoupl)[37][im];
      spinzerohiggs_anomcoupl_.gh2z4_prime7[im] = (Hcoupling->H2zzcoupl)[38][im];
      //
      spinzerohiggs_anomcoupl_.gh2zgs1_prime2[im] = (Hcoupling->H2zzcoupl)[30][im];
      spinzerohiggs_anomcoupl_.gh2zgs2[im] = (Hcoupling->H2zzcoupl)[4][im];
      spinzerohiggs_anomcoupl_.gh2zgs3[im] = (Hcoupling->H2zzcoupl)[5][im];
      spinzerohiggs_anomcoupl_.gh2zgs4[im] = (Hcoupling->H2zzcoupl)[6][im];
      spinzerohiggs_anomcoupl_.gh2gsgs2[im] = (Hcoupling->H2zzcoupl)[7][im];
      spinzerohiggs_anomcoupl_.gh2gsgs3[im] = (Hcoupling->H2zzcoupl)[8][im];
      spinzerohiggs_anomcoupl_.gh2gsgs4[im] = (Hcoupling->H2zzcoupl)[9][im];
    }
    //
    if (spinzerohiggs_anomcoupl_.distinguish_HWWcouplings){
      //
      spinzerohiggs_anomcoupl_.c2w_q1sq = (Hcoupling->H2wwCLambda_qsq)[0];
      spinzerohiggs_anomcoupl_.Lambda2_w11 = (Hcoupling->H2wwLambda_qsq)[0][0];
      spinzerohiggs_anomcoupl_.Lambda2_w21 = (Hcoupling->H2wwLambda_qsq)[1][0];
      spinzerohiggs_anomcoupl_.Lambda2_w31 = (Hcoupling->H2wwLambda_qsq)[2][0];
      spinzerohiggs_anomcoupl_.Lambda2_w41 = (Hcoupling->H2wwLambda_qsq)[3][0];
      spinzerohiggs_anomcoupl_.c2w_q2sq = (Hcoupling->H2wwCLambda_qsq)[1];
      spinzerohiggs_anomcoupl_.Lambda2_w12 = (Hcoupling->H2wwLambda_qsq)[0][1];
      spinzerohiggs_anomcoupl_.Lambda2_w22 = (Hcoupling->H2wwLambda_qsq)[1][1];
      spinzerohiggs_anomcoupl_.Lambda2_w32 = (Hcoupling->H2wwLambda_qsq)[2][1];
      spinzerohiggs_anomcoupl_.Lambda2_w42 = (Hcoupling->H2wwLambda_qsq)[3][1];
      spinzerohiggs_anomcoupl_.c2w_q12sq = (Hcoupling->H2wwCLambda_qsq)[2];
      spinzerohiggs_anomcoupl_.Lambda2_w10 = (Hcoupling->H2wwLambda_qsq)[0][2];
      spinzerohiggs_anomcoupl_.Lambda2_w20 = (Hcoupling->H2wwLambda_qsq)[1][2];
      spinzerohiggs_anomcoupl_.Lambda2_w30 = (Hcoupling->H2wwLambda_qsq)[2][2];
      spinzerohiggs_anomcoupl_.Lambda2_w40 = (Hcoupling->H2wwLambda_qsq)[3][2];
      //
      for (int im=0; im<2; im++){
        spinzerohiggs_anomcoupl_.gh2w1[im] = (Hcouplings->H2wwcoupl)[0][im];
        spinzerohiggs_anomcoupl_.gh2w2[im] = (Hcouplings->H2wwcoupl)[1][im];
        spinzerohiggs_anomcoupl_.gh2w3[im] = (Hcouplings->H2wwcoupl)[2][im];
        spinzerohiggs_anomcoupl_.gh2w4[im] = (Hcouplings->H2wwcoupl)[3][im];
        spinzerohiggs_anomcoupl_.gh2w1_prime[im] = (Hcouplings->H2wwcoupl)[10][im];
        spinzerohiggs_anomcoupl_.gh2w1_prime2[im] = (Hcouplings->H2wwcoupl)[11][im];
        spinzerohiggs_anomcoupl_.gh2w1_prime3[im] = (Hcouplings->H2wwcoupl)[12][im];
        spinzerohiggs_anomcoupl_.gh2w1_prime4[im] = (Hcouplings->H2wwcoupl)[13][im];
        spinzerohiggs_anomcoupl_.gh2w1_prime5[im] = (Hcouplings->H2wwcoupl)[14][im];
        spinzerohiggs_anomcoupl_.gh2w2_prime[im] = (Hcouplings->H2wwcoupl)[15][im];
        spinzerohiggs_anomcoupl_.gh2w2_prime2[im] = (Hcouplings->H2wwcoupl)[16][im];
        spinzerohiggs_anomcoupl_.gh2w2_prime3[im] = (Hcouplings->H2wwcoupl)[17][im];
        spinzerohiggs_anomcoupl_.gh2w2_prime4[im] = (Hcouplings->H2wwcoupl)[18][im];
        spinzerohiggs_anomcoupl_.gh2w2_prime5[im] = (Hcouplings->H2wwcoupl)[19][im];
        spinzerohiggs_anomcoupl_.gh2w3_prime[im] = (Hcouplings->H2wwcoupl)[20][im];
        spinzerohiggs_anomcoupl_.gh2w3_prime2[im] = (Hcouplings->H2wwcoupl)[21][im];
        spinzerohiggs_anomcoupl_.gh2w3_prime3[im] = (Hcouplings->H2wwcoupl)[22][im];
        spinzerohiggs_anomcoupl_.gh2w3_prime4[im] = (Hcouplings->H2wwcoupl)[23][im];
        spinzerohiggs_anomcoupl_.gh2w3_prime5[im] = (Hcouplings->H2wwcoupl)[24][im];
        spinzerohiggs_anomcoupl_.gh2w4_prime[im] = (Hcouplings->H2wwcoupl)[25][im];
        spinzerohiggs_anomcoupl_.gh2w4_prime2[im] = (Hcouplings->H2wwcoupl)[26][im];
        spinzerohiggs_anomcoupl_.gh2w4_prime3[im] = (Hcouplings->H2wwcoupl)[27][im];
        spinzerohiggs_anomcoupl_.gh2w4_prime4[im] = (Hcouplings->H2wwcoupl)[28][im];
        spinzerohiggs_anomcoupl_.gh2w4_prime5[im] = (Hcouplings->H2wwcoupl)[29][im];
        spinzerohiggs_anomcoupl_.gh2w1_prime6[im] = (Hcouplings->H2wwcoupl)[31][im];
        spinzerohiggs_anomcoupl_.gh2w1_prime7[im] = (Hcouplings->H2wwcoupl)[32][im];
        spinzerohiggs_anomcoupl_.gh2w2_prime6[im] = (Hcouplings->H2wwcoupl)[33][im];
        spinzerohiggs_anomcoupl_.gh2w2_prime7[im] = (Hcouplings->H2wwcoupl)[34][im];
        spinzerohiggs_anomcoupl_.gh2w3_prime6[im] = (Hcouplings->H2wwcoupl)[35][im];
        spinzerohiggs_anomcoupl_.gh2w3_prime7[im] = (Hcouplings->H2wwcoupl)[36][im];
        spinzerohiggs_anomcoupl_.gh2w4_prime6[im] = (Hcouplings->H2wwcoupl)[37][im];
        spinzerohiggs_anomcoupl_.gh2w4_prime7[im] = (Hcouplings->H2wwcoupl)[38][im];
      }
    }
    else{
      //
      spinzerohiggs_anomcoupl_.c2w_q1sq = (Hcoupling->H2zzCLambda_qsq)[0];
      spinzerohiggs_anomcoupl_.Lambda2_w11 = (Hcoupling->H2zzLambda_qsq)[0][0];
      spinzerohiggs_anomcoupl_.Lambda2_w21 = (Hcoupling->H2zzLambda_qsq)[1][0];
      spinzerohiggs_anomcoupl_.Lambda2_w31 = (Hcoupling->H2zzLambda_qsq)[2][0];
      spinzerohiggs_anomcoupl_.Lambda2_w41 = (Hcoupling->H2zzLambda_qsq)[3][0];
      spinzerohiggs_anomcoupl_.c2w_q2sq = (Hcoupling->H2zzCLambda_qsq)[1];
      spinzerohiggs_anomcoupl_.Lambda2_w12 = (Hcoupling->H2zzLambda_qsq)[0][1];
      spinzerohiggs_anomcoupl_.Lambda2_w22 = (Hcoupling->H2zzLambda_qsq)[1][1];
      spinzerohiggs_anomcoupl_.Lambda2_w32 = (Hcoupling->H2zzLambda_qsq)[2][1];
      spinzerohiggs_anomcoupl_.Lambda2_w42 = (Hcoupling->H2zzLambda_qsq)[3][1];
      spinzerohiggs_anomcoupl_.c2w_q12sq = (Hcoupling->H2zzCLambda_qsq)[2];
      spinzerohiggs_anomcoupl_.Lambda2_w10 = (Hcoupling->H2zzLambda_qsq)[0][2];
      spinzerohiggs_anomcoupl_.Lambda2_w20 = (Hcoupling->H2zzLambda_qsq)[1][2];
      spinzerohiggs_anomcoupl_.Lambda2_w30 = (Hcoupling->H2zzLambda_qsq)[2][2];
      spinzerohiggs_anomcoupl_.Lambda2_w40 = (Hcoupling->H2zzLambda_qsq)[3][2];
      //
      for (int im=0; im<2; im++){
        spinzerohiggs_anomcoupl_.gh2w1[im] = (Hcouplings->H2zzcoupl)[0][im];
        spinzerohiggs_anomcoupl_.gh2w2[im] = (Hcouplings->H2zzcoupl)[1][im];
        spinzerohiggs_anomcoupl_.gh2w3[im] = (Hcouplings->H2zzcoupl)[2][im];
        spinzerohiggs_anomcoupl_.gh2w4[im] = (Hcouplings->H2zzcoupl)[3][im];
        spinzerohiggs_anomcoupl_.gh2w1_prime[im] = (Hcouplings->H2zzcoupl)[10][im];
        spinzerohiggs_anomcoupl_.gh2w1_prime2[im] = (Hcouplings->H2zzcoupl)[11][im];
        spinzerohiggs_anomcoupl_.gh2w1_prime3[im] = (Hcouplings->H2zzcoupl)[12][im];
        spinzerohiggs_anomcoupl_.gh2w1_prime4[im] = (Hcouplings->H2zzcoupl)[13][im];
        spinzerohiggs_anomcoupl_.gh2w1_prime5[im] = (Hcouplings->H2zzcoupl)[14][im];
        spinzerohiggs_anomcoupl_.gh2w2_prime[im] = (Hcouplings->H2zzcoupl)[15][im];
        spinzerohiggs_anomcoupl_.gh2w2_prime2[im] = (Hcouplings->H2zzcoupl)[16][im];
        spinzerohiggs_anomcoupl_.gh2w2_prime3[im] = (Hcouplings->H2zzcoupl)[17][im];
        spinzerohiggs_anomcoupl_.gh2w2_prime4[im] = (Hcouplings->H2zzcoupl)[18][im];
        spinzerohiggs_anomcoupl_.gh2w2_prime5[im] = (Hcouplings->H2zzcoupl)[19][im];
        spinzerohiggs_anomcoupl_.gh2w3_prime[im] = (Hcouplings->H2zzcoupl)[20][im];
        spinzerohiggs_anomcoupl_.gh2w3_prime2[im] = (Hcouplings->H2zzcoupl)[21][im];
        spinzerohiggs_anomcoupl_.gh2w3_prime3[im] = (Hcouplings->H2zzcoupl)[22][im];
        spinzerohiggs_anomcoupl_.gh2w3_prime4[im] = (Hcouplings->H2zzcoupl)[23][im];
        spinzerohiggs_anomcoupl_.gh2w3_prime5[im] = (Hcouplings->H2zzcoupl)[24][im];
        spinzerohiggs_anomcoupl_.gh2w4_prime[im] = (Hcouplings->H2zzcoupl)[25][im];
        spinzerohiggs_anomcoupl_.gh2w4_prime2[im] = (Hcouplings->H2zzcoupl)[26][im];
        spinzerohiggs_anomcoupl_.gh2w4_prime3[im] = (Hcouplings->H2zzcoupl)[27][im];
        spinzerohiggs_anomcoupl_.gh2w4_prime4[im] = (Hcouplings->H2zzcoupl)[28][im];
        spinzerohiggs_anomcoupl_.gh2w4_prime5[im] = (Hcouplings->H2zzcoupl)[29][im];
        spinzerohiggs_anomcoupl_.gh2w1_prime6[im] = (Hcouplings->H2zzcoupl)[31][im];
        spinzerohiggs_anomcoupl_.gh2w1_prime7[im] = (Hcouplings->H2zzcoupl)[32][im];
        spinzerohiggs_anomcoupl_.gh2w2_prime6[im] = (Hcouplings->H2zzcoupl)[33][im];
        spinzerohiggs_anomcoupl_.gh2w2_prime7[im] = (Hcouplings->H2zzcoupl)[34][im];
        spinzerohiggs_anomcoupl_.gh2w3_prime6[im] = (Hcouplings->H2zzcoupl)[35][im];
        spinzerohiggs_anomcoupl_.gh2w3_prime7[im] = (Hcouplings->H2zzcoupl)[36][im];
        spinzerohiggs_anomcoupl_.gh2w4_prime6[im] = (Hcouplings->H2zzcoupl)[37][im];
        spinzerohiggs_anomcoupl_.gh2w4_prime7[im] = (Hcouplings->H2zzcoupl)[38][im];
      }
    }
    /***** END SECOND RESONANCE *****/
  }
}
void TUtil::SetJHUGenSpinZeroVVCouplings(double Hvvcoupl[SIZE_HVV][2], int Hvvcoupl_cqsq[3], double HvvLambda_qsq[4][3], bool useWWcoupl){
  const double GeV = 1./100.;
  int iWWcoupl = (useWWcoupl ? 1 : 0);
  for (int c=0; c<4; c++){ for (int k=0; k<3; k++) HvvLambda_qsq[c][k] *= GeV; } // GeV units in JHUGen
  __modjhugenmela_MOD_setspinzerovvcouplings(Hvvcoupl, Hvvcoupl_cqsq, HvvLambda_qsq, &iWWcoupl);
}
void TUtil::SetJHUGenSpinZeroGGCouplings(double Hggcoupl[SIZE_HGG][2]){ __modjhugenmela_MOD_setspinzeroggcouplings(Hggcoupl); }
void TUtil::SetJHUGenSpinZeroQQCouplings(double Hqqcoupl[SIZE_HQQ][2]){ __modjhugenmela_MOD_setspinzeroqqcouplings(Hqqcoupl); }
void TUtil::SetJHUGenSpinOneCouplings(double Zqqcoupl[SIZE_ZQQ][2], double Zvvcoupl[SIZE_ZVV][2]){ __modjhugenmela_MOD_setspinonecouplings(Zqqcoupl, Zvvcoupl); }
void TUtil::SetJHUGenSpinTwoCouplings(double Gacoupl[SIZE_GGG][2], double Gbcoupl[SIZE_GVV][2], double qLeftRightcoupl[SIZE_GQQ][2]){ __modjhugenmela_MOD_setspintwocouplings(Gacoupl, Gbcoupl, qLeftRightcoupl); }

bool TUtil::MCFM_masscuts(double s[][mxpart], TVar::Process process){
  double minZmassSqr=10*10;
  if (
    process==TVar::bkgZZ
    &&
    (s[2][3]<minZmassSqr || s[4][5]<minZmassSqr)
    ) return true;
  return false;
}
bool TUtil::MCFM_smalls(double s[][mxpart], int npart){

  // Reject event if any s(i,j) is too small
  // cutoff is defined in technical.Dat

  if (
    npart == 3 &&
    (
    (-s[5-1][1-1]< cutoff_.cutoff)  //gamma p1
    || (-s[5-1][2-1]< cutoff_.cutoff)  //gamma p2
    || (-s[4-1][1-1]< cutoff_.cutoff)  //e+    p1
    || (-s[4-1][2-1]< cutoff_.cutoff)  //e-    p2
    || (-s[3-1][1-1]< cutoff_.cutoff)  //nu    p1
    || (-s[3-1][2-1]< cutoff_.cutoff)  //nu    p2
    || (+s[5-1][4-1]< cutoff_.cutoff)  //gamma e+
    || (+s[5-1][3-1]< cutoff_.cutoff)  //gamma nu
    || (+s[4-1][3-1]< cutoff_.cutoff)  //e+    nu
    )
    )
    return true;

  else if (
    npart == 4 &&
    (
    (-s[5-1][1-1]< cutoff_.cutoff)  //e-    p1
    || (-s[5-1][2-1]< cutoff_.cutoff)  //e-    p2
    || (-s[6-1][1-1]< cutoff_.cutoff)  //nb    p1
    || (-s[6-1][2-1]< cutoff_.cutoff)  //nb    p2
    || (+s[6-1][5-1]< cutoff_.cutoff)  //e-    nb
    )

    )

    return true;

  return false;
}

//Make sure
// 1. tot Energy Sum < 2EBEAM
// 2. PartonEnergy Fraction minimum<x0,x1<1
// 3. number of final state particle is defined
//
double TUtil::SumMatrixElementPDF(
  TVar::Process process, TVar::Production production, TVar::MatrixElement matrixElement,
  event_scales_type* event_scales, MelaIO* RcdME,
  double EBEAM,
  double coupling[SIZE_HVV_FREENORM], // This last argument is unfortunately the simplest way to pass these couplings
  TVar::VerbosityLevel verbosity
  ){

  int partIncCode=TVar::kNoAssociated; // Do not use associated particles in the pT=0 frame boost
  int nRequested_AssociatedJets = 0;
  if (
    ((process==TVar::bkgZZ_SMHiggs || process==TVar::HSMHiggs || process==TVar::bkgZZ) && production==TVar::JJVBF)
    ){ // Use asociated jets in the pT=0 frame boost
    partIncCode=TVar::kUseAssociated_Jets;
    nRequested_AssociatedJets = 2;
  }
  simple_event_record mela_event;
  mela_event.AssociationCode=partIncCode;
  mela_event.nRequested_AssociatedJets=nRequested_AssociatedJets;
  GetBoostedParticleVectors(
    RcdME->melaCand,
    mela_event
    );

  double xx[2]={ 0 };
  if (!CheckPartonMomFraction(mela_event.pMothers.at(0).second, mela_event.pMothers.at(1).second, xx, EBEAM, TVar::ERROR)) return 0;
  if (xx[0]==0. || xx[1]==0. || EBEAM==0) return 0;

  TLorentzVector MomStore[mxpart];
  for (int i = 0; i < mxpart; i++) MomStore[i].SetXYZT(0, 0, 0, 0);

  int NPart=npart_.npart+2;
  double p4[4][mxpart] = { { 0 } };
  double s[mxpart][mxpart] = { { 0 } };
  double msq[nmsq][nmsq];
  double msqjk=0;
  int channeltoggle=0;

  //Convert TLorentzVector into 4xNPart Matrix
  //reverse sign of incident partons
  for(int ipar=0;ipar<2;ipar++){
    if(mela_event.pMothers.at(ipar).second.T()>0.){
      p4[0][ipar] = -mela_event.pMothers.at(ipar).second.X();
      p4[1][ipar] = -mela_event.pMothers.at(ipar).second.Y();
      p4[2][ipar] = -mela_event.pMothers.at(ipar).second.Z();
      p4[3][ipar] = -mela_event.pMothers.at(ipar).second.T();
      MomStore[ipar] = mela_event.pMothers.at(ipar).second;
    }
    else{
      p4[0][ipar] = mela_event.pMothers.at(ipar).second.X();
      p4[1][ipar] = mela_event.pMothers.at(ipar).second.Y();
      p4[2][ipar] = mela_event.pMothers.at(ipar).second.Z();
      p4[3][ipar] = mela_event.pMothers.at(ipar).second.T();
      MomStore[ipar] = -mela_event.pMothers.at(ipar).second;
    }
  }

  //initialize decayed particles
  for (int ipar=2; ipar<min(NPart, mela_event.pDaughters.size()+mela_event.pAssociated.size()+2); ipar++){
    TLorentzVector* momTmp;
    if (ipar<mela_event.pDaughters.size()+2) momTmp=&(mela_event.pDaughters.at(ipar-2).second);
    else momTmp=&(mela_event.pAssociated.at(ipar-2).second);
    p4[0][ipar] = momTmp->X();
    p4[1][ipar] = momTmp->Y();
    p4[2][ipar] = momTmp->Z();
    p4[3][ipar] = momTmp->T();
    MomStore[ipar]=*momTmp;
  }

  //for (int i=0; i<NPart; i++) cout << "p["<<i<<"] (Px, Py, Pz, E):\t" << p4[0][i] << '\t' << p4[1][i] << '\t' << p4[2][i] << '\t' << p4[3][i] << endl;

  double defaultRenScale = scale_.scale;
  double defaultFacScale = facscale_.facscale;
//  cout << "Default scales: " << defaultRenScale << '\t' << defaultFacScale << endl;
  int defaultNloop = nlooprun_.nlooprun;
  int defaultNflav = nflav_.nflav;
  string defaultPdflabel = pdlabel_.pdlabel;
  double renQ = InterpretScaleScheme(production, matrixElement, event_scales->renomalizationScheme, MomStore);
//  cout << "renQ: " << renQ << " x " << event_scales->ren_scale_factor << endl;
  double facQ = InterpretScaleScheme(production, matrixElement, event_scales->factorizationScheme, MomStore);
//  cout << "facQ: " << facQ << " x " << event_scales->fac_scale_factor << endl;
  SetAlphaS(renQ, facQ, event_scales->ren_scale_factor, event_scales->fac_scale_factor, 1, 5, "cteq6_l"); // Set AlphaS(|Q|/2, mynloop, mynflav, mypartonPDF) for MCFM ME-related calculations

  //calculate invariant masses between partons/final state particles
  for (int jdx=0; jdx< NPart; jdx++){
    s[jdx][jdx]=0;
    for (int kdx=jdx+1; kdx<NPart; kdx++){
      s[jdx][kdx]=2*(p4[3][jdx]*p4[3][kdx]-p4[2][jdx]*p4[2][kdx]-p4[1][jdx]*p4[1][kdx]-p4[0][jdx]*p4[0][kdx]);
      s[kdx][jdx]=s[jdx][kdx];
    }
  }

  bool passMassCuts=true;
  if (passMassCuts){
    if ((production == TVar::ZZINDEPENDENT || production == TVar::ZZQQB) && process == TVar::bkgZZ) qqb_zz_(p4[0], msq[0]);
    if (production == TVar::ZZQQB_STU && process == TVar::bkgZZ){
      channeltoggle=0;
      qqb_zz_stu_(p4[0], msq[0], &channeltoggle);
    }
    if (production == TVar::ZZQQB_S && process == TVar::bkgZZ){
      channeltoggle=1;
      qqb_zz_stu_(p4[0], msq[0], &channeltoggle);
    }
    if (production == TVar::ZZQQB_TU && process == TVar::bkgZZ){
      channeltoggle=2;
      qqb_zz_stu_(p4[0], msq[0], &channeltoggle);
    }
    //if( process==TVar::HZZ_4l)     qqb_hzz_(p4[0],msq[0]);
    // the subroutine for the calculations including the interfenrence
    // ME =  sig + inter (sign, bkg)
    // 1161 '  f(p1)+f(p2) --> H(--> Z^0(mu^-(p3)+mu^+(p4)) + Z^0(e^-(p5)+e^+(p6)) [including gg->ZZ intf.]' 'L'
    if (process==TVar::bkgZZ_SMHiggs && matrixElement==TVar::JHUGen) gg_zz_int_freenorm_(p4[0], coupling, msq[0]); // |ggZZ + ggHZZ|**2 MCFM 6.6 version
    if (process==TVar::bkgZZ_SMHiggs && matrixElement==TVar::MCFM) gg_zz_all_(p4[0], msq[0]); // |ggZZ + ggHZZ|**2
    if (process==TVar::HSMHiggs && production == TVar::ZZGG) gg_hzz_tb_(p4[0], msq[0]); // |ggHZZ|**2
    if (process==TVar::bkgZZ && production==TVar::ZZGG) gg_zz_(p4[0], &msq[5][5]); // |ggZZ|**2
    if ((process==TVar::bkgZZ_SMHiggs || process==TVar::HSMHiggs || process==TVar::bkgZZ) && production==TVar::JJVBF) qq_zzqq_(p4[0], msq[0]); // VBF MCFM SBI, S or B
/*
    // Below code sums over all production parton flavors according to PDF
    // This is disabled as we are not using the intial production information
    // the below code is fine for the particle produced by single flavor of incoming partons
    for(int ii=0;ii<nmsq;ii++){
      for(int jj=0;jj<nmsq;jj++){

        //2-D matrix is transposed in fortran
        //msq[ parton2 ] [ parton1 ]
        //flavor_msq[jj][ii] = fx1[ii]*fx2[jj]*msq[jj][ii];

        flavor_msq[jj][ii] = msq[jj][ii];
        //cout<<jj<<ii<<"="<<msq[jj][ii]<<"  ";
        msqjk+=flavor_msq[jj][ii];
      }//ii
    //cout<<"\n";
    }//jj
    // by default assume only gg productions
    // FOTRAN convention -5    -4   -3  -2    -1  0 1 2 3 4 5
    //     parton flavor bbar cbar sbar ubar dbar g d u s c b
    // C++ convention     0     1   2    3    4   5 6 7 8 9 10
    //
    */
    double msqjk_sum = SumMEPDF(MomStore[0], MomStore[1], msq, RcdME, EBEAM, verbosity);
    if (process==TVar::bkgZZ && (production == TVar::ZZQQB || production == TVar::ZZQQB_STU || production == TVar::ZZQQB_S || production == TVar::ZZQQB_TU || production ==TVar::ZZINDEPENDENT)) msqjk = msq[3][7] + msq[7][3]; // all of the unweighted MEs are the same
    else if ((process==TVar::bkgZZ_SMHiggs || process==TVar::HSMHiggs || process==TVar::bkgZZ) && production==TVar::JJVBF) msqjk = msqjk_sum; // MCFM VVH sum
    else msqjk = msq[5][5]; // gg-only

    //flux=fbGeV2/(8.*xx[0]*xx[1]*EBEAM*EBEAM);
  }

  if (msqjk != msqjk){
    cout << "SumMatrixPDF: "<< TVar::ProcessName(process) << " msqjk="  << msqjk << endl;
    msqjk=0;
  }

//  cout << "Before reset: " << scale_.scale << '\t' << facscale_.facscale << endl;
  SetAlphaS(defaultRenScale, defaultFacScale, 1., 1., defaultNloop, defaultNflav, defaultPdflabel); // Protection for other probabilities
//  cout << "Default scale reset: " << scale_.scale << '\t' << facscale_.facscale << endl;
  return msqjk;
}


double TUtil::JHUGenMatEl(
  TVar::Process process, TVar::Production production, TVar::MatrixElement matrixElement,
  event_scales_type* event_scales, MelaIO* RcdME,
  double EBEAM,
  TVar::VerbosityLevel verbosity
  ){
  const double GeV=1./100.; // JHUGen mom. scale factor
  double MatElSq=0; // Return value

  if (matrixElement!=TVar::JHUGen){ cerr << "TUtil::JHUGenMatEl: Non-JHUGen MEs are not supported" << endl; return MatElSq; }
  bool isSpinZero = (
    process == TVar::HSMHiggs
    || process == TVar::H0minus
    || process == TVar::H0hplus
    || process == TVar::H0_g1prime2
    || process== TVar::H0_Zgs
    || process ==TVar::H0_gsgs
    || process ==TVar::H0_Zgs_PS
    || process ==TVar::H0_gsgs_PS
    || process ==TVar::H0_Zgsg1prime2
    || process == TVar::SelfDefine_spin0
    );
  bool isSpinOne = (
    process == TVar::H1minus
    || process == TVar::H1plus
    || process == TVar::SelfDefine_spin1
    );
  bool isSpinTwo = (
    process == TVar::H2_g1g5
    || process == TVar::H2_g1
    || process == TVar::H2_g8
    || process == TVar::H2_g4
    || process == TVar::H2_g5
    || process == TVar::H2_g2
    || process == TVar::H2_g3
    || process == TVar::H2_g6
    || process == TVar::H2_g7
    || process == TVar::H2_g9
    || process == TVar::H2_g10
    || process == TVar::SelfDefine_spin2
    );
  if (!(isSpinZero || isSpinOne || isSpinTwo)){ cerr << "TUtil::JHUGenMatEl: Process " << process << " not supported." << endl; return MatElSq; }

  double msq[nmsq][nmsq]={ 0 }; // ME**2[parton2][parton1] for each incoming parton 1 and 2, used in RcdME
  int MYIDUP_tmp[4]={ 0 }; // Initial assignment array, unconverted. 0==Unassigned
  vector<int> idarray[4]; // All possible ids for each daughter based on the value of MYIDUP_tmp[0:3] and the desired V ids taken from mela_event.intermediateVid.at(0:1)
  int MYIDUP[4]={ 0 }; // "Working" assignment, converted
  int idfirst[2]={ 0 }; // Used to set DecayMode1, = MYIDUP[0:1]
  int idsecond[2]={ 0 }; // Used to set DecayMode2, = MYIDUP[2:3]
  double p4[6][4]={ { 0 } }; // Mom (*GeV) to pass into JHUGen
  TLorentzVector MomStore[mxpart]; // Mom (in natural units) to compute alphaS (FIXME: SETALPHAS)
  for (int i = 0; i < mxpart; i++) MomStore[i].SetXYZT(0, 0, 0, 0);

  // Notice that partIncCode is specific for this subroutine
  int partIncCode=TVar::kNoAssociated; // Do not use associated particles in the pT=0 frame boost
  simple_event_record mela_event;
  mela_event.AssociationCode=partIncCode;
  GetBoostedParticleVectors(
    RcdME->melaCand,
    mela_event
    );
  if (mela_event.pDaughters.size()<2 || mela_event.intermediateVid.size()!=2){
    cerr << "TUtil::JHUGenMatEl: Number of daughters " << mela_event.pDaughters.size() << " or number of intermediate Vs " << mela_event.intermediateVid.size() << " not supported!" << endl;
    return MatElSq;
  }

  // p(i,0:3) = (E(i),px(i),py(i),pz(i))
  // i=0,1: g1,g2 or q1, qb2 (outgoing convention)
  // i=2,3: correspond to MY_IDUP(0),MY_IDUP(1)
  // i=4,5: correspond to MY_IDUP(2),MY_IDUP(3)
  for (int ipar=0; ipar<2; ipar++){
    if (mela_event.pMothers.at(ipar).second.T()>0.){
      p4[ipar][0] = -mela_event.pMothers.at(ipar).second.T()*GeV;
      p4[ipar][1] = -mela_event.pMothers.at(ipar).second.X()*GeV;
      p4[ipar][2] = -mela_event.pMothers.at(ipar).second.Y()*GeV;
      p4[ipar][3] = -mela_event.pMothers.at(ipar).second.Z()*GeV;
      MomStore[ipar] = mela_event.pMothers.at(ipar).second;
    }
    else{
      p4[ipar][0] = mela_event.pMothers.at(ipar).second.T()*GeV;
      p4[ipar][1] = mela_event.pMothers.at(ipar).second.X()*GeV;
      p4[ipar][2] = mela_event.pMothers.at(ipar).second.Y()*GeV;
      p4[ipar][3] = mela_event.pMothers.at(ipar).second.Z()*GeV;
      MomStore[ipar] = -mela_event.pMothers.at(ipar).second;
    }
    // From Markus:
    // Note that the momentum no.2, p(1:4, 2), is a dummy which is not used in case production==TVar::ZZINDEPENDENT.
    if (ipar==1 && production==TVar::ZZINDEPENDENT){ for (int ix=0; ix<4; ix++){ p4[0][ix] += p4[ipar][ix]; p4[ipar][ix]=0.; } }
  }
  //initialize decayed particles
  if (mela_event.pDaughters.size()==2){
    for (int ipar=0; ipar<mela_event.pDaughters.size(); ipar++){
      TLorentzVector* momTmp = &(mela_event.pDaughters.at(ipar).second);
      int* idtmp = &(mela_event.pDaughters.at(ipar).first);

      int arrindex = ipar;
      if (PDGHelpers::isAPhoton(*idtmp) && ipar==1) arrindex=2; // In GG
      if (!PDGHelpers::isAnUnknownJet(*idtmp)) MYIDUP_tmp[arrindex] = *idtmp;
      else MYIDUP_tmp[ipar] = 0;
      p4[arrindex+2][0] = momTmp->T()*GeV;
      p4[arrindex+2][1] = momTmp->X()*GeV;
      p4[arrindex+2][2] = momTmp->Y()*GeV;
      p4[arrindex+2][3] = momTmp->Z()*GeV;
      MomStore[arrindex+2] = *momTmp;
    }
  }
  else{
    for (int ipar=0; ipar<4; ipar++){
      if (ipar<mela_event.pDaughters.size()){
        TLorentzVector* momTmp = &(mela_event.pDaughters.at(ipar).second);
        int* idtmp = &(mela_event.pDaughters.at(ipar).first);

        if (!PDGHelpers::isAnUnknownJet(*idtmp)) MYIDUP_tmp[ipar] = *idtmp;
        else MYIDUP_tmp[ipar] = 0;
        p4[ipar+2][0] = momTmp->T()*GeV;
        p4[ipar+2][1] = momTmp->X()*GeV;
        p4[ipar+2][2] = momTmp->Y()*GeV;
        p4[ipar+2][3] = momTmp->Z()*GeV;
        MomStore[ipar+2] = *momTmp;
      }
      else MYIDUP_tmp[ipar] = -9000; // No need to set p4, which is already 0 by initialization
      // __modparameters_MOD_not_a_particle__?
    }
    cout << "MYIDUP_tmp[" << ipar << "]=" << MYIDUP_tmp[ipar] << endl;
  }

  for (int idau=0; idau<4; idau++) cout << "MYIDUP_tmp[" << idau << "]=" << MYIDUP_tmp[idau] << endl;

  // FIXME: SETALPHAS DOES NOT MODIFY JHUGENMELA
  // Set alphas
  double defaultRenScale = scale_.scale;
  double defaultFacScale = facscale_.facscale;
  //  cout << "Default scales: " << defaultRenScale << '\t' << defaultFacScale << endl;
  int defaultNloop = nlooprun_.nlooprun;
  int defaultNflav = nflav_.nflav;
  string defaultPdflabel = pdlabel_.pdlabel;
  double renQ = InterpretScaleScheme(production, matrixElement, event_scales->renomalizationScheme, MomStore);
  //  cout << "renQ: " << renQ << " x " << event_scales->ren_scale_factor << endl;
  double facQ = InterpretScaleScheme(production, matrixElement, event_scales->factorizationScheme, MomStore);
  //  cout << "facQ: " << facQ << " x " << event_scales->fac_scale_factor << endl;
  SetAlphaS(renQ, facQ, event_scales->ren_scale_factor, event_scales->fac_scale_factor, 1, 5, "cteq6_l"); // Set AlphaS(|Q|/2, mynloop, mynflav, mypartonPDF) for MCFM ME-related calculations

  for (int idau=0; idau<4; idau++){
    // If id is known, just assign it.
    if (MYIDUP_tmp[idau]!=0) idarray[idau].push_back(MYIDUP_tmp[idau]);
    else{
      if (idau<2 && PDGHelpers::isAZBoson(mela_event.intermediateVid.at(0))){ // Z->ffb
        if (MYIDUP_tmp[1-idau]!=0) idarray[idau].push_back(-MYIDUP_tmp[1-idau]); // Z->f+unknown (I don't know how this could happen, but cover this case as well)
        else if(idau==0){
          for (int iquark=1; iquark<=5; iquark++){ // iquark in PDG convention; 1/3/5 are dn-type
            idarray[idau].push_back(iquark);
            idarray[1-idau].push_back(-iquark);
          }
        }
      }
      else if (idau<2 && PDGHelpers::isAWBoson(mela_event.intermediateVid.at(0))){ // (W+)->ffb
        if (MYIDUP_tmp[1-idau]!=0){ // One quark is known
          if (PDGHelpers::isUpTypeQuark(MYIDUP_tmp[1-idau])){
            int id_dn = TMath::Sign(1, -MYIDUP_tmp[1-idau]);
            int id_st = TMath::Sign(3, -MYIDUP_tmp[1-idau]);
            int id_bt = TMath::Sign(5, -MYIDUP_tmp[1-idau]);
            idarray[idau].push_back(id_dn);
            idarray[idau].push_back(id_st);
            idarray[idau].push_back(id_bt);
          }
          else if (PDGHelpers::isDownTypeQuark(MYIDUP_tmp[1-idau])){
            int id_up = TMath::Sign(2, -MYIDUP_tmp[1-idau]);
            int id_ch = TMath::Sign(4, -MYIDUP_tmp[1-idau]);
            idarray[idau].push_back(id_up);
            idarray[idau].push_back(id_ch);
          }
        }
        else if (idau==0){ // Both quarks unknown
          for (int iquark=1; iquark<=5; iquark++){ // iquark in PDG convention; 1/3/5 are dn-type
            if (PDGHelpers::isUpTypeQuark(iquark)) idarray[idau].push_back(iquark); // Assign u to f
            else idarray[1-idau].push_back(-iquark); // Assign db to fb
          }
        }
      }
      if (idau>=2 && PDGHelpers::isAZBoson(mela_event.intermediateVid.at(1))){ // Z->ffb
        if (MYIDUP_tmp[5-idau]!=0) idarray[idau].push_back(-MYIDUP_tmp[5-idau]); // Z->f+unknown (I don't know how this could happen, but cover this case as well)
        else if (idau==2){
          for (int iquark=1; iquark<=5; iquark++){ // iquark in PDG convention; 1/3/5 are dn-type
            idarray[idau].push_back(iquark);
            idarray[5-idau].push_back(-iquark);
          }
        }
      }
      else if (idau>=2 && PDGHelpers::isAWBoson(mela_event.intermediateVid.at(1))){ // (W-)->ffb
        if (MYIDUP_tmp[5-idau]!=0){ // One quark is known
          if (PDGHelpers::isUpTypeQuark(MYIDUP_tmp[5-idau])){
            int id_dn = TMath::Sign(1, -MYIDUP_tmp[5-idau]);
            int id_st = TMath::Sign(3, -MYIDUP_tmp[5-idau]);
            int id_bt = TMath::Sign(5, -MYIDUP_tmp[5-idau]);
            idarray[idau].push_back(id_dn);
            idarray[idau].push_back(id_st);
            idarray[idau].push_back(id_bt);
          }
          if (PDGHelpers::isDownTypeQuark(MYIDUP_tmp[5-idau])){
            int id_up = TMath::Sign(2, -MYIDUP_tmp[5-idau]);
            int id_ch = TMath::Sign(4, -MYIDUP_tmp[5-idau]);
            idarray[idau].push_back(id_up);
            idarray[idau].push_back(id_ch);
          }
        }
        else if (idau==2){ // Both quarks unknown
          for (int iquark=-5; iquark<=-1; iquark++){ // iquark in PDG convention; 1/3/5 are dn-type
            if (PDGHelpers::isDownTypeQuark(iquark)) idarray[idau].push_back(-iquark); // Assign d to f
            else idarray[5-idau].push_back(iquark); // Assign ub to fb
          }
        }
      }
    }
  }

  int nNonZero=0;
  for (int d1=0; d1<idarray[0].size(); d1++){
    for (int d2=0; d2<idarray[1].size(); d2++){
      for (int d3=0; d1<idarray[2].size(); d3++){
        for (int d4=0; d1<idarray[3].size(); d4++){
          if (
            (PDGHelpers::isAZBoson(mela_event.intermediateVid.at(0)) && idarray[0].at(d1)!=-idarray[1].at(d2))
            ||
            (PDGHelpers::isAZBoson(mela_event.intermediateVid.at(1)) && idarray[2].at(d3)!=-idarray[2].at(d4))
            ) continue;
          // Convert the particles
          MYIDUP[0] = convertLHEreverse(&(idarray[0].at(d1)));
          MYIDUP[1] = convertLHEreverse(&(idarray[1].at(d2)));
          MYIDUP[2] = convertLHEreverse(&(idarray[2].at(d3)));
          MYIDUP[3] = convertLHEreverse(&(idarray[3].at(d4)));
          // TEST: Print particles
          cout << "TUtil::JHUGenMatEl: d1-d4=\t" << d1 << '\t' << d2 << '\t' << d3 << '\t' << d4 << endl;
          for (int idau=0; idau<4; idau++) cout << "MYIDUP[" << idau << "]=" << MYIDUP[idau] << endl;
          // Determine M_V and Ga_V in JHUGen, needed for g1 vs everything else.
          for (int ip=0; ip<2; ip++){ idfirst[ip]=MYIDUP[ip]; idsecond[ip]=MYIDUP[ip+2]; }
          __modjhugenmela_MOD_setdecaymodes(idfirst, idsecond); // Set M_V and Ga_V in JHUGen

          double MatElTmp=0.;
          int idIntermediate[2]={ mela_event.intermediateVid.at(0), mela_event.intermediateVid.at(1) };
          if (
            (PDGHelpers::isAWBoson(idIntermediate[0]) && PDGHelpers::isAWBoson(idIntermediate[1]))
            ||
            (
            (PDGHelpers::isAZBoson(idIntermediate[0]) || PDGHelpers::isAPhoton(idIntermediate[0]))
            &&
            (PDGHelpers::isAZBoson(idIntermediate[1]) || PDGHelpers::isAPhoton(idIntermediate[1]))
            )
            ){
            if (production == TVar::ZZGG){
              if (isSpinZero){
                //__modkinematics_MOD_evalalphas();
                __modhiggs_MOD_evalamp_gg_h_vv(p4, MYIDUP, &MatElTmp);
              }
              else if (isSpinTwo) __modgraviton_MOD_evalamp_gg_g_vv(p4, MYIDUP, &MatElTmp);
            }
            else if (production == TVar::ZZQQB){
              if (isSpinOne) __modzprime_MOD_evalamp_qqb_zprime_vv(p4, MYIDUP, &MatElTmp);
              else if (isSpinTwo) __modgraviton_MOD_evalamp_qqb_g_vv(p4, MYIDUP, &MatElTmp);
            }
            else if (production == TVar::ZZINDEPENDENT){
              if (isSpinZero) __modhiggs_MOD_evalamp_h_vv(p4, MYIDUP, &MatElTmp);
              else if (isSpinOne) __modzprime_MOD_evalamp_zprime_vv(p4, MYIDUP, &MatElTmp);
              else if (isSpinTwo) __modgraviton_MOD_evalamp_g_vv(p4, MYIDUP, &MatElTmp);
            }
            // Add CKM elements since they are not included
            if (PDGHelpers::isAWBoson(idIntermediate[0])) MatElTmp *= pow(__modparameters_MOD_ckm(&(idarray[0].at(d1)), &(idarray[1].at(d2)))/__modparameters_MOD_scalefactor(&(idarray[0].at(d1)), &(idarray[1].at(d2))), 2);
            if (PDGHelpers::isAWBoson(idIntermediate[1])) MatElTmp *= pow(__modparameters_MOD_ckm(&(idarray[2].at(d3)), &(idarray[3].at(d4)))/__modparameters_MOD_scalefactor(&(idarray[2].at(d3)), &(idarray[3].at(d4))), 2);
          }

          cout << "TUtil::JHUGenMatEl: MatElTmp = " << MatElTmp << endl;
          MatElSq += MatElTmp;
          if (MatElTmp>0.) nNonZero++;
        }
      }
    }
  }
  if (nNonZero>0) MatElSq /= ((double)nNonZero);
  cout << "TUtil::JHUGenMatEl: Number of matrix element instances computed: " << nNonZero << endl;
  // This constant is needed to account for the different units used in
  // JHUGen compared to the MCFM
  double constant = 1.45e-8;
  if (isSpinZero && production!=TVar::ZZINDEPENDENT) constant = 4.46162946e-4;
  MatElSq *= constant;
  // == 1.45e-8/pow(0.13229060/(3.*3.141592653589793238462643383279502884197*2.4621845810181631), 2), constant was 1.45e-8 before with alpha_s=0.13229060, vev=2.4621845810181631/GeV

  // Set RcdME information for ME and parton distributions, taking into account the mothers if id!=0 (i.e. if not unknown).
  if (mela_event.pMothers.at(0).first!=0 && mela_event.pMothers.at(1).first!=0){
    int ix=0, iy=0;
    if (abs(mela_event.pMothers.at(0).first)<6) ix=mela_event.pMothers.at(0).first;
    if (abs(mela_event.pMothers.at(1).first)<6) iy=mela_event.pMothers.at(1).first;
    msq[iy][ix]=MatElSq; // Note that SumMEPdf receives a transposed msq
  }
  else{
    if (production == TVar::ZZGG || production == TVar::ZZINDEPENDENT) msq[5][5]=MatElSq;
    else if (production == TVar::ZZQQB){
      for (int ix=0; ix<5; ix++){ msq[ix][10-ix]=MatElSq; msq[10-ix][ix]=MatElSq; }
    }
  }
  if (production!=TVar::ZZINDEPENDENT) SumMEPDF(MomStore[0], MomStore[1], msq, RcdME, EBEAM, verbosity);
  else{ // If production is ZZINDEPENDENT, only set gg index with fx1,2[g,g]=1.
    double fx_dummy[nmsq]={ 0 }; fx_dummy[5]=1.;
    RcdME->setPartonWeights(fx_dummy, fx_dummy);
    RcdME->setMEArray(msq, true);
    RcdME->computeWeightedMEArray();
  }

/*
  for (int i=0; i<6; i++) cout << "p4["<<i<<"] (Px, Py, Pz, E):\t" << p4[i][1]*100. << '\t' << p4[i][2]*100. << '\t' << p4[i][3]*100. << '\t' << p4[i][0]*100. << endl;
  cout << "MatElSq=" << endl;
*/

  // Reset alphas
  SetAlphaS(defaultRenScale, defaultFacScale, 1., 1., defaultNloop, defaultNflav, defaultPdflabel); // Protection for other probabilities
  //cout << "TUtil::MatElSq = " << MatElSq << endl;
  return MatElSq;
}

double TUtil::HJJMatEl(
  TVar::Process process, TVar::Production production, TVar::MatrixElement matrixElement,
  event_scales_type* event_scales, MelaIO* RcdME,
  double EBEAM,
  TVar::VerbosityLevel verbosity
  ){
  const double GeV=1./100.; // JHUGen mom. scale factor
  double sum_msqjk = 0;
  // by default assume only gg productions
  // FOTRAN convention -5    -4   -3  -2    -1  0 1 2 3 4 5
  //     parton flavor bbar cbar sbar ubar dbar g d u s c b
  // C++ convention     0     1   2    3    4   5 6 7 8 9 10
  //2-D matrix is reversed in fortran
  // msq[ parton2 ] [ parton1 ]
  //      flavor_msq[jj][ii] = fx1[ii]*fx2[jj]*msq[jj][ii];
  double MatElsq[nmsq][nmsq]={ { 0 } };
  double MatElsq_tmp[nmsq][nmsq]={ { 0 } };

  if (matrixElement!=TVar::JHUGen){ cerr << "TUtil::HJJMatEl: Non-JHUGen MEs are not supported" << endl; return sum_msqjk; }
  if (!(production==TVar::JJGG || production==TVar::JJVBF || production==TVar::JH)){ cerr << "TUtil::HJJMatEl: Production is not supported!" << endl; return sum_msqjk; }

  // Notice that partIncCode is specific for this subroutine
  int nRequested_AssociatedJets=2;
  if (production == TVar::JH) nRequested_AssociatedJets=1;
  int partIncCode=TVar::kUseAssociated_Jets; // Only use associated partons in the pT=0 frame boost
  simple_event_record mela_event;
  mela_event.AssociationCode=partIncCode;
  mela_event.nRequested_AssociatedJets=nRequested_AssociatedJets;
  GetBoostedParticleVectors(
    RcdME->melaCand,
    mela_event
    );
  if (mela_event.pAssociated.size()==0){ cerr << "TUtil::HJJMatEl: Number of associated particles is 0!" << endl; return sum_msqjk; }

  int MYIDUP_tmp[4]={ 0 }; // "Incoming" partons 1, 2, "outgoing" partons 3, 4
  double p4[5][4]={ { 0 } };
  double pOneJet[4][4] ={ { 0 } }; // For HJ
  TLorentzVector MomStore[mxpart];
  for (int i = 0; i < mxpart; i++) MomStore[i].SetXYZT(0, 0, 0, 0);

  // p4(0:4,i) = (E(i),px(i),py(i),pz(i))
  // i=0,1: g1,g2 or q1, qb2 (outgoing convention)
  // i=2,3: J1, J2 in outgoing convention when possible
  // i=4: H
  for (int ipar=0; ipar<2; ipar++){
    TLorentzVector* momTmp = &(mela_event.pMothers.at(ipar).second);
    int* idtmp = &(mela_event.pMothers.at(ipar).first);
    if (!PDGHelpers::isAnUnknownJet(*idtmp)) MYIDUP_tmp[ipar] = *idtmp;
    else MYIDUP_tmp[ipar] = 0;
    if (momTmp->T()>0.){
      p4[ipar][0] = momTmp->T()*GeV;
      p4[ipar][1] = momTmp->X()*GeV;
      p4[ipar][2] = momTmp->Y()*GeV;
      p4[ipar][3] = momTmp->Z()*GeV;
      MomStore[ipar] = (*momTmp);
    }
    else{
      p4[ipar][0] = -momTmp->T()*GeV;
      p4[ipar][1] = -momTmp->X()*GeV;
      p4[ipar][2] = -momTmp->Y()*GeV;
      p4[ipar][3] = -momTmp->Z()*GeV;
      MomStore[ipar] = -(*momTmp);
      MYIDUP_tmp[ipar] = -MYIDUP_tmp[ipar];
    }
  }
  for (int ipar=0; ipar<2; ipar++){
    if (ipar<mela_event.pAssociated.size()){
      TLorentzVector* momTmp = &(mela_event.pAssociated.at(ipar).second);
      int* idtmp = &(mela_event.pAssociated.at(ipar).first);
      if (!PDGHelpers::isAnUnknownJet(*idtmp)) MYIDUP_tmp[ipar+2] = *idtmp;
      else MYIDUP_tmp[ipar+2] = 0;
      p4[ipar+2][0] = momTmp->T()*GeV;
      p4[ipar+2][1] = momTmp->X()*GeV;
      p4[ipar+2][2] = momTmp->Y()*GeV;
      p4[ipar+2][3] = momTmp->Z()*GeV;
      MomStore[ipar+6] = (*momTmp); // i==(2, 3, 4) is (J1, J2, H), recorded as MomStore (I1, I2, 0, 0, 0, H, J1, J2)
    }
    else MYIDUP_tmp[ipar+2] = -9000; // No need to set p4, which is already 0 by initialization
  }
  for (int ipar=0; ipar<mela_event.pDaughters.size(); ipar++){
    TLorentzVector* momTmp = &(mela_event.pDaughters.at(ipar).second);
    p4[4][0] += momTmp->T()*GeV;
    p4[4][1] += momTmp->X()*GeV;
    p4[4][2] += momTmp->Y()*GeV;
    p4[4][3] += momTmp->Z()*GeV;
    MomStore[5] = MomStore[5] + (*momTmp); // i==(2, 3, 4) is (J1, J2, H), recorded as MomStore (I1, I2, 0, 0, 0, H, J1, J2)
  }
  // Momenta for HJ
  for (int i = 0; i < 4; i++) {
    if (i<3){ for (int j = 0; j < 4; j++) pOneJet[i][j] = p4[i][j]; } // p1 p2 J1
    else{ for (int j = 0; j < 4; j++) pOneJet[i][j] = p4[i+1][j]; } // H
  }
  if (verbosity >= TVar::DEBUG){
    for (int i=0; i<5; i++) cout << "p["<<i<<"] (Px, Py, Pz, E):\t" << p4[i][1]/GeV << '\t' << p4[i][2]/GeV << '\t' << p4[i][3]/GeV << '\t' << p4[i][0]/GeV << endl;
  }

  double defaultRenScale = scale_.scale;
  double defaultFacScale = facscale_.facscale;
  //cout << "Default scales: " << defaultRenScale << '\t' << defaultFacScale << endl;
  int defaultNloop = nlooprun_.nlooprun;
  int defaultNflav = nflav_.nflav;
  string defaultPdflabel = pdlabel_.pdlabel;
  double renQ = InterpretScaleScheme(production, matrixElement, event_scales->renomalizationScheme, MomStore);
  //cout << "renQ: " << renQ << " x " << event_scales->ren_scale_factor << endl;
  double facQ = InterpretScaleScheme(production, matrixElement, event_scales->factorizationScheme, MomStore);
  //cout << "facQ: " << facQ << " x " << event_scales->fac_scale_factor << endl;
  SetAlphaS(renQ, facQ, event_scales->ren_scale_factor, event_scales->fac_scale_factor, 1, 5, "cteq6_l"); // Set AlphaS(|Q|/2, mynloop, mynflav, mypartonPDF) for MCFM ME-related calculations

  // NOTE ON CHANNEL HASHES:
  // THEY ONLY RETURN ISEL>=JSEL CASES. ISEL<JSEL NEEDS TO BE DONE MANUALLY.
  if (production == TVar::JH){ // Computation is already for all possible qqb/qg/qbg/gg, and incoming q, qb and g flavor have 1-1 correspondence to the outgoing jet flavor.
    __modhiggsj_MOD_evalamp_hj(pOneJet, MatElsq_tmp);
    for (int isel=-5; isel<=5; isel++){
      if (MYIDUP_tmp[0]!=0 && !((PDGHelpers::isAGluon(MYIDUP_tmp[0]) && isel==0) || MYIDUP_tmp[0]==isel)) continue;
      for (int jsel=-5; jsel<=5; jsel++){
        if (MYIDUP_tmp[1]!=0 && !((PDGHelpers::isAGluon(MYIDUP_tmp[1]) && jsel==0) || MYIDUP_tmp[1]==jsel)) continue;
        int rsel;
        if (isel!=0 && jsel!=0) rsel=0; // Covers qqb->Hg
        else if (isel==0) rsel=jsel; // Covers gg->Hg, gq->Hq, gqb->Hqb
        else rsel=isel; // Covers qg->Hq, qbg->Hqb
        if (MYIDUP_tmp[2]!=0 && !((PDGHelpers::isAGluon(MYIDUP_tmp[2]) && rsel==0) || MYIDUP_tmp[2]==rsel)) continue;
        MatElsq[jsel+5][isel+5] = MatElsq_tmp[jsel+5][isel+5]; // Assign only those that match gen. info, if present at all.
      }
    }
  }
  else if (production==TVar::JJGG){
    int ijsel[3][121];
    int nijchannels=77;
    __modhiggsjj_MOD_get_hjjchannelhash_nosplit(ijsel, &nijchannels);
    for (int ic=0; ic<nijchannels; ic++){
      // Emulate EvalWeighted_HJJ_test
      int isel = ijsel[0][ic];
      int jsel = ijsel[1][ic];
      int code = ijsel[2][ic];
      // Default assignments
      int rsel=isel;
      int ssel=jsel;
      if (
        (MYIDUP_tmp[0]==0 || ((PDGHelpers::isAGluon(MYIDUP_tmp[0]) && isel==0) || MYIDUP_tmp[0]==isel))
        &&
        (MYIDUP_tmp[1]==0 || ((PDGHelpers::isAGluon(MYIDUP_tmp[1]) && jsel==0) || MYIDUP_tmp[1]==jsel))
        ){ // Do it this way to be able to swap isel and jsel later

        if (isel==0 && jsel==0){ // gg->?
          if (code==2){ // gg->qqb
            // Only compute u-ub. The amplitude is multiplied by nf=5
            rsel=1;
            ssel=-1;
            double avgfac=1.; if (MYIDUP_tmp[2]==0 && MYIDUP_tmp[3]==0) avgfac=0.5;
            if (
              (MYIDUP_tmp[2]==0 || (PDGHelpers::isAQuark(MYIDUP_tmp[2]) && MYIDUP_tmp[2]>0))
              &&
              (MYIDUP_tmp[3]==0 || (PDGHelpers::isAQuark(MYIDUP_tmp[3]) && MYIDUP_tmp[3]<0))
              ){
              __modhiggsjj_MOD_evalamp_sbfh_unsymm_sa_select_exact(p4, &isel, &jsel, &rsel, &ssel, MatElsq_tmp);
              MatElsq[jsel+5][isel+5] += MatElsq_tmp[jsel+5][isel+5]*avgfac; // Assign only those that match gen. info, if present at all.
            }
            if (
              (MYIDUP_tmp[2]==0 || (PDGHelpers::isAQuark(MYIDUP_tmp[2]) && MYIDUP_tmp[2]<0))
              &&
              (MYIDUP_tmp[3]==0 || (PDGHelpers::isAQuark(MYIDUP_tmp[3]) && MYIDUP_tmp[3]>0))
              ){
              __modhiggsjj_MOD_evalamp_sbfh_unsymm_sa_select_exact(p4, &isel, &jsel, &ssel, &rsel, MatElsq_tmp);
              MatElsq[jsel+5][isel+5] += MatElsq_tmp[jsel+5][isel+5]*avgfac; // Assign only those that match gen. info, if present at all.
            }
          }
          else{ // gg->gg
            // rsel=ssel=g already
            if (
              (MYIDUP_tmp[2]==0 || (PDGHelpers::isAGluon(MYIDUP_tmp[2]) && rsel==0))
              &&
              (MYIDUP_tmp[3]==0 || (PDGHelpers::isAGluon(MYIDUP_tmp[3]) && ssel==0))
              ){
              __modhiggsjj_MOD_evalamp_sbfh_unsymm_sa_select_exact(p4, &isel, &jsel, &rsel, &ssel, MatElsq_tmp);
              MatElsq[jsel+5][isel+5] += MatElsq_tmp[jsel+5][isel+5]; // Assign only those that match gen. info, if present at all.
            }
          }
        }
        else if (isel==0 || jsel==0){ // qg/qbg/gq/gqb->qg/qbg/gq/gqb
          double avgfac=1.; if (MYIDUP_tmp[2]==0 && MYIDUP_tmp[3]==0) avgfac=0.5;
          if (
            (MYIDUP_tmp[2]==0 || ((PDGHelpers::isAGluon(MYIDUP_tmp[2]) && rsel==0) || MYIDUP_tmp[2]==rsel))
            &&
            (MYIDUP_tmp[3]==0 || ((PDGHelpers::isAGluon(MYIDUP_tmp[3]) && ssel==0) || MYIDUP_tmp[3]==ssel))
            ){
            __modhiggsjj_MOD_evalamp_sbfh_unsymm_sa_select_exact(p4, &isel, &jsel, &rsel, &ssel, MatElsq_tmp);
            MatElsq[jsel+5][isel+5] += MatElsq_tmp[jsel+5][isel+5]*avgfac; // Assign only those that match gen. info, if present at all.
          }
          if (
            (MYIDUP_tmp[2]==0 || ((PDGHelpers::isAGluon(MYIDUP_tmp[2]) && ssel==0) || MYIDUP_tmp[2]==ssel))
            &&
            (MYIDUP_tmp[3]==0 || ((PDGHelpers::isAGluon(MYIDUP_tmp[3]) && rsel==0) || MYIDUP_tmp[3]==rsel))
            ){
            __modhiggsjj_MOD_evalamp_sbfh_unsymm_sa_select_exact(p4, &isel, &jsel, &ssel, &rsel, MatElsq_tmp);
            MatElsq[jsel+5][isel+5] += MatElsq_tmp[jsel+5][isel+5]*avgfac; // Assign only those that match gen. info, if present at all.
          }
        }
        else if ((isel>0 && jsel<0) || (isel<0 && jsel>0)){ // qQb/qbQ->?
          if (code==1){ // qqb/qbq->gg
            rsel=0; ssel=0;
            if (
              (MYIDUP_tmp[2]==0 || PDGHelpers::isAGluon(MYIDUP_tmp[2]))
              &&
              (MYIDUP_tmp[3]==0 || PDGHelpers::isAGluon(MYIDUP_tmp[3]))
              ){
              __modhiggsjj_MOD_evalamp_sbfh_unsymm_sa_select_exact(p4, &isel, &jsel, &rsel, &ssel, MatElsq_tmp);
              MatElsq[jsel+5][isel+5] += MatElsq_tmp[jsel+5][isel+5]; // Assign only those that match gen. info, if present at all.
            }
          }
          else if (code==2){ // qQb/qbQ->qQb/qbQ
            double avgfac=1.; if (MYIDUP_tmp[2]==0 && MYIDUP_tmp[3]==0) avgfac=0.5;
            if (
              (MYIDUP_tmp[2]==0 || MYIDUP_tmp[2]==rsel)
              &&
              (MYIDUP_tmp[3]==0 || MYIDUP_tmp[3]==ssel)
              ){
              __modhiggsjj_MOD_evalamp_sbfh_unsymm_sa_select_exact(p4, &isel, &jsel, &rsel, &ssel, MatElsq_tmp);
              MatElsq[jsel+5][isel+5] += MatElsq_tmp[jsel+5][isel+5]*avgfac; // Assign only those that match gen. info, if present at all.
            }
            if (
              (MYIDUP_tmp[2]==0 || MYIDUP_tmp[2]==ssel)
              &&
              (MYIDUP_tmp[3]==0 || MYIDUP_tmp[3]==rsel)
              ){
              __modhiggsjj_MOD_evalamp_sbfh_unsymm_sa_select_exact(p4, &isel, &jsel, &ssel, &rsel, MatElsq_tmp);
              MatElsq[jsel+5][isel+5] += MatElsq_tmp[jsel+5][isel+5]*avgfac; // Assign only those that match gen. info, if present at all.
            }
          }
          else{ // qqb->QQb
            if (abs(isel)!=1){ rsel=1; ssel=-1; } // Make sure rsel, ssel are not of same flavor as isel, jsel
            else{ rsel=2; ssel=-2; }
            // The amplitude is aready multiplied by nf-1, so no need to calculate everything (nf-1) times.
            double avgfac=1.; if (MYIDUP_tmp[2]==0 && MYIDUP_tmp[3]==0) avgfac=0.5;
            if (
              (MYIDUP_tmp[2]==0 || MYIDUP_tmp[2]==rsel)
              &&
              (MYIDUP_tmp[3]==0 || MYIDUP_tmp[3]==ssel)
              ){
              __modhiggsjj_MOD_evalamp_sbfh_unsymm_sa_select_exact(p4, &isel, &jsel, &rsel, &ssel, MatElsq_tmp);
              MatElsq[jsel+5][isel+5] += MatElsq_tmp[jsel+5][isel+5]*avgfac; // Assign only those that match gen. info, if present at all.
            }
            if (
              (MYIDUP_tmp[2]==0 || MYIDUP_tmp[2]==ssel)
              &&
              (MYIDUP_tmp[3]==0 || MYIDUP_tmp[3]==rsel)
              ){
              __modhiggsjj_MOD_evalamp_sbfh_unsymm_sa_select_exact(p4, &isel, &jsel, &ssel, &rsel, MatElsq_tmp);
              MatElsq[jsel+5][isel+5] += MatElsq_tmp[jsel+5][isel+5]*avgfac; // Assign only those that match gen. info, if present at all.
            }
          }
        }
        else{
          double avgfac=1.; if (MYIDUP_tmp[2]==0 && MYIDUP_tmp[3]==0 && rsel!=ssel) avgfac=0.5;
          if (
            (MYIDUP_tmp[2]==0 || MYIDUP_tmp[2]==rsel)
            &&
            (MYIDUP_tmp[3]==0 || MYIDUP_tmp[3]==ssel)
            ){
            __modhiggsjj_MOD_evalamp_sbfh_unsymm_sa_select_exact(p4, &isel, &jsel, &rsel, &ssel, MatElsq_tmp);
            MatElsq[jsel+5][isel+5] += MatElsq_tmp[jsel+5][isel+5]*avgfac; // Assign only those that match gen. info, if present at all.
          }
          if (
            rsel!=ssel
            &&
            (MYIDUP_tmp[2]==0 || MYIDUP_tmp[2]==ssel)
            &&
            (MYIDUP_tmp[3]==0 || MYIDUP_tmp[3]==rsel)
            ){
            __modhiggsjj_MOD_evalamp_sbfh_unsymm_sa_select_exact(p4, &isel, &jsel, &ssel, &rsel, MatElsq_tmp);
            MatElsq[jsel+5][isel+5] += MatElsq_tmp[jsel+5][isel+5]*avgfac; // Assign only those that match gen. info, if present at all.
          }
        }
      } // End unswapped isel>=jsel cases
      if (isel==jsel) continue;
      isel = ijsel[1][ic];
      jsel = ijsel[0][ic];
      // Reset to default assignments
      rsel=isel;
      ssel=jsel;
      if (
        (MYIDUP_tmp[0]==0 || ((PDGHelpers::isAGluon(MYIDUP_tmp[0]) && isel==0) || MYIDUP_tmp[0]==isel))
        &&
        (MYIDUP_tmp[1]==0 || ((PDGHelpers::isAGluon(MYIDUP_tmp[1]) && jsel==0) || MYIDUP_tmp[1]==jsel))
        ){
        // isel==jsel==0 is already eliminated by isel!=jsel condition
        if (isel==0 || jsel==0){ // qg/qbg/gq/gqb->qg/qbg/gq/gqb
          double avgfac=1.; if (MYIDUP_tmp[2]==0 && MYIDUP_tmp[3]==0) avgfac=0.5;
          if (
            (MYIDUP_tmp[2]==0 || ((PDGHelpers::isAGluon(MYIDUP_tmp[2]) && rsel==0) || MYIDUP_tmp[2]==rsel))
            &&
            (MYIDUP_tmp[3]==0 || ((PDGHelpers::isAGluon(MYIDUP_tmp[3]) && ssel==0) || MYIDUP_tmp[3]==ssel))
            ){
            __modhiggsjj_MOD_evalamp_sbfh_unsymm_sa_select_exact(p4, &isel, &jsel, &rsel, &ssel, MatElsq_tmp);
            MatElsq[jsel+5][isel+5] += MatElsq_tmp[jsel+5][isel+5]*avgfac; // Assign only those that match gen. info, if present at all.
          }
          if (
            (MYIDUP_tmp[2]==0 || ((PDGHelpers::isAGluon(MYIDUP_tmp[2]) && ssel==0) || MYIDUP_tmp[2]==ssel))
            &&
            (MYIDUP_tmp[3]==0 || ((PDGHelpers::isAGluon(MYIDUP_tmp[3]) && rsel==0) || MYIDUP_tmp[3]==rsel))
            ){
            __modhiggsjj_MOD_evalamp_sbfh_unsymm_sa_select_exact(p4, &isel, &jsel, &ssel, &rsel, MatElsq_tmp);
            MatElsq[jsel+5][isel+5] += MatElsq_tmp[jsel+5][isel+5]*avgfac; // Assign only those that match gen. info, if present at all.
          }
        }
        else if ((isel>0 && jsel<0) || (isel<0 && jsel>0)){ // qQb/qbQ->?
          if (code==1){ // qqb/qbq->gg
            rsel=0; ssel=0;
            if (
              (MYIDUP_tmp[2]==0 || PDGHelpers::isAGluon(MYIDUP_tmp[2]))
              &&
              (MYIDUP_tmp[3]==0 || PDGHelpers::isAGluon(MYIDUP_tmp[3]))
              ){
              __modhiggsjj_MOD_evalamp_sbfh_unsymm_sa_select_exact(p4, &isel, &jsel, &rsel, &ssel, MatElsq_tmp);
              MatElsq[jsel+5][isel+5] += MatElsq_tmp[jsel+5][isel+5]; // Assign only those that match gen. info, if present at all.
            }
          }
          else if (code==2){ // qQb/qbQ->qQb/qbQ
            double avgfac=1.; if (MYIDUP_tmp[2]==0 && MYIDUP_tmp[3]==0) avgfac=0.5;
            if (
              (MYIDUP_tmp[2]==0 || MYIDUP_tmp[2]==rsel)
              &&
              (MYIDUP_tmp[3]==0 || MYIDUP_tmp[3]==ssel)
              ){
              __modhiggsjj_MOD_evalamp_sbfh_unsymm_sa_select_exact(p4, &isel, &jsel, &rsel, &ssel, MatElsq_tmp);
              MatElsq[jsel+5][isel+5] += MatElsq_tmp[jsel+5][isel+5]*avgfac; // Assign only those that match gen. info, if present at all.
            }
            if (
              (MYIDUP_tmp[2]==0 || MYIDUP_tmp[2]==ssel)
              &&
              (MYIDUP_tmp[3]==0 || MYIDUP_tmp[3]==rsel)
              ){
              __modhiggsjj_MOD_evalamp_sbfh_unsymm_sa_select_exact(p4, &isel, &jsel, &ssel, &rsel, MatElsq_tmp);
              MatElsq[jsel+5][isel+5] += MatElsq_tmp[jsel+5][isel+5]*avgfac; // Assign only those that match gen. info, if present at all.
            }
          }
          else{ // qqb->QQb
            if (abs(isel)!=1){ rsel=1; ssel=-1; } // Make sure rsel, ssel are not of same flavor as isel, jsel
            else{ rsel=2; ssel=-2; }
            // The amplitude is aready multiplied by nf-1, so no need to calculate everything (nf-1) times.
            double avgfac=1.; if (MYIDUP_tmp[2]==0 && MYIDUP_tmp[3]==0) avgfac=0.5;
            if (
              (MYIDUP_tmp[2]==0 || MYIDUP_tmp[2]==rsel)
              &&
              (MYIDUP_tmp[3]==0 || MYIDUP_tmp[3]==ssel)
              ){
              __modhiggsjj_MOD_evalamp_sbfh_unsymm_sa_select_exact(p4, &isel, &jsel, &rsel, &ssel, MatElsq_tmp);
              MatElsq[jsel+5][isel+5] += MatElsq_tmp[jsel+5][isel+5]*avgfac; // Assign only those that match gen. info, if present at all.
            }
            if (
              (MYIDUP_tmp[2]==0 || MYIDUP_tmp[2]==ssel)
              &&
              (MYIDUP_tmp[3]==0 || MYIDUP_tmp[3]==rsel)
              ){
              __modhiggsjj_MOD_evalamp_sbfh_unsymm_sa_select_exact(p4, &isel, &jsel, &ssel, &rsel, MatElsq_tmp);
              MatElsq[jsel+5][isel+5] += MatElsq_tmp[jsel+5][isel+5]*avgfac; // Assign only those that match gen. info, if present at all.
            }
          }
        }
        else{
          double avgfac=1.; if (MYIDUP_tmp[2]==0 && MYIDUP_tmp[3]==0 && rsel!=ssel) avgfac=0.5;
          if (
            (MYIDUP_tmp[2]==0 || MYIDUP_tmp[2]==rsel)
            &&
            (MYIDUP_tmp[3]==0 || MYIDUP_tmp[3]==ssel)
            ){
            __modhiggsjj_MOD_evalamp_sbfh_unsymm_sa_select_exact(p4, &isel, &jsel, &rsel, &ssel, MatElsq_tmp);
            MatElsq[jsel+5][isel+5] += MatElsq_tmp[jsel+5][isel+5]*avgfac; // Assign only those that match gen. info, if present at all.
          }
          if (
            rsel!=ssel
            &&
            (MYIDUP_tmp[2]==0 || MYIDUP_tmp[2]==ssel)
            &&
            (MYIDUP_tmp[3]==0 || MYIDUP_tmp[3]==rsel)
            ){
            __modhiggsjj_MOD_evalamp_sbfh_unsymm_sa_select_exact(p4, &isel, &jsel, &ssel, &rsel, MatElsq_tmp);
            MatElsq[jsel+5][isel+5] += MatElsq_tmp[jsel+5][isel+5]*avgfac; // Assign only those that match gen. info, if present at all.
          }
        }
      } // End swapped isel<jsel cases
    } // End loop over ic<nijchannels
  } // End production==TVar::JJGG
  else if (production==TVar::JJVBF){
    int ijsel[3][121];
    int nijchannels=68;
    __modhiggsjj_MOD_get_vbfchannelhash_nosplit(ijsel, &nijchannels);
    for (int ic=0; ic<nijchannels; ic++){
      // Emulate EvalWeighted_HJJ_test
      int isel = ijsel[0][ic];
      int jsel = ijsel[1][ic];
      int code = ijsel[2][ic];
      // Default assignments
      int rsel=isel;
      int ssel=jsel;
      if (
        (MYIDUP_tmp[0]==0 || ((PDGHelpers::isAGluon(MYIDUP_tmp[0]) && isel==0) || MYIDUP_tmp[0]==isel))
        &&
        (MYIDUP_tmp[1]==0 || ((PDGHelpers::isAGluon(MYIDUP_tmp[1]) && jsel==0) || MYIDUP_tmp[1]==jsel))
        ){ // Do it this way to be able to swap isel and jsel later
        if (code==1){ // Only ZZ->H possible
          // rsel=isel and ssel=jsel already
          double avgfac=1.; if (MYIDUP_tmp[2]==0 && MYIDUP_tmp[3]==0 && rsel!=ssel) avgfac=0.5;
          if (
            (MYIDUP_tmp[2]==0 || MYIDUP_tmp[2]==rsel)
            &&
            (MYIDUP_tmp[3]==0 || MYIDUP_tmp[3]==ssel)
            ){
            __modhiggsjj_MOD_evalamp_wbfh_unsymm_sa_select_exact(p4, &isel, &jsel, &rsel, &ssel, MatElsq_tmp);
            MatElsq[jsel+5][isel+5] += MatElsq_tmp[jsel+5][isel+5]*avgfac; // Assign only those that match gen. info, if present at all.
          }
          if (
            rsel!=ssel
            &&
            (MYIDUP_tmp[2]==0 || MYIDUP_tmp[2]==ssel)
            &&
            (MYIDUP_tmp[3]==0 || MYIDUP_tmp[3]==rsel)
            ){
            __modhiggsjj_MOD_evalamp_wbfh_unsymm_sa_select_exact(p4, &isel, &jsel, &ssel, &rsel, MatElsq_tmp);
            MatElsq[jsel+5][isel+5] += MatElsq_tmp[jsel+5][isel+5]*avgfac; // Assign only those that match gen. info, if present at all.
          }
        }
        else{ // code==0 means WW->H is also possible with no interference to ZZ->H, and code==2 means the same except interference to ZZ->H could occur for some outgoing quark flavors.
          vector<int> possible_rsel;
          vector<int> possible_ssel;
          if (PDGHelpers::isUpTypeQuark(isel)){ possible_rsel.push_back(1); possible_rsel.push_back(3); possible_rsel.push_back(5); }
          else if (PDGHelpers::isDownTypeQuark(isel)){ possible_rsel.push_back(2); possible_rsel.push_back(4); }
          if (PDGHelpers::isUpTypeQuark(jsel)){ possible_ssel.push_back(1); possible_ssel.push_back(3); possible_ssel.push_back(5); }
          else if (PDGHelpers::isDownTypeQuark(jsel)){ possible_ssel.push_back(2); possible_ssel.push_back(4); }
          for (unsigned int ix=0; ix<possible_rsel.size(); ix++){
            for (unsigned int iy=0; iy<possible_ssel.size(); iy++){
              rsel=possible_rsel.at(ix)*TMath::Sign(1, isel);
              ssel=possible_ssel.at(ix)*TMath::Sign(1, jsel);
              double avgfac=1.; if (MYIDUP_tmp[2]==0 && MYIDUP_tmp[3]==0 && rsel!=ssel) avgfac=0.5;
              if (
                (MYIDUP_tmp[2]==0 || MYIDUP_tmp[2]==rsel)
                &&
                (MYIDUP_tmp[3]==0 || MYIDUP_tmp[3]==ssel)
                ){
                __modhiggsjj_MOD_evalamp_wbfh_unsymm_sa_select_exact(p4, &isel, &jsel, &rsel, &ssel, MatElsq_tmp);
                MatElsq[jsel+5][isel+5] += MatElsq_tmp[jsel+5][isel+5]*avgfac; // Assign only those that match gen. info, if present at all.
              }
              if (
                rsel!=ssel
                &&
                (MYIDUP_tmp[2]==0 || MYIDUP_tmp[2]==ssel)
                &&
                (MYIDUP_tmp[3]==0 || MYIDUP_tmp[3]==rsel)
                ){
                __modhiggsjj_MOD_evalamp_wbfh_unsymm_sa_select_exact(p4, &isel, &jsel, &ssel, &rsel, MatElsq_tmp);
                MatElsq[jsel+5][isel+5] += MatElsq_tmp[jsel+5][isel+5]*avgfac; // Assign only those that match gen. info, if present at all.
              }
            }
          }
        }
      } // End unswapped isel>=jsel cases
      if (isel==jsel) continue;
      isel = ijsel[1][ic];
      jsel = ijsel[0][ic];
      // Reset to default assignments
      rsel=isel;
      ssel=jsel;
      if (
        (MYIDUP_tmp[0]==0 || ((PDGHelpers::isAGluon(MYIDUP_tmp[0]) && isel==0) || MYIDUP_tmp[0]==isel))
        &&
        (MYIDUP_tmp[1]==0 || ((PDGHelpers::isAGluon(MYIDUP_tmp[1]) && jsel==0) || MYIDUP_tmp[1]==jsel))
        ){
        // isel==jsel==0 is already eliminated by isel!=jsel condition
        if (code==1){ // Only ZZ->H possible
          // rsel=isel and ssel=jsel already
          double avgfac=1.; if (MYIDUP_tmp[2]==0 && MYIDUP_tmp[3]==0 && rsel!=ssel) avgfac=0.5;
          if (
            (MYIDUP_tmp[2]==0 || MYIDUP_tmp[2]==rsel)
            &&
            (MYIDUP_tmp[3]==0 || MYIDUP_tmp[3]==ssel)
            ){
            __modhiggsjj_MOD_evalamp_wbfh_unsymm_sa_select_exact(p4, &isel, &jsel, &rsel, &ssel, MatElsq_tmp);
            MatElsq[jsel+5][isel+5] += MatElsq_tmp[jsel+5][isel+5]*avgfac; // Assign only those that match gen. info, if present at all.
          }
          if (
            rsel!=ssel
            &&
            (MYIDUP_tmp[2]==0 || MYIDUP_tmp[2]==ssel)
            &&
            (MYIDUP_tmp[3]==0 || MYIDUP_tmp[3]==rsel)
            ){
            __modhiggsjj_MOD_evalamp_wbfh_unsymm_sa_select_exact(p4, &isel, &jsel, &ssel, &rsel, MatElsq_tmp);
            MatElsq[jsel+5][isel+5] += MatElsq_tmp[jsel+5][isel+5]*avgfac; // Assign only those that match gen. info, if present at all.
          }
        }
        else{ // code==0 means WW->H is also possible with no interference to ZZ->H, and code==2 means the same except interference to ZZ->H could occur for some outgoing quark flavors.
          vector<int> possible_rsel;
          vector<int> possible_ssel;
          if (PDGHelpers::isUpTypeQuark(isel)){ possible_rsel.push_back(1); possible_rsel.push_back(3); possible_rsel.push_back(5); }
          else if (PDGHelpers::isDownTypeQuark(isel)){ possible_rsel.push_back(2); possible_rsel.push_back(4); }
          if (PDGHelpers::isUpTypeQuark(jsel)){ possible_ssel.push_back(1); possible_ssel.push_back(3); possible_ssel.push_back(5); }
          else if (PDGHelpers::isDownTypeQuark(jsel)){ possible_ssel.push_back(2); possible_ssel.push_back(4); }
          for (unsigned int ix=0; ix<possible_rsel.size(); ix++){
            for (unsigned int iy=0; iy<possible_ssel.size(); iy++){
              rsel=possible_rsel.at(ix)*TMath::Sign(1, isel);
              ssel=possible_ssel.at(ix)*TMath::Sign(1, jsel);
              double avgfac=1.; if (MYIDUP_tmp[2]==0 && MYIDUP_tmp[3]==0 && rsel!=ssel) avgfac=0.5;
              if (
                (MYIDUP_tmp[2]==0 || MYIDUP_tmp[2]==rsel)
                &&
                (MYIDUP_tmp[3]==0 || MYIDUP_tmp[3]==ssel)
                ){
                __modhiggsjj_MOD_evalamp_wbfh_unsymm_sa_select_exact(p4, &isel, &jsel, &rsel, &ssel, MatElsq_tmp);
                MatElsq[jsel+5][isel+5] += MatElsq_tmp[jsel+5][isel+5]*avgfac; // Assign only those that match gen. info, if present at all.
              }
              if (
                rsel!=ssel
                &&
                (MYIDUP_tmp[2]==0 || MYIDUP_tmp[2]==ssel)
                &&
                (MYIDUP_tmp[3]==0 || MYIDUP_tmp[3]==rsel)
                ){
                __modhiggsjj_MOD_evalamp_wbfh_unsymm_sa_select_exact(p4, &isel, &jsel, &ssel, &rsel, MatElsq_tmp);
                MatElsq[jsel+5][isel+5] += MatElsq_tmp[jsel+5][isel+5]*avgfac; // Assign only those that match gen. info, if present at all.
              }
            }
          }
        }
      } // End swapped isel<jsel cases
    } // End loop over ic<nijchannels
  } // End production==TVar::JJVBF

  //    FOTRAN convention    -5    -4   -3   -2   -1    0   1   2   3  4  5
  //     parton flavor      bbar  cbar  sbar ubar dbar  g   d   u   s  c  b
  //      C++ convention     0      1    2    3    4    5   6   7   8  9  10
  for (int ii = 0; ii < nmsq; ii++){ for (int jj = 0; jj < nmsq; jj++){ if (verbosity >= TVar::DEBUG) cout<< "MatElsq: " << ii-5 << " " << jj-5 << " " << MatElsq[jj][ii] << endl; } }

  sum_msqjk = SumMEPDF(MomStore[0], MomStore[1], MatElsq, RcdME, EBEAM, verbosity);

  //cout << "Before reset: " << scale_.scale << '\t' << facscale_.facscale << endl;
  SetAlphaS(defaultRenScale, defaultFacScale, 1., 1., defaultNloop, defaultNflav, defaultPdflabel); // Protection for other probabilities
  //cout << "Default scale reset: " << scale_.scale << '\t' << facscale_.facscale << endl;
  return sum_msqjk;
}

double TUtil::VHiggsMatEl(
  TVar::Process process, TVar::Production production, TVar::MatrixElement matrixElement,
  event_scales_type* event_scales, MelaIO* RcdME,
  double EBEAM,
  TVar::VerbosityLevel verbosity
  ){
  const double GeV=1./100.; // JHUGen mom. scale factor
  double sum_msqjk = 0;
  // by default assume only gg productions
  // FOTRAN convention -5    -4   -3  -2    -1  0 1 2 3 4 5
  //     parton flavor bbar cbar sbar ubar dbar g d u s c b
  // C++ convention     0     1   2    3    4   5 6 7 8 9 10
  //2-D matrix is reversed in fortran
  // msq[ parton2 ] [ parton1 ]
  //      flavor_msq[jj][ii] = fx1[ii]*fx2[jj]*msq[jj][ii];
  double MatElsq[nmsq][nmsq]={ { 0 } };

  if (matrixElement!=TVar::JHUGen){ cerr << "TUtil::VHiggsMatEl: Non-JHUGen MEs are not supported" << endl; return sum_msqjk; }
  if (!(production == TVar::Lep_ZH || production == TVar::Lep_WH || production == TVar::Had_ZH || production == TVar::Had_WH)){ cerr << "TUtil::VHiggsMatEl: Production is not supported!" << endl; return sum_msqjk; }

  int nRequested_AssociatedJets=0;
  int nRequested_AssociatedLeptons=0;
  int AssociationVCompatibility=0;
  int partIncCode;
  if (production == TVar::Had_ZH || production == TVar::Had_WH){ // Only use associated partons
    partIncCode=TVar::kUseAssociated_Jets;
    nRequested_AssociatedJets=2;
  }
  else if (production == TVar::Lep_ZH || production == TVar::Lep_WH){ // Only use associated leptons(+)neutrinos
    partIncCode=TVar::kUseAssociated_Leptons;
    nRequested_AssociatedLeptons=2;
  }
  if (production==TVar::Lep_WH || production==TVar::Had_WH)) AssociationVCompatibility=24;
  else if (production==TVar::Lep_ZH || production==TVar::Had_ZH)) AssociationVCompatibility=23;
  simple_event_record mela_event;
  mela_event.AssociationCode=partIncCode;
  mela_event.AssociationVCompatibility=AssociationVCompatibility;
  mela_event.nRequested_AssociatedJets=nRequested_AssociatedJets;
  mela_event.nRequested_AssociatedLeptons=nRequested_AssociatedLeptons;
  GetBoostedParticleVectors(
    RcdME->melaCand,
    mela_event,
    partIncCode
    );
  if (mela_event.pAssociated.size()<2){ cerr << "TUtil::VHiggsMatEl: Number of associated particles is 0!" << endl; return sum_msqjk; }

  int MYIDUP_prod[4]={ 0 }; // "Incoming" partons 1, 2, "outgoing" partons 3, 4
  int MYIDUP_dec[2]={ 0 }; // "Outgoing" partons 1, 2
  double p4[9][4] ={ { 0 } };
  double helicities[9] ={ 0 };
  int vh_ids[9] ={ 0 };
  TLorentzVector MomStore[mxpart];
  for (int i = 0; i < mxpart; i++) MomStore[i].SetXYZT(0, 0, 0, 0);

  // p4(0:8,i) = (E(i),px(i),py(i),pz(i))
  // i=0,1: q1, qb2 (outgoing convention)
  // i=2,3: V*, V
  // i=4: H
  // i=5,6: f, fb from V
  // i=7,8: b, bb from H
  for (int ipar=0; ipar<2; ipar++){
    TLorentzVector* momTmp = &(mela_event.pMothers.at(ipar).second);
    int* idtmp = &(mela_event.pMothers.at(ipar).first);
    if (!PDGHelpers::isAnUnknownJet(*idtmp)) MYIDUP_prod[ipar] = *idtmp;
    else MYIDUP_prod[ipar] = 0;
    if (momTmp->T()>0.){
      p4[ipar][0] = momTmp->T()*GeV;
      p4[ipar][1] = momTmp->X()*GeV;
      p4[ipar][2] = momTmp->Y()*GeV;
      p4[ipar][3] = momTmp->Z()*GeV;
      MomStore[ipar] = (*momTmp);
    }
    else{
      p4[ipar][0] = -momTmp->T()*GeV;
      p4[ipar][1] = -momTmp->X()*GeV;
      p4[ipar][2] = -momTmp->Y()*GeV;
      p4[ipar][3] = -momTmp->Z()*GeV;
      MomStore[ipar] = -(*momTmp);
      MYIDUP_prod[ipar] = -MYIDUP_prod[ipar];
    }
  }
  // Associated particles
  for (int ipar=0; ipar<2; ipar++){
    TLorentzVector* momTmp = &(mela_event.pAssociated.at(ipar).second);
    int* idtmp = &(mela_event.pAssociated.at(ipar).first);
    if (!PDGHelpers::isAnUnknownJet(*idtmp)) MYIDUP_prod[ipar+2] = *idtmp;
    else MYIDUP_prod[ipar+2] = 0;
    p4[ipar+5][0] = momTmp->T()*GeV;
    p4[ipar+5][1] = momTmp->X()*GeV;
    p4[ipar+5][2] = momTmp->Y()*GeV;
    p4[ipar+5][3] = momTmp->Z()*GeV;
    MomStore[ipar+6] = (*momTmp);
  }

  if (PDGHelpers::isAGluon(MYIDUP_prod[0]) || PDGHelpers::isAGluon(MYIDUP_prod[1])){ cerr << "TUtil::VHiggsMatEl: Initial state gluons are not permitted!" << endl; return sum_msqjk; }
  if (PDGHelpers::isAGluon(MYIDUP_prod[2]) || PDGHelpers::isAGluon(MYIDUP_prod[3])){ cerr << "TUtil::VHiggsMatEl: Final state gluons are not permitted!" << endl; return sum_msqjk; }

  // Decay V/f ids
  for (int iv=0; iv<2; iv++){
    int idtmp = mela_event.intermediateVid.at(iv);
    if (!PDGHelpers::isAnUnknownJet(idtmp)) MYIDUP_dec[iv] = idtmp;
    else MYIDUP_dec[iv] = 0;
  }
  // Decay daughters
  for (int ipar=0; ipar<mela_event.pDaughters.size(); ipar++){
    TLorentzVector* momTmp = &(mela_event.pDaughters.at(ipar).second);
    if (mela_event.pDaughters.size()==1){
      p4[ipar+7][0] = momTmp->T()*GeV;
      p4[ipar+7][1] = momTmp->X()*GeV;
      p4[ipar+7][2] = momTmp->Y()*GeV;
      p4[ipar+7][3] = momTmp->Z()*GeV;
      MomStore[5] = (*momTmp); // 5
    }
    else if (mela_event.pDaughters.size()==2){
      p4[ipar+7][0] = momTmp->T()*GeV;
      p4[ipar+7][1] = momTmp->X()*GeV;
      p4[ipar+7][2] = momTmp->Y()*GeV;
      p4[ipar+7][3] = momTmp->Z()*GeV;
      MomStore[2*ipar+2] = (*momTmp); // 2,4
    }
    else if (mela_event.pDaughters.size()==3){
      if (ipar<2){
        p4[7][0] += momTmp->T()*GeV;
        p4[7][1] += momTmp->X()*GeV;
        p4[7][2] += momTmp->Y()*GeV;
        p4[7][3] += momTmp->Z()*GeV;
      }
      else{
        p4[8][0] = momTmp->T()*GeV;
        p4[8][1] = momTmp->X()*GeV;
        p4[8][2] = momTmp->Y()*GeV;
        p4[8][3] = momTmp->Z()*GeV;
      }
      MomStore[ipar+2] = (*momTmp); // 2,3,4
    }
    else if (mela_event.pDaughters.size()==4){
      if (ipar<2){
        p4[7][0] += momTmp->T()*GeV;
        p4[7][1] += momTmp->X()*GeV;
        p4[7][2] += momTmp->Y()*GeV;
        p4[7][3] += momTmp->Z()*GeV;
      }
      else{
        p4[8][0] += momTmp->T()*GeV;
        p4[8][1] += momTmp->X()*GeV;
        p4[8][2] += momTmp->Y()*GeV;
        p4[8][3] += momTmp->Z()*GeV;
      }
      MomStore[ipar+2] = (*momTmp); // 2,3,4,5
    }
    else{ // Should never happen
      p4[7][0] += momTmp->T()*GeV;
      p4[7][1] += momTmp->X()*GeV;
      p4[7][2] += momTmp->Y()*GeV;
      p4[7][3] += momTmp->Z()*GeV;
      MomStore[5] = MomStore[5] + (*momTmp);
    }
  }
  for (int ix=0; ix<4; ix++){
    p4[3][ix] = p4[5][ix] + p4[6][ix];
    p4[4][ix] = p4[7][ix] + p4[8][ix];
    p4[2][ix] = p4[3][ix] + p4[4][ix];
  }

  vh_ids[4] = 25;
  if (production==TVar::Lep_ZH || production==TVar::Had_ZH) vh_ids[2] = 23;
  else if (production==TVar::Lep_WH || production==TVar::Had_WH) vh_ids[2] = 24;
  vh_ids[3] = vh_ids[2];

  // H->ffb decay is turned off, so no need to loop over helicities[7]=helicities[8]=+-1
  vh_ids[7] = 5; helicities[7] = 1;
  vh_ids[8] = 5; helicities[8] = 1;

  if ( verbosity >= TVar::DEBUG ) {
    for (int i=0; i<9; i++) cout << "p4[0] = "  << p4[i][0] << ", " <<  p4[i][1] << ", "  <<  p4[i][2] << ", "  <<  p4[i][3] << "\n";
    for (int i=0; i<9; i++) cout << "id(" << i << ") = "  << vh_ids[i] << endl;
  }

  double defaultRenScale = scale_.scale;
  double defaultFacScale = facscale_.facscale;
  //cout << "Default scales: " << defaultRenScale << '\t' << defaultFacScale << endl;
  int defaultNloop = nlooprun_.nlooprun;
  int defaultNflav = nflav_.nflav;
  string defaultPdflabel = pdlabel_.pdlabel;
  double renQ = InterpretScaleScheme(production, matrixElement, event_scales->renomalizationScheme, MomStore);
  //cout << "renQ: " << renQ << " x " << event_scales->ren_scale_factor << endl;
  double facQ = InterpretScaleScheme(production, matrixElement, event_scales->factorizationScheme, MomStore);
  //cout << "facQ: " << facQ << " x " << event_scales->fac_scale_factor << endl;
  SetAlphaS(renQ, facQ, event_scales->ren_scale_factor, event_scales->fac_scale_factor, 1, 5, "cteq6_l"); // Set AlphaS(|Q|/2, mynloop, mynflav, mypartonPDF) for MCFM ME-related calculations

  const double allowed_helicities[2] = { -1, 1 }; // L,R
  double sumME=0;
  for (int h01 = 0; h01 < 2; h01++){
    helicities[0] = allowed_helicities[h01];
    helicities[1] = -helicities[0];
    for (int h56 = 0; h56 < 2; h56++){
      helicities[5] = allowed_helicities[h56];
      helicities[6] = -helicities[5];
      for (int incoming1 = -nf; incoming1 <= nf; incoming1++){
        if (incoming1==0) continue;
        if (production==TVar::Lep_ZH || production==TVar::Had_ZH){
          vh_ids[0] = incoming1;
          vh_ids[1] = -incoming1;
          if (
            (MYIDUP_prod[0]!=0 && MYIDUP_prod[0]!=vh_ids[0])
            ||
            (MYIDUP_prod[1]!=0 && MYIDUP_prod[1]!=vh_ids[1])
            ) continue;

          if (production==TVar::Had_ZH){
            for (int outgoing1=-nf; outgoing1<=nf; outgoing1++){
              vh_ids[5] = outgoing1;
              vh_ids[6] = -outgoing1;
              if (
                (MYIDUP_prod[2]!=0 && MYIDUP_prod[2]!=vh_ids[5])
                ||
                (MYIDUP_prod[3]!=0 && MYIDUP_prod[3]!=vh_ids[6])
                ) continue;

              double msq=0;
              __modvhiggs_MOD_evalamp_vhiggs(vh_ids, helicities, p4, &msq);
              MatElsq[vh_ids[0]+5][vh_ids[1]+5] += msq * 0.25; // Average over initial states with helicities +-1 only
            } // End loop over outgoing1=-outgoing2
          }
          else{
            vh_ids[5] = MYIDUP_prod[2];
            vh_ids[6] = MYIDUP_prod[3];

            double msq=0;
            __modvhiggs_MOD_evalamp_vhiggs(vh_ids, helicities, p4, &msq);
            MatElsq[vh_ids[0]+5][vh_ids[1]+5] += msq * 0.25; // Average over initial states with helicities +-1 only
          } // End check on Had vs Lep
        } // End ZH case
        else if (production==TVar::Lep_WH || production==TVar::Had_WH){
          vh_ids[0] = incoming1;
          if (MYIDUP_prod[0]!=0 && MYIDUP_prod[0]!=vh_ids[0]) continue;

          for (int incoming2 = -nf; incoming2 < 0; incoming2++){
            if (abs(incoming2)==abs(incoming1) || TMath::Sign(1, incoming1)==TMath::Sign(1, incoming2)) continue;

            vh_ids[1] = incoming2;
            if (MYIDUP_prod[1]!=0 && MYIDUP_prod[1]!=vh_ids[1]) continue;

            if (production==TVar::Had_WH){
              for (int outgoing1=-nf; outgoing1<=nf; outgoing1++){
                for (int outgoing2=-nf; outgoing2<=nf; outgoing2++){
                  if (abs(outgoing2)==abs(outgoing1) || TMath::Sign(1, outgoing1)==TMath::Sign(1, outgoing2)) continue;

                  // Determine whether the decay is from a W+ or a W-
                  if (
                    (PDGHelpers::isUpTypeQuark(outgoing1) && outgoing1>0)
                    ||
                    (PDGHelpers::isUpTypeQuark(outgoing2) && outgoing2>0)
                    ) vh_ids[2]=24;
                  else vh_ids[2]=-24;

                  vh_ids[5] = outgoing1;
                  vh_ids[6] = outgoing2;
                  if (
                    (MYIDUP_prod[2]!=0 && MYIDUP_prod[2]!=vh_ids[5])
                    ||
                    (MYIDUP_prod[3]!=0 && MYIDUP_prod[3]!=vh_ids[6])
                    ) continue;

                  double msq=0;
                  __modvhiggs_MOD_evalamp_vhiggs(vh_ids, helicities, p4, &msq);
                  MatElsq[vh_ids[0]+5][vh_ids[1]+5] += msq * 0.25; // Average over initial states with helicities +-1 only
                } // End loop over outgoing2
              } // End loop over outgoing1
            }
            else{
              // Determine whether the decay is from a W+ or a W-
              if (
                (PDGHelpers::isANeutrino(outgoing1) && outgoing1>0)
                ||
                (PDGHelpers::isANeutrino(outgoing2) && outgoing2>0)
                ) vh_ids[2]=24;
              else vh_ids[2]=-24;

              vh_ids[5] = MYIDUP_prod[2];
              vh_ids[6] = MYIDUP_prod[3];

              double msq=0;
              __modvhiggs_MOD_evalamp_vhiggs(vh_ids, helicities, p4, &msq);
              MatElsq[vh_ids[0]+5][vh_ids[1]+5] += msq * 0.25; // Average over initial states with helicities +-1 only
            } // End check on Had vs Lep
          } // End loop over incoming2
        } // End WH case
      } // End loop over incoming1
    } // End loop over h56
  } // End loop over h01

  sum_msqjk = SumMEPDF(MomStore[0], MomStore[1], MatElsq, RcdME, EBEAM, verbosity);

  //cout << "Before reset: " << scale_.scale << '\t' << facscale_.facscale << endl;
  SetAlphaS(defaultRenScale, defaultFacScale, 1., 1., defaultNloop, defaultNflav, defaultPdflabel); // Protection for other probabilities
  //cout << "Default scale reset: " << scale_.scale << '\t' << facscale_.facscale << endl;
  return sum_msqjk;
}


// ttH
double TUtil::TTHiggsMatEl(
  TVar::Process process, TVar::Production production, TVar::MatrixElement matrixElement,
  event_scales_type* event_scales, MelaIO* RcdME,
  double EBEAM,
  int topDecay, int topProcess,
  TVar::VerbosityLevel verbosity
  ){
  const double GeV=1./100.; // JHUGen mom. scale factor
  double sum_msqjk = 0;
  double MatElsq[nmsq][nmsq]={ { 0 } };
  double MatElsq_tmp[nmsq][nmsq]={ { 0 } };

  if (matrixElement!=TVar::JHUGen){ cerr << "TUtil::TTHiggsMatEl: Non-JHUGen MEs are not supported." << endl; return sum_msqjk; }
  if (production!=TVar::ttH){ cerr << "TUtil::TTHiggsMatEl: Only ttH is supported." << endl; return sum_msqjk; }

  int partIncCode;
  int nRequested_Tops=1;
  int nRequested_Antitops=1;
  if (topDecay>0) partIncCode=TVar::kUseAssociated_UnstableTops; // Look for unstable tops
  else partIncCode=TVar::kUseAssociated_StableTops; // Look for unstable tops

  simple_event_record mela_event;
  mela_event.AssociationCode=partIncCode;
  mela_event.nRequested_Tops=nRequested_Tops;
  mela_event.nRequested_Antitops=nRequested_Antitops;
  GetBoostedParticleVectors(
    RcdME->melaCand,
    mela_event
    );


  if (topDecay>0 && mela_event.pTopDaughters.size()<1 && mela_event.pAntitopDaughters.size()<1){
    cerr
      << "TUtil::TTHiggsMatEl: Number of set of top daughters (" << mela_event.pTopDaughters.size() << ")"
      << "and number of set of antitop daughters (" << mela_event.pAntitopDaughters.size() << ")"
      <<" in ttH process is not 1!" << endl;
    return sum_msqjk;
  }
  else if (topDecay>0 && mela_event.pTopDaughters.at(0).size()!=3 && mela_event.pAntitopDaughters.at(0).size()!=3){
    cerr
      << "TUtil::TTHiggsMatEl: Number of top daughters (" << mela_event.pTopDaughters.at(0).size() << ")"
      << "and number of antitop daughters (" << mela_event.pAntitopDaughters.at(0).size() << ")"
      <<" in ttH process is not 3!" << endl;
    return sum_msqjk;
  }
  else if (topDecay==0 && mela_event.pStableTops.size()<1 && mela_event.pStableAntitops.size()<1){
    cerr
      << "TUtil::TTHiggsMatEl: Number of stable tops (" << mela_event.pStableTops.size() << ")"
      << "and number of sstable antitops (" << mela_event.pStableAntitops.size() << ")"
      <<" in ttH process is not 1!" << endl;
    return sum_msqjk;
  }

  SimpleParticleCollection_t topDaughters;
  SimpleParticleCollection_t topbarDaughters;
  bool isUnknown[2]; isUnknown[0]=true; isUnknown[1]=true;

  if (topDecay>0){
    // Daughters are assumed to have been ordered as b, Wf, Wfb already.
    for (int itd=0; itd<mela_event.pTopDaughters.at(0).size(); itd++) topDaughter.push_back(mela_event.pTopDaughters.at(0).at(itd));
    for (int itd=0; itd<mela_event.pAntitopDaughters.at(0).size(); itd++) topbarDaughters.push_back(mela_event.pAntitopDaughters.at(0).at(itd));
  }
  else{
    for (int itop=0; itop<mela_event.pStableTops.size(); itop++) topDaughter.push_back(mela_event.pStableTops.at(itop));
    for (int itop=0; itop<mela_event.pStableAntitops.size(); itop++) topbarDaughters.push_back(mela_event.pStableAntitops.at(itop));
  }
  // Check if either top is definitely identified
  for (int itd=0; itd<topDaughter.size(); itd++){ if (topDaughter.at(itd).first!=0){ isUnknown[0]=false; break; } }
  for (int itd=0; itd<topbarDaughters.size(); itd++){ if (topbarDaughters.at(itd).first!=0){ isUnknown[1]=false; break; } }

  // Start assigning the momenta
  // 0,1: p1 p2
  // 2-4: H,tb,t
  // 5-8: bb,W-,f,fb
  // 9-12: b,W+,fb,f
  double p4[13][4]={ { 0 } };
  MYIDUP_prod[2]={ 0 };
  TLorentzVector MomStore[mxpart];
  for (int i = 0; i < mxpart; i++) MomStore[i].SetXYZT(0, 0, 0, 0);
  for (int ipar=0; ipar<2; ipar++){
    TLorentzVector* momTmp = &(mela_event.pMothers.at(ipar).second);
    int* idtmp = &(mela_event.pMothers.at(ipar).first);
    if (!PDGHelpers::isAnUnknownJet(*idtmp)) MYIDUP_prod[ipar] = *idtmp;
    else MYIDUP_prod[ipar] = 0;
    if (momTmp->T()>0.){
      p4[ipar][0] = -momTmp->T()*GeV;
      p4[ipar][1] = -momTmp->X()*GeV;
      p4[ipar][2] = -momTmp->Y()*GeV;
      p4[ipar][3] = -momTmp->Z()*GeV;
      MomStore[ipar] = (*momTmp);
    }
    else{
      p4[ipar][0] = momTmp->T()*GeV;
      p4[ipar][1] = momTmp->X()*GeV;
      p4[ipar][2] = momTmp->Y()*GeV;
      p4[ipar][3] = momTmp->Z()*GeV;
      MomStore[ipar] = -(*momTmp);
      MYIDUP_prod[ipar] = -MYIDUP_prod[ipar];
    }
  }

  // Assign top momenta
  for (unsigned int ipar=0; ipar<topDaughter.size(); ipar++){
    TLorentzVector* momTmp = &(topDaughter.at(ipar).second);
    if (topDaughter.size()==1){
      p4[4][0] = momTmp->T()*GeV;
      p4[4][1] = momTmp->X()*GeV;
      p4[4][2] = momTmp->Y()*GeV;
      p4[4][3] = momTmp->Z()*GeV;
    }
    // size==3
    else if (ipar==0){ // b
      p4[9][0] = momTmp->T()*GeV;
      p4[9][1] = momTmp->X()*GeV;
      p4[9][2] = momTmp->Y()*GeV;
      p4[9][3] = momTmp->Z()*GeV;
    }
    else if (ipar==2){ // Wfb
      p4[11][0] = momTmp->T()*GeV;
      p4[11][1] = momTmp->X()*GeV;
      p4[11][2] = momTmp->Y()*GeV;
      p4[11][3] = momTmp->Z()*GeV;
      p4[10][0] += p4[11][0];
      p4[10][1] += p4[11][1];
      p4[10][2] += p4[11][2];
      p4[10][3] += p4[11][3];
    }
    else/* if (ipar==1)*/{ // Wf
      p4[12][0] = momTmp->T()*GeV;
      p4[12][1] = momTmp->X()*GeV;
      p4[12][2] = momTmp->Y()*GeV;
      p4[12][3] = momTmp->Z()*GeV;
      p4[10][0] += p4[12][0];
      p4[10][1] += p4[12][1];
      p4[10][2] += p4[12][2];
      p4[10][3] += p4[12][3];
    }
    MomStore[6] = MomStore[6] + (*momTmp); // MomStore (I1, I2, 0, 0, 0, H, J1, J2)
  }
  if (topDaughter.size()!=1){ for (int ix=0; ix<4; ix++){ for (int ip=9; ip<=10; ip++) p4[4][ix] = p4[ip][ix]; } }

  // Assign antitop momenta
  for (unsigned int ipar=0; ipar<antitopDaughter.size(); ipar++){
    TLorentzVector* momTmp = &(antitopDaughter.at(ipar).second);
    if (antitopDaughter.size()==1){
      p4[3][0] = momTmp->T()*GeV;
      p4[3][1] = momTmp->X()*GeV;
      p4[3][2] = momTmp->Y()*GeV;
      p4[3][3] = momTmp->Z()*GeV;
    }
    // size==3
    else if (ipar==0){ // bb
      p4[5][0] = momTmp->T()*GeV;
      p4[5][1] = momTmp->X()*GeV;
      p4[5][2] = momTmp->Y()*GeV;
      p4[5][3] = momTmp->Z()*GeV;
    }
    else if (ipar==1){ // Wf
      p4[7][0] = momTmp->T()*GeV;
      p4[7][1] = momTmp->X()*GeV;
      p4[7][2] = momTmp->Y()*GeV;
      p4[7][3] = momTmp->Z()*GeV;
      p4[6][0] += p4[7][0];
      p4[6][1] += p4[7][1];
      p4[6][2] += p4[7][2];
      p4[6][3] += p4[7][3];
    }
    else/* if (ipar==1)*/{ // Wfb
      p4[8][0] = momTmp->T()*GeV;
      p4[8][1] = momTmp->X()*GeV;
      p4[8][2] = momTmp->Y()*GeV;
      p4[8][3] = momTmp->Z()*GeV;
      p4[6][0] += p4[8][0];
      p4[6][1] += p4[8][1];
      p4[6][2] += p4[8][2];
      p4[6][3] += p4[8][3];
    }
    MomStore[7] = MomStore[7] + (*momTmp); // MomStore (I1, I2, 0, 0, 0, H, J1, J2)
  }
  if (antitopDaughter.size()!=1){ for (int ix=0; ix<4; ix++){ for (int ip=5; ip<=6; ip++) p4[3][ix] = p4[ip][ix]; } }

  for (int ipar=0; ipar<mela_event.pDaughters.size(); ipar++){
    TLorentzVector* momTmp = &(mela_event.pDaughters.at(ipar).second);
    p4[2][0] += momTmp->T()*GeV;
    p4[2][1] += momTmp->X()*GeV;
    p4[2][2] += momTmp->Y()*GeV;
    p4[2][3] += momTmp->Z()*GeV;
    MomStore[5] = MomStore[5] + (*momTmp); // i==(2, 3, 4) is (J1, J2, H), recorded as MomStore (I1, I2, 0, 0, 0, H, J1, J2)
  }

  if (verbosity >= TVar::DEBUG){
    for (int ii=0; ii<13; ii++){ cout << "p4[" << ii << "] = "; for (int jj=0; jj<4; jj++) cout << p4[ii][jj]/GeV << '\t'; cout << endl; }
  }

  double defaultRenScale = scale_.scale;
  double defaultFacScale = facscale_.facscale;
  //cout << "Default scales: " << defaultRenScale << '\t' << defaultFacScale << endl;
  int defaultNloop = nlooprun_.nlooprun;
  int defaultNflav = nflav_.nflav;
  string defaultPdflabel = pdlabel_.pdlabel;
  double renQ = InterpretScaleScheme(production, matrixElement, event_scales->renomalizationScheme, MomStore);
  //cout << "renQ: " << renQ << " x " << event_scales->ren_scale_factor << endl;
  double facQ = InterpretScaleScheme(production, matrixElement, event_scales->factorizationScheme, MomStore);
  //cout << "facQ: " << facQ << " x " << event_scales->fac_scale_factor << endl;
  SetAlphaS(renQ, facQ, event_scales->ren_scale_factor, event_scales->fac_scale_factor, 1, 5, "cteq6_l"); // Set AlphaS(|Q|/2, mynloop, mynflav, mypartonPDF) for MCFM ME-related calculations

  __modjhugenmela_MOD_settopdecays(&topDecay);
  __modttbh_MOD_evalxsec_pp_ttbh(p4, &topProcess, MatElsq);
  if (isUnknown[0] && isUnknown[1]){
    Swap_Momenta(p4[3], p4[4]);
    Swap_Momenta(p4[5], p4[9]);
    Swap_Momenta(p4[6], p4[10]);
    Swap_Momenta(p4[7], p4[12]);
    Swap_Momenta(p4[8], p4[11]);
    __modttbh_MOD_evalxsec_pp_ttbh(p4, &topProcess, MatElsq_tmp);
    for (int ix=0; ix<11; ix++){ for (int iy=0; iy<11; iy++) MatElsq[iy][ix] = (MatElsq[iy][ix]+MatElsq_tmp[iy][ix])/2.; }
  }
  int defaultTopDecay=-1;
  __modjhugenmela_MOD_settopdecays(&defaultTopDecay); // reset top decay

  sum_msqjk = SumMEPDF(MomStore[0], MomStore[1], MatElsq, RcdME, EBEAM, verbosity);

  //cout << "Before reset: " << scale_.scale << '\t' << facscale_.facscale << endl;
  SetAlphaS(defaultRenScale, defaultFacScale, 1., 1., defaultNloop, defaultNflav, defaultPdflabel); // Protection for other probabilities
  //cout << "Default scale reset: " << scale_.scale << '\t' << facscale_.facscale << endl;
  return sum_msqjk;
}

// bbH
double TUtil::BBHiggsMatEl(
  TVar::Process process, TVar::Production production, TVar::MatrixElement matrixElement,
  event_scales_type* event_scales, MelaIO* RcdME,
  double EBEAM,
  int botProcess,
  TVar::VerbosityLevel verbosity
  ){
  const double GeV=1./100.; // JHUGen mom. scale factor
  double sum_msqjk = 0;
  double MatElsq[nmsq][nmsq]={ { 0 } };

  if (matrixElement!=TVar::JHUGen){ cerr << "TUtil::BBHiggsMatEl: Non-JHUGen MEs are not supported." << endl; return sum_msqjk; }
  if (production!=TVar::bbH){ cerr << "TUtil::BBHiggsMatEl: Only bbH is supported." << endl; return sum_msqjk; }

  int partIncCode=TVar::kUseAssociated_Jets; // Look for jets
  int nRequested_AssociatedJets=2;

  simple_event_record mela_event;
  mela_event.AssociationCode=partIncCode;
  mela_event.kUseAssociated_Jets=nRequested_Tops;
  GetBoostedParticleVectors(
    RcdME->melaCand,
    mela_event
    );

  if (mela_event.pAssociated.size()<2){
    cerr
      << "TUtil::BBHiggsMatEl: Number of stable bs (" << mela_event.pAssociated.size() << ")"
      <<" in bbH process is not 2!" << endl;
    return sum_msqjk;
  }

  SimpleParticleCollection_t topDaughters;
  SimpleParticleCollection_t topbarDaughters;
  bool isUnknown[2]; isUnknown[0]=false; isUnknown[1]=false;

  if (pAssociated.at(0).first>=0){ topDaughter.push_back(mela_event.pAssociated.at(0)); isUnknown[0]=(PDGHelpers::isAnUnknownJet(pAssociated.at(0).first)); }
  else antitopDaughter.push_back(mela_event.pAssociated.at(0));
  if (pAssociated.at(1).first<=0){ antitopDaughter.push_back(mela_event.pAssociated.at(1)); isUnknown[1]=(PDGHelpers::isAnUnknownJet(pAssociated.at(1).first)); }
  else topDaughter.push_back(mela_event.pAssociated.at(1));

  // Start assigning the momenta
  // 0,1: p1 p2
  // 2-4: H,tb,t
  // 5-8: bb,W-,f,fb
  // 9-12: b,W+,fb,f
  double p4[13][4]={ { 0 } };
  MYIDUP_prod[2]={ 0 };
  TLorentzVector MomStore[mxpart];
  for (int i = 0; i < mxpart; i++) MomStore[i].SetXYZT(0, 0, 0, 0);
  for (int ipar=0; ipar<2; ipar++){
    TLorentzVector* momTmp = &(mela_event.pMothers.at(ipar).second);
    int* idtmp = &(mela_event.pMothers.at(ipar).first);
    if (!PDGHelpers::isAnUnknownJet(*idtmp)) MYIDUP_prod[ipar] = *idtmp;
    else MYIDUP_prod[ipar] = 0;
    if (momTmp->T()>0.){
      p4[ipar][0] = -momTmp->T()*GeV;
      p4[ipar][1] = -momTmp->X()*GeV;
      p4[ipar][2] = -momTmp->Y()*GeV;
      p4[ipar][3] = -momTmp->Z()*GeV;
      MomStore[ipar] = (*momTmp);
    }
    else{
      p4[ipar][0] = momTmp->T()*GeV;
      p4[ipar][1] = momTmp->X()*GeV;
      p4[ipar][2] = momTmp->Y()*GeV;
      p4[ipar][3] = momTmp->Z()*GeV;
      MomStore[ipar] = -(*momTmp);
      MYIDUP_prod[ipar] = -MYIDUP_prod[ipar];
    }
  }

  // Assign b momenta
  for (unsigned int ipar=0; ipar<topDaughter.size(); ipar++){
    TLorentzVector* momTmp = &(topDaughter.at(ipar).second);
    p4[4][0] = momTmp->T()*GeV;
    p4[4][1] = momTmp->X()*GeV;
    p4[4][2] = momTmp->Y()*GeV;
    p4[4][3] = momTmp->Z()*GeV;
    MomStore[6] = MomStore[6] + (*momTmp); // MomStore (I1, I2, 0, 0, 0, H, J1, J2)
  }

  // Assign bb momenta
  for (unsigned int ipar=0; ipar<antitopDaughter.size(); ipar++){
    TLorentzVector* momTmp = &(antitopDaughter.at(ipar).second);
    p4[3][0] = momTmp->T()*GeV;
    p4[3][1] = momTmp->X()*GeV;
    p4[3][2] = momTmp->Y()*GeV;
    p4[3][3] = momTmp->Z()*GeV;
    MomStore[7] = MomStore[7] + (*momTmp); // MomStore (I1, I2, 0, 0, 0, H, J1, J2)
  }

  for (int ipar=0; ipar<mela_event.pDaughters.size(); ipar++){
    TLorentzVector* momTmp = &(mela_event.pDaughters.at(ipar).second);
    p4[2][0] += momTmp->T()*GeV;
    p4[2][1] += momTmp->X()*GeV;
    p4[2][2] += momTmp->Y()*GeV;
    p4[2][3] += momTmp->Z()*GeV;
    MomStore[5] = MomStore[5] + (*momTmp); // i==(2, 3, 4) is (J1, J2, H), recorded as MomStore (I1, I2, 0, 0, 0, H, J1, J2)
  }

  if (verbosity >= TVar::DEBUG){
    for (int ii=0; ii<13; ii++){ cout << "p4[" << ii << "] = "; for (int jj=0; jj<4; jj++) cout << p4[ii][jj]/GeV << '\t'; cout << endl; }
  }

  double defaultRenScale = scale_.scale;
  double defaultFacScale = facscale_.facscale;
  //cout << "Default scales: " << defaultRenScale << '\t' << defaultFacScale << endl;
  int defaultNloop = nlooprun_.nlooprun;
  int defaultNflav = nflav_.nflav;
  string defaultPdflabel = pdlabel_.pdlabel;
  double renQ = InterpretScaleScheme(production, matrixElement, event_scales->renomalizationScheme, MomStore);
  //cout << "renQ: " << renQ << " x " << event_scales->ren_scale_factor << endl;
  double facQ = InterpretScaleScheme(production, matrixElement, event_scales->factorizationScheme, MomStore);
  //cout << "facQ: " << facQ << " x " << event_scales->fac_scale_factor << endl;
  SetAlphaS(renQ, facQ, event_scales->ren_scale_factor, event_scales->fac_scale_factor, 1, 5, "cteq6_l"); // Set AlphaS(|Q|/2, mynloop, mynflav, mypartonPDF) for MCFM ME-related calculations

  __modbbbh_MOD_evalxsec_pp_bbbh(p4, &botProcess, MatElsq);
  if (isUnknown[0] && isUnknown[1]){
    Swap_Momenta(p4[3], p4[4]);
    __modttbh_MOD_evalxsec_pp_bbbh(p4, &botProcess, MatElsq_tmp);
    for (int ix=0; ix<11; ix++){ for (int iy=0; iy<11; iy++) MatElsq[iy][ix] = (MatElsq[iy][ix]+MatElsq_tmp[iy][ix])/2.; }
  }
  sum_msqjk = SumMEPDF(MomStore[0], MomStore[1], MatElsq, RcdME, EBEAM, verbosity);

  //cout << "Before reset: " << scale_.scale << '\t' << facscale_.facscale << endl;
  SetAlphaS(defaultRenScale, defaultFacScale, 1., 1., defaultNloop, defaultNflav, defaultPdflabel); // Protection for other probabilities
  //cout << "Default scale reset: " << scale_.scale << '\t' << facscale_.facscale << endl;
  return sum_msqjk;
}

void TUtil::Swap_Momenta(double(&p)[4], double(&q)[4]){
  for (int ix=0; ix<4; ix++){
    double tmp=p[ix];
    p[ix]=q[ix];
    q[ix]=tmp;
  }
}

// CheckPartonMomFraction computes xx[0:1] based on p0, p1
bool TUtil::CheckPartonMomFraction(const TLorentzVector p0, const TLorentzVector p1, double xx[2], double EBEAM, TVar::VerbosityLevel verbosity){
  //Make sure parton Level Energy fraction is [0,1]
  //phase space function already makes sure the parton energy fraction between [min,1]
  //  x0 EBeam =>   <= -x1 EBeam
  double sysPz=p0.Pz()    + p1.Pz();
  double sysE =p0.Energy()+ p1.Energy();

  //Ignore the Pt doesn't make significant effect
  //double sysPt_sqr=sysPx*sysPx+sysPy*sysPy;
  //if(sysPt_sqr>=1.0E-10)  sysE=TMath::Sqrt(sysE*sysE-sysPt_sqr);
  xx[0]=(sysE+sysPz)/EBEAM/2.;
  xx[1]=(sysE-sysPz)/EBEAM/2.;
  if (verbosity >= TVar::DEBUG) cout << "xx[0]: " << xx[0] << ", xx[1] = " << xx[1] << '\n';

  if (
    xx[0] > 1.0 || xx[0]<=xmin_.xmin
    ||
    xx[1] > 1.0 || xx[1]<=xmin_.xmin
    ) return false;
  else return true;
}
// ComputePDF does the PDF computation
void TUtil::ComputePDF(const TLorentzVector p0, const TLorentzVector p1, double fx1[nmsq], double fx2[nmsq], double EBEAM, TVar::VerbosityLevel verbosity){
  double xx[2]={ 0 };
  bool passPartonErgFrac=CheckPartonMomFraction(p0, p1, xx, verbosity, EBEAM);
  if (passPartonErgFrac){
    //Calculate Pdf
    //Parton Density Function is always evalualted at pT=0 frame
    //Always pass address through fortran function
    fdist_(&density_.ih1, &xx[0], &facscale_.facscale, fx1);
    fdist_(&density_.ih2, &xx[1], &facscale_.facscale, fx2);
  }
}
// SumMEPDF sums over all production parton flavors according to PDF and calls ComputePDF
double TUtil::SumMEPDF(const TLorentzVector p0, const TLorentzVector p1, double msq[nmsq][nmsq], MelaIO* RcdME, double EBEAM, TVar::VerbosityLevel verbosity){
  double fx1[nmsq]={ 0 };
  double fx2[nmsq]={ 0 };
  //double wgt_msq[nmsq][nmsq]={ { 0 } };

  ComputePDF(p0, p1, fx1, fx2, EBEAM, verbosity);
  RcdME->setPartonWeights(fx1, fx2);
  RcdME->setMEArray(msq,true);
  RcdME->computeWeightedMEArray();
  //RcdME->getWeightedMEArray(wgt_msq);
  return RcdME->getSumME();
}

// GetBoostedParticleVectors decomposes the MELACandidate object melaCand into mothers, daughters, associateds and tops
// and boosts them to the pT=0 frame for the particular mela_event.AssociationCode, which is a product of TVar::kUseAssociated* (prime numbers).
// If no mothers are present, it assigns mothers as Z>0, Z<0. If they are present, it orders them as "incoming" q-qbar / g-qbar / q-g / g-g
// Associated particles passed are different based on the code. If code==1, no associated particles are passed.
// mela_event.intermediateVids are needed to keep track of the decay mode. TVar::Process or TVar::Production do not keep track of the V/f decay modes.
// This is the major replacement functionality of the TVar lepton flavors.
void TUtil::GetBoostedParticleVectors(
  MELACandidate* melaCand,
  simple_event_record& mela_event
  ){
  // This is the beginning of one long function.

  int code = mela_event.AssociationCode;
  int aVhypo = mela_event.AssociationVCompatibility;

  pair<vector<int>, vector<TLorentzVector>> daughters;
  vector<int> idVstar;
  if (melaCand->getNDaughters()==0){
    // Undecayed Higgs has V1=H, V2=empty, no sortedDaughters!
    daughters.push_back(
      pair<vector<int>, vector<TLorentzVector>>(melaCand->id, melaCand->p4)
      );
    idVstar.push_back(melaCand->id);
  }
  else{
    // H->ffb has V1=f->f, V2=fb->fb
    // H->GG has V1=G->G, V2=G->G
    // H->ZG has V1=Z->ffb, V2=G->G
    // Everything else is as expected.
    for (int iv=0; iv<2; iv++){ // 2 Vs are guaranteed in MELACandidate.
      MELAParticle* Vdau = melaCand->getSortedV(iv);
      if (Vdau!=0){
        int idtmp = Vdau->id;
        for (int ivd=0; ivd<Vdau->getNDaughters(); ivd++){
          MELAParticle* Vdau_i = Vdau->getDaughter(ivd);
          if (Vdau_i!=0 && Vdau_i->passSelection) daughters.push_back(SimpleParticle_t(Vdau_i->id, Vdau_i->p4));
        }
        if (idtmp!=0 || Vdau->getNDaughters()>0){ // Avoid "empty" intermediate Vs of the MELACandidate object
          if (Vdau->getNDaughters()>=2 && PDGHelpers::isAPhoton(idtmp)) idtmp=23; // Special case to avoid V->2f with HVVmass==Zeromass setting (could happen by mistake)
          idVstar.push_back(idtmp);
        }
      }
    }
  }

  /***** ASSOCIATED PARTICLES *****/
  int nsatisfied_jets=0;
  int nsatisfied_lnus=0;
  int nsatisfied_gammas=0;
  vector<MELAParticle*> candidateVs; // Used if aVhypo!=0
  SimpleParticle_t associated;
  if (aVhypo!=0){

    for (int iv=2; iv<melaCand->getNSortedVs(); iv++){ // Loop over associated Vs
      MELAParticle* Vdau = melaCand->getSortedV(iv);
      if (Vdau!=0){
        bool doAdd=false;
        int idV = Vdau->id;
        if ((abs(idV)==aVhypo || idV==0) && Vdau->getNDaughters()>0){ // If the V is unknown or compatible with the requested hypothesis
          doAdd=true;
          for (int ivd=0; ivd<Vdau->getNDaughters(); ivd++){ // Loop over the daughters of V
            MELAParticle* Vdau_i = Vdau->getDaughter(ivd);
            if (Vdau_i==0){ doAdd=false; break; }
            else if (
              (mela_event.nRequested_AssociatedLeptons==0 && (PDGHelpers::isALepton(Vdau_i->id) || PDGHelpers::isANeutrino(Vdau_i->id)))
              ||
              (mela_event.nRequested_AssociatedJets==0 && PDGHelpers::isAJet(Vdau_i->id))
              ){
              doAdd=false; break;
            }
          }
        }
        if (doAdd) candidateVs.push_back(Vdau);
      }
    } // End loop over associated Vs

    cout << "TUtil::GetBoostedParticleVectors: candidateVs size = " << candidateVs.size() << endl;

    // Pick however many candidates necessary to fill up the requested number of jets or lepton(+)neutrinos
    for (unsigned int iv=0; iv<candidateVs.size(); iv++){
      MELAParticle* Vdau = candidateVs.at(iv);
      for (int ivd=0; ivd<Vdau->getNDaughters(); ivd++){ // Loop over the daughters of V
        MELAParticle* part = Vdau->getDaughter(ivd);
        if (
          part->passSelection
          &&
          (PDGHelpers::isALepton(part->id) || PDGHelpers::isANeutrino(part->id))
          &&
          nsatisfied_lnus<mela_event.nRequested_AssociatedLeptons
          ){
          nsatisfied_lnus++;
          associated.push_back(SimpleParticle_t(part->id, part->p4));
        }
        else if (
          part->passSelection
          &&
          PDGHelpers::isAJet(part->id)
          &&
          nsatisfied_jets<mela_event.nRequested_AssociatedJets
          ){
          nsatisfied_jets++;
          associated.push_back(SimpleParticle_t(part->id, part->p4));
        }
      }
    }

  }
  else{ // Could be split to aVhypo==0 and aVhypo<0 if associated V+jets is needed

    if (code%TVar::kUseAssociated_Leptons==0){
      for (int ip=0; ip<melaCand->getNAssociatedLeptons(); ip++){
        MELAParticle* part = melaCand->getAssociatedLepton(ip);
        if (part!=0 && part->passSelection && nsatisfied_lnus<mela_event.nRequested_AssociatedLeptons){
          nsatisfied_lnus++;
          associated.push_back(SimpleParticle_t(part->id, part->p4));
        }
      }
    }
    if (code%TVar::kUseAssociated_Jets==0){
      for (int ip=0; ip<melaCand->getNAssociatedJets(); ip++){
        MELAParticle* part = melaCand->getAssociatedJet(ip);
        if (part!=0 && part->passSelection && nsatisfied_jets<mela_event.nRequested_AssociatedJets){
          nsatisfied_jets++;
          associated.push_back(SimpleParticle_t(part->id, part->p4));
        }
      }
    }

  } // End if(aVhypo!=0)-else statement

  if (code%TVar::kUseAssociated_Photons==0){
    for (int ip=0; ip<melaCand->getNAssociatedPhotons(); ip++){
      MELAParticle* part = melaCand->getAssociatedPhoton(ip);
      if (part!=0 && part->passSelection && nsatisfied_gammas<mela_event.nRequested_AssociatedPhotons){
        nsatisfied_gammas++;
        associated.push_back(SimpleParticle_t(part->id, part->p4));
      }
    }
  }
  /***** END ASSOCIATED PARTICLES *****/


  /***** ASSOCIATED TOP OBJECTS *****/
  int nsatisfied_tops=0;
  int nsatisfied_antitops=0;
  vector<SimpleParticleCollection_t> topDaughters;
  vector<SimpleParticleCollection_t> antitopDaughters;
  SimpleParticleCollection_t stableTops;
  SimpleParticleCollection_t stableAntiTops;

  vector<MELATopCandidate*> tops;
  vector<MELATopCandidate*> topbars;
  vector<MELATopCandidate*> unknowntops;
  if (code%TVar::kUseAssociated_StableTops==0 && code%TVar::kUseAssociated_UnstableTops==0){ cerr << "TUtil::GetBoostedParticleVectors: Stable and unstable tops ar not supported at the same time!"  << endl; }
  else if (code%TVar::kUseAssociated_StableTops==0 || code%TVar::kUseAssociated_UnstableTops==0){

    for (int itop=0; itop<melaCand->getNAssociatedTops(); itop++){
      MELATopCandidate* theTop = melaCand->getAssociatedTop(ip);
      if (theTop!=0 && theTop->passSelection){
        vector<MELATopCandidate*>* particleArray;
        if (theTop->id==6) particleArray = &tops;
        else if (theTop->id==-6) particleArray = &topbars;
        else particleArray = &unknowntops;
        if (
          (code%TVar::kUseAssociated_StableTops==0)
          ||
          (theTop->getNDaughters()==3 && code%TVar::kUseAssociated_UnstableTops==0)
          ) particleArray->push_back((MELATopCandidate*)theTop);
      }
    }

    // Fill the stable/unstable top arrays
    for (unsigned int itop=0; itop<tops.size(); itop++){
      MELATopCandidate* theTop = tops.at(itop); // Checks are already performed
      if (code%TVar::kUseAssociated_StableTops==0 && nsatisfied_tops<mela_event.nRequested_Tops){ // Case with no daughters needed
        nsatisfied_tops++;
        stableTops.push_back(SimpleParticle_t(theTop->id, theTop->p4));
      }
      else if (code%TVar::kUseAssociated_UnstableTops==0 && theTop->getNDaughters()==3 && nsatisfied_tops<mela_event.nRequested_Tops){ // Case with daughters needed
        nsatisfied_tops++;
        SimpleParticleCollection_t vdaughters;

        MELAParticle* bottom = theTop->getLightQuark();
        MELAParticle* Wf = theTop->getWFermion();
        MELAParticle* Wfb = theTop->getWAntifermion();
        if (bottom!=0) vdaughters.push_back(SimpleParticle_t(bottom->id, bottom->p4));
        if (Wf!=0) vdaughters.push_back(SimpleParticle_t(Wf->id, Wf->p4));
        if (Wfb!=0) vdaughters.push_back(SimpleParticle_t(Wfb->id, Wfb->p4));

        if (vdaughters.size()==3) topDaughters.push_back(vdaughters);
      }
    }
    for (unsigned int itop=0; itop<topbars.size(); itop++){
      MELATopCandidate* theTop = topbars.at(itop);
      if (code%TVar::kUseAssociated_StableTops==0 && nsatisfied_antitops<mela_event.nRequested_Antitops){ // Case with no daughters needed
        nsatisfied_antitops++;
        stableAntiTops.push_back(SimpleParticle_t(theTop->id, theTop->p4));
      }
      else if (code%TVar::kUseAssociated_UnstableTops==0 && nsatisfied_antitops<mela_event.nRequested_Antitops){ // Case with daughters needed
        nsatisfied_antitops++;
        SimpleParticleCollection_t vdaughters;

        MELAParticle* bottom = theTop->getLightQuark();
        MELAParticle* Wf = theTop->getWFermion();
        MELAParticle* Wfb = theTop->getWAntifermion();
        if (bottom!=0) vdaughters.push_back(SimpleParticle_t(bottom->id, bottom->p4));
        if (Wf!=0) vdaughters.push_back(SimpleParticle_t(Wf->id, Wf->p4));
        if (Wfb!=0) vdaughters.push_back(SimpleParticle_t(Wfb->id, Wfb->p4));

        if (vdaughters.size()==3) antitopDaughters.push_back(vdaughters);
      }
      else break;
    }
    if (unknowntops.size()!=0){ // Loops over te unknown-id tops
      // Fill tops, then antitops from the unknown tops
      for (unsigned int itop=0; itop<unknowntops.size(); itop++){
        MELATopCandidate* theTop = unknowntops.at(itop);
        // t, then tb cases with no daughters needed
        if (code%TVar::kUseAssociated_StableTops==0 && nsatisfied_tops<mela_event.nRequested_Tops){
          nsatisfied_tops++;
          stableTops.push_back(SimpleParticle_t(theTop->id, theTop->p4));
        }
        else if (code%TVar::kUseAssociated_StableTops==0 && nsatisfied_antitops<mela_event.nRequested_Antitops){
          nsatisfied_antitops++;
          stableAntiTops.push_back(SimpleParticle_t(theTop->id, theTop->p4));
        }
        // t, then tb cases with daughters needed
        else if (code%TVar::kUseAssociated_UnstableTops==0 && nsatisfied_tops<mela_event.nRequested_Tops){
          nsatisfied_tops++;
          SimpleParticleCollection_t vdaughters;

          MELAParticle* bottom = theTop->getLightQuark();
          MELAParticle* Wf = theTop->getWFermion();
          MELAParticle* Wfb = theTop->getWAntifermion();
          if (bottom!=0) vdaughters.push_back(SimpleParticle_t(bottom->id, bottom->p4));
          if (Wf!=0) vdaughters.push_back(SimpleParticle_t(Wf->id, Wf->p4));
          if (Wfb!=0) vdaughters.push_back(SimpleParticle_t(Wfb->id, Wfb->p4));

          if (vdaughters.size()==3) topDaughters.push_back(vdaughters);
        }
        else if (code%TVar::kUseAssociated_UnstableTops==0 && nsatisfied_antitops<mela_event.nRequested_Antitops){
          nsatisfied_antitops++;
          SimpleParticleCollection_t vdaughters;

          MELAParticle* bottom = theTop->getLightQuark();
          MELAParticle* Wf = theTop->getWFermion();
          MELAParticle* Wfb = theTop->getWAntifermion();
          if (bottom!=0) vdaughters.push_back(SimpleParticle_t(bottom->id, bottom->p4));
          if (Wf!=0) vdaughters.push_back(SimpleParticle_t(Wf->id, Wf->p4));
          if (Wfb!=0) vdaughters.push_back(SimpleParticle_t(Wfb->id, Wfb->p4));

          if (vdaughters.size()==3) antitopDaughters.push_back(vdaughters);
        }
      }
    }

  }
  /***** END ASSOCIATED TOP OBJECTS *****/


  /***** BOOSTS TO THE CORRECT PT=0 FRAME *****/
  // Gather all final state particles collected for this frame
  TLorentzVector pTotal(0, 0, 0, 0);
  for (unsigned int ip=0; ip<daughters.size(); ip++) pTotal = pTotal + daughters.at(ip).second;
  for (unsigned int ip=0; ip<associated.size(); ip++) pTotal = pTotal + associated.at(ip).second;
  for (unsigned int ip=0; ip<stableTops.size(); ip++) pTotal = pTotal + stableTops.at(ip).second;
  for (unsigned int ip=0; ip<stableAntiTops.size(); ip++) pTotal = pTotal + stableAntiTops.at(ip).second;
  for (unsigned int ip=0; ip<topDaughters.size(); ip++){ for (unsigned int ipd=0; ipd<topDaughters.at(ip).size(); ipd++) pTotal = pTotal + topDaughters.at(ip).at(ipd).second; }
  for (unsigned int ip=0; ip<antitopDaughters.size(); ip++){ for (unsigned int ipd=0; ipd<antitopDaughters.at(ip).size(); ipd++) pTotal = pTotal + antitopDaughters.at(ip).at(ipd).second; }

  // Get the boost vector and boost all final state particles
  double qX = pTotal.X();
  double qY = pTotal.Y();
  double qE = pTotal.T();;
  if ((qX*qX+qY*qY)>0.){
    TVector3 boostV(-qX/qE, -qY/qE, 0.);
    for (unsigned int ip=0; ip<daughters.size(); ip++) daughters.at(ip).second.Boost(boostV);
    for (unsigned int ip=0; ip<associated.size(); ip++) associated.at(ip).second.Boost(boostV);
    for (unsigned int ip=0; ip<stableTops.size(); ip++) stableTops.at(ip).second.Boost(boostV);
    for (unsigned int ip=0; ip<stableAntiTops.size(); ip++) stableAntiTops.at(ip).second.Boost(boostV);
    for (unsigned int ip=0; ip<topDaughters.size(); ip++){ for (unsigned int ipd=0; ipd<topDaughters.at(ip).size(); ipd++) topDaughters.at(ip).at(ipd).second.Boost(boostV); }
    for (unsigned int ip=0; ip<antitopDaughters.size(); ip++){ for (unsigned int ipd=0; ipd<antitopDaughters.at(ip).size(); ipd++) antitopDaughters.at(ip).at(ipd).second.Boost(boostV); }
    pTotal.Boost(boostV);
  }

  // Mothers need special treatment:
  // In case they are undefined, mother id is unknown, and mothers are ordered as M1==(pz>0) and M2==(pz<0).
  // In case they are defined, mothers are first matched to the assumed momenta through their pz's.
  // They are ordered by M1==f, M2==fb afterward if this is possible.
  // Notice that the momenta of the mother objects are not actually used. If the event is truly a gen. particle, everything would be in the pT=0 frame to begin with, and everybody is happy in reweighting.
  double sysPz= pTotal.Z();
  double sysE = pTotal.T();
  double pz0 = (sysE+sysPz)/2.;
  double pz1 = -(sysE-sysPz)/2.;
  int motherId[2]={ 0, 0 };
  if (melaCand->getNMothers()==2){
    // Match the assumed M1==(pz>0) and M2==(pz<0) to the pz of the actual mothers.
    if ((melaCand->getMother(0)->p4).Z()<0.){ // Swap by Z to get the correct momentum matching default M1, M2
      double pztmp = pz1;
      pz1=pz0;
      pz0=pztmp;
    }
    if (motherId[0]<0){ // Swap the ids of mothers and their "assumed" momenta to achieve ordering as "incoming" q-qbar
      double pztmp = pz1;
      pz1=pz0;
      pz0=pztmp;
      for (int ip=0; ip<2; ip++) motherId[1-ip]=melaCand->getMother(ip)->id;
    }
    else{ // No swap of mothers necessary, just get their ids
      for (int ip=0; ip<2; ip++) motherId[ip]=melaCand->getMother(ip)->id;
    }
  }
  TLorentzVector pM[2];
  pM[0].SetXYZT(0., 0., pz0, fabs(pz0));
  pM[1].SetXYZT(0., 0., pz1, fabs(pz1));

  // Fill the ids of the V intermediates to the candidate daughters
  mela_event.intermediateVid.clear();
  for (unsigned int ip=0; ip<idVstar.size(); ip++) mela_event.intermediateVid.push_back(idVstar.at(ip));
  // Fill the mothers
  mela_event.pMothers.clear();
  for (unsigned int ip=0; ip<2; ip++){ mela_event.pMothers.push_back(SimpleParticle_t(motherId[ip], pM[ip])); }
  // Fill the daughters
  mela_event.pDaughters.clear();
  for (unsigned int ip=0; ip<daughters.size(); ip++) mela_event.pDaughters.push_back(daughters.at(ip));
  // Fill the associated particles
  mela_event.pAssociated.clear();
  for (unsigned int ip=0; ip<associated.size(); ip++) mela_event.pAssociated.push_back(associated.at(ip));
  // Fill the stable tops and antitops
  mela_event.pStableTops.clear();
  for (unsigned int ip=0; ip<stableTops.size(); ip++) mela_event.pStableTops.push_back(stableTops.at(ip));
  mela_event.pStableAntitops.clear();
  for (unsigned int ip=0; ip<stableAntitops.size(); ip++) mela_event.pStableAntitops.push_back(stableAntitops.at(ip));
  // Fill the daughters of unstable tops and antitops
  mela_event.pTopDaughters.clear();
  for (unsigned int ip=0; ip<topDaughters.size(); ip++) mela_event.pTopDaughters.push_back(topDaughters.at(ip));
  mela_event.pAntitopDaughters.clear();
  for (unsigned int ip=0; ip<antitopDaughters.size(); ip++) mela_event.pAntitopDaughters.push_back(antitopDaughters.at(ip));

  // This is the end of one long function.
}

// Convert vectors of simple particles to MELAParticles and create a MELACandidate
// The output lists could be members of TEvtProb directly.
MELACandidate* TUtil::ConvertVectorFormat(
  // Inputs
  SimpleParticleCollection_t* pDaughters,
  SimpleParticleCollection_t* pAssociated,
  SimpleParticleCollection_t* pMothers,
  bool isGen,
  // Outputs
  std::vector<MELAParticle*>* particleList,
  std::vector<MELACandidate*>* candList
  ){
  MELACandidate* cand=0;

  if (pDaughters==0){ cerr << "TUtil::ConvertVectorFormat: No daughters!" << endl; return cand; }
  else if (pDaughters->size()==0){ cerr << "TUtil::ConvertVectorFormat: Daughter size==0!" << endl; return cand; }
  else if (pDaughters->size()>4){ cerr << "TUtil::ConvertVectorFormat: Daughter size " << pDaughters->size() << ">4 is not supported!" << endl; return cand; }
  if (pMothers!=0 && pMothers->size()!=2){ cerr << "TUtil::ConvertVectorFormat: Mothers momentum size (" << pMothers->size() << ") has to have had been 2! Continuing by omitting mothers." << endl; /*return cand;*/ }

  std::vector<MELAParticle>* daughters;
  std::vector<MELAParticle>* aparticles;
  std::vector<MELAParticle>* mothers;

  for (unsigned int ip=0; ip<pDaughters->size(); ip++){
    MELAParticle* onePart = new MELAParticle((pDaughters->at(ip)).first, (pDaughters->at(ip)).second);
    onePart->setGenStatus(1); // Final state status
    particleList->push_back(onePart);
    daughters.push_back(onePart);
  }
  if (pAssociated!=0){
    for (unsigned int ip=0; ip<pAssociated->size(); ip++){
      MELAParticle* onePart = new MELAParticle((pAssociated->at(ip)).first, (pAssociated->at(ip)).second);
      onePart->setGenStatus(1); // Final state status
      particleList->push_back(onePart);
      aparticles.push_back(onePart);
    }
  }
  if (pMothers!=0 && pMothers->size()==2){
    for (unsigned int ip=0; ip<pMothers->size(); ip++){
      MELAParticle* onePart = new MELAParticle((pMothers->at(ip)).first, (pMothers->at(ip)).second);
      onePart->setGenStatus(-1); // Mother status
      particleList->push_back(onePart);
      mothers.push_back(onePart);
    }
  }

  /***** Adaptation of LHEAnalyzer::Event::constructVVCandidates *****/
  /*
  The assumption is that the daughters make sense for either ffb, gamgam, Zgam, ZZ or WW.
  No checking is done on whether particle-antiparticle pairing is correct when necessary.
  If not, you will get a seg. fault!
  */
  // Undecayed Higgs
  if (daughters.size()==1) cand = new MELACandidate(25, (daughters.at(0))->p4); // No sorting!
  // GG / ff final states
  else if (daughters.size()==2){
    MELAParticle* F1 = daughters.at(0);
    MELAParticle* F2 = daughters.at(1);
    TLorentzVector pH = F1->p4+F2->p4;
    cand = new MELACandidate(25, pH);
    cand->addDaughter(F1);
    cand->addDaughter(F2);
    double defaultHVVmass = PDGHelpers::HVVmass;
    PDGHelpers::setHVVmass(Zeromass);
    cand->sortDaughters();
    PDGHelpers::setHVVmass(defaultHVVmass);
  }
  // ZG / WG
  else if (daughters.size()==3){
    MELAParticle* F1 = daughters.at(0);
    MELAParticle* F2 = daughters.at(1);
    MELAParticle* gamma = daughters.at(2);
    if (PDGHelpers::isAPhoton(F1->id)){
      MELAParticle* tmp = F1;
      F1 = gamma;
      gamma = tmp;
    }
    else if (PDGHelpers::isAPhoton(F2->id)){
      MELAParticle* tmp = F2;
      F2 = gamma;
      gamma = tmp;
    }
    TLorentzVector pH = F1->p4+F2->p4+gamma->p4;
    double charge = F1->charge()+F2->charge()+gamma->charge();
    cand = new MELACandidate(25, pH);
    cand->addDaughter(F1);
    cand->addDaughter(F2);
    cand->addDaughter(gamma);
    double defaultHVVmass = PDGHelpers::HVVmass;
    if (fabs(charge)<0.01) PDGHelpers::setHVVmass(PDGHelpers::Zmass); // ZG
    else PDGHelpers::setHVVmass(PDGHelpers::Wmass); // WG,GW (?), un-tested
    cand->sortDaughters();
    PDGHelpers::setHVVmass(defaultHVVmass);
  }
  // ZZ / WW / ZW
  else/* if (daughters.size()==4)*/{
    TLorentzVector pH(0, 0, 0, 0);
    double charge = 0.;
    for (int ip=0; ip<4; ip++){ pH = pH + (daughters.at(ip))->p4; charge += (daughters.at(ip))->charge(); }
    cand = new MELACandidate(25, pH);
    for (int ip=0; ip<4; ip++) cand->addDaughter(daughters.at(ip));
    double defaultHVVmass = PDGHelpers::HVVmass;
    if (fabs(charge)>0.01) PDGHelpers::setHVVmass(PDGHelpers::Wmass); // WZ,ZW (?), un-tested
    cand->sortDaughters();
    PDGHelpers::setHVVmass(defaultHVVmass);
  }
  /***** Adaptation of LHEAnalyzer::Event::addVVCandidateMother *****/
  if (mothers.size()>0){ // ==2
    for (int ip=0; ip<mothers.size(); ip++) cand->addMother(mothers.at(ip));
    if (isGen) cand->setGenStatus(-1); // Candidate is a gen. particle!
  }
  /***** Adaptation of LHEAnalyzer::Event::addVVCandidateAppendages *****/
  if (aparticles.size()>0){ // ==2
    for (int ip=0; ip<aparticles.size(); ip++){
      const int partId = (aparticles.at(ip))->id;
      if (PDGHelpers::isALepton(partId)) cand->addAssociatedLeptons(aparticles.at(ip));
      else if (PDGHelpers::isANeutrino(partId)) cand->addAssociatedNeutrinos(aparticles.at(ip)); // Be careful: Neutrinos are neutrinos, but also "leptons" in MELACandidate!
      else if (PDGHelpers::isAPhoton(partId)) cand->addAssociatedPhotons(aparticles.at(ip));
      else if (PDGHelpers::isAJet(partId)) cand->addAssociatedJets(aparticles.at(ip));
    }
    cand->addAssociatedVs(); // For the VH topology
  }
  candList->push_back(cand);
  return cand;
}

// Convert the vector of top daughters (as simple particles) to MELAParticles and create a MELATopCandidate
// The output lists could be members of TEvtProb directly.
MELATopCandidate* TUtil::ConvertTopCandidate(
  // Input
  SimpleParticleCollection_t* TopDaughters,
  // Outputs
  std::vector<MELAParticle*>* particleList,
  std::vector<MELATopCandidate*>* topCandList
  ){
  MELATopCandidate* cand=0;

  if (TopDaughters==0){ cerr << "TUtil::ConvertTopCandidate: No daughters!" << endl; return cand; }
  else if (TopDaughters->size()==0){ cerr << "TUtil::ConvertTopCandidate: Daughter size==0!" << endl; return cand; }
  else if (!(TopDaughters->size()==1 || TopDaughters->size()==3)){ cerr << "TUtil::ConvertVectorFormat: Daughter size " << TopDaughters->size() << "!=1 or 3 is not supported!" << endl; return cand; }

  if (TopDaughters->size()==1){
    if (abs((TopDaughters->at(0)).first)==6 || (TopDaughters->at(0)).first==0){
      cand = new MELATopCandidate((TopDaughters->at(0)).first, (TopDaughters->at(0)).second);
      topCandList->push_back(cand);
    }
  }
  else if (TopDaughters->size()==3){
    MELAParticle* bottom = new MELAParticle((TopDaughters->at(0)).first, (TopDaughters->at(0)).second);
    MELAParticle* Wf = new MELAParticle((TopDaughters->at(0)).first, (TopDaughters->at(0)).second);
    MELAParticle* Wfb = new MELAParticle((TopDaughters->at(0)).first, (TopDaughters->at(0)).second);

    if (Wf->id<0 || Wfb->id>0){
      MELAParticle* parttmp = Wf;
      Wf=Wfb;
      Wfb=parttmp;
    }

    particleList->push_back(bottom);
    particleList->push_back(Wf);
    particleList->push_back(Wfb);

    cand = new MELATopCandidate(bottom, Wf, Wfb);
    topCandList->push_back(cand);
  }
  return cand;
}


