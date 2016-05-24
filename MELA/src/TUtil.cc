#include "ZZMatrixElement/MELA/interface/TUtil.hh"
#include <iostream>
#include <cstdio>
#include <cmath>
#include <TMath.h>

using namespace std;

void SetEwkCouplingParameters(){

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


double InterpretScaleScheme(TVar::Production production, TVar::MatrixElement matrixElement, TVar::EventScaleScheme scheme, TLorentzVector p[mxpart]){
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


void SetAlphaS(double Q_ren, double Q_fac, double multiplier_ren, double multiplier_fac, int mynloop, int mynflav, string mypartons){
  bool hasReset=false;
  if (multiplier_ren<=0 || multiplier_fac<=0){
    cout << "Invalid scale multipliers" << endl;
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
    cout << "Only default :: cteq6l1 or cteq6_l are supported. Modify mela.cc symlinks, put the pdf table into data/Pdfdata and retry. Setting to Default..." << endl;
    mypartons = "Default";
  }
  // From pdfwrapper_linux.f:
  if (mypartons.compare("cteq6_l")==0) couple_.amz = 0.118;
	else if(mypartons.compare("cteq6l1")==0 || mypartons.compare("Default")==0) couple_.amz = 0.130;
	else couple_.amz = 0.118; // Add pdf as appropriate

// For proper pdfwrapper_linux.f execution (alpha_s computation does not use pdf but examines the pdf name to initialize amz.)
	if(!nflav_is_same){
		nflav_.nflav = mynflav;

		if(mypartons.compare("Default")!=0) sprintf(pdlabel_.pdlabel,"%s",mypartons.c_str());
		else sprintf(pdlabel_.pdlabel,"%s","cteq6l1"); // Default pdf is cteq6l1
		coupling2_();
	}
	else{
		qcdcouple_.as = alphas_(&(scale_.scale),&(couple_.amz),&(nlooprun_.nlooprun));
	}

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

void InitJHUGenMELA(const char* pathtoPDFSet, int PDFMember){
  char path_pdf_c[200];
  sprintf(path_pdf_c, "%s", pathtoPDFSet);
  int pathpdfLength = strlen(path_pdf_c);
  __modjhugen_MOD_initfirsttime(path_pdf_c, &pathpdfLength, &PDFMember);
}
void SetJHUGenHiggsMassWidth(double MReso, double GaReso){
  MReso /= 100.; // GeV units in JHUGen
  GaReso /= 100.; // GeV units in JHUGen
  __modjhugenmela_MOD_sethiggsmasswidth(&MReso, &GaReso);
}
void SetJHUGenDistinguishWWCouplings(bool doAllow){
  int iAllow = (doAllow ? 1 : 0);
  __modjhugenmela_MOD_setdistinguishwwcouplingsflag(&iAllow);
}
void SetMCFMSpinZeroVVCouplings(bool useBSM, double Hvvcoupl[SIZE_HVV][2], double Hwwcoupl[SIZE_HVV][2]){
  if (!useBSM){
    spinzerohiggs_anomcoupl_.AllowAnomalousCouplings = false;
    spinzerohiggs_anomcoupl_.ghz1[0] =  1;
    spinzerohiggs_anomcoupl_.ghz2[0] =  0;
    spinzerohiggs_anomcoupl_.ghz3[0] =  0;
    spinzerohiggs_anomcoupl_.ghz4[0] =  0;
    /*
    spinzerohiggs_anomcoupl_.ghzgs2[0] = 0;
    spinzerohiggs_anomcoupl_.ghzgs3[0] = 0;
    spinzerohiggs_anomcoupl_.ghzgs4[0] = 0;
    spinzerohiggs_anomcoupl_.ghgsgs2[0] = 0;
    spinzerohiggs_anomcoupl_.ghgsgs3[0] = 0;
    spinzerohiggs_anomcoupl_.ghgsgs4[0] = 0;
    */
    spinzerohiggs_anomcoupl_.ghz1_prime[0] = 0;
    spinzerohiggs_anomcoupl_.ghz1_prime2[0] = 0;
    spinzerohiggs_anomcoupl_.ghz1_prime3[0] = 0;
    spinzerohiggs_anomcoupl_.ghz1_prime4[0] = 0;
    spinzerohiggs_anomcoupl_.ghz1_prime5[0] = 0;
    spinzerohiggs_anomcoupl_.ghz1_prime6[0] = 0;
    spinzerohiggs_anomcoupl_.ghz1_prime7[0] = 0;
    spinzerohiggs_anomcoupl_.ghz2_prime[0] = 0;
    spinzerohiggs_anomcoupl_.ghz2_prime2[0] = 0;
    spinzerohiggs_anomcoupl_.ghz2_prime3[0] = 0;
    spinzerohiggs_anomcoupl_.ghz2_prime4[0] = 0;
    spinzerohiggs_anomcoupl_.ghz2_prime5[0] = 0;
    spinzerohiggs_anomcoupl_.ghz2_prime6[0] = 0;
    spinzerohiggs_anomcoupl_.ghz2_prime7[0] = 0;
    spinzerohiggs_anomcoupl_.ghz3_prime[0] = 0;
    spinzerohiggs_anomcoupl_.ghz3_prime2[0] = 0;
    spinzerohiggs_anomcoupl_.ghz3_prime3[0] = 0;
    spinzerohiggs_anomcoupl_.ghz3_prime4[0] = 0;
    spinzerohiggs_anomcoupl_.ghz3_prime5[0] = 0;
    spinzerohiggs_anomcoupl_.ghz3_prime6[0] = 0;
    spinzerohiggs_anomcoupl_.ghz3_prime7[0] = 0;
    spinzerohiggs_anomcoupl_.ghz4_prime[0] = 0;
    spinzerohiggs_anomcoupl_.ghz4_prime2[0] = 0;
    spinzerohiggs_anomcoupl_.ghz4_prime3[0] = 0;
    spinzerohiggs_anomcoupl_.ghz4_prime4[0] = 0;
    spinzerohiggs_anomcoupl_.ghz4_prime5[0] = 0;
    spinzerohiggs_anomcoupl_.ghz4_prime6[0] = 0;
    spinzerohiggs_anomcoupl_.ghz4_prime7[0] = 0;
    //spinzerohiggs_anomcoupl_.ghzgs1_prime2[0] = 0;

    spinzerohiggs_anomcoupl_.ghz1[1] =  0;
    spinzerohiggs_anomcoupl_.ghz2[1] =  0;
    spinzerohiggs_anomcoupl_.ghz3[1] =  0;
    spinzerohiggs_anomcoupl_.ghz4[1] =  0;
    /*
    spinzerohiggs_anomcoupl_.ghzgs2[1] = 0;
    spinzerohiggs_anomcoupl_.ghzgs3[1] = 0;
    spinzerohiggs_anomcoupl_.ghzgs4[1] = 0;
    spinzerohiggs_anomcoupl_.ghgsgs2[1] = 0;
    spinzerohiggs_anomcoupl_.ghgsgs3[1] = 0;
    spinzerohiggs_anomcoupl_.ghgsgs4[1] = 0;
    */
    spinzerohiggs_anomcoupl_.ghz1_prime[1] = 0;
    spinzerohiggs_anomcoupl_.ghz1_prime2[1] = 0;
    spinzerohiggs_anomcoupl_.ghz1_prime3[1] = 0;
    spinzerohiggs_anomcoupl_.ghz1_prime4[1] = 0;
    spinzerohiggs_anomcoupl_.ghz1_prime5[1] = 0;
    spinzerohiggs_anomcoupl_.ghz1_prime6[1] = 0;
    spinzerohiggs_anomcoupl_.ghz1_prime7[1] = 0;
    spinzerohiggs_anomcoupl_.ghz2_prime[1] = 0;
    spinzerohiggs_anomcoupl_.ghz2_prime2[1] = 0;
    spinzerohiggs_anomcoupl_.ghz2_prime3[1] = 0;
    spinzerohiggs_anomcoupl_.ghz2_prime4[1] = 0;
    spinzerohiggs_anomcoupl_.ghz2_prime5[1] = 0;
    spinzerohiggs_anomcoupl_.ghz2_prime6[1] = 0;
    spinzerohiggs_anomcoupl_.ghz2_prime7[1] = 0;
    spinzerohiggs_anomcoupl_.ghz3_prime[1] = 0;
    spinzerohiggs_anomcoupl_.ghz3_prime2[1] = 0;
    spinzerohiggs_anomcoupl_.ghz3_prime3[1] = 0;
    spinzerohiggs_anomcoupl_.ghz3_prime4[1] = 0;
    spinzerohiggs_anomcoupl_.ghz3_prime5[1] = 0;
    spinzerohiggs_anomcoupl_.ghz3_prime6[1] = 0;
    spinzerohiggs_anomcoupl_.ghz3_prime7[1] = 0;
    spinzerohiggs_anomcoupl_.ghz4_prime[1] = 0;
    spinzerohiggs_anomcoupl_.ghz4_prime2[1] = 0;
    spinzerohiggs_anomcoupl_.ghz4_prime3[1] = 0;
    spinzerohiggs_anomcoupl_.ghz4_prime4[1] = 0;
    spinzerohiggs_anomcoupl_.ghz4_prime5[1] = 0;
    spinzerohiggs_anomcoupl_.ghz4_prime6[1] = 0;
    spinzerohiggs_anomcoupl_.ghz4_prime7[1] = 0;
    //spinzerohiggs_anomcoupl_.ghzgs1_prime2[1] = 0;
    //
    spinzerohiggs_anomcoupl_.ghw1[0] =  1;
    spinzerohiggs_anomcoupl_.ghw2[0] =  0;
    spinzerohiggs_anomcoupl_.ghw3[0] =  0;
    spinzerohiggs_anomcoupl_.ghw4[0] =  0;
    /*
    spinzerohiggs_anomcoupl_.ghwgs2[0] = 0;
    spinzerohiggs_anomcoupl_.ghwgs3[0] = 0;
    spinzerohiggs_anomcoupl_.ghwgs4[0] = 0;
    //spinzerohiggs_anomcoupl_.ghgsgs2[0] = 0;
    //spinzerohiggs_anomcoupl_.ghgsgs3[0] = 0;
    //spinzerohiggs_anomcoupl_.ghgsgs4[0] = 0;
    */
    spinzerohiggs_anomcoupl_.ghw1_prime[0] = 0;
    spinzerohiggs_anomcoupl_.ghw1_prime2[0] = 0;
    spinzerohiggs_anomcoupl_.ghw1_prime3[0] = 0;
    spinzerohiggs_anomcoupl_.ghw1_prime4[0] = 0;
    spinzerohiggs_anomcoupl_.ghw1_prime5[0] = 0;
    spinzerohiggs_anomcoupl_.ghw1_prime6[0] = 0;
    spinzerohiggs_anomcoupl_.ghw1_prime7[0] = 0;
    spinzerohiggs_anomcoupl_.ghw2_prime[0] = 0;
    spinzerohiggs_anomcoupl_.ghw2_prime2[0] = 0;
    spinzerohiggs_anomcoupl_.ghw2_prime3[0] = 0;
    spinzerohiggs_anomcoupl_.ghw2_prime4[0] = 0;
    spinzerohiggs_anomcoupl_.ghw2_prime5[0] = 0;
    spinzerohiggs_anomcoupl_.ghw2_prime6[0] = 0;
    spinzerohiggs_anomcoupl_.ghw2_prime7[0] = 0;
    spinzerohiggs_anomcoupl_.ghw3_prime[0] = 0;
    spinzerohiggs_anomcoupl_.ghw3_prime2[0] = 0;
    spinzerohiggs_anomcoupl_.ghw3_prime3[0] = 0;
    spinzerohiggs_anomcoupl_.ghw3_prime4[0] = 0;
    spinzerohiggs_anomcoupl_.ghw3_prime5[0] = 0;
    spinzerohiggs_anomcoupl_.ghw3_prime6[0] = 0;
    spinzerohiggs_anomcoupl_.ghw3_prime7[0] = 0;
    spinzerohiggs_anomcoupl_.ghw4_prime[0] = 0;
    spinzerohiggs_anomcoupl_.ghw4_prime2[0] = 0;
    spinzerohiggs_anomcoupl_.ghw4_prime3[0] = 0;
    spinzerohiggs_anomcoupl_.ghw4_prime4[0] = 0;
    spinzerohiggs_anomcoupl_.ghw4_prime5[0] = 0;
    spinzerohiggs_anomcoupl_.ghw4_prime6[0] = 0;
    spinzerohiggs_anomcoupl_.ghw4_prime7[0] = 0;
    //spinzerohiggs_anomcoupl_.ghwgs1_prime2[0] = 0;

    spinzerohiggs_anomcoupl_.ghw1[1] =  0;
    spinzerohiggs_anomcoupl_.ghw2[1] =  0;
    spinzerohiggs_anomcoupl_.ghw3[1] =  0;
    spinzerohiggs_anomcoupl_.ghw4[1] =  0;
    /*
    spinzerohiggs_anomcoupl_.ghwgs2[1] = 0;
    spinzerohiggs_anomcoupl_.ghwgs3[1] = 0;
    spinzerohiggs_anomcoupl_.ghwgs4[1] = 0;
    //spinzerohiggs_anomcoupl_.ghgsgs2[1] = 0;
    //spinzerohiggs_anomcoupl_.ghgsgs3[1] = 0;
    //spinzerohiggs_anomcoupl_.ghgsgs4[1] = 0;
    */
    spinzerohiggs_anomcoupl_.ghw1_prime[1] = 0;
    spinzerohiggs_anomcoupl_.ghw1_prime2[1] = 0;
    spinzerohiggs_anomcoupl_.ghw1_prime3[1] = 0;
    spinzerohiggs_anomcoupl_.ghw1_prime4[1] = 0;
    spinzerohiggs_anomcoupl_.ghw1_prime5[1] = 0;
    spinzerohiggs_anomcoupl_.ghw1_prime6[1] = 0;
    spinzerohiggs_anomcoupl_.ghw1_prime7[1] = 0;
    spinzerohiggs_anomcoupl_.ghw2_prime[1] = 0;
    spinzerohiggs_anomcoupl_.ghw2_prime2[1] = 0;
    spinzerohiggs_anomcoupl_.ghw2_prime3[1] = 0;
    spinzerohiggs_anomcoupl_.ghw2_prime4[1] = 0;
    spinzerohiggs_anomcoupl_.ghw2_prime5[1] = 0;
    spinzerohiggs_anomcoupl_.ghw2_prime6[1] = 0;
    spinzerohiggs_anomcoupl_.ghw2_prime7[1] = 0;
    spinzerohiggs_anomcoupl_.ghw3_prime[1] = 0;
    spinzerohiggs_anomcoupl_.ghw3_prime2[1] = 0;
    spinzerohiggs_anomcoupl_.ghw3_prime3[1] = 0;
    spinzerohiggs_anomcoupl_.ghw3_prime4[1] = 0;
    spinzerohiggs_anomcoupl_.ghw3_prime5[1] = 0;
    spinzerohiggs_anomcoupl_.ghw3_prime6[1] = 0;
    spinzerohiggs_anomcoupl_.ghw3_prime7[1] = 0;
    spinzerohiggs_anomcoupl_.ghw4_prime[1] = 0;
    spinzerohiggs_anomcoupl_.ghw4_prime2[1] = 0;
    spinzerohiggs_anomcoupl_.ghw4_prime3[1] = 0;
    spinzerohiggs_anomcoupl_.ghw4_prime4[1] = 0;
    spinzerohiggs_anomcoupl_.ghw4_prime5[1] = 0;
    spinzerohiggs_anomcoupl_.ghw4_prime6[1] = 0;
    spinzerohiggs_anomcoupl_.ghw4_prime7[1] = 0;
    //spinzerohiggs_anomcoupl_.ghwgs1_prime2[1] = 0;
    //
  }
  else{
    spinzerohiggs_anomcoupl_.AllowAnomalousCouplings=true;
    spinzerohiggs_anomcoupl_.ghz1[0] =  Hvvcoupl[0][0];
    spinzerohiggs_anomcoupl_.ghz2[0] =  Hvvcoupl[1][0];
    spinzerohiggs_anomcoupl_.ghz3[0] =  Hvvcoupl[2][0];
    spinzerohiggs_anomcoupl_.ghz4[0] =  Hvvcoupl[3][0];
    /*
    spinzerohiggs_anomcoupl_.ghzgs2[0] = Hvvcoupl[4][0];
    spinzerohiggs_anomcoupl_.ghzgs3[0] = Hvvcoupl[5][0];
    spinzerohiggs_anomcoupl_.ghzgs4[0] = Hvvcoupl[6][0];
    spinzerohiggs_anomcoupl_.ghgsgs2[0] = Hvvcoupl[7][0];
    spinzerohiggs_anomcoupl_.ghgsgs3[0] = Hvvcoupl[8][0];
    spinzerohiggs_anomcoupl_.ghgsgs4[0] = Hvvcoupl[9][0];
    */
    spinzerohiggs_anomcoupl_.ghz1_prime[0] = Hvvcoupl[10][0];
    spinzerohiggs_anomcoupl_.ghz1_prime2[0] = Hvvcoupl[11][0];
    spinzerohiggs_anomcoupl_.ghz1_prime3[0] = Hvvcoupl[12][0];
    spinzerohiggs_anomcoupl_.ghz1_prime4[0] = Hvvcoupl[13][0];
    spinzerohiggs_anomcoupl_.ghz1_prime5[0] = Hvvcoupl[14][0];
    spinzerohiggs_anomcoupl_.ghz2_prime[0] = Hvvcoupl[15][0];
    spinzerohiggs_anomcoupl_.ghz2_prime2[0] = Hvvcoupl[16][0];
    spinzerohiggs_anomcoupl_.ghz2_prime3[0] = Hvvcoupl[17][0];
    spinzerohiggs_anomcoupl_.ghz2_prime4[0] = Hvvcoupl[18][0];
    spinzerohiggs_anomcoupl_.ghz2_prime5[0] = Hvvcoupl[19][0];
    spinzerohiggs_anomcoupl_.ghz3_prime[0] = Hvvcoupl[20][0];
    spinzerohiggs_anomcoupl_.ghz3_prime2[0] = Hvvcoupl[21][0];
    spinzerohiggs_anomcoupl_.ghz3_prime3[0] = Hvvcoupl[22][0];
    spinzerohiggs_anomcoupl_.ghz3_prime4[0] = Hvvcoupl[23][0];
    spinzerohiggs_anomcoupl_.ghz3_prime5[0] = Hvvcoupl[24][0];
    spinzerohiggs_anomcoupl_.ghz4_prime[0] = Hvvcoupl[25][0];
    spinzerohiggs_anomcoupl_.ghz4_prime2[0] = Hvvcoupl[26][0];
    spinzerohiggs_anomcoupl_.ghz4_prime3[0] = Hvvcoupl[27][0];
    spinzerohiggs_anomcoupl_.ghz4_prime4[0] = Hvvcoupl[28][0];
    spinzerohiggs_anomcoupl_.ghz4_prime5[0] = Hvvcoupl[29][0];
    //spinzerohiggs_anomcoupl_.ghzgs1_prime2[0] = Hvvcoupl[30][0];
    spinzerohiggs_anomcoupl_.ghz1_prime6[0] = Hvvcoupl[31][0];
    spinzerohiggs_anomcoupl_.ghz1_prime7[0] = Hvvcoupl[32][0];
    spinzerohiggs_anomcoupl_.ghz2_prime6[0] = Hvvcoupl[33][0];
    spinzerohiggs_anomcoupl_.ghz2_prime7[0] = Hvvcoupl[34][0];
    spinzerohiggs_anomcoupl_.ghz3_prime6[0] = Hvvcoupl[35][0];
    spinzerohiggs_anomcoupl_.ghz3_prime7[0] = Hvvcoupl[36][0];
    spinzerohiggs_anomcoupl_.ghz4_prime6[0] = Hvvcoupl[37][0];
    spinzerohiggs_anomcoupl_.ghz4_prime7[0] = Hvvcoupl[38][0];

    spinzerohiggs_anomcoupl_.ghz1[1] =  Hvvcoupl[0][1];
    spinzerohiggs_anomcoupl_.ghz2[1] =  Hvvcoupl[1][1];
    spinzerohiggs_anomcoupl_.ghz3[1] =  Hvvcoupl[2][1];
    spinzerohiggs_anomcoupl_.ghz4[1] =  Hvvcoupl[3][1];
    /*
    spinzerohiggs_anomcoupl_.ghzgs2[1] = Hvvcoupl[4][1];
    spinzerohiggs_anomcoupl_.ghzgs3[1] = Hvvcoupl[5][1];
    spinzerohiggs_anomcoupl_.ghzgs4[1] = Hvvcoupl[6][1];
    spinzerohiggs_anomcoupl_.ghgsgs2[1] = Hvvcoupl[7][1];
    spinzerohiggs_anomcoupl_.ghgsgs3[1] = Hvvcoupl[8][1];
    spinzerohiggs_anomcoupl_.ghgsgs4[1] = Hvvcoupl[9][1];
    */
    spinzerohiggs_anomcoupl_.ghz1_prime[1] = Hvvcoupl[10][1];
    spinzerohiggs_anomcoupl_.ghz1_prime2[1] = Hvvcoupl[11][1];
    spinzerohiggs_anomcoupl_.ghz1_prime3[1] = Hvvcoupl[12][1];
    spinzerohiggs_anomcoupl_.ghz1_prime4[1] = Hvvcoupl[13][1];
    spinzerohiggs_anomcoupl_.ghz1_prime5[1] = Hvvcoupl[14][1];
    spinzerohiggs_anomcoupl_.ghz2_prime[1] = Hvvcoupl[15][1];
    spinzerohiggs_anomcoupl_.ghz2_prime2[1] = Hvvcoupl[16][1];
    spinzerohiggs_anomcoupl_.ghz2_prime3[1] = Hvvcoupl[17][1];
    spinzerohiggs_anomcoupl_.ghz2_prime4[1] = Hvvcoupl[18][1];
    spinzerohiggs_anomcoupl_.ghz2_prime5[1] = Hvvcoupl[19][1];
    spinzerohiggs_anomcoupl_.ghz3_prime[1] = Hvvcoupl[20][1];
    spinzerohiggs_anomcoupl_.ghz3_prime2[1] = Hvvcoupl[21][1];
    spinzerohiggs_anomcoupl_.ghz3_prime3[1] = Hvvcoupl[22][1];
    spinzerohiggs_anomcoupl_.ghz3_prime4[1] = Hvvcoupl[23][1];
    spinzerohiggs_anomcoupl_.ghz3_prime5[1] = Hvvcoupl[24][1];
    spinzerohiggs_anomcoupl_.ghz4_prime[1] = Hvvcoupl[25][1];
    spinzerohiggs_anomcoupl_.ghz4_prime2[1] = Hvvcoupl[26][1];
    spinzerohiggs_anomcoupl_.ghz4_prime3[1] = Hvvcoupl[27][1];
    spinzerohiggs_anomcoupl_.ghz4_prime4[1] = Hvvcoupl[28][1];
    spinzerohiggs_anomcoupl_.ghz4_prime5[1] = Hvvcoupl[29][1];
    //spinzerohiggs_anomcoupl_.ghzgs1_prime2[1] = Hvvcoupl[30][1];
    spinzerohiggs_anomcoupl_.ghz1_prime6[1] = Hvvcoupl[31][1];
    spinzerohiggs_anomcoupl_.ghz1_prime7[1] = Hvvcoupl[32][1];
    spinzerohiggs_anomcoupl_.ghz2_prime6[1] = Hvvcoupl[33][1];
    spinzerohiggs_anomcoupl_.ghz2_prime7[1] = Hvvcoupl[34][1];
    spinzerohiggs_anomcoupl_.ghz3_prime6[1] = Hvvcoupl[35][1];
    spinzerohiggs_anomcoupl_.ghz3_prime7[1] = Hvvcoupl[36][1];
    spinzerohiggs_anomcoupl_.ghz4_prime6[1] = Hvvcoupl[37][1];
    spinzerohiggs_anomcoupl_.ghz4_prime7[1] = Hvvcoupl[38][1];
    //
    spinzerohiggs_anomcoupl_.ghw1[0] =  Hwwcoupl[0][0];
    spinzerohiggs_anomcoupl_.ghw2[0] =  Hwwcoupl[1][0];
    spinzerohiggs_anomcoupl_.ghw3[0] =  Hwwcoupl[2][0];
    spinzerohiggs_anomcoupl_.ghw4[0] =  Hwwcoupl[3][0];
    /*
    spinzerohiggs_anomcoupl_.ghwgs2[0] = Hwwcoupl[4][0];
    spinzerohiggs_anomcoupl_.ghwgs3[0] = Hwwcoupl[5][0];
    spinzerohiggs_anomcoupl_.ghwgs4[0] = Hwwcoupl[6][0];
    //spinzerohiggs_anomcoupl_.ghgsgs2[0] = Hwwcoupl[7][0];
    //spinzerohiggs_anomcoupl_.ghgsgs3[0] = Hwwcoupl[8][0];
    //spinzerohiggs_anomcoupl_.ghgsgs4[0] = Hwwcoupl[9][0];
    */
    spinzerohiggs_anomcoupl_.ghw1_prime[0] = Hwwcoupl[10][0];
    spinzerohiggs_anomcoupl_.ghw1_prime2[0] = Hwwcoupl[11][0];
    spinzerohiggs_anomcoupl_.ghw1_prime3[0] = Hwwcoupl[12][0];
    spinzerohiggs_anomcoupl_.ghw1_prime4[0] = Hwwcoupl[13][0];
    spinzerohiggs_anomcoupl_.ghw1_prime5[0] = Hwwcoupl[14][0];
    spinzerohiggs_anomcoupl_.ghw2_prime[0] = Hwwcoupl[15][0];
    spinzerohiggs_anomcoupl_.ghw2_prime2[0] = Hwwcoupl[16][0];
    spinzerohiggs_anomcoupl_.ghw2_prime3[0] = Hwwcoupl[17][0];
    spinzerohiggs_anomcoupl_.ghw2_prime4[0] = Hwwcoupl[18][0];
    spinzerohiggs_anomcoupl_.ghw2_prime5[0] = Hwwcoupl[19][0];
    spinzerohiggs_anomcoupl_.ghw3_prime[0] = Hwwcoupl[20][0];
    spinzerohiggs_anomcoupl_.ghw3_prime2[0] = Hwwcoupl[21][0];
    spinzerohiggs_anomcoupl_.ghw3_prime3[0] = Hwwcoupl[22][0];
    spinzerohiggs_anomcoupl_.ghw3_prime4[0] = Hwwcoupl[23][0];
    spinzerohiggs_anomcoupl_.ghw3_prime5[0] = Hwwcoupl[24][0];
    spinzerohiggs_anomcoupl_.ghw4_prime[0] = Hwwcoupl[25][0];
    spinzerohiggs_anomcoupl_.ghw4_prime2[0] = Hwwcoupl[26][0];
    spinzerohiggs_anomcoupl_.ghw4_prime3[0] = Hwwcoupl[27][0];
    spinzerohiggs_anomcoupl_.ghw4_prime4[0] = Hwwcoupl[28][0];
    spinzerohiggs_anomcoupl_.ghw4_prime5[0] = Hwwcoupl[29][0];
    //spinzerohiggs_anomcoupl_.ghwgs1_prime2[0] = Hwwcoupl[30][0];
    spinzerohiggs_anomcoupl_.ghw1_prime6[0] = Hwwcoupl[31][0];
    spinzerohiggs_anomcoupl_.ghw1_prime7[0] = Hwwcoupl[32][0];
    spinzerohiggs_anomcoupl_.ghw2_prime6[0] = Hwwcoupl[33][0];
    spinzerohiggs_anomcoupl_.ghw2_prime7[0] = Hwwcoupl[34][0];
    spinzerohiggs_anomcoupl_.ghw3_prime6[0] = Hwwcoupl[35][0];
    spinzerohiggs_anomcoupl_.ghw3_prime7[0] = Hwwcoupl[36][0];
    spinzerohiggs_anomcoupl_.ghw4_prime6[0] = Hwwcoupl[37][0];
    spinzerohiggs_anomcoupl_.ghw4_prime7[0] = Hwwcoupl[38][0];

    spinzerohiggs_anomcoupl_.ghw1[1] =  Hwwcoupl[0][1];
    spinzerohiggs_anomcoupl_.ghw2[1] =  Hwwcoupl[1][1];
    spinzerohiggs_anomcoupl_.ghw3[1] =  Hwwcoupl[2][1];
    spinzerohiggs_anomcoupl_.ghw4[1] =  Hwwcoupl[3][1];
    /*
    spinzerohiggs_anomcoupl_.ghwgs2[1] = Hwwcoupl[4][1];
    spinzerohiggs_anomcoupl_.ghwgs3[1] = Hwwcoupl[5][1];
    spinzerohiggs_anomcoupl_.ghwgs4[1] = Hwwcoupl[6][1];
    //spinzerohiggs_anomcoupl_.ghgsgs2[1] = Hwwcoupl[7][1];
    //spinzerohiggs_anomcoupl_.ghgsgs3[1] = Hwwcoupl[8][1];
    //spinzerohiggs_anomcoupl_.ghgsgs4[1] = Hwwcoupl[9][1];
    */
    spinzerohiggs_anomcoupl_.ghw1_prime[1] = Hwwcoupl[10][1];
    spinzerohiggs_anomcoupl_.ghw1_prime2[1] = Hwwcoupl[11][1];
    spinzerohiggs_anomcoupl_.ghw1_prime3[1] = Hwwcoupl[12][1];
    spinzerohiggs_anomcoupl_.ghw1_prime4[1] = Hwwcoupl[13][1];
    spinzerohiggs_anomcoupl_.ghw1_prime5[1] = Hwwcoupl[14][1];
    spinzerohiggs_anomcoupl_.ghw2_prime[1] = Hwwcoupl[15][1];
    spinzerohiggs_anomcoupl_.ghw2_prime2[1] = Hwwcoupl[16][1];
    spinzerohiggs_anomcoupl_.ghw2_prime3[1] = Hwwcoupl[17][1];
    spinzerohiggs_anomcoupl_.ghw2_prime4[1] = Hwwcoupl[18][1];
    spinzerohiggs_anomcoupl_.ghw2_prime5[1] = Hwwcoupl[19][1];
    spinzerohiggs_anomcoupl_.ghw3_prime[1] = Hwwcoupl[20][1];
    spinzerohiggs_anomcoupl_.ghw3_prime2[1] = Hwwcoupl[21][1];
    spinzerohiggs_anomcoupl_.ghw3_prime3[1] = Hwwcoupl[22][1];
    spinzerohiggs_anomcoupl_.ghw3_prime4[1] = Hwwcoupl[23][1];
    spinzerohiggs_anomcoupl_.ghw3_prime5[1] = Hwwcoupl[24][1];
    spinzerohiggs_anomcoupl_.ghw4_prime[1] = Hwwcoupl[25][1];
    spinzerohiggs_anomcoupl_.ghw4_prime2[1] = Hwwcoupl[26][1];
    spinzerohiggs_anomcoupl_.ghw4_prime3[1] = Hwwcoupl[27][1];
    spinzerohiggs_anomcoupl_.ghw4_prime4[1] = Hwwcoupl[28][1];
    spinzerohiggs_anomcoupl_.ghw4_prime5[1] = Hwwcoupl[29][1];
    //spinzerohiggs_anomcoupl_.ghwgs1_prime2[1] = Hwwcoupl[30][1];
    spinzerohiggs_anomcoupl_.ghw1_prime6[1] = Hwwcoupl[31][1];
    spinzerohiggs_anomcoupl_.ghw1_prime7[1] = Hwwcoupl[32][1];
    spinzerohiggs_anomcoupl_.ghw2_prime6[1] = Hwwcoupl[33][1];
    spinzerohiggs_anomcoupl_.ghw2_prime7[1] = Hwwcoupl[34][1];
    spinzerohiggs_anomcoupl_.ghw3_prime6[1] = Hwwcoupl[35][1];
    spinzerohiggs_anomcoupl_.ghw3_prime7[1] = Hwwcoupl[36][1];
    spinzerohiggs_anomcoupl_.ghw4_prime6[1] = Hwwcoupl[37][1];
    spinzerohiggs_anomcoupl_.ghw4_prime7[1] = Hwwcoupl[38][1];
    //
  }
}
void SetJHUGenSpinZeroVVCouplings(double Hvvcoupl[SIZE_HVV][2], int Hvvcoupl_cqsq[3], double HvvLambda_qsq[4][3], bool useWWcoupl){
  int iWWcoupl = (useWWcoupl ? 1 : 0);
  for (int c=0; c<4; c++){ for (int k=0; k<3; k++) HvvLambda_qsq[c][k] /= 100.; } // GeV units in JHUGen
  __modjhugenmela_MOD_setspinzerovvcouplings(Hvvcoupl, Hvvcoupl_cqsq, HvvLambda_qsq, &iWWcoupl);
}
void SetJHUGenSpinZeroVVCouplings_NoGamma(double Hvvcoupl[SIZE_HVV_VBF][2], int Hvvcoupl_cqsq[3], double HvvLambda_qsq[4][3], bool useWWcoupl){
  int iWWcoupl = (useWWcoupl ? 1 : 0);
  for (int c=0; c<4; c++){ for (int k=0; k<3; k++) HvvLambda_qsq[c][k] /= 100.; } // GeV units in JHUGen
  __modjhugenmela_MOD_setspinzerovvcouplings_nogamma(Hvvcoupl, Hvvcoupl_cqsq, HvvLambda_qsq, &iWWcoupl);
}
void SetJHUGenSpinZeroGGCouplings(double Hggcoupl[SIZE_HGG][2]){ __modjhugenmela_MOD_setspinzeroggcouplings(Hggcoupl); }
void SetJHUGenSpinZeroQQCouplings(double Hqqcoupl[SIZE_HQQ][2]){ __modjhugenmela_MOD_setspinzeroqqcouplings(Hqqcoupl); }
void SetJHUGenSpinOneCouplings(double Zqqcoupl[SIZE_ZQQ][2], double Zvvcoupl[SIZE_ZVV][2]){ __modjhugenmela_MOD_setspinonecouplings(Zqqcoupl, Zvvcoupl); }
void SetJHUGenSpinTwoCouplings(double Gacoupl[SIZE_GGG][2], double Gbcoupl[SIZE_GVV][2], double qLeftRightcoupl[SIZE_GQQ][2]){ __modjhugenmela_MOD_setspintwocouplings(Gacoupl, Gbcoupl, qLeftRightcoupl); }


bool MCFM_chooser(TVar::Process process, TVar::Production production, TVar::LeptonInterference leptonInterf, vector<int> id_dau, vector<int> id_associated){
  unsigned int ndau = id_dau.size();
  unsigned int naparts = id_associated.size();
  unsigned int najets = 0;
  unsigned int naneutrinos = 0;
  unsigned int naleps = 0;
  unsigned int naferms = 0;
  unsigned int naphotons = 0;
  unsigned int nainvalid = 0;
  for (unsigned int ap=0; ap<naparts; ap++){
    if (isALepton(id_associated.at(ap))) naleps++;
    else if (isANeutrino(id_associated.at(ap))) naneutrinos++;
    else if (isAPhoton(id_associated.at(ap))) naphotons++;
    else if (isAJet(id_associated.at(ap))) najets++;
    else nainvalid++;
  }
  naferms = naleps + naneutrinos + najets;
  if (nainvalid>0) cerr << "TUtil::MCFM_chooser: nainvalid=" << nainvalid << ">0! This should not have happened!" << endl;
  naparts -= nainvalid;
  bool definiteInterf = false;
  if (ndau>3) definiteInterf = (
    id_dau.at(0)==id_dau.at(2) && id_dau.at(1)==id_dau.at(3)
    &&
    !isAnUnknownJet(id_dau.at(0)) && !isAnUnknownJet(id_dau.at(2))
    &&
    !isInvalid(id_dau.at(0)) && !isInvalid(id_dau.at(2))
    );

  // VV->4f
  if (
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
      //		vsymfact_.vsymfact=0.25; // MELA FACTOR (Same 0.25 in MCFM 6.7)
      interference_.interference=true;
    }

  }
  else if (
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
  else if (
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

bool My_masscuts(double s[][mxpart],TVar::Process process){

 double minZmassSqr=10*10;

 if(process==TVar::bkgZZ){
   if(s[2][3]< minZmassSqr) return true;
   if(s[4][5]< minZmassSqr) return true;
 }
 return false;	 

}


bool My_smalls(double s[][mxpart],int npart){

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
double SumMatrixElementPDF(
  TVar::Process process, TVar::Production production, TVar::MatrixElement matrixElement,
  event_scales_type* event_scales, MelaIO* RcdME,
  double* flux, double EBEAM, double coupling[SIZE_HVV_FREENORM]
  ){

  int partIncCode=TVar::kNoAssociated; // Do not use associated particles in the pT=0 frame boost
  if (
    ((process==TVar::bkgZZ_SMHiggs || process==TVar::HSMHiggs || process==TVar::bkgZZ) && production==TVar::JJVBF)
    ) partIncCode=TVar::kUseAssociated_Jets; // Use asociated jets in the pT=0 frame boost

  simple_event_record mela_event;
  GetBoostedParticleVectors(
    RcdME->melaCand,
    mela_event.pDaughters,
    mela_event.pAssociated,
    mela_Event.pMothers,
    mela_event.intermediateVid,
    partIncCode
    );

  double xx[2]={ 0 };
  if (!CheckPartonMomFraction(mela_event.pMothers.at(0).second, mela_event.pMothers.at(1).second, xx, TVar::ERROR, EBEAM)) return 0;

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
    double msqjk_sum = SumMEPDF(MomStore[0], MomStore[1], msq, RcdME, TVar::ERROR, EBEAM);
    if (process==TVar::bkgZZ && (production == TVar::ZZQQB || production == TVar::ZZQQB_STU || production == TVar::ZZQQB_S || production == TVar::ZZQQB_TU || production ==TVar::ZZINDEPENDENT)) msqjk = msq[3][7] + msq[7][3]; // all of the unweighted MEs are the same
    else if ((process==TVar::bkgZZ_SMHiggs || process==TVar::HSMHiggs || process==TVar::bkgZZ) && production==TVar::JJVBF) msqjk = msqjk_sum; // MCFM VVH sum
    else msqjk = msq[5][5]; // gg-only

    (*flux)=fbGeV2/(8.*xx[0]*xx[1]*EBEAM*EBEAM);
  }

  if (msqjk != msqjk || flux!=flux){
    cout << "SumMatrixPDF: "<< TVar::ProcessName(process) << " msqjk="  << msqjk << " flux="<< *flux << endl;
    msqjk=0;
    *flux=0;
  }

//  cout << "Before reset: " << scale_.scale << '\t' << facscale_.facscale << endl;
  SetAlphaS(defaultRenScale, defaultFacScale, 1., 1., defaultNloop, defaultNflav, defaultPdflabel); // Protection for other probabilities
//  cout << "Default scale reset: " << scale_.scale << '\t' << facscale_.facscale << endl;
  return msqjk;
}


double JHUGenMatEl(
  TVar::Process process, TVar::Production production, TVar::MatrixElement matrixElement,
  event_scales_type* event_scales, MelaIO* RcdME,
  double EBEAM
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
  GetBoostedParticleVectors(
    RcdME->melaCand,
    mela_event.pDaughters,
    mela_event.pAssociated,
    mela_Event.pMothers,
    mela_event.intermediateVid,
    partIncCode
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
    // If is is known, just assign it.
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
      for (int ix=0; ix<5; ix++){ msq[ix][ix+6]=MatElSq; msq[ix+6][ix]=MatElSq; }
    }
  }
  if (production!=TVar::ZZINDEPENDENT) SumMEPDF(MomStore[0], MomStore[1], msq, RcdME, TVar::ERROR, EBEAM);
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

double HJJMatEl(
  TVar::Process process, TVar::Production production, TVar::MatrixElement matrixElement,
  event_scales_type* event_scales,
  MelaIO* RcdME,
  TVar::VerbosityLevel verbosity,
  double EBEAM
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
  if (!(production == TVar::JJGG || production == TVar::JJVBF || production == TVar::JH)){ cerr << "TUtil::HJJMatEl: Production is not supported!" << endl; return sum_msqjk; }

  // Notice that partIncCode is specific for this subroutine
  int partIncCode=TVar::kUseAssociated_Jets; // Only use associated partons in the pT=0 frame boost
  simple_event_record mela_event;
  GetBoostedParticleVectors(
    RcdME->melaCand,
    mela_event.pDaughters,
    mela_event.pAssociated,
    mela_Event.pMothers,
    mela_event.intermediateVid,
    partIncCode
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
        MatElsq[jsel][isel] = MatElsq_tmp[jsel][isel]; // Assign only those that match gen. info, if present at all.
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
              MatElsq[jsel][isel] += MatElsq_tmp[jsel][isel]*avgfac; // Assign only those that match gen. info, if present at all.
            }
            if (
              (MYIDUP_tmp[2]==0 || (PDGHelpers::isAQuark(MYIDUP_tmp[2]) && MYIDUP_tmp[2]<0))
              &&
              (MYIDUP_tmp[3]==0 || (PDGHelpers::isAQuark(MYIDUP_tmp[3]) && MYIDUP_tmp[3]>0))
              ){
              __modhiggsjj_MOD_evalamp_sbfh_unsymm_sa_select_exact(p4, &isel, &jsel, &ssel, &rsel, MatElsq_tmp);
              MatElsq[jsel][isel] += MatElsq_tmp[jsel][isel]*avgfac; // Assign only those that match gen. info, if present at all.
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
              MatElsq[jsel][isel] += MatElsq_tmp[jsel][isel]; // Assign only those that match gen. info, if present at all.
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
            MatElsq[jsel][isel] += MatElsq_tmp[jsel][isel]*avgfac; // Assign only those that match gen. info, if present at all.
          }
          if (
            (MYIDUP_tmp[2]==0 || ((PDGHelpers::isAGluon(MYIDUP_tmp[2]) && ssel==0) || MYIDUP_tmp[2]==ssel))
            &&
            (MYIDUP_tmp[3]==0 || ((PDGHelpers::isAGluon(MYIDUP_tmp[3]) && rsel==0) || MYIDUP_tmp[3]==rsel))
            ){
            __modhiggsjj_MOD_evalamp_sbfh_unsymm_sa_select_exact(p4, &isel, &jsel, &ssel, &rsel, MatElsq_tmp);
            MatElsq[jsel][isel] += MatElsq_tmp[jsel][isel]*avgfac; // Assign only those that match gen. info, if present at all.
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
              MatElsq[jsel][isel] += MatElsq_tmp[jsel][isel]; // Assign only those that match gen. info, if present at all.
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
              MatElsq[jsel][isel] += MatElsq_tmp[jsel][isel]*avgfac; // Assign only those that match gen. info, if present at all.
            }
            if (
              (MYIDUP_tmp[2]==0 || MYIDUP_tmp[2]==ssel)
              &&
              (MYIDUP_tmp[3]==0 || MYIDUP_tmp[3]==rsel)
              ){
              __modhiggsjj_MOD_evalamp_sbfh_unsymm_sa_select_exact(p4, &isel, &jsel, &ssel, &rsel, MatElsq_tmp);
              MatElsq[jsel][isel] += MatElsq_tmp[jsel][isel]*avgfac; // Assign only those that match gen. info, if present at all.
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
              MatElsq[jsel][isel] += MatElsq_tmp[jsel][isel]*avgfac; // Assign only those that match gen. info, if present at all.
            }
            if (
              (MYIDUP_tmp[2]==0 || MYIDUP_tmp[2]==ssel)
              &&
              (MYIDUP_tmp[3]==0 || MYIDUP_tmp[3]==rsel)
              ){
              __modhiggsjj_MOD_evalamp_sbfh_unsymm_sa_select_exact(p4, &isel, &jsel, &ssel, &rsel, MatElsq_tmp);
              MatElsq[jsel][isel] += MatElsq_tmp[jsel][isel]*avgfac; // Assign only those that match gen. info, if present at all.
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
            MatElsq[jsel][isel] += MatElsq_tmp[jsel][isel]*avgfac; // Assign only those that match gen. info, if present at all.
          }
          if (
            rsel!=ssel
            &&
            (MYIDUP_tmp[2]==0 || MYIDUP_tmp[2]==ssel)
            &&
            (MYIDUP_tmp[3]==0 || MYIDUP_tmp[3]==rsel)
            ){
            __modhiggsjj_MOD_evalamp_sbfh_unsymm_sa_select_exact(p4, &isel, &jsel, &ssel, &rsel, MatElsq_tmp);
            MatElsq[jsel][isel] += MatElsq_tmp[jsel][isel]*avgfac; // Assign only those that match gen. info, if present at all.
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
            MatElsq[jsel][isel] += MatElsq_tmp[jsel][isel]*avgfac; // Assign only those that match gen. info, if present at all.
          }
          if (
            (MYIDUP_tmp[2]==0 || ((PDGHelpers::isAGluon(MYIDUP_tmp[2]) && ssel==0) || MYIDUP_tmp[2]==ssel))
            &&
            (MYIDUP_tmp[3]==0 || ((PDGHelpers::isAGluon(MYIDUP_tmp[3]) && rsel==0) || MYIDUP_tmp[3]==rsel))
            ){
            __modhiggsjj_MOD_evalamp_sbfh_unsymm_sa_select_exact(p4, &isel, &jsel, &ssel, &rsel, MatElsq_tmp);
            MatElsq[jsel][isel] += MatElsq_tmp[jsel][isel]*avgfac; // Assign only those that match gen. info, if present at all.
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
              MatElsq[jsel][isel] += MatElsq_tmp[jsel][isel]; // Assign only those that match gen. info, if present at all.
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
              MatElsq[jsel][isel] += MatElsq_tmp[jsel][isel]*avgfac; // Assign only those that match gen. info, if present at all.
            }
            if (
              (MYIDUP_tmp[2]==0 || MYIDUP_tmp[2]==ssel)
              &&
              (MYIDUP_tmp[3]==0 || MYIDUP_tmp[3]==rsel)
              ){
              __modhiggsjj_MOD_evalamp_sbfh_unsymm_sa_select_exact(p4, &isel, &jsel, &ssel, &rsel, MatElsq_tmp);
              MatElsq[jsel][isel] += MatElsq_tmp[jsel][isel]*avgfac; // Assign only those that match gen. info, if present at all.
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
              MatElsq[jsel][isel] += MatElsq_tmp[jsel][isel]*avgfac; // Assign only those that match gen. info, if present at all.
            }
            if (
              (MYIDUP_tmp[2]==0 || MYIDUP_tmp[2]==ssel)
              &&
              (MYIDUP_tmp[3]==0 || MYIDUP_tmp[3]==rsel)
              ){
              __modhiggsjj_MOD_evalamp_sbfh_unsymm_sa_select_exact(p4, &isel, &jsel, &ssel, &rsel, MatElsq_tmp);
              MatElsq[jsel][isel] += MatElsq_tmp[jsel][isel]*avgfac; // Assign only those that match gen. info, if present at all.
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
            MatElsq[jsel][isel] += MatElsq_tmp[jsel][isel]*avgfac; // Assign only those that match gen. info, if present at all.
          }
          if (
            rsel!=ssel
            &&
            (MYIDUP_tmp[2]==0 || MYIDUP_tmp[2]==ssel)
            &&
            (MYIDUP_tmp[3]==0 || MYIDUP_tmp[3]==rsel)
            ){
            __modhiggsjj_MOD_evalamp_sbfh_unsymm_sa_select_exact(p4, &isel, &jsel, &ssel, &rsel, MatElsq_tmp);
            MatElsq[jsel][isel] += MatElsq_tmp[jsel][isel]*avgfac; // Assign only those that match gen. info, if present at all.
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
            MatElsq[jsel][isel] += MatElsq_tmp[jsel][isel]*avgfac; // Assign only those that match gen. info, if present at all.
          }
          if (
            rsel!=ssel
            &&
            (MYIDUP_tmp[2]==0 || MYIDUP_tmp[2]==ssel)
            &&
            (MYIDUP_tmp[3]==0 || MYIDUP_tmp[3]==rsel)
            ){
            __modhiggsjj_MOD_evalamp_wbfh_unsymm_sa_select_exact(p4, &isel, &jsel, &ssel, &rsel, MatElsq_tmp);
            MatElsq[jsel][isel] += MatElsq_tmp[jsel][isel]*avgfac; // Assign only those that match gen. info, if present at all.
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
                MatElsq[jsel][isel] += MatElsq_tmp[jsel][isel]*avgfac; // Assign only those that match gen. info, if present at all.
              }
              if (
                rsel!=ssel
                &&
                (MYIDUP_tmp[2]==0 || MYIDUP_tmp[2]==ssel)
                &&
                (MYIDUP_tmp[3]==0 || MYIDUP_tmp[3]==rsel)
                ){
                __modhiggsjj_MOD_evalamp_wbfh_unsymm_sa_select_exact(p4, &isel, &jsel, &ssel, &rsel, MatElsq_tmp);
                MatElsq[jsel][isel] += MatElsq_tmp[jsel][isel]*avgfac; // Assign only those that match gen. info, if present at all.
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
            MatElsq[jsel][isel] += MatElsq_tmp[jsel][isel]*avgfac; // Assign only those that match gen. info, if present at all.
          }
          if (
            rsel!=ssel
            &&
            (MYIDUP_tmp[2]==0 || MYIDUP_tmp[2]==ssel)
            &&
            (MYIDUP_tmp[3]==0 || MYIDUP_tmp[3]==rsel)
            ){
            __modhiggsjj_MOD_evalamp_wbfh_unsymm_sa_select_exact(p4, &isel, &jsel, &ssel, &rsel, MatElsq_tmp);
            MatElsq[jsel][isel] += MatElsq_tmp[jsel][isel]*avgfac; // Assign only those that match gen. info, if present at all.
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
                MatElsq[jsel][isel] += MatElsq_tmp[jsel][isel]*avgfac; // Assign only those that match gen. info, if present at all.
              }
              if (
                rsel!=ssel
                &&
                (MYIDUP_tmp[2]==0 || MYIDUP_tmp[2]==ssel)
                &&
                (MYIDUP_tmp[3]==0 || MYIDUP_tmp[3]==rsel)
                ){
                __modhiggsjj_MOD_evalamp_wbfh_unsymm_sa_select_exact(p4, &isel, &jsel, &ssel, &rsel, MatElsq_tmp);
                MatElsq[jsel][isel] += MatElsq_tmp[jsel][isel]*avgfac; // Assign only those that match gen. info, if present at all.
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

  sum_msqjk = SumMEPDF(MomStore[0], MomStore[1], MatElsq, RcdME, verbosity, EBEAM);

  //cout << "Before reset: " << scale_.scale << '\t' << facscale_.facscale << endl;
  SetAlphaS(defaultRenScale, defaultFacScale, 1., 1., defaultNloop, defaultNflav, defaultPdflabel); // Protection for other probabilities
  //cout << "Default scale reset: " << scale_.scale << '\t' << facscale_.facscale << endl;
  return sum_msqjk;
}

double VHiggsMatEl(
  TVar::Process process, TVar::Production production, TVar::MatrixElement matrixElement,
  event_scales_type* event_scales, MelaIO* RcdME,
  TVar::VerbosityLevel verbosity,
  double EBEAM
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

  // Notice that partIncCode is specific for this subroutine
  int partIncCode;
  if (production == TVar::Had_ZH || production == TVar::Had_WH) partIncCode=TVar::kUseAssociated_Jets; // Only use associated partons
  else if (production == TVar::Lep_ZH || production == TVar::Lep_WH) partIncCode=TVar::kUseAssociated_Leptons; // Only use associated leptons(+)neutrinos
  simple_event_record mela_event;
  GetBoostedParticleVectors(
    RcdME->melaCand,
    mela_event.pDaughters,
    mela_event.pAssociated,
    mela_Event.pMothers,
    mela_event.intermediateVid,
    partIncCode
    );
  if (mela_event.pAssociated.size()!=2){ cerr << "TUtil::VHiggsMatEl: Number of associated particles is 0!" << endl; return sum_msqjk; }

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
        } // End ZH case
        else if (production==TVar::Lep_WH || production==TVar::Had_WH){
          vh_ids[0] = incoming1;
          if (MYIDUP_prod[0]!=0 && MYIDUP_prod[0]!=vh_ids[0]) continue;

          for (int incoming2 = -nf; incoming2 < 0; incoming2++){
            if (abs(incoming2)==abs(incoming1) || TMath::Sign(1, incoming1)==TMath::Sign(1, incoming2)) continue;

            vh_ids[1] = incoming2;
            if (MYIDUP_prod[1]!=0 && MYIDUP_prod[1]!=vh_ids[1]) continue;

            for (int outgoing1=-nf; outgoing1<=nf; outgoing1++){
              for (int outgoing2=-nf; outgoing2<=nf; outgoing2++){
                if (abs(outgoing2)==abs(outgoing1) || TMath::Sign(1, outgoing1)==TMath::Sign(1, outgoing2)) continue;

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
          } // End loop over incoming2
        } // End WH case
      } // End loop over incoming1
    } // End loop over h56
  } // End loop over h01

  sum_msqjk = SumMEPDF(MomStore[0], MomStore[1], MatElsq, RcdME, verbosity, EBEAM);

  //cout << "Before reset: " << scale_.scale << '\t' << facscale_.facscale << endl;
  SetAlphaS(defaultRenScale, defaultFacScale, 1., 1., defaultNloop, defaultNflav, defaultPdflabel); // Protection for other probabilities
  //cout << "Default scale reset: " << scale_.scale << '\t' << facscale_.facscale << endl;
  return sum_msqjk;
}


// LEFT HERE: NEED TO MERGE OLD AND NEW FORTRAN CODE, MOST NOTABLY ADDING EVALXSEC_PP*
double TTHiggsMatEl(
  TVar::Process process, TVar::Production production, TVar::MatrixElement matrixElement,
  event_scales_type* event_scales, MelaIO* RcdME,
  TVar::VerbosityLevel verbosity,
  double EBEAM,
  int topDecay, int topProcess
  ){
  double sum_msqjk=0;
  double MatElsq[nmsq][nmsq]={ { 0 } };
  double p4[13][4]={ { 0 } };

  if (matrixElement!=TVar::JHUGen){ cerr << "TUtil::TTHiggsMatEl: Non-JHUGen MEs are not supported." << endl; return sum_msqjk; }
  if (production!=TVar::ttH && production!=TVar::bbH){ cerr << "TUtil::TTHiggsMatEl: Only ttH or bbH are supported." << endl; return sum_msqjk; }

  // LEFT HERE

  for (int i = 0; i < 2; i++){
    p4[i][0] = -p[i].Energy()/100.;
    p4[i][1] = -p[i].Px()/100.;
    p4[i][2] = -p[i].Py()/100.;
    p4[i][3] = -p[i].Pz()/100.;
  }
  // 2-10 are what?
  for (int i = 2; i < 11; i++){
    p4[i][0] = p[i].Energy()/100.;
    p4[i][1] = p[i].Px()/100.;
    p4[i][2] = p[i].Py()/100.;
    p4[i][3] = p[i].Pz()/100.;
  }

  if (verbosity >= TVar::DEBUG){
    for (int i=0; i<11; i++){
      cout << "p4[" << i << "] = ";
      for (int jj=0; jj<4; jj++) cout << p4[i][jj] << '\t';
      cout << endl;
    }
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


  if (production == TVar::ttH)     __modttbh_MOD_evalxsec_pp_ttbh(p4, &topProcess, MatElsq);
  else if (production ==TVar::bbH) __modttbh_MOD_evalxsec_pp_bbbh(p4, &topProcess, MatElsq);


  sum_msqjk = SumMEPDF(MomStore[0], MomStore[1], MatElsq, RcdME, verbosity, EBEAM);

  //cout << "Before reset: " << scale_.scale << '\t' << facscale_.facscale << endl;
  SetAlphaS(defaultRenScale, defaultFacScale, 1., 1., defaultNloop, defaultNflav, defaultPdflabel); // Protection for other probabilities
  //cout << "Default scale reset: " << scale_.scale << '\t' << facscale_.facscale << endl;
  return sum_msqjk;
}


// CheckPartonMomFraction computes xx[0:1] based on p0, p1
bool CheckPartonMomFraction(const TLorentzVector p0, const TLorentzVector p1, double xx[2], TVar::VerbosityLevel verbosity, double EBEAM){
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
void ComputePDF(const TLorentzVector p0, const TLorentzVector p1, double fx1[nmsq], double fx2[nmsq], TVar::VerbosityLevel verbosity, double EBEAM){
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
double SumMEPDF(const TLorentzVector p0, const TLorentzVector p1, double msq[nmsq][nmsq], MelaIO* RcdME, TVar::VerbosityLevel verbosity, double EBEAM){
  double fx1[nmsq]={ 0 };
  double fx2[nmsq]={ 0 };
  //double wgt_msq[nmsq][nmsq]={ { 0 } };

  ComputePDF(p0, p1, fx1, fx2, verbosity, EBEAM);
  RcdME->setPartonWeights(fx1, fx2);
  RcdME->setMEArray(msq,true);
  RcdME->computeWeightedMEArray();
  //RcdME->getWeightedMEArray(wgt_msq);
  return RcdME->getSumME();
}

// GetBoostedParticleVectors decomposes the MELACandidate object melaCand into mothers, daughters, and associated particles in the pT=0 frame for the particular the useAssociatedCode, which is a product of TVar::kUseAssociated* (prime numbers)
// If no mothers are present, it assigns mothers as Z>0, Z<0. If they are present, it orders them as "incoming" qbar-q / qbar-g / g-q / g-g
// Associated particles passed are different based on the useAssociatedCode!=1. If useAssociatedCode==1, no associated particles are passed.
// intermediateVids are needed to keep track of the decay mode. TVar::Process or TVar::Production do not keep track of V/f decay modes.
void GetBoostedParticleVectors(
  MELACandidate* melaCand,
  pair<vector<int>, vector<TLorentzVector>>& pDaughters,
  pair<vector<int>, vector<TLorentzVector>>& pAssociated,
  pair<vector<int>, vector<TLorentzVector>>& pMothers,
  vector<int>& intermediateVid,
  int useAssociatedCode
  ){
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
          if (Vdau_i!=0) daughters.push_back(
            pair<vector<int>, vector<TLorentzVector>>(Vdau_i->id, Vdau_i->p4)
            );
        }
        if (idtmp!=0 || Vdau->getNDaughters()>0){ // Avoid "empty" intermediate Vs of the MELACandidate object
          if (Vdau->getNDaughters()>=2 && PDGHelpers::isAPhoton(idtmp)) idtmp=23; // Special case to avoid V->2f with HVVmass==Zeromass setting (could happen by mistake)
          idVstar.push_back(idtmp);
        }
      }
    }
  }

  pair<vector<int>, vector<TLorentzVector>> associated;
  if (code%TVar::kUseAssociated_Leptons){
    for (int ip=0; ip<melaCand->getNAssociatedLeptons(); ip++){
      if (melaCand->getAssociatedLepton(ip)!=0) associated.push_back(
        pair<vector<int>, vector<TLorentzVector>>(melaCand->getAssociatedLepton(ip)->id, melaCand->getAssociatedLepton(ip)->p4)
        );
    }
  }
  if (code%TVar::kUseAssociated_Photons){
    for (int ip=0; ip<melaCand->getNAssociatedPhotons(); ip++){
      if (melaCand->getAssociatedPhoton(ip)!=0) associated.push_back(
        pair<vector<int>, vector<TLorentzVector>>(melaCand->getAssociatedPhoton(ip)->id, melaCand->getAssociatedPhoton(ip)->p4)
        );
    }
  }
  if (code%TVar::kUseAssociated_Jets){
    for (int ip=0; ip<melaCand->getNAssociatedJets(); ip++){
      if (melaCand->getAssociatedJet(ip)!=0) associated.push_back(
        pair<vector<int>, vector<TLorentzVector>>(melaCand->getAssociatedJet(ip)->id, melaCand->getAssociatedJet(ip)->p4)
        );
    }
  }

  TLorentzVector pTotal(0, 0, 0, 0);
  for (unsigned int ip=0; ip<daughters.size(); ip++) pTotal = pTotal + daughters.at(ip).second;
  if (useAssociated){ for (unsigned int ip=0; ip<associated.size(); ip++) pTotal = pTotal + associated.at(ip).second; }

  double qX = pTotal.X();
  double qY = pTotal.Y();
  double qE = pTotal.T();;
  if ((qX*qX+qY*qY)>0.){
    TVector3 boostV(-qX/qE, -qY/qE, 0.);
    for (unsigned int ip=0; ip<daughters.size(); ip++) daughters.at(ip).second.Boost(boostV);
    for (unsigned int ip=0; ip<associated.size(); ip++) associated.at(ip).second.Boost(boostV);
    pTotal.Boost(boostV);
  }

  double sysPz= pTotal.Z();
  double sysE = pTotal.T();
  double pz0 = (sysE+sysPz)/2.;
  double pz1 = -(sysE-sysPz)/2.;
  int motherId[2]={ 0, 0 };
  if (melaCand->getNMothers()==2){
    if ((melaCand->getMother(0)->p4).Z()<0.){ // Swap by Z to get the correct momentum matching default M1, M2
      double pztmp = pz1;
      pz1=pz0;
      pz0=pztmp;
    }
    for (int ip=0; ip<2; ip++) motherId[ip]=melaCand->getMother(ip)->id;
    if (motherId[0]<0){ // Swap to achieve ordering as "incoming" q-qbar
      double pztmp = pz1;
      pz1=pz0;
      pz0=pztmp;
      for (int ip=0; ip<2; ip++) motherId[1-ip]=melaCand->getMother(ip)->id;
    }
  }
  TLorentzVector pM[2];
  pM[0].SetPxPyPzE(0., 0., pz0, TMath::Abs(pz0));
  pM[1].SetPxPyPzE(0., 0., pz1, TMath::Abs(pz1));

  pDaughters.clear();
  for (unsigned int ip=0; ip<daughters.size(); ip++) pDaughters.push_back(daughters.at(ip));
  intermediateVid.clear();
  for (unsigned int ip=0; ip<idVstar.size(); ip++) intermediateVid.push_back(idVstar.at(ip));
  pAssociated.clear();
  for (unsigned int ip=0; ip<associated.size(); ip++) pAssociated.push_back(associated.at(ip));
  pMothers.clear();
  for (unsigned int ip=0; ip<2; ip++){
    int idmother = 0; // In case it is unknown.
    if (melaCand->getNMothers()==2) idmother = motherId[ip];
    pMothers.push_back(
      pair<vector<int>, vector<TLorentzVector>>(idmother, pM[ip])
      );
  }
}


