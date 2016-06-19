#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include <fstream>
#include <cstdlib>
#include <iomanip>
#include "TMath.h"
#include "TLorentzVector.h"
#include "TLorentzRotation.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TString.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include <ZZMatrixElement/MELA/interface/Mela.h>


using namespace RooFit;
using namespace std;


void testME_Dec_MCFM_Ping(int flavor=2, int useMothers=0){
  int erg_tev=13;
  float mPOLE=125.;
  float wPOLE=4.07e-3;

  TVar::VerbosityLevel verbosity = TVar::DEBUG;
  Mela mela(erg_tev, mPOLE, verbosity);
  if (verbosity>=TVar::DEBUG) cout << "Mela is initialized" << endl;
  //mela.resetMCFM_EWKParameters(1.16639E-05, 1./128., 79.9549392, 91.1876, 0.23119);
  if (verbosity>=TVar::DEBUG) cout << "Mela candidate decay mode initializing" << endl;
  mela.setCandidateDecayMode(TVar::CandidateDecay_ZZ);
  if (verbosity>=TVar::DEBUG) cout << "Mela candidate decay mode initialized" << endl;

  float GenLep1Id=0, GenLep2Id=0, GenLep3Id=0, GenLep4Id=0;
  /*
  const int nEntries = 3;
  double l1_array[nEntries][4] = {
  {1365.4973807340846,        10.289826593755228,        25.205694382277809,       -1365.2259480507332},
  {238.65751023078761,        9.2808858562825005,        15.827726043466324,       -237.95116187061188},
  {101.52463181523598,        27.359569630718468,      -0.90299073100241323,       -97.764458892691749}
  };
  double l2_array[nEntries][4] = {
  {22.786181013986834,      -0.15136300982222117,      -0.90077551414353962,       -22.767866345236371},
  {101.67043553688544,        2.1169375132239789,       0.77953005873937187,       -101.64540506443268},
  {24.717634703436786,       -1.1722249478288802,       -5.9599387484197646,       -23.959684558009428}
  };
  double l3_array[nEntries][4] = {
  {1895.7562628816693,        25.837804322120054,       -28.821887970086259,       -1895.3610513294620},
  {317.81904277258536,        2.5882005498984775,        21.352807448987718,       -317.09037005377883},
  {180.10885677707822,       -6.7240759244122792,        35.742176497019194,       -176.39865053838915}
  };
  double l4_array[nEntries][4] = {
  {471.71918486784784,       -35.976267906053060,        4.5169691019519895,       -470.32360615864354},
  {95.655512770627581,       -13.986023919404957,       -37.960063551193414,       -86.679881365440792},
  {49.137252081251319,       -19.463268758477309,       -28.879247017597017,       -34.664676589120688}
  };
  */
  const int nEntries = 6;
  double l1_array[nEntries][4] ={
    { 51.374202, 25.924766, 12.290178, 42.616376 },
    { 51.374202, 25.924766, 12.290178, 42.616376 },
    { 1365.4973807340846, 10.289826593755228, 25.205694382277809, -1365.2259480507332 }, // Massless
    { 1365.4973807340846, 10.289826593755228, 25.205694382277809, -1365.2259480507332 }, // Massless
    { 1365.4973848483, 10.289826593755228, 25.205694382277809, -1365.2259480507332 }, // Muon via E
    { 1365.4973848483, 10.289826593755228, 25.205694382277809, -1365.2259480507332 } // Muon via E
  };
  double l2_array[nEntries][4] ={
    { 271.875752, 70.427173, -11.138146, 261.769598 },
    { 21.481452, 9.489680, -9.336587, 16.858699 },
    { 22.786181013986834, -0.15136300982222117, -0.90077551414353962, -22.767866345236371 },
    { 471.71918486784784, -35.976267906053060, 4.5169691019519895, -470.32360615864354 },
    { 22.7864275656, -0.15136300982222117, -0.90077551414353962, -22.767866345236371 },
    { 471.7191967775, -35.976267906053060, 4.5169691019519895, -470.32360615864354 }
  };
  double l3_array[nEntries][4] ={
    { 75.823478, -16.640412, 23.246999, 70.227220 },
    { 75.823478, -16.640412, 23.246999, 70.227220 },
    { 1895.7562628816693, 25.837804322120054, -28.821887970086259, -1895.3610513294620 },
    { 1895.7562628816693, 25.837804322120054, -28.821887970086259, -1895.3610513294620 },
    { 1895.7562658451, 25.837804322120054, -28.821887970086259, -1895.3610513294620 },
    { 1895.7562658451, 25.837804322120054, -28.821887970086259, -1895.3610513294620 }
  };
  double l4_array[nEntries][4] ={
    { 21.481452, 9.489680, -9.336587, 16.858699 },
    { 271.875752, 70.427173, -11.138146, 261.769598 },
    { 471.71918486784784, -35.976267906053060, 4.5169691019519895, -470.32360615864354 },
    { 22.786181013986834, -0.15136300982222117, -0.90077551414353962, -22.767866345236371 },
    { 471.7191967775, -35.976267906053060, 4.5169691019519895, -470.32360615864354 },
    { 22.7864275656, -0.15136300982222117, -0.90077551414353962, -22.767866345236371 }
  };

  for (int ev = 2; ev < 3; ev++){
    if (flavor == 2){
      GenLep1Id=13;
      GenLep2Id=-13;
      GenLep3Id=11;
      GenLep4Id=-11;
    }
    else if (flavor == 1){
      GenLep1Id=11;
      GenLep2Id=-11;
      GenLep3Id=11;
      GenLep4Id=-11;
    }
    else if (flavor == 0){
      GenLep1Id=13;
      GenLep2Id=-13;
      GenLep3Id=13;
      GenLep4Id=-13;
    }
    else if (flavor == 3){
      GenLep1Id=14;
      GenLep2Id=-14;
      GenLep3Id=13;
      GenLep4Id=-13;
    }
    else if (flavor == 4){
      GenLep1Id=0;
      GenLep2Id=-0;
      GenLep3Id=1;
      GenLep4Id=-1;
    }
    else if (flavor == 5){
      GenLep1Id=1;
      GenLep2Id=-1;
      GenLep3Id=2;
      GenLep4Id=-2;
    }

    int idOrdered[4] ={ static_cast<int>(GenLep1Id), static_cast<int>(GenLep2Id), static_cast<int>(GenLep3Id), static_cast<int>(GenLep4Id) };
    int idOrdered_WW[4] ={ 11, -12, -11, 12 };
    TLorentzVector pOrdered[4];
    pOrdered[0].SetXYZT(l1_array[ev][1], l1_array[ev][2], l1_array[ev][3], l1_array[ev][0]);
    pOrdered[1].SetXYZT(l2_array[ev][1], l2_array[ev][2], l2_array[ev][3], l2_array[ev][0]);
    pOrdered[2].SetXYZT(l3_array[ev][1], l3_array[ev][2], l3_array[ev][3], l3_array[ev][0]);
    pOrdered[3].SetXYZT(l4_array[ev][1], l4_array[ev][2], l4_array[ev][3], l4_array[ev][0]);
    SimpleParticleCollection_t daughters_WW;
    for (unsigned int idau=0; idau<4; idau++) daughters_WW.push_back(SimpleParticle_t(idOrdered_WW[idau], pOrdered[idau]));
    SimpleParticleCollection_t daughters_WWasZZ;
    for (unsigned int iv=0; iv<2; iv++){ for (int ivd=0; ivd<2; ivd++) daughters_WWasZZ.push_back(SimpleParticle_t(idOrdered_WW[iv+2*ivd], pOrdered[iv+2*ivd])); }
    SimpleParticleCollection_t daughters_ZZ;
    for (unsigned int idau=0; idau<4; idau++) daughters_ZZ.push_back(SimpleParticle_t(idOrdered[idau], pOrdered[idau]));

    TLorentzVector pOrdered_ZG[3];
    pOrdered_ZG[0]=pOrdered[0];
    pOrdered_ZG[1]=pOrdered[1];
    pOrdered_ZG[2]=pOrdered[2]+pOrdered[3];
    SimpleParticleCollection_t daughters_ZG;
    for (unsigned int idau=0; idau<2; idau++) daughters_ZG.push_back(SimpleParticle_t(idOrdered[idau], pOrdered_ZG[idau]));
    for (unsigned int idau=2; idau<3; idau++) daughters_ZG.push_back(SimpleParticle_t(22, pOrdered_ZG[idau]));

    // Some gymnastics to get pseudo-mothers for these events
    TLorentzVector pTotal = pOrdered[0]+pOrdered[1]+pOrdered[2]+pOrdered[3];
    TLorentzVector pTotal_dummy=pTotal;
    TLorentzVector pTotal_perp(pTotal.X(), pTotal.Y(), 0, pTotal.T());
    pTotal_dummy.Boost(-pTotal_perp.BoostVector());
    TLorentzVector pMothers[2];
    pMothers[0].SetXYZT(0, 0, (pTotal_dummy.Z()+pTotal_dummy.T())/2., (pTotal_dummy.Z()+pTotal_dummy.T())/2.);
    pMothers[1].SetXYZT(0, 0, (pTotal_dummy.Z()-pTotal_dummy.T())/2., (-pTotal_dummy.Z()+pTotal_dummy.T())/2.);
    for (int im=0; im<2; im++) pMothers[im].Boost(pTotal_perp.BoostVector());
    SimpleParticleCollection_t mothers_QQB;
    for (unsigned int im=0; im<2; im++) mothers_QQB.push_back(SimpleParticle_t(1-2*im, pMothers[im]));
    SimpleParticleCollection_t mothers_GG;
    for (unsigned int im=0; im<2; im++) mothers_GG.push_back(SimpleParticle_t(21, pMothers[im]));
    SimpleParticleCollection_t* mothersPtr=0;
    if (useMothers==1) mothersPtr=&mothers_GG;
    else if (useMothers==2) mothersPtr=&mothers_QQB;

    TLorentzVector pOrdered_GG[3];
    pOrdered_GG[0]=pOrdered[0]+pOrdered[1];
    pOrdered_GG[1]=pOrdered[2]+pOrdered[3];
    SimpleParticleCollection_t daughters_GG;
    for (unsigned int idau=0; idau<2; idau++) daughters_GG.push_back(SimpleParticle_t(22, pOrdered_GG[idau]));

    mela.setCandidateDecayMode(TVar::CandidateDecay_WW);
    mela.setInputEvent(&daughters_WW, (SimpleParticleCollection_t*)0, mothersPtr, (mothersPtr!=0));
    mela.setCandidateDecayMode(TVar::CandidateDecay_ZZ);
    mela.setInputEvent(&daughters_WWasZZ, (SimpleParticleCollection_t*)0, mothersPtr, (mothersPtr!=0));
    mela.setCandidateDecayMode(TVar::CandidateDecay_ZZ);
    mela.setInputEvent(&daughters_ZZ, (SimpleParticleCollection_t*)0, mothersPtr, (mothersPtr!=0));
    mela.setInputEvent(&daughters_ZG, (SimpleParticleCollection_t*)0, mothersPtr, (mothersPtr!=0));
    mela.setCandidateDecayMode(TVar::CandidateDecay_GG);
    mela.setInputEvent(&daughters_GG, (SimpleParticleCollection_t*)0, mothersPtr, (mothersPtr!=0));

    int cindex;

    /***** WW *****/
    cindex=0;
    mela.setCurrentCandidateFromIndex(cindex);

    float pVAMCFM_qqWW_bkg;
    mela.setProcess(TVar::bkgWW, TVar::MCFM, TVar::ZZQQB);
    mela.computeP(pVAMCFM_qqWW_bkg, false);
    cout << "pVAMCFM_qqWW_bkg: " << pVAMCFM_qqWW_bkg << '\n' << endl;

    float pVAMCFM_ggVV_total;
    mela.setProcess(TVar::bkgWWZZ_SMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.computeP(pVAMCFM_ggVV_total, false);
    cout << "pVAMCFM_ggVV_total: " << pVAMCFM_ggVV_total << '\n' << endl;
    float pVAMCFM_ggVV_bkg;
    mela.setProcess(TVar::bkgWWZZ, TVar::MCFM, TVar::ZZGG);
    mela.computeP(pVAMCFM_ggVV_bkg, false);
    cout << "pVAMCFM_ggVV_bkg: " << pVAMCFM_ggVV_bkg << '\n' << endl;
    float pVAMCFM_ggVV_sig;
    mela.setProcess(TVar::HSMHiggs_WWZZ, TVar::MCFM, TVar::ZZGG);
    mela.computeP(pVAMCFM_ggVV_sig, false);
    cout << "pVAMCFM_ggVV_sig: " << pVAMCFM_ggVV_sig << '\n' << endl;

    float pVAMCFM_ggWW_total;
    mela.setProcess(TVar::bkgWW_SMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.computeP(pVAMCFM_ggWW_total, false);
    cout << "pVAMCFM_ggWW_total: " << pVAMCFM_ggWW_total << '\n' << endl;
    float pVAMCFM_ggWW_bkg;
    mela.setProcess(TVar::bkgWW, TVar::MCFM, TVar::ZZGG);
    mela.computeP(pVAMCFM_ggWW_bkg, false);
    cout << "pVAMCFM_ggWW_bkg: " << pVAMCFM_ggWW_bkg << '\n' << endl;
    float pVAMCFM_ggWW_sig;
    mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.computeP(pVAMCFM_ggWW_sig, false);
    cout << "pVAMCFM_ggWW_sig: " << pVAMCFM_ggWW_sig << '\n' << endl;

    /***** WW (as ZZ) *****/
    cindex=1;
    mela.setCurrentCandidateFromIndex(cindex);

    float pVAMCFM_ggVV_fromZZ_total;
    mela.setProcess(TVar::bkgWWZZ_SMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.computeP(pVAMCFM_ggVV_fromZZ_total, false);
    cout << "pVAMCFM_ggVV_fromZZ_total: " << pVAMCFM_ggVV_fromZZ_total << '\n' << endl;
    float pVAMCFM_ggVV_fromZZ_bkg;
    mela.setProcess(TVar::bkgWWZZ, TVar::MCFM, TVar::ZZGG);
    mela.computeP(pVAMCFM_ggVV_fromZZ_bkg, false);
    cout << "pVAMCFM_ggVV_fromZZ_bkg: " << pVAMCFM_ggVV_fromZZ_bkg << '\n' << endl;
    float pVAMCFM_ggVV_fromZZ_sig;
    mela.setProcess(TVar::HSMHiggs_WWZZ, TVar::MCFM, TVar::ZZGG);
    mela.computeP(pVAMCFM_ggVV_fromZZ_sig, false);
    cout << "pVAMCFM_ggVV_fromZZ_sig: " << pVAMCFM_ggVV_fromZZ_sig << '\n' << endl;

    float pVAMCFM_ggZZ_total;
    mela.setProcess(TVar::bkgZZ_SMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.computeP(pVAMCFM_ggZZ_total, false);
    cout << "pVAMCFM_ggZZ_total: " << pVAMCFM_ggZZ_total << '\n' << endl;
    float pVAMCFM_ggZZ_bkg;
    mela.setProcess(TVar::bkgZZ, TVar::MCFM, TVar::ZZGG);
    mela.computeP(pVAMCFM_ggZZ_bkg, false);
    cout << "pVAMCFM_ggZZ_bkg: " << pVAMCFM_ggZZ_bkg << '\n' << endl;
    float pVAMCFM_ggWWasZZ_sig;
    mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.computeP(pVAMCFM_ggWWasZZ_sig, false);
    cout << "pVAMCFM_ggWWasZZ_sig: " << pVAMCFM_ggWWasZZ_sig << '\n' << endl;

    /***** ZZ *****/
    cindex=2;
    mela.setCurrentCandidateFromIndex(cindex);

    float pVAMCFM_qqZZ_bkg;
    mela.setProcess(TVar::bkgZZ, TVar::MCFM, TVar::ZZQQB);
    mela.computeP(pVAMCFM_qqZZ_bkg, false);
    cout << "pVAMCFM_qqZZ_bkg: " << pVAMCFM_qqZZ_bkg << '\n' << endl;

    float pVAMCFM_qqZJJ_bkg;
    mela.setProcess(TVar::bkgZJJ, TVar::MCFM, TVar::ZZQQB);
    mela.computeP(pVAMCFM_qqZJJ_bkg, false);
    cout << "pVAMCFM_qqZJJ_bkg: " << pVAMCFM_qqZJJ_bkg << '\n' << endl;

    float pVAMCFM_ggZZ_sig;
    mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
    for (int ii = 0; ii < SIZE_HVV; ii++){ for (int jj = 0; jj < 2; jj++)   mela.selfDHzzcoupl[0][ii][jj] = 0; }
    mela.setMelaHiggsWidth(wPOLE);
    mela.setMelaLeptonInterference(TVar::InterfOn);
    mela.computeP(pVAMCFM_ggZZ_sig, false);
    cout << "pVAMCFM_ggZZ_sig: " << pVAMCFM_ggZZ_sig << '\n' << endl;

    float pVAMCFM_ZZ_sig_selfDg1;
    mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
    for (int ii = 0; ii < SIZE_HVV; ii++){ for (int jj = 0; jj < 2; jj++)   mela.selfDHzzcoupl[0][ii][jj] = 0; }
    mela.selfDHzzcoupl[0][0][0]=1.0;
    mela.setMelaHiggsWidth(wPOLE);
    mela.setMelaLeptonInterference(TVar::InterfOn);
    mela.computeP(pVAMCFM_ZZ_sig_selfDg1, false);
    cout << "pVAMCFM_ggZZ_sig_selfDg1: " << pVAMCFM_ZZ_sig_selfDg1 << '\n' << endl;

    /***** ZG *****/
    cindex=3;
    mela.setCurrentCandidateFromIndex(cindex);
    float pVAMCFM_ZG_bkg;
    mela.setProcess(TVar::bkgZGamma, TVar::MCFM, TVar::ZZQQB);
    mela.computeP(pVAMCFM_ZG_bkg, false);
    cout << "pVAMCFM_qqZG_bkg: " << pVAMCFM_ZG_bkg << '\n' << endl;

    if (verbosity>=TVar::DEBUG){
      cout << "Removing Mela candidates\nSummary:" << endl;
      for (int ic=0; ic<mela.getNCandidates(); ic++){
        cout << "Candidate " << ic << endl;
        mela.setCurrentCandidateFromIndex(ic);
        TUtil::PrintCandidateSummary(mela.getCurrentCandidate());
      }
      cout << endl;
    }
    mela.resetInputEvent();
    cout << "Removed..." << endl;
  }
}
void testME_ProdDec_MCFM_Ping(int flavor=2){
  int erg_tev=13;
  float mPOLE=125.;
  float wPOLE=4.07e-3;

  TVar::VerbosityLevel verbosity = TVar::DEBUG;
  Mela mela(erg_tev, mPOLE, verbosity);
  if (verbosity>=TVar::DEBUG) cout << "Mela is initialized" << endl;
  //mela.resetMCFM_EWKParameters(1.16639E-05, 1./128., 79.9549392, 91.1876, 0.23119);
  if (verbosity>=TVar::DEBUG) cout << "Mela candidate decay mode initializing" << endl;
  mela.setCandidateDecayMode(TVar::CandidateDecay_ZZ);
  if (verbosity>=TVar::DEBUG) cout << "Mela candidate decay mode initialized" << endl;

  float GenLep1Id=0, GenLep2Id=0, GenLep3Id=0, GenLep4Id=0;
  const int nEntries = 6;
  double a1_array[nEntries][4] ={
    { 1365.4973807340846, 10.289826593755228, 25.205694382277809, -1365.2259480507332 },
    { 238.65751023078761, 9.2808858562825005, 15.827726043466324, -237.95116187061188 },
    { 101.52463181523598, 27.359569630718468, -0.90299073100241323, -97.764458892691749 },
    { 22.786181013986834, -0.15136300982222117, -0.90077551414353962, -22.767866345236371 },
    { 101.67043553688544, 2.1169375132239789, 0.77953005873937187, -101.64540506443268 },
    { 24.717634703436786, -1.1722249478288802, -5.9599387484197646, -23.959684558009428 }
  };
  double a2_array[nEntries][4] ={
    { 1895.7562628816693, 25.837804322120054, -28.821887970086259, -1895.3610513294620 },
    { 317.81904277258536, 2.5882005498984775, 21.352807448987718, -317.09037005377883 },
    { 180.10885677707822, -6.7240759244122792, 35.742176497019194, -176.39865053838915 },
    { 471.71918486784784, -35.976267906053060, 4.5169691019519895, -470.32360615864354 },
    { 95.655512770627581, -13.986023919404957, -37.960063551193414, -86.679881365440792 },
    { 49.137252081251319, -19.463268758477309, -28.879247017597017, -34.664676589120688 }
  };
  double l1_array[nEntries][4] ={
    { 51.374202, 25.924766, 12.290178, 42.616376 },
    { 51.374202, 25.924766, 12.290178, 42.616376 },
    { 1365.4973807340846, 10.289826593755228, 25.205694382277809, -1365.2259480507332 }, // Massless
    { 1365.4973807340846, 10.289826593755228, 25.205694382277809, -1365.2259480507332 }, // Massless
    { 1365.4973848483, 10.289826593755228, 25.205694382277809, -1365.2259480507332 }, // Muon via E
    { 1365.4973848483, 10.289826593755228, 25.205694382277809, -1365.2259480507332 } // Muon via E
  };
  double l2_array[nEntries][4] ={
    { 271.875752, 70.427173, -11.138146, 261.769598 },
    { 21.481452, 9.489680, -9.336587, 16.858699 },
    { 22.786181013986834, -0.15136300982222117, -0.90077551414353962, -22.767866345236371 },
    { 471.71918486784784, -35.976267906053060, 4.5169691019519895, -470.32360615864354 },
    { 22.7864275656, -0.15136300982222117, -0.90077551414353962, -22.767866345236371 },
    { 471.7191967775, -35.976267906053060, 4.5169691019519895, -470.32360615864354 }
  };
  double l3_array[nEntries][4] ={
    { 75.823478, -16.640412, 23.246999, 70.227220 },
    { 75.823478, -16.640412, 23.246999, 70.227220 },
    { 1895.7562628816693, 25.837804322120054, -28.821887970086259, -1895.3610513294620 },
    { 1895.7562628816693, 25.837804322120054, -28.821887970086259, -1895.3610513294620 },
    { 1895.7562658451, 25.837804322120054, -28.821887970086259, -1895.3610513294620 },
    { 1895.7562658451, 25.837804322120054, -28.821887970086259, -1895.3610513294620 }
  };
  double l4_array[nEntries][4] ={
    { 21.481452, 9.489680, -9.336587, 16.858699 },
    { 271.875752, 70.427173, -11.138146, 261.769598 },
    { 471.71918486784784, -35.976267906053060, 4.5169691019519895, -470.32360615864354 },
    { 22.786181013986834, -0.15136300982222117, -0.90077551414353962, -22.767866345236371 },
    { 471.7191967775, -35.976267906053060, 4.5169691019519895, -470.32360615864354 },
    { 22.7864275656, -0.15136300982222117, -0.90077551414353962, -22.767866345236371 }
  };

  if (flavor == 2){
    GenLep1Id=13;
    GenLep2Id=-13;
    GenLep3Id=11;
    GenLep4Id=-11;
  }
  else if (flavor == 1){
    GenLep1Id=11;
    GenLep2Id=-11;
    GenLep3Id=11;
    GenLep4Id=-11;
  }
  else if (flavor == 0){
    GenLep1Id=13;
    GenLep2Id=-13;
    GenLep3Id=13;
    GenLep4Id=-13;
  }
  else if (flavor == 3){
    GenLep1Id=14;
    GenLep2Id=-14;
    GenLep3Id=13;
    GenLep4Id=-13;
  }
  else if (flavor == 4){
    GenLep1Id=0;
    GenLep2Id=-0;
    GenLep3Id=1;
    GenLep4Id=-1;
  }
  else if (flavor == 5){
    GenLep1Id=1;
    GenLep2Id=-1;
    GenLep3Id=2;
    GenLep4Id=-2;
  }

  for (int ev = 0; ev < 1; ev++){
    SimpleParticleCollection_t aparticles;
    TLorentzVector pAPart[4];
    pAPart[0].SetXYZT(a1_array[ev][1], a1_array[ev][2], a1_array[ev][3], a1_array[ev][0]);
    pAPart[1].SetXYZT(a2_array[ev][1], a2_array[ev][2], a2_array[ev][3], a2_array[ev][0]);
    for (unsigned int iap=0; iap<2; iap++) aparticles.push_back(SimpleParticle_t(1-2*iap, pAPart[iap]));
    for (unsigned int iap=0; iap<2; iap++) aparticles.push_back(SimpleParticle_t((1-2*iap)*13, pAPart[iap]));

    int idOrdered[4] ={ static_cast<int>(GenLep1Id), static_cast<int>(GenLep2Id), static_cast<int>(GenLep3Id), static_cast<int>(GenLep4Id) };
    int idOrdered_WW[4] ={ 11, -12, -11, 12 };
    TLorentzVector pOrdered[4];
    pOrdered[0].SetXYZT(l1_array[ev][1], l1_array[ev][2], l1_array[ev][3], l1_array[ev][0]);
    pOrdered[1].SetXYZT(l2_array[ev][1], l2_array[ev][2], l2_array[ev][3], l2_array[ev][0]);
    pOrdered[2].SetXYZT(l3_array[ev][1], l3_array[ev][2], l3_array[ev][3], l3_array[ev][0]);
    pOrdered[3].SetXYZT(l4_array[ev][1], l4_array[ev][2], l4_array[ev][3], l4_array[ev][0]);
    SimpleParticleCollection_t daughters_WW;
    for (unsigned int idau=0; idau<4; idau++) daughters_WW.push_back(SimpleParticle_t(idOrdered_WW[idau], pOrdered[idau]));
    SimpleParticleCollection_t daughters_WWasZZ;
    for (unsigned int iv=0; iv<2; iv++){ for (int ivd=0; ivd<2; ivd++) daughters_WWasZZ.push_back(SimpleParticle_t(idOrdered_WW[iv+2*ivd], pOrdered[iv+2*ivd])); }
    SimpleParticleCollection_t daughters_ZZ;
    for (unsigned int idau=0; idau<4; idau++) daughters_ZZ.push_back(SimpleParticle_t(idOrdered[idau], pOrdered[idau]));

    mela.setCandidateDecayMode(TVar::CandidateDecay_WW);
    mela.setInputEvent(&daughters_WW, &aparticles, (SimpleParticleCollection_t*)0, false);
    mela.setCandidateDecayMode(TVar::CandidateDecay_ZZ);
    mela.setInputEvent(&daughters_WWasZZ, &aparticles, (SimpleParticleCollection_t*)0, false);
    mela.setCandidateDecayMode(TVar::CandidateDecay_ZZ);
    mela.setInputEvent(&daughters_ZZ, &aparticles, (SimpleParticleCollection_t*)0, false);


    int cindex;

    /***** WW *****/
    cindex=0;
    mela.setCurrentCandidateFromIndex(cindex);

    float pVAMCFM_wbfVV_total;
    mela.setProcess(TVar::bkgWWZZ_SMHiggs, TVar::MCFM, TVar::JJVBF);
    mela.computeProdDecP(pVAMCFM_wbfVV_total, true);
    cout << "pVAMCFM_wbfVV_total: " << pVAMCFM_wbfVV_total << '\n' << endl;
    float pVAMCFM_wbfVV_bkg;
    mela.setProcess(TVar::bkgWWZZ, TVar::MCFM, TVar::JJVBF);
    mela.computeProdDecP(pVAMCFM_wbfVV_bkg, true);
    cout << "pVAMCFM_wbfVV_bkg: " << pVAMCFM_wbfVV_bkg << '\n' << endl;
    float pVAMCFM_wbfVV_sig;
    mela.setProcess(TVar::HSMHiggs_WWZZ, TVar::MCFM, TVar::JJVBF);
    mela.computeProdDecP(pVAMCFM_wbfVV_sig, true);
    cout << "pVAMCFM_wbfVV_sig: " << pVAMCFM_wbfVV_sig << '\n' << endl;

    float pVAMCFM_wbfWW_total;
    mela.setProcess(TVar::bkgWW_SMHiggs, TVar::MCFM, TVar::JJVBF);
    mela.computeProdDecP(pVAMCFM_wbfWW_total, true);
    cout << "pVAMCFM_wbfWW_total: " << pVAMCFM_wbfWW_total << '\n' << endl;
    float pVAMCFM_wbfWW_bkg;
    mela.setProcess(TVar::bkgWW, TVar::MCFM, TVar::JJVBF);
    mela.computeProdDecP(pVAMCFM_wbfWW_bkg, true);
    cout << "pVAMCFM_wbfWW_bkg: " << pVAMCFM_wbfWW_bkg << '\n' << endl;
    float pVAMCFM_wbfWW_sig;
    mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::JJVBF);
    mela.computeProdDecP(pVAMCFM_wbfWW_sig, true);
    cout << "pVAMCFM_wbfWW_sig: " << pVAMCFM_wbfWW_sig << '\n' << endl;

    /***** WW (as ZZ) *****/
    cindex=1;
    mela.setCurrentCandidateFromIndex(cindex);
    float pVAMCFM_wbfZZ_total;
    mela.setProcess(TVar::bkgZZ_SMHiggs, TVar::MCFM, TVar::JJVBF);
    mela.computeProdDecP(pVAMCFM_wbfZZ_total, true);
    cout << "pVAMCFM_wbfZZ_total: " << pVAMCFM_wbfZZ_total << '\n' << endl;
    float pVAMCFM_wbfZZ_bkg;
    mela.setProcess(TVar::bkgZZ, TVar::MCFM, TVar::JJVBF);
    mela.computeProdDecP(pVAMCFM_wbfZZ_bkg, true);
    cout << "pVAMCFM_wbfZZ_bkg: " << pVAMCFM_wbfZZ_bkg << '\n' << endl;
    float pVAMCFM_wbfWWasZZ_sig;
    mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::JJVBF);
    mela.computeProdDecP(pVAMCFM_wbfWWasZZ_sig, true);
    cout << "pVAMCFM_wbfWWasZZ_sig: " << pVAMCFM_wbfWWasZZ_sig << '\n' << endl;

    /***** ZZ *****/
    cindex=2;
    mela.setCurrentCandidateFromIndex(cindex);
    float pVAMCFM_wbfZZ_sig;
    mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::JJVBF);
    mela.computeProdDecP(pVAMCFM_wbfZZ_sig, true);
    cout << "pVAMCFM_wbfZZ_sig: " << pVAMCFM_wbfZZ_sig << '\n' << endl;

    if (verbosity>=TVar::DEBUG){
      cout << "Removing Mela candidates\nSummary:" << endl;
      for (int ic=0; ic<mela.getNCandidates(); ic++){
        cout << "Candidate " << ic << endl;
        mela.setCurrentCandidateFromIndex(ic);
        TUtil::PrintCandidateSummary(mela.getCurrentCandidate());
      }
      cout << endl;
    }
    mela.resetInputEvent();
    cout << "Removed..." << endl;
  }
}


void testME_Dec_JHUGenMCFM_Ping(int flavor=2){
  int erg_tev=13;
  float mPOLE=125.;
  float wPOLE=4.07e-3;

  TVar::VerbosityLevel verbosity = TVar::DEBUG;
  Mela mela(erg_tev, mPOLE, verbosity);
  if (verbosity>=TVar::DEBUG) cout << "Mela is initialized" << endl;
  mela.resetMCFM_EWKParameters(1.16639E-05, 1./128., 80.399, 91.1876, 0.23119);

  float GenLep1Id=0, GenLep2Id=0, GenLep3Id=0, GenLep4Id=0;
  /*
  const int nEntries = 3;
  double l1_array[nEntries][4] = {
  {1365.4973807340846,        10.289826593755228,        25.205694382277809,       -1365.2259480507332},
  {238.65751023078761,        9.2808858562825005,        15.827726043466324,       -237.95116187061188},
  {101.52463181523598,        27.359569630718468,      -0.90299073100241323,       -97.764458892691749}
  };
  double l2_array[nEntries][4] = {
  {22.786181013986834,      -0.15136300982222117,      -0.90077551414353962,       -22.767866345236371},
  {101.67043553688544,        2.1169375132239789,       0.77953005873937187,       -101.64540506443268},
  {24.717634703436786,       -1.1722249478288802,       -5.9599387484197646,       -23.959684558009428}
  };
  double l3_array[nEntries][4] = {
  {1895.7562628816693,        25.837804322120054,       -28.821887970086259,       -1895.3610513294620},
  {317.81904277258536,        2.5882005498984775,        21.352807448987718,       -317.09037005377883},
  {180.10885677707822,       -6.7240759244122792,        35.742176497019194,       -176.39865053838915}
  };
  double l4_array[nEntries][4] = {
  {471.71918486784784,       -35.976267906053060,        4.5169691019519895,       -470.32360615864354},
  {95.655512770627581,       -13.986023919404957,       -37.960063551193414,       -86.679881365440792},
  {49.137252081251319,       -19.463268758477309,       -28.879247017597017,       -34.664676589120688}
  };
  */
  const int nEntries = 6;
  double l1_array[nEntries][4] ={
    { 51.374202, 25.924766, 12.290178, 42.616376 },
    { 51.374202, 25.924766, 12.290178, 42.616376 },
    { 1365.4973807340846, 10.289826593755228, 25.205694382277809, -1365.2259480507332 }, // Massless
    { 1365.4973807340846, 10.289826593755228, 25.205694382277809, -1365.2259480507332 }, // Massless
    { 1365.4973848483, 10.289826593755228, 25.205694382277809, -1365.2259480507332 }, // Muon via E
    { 1365.4973848483, 10.289826593755228, 25.205694382277809, -1365.2259480507332 } // Muon via E
  };
  double l2_array[nEntries][4] ={
    { 271.875752, 70.427173, -11.138146, 261.769598 },
    { 21.481452, 9.489680, -9.336587, 16.858699 },
    { 22.786181013986834, -0.15136300982222117, -0.90077551414353962, -22.767866345236371 },
    { 471.71918486784784, -35.976267906053060, 4.5169691019519895, -470.32360615864354 },
    { 22.7864275656, -0.15136300982222117, -0.90077551414353962, -22.767866345236371 },
    { 471.7191967775, -35.976267906053060, 4.5169691019519895, -470.32360615864354 }
  };
  double l3_array[nEntries][4] ={
    { 75.823478, -16.640412, 23.246999, 70.227220 },
    { 75.823478, -16.640412, 23.246999, 70.227220 },
    { 1895.7562628816693, 25.837804322120054, -28.821887970086259, -1895.3610513294620 },
    { 1895.7562628816693, 25.837804322120054, -28.821887970086259, -1895.3610513294620 },
    { 1895.7562658451, 25.837804322120054, -28.821887970086259, -1895.3610513294620 },
    { 1895.7562658451, 25.837804322120054, -28.821887970086259, -1895.3610513294620 }
  };
  double l4_array[nEntries][4] ={
    { 21.481452, 9.489680, -9.336587, 16.858699 },
    { 271.875752, 70.427173, -11.138146, 261.769598 },
    { 471.71918486784784, -35.976267906053060, 4.5169691019519895, -470.32360615864354 },
    { 22.786181013986834, -0.15136300982222117, -0.90077551414353962, -22.767866345236371 },
    { 471.7191967775, -35.976267906053060, 4.5169691019519895, -470.32360615864354 },
    { 22.7864275656, -0.15136300982222117, -0.90077551414353962, -22.767866345236371 }
  };

  for (int ev = 0; ev < 1; ev++){
    if (flavor == 2){
      GenLep1Id=13;
      GenLep2Id=-13;
      GenLep3Id=11;
      GenLep4Id=-11;
    }
    else if (flavor == 1){
      GenLep1Id=11;
      GenLep2Id=-11;
      GenLep3Id=11;
      GenLep4Id=-11;
    }
    else if (flavor == 0){
      GenLep1Id=13;
      GenLep2Id=-13;
      GenLep3Id=13;
      GenLep4Id=-13;
    }
    else if (flavor == 3){
      GenLep1Id=14;
      GenLep2Id=-14;
      GenLep3Id=13;
      GenLep4Id=-13;
    }
    else if (flavor == 4){
      GenLep1Id=0;
      GenLep2Id=-0;
      GenLep3Id=1;
      GenLep4Id=-1;
    }
    else if (flavor == 5){
      GenLep1Id=1;
      GenLep2Id=-1;
      GenLep3Id=2;
      GenLep4Id=-2;
    }

    int idOrdered[4] ={ static_cast<int>(GenLep1Id), static_cast<int>(GenLep2Id), static_cast<int>(GenLep3Id), static_cast<int>(GenLep4Id) };
    int idOrdered_WW[4] ={ 11, -12, -11, 12 };
    TLorentzVector pOrdered[4];
    pOrdered[0].SetXYZT(l1_array[ev][1], l1_array[ev][2], l1_array[ev][3], l1_array[ev][0]);
    pOrdered[1].SetXYZT(l2_array[ev][1], l2_array[ev][2], l2_array[ev][3], l2_array[ev][0]);
    pOrdered[2].SetXYZT(l3_array[ev][1], l3_array[ev][2], l3_array[ev][3], l3_array[ev][0]);
    pOrdered[3].SetXYZT(l4_array[ev][1], l4_array[ev][2], l4_array[ev][3], l4_array[ev][0]);
    SimpleParticleCollection_t daughters_WW;
    for (unsigned int idau=0; idau<4; idau++) daughters_WW.push_back(SimpleParticle_t(idOrdered_WW[idau], pOrdered[idau]));
    SimpleParticleCollection_t daughters_WWasZZ;
    for (unsigned int iv=0; iv<2; iv++){ for (int ivd=0; ivd<2; ivd++) daughters_WWasZZ.push_back(SimpleParticle_t(idOrdered_WW[iv+2*ivd], pOrdered[iv+2*ivd])); }
    SimpleParticleCollection_t daughters_ZZ;
    for (unsigned int idau=0; idau<4; idau++) daughters_ZZ.push_back(SimpleParticle_t(idOrdered[idau], pOrdered[idau]));

    TLorentzVector pOrdered_ZG[3];
    pOrdered_ZG[0]=pOrdered[0];
    pOrdered_ZG[1]=pOrdered[1];
    pOrdered_ZG[2]=pOrdered[2]+pOrdered[3];
    SimpleParticleCollection_t daughters_ZG;
    for (unsigned int idau=0; idau<2; idau++) daughters_ZG.push_back(SimpleParticle_t(idOrdered[idau], pOrdered_ZG[idau]));
    for (unsigned int idau=2; idau<3; idau++) daughters_ZG.push_back(SimpleParticle_t(22, pOrdered_ZG[idau]));

    TLorentzVector pOrdered_GG[3];
    pOrdered_GG[0]=pOrdered[0]+pOrdered[1];
    pOrdered_GG[1]=pOrdered[2]+pOrdered[3];
    SimpleParticleCollection_t daughters_GG;
    for (unsigned int idau=0; idau<2; idau++) daughters_GG.push_back(SimpleParticle_t(22, pOrdered_GG[idau]));

    mela.setCandidateDecayMode(TVar::CandidateDecay_WW);
    mela.setInputEvent(&daughters_WW, (SimpleParticleCollection_t*)0, (SimpleParticleCollection_t*)0, false);
    mela.setCandidateDecayMode(TVar::CandidateDecay_ZZ);
    mela.setInputEvent(&daughters_ZZ, (SimpleParticleCollection_t*)0, (SimpleParticleCollection_t*)0, false);
    mela.setInputEvent(&daughters_ZG, (SimpleParticleCollection_t*)0, (SimpleParticleCollection_t*)0, false);
    mela.setCandidateDecayMode(TVar::CandidateDecay_GG);
    mela.setInputEvent(&daughters_GG, (SimpleParticleCollection_t*)0, (SimpleParticleCollection_t*)0, false);

    int cindex;

    /***** WW *****/
    cindex=0;
    mela.setCurrentCandidateFromIndex(cindex);

    // JHUGen ME
    float pVAJHUGen_ggWW_SM_sig;
    mela.setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::ZZGG);
    mela.computeP(pVAJHUGen_ggWW_SM_sig, true);
    cout << "pVAJHUGen_ggWW_SM_sig: " << pVAJHUGen_ggWW_SM_sig << '\n' << endl;

    float pVAJHUGen_ggWW_0pm_sig;
    mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
    mela.selfDHggcoupl[0][0]=1.;
    mela.selfDHzzcoupl[0][0][0]=1.;
    mela.computeP(pVAJHUGen_ggWW_0pm_sig, true);
    cout << "pVAJHUGen_ggWW_0pm_sig: " << pVAJHUGen_ggWW_0pm_sig << '\n' << endl;

    float pVAJHUGen_ggWW_PS_sig;
    mela.setProcess(TVar::H0minus, TVar::JHUGen, TVar::ZZGG);
    mela.computeP(pVAJHUGen_ggWW_PS_sig, true);
    cout << "pVAJHUGen_ggWW_PS_sig: " << pVAJHUGen_ggWW_PS_sig << '\n' << endl;

    float pVAJHUGen_ggWW_0m_sig;
    mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
    mela.selfDHggcoupl[0][0]=1.;
    mela.selfDHzzcoupl[0][3][0]=1.;
    mela.computeP(pVAJHUGen_ggWW_0m_sig, true);
    cout << "pVAJHUGen_ggWW_0m_sig: " << pVAJHUGen_ggWW_0m_sig << '\n' << endl;

    // MCFM ME
    float pVAMCFM_ggWW_SM_total;
    mela.setProcess(TVar::bkgWW_SMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.computeP(pVAMCFM_ggWW_SM_total, true);
    cout << "pVAMCFM_ggWW_SM_total: " << pVAMCFM_ggWW_SM_total << '\n' << endl;

    float pVAMCFM_ggWW_0pm_total;
    mela.setProcess(TVar::bkgWW_SMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.selfDHggcoupl[0][0]=1.;
    mela.selfDHzzcoupl[0][0][0]=1.;
    mela.computeP(pVAMCFM_ggWW_0pm_total, true);
    cout << "pVAMCFM_ggWW_0pm_total: " << pVAMCFM_ggWW_0pm_total << '\n' << endl;

    float pVAMCFM_ggWW_SM_sig;
    mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.computeP(pVAMCFM_ggWW_SM_sig, true);
    cout << "pVAMCFM_ggWW_SM_sig: " << pVAMCFM_ggWW_SM_sig << '\n' << endl;

    float pVAMCFM_ggWW_0pm_sig;
    mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.selfDHggcoupl[0][0]=1.;
    mela.selfDHzzcoupl[0][0][0]=1.;
    mela.computeP(pVAMCFM_ggWW_0pm_sig, true);
    cout << "pVAMCFM_ggWW_0pm_sig: " << pVAMCFM_ggWW_0pm_sig << '\n' << endl;

    float pVAMCFM_ggWW_0m_sig;
    mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.selfDHggcoupl[0][0]=1.;
    mela.selfDHzzcoupl[0][3][0]=1.;
    mela.computeP(pVAMCFM_ggWW_0m_sig, true);
    cout << "pVAMCFM_ggWW_0m_sig: " << pVAMCFM_ggWW_0m_sig << '\n' << endl;

    // MCFM HM ME
    float pVAMCFM_HM_ggWW_SM_total;
    mela.setProcess(TVar::bkgWW_SMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.setMelaHiggsMassWidth(-1, 0, 0); mela.setMelaHiggsMassWidth(mPOLE, -1, 1);
    mela.computeP(pVAMCFM_HM_ggWW_SM_total, true);
    cout << "pVAMCFM_HM_ggWW_SM_total: " << pVAMCFM_HM_ggWW_SM_total << '\n' << endl;

    float pVAMCFM_HM_ggWW_0pm_total;
    mela.setProcess(TVar::bkgWW_SMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.setMelaHiggsMassWidth(-1, 0, 0); mela.setMelaHiggsMassWidth(mPOLE, -1, 1);
    mela.selfDHggcoupl[0][0]=1.;
    mela.selfDHzzcoupl[1][0][0]=1.;
    mela.computeP(pVAMCFM_HM_ggWW_0pm_total, true);
    cout << "pVAMCFM_HM_ggWW_0pm_total: " << pVAMCFM_HM_ggWW_0pm_total << '\n' << endl;

    float pVAMCFM_HM_ggWW_SM_sig;
    mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.setMelaHiggsMassWidth(-1, 0, 0); mela.setMelaHiggsMassWidth(mPOLE, -1, 1);
    mela.computeP(pVAMCFM_HM_ggWW_SM_sig, true);
    cout << "pVAMCFM_HM_ggWW_SM_sig: " << pVAMCFM_HM_ggWW_SM_sig << '\n' << endl;

    float pVAMCFM_HM_ggWW_0pm_sig;
    mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.setMelaHiggsMassWidth(-1, 0, 0); mela.setMelaHiggsMassWidth(mPOLE, -1, 1);
    mela.selfDHggcoupl[0][0]=1.;
    mela.selfDHzzcoupl[1][0][0]=1.;
    mela.computeP(pVAMCFM_HM_ggWW_0pm_sig, true);
    cout << "pVAMCFM_HM_ggWW_0pm_sig: " << pVAMCFM_HM_ggWW_0pm_sig << '\n' << endl;

    float pVAMCFM_HM_ggWW_0m_sig;
    mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.setMelaHiggsMassWidth(-1, 0, 0); mela.setMelaHiggsMassWidth(mPOLE, -1, 1);
    mela.selfDHggcoupl[0][0]=1.;
    mela.selfDHzzcoupl[1][3][0]=1.;
    mela.computeP(pVAMCFM_HM_ggWW_0m_sig, true);
    cout << "pVAMCFM_HM_ggWW_0m_sig: " << pVAMCFM_HM_ggWW_0m_sig << '\n' << endl;

    /***** ZZ *****/
    cindex=1;
    mela.setCurrentCandidateFromIndex(cindex);

    // JHUGen ME
    float pVAJHUGen_ggZZ_SM_sig;
    mela.setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::ZZGG);
    mela.computeP(pVAJHUGen_ggZZ_SM_sig, true);
    cout << "pVAJHUGen_ggZZ_SM_sig: " << pVAJHUGen_ggZZ_SM_sig << '\n' << endl;

    float pVAJHUGen_ggZZ_0pm_sig;
    mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
    mela.selfDHggcoupl[0][0]=1.;
    mela.selfDHzzcoupl[0][0][0]=1.;
    mela.computeP(pVAJHUGen_ggZZ_0pm_sig, true);
    cout << "pVAJHUGen_ggZZ_0pm_sig: " << pVAJHUGen_ggZZ_0pm_sig << '\n' << endl;

    float pVAJHUGen_ggZZ_PS_sig;
    mela.setProcess(TVar::H0minus, TVar::JHUGen, TVar::ZZGG);
    mela.computeP(pVAJHUGen_ggZZ_PS_sig, true);
    cout << "pVAJHUGen_ggZZ_PS_sig: " << pVAJHUGen_ggZZ_PS_sig << '\n' << endl;

    float pVAJHUGen_ggZZ_0m_sig;
    mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
    mela.selfDHggcoupl[0][0]=1.;
    mela.selfDHzzcoupl[0][3][0]=1.;
    mela.computeP(pVAJHUGen_ggZZ_0m_sig, true);
    cout << "pVAJHUGen_ggZZ_0m_sig: " << pVAJHUGen_ggZZ_0m_sig << '\n' << endl;

    // MCFM ME
    float pVAMCFM_ggZZ_SM_total;
    mela.setProcess(TVar::bkgZZ_SMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.computeP(pVAMCFM_ggZZ_SM_total, true);
    cout << "pVAMCFM_ggZZ_SM_total: " << pVAMCFM_ggZZ_SM_total << '\n' << endl;

    float pVAMCFM_ggZZ_0pm_total;
    mela.setProcess(TVar::bkgZZ_SMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.selfDHggcoupl[0][0]=1.;
    mela.selfDHzzcoupl[0][0][0]=1.;
    mela.computeP(pVAMCFM_ggZZ_0pm_total, true);
    cout << "pVAMCFM_ggZZ_0pm_total: " << pVAMCFM_ggZZ_0pm_total << '\n' << endl;

    float pVAMCFM_ggZZ_SM_sig;
    mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.computeP(pVAMCFM_ggZZ_SM_sig, true);
    cout << "pVAMCFM_ggZZ_SM_sig: " << pVAMCFM_ggZZ_SM_sig << '\n' << endl;

    float pVAMCFM_ggZZ_0pm_sig;
    mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.selfDHggcoupl[0][0]=1.;
    mela.selfDHzzcoupl[0][0][0]=1.;
    mela.computeP(pVAMCFM_ggZZ_0pm_sig, true);
    cout << "pVAMCFM_ggZZ_0pm_sig: " << pVAMCFM_ggZZ_0pm_sig << '\n' << endl;

    float pVAMCFM_ggZZ_0m_sig;
    mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.selfDHggcoupl[0][0]=1.;
    mela.selfDHzzcoupl[0][3][0]=1.;
    mela.computeP(pVAMCFM_ggZZ_0m_sig, true);
    cout << "pVAMCFM_ggZZ_0m_sig: " << pVAMCFM_ggZZ_0m_sig << '\n' << endl;

    // MCFM HM ME
    float pVAMCFM_HM_ggZZ_SM_total;
    mela.setProcess(TVar::bkgZZ_SMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.setMelaHiggsMassWidth(-1, 0, 0); mela.setMelaHiggsMassWidth(mPOLE, -1, 1);
    mela.computeP(pVAMCFM_HM_ggZZ_SM_total, true);
    cout << "pVAMCFM_HM_ggZZ_SM_total: " << pVAMCFM_HM_ggZZ_SM_total << '\n' << endl;

    float pVAMCFM_HM_ggZZ_0pm_total;
    mela.setProcess(TVar::bkgZZ_SMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.setMelaHiggsMassWidth(-1, 0, 0); mela.setMelaHiggsMassWidth(mPOLE, -1, 1);
    mela.selfDHggcoupl[0][0]=1.;
    mela.selfDHzzcoupl[1][0][0]=1.;
    mela.computeP(pVAMCFM_HM_ggZZ_0pm_total, true);
    cout << "pVAMCFM_HM_ggZZ_0pm_total: " << pVAMCFM_HM_ggZZ_0pm_total << '\n' << endl;

    float pVAMCFM_HM_ggZZ_SM_sig;
    mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.setMelaHiggsMassWidth(-1, 0, 0); mela.setMelaHiggsMassWidth(mPOLE, -1, 1);
    mela.computeP(pVAMCFM_HM_ggZZ_SM_sig, true);
    cout << "pVAMCFM_HM_ggZZ_SM_sig: " << pVAMCFM_HM_ggZZ_SM_sig << '\n' << endl;

    float pVAMCFM_HM_ggZZ_0pm_sig;
    mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.setMelaHiggsMassWidth(-1, 0, 0); mela.setMelaHiggsMassWidth(mPOLE, -1, 1);
    mela.selfDHggcoupl[0][0]=1.;
    mela.selfDHzzcoupl[1][0][0]=1.;
    mela.computeP(pVAMCFM_HM_ggZZ_0pm_sig, true);
    cout << "pVAMCFM_HM_ggZZ_0pm_sig: " << pVAMCFM_HM_ggZZ_0pm_sig << '\n' << endl;

    float pVAMCFM_HM_ggZZ_0m_sig;
    mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.setMelaHiggsMassWidth(-1, 0, 0); mela.setMelaHiggsMassWidth(mPOLE, -1, 1);
    mela.selfDHggcoupl[0][0]=1.;
    mela.selfDHzzcoupl[1][3][0]=1.;
    mela.computeP(pVAMCFM_HM_ggZZ_0m_sig, true);
    cout << "pVAMCFM_HM_ggZZ_0m_sig: " << pVAMCFM_HM_ggZZ_0m_sig << '\n' << endl;

    /***** ZG *****/
    cindex=2;
    mela.setCurrentCandidateFromIndex(cindex);
    float pVAJHUGen_ggZG_SM_sig;
    mela.setProcess(TVar::H0_Zgs, TVar::JHUGen, TVar::ZZGG);
    mela.computeP(pVAJHUGen_ggZG_SM_sig, true);
    cout << "pVAJHUGen_ggZG_SM_sig: " << pVAJHUGen_ggZG_SM_sig << '\n' << endl;

    float pVAJHUGen_ggZG_0pm_sig;
    mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
    mela.selfDHggcoupl[0][0]=1.;
    mela.selfDHzzcoupl[0][4][0]=1.;
    mela.computeP(pVAJHUGen_ggZG_0pm_sig, true);
    cout << "pVAJHUGen_ggZG_0pm_sig: " << pVAJHUGen_ggZG_0pm_sig << '\n' << endl;

    float pVAJHUGen_ggZG_PS_sig;
    mela.setProcess(TVar::H0_Zgs_PS, TVar::JHUGen, TVar::ZZGG);
    mela.computeP(pVAJHUGen_ggZG_PS_sig, true);
    cout << "pVAJHUGen_ggZG_PS_sig: " << pVAJHUGen_ggZG_PS_sig << '\n' << endl;

    float pVAJHUGen_ggZG_0m_sig;
    mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
    mela.selfDHggcoupl[0][0]=1.;
    mela.selfDHzzcoupl[0][6][0]=1.;
    mela.computeP(pVAJHUGen_ggZG_0m_sig, true);
    cout << "pVAJHUGen_ggZG_0m_sig: " << pVAJHUGen_ggZG_0m_sig << '\n' << endl;

    /***** GG *****/
    cindex=3;
    mela.setCurrentCandidateFromIndex(cindex);
    float pVAJHUGen_ggGG_SM_sig;
    mela.setProcess(TVar::H0_gsgs, TVar::JHUGen, TVar::ZZGG);
    mela.computeP(pVAJHUGen_ggGG_SM_sig, true);
    cout << "pVAJHUGen_ggGG_SM_sig: " << pVAJHUGen_ggGG_SM_sig << '\n' << endl;

    float pVAJHUGen_ggGG_0pm_sig;
    mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
    mela.selfDHggcoupl[0][0]=1.;
    mela.selfDHzzcoupl[0][7][0]=1.;
    mela.computeP(pVAJHUGen_ggGG_0pm_sig, true);
    cout << "pVAJHUGen_ggGG_0pm_sig: " << pVAJHUGen_ggGG_0pm_sig << '\n' << endl;

    float pVAJHUGen_ggGG_PS_sig;
    mela.setProcess(TVar::H0_gsgs_PS, TVar::JHUGen, TVar::ZZGG);
    mela.computeP(pVAJHUGen_ggGG_PS_sig, true);
    cout << "pVAJHUGen_ggGG_PS_sig: " << pVAJHUGen_ggGG_PS_sig << '\n' << endl;

    float pVAJHUGen_ggGG_0m_sig;
    mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
    mela.selfDHggcoupl[0][0]=1.;
    mela.selfDHzzcoupl[0][9][0]=1.;
    mela.computeP(pVAJHUGen_ggGG_0m_sig, true);
    cout << "pVAJHUGen_ggGG_0m_sig: " << pVAJHUGen_ggGG_0m_sig << '\n' << endl;


    if (verbosity>=TVar::DEBUG){
      cout << "Removing Mela candidates\nSummary:" << endl;
      for (int ic=0; ic<mela.getNCandidates(); ic++){
        cout << "Candidate " << ic << endl;
        mela.setCurrentCandidateFromIndex(ic);
        TUtil::PrintCandidateSummary(mela.getCurrentCandidate());
      }
      cout << endl;
    }
    mela.resetInputEvent();
    cout << "Removed..." << endl;
  }
}

void testME_Dec_FullSim(int flavor=2){
  int erg_tev=8;
  float mPOLE=125.6;
  float wPOLE=4.07e-3;
  TString TREE_NAME = "SelectedTree";

  TVar::VerbosityLevel verbosity = TVar::ERROR;
  Mela mela(erg_tev, mPOLE, verbosity);
  if (verbosity>=TVar::DEBUG) cout << "Mela is initialized" << endl;
  mela.resetMCFM_EWKParameters(1.16639E-05, 1./128., 80.399, 91.1876, 0.23119);

  TString cinput_main = "/scratch0/hep/ianderso/CJLST/140519/PRODFSR_8TeV";
  TFile* finput = new TFile(Form("%s/%s/HZZ4lTree_powheg15jhuGenV3-0PMH125.6.root", cinput_main.Data(), (flavor>=2 ? "2mu2e" : "4e")), "read");
  TFile* foutput = new TFile(Form("HZZ4lTree_powheg15jhuGenV3-0PMH125.6_%s_OriginalMEv2ValidationTestOnly.root", (flavor>=2 ? "2mu2e" : "4e")), "recreate");

  float p0plus_VAJHU;
  float p0hplus_VAJHU;
  float p0minus_VAJHU;
  float p0_g1prime2_VAJHU;
  float pg1g1prime2_VAJHU;
  //float pg1g1prime2_pi2_VAJHU;
  float pg1g2_VAJHU;
  float pg1g2_pi2_VAJHU;
  float pg1g4_VAJHU;
  float pg1g4_pi2_VAJHU;

  float p0plus_VAJHU_NEW;
  float p0hplus_VAJHU_NEW;
  float p0minus_VAJHU_NEW;
  float p0_g1prime2_VAJHU_NEW;
  float pg1g1prime2_VAJHU_NEW;
  float pg1g1prime2_pi2_VAJHU_NEW;
  float pg1g2_VAJHU_NEW;
  float pg1g2_pi2_VAJHU_NEW;
  float pg1g4_VAJHU_NEW;
  float pg1g4_pi2_VAJHU_NEW;

  float p0hplus_VAMCFM_NEW;
  float p0minus_VAMCFM_NEW;
  float p0_g1prime2_VAMCFM_NEW;
  float pg1g1prime2_VAMCFM_NEW;
  float pg1g1prime2_pi2_VAMCFM_NEW;
  float pg1g2_VAMCFM_NEW;
  float pg1g2_pi2_VAMCFM_NEW;
  float pg1g4_VAMCFM_NEW;
  float pg1g4_pi2_VAMCFM_NEW;

  float p0hplus_VAMCFM_ratio_NEW;
  float p0minus_VAMCFM_ratio_NEW;
  float p0_g1prime2_VAMCFM_ratio_NEW;
  float pg1g1prime2_VAMCFM_ratio_NEW;
  float pg1g1prime2_pi2_VAMCFM_ratio_NEW;
  float pg1g2_VAMCFM_ratio_NEW;
  float pg1g2_pi2_VAMCFM_ratio_NEW;
  float pg1g4_VAMCFM_ratio_NEW;
  float pg1g4_pi2_VAMCFM_ratio_NEW;

  float p0hplus_VAJHU_ratio_NEW;
  float p0minus_VAJHU_ratio_NEW;
  float p0_g1prime2_VAJHU_ratio_NEW;
  float pg1g1prime2_VAJHU_ratio_NEW;
  float pg1g1prime2_pi2_VAJHU_ratio_NEW;
  float pg1g2_VAJHU_ratio_NEW;
  float pg1g2_pi2_VAJHU_ratio_NEW;
  float pg1g4_VAJHU_ratio_NEW;
  float pg1g4_pi2_VAJHU_ratio_NEW;

  float p0plus_VAMCFM, ggzz_p0plus_VAMCFM;
  float bkg_VAMCFM, ggzz_VAMCFM;
  float p0plus_VAMCFM_NEW, ggzz_p0plus_VAMCFM_NEW, p0plus_VAMCFM_NEW_BSMOn;
  float bkg_VAMCFM_NEW, ggzz_VAMCFM_NEW;
  //float bkg_VAMCFM_STU, bkg_VAMCFM_TU, bkg_VAMCFM_S;

  float mzz = 126.;
  float m1 = 91.471450;
  float m2 = 12.139782;
  float h1 = 0.2682896;
  float h2 = 0.1679779;
  float phi = 1.5969792;
  float hs = -0.727181;
  float phi1 = 1.8828257;

  TTree* tree = (TTree*)finput->Get(TREE_NAME);
  tree->SetBranchAddress("ZZMass", &mzz);
  tree->SetBranchAddress("Z1Mass", &m1);
  tree->SetBranchAddress("Z2Mass", &m2);
  tree->SetBranchAddress("helcosthetaZ1", &h1);
  tree->SetBranchAddress("helcosthetaZ2", &h2);
  tree->SetBranchAddress("helphi", &phi);
  tree->SetBranchAddress("costhetastar", &hs);
  tree->SetBranchAddress("phistarZ1", &phi1);
  tree->SetBranchAddress("p0plus_VAJHU", &p0plus_VAJHU);
  tree->SetBranchAddress("p0hplus_VAJHU", &p0hplus_VAJHU);
  tree->SetBranchAddress("p0minus_VAJHU", &p0minus_VAJHU);
  tree->SetBranchAddress("p0_g1prime2_VAJHU", &p0_g1prime2_VAJHU);
  tree->SetBranchAddress("pg1g2_VAJHU", &pg1g2_VAJHU);
  tree->SetBranchAddress("pg1g2_pi2_VAJHU", &pg1g2_pi2_VAJHU);
  tree->SetBranchAddress("pg1g4_VAJHU", &pg1g4_VAJHU);
  tree->SetBranchAddress("pg1g4_pi2_VAJHU", &pg1g4_pi2_VAJHU);
  tree->SetBranchAddress("pg1g1prime2_VAJHU", &pg1g1prime2_VAJHU);
  tree->SetBranchAddress("bkg_VAMCFM", &bkg_VAMCFM);
  tree->SetBranchAddress("ggzz_VAMCFM", &ggzz_VAMCFM);
  tree->SetBranchAddress("p0plus_VAMCFM", &p0plus_VAMCFM);
  tree->SetBranchAddress("ggzz_p0plus_VAMCFM", &ggzz_p0plus_VAMCFM);

  TTree* newtree = new TTree("TestTree", "");
  newtree->Branch("ZZMass", &mzz);
  newtree->Branch("Z1Mass", &m1);
  newtree->Branch("Z2Mass", &m2);
  newtree->Branch("helcosthetaZ1", &h1);
  newtree->Branch("helcosthetaZ2", &h2);
  newtree->Branch("helphi", &phi);
  newtree->Branch("costhetastar", &hs);
  newtree->Branch("phistarZ1", &phi1);

  newtree->Branch("p0plus_VAJHU", &p0plus_VAJHU);
  newtree->Branch("p0hplus_VAJHU", &p0hplus_VAJHU);
  newtree->Branch("p0minus_VAJHU", &p0minus_VAJHU);
  newtree->Branch("p0_g1prime2_VAJHU", &p0_g1prime2_VAJHU);
  newtree->Branch("pg1g2_VAJHU", &pg1g2_VAJHU);
  newtree->Branch("pg1g2_pi2_VAJHU", &pg1g2_pi2_VAJHU);
  newtree->Branch("pg1g4_VAJHU", &pg1g4_VAJHU);
  newtree->Branch("pg1g4_pi2_VAJHU", &pg1g4_pi2_VAJHU);
  newtree->Branch("pg1g1prime2_VAJHU", &pg1g1prime2_VAJHU);

  newtree->Branch("p0plus_VAJHU_NEW", &p0plus_VAJHU_NEW);
  newtree->Branch("p0hplus_VAJHU_NEW", &p0hplus_VAJHU_NEW);
  newtree->Branch("p0minus_VAJHU_NEW", &p0minus_VAJHU_NEW);
  newtree->Branch("p0_g1prime2_VAJHU_NEW", &p0_g1prime2_VAJHU_NEW);
  newtree->Branch("pg1g2_VAJHU_NEW", &pg1g2_VAJHU_NEW);
  newtree->Branch("pg1g2_pi2_VAJHU_NEW", &pg1g2_pi2_VAJHU_NEW);
  newtree->Branch("pg1g4_VAJHU_NEW", &pg1g4_VAJHU_NEW);
  newtree->Branch("pg1g4_pi2_VAJHU_NEW", &pg1g4_pi2_VAJHU_NEW);
  newtree->Branch("pg1g1prime2_VAJHU_NEW", &pg1g1prime2_VAJHU_NEW);
  newtree->Branch("pg1g1prime2_pi2_VAJHU_NEW", &pg1g1prime2_pi2_VAJHU_NEW);


  newtree->Branch("bkg_VAMCFM", &bkg_VAMCFM);
  newtree->Branch("ggzz_VAMCFM", &ggzz_VAMCFM);
  newtree->Branch("p0plus_VAMCFM", &p0plus_VAMCFM);
  newtree->Branch("ggzz_p0plus_VAMCFM", &ggzz_p0plus_VAMCFM);
  newtree->Branch("bkg_VAMCFM_NEW", &bkg_VAMCFM_NEW);
  newtree->Branch("ggzz_VAMCFM_NEW", &ggzz_VAMCFM_NEW);
  newtree->Branch("p0plus_VAMCFM_NEW", &p0plus_VAMCFM_NEW);
  newtree->Branch("ggzz_p0plus_VAMCFM_NEW", &ggzz_p0plus_VAMCFM_NEW);
  newtree->Branch("p0plus_VAMCFM_NEW_BSMOn", &p0plus_VAMCFM_NEW_BSMOn);

  newtree->Branch("p0hplus_VAMCFM_NEW", &p0hplus_VAMCFM_NEW);
  newtree->Branch("p0minus_VAMCFM_NEW", &p0minus_VAMCFM_NEW);
  newtree->Branch("p0_g1prime2_VAMCFM_NEW", &p0_g1prime2_VAMCFM_NEW);
  newtree->Branch("pg1g2_VAMCFM_NEW", &pg1g2_VAMCFM_NEW);
  newtree->Branch("pg1g2_pi2_VAMCFM_NEW", &pg1g2_pi2_VAMCFM_NEW);
  newtree->Branch("pg1g4_VAMCFM_NEW", &pg1g4_VAMCFM_NEW);
  newtree->Branch("pg1g4_pi2_VAMCFM_NEW", &pg1g4_pi2_VAMCFM_NEW);
  newtree->Branch("pg1g1prime2_VAMCFM_NEW", &pg1g1prime2_VAMCFM_NEW);
  newtree->Branch("pg1g1prime2_pi2_VAMCFM_NEW", &pg1g1prime2_pi2_VAMCFM_NEW);

  newtree->Branch("p0hplus_VAMCFM_ratio_NEW", &p0hplus_VAMCFM_ratio_NEW);
  newtree->Branch("p0minus_VAMCFM_ratio_NEW", &p0minus_VAMCFM_ratio_NEW);
  newtree->Branch("p0_g1prime2_VAMCFM_ratio_NEW", &p0_g1prime2_VAMCFM_ratio_NEW);
  newtree->Branch("pg1g2_VAMCFM_ratio_NEW", &pg1g2_VAMCFM_ratio_NEW);
  newtree->Branch("pg1g2_pi2_VAMCFM_ratio_NEW", &pg1g2_pi2_VAMCFM_ratio_NEW);
  newtree->Branch("pg1g4_VAMCFM_ratio_NEW", &pg1g4_VAMCFM_ratio_NEW);
  newtree->Branch("pg1g4_pi2_VAMCFM_ratio_NEW", &pg1g4_pi2_VAMCFM_ratio_NEW);
  newtree->Branch("pg1g1prime2_VAMCFM_ratio_NEW", &pg1g1prime2_VAMCFM_ratio_NEW);
  newtree->Branch("pg1g1prime2_pi2_VAMCFM_ratio_NEW", &pg1g1prime2_pi2_VAMCFM_ratio_NEW);

  newtree->Branch("p0hplus_VAJHU_ratio_NEW", &p0hplus_VAJHU_ratio_NEW);
  newtree->Branch("p0minus_VAJHU_ratio_NEW", &p0minus_VAJHU_ratio_NEW);
  newtree->Branch("p0_g1prime2_VAJHU_ratio_NEW", &p0_g1prime2_VAJHU_ratio_NEW);
  newtree->Branch("pg1g2_VAJHU_ratio_NEW", &pg1g2_VAJHU_ratio_NEW);
  newtree->Branch("pg1g2_pi2_VAJHU_ratio_NEW", &pg1g2_pi2_VAJHU_ratio_NEW);
  newtree->Branch("pg1g4_VAJHU_ratio_NEW", &pg1g4_VAJHU_ratio_NEW);
  newtree->Branch("pg1g4_pi2_VAJHU_ratio_NEW", &pg1g4_pi2_VAJHU_ratio_NEW);
  newtree->Branch("pg1g1prime2_VAJHU_ratio_NEW", &pg1g1prime2_VAJHU_ratio_NEW);
  newtree->Branch("pg1g1prime2_pi2_VAJHU_ratio_NEW", &pg1g1prime2_pi2_VAJHU_ratio_NEW);

  float GenLep1Id=0, GenLep2Id=0, GenLep3Id=0, GenLep4Id=0;
  if (flavor == 2){
    GenLep1Id=13;
    GenLep2Id=-13;
    GenLep3Id=11;
    GenLep4Id=-11;
  }
  else if (flavor == 1){
    GenLep1Id=11;
    GenLep2Id=-11;
    GenLep3Id=11;
    GenLep4Id=-11;
  }
  else if (flavor == 0){
    GenLep1Id=13;
    GenLep2Id=-13;
    GenLep3Id=13;
    GenLep4Id=-13;
  }
  else if (flavor == 3){
    GenLep1Id=14;
    GenLep2Id=-14;
    GenLep3Id=13;
    GenLep4Id=-13;
  }
  else if (flavor == 4){
    GenLep1Id=0;
    GenLep2Id=-0;
    GenLep3Id=1;
    GenLep4Id=-1;
  }
  else if (flavor == 5){
    GenLep1Id=1;
    GenLep2Id=-1;
    GenLep3Id=2;
    GenLep4Id=-2;
  }
  mela.setCandidateDecayMode(TVar::CandidateDecay_ZZ);

  for (int ev = 0; ev < min(1000, (int)tree->GetEntries()); ev++){
    tree->GetEntry(ev);

    int idOrdered[4] ={ static_cast<int>(GenLep1Id), static_cast<int>(GenLep2Id), static_cast<int>(GenLep3Id), static_cast<int>(GenLep4Id) };
    TLorentzVector pOrdered[4];
    std::vector<TLorentzVector> daus = mela.calculate4Momentum(mzz, m1, m2, acos(hs), acos(h1), acos(h2), phi1, phi);
    for (int ip=0; ip<min(4, (int)daus.size()); ip++) pOrdered[ip]=daus.at(ip);
    SimpleParticleCollection_t daughters_ZZ;
    for (unsigned int idau=0; idau<4; idau++) daughters_ZZ.push_back(SimpleParticle_t(idOrdered[idau], pOrdered[idau]));
    mela.setInputEvent(&daughters_ZZ, (SimpleParticleCollection_t*)0, (SimpleParticleCollection_t*)0, false);

    /***** ZZ *****/
    if (verbosity>=TVar::DEBUG) cout << "Computing bkg_VAMCFM_NEW" << endl;
    mela.setProcess(TVar::bkgZZ, TVar::MCFM, TVar::ZZQQB);
    mela.computeP(bkg_VAMCFM_NEW, false);

    if (verbosity>=TVar::DEBUG) cout << "Computing ggzz_VAMCFM_NEW" << endl;
    mela.setProcess(TVar::bkgZZ, TVar::MCFM, TVar::ZZGG);
    //mela.setMelaLeptonInterference(TVar::InterfOff);
    mela.computeP(ggzz_VAMCFM_NEW, false);

    if (verbosity>=TVar::DEBUG) cout << "Computing p0plus_VAMCFM_NEW" << endl;
    mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
    //mela.setMelaLeptonInterference(TVar::InterfOff);
    mela.computeP(p0plus_VAMCFM_NEW, false);

    if (verbosity>=TVar::DEBUG) cout << "Computing ggzz_p0plus_VAMCFM_NEW" << endl;
    mela.setProcess(TVar::bkgZZ_SMHiggs, TVar::MCFM, TVar::ZZGG);
    //mela.setMelaLeptonInterference(TVar::InterfOff);
    mela.computeP(ggzz_p0plus_VAMCFM_NEW, false);

    if (verbosity>=TVar::DEBUG) cout << "Computing p0plus_VAMCFM_NEW_BSMOn" << endl;
    mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.selfDHzzcoupl[0][0][0]=1;
    mela.computeP(p0plus_VAMCFM_NEW_BSMOn, false);


    mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);

    mela.selfDHggcoupl[0][0]=1;
    mela.selfDHzzcoupl[0][0][0]=1;
    mela.computeP(p0plus_VAJHU_NEW, false);

    mela.selfDHggcoupl[0][0]=1;
    mela.selfDHzzcoupl[0][1][0]=1.638;
    mela.computeP(p0hplus_VAJHU_NEW, false);

    mela.selfDHggcoupl[0][0]=1;
    mela.selfDHzzcoupl[0][3][0]=2.521;
    mela.computeP(p0minus_VAJHU_NEW, false);

    mela.selfDHggcoupl[0][0]=1;
    mela.selfDHzzcoupl[0][11][0]=-12046.01;
    mela.computeP(p0_g1prime2_VAJHU_NEW, false);

    mela.selfDHggcoupl[0][0]=1;
    mela.selfDHzzcoupl[0][0][0]=1;
    mela.selfDHzzcoupl[0][1][0]=1.638;
    mela.computeP(pg1g2_VAJHU_NEW, false);
    pg1g2_VAJHU_NEW -= (p0plus_VAJHU_NEW + p0hplus_VAJHU_NEW);

    mela.selfDHggcoupl[0][0]=1;
    mela.selfDHzzcoupl[0][0][0]=1;
    mela.selfDHzzcoupl[0][1][1]=1.638;
    mela.computeP(pg1g2_pi2_VAJHU_NEW, false);
    pg1g2_pi2_VAJHU_NEW -= (p0plus_VAJHU_NEW + p0hplus_VAJHU_NEW);

    mela.selfDHggcoupl[0][0]=1;
    mela.selfDHzzcoupl[0][0][0]=1;
    mela.selfDHzzcoupl[0][3][0]=2.521;
    mela.computeP(pg1g4_VAJHU_NEW, false);
    pg1g4_VAJHU_NEW -= (p0plus_VAJHU_NEW + p0minus_VAJHU_NEW);

    mela.selfDHggcoupl[0][0]=1;
    mela.selfDHzzcoupl[0][0][0]=1;
    mela.selfDHzzcoupl[0][3][1]=2.521;
    mela.computeP(pg1g4_pi2_VAJHU_NEW, false);
    pg1g4_pi2_VAJHU_NEW -= (p0plus_VAJHU_NEW + p0minus_VAJHU_NEW);

    mela.selfDHggcoupl[0][0]=1;
    mela.selfDHzzcoupl[0][0][0]=1;
    mela.selfDHzzcoupl[0][11][0]=12046.01;
    mela.computeP(pg1g1prime2_VAJHU_NEW, false);
    pg1g1prime2_VAJHU_NEW -= (p0plus_VAJHU_NEW + p0_g1prime2_VAJHU_NEW);

    mela.selfDHggcoupl[0][0]=1;
    mela.selfDHzzcoupl[0][0][0]=1;
    mela.selfDHzzcoupl[0][11][1]=12046.01;
    mela.computeP(pg1g1prime2_pi2_VAJHU_NEW, false);
    pg1g1prime2_pi2_VAJHU_NEW -= (p0plus_VAJHU_NEW + p0_g1prime2_VAJHU_NEW);


    mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);

    mela.selfDHzzcoupl[0][0][0]=1;
    mela.computeP(p0plus_VAMCFM_NEW, false);

    mela.selfDHzzcoupl[0][1][0]=1.638;
    mela.computeP(p0hplus_VAMCFM_NEW, false);

    mela.selfDHzzcoupl[0][3][0]=2.521;
    mela.computeP(p0minus_VAMCFM_NEW, false);

    mela.selfDHzzcoupl[0][11][0]=-12046.01;
    mela.computeP(p0_g1prime2_VAMCFM_NEW, false);

    mela.selfDHzzcoupl[0][0][0]=1;
    mela.selfDHzzcoupl[0][1][0]=1.638;
    mela.computeP(pg1g2_VAMCFM_NEW, false);
    pg1g2_VAMCFM_NEW -= (p0plus_VAMCFM_NEW + p0hplus_VAMCFM_NEW);

    mela.selfDHzzcoupl[0][0][0]=1;
    mela.selfDHzzcoupl[0][1][1]=1.638;
    mela.computeP(pg1g2_pi2_VAMCFM_NEW, false);
    pg1g2_pi2_VAMCFM_NEW -= (p0plus_VAMCFM_NEW + p0hplus_VAMCFM_NEW);

    mela.selfDHzzcoupl[0][0][0]=1;
    mela.selfDHzzcoupl[0][3][0]=2.521;
    mela.computeP(pg1g4_VAMCFM_NEW, false);
    pg1g4_VAMCFM_NEW -= (p0plus_VAMCFM_NEW + p0minus_VAMCFM_NEW);

    mela.selfDHzzcoupl[0][0][0]=1;
    mela.selfDHzzcoupl[0][3][1]=2.521;
    mela.computeP(pg1g4_pi2_VAMCFM_NEW, false);
    pg1g4_pi2_VAMCFM_NEW -= (p0plus_VAMCFM_NEW + p0minus_VAMCFM_NEW);

    mela.selfDHzzcoupl[0][0][0]=1;
    mela.selfDHzzcoupl[0][11][0]=12046.01;
    mela.computeP(pg1g1prime2_VAMCFM_NEW, false);
    pg1g1prime2_VAMCFM_NEW -= (p0plus_VAMCFM_NEW + p0_g1prime2_VAMCFM_NEW);

    mela.selfDHzzcoupl[0][0][0]=1;
    mela.selfDHzzcoupl[0][11][1]=12046.01;
    mela.computeP(pg1g1prime2_pi2_VAMCFM_NEW, false);
    pg1g1prime2_pi2_VAMCFM_NEW -= (p0plus_VAMCFM_NEW + p0_g1prime2_VAMCFM_NEW);


    p0hplus_VAJHU_ratio_NEW = p0hplus_VAJHU_NEW/p0plus_VAJHU_NEW;
    p0minus_VAJHU_ratio_NEW = p0minus_VAJHU_NEW/p0plus_VAJHU_NEW;
    p0_g1prime2_VAJHU_ratio_NEW = p0_g1prime2_VAJHU_NEW/p0plus_VAJHU_NEW;
    pg1g2_VAJHU_ratio_NEW = pg1g2_VAJHU_NEW/p0plus_VAJHU_NEW;
    pg1g4_VAJHU_ratio_NEW = pg1g4_VAJHU_NEW/p0plus_VAJHU_NEW;
    pg1g1prime2_VAJHU_ratio_NEW = pg1g1prime2_VAJHU_NEW/p0plus_VAJHU_NEW;
    pg1g2_pi2_VAJHU_ratio_NEW = pg1g2_pi2_VAJHU_NEW/p0plus_VAJHU_NEW;
    pg1g4_pi2_VAJHU_ratio_NEW = pg1g4_pi2_VAJHU_NEW/p0plus_VAJHU_NEW;
    pg1g1prime2_pi2_VAJHU_ratio_NEW = pg1g1prime2_pi2_VAJHU_NEW/p0plus_VAJHU_NEW;

    p0hplus_VAMCFM_ratio_NEW = p0hplus_VAMCFM_NEW/p0plus_VAMCFM_NEW;
    p0minus_VAMCFM_ratio_NEW = p0minus_VAMCFM_NEW/p0plus_VAMCFM_NEW;
    p0_g1prime2_VAMCFM_ratio_NEW = p0_g1prime2_VAMCFM_NEW/p0plus_VAMCFM_NEW;
    pg1g2_VAMCFM_ratio_NEW = pg1g2_VAMCFM_NEW/p0plus_VAMCFM_NEW;
    pg1g4_VAMCFM_ratio_NEW = pg1g4_VAMCFM_NEW/p0plus_VAMCFM_NEW;
    pg1g1prime2_VAMCFM_ratio_NEW = pg1g1prime2_VAMCFM_NEW/p0plus_VAMCFM_NEW;
    pg1g2_pi2_VAMCFM_ratio_NEW = pg1g2_pi2_VAMCFM_NEW/p0plus_VAMCFM_NEW;
    pg1g4_pi2_VAMCFM_ratio_NEW = pg1g4_pi2_VAMCFM_NEW/p0plus_VAMCFM_NEW;
    pg1g1prime2_pi2_VAMCFM_ratio_NEW = pg1g1prime2_pi2_VAMCFM_NEW/p0plus_VAMCFM_NEW;


    mela.resetInputEvent();
    newtree->Fill();
  }


  foutput->WriteTObject(newtree);
  foutput->Close();
  finput->Close();
}

void testME_ProdP_VBFHJJ_FullSim(int flavor=2){
  int erg_tev=8;
  float mPOLE=125.6;
  float wPOLE=4.07e-3;
  TString TREE_NAME = "SelectedTree";

  TVar::VerbosityLevel verbosity = TVar::DEBUG_MECHECK;
  Mela mela(erg_tev, mPOLE, verbosity);

  TString cinput_main = "/scratch0/hep/ianderso/CJLST/140519/PRODFSR_8TeV";
  TFile* finput = new TFile(Form("%s/%s/HZZ4lTree_ZZTo%s.root", cinput_main.Data(), (flavor>=2 ? "2mu2e" : "4e"), (flavor==2 ? "2e2mu" : "4e")), "read");
  TFile* foutput = new TFile(Form("HZZ4lTree_ZZTo%s_vbfMELATest.root", (flavor>=2 ? "2e2mu" : "4e")), "recreate");

  float phjj_VAJHU_old;
  float pvbf_VAJHU_old;
  float phjj_VAJHU_old_NEW;
  float pvbf_VAJHU_old_NEW;
  float phjj_VAJHU_old_NEW_selfD;
  float pvbf_VAJHU_old_NEW_selfD;
  float phjj0minus_VAJHU_old_NEW;
  float pvbf0minus_VAJHU_old_NEW;
  float phjj0minus_VAJHU_old_NEW_selfD;
  float pvbf0minus_VAJHU_old_NEW_selfD;

  //float jet1Pt, jet2Pt;
  //float jet1px, jet1py, jet1pz, jet1E;
  //float jet2px, jet2py, jet2pz, jet2E;
  //float ZZPx, ZZPy, ZZPz, ZZE, dR;
  short NJets30;
  std::vector<double>* JetPt=0;
  std::vector<double>* JetEta=0;
  std::vector<double>* JetPhi=0;
  std::vector<double>* JetMass=0;
  std::vector<double> myJetPt;
  std::vector<double> myJetEta;
  std::vector<double> myJetPhi;
  std::vector<double> myJetMass;
  TBranch* bJetPt=0;
  TBranch* bJetEta=0;
  TBranch* bJetPhi=0;
  TBranch* bJetMass=0;

  float mzz = 126.;
  float m1 = 91.471450;
  float m2 = 12.139782;
  float h1 = 0.2682896;
  float h2 = 0.1679779;
  float phi = 1.5969792;
  float hs = -0.727181;
  float phi1 = 1.8828257;
  float ZZPt, ZZPhi, ZZEta;

  TTree* tree = (TTree*)finput->Get(TREE_NAME);
  tree->SetBranchAddress("NJets30", &NJets30);
  tree->SetBranchAddress("JetPt", &JetPt, &bJetPt);
  tree->SetBranchAddress("JetEta", &JetEta, &bJetEta);
  tree->SetBranchAddress("JetPhi", &JetPhi, &bJetPhi);
  tree->SetBranchAddress("JetMass", &JetMass, &bJetMass);
  tree->SetBranchAddress("phjj_VAJHU_old", &phjj_VAJHU_old);
  tree->SetBranchAddress("pvbf_VAJHU_old", &pvbf_VAJHU_old);
  tree->SetBranchAddress("ZZMass", &mzz);
  tree->SetBranchAddress("ZZPt", &ZZPt);
  tree->SetBranchAddress("ZZEta", &ZZEta);
  tree->SetBranchAddress("ZZPhi", &ZZPhi);
  tree->SetBranchAddress("Z1Mass", &m1);
  tree->SetBranchAddress("Z2Mass", &m2);
  tree->SetBranchAddress("helcosthetaZ1", &h1);
  tree->SetBranchAddress("helcosthetaZ2", &h2);
  tree->SetBranchAddress("helphi", &phi);
  tree->SetBranchAddress("costhetastar", &hs);
  tree->SetBranchAddress("phistarZ1", &phi1);

  TTree* newtree = new TTree("TestTree", "");
  newtree->Branch("phjj_VAJHU_old", &phjj_VAJHU_old);
  newtree->Branch("pvbf_VAJHU_old", &pvbf_VAJHU_old);
  newtree->Branch("phjj_VAJHU_old_NEW", &phjj_VAJHU_old_NEW);
  newtree->Branch("pvbf_VAJHU_old_NEW", &pvbf_VAJHU_old_NEW);
  newtree->Branch("phjj0minus_VAJHU_old_NEW", &phjj0minus_VAJHU_old_NEW);
  newtree->Branch("pvbf0minus_VAJHU_old_NEW", &pvbf0minus_VAJHU_old_NEW);
  newtree->Branch("phjj_VAJHU_old_NEW_selfD", &phjj_VAJHU_old_NEW_selfD);
  newtree->Branch("pvbf_VAJHU_old_NEW_selfD", &pvbf_VAJHU_old_NEW_selfD);
  newtree->Branch("phjj0minus_VAJHU_old_NEW_selfD", &phjj0minus_VAJHU_old_NEW_selfD);
  newtree->Branch("pvbf0minus_VAJHU_old_NEW_selfD", &pvbf0minus_VAJHU_old_NEW_selfD);
  newtree->Branch("ZZMass", &mzz);

  float GenLep1Id=0, GenLep2Id=0, GenLep3Id=0, GenLep4Id=0;
  if (flavor == 2){
    GenLep1Id=13;
    GenLep2Id=-13;
    GenLep3Id=11;
    GenLep4Id=-11;
  }
  else if (flavor == 1){
    GenLep1Id=11;
    GenLep2Id=-11;
    GenLep3Id=11;
    GenLep4Id=-11;
  }
  else if (flavor == 0){
    GenLep1Id=13;
    GenLep2Id=-13;
    GenLep3Id=13;
    GenLep4Id=-13;
  }
  else if (flavor == 3){
    GenLep1Id=14;
    GenLep2Id=-14;
    GenLep3Id=13;
    GenLep4Id=-13;
  }
  else if (flavor == 4){
    GenLep1Id=0;
    GenLep2Id=-0;
    GenLep3Id=1;
    GenLep4Id=-1;
  }
  else if (flavor == 5){
    GenLep1Id=1;
    GenLep2Id=-1;
    GenLep3Id=2;
    GenLep4Id=-2;
  }
  mela.setCandidateDecayMode(TVar::CandidateDecay_ZZ);
  int idOrdered[4] ={ static_cast<int>(GenLep1Id), static_cast<int>(GenLep2Id), static_cast<int>(GenLep3Id), static_cast<int>(GenLep4Id) };

  int nEntries = tree->GetEntries();
  int recorded=0;
  for (int ev = 0; ev < nEntries; ev++){
    if (recorded>=(verbosity>=TVar::DEBUG ? 1 : 1000)) break;
    tree->GetEntry(ev);

    if (JetPt->size()>=2 && NJets30>=2){
      TLorentzVector jet1(0, 0, 1e-3, 1e-3), jet2(0, 0, 1e-3, 1e-3), higgs(0, 0, 0, 0);
      jet1.SetPtEtaPhiM(JetPt->at(0), JetEta->at(0), JetPhi->at(0), JetMass->at(0));
      jet2.SetPtEtaPhiM(JetPt->at(1), JetEta->at(1), JetPhi->at(1), JetMass->at(1));
      higgs.SetPtEtaPhiM(ZZPt, ZZEta, ZZPhi, mzz);
      TVector3 boostH = higgs.BoostVector();

      SimpleParticleCollection_t associated;
      associated.push_back(SimpleParticle_t(0, jet1));
      associated.push_back(SimpleParticle_t(0, jet2));

      TLorentzVector pDaughters[4];
      std::vector<TLorentzVector> daus = mela.calculate4Momentum(mzz, m1, m2, acos(hs), acos(h1), acos(h2), phi1, phi);
      for (int ip=0; ip<min(4, (int)daus.size()); ip++){ pDaughters[ip]=daus.at(ip); pDaughters[ip].Boost(boostH); }
      SimpleParticleCollection_t daughters_ZZ;
      for (unsigned int idau=0; idau<4; idau++) daughters_ZZ.push_back(SimpleParticle_t(idOrdered[idau], pDaughters[idau]));
      mela.setInputEvent(&daughters_ZZ, &associated, (SimpleParticleCollection_t*)0, false);


      mela.setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::JJVBF);
      mela.computeProdP(pvbf_VAJHU_old_NEW, false);
/*
      mela.setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::JJQCD);
      mela.computeProdP(phjj_VAJHU_old_NEW, false);

      mela.setProcess(TVar::H0minus, TVar::JHUGen, TVar::JJVBF);
      mela.computeProdP(pvbf0minus_VAJHU_old_NEW, false);

      mela.setProcess(TVar::H0minus, TVar::JHUGen, TVar::JJQCD);
      mela.computeProdP(phjj0minus_VAJHU_old_NEW, false);


      mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::JJVBF);
      mela.selfDHzzcoupl[0][0][0]=1;
      mela.computeProdP(pvbf_VAJHU_old_NEW_selfD, false);

      mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::JJQCD);
      mela.selfDHggcoupl[0][0]=1;
      mela.computeProdP(phjj_VAJHU_old_NEW_selfD, false);


      mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::JJVBF);
      mela.selfDHzzcoupl[0][3][0]=1;
      mela.computeProdP(pvbf0minus_VAJHU_old_NEW_selfD, false);

      mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::JJQCD);
      mela.selfDHggcoupl[2][0]=1;
      mela.computeProdP(phjj0minus_VAJHU_old_NEW_selfD, false);
      */
      newtree->Fill();
      recorded++;
      mela.resetInputEvent();
    }
  }


  foutput->WriteTObject(newtree);
  delete newtree;
  foutput->Close();
  finput->Close();
}

void testME_ProdP_VH_FullSim(){
  int erg_tev=8;
  float mPOLE=125.6;
  float wPOLE=4.07e-3;
  TString TREE_NAME = "SelectedTree";

  TVar::VerbosityLevel verbosity = TVar::DEBUG;
  Mela mela(erg_tev, mPOLE, verbosity);

  TString cinput_main = "/scratch0/hep/ianderso/CJLST/140519/PRODFSR_8TeV";
  TFile* finput = new TFile(Form("%s/%s/HZZ4lTree_ZZTo%s.root", cinput_main.Data(), "2mu2e", "2e2mu"), "read");
  TFile* foutput = new TFile(Form("HZZ4lTree_ZZTo%s_vhMELATest.root", "2e2mu"), "recreate");

  float pwh_leptonic_VAJHU_old_NEW;
  float pzh_leptonic_VAJHU_old_NEW;
  float pwh_leptonic_VAJHU_old_NEW_selfD;
  float pzh_leptonic_VAJHU_old_NEW_selfD;
  float pwh0minus_leptonic_VAJHU_old_NEW;
  float pzh0minus_leptonic_VAJHU_old_NEW;
  float pwh0minus_leptonic_VAJHU_old_NEW_selfD;
  float pzh0minus_leptonic_VAJHU_old_NEW_selfD;
  float pwh_hadronic_VAJHU_old_NEW;
  float pzh_hadronic_VAJHU_old_NEW;
  float pwh_hadronic_VAJHU_old_NEW_selfD;
  float pzh_hadronic_VAJHU_old_NEW_selfD;
  float pwh0minus_hadronic_VAJHU_old_NEW;
  float pzh0minus_hadronic_VAJHU_old_NEW;
  float pwh0minus_hadronic_VAJHU_old_NEW_selfD;
  float pzh0minus_hadronic_VAJHU_old_NEW_selfD;

  //float jet1Pt, jet2Pt;
  //float jet1px, jet1py, jet1pz, jet1E;
  //float jet2px, jet2py, jet2pz, jet2E;
  //float ZZPx, ZZPy, ZZPz, ZZE, dR;
  short NJets30;
  std::vector<double>* JetPt=0;
  std::vector<double>* JetEta=0;
  std::vector<double>* JetPhi=0;
  std::vector<double>* JetMass=0;
  std::vector<double> myJetPt;
  std::vector<double> myJetEta;
  std::vector<double> myJetPhi;
  std::vector<double> myJetMass;
  TBranch* bJetPt=0;
  TBranch* bJetEta=0;
  TBranch* bJetPhi=0;
  TBranch* bJetMass=0;

  float mzz = 126.;
  float m1 = 91.471450;
  float m2 = 12.139782;
  float h1 = 0.2682896;
  float h2 = 0.1679779;
  float phi = 1.5969792;
  float hs = -0.727181;
  float phi1 = 1.8828257;
  float ZZPt, ZZPhi, ZZEta;

  TTree* tree = (TTree*)finput->Get(TREE_NAME);
  tree->SetBranchAddress("NJets30", &NJets30);
  tree->SetBranchAddress("JetPt", &JetPt, &bJetPt);
  tree->SetBranchAddress("JetEta", &JetEta, &bJetEta);
  tree->SetBranchAddress("JetPhi", &JetPhi, &bJetPhi);
  tree->SetBranchAddress("JetMass", &JetMass, &bJetMass);
  tree->SetBranchAddress("ZZMass", &mzz);
  tree->SetBranchAddress("ZZPt", &ZZPt);
  tree->SetBranchAddress("ZZEta", &ZZEta);
  tree->SetBranchAddress("ZZPhi", &ZZPhi);
  tree->SetBranchAddress("Z1Mass", &m1);
  tree->SetBranchAddress("Z2Mass", &m2);
  tree->SetBranchAddress("helcosthetaZ1", &h1);
  tree->SetBranchAddress("helcosthetaZ2", &h2);
  tree->SetBranchAddress("helphi", &phi);
  tree->SetBranchAddress("costhetastar", &hs);
  tree->SetBranchAddress("phistarZ1", &phi1);

  TTree* newtree = new TTree("TestTree", "");
  newtree->Branch("pwh_leptonic_VAJHU_old_NEW", &pwh_leptonic_VAJHU_old_NEW);
  newtree->Branch("pzh_leptonic_VAJHU_old_NEW", &pzh_leptonic_VAJHU_old_NEW);
  newtree->Branch("pwh0minus_leptonic_VAJHU_old_NEW", &pwh0minus_leptonic_VAJHU_old_NEW);
  newtree->Branch("pzh0minus_leptonic_VAJHU_old_NEW", &pzh0minus_leptonic_VAJHU_old_NEW);
  newtree->Branch("pwh_leptonic_VAJHU_old_NEW_selfD", &pwh_leptonic_VAJHU_old_NEW_selfD);
  newtree->Branch("pzh_leptonic_VAJHU_old_NEW_selfD", &pzh_leptonic_VAJHU_old_NEW_selfD);
  newtree->Branch("pwh0minus_leptonic_VAJHU_old_NEW_selfD", &pwh0minus_leptonic_VAJHU_old_NEW_selfD);
  newtree->Branch("pzh0minus_leptonic_VAJHU_old_NEW_selfD", &pzh0minus_leptonic_VAJHU_old_NEW_selfD);
  newtree->Branch("pwh_hadronic_VAJHU_old_NEW", &pwh_hadronic_VAJHU_old_NEW);
  newtree->Branch("pzh_hadronic_VAJHU_old_NEW", &pzh_hadronic_VAJHU_old_NEW);
  newtree->Branch("pwh0minus_hadronic_VAJHU_old_NEW", &pwh0minus_hadronic_VAJHU_old_NEW);
  newtree->Branch("pzh0minus_hadronic_VAJHU_old_NEW", &pzh0minus_hadronic_VAJHU_old_NEW);
  newtree->Branch("pwh_hadronic_VAJHU_old_NEW_selfD", &pwh_hadronic_VAJHU_old_NEW_selfD);
  newtree->Branch("pzh_hadronic_VAJHU_old_NEW_selfD", &pzh_hadronic_VAJHU_old_NEW_selfD);
  newtree->Branch("pwh0minus_hadronic_VAJHU_old_NEW_selfD", &pwh0minus_hadronic_VAJHU_old_NEW_selfD);
  newtree->Branch("pzh0minus_hadronic_VAJHU_old_NEW_selfD", &pzh0minus_hadronic_VAJHU_old_NEW_selfD);
  newtree->Branch("ZZMass", &mzz);

  float GenLep1Id=0, GenLep2Id=0, GenLep3Id=0, GenLep4Id=0;
    GenLep1Id=13;
    GenLep2Id=-13;
    GenLep3Id=14;
    GenLep4Id=-14;
  mela.setCandidateDecayMode(TVar::CandidateDecay_ZZ);
  int idOrdered[4] ={ static_cast<int>(GenLep1Id), static_cast<int>(GenLep2Id), static_cast<int>(GenLep3Id), static_cast<int>(GenLep4Id) };

  int nEntries = tree->GetEntries();
  int recorded=0;
  for (int ev = 0; ev < nEntries; ev++){
    if (recorded>=(verbosity>=TVar::DEBUG ? 1 : 1000)) break;
    tree->GetEntry(ev);

    if (JetPt->size()>=2 && NJets30>=2){
      TLorentzVector jet1(0, 0, 1e-3, 1e-3), jet2(0, 0, 1e-3, 1e-3), higgs(0, 0, 0, 0);
      jet1.SetPtEtaPhiM(JetPt->at(0), JetEta->at(0), JetPhi->at(0), JetMass->at(0));
      jet2.SetPtEtaPhiM(JetPt->at(1), JetEta->at(1), JetPhi->at(1), JetMass->at(1));
      higgs.SetPtEtaPhiM(ZZPt, ZZEta, ZZPhi, mzz);
      TVector3 boostH = higgs.BoostVector();

      TLorentzVector pDaughters[4];
      std::vector<TLorentzVector> daus = mela.calculate4Momentum(mzz, m1, m2, acos(hs), acos(h1), acos(h2), phi1, phi);
      for (int ip=0; ip<min(4, (int)daus.size()); ip++){ pDaughters[ip]=daus.at(ip); pDaughters[ip].Boost(boostH); }

      SimpleParticleCollection_t associated;
      associated.push_back(SimpleParticle_t(0, jet1));
      associated.push_back(SimpleParticle_t(0, jet2));
      for (unsigned int idau=0; idau<4; idau++) associated.push_back(SimpleParticle_t(idOrdered[idau], pDaughters[idau])); // Let's put some extra leptons

      SimpleParticleCollection_t daughters_ZZ;
      for (unsigned int idau=0; idau<4; idau++) daughters_ZZ.push_back(SimpleParticle_t(idOrdered[idau], pDaughters[idau]));
      mela.setInputEvent(&daughters_ZZ, &associated, (SimpleParticleCollection_t*)0, false);

      if (verbosity>=TVar::DEBUG){
        cout << "Mela candidates summary:" << endl;
        for (int ic=0; ic<mela.getNCandidates(); ic++){
          cout << "Candidate " << ic << endl;
          mela.setCurrentCandidateFromIndex(ic);
          TUtil::PrintCandidateSummary(mela.getCurrentCandidate());
        }
        cout << endl;
      }

      mela.setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::Lep_ZH);
      mela.computeProdP_VH(pzh_leptonic_VAJHU_old_NEW, false, false);

      mela.setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::Lep_WH);
      mela.computeProdP_VH(pwh_leptonic_VAJHU_old_NEW, false, false);

      mela.setProcess(TVar::H0minus, TVar::JHUGen, TVar::Lep_ZH);
      mela.computeProdP_VH(pzh0minus_leptonic_VAJHU_old_NEW, false, false);

      mela.setProcess(TVar::H0minus, TVar::JHUGen, TVar::Lep_WH);
      mela.computeProdP_VH(pwh0minus_leptonic_VAJHU_old_NEW, false, false);


      mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::Lep_ZH);
      mela.selfDHzzcoupl[0][0][0]=1;
      mela.computeProdP_VH(pzh_leptonic_VAJHU_old_NEW_selfD, false, false);

      mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::Lep_WH);
      mela.selfDHzzcoupl[0][0]=1;
      mela.computeProdP_VH(pwh_leptonic_VAJHU_old_NEW_selfD, false, false);


      mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::Lep_ZH);
      mela.selfDHzzcoupl[0][3][0]=1;
      mela.computeProdP_VH(pzh0minus_leptonic_VAJHU_old_NEW_selfD, false, false);

      mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::Lep_WH);
      mela.selfDHzzcoupl[2][0]=1;
      mela.computeProdP_VH(pwh0minus_leptonic_VAJHU_old_NEW_selfD, false, false);

      mela.setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::Had_ZH);
      mela.computeProdP_VH(pzh_hadronic_VAJHU_old_NEW, false, false);

      mela.setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::Had_WH);
      mela.computeProdP_VH(pwh_hadronic_VAJHU_old_NEW, false, false);

      mela.setProcess(TVar::H0minus, TVar::JHUGen, TVar::Had_ZH);
      mela.computeProdP_VH(pzh0minus_hadronic_VAJHU_old_NEW, false, false);

      mela.setProcess(TVar::H0minus, TVar::JHUGen, TVar::Had_WH);
      mela.computeProdP_VH(pwh0minus_hadronic_VAJHU_old_NEW, false, false);


      mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::Had_ZH);
      mela.selfDHzzcoupl[0][0][0]=1;
      mela.computeProdP_VH(pzh_hadronic_VAJHU_old_NEW_selfD, false, false);

      mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::Had_WH);
      mela.selfDHzzcoupl[0][0]=1;
      mela.computeProdP_VH(pwh_hadronic_VAJHU_old_NEW_selfD, false, false);


      mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::Had_ZH);
      mela.selfDHzzcoupl[0][3][0]=1;
      mela.computeProdP_VH(pzh0minus_hadronic_VAJHU_old_NEW_selfD, false, false);

      mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::Had_WH);
      mela.selfDHzzcoupl[2][0]=1;
      mela.computeProdP_VH(pwh0minus_hadronic_VAJHU_old_NEW_selfD, false, false);

      newtree->Fill();
      recorded++;
      mela.resetInputEvent();
    }
  }


  foutput->WriteTObject(newtree);
  delete newtree;
  foutput->Close();
  finput->Close();
}

void testME_ProdP_TTHBBH_FullSim(){
  int erg_tev=8;
  float mPOLE=125.6;
  float wPOLE=4.07e-3;
  TString TREE_NAME = "SelectedTree";

  TVar::VerbosityLevel verbosity = TVar::ERROR;
  Mela mela(erg_tev, mPOLE, verbosity);

  TString cinput_main = "/scratch0/hep/ianderso/CJLST/140519/PRODFSR_8TeV";
  TFile* finput = new TFile(Form("%s/%s/HZZ4lTree_ZZTo%s.root", cinput_main.Data(), "2mu2e", "2e2mu"), "read");
  TFile* foutput = new TFile(Form("HZZ4lTree_ZZTo%s_ttHbbHMELATest.root", "2e2mu"), "recreate");

  float pbbh_VAJHU_old_NEW;
  float ptth_VAJHU_old_NEW;
  float pbbh_VAJHU_old_NEW_selfD;
  float ptth_VAJHU_old_NEW_selfD;
  float pbbh0minus_VAJHU_old_NEW;
  float ptth0minus_VAJHU_old_NEW;
  float pbbh0minus_VAJHU_old_NEW_selfD;
  float ptth0minus_VAJHU_old_NEW_selfD;

  short NJets30;
  std::vector<double>* JetPt=0;
  std::vector<double>* JetEta=0;
  std::vector<double>* JetPhi=0;
  std::vector<double>* JetMass=0;
  std::vector<double> myJetPt;
  std::vector<double> myJetEta;
  std::vector<double> myJetPhi;
  std::vector<double> myJetMass;
  TBranch* bJetPt=0;
  TBranch* bJetEta=0;
  TBranch* bJetPhi=0;
  TBranch* bJetMass=0;

  float mzz = 126.;
  float m1 = 91.471450;
  float m2 = 12.139782;
  float h1 = 0.2682896;
  float h2 = 0.1679779;
  float phi = 1.5969792;
  float hs = -0.727181;
  float phi1 = 1.8828257;
  float ZZPt, ZZPhi, ZZEta;

  TTree* tree = (TTree*)finput->Get(TREE_NAME);
  tree->SetBranchAddress("NJets30", &NJets30);
  tree->SetBranchAddress("JetPt", &JetPt, &bJetPt);
  tree->SetBranchAddress("JetEta", &JetEta, &bJetEta);
  tree->SetBranchAddress("JetPhi", &JetPhi, &bJetPhi);
  tree->SetBranchAddress("JetMass", &JetMass, &bJetMass);
  tree->SetBranchAddress("ZZMass", &mzz);
  tree->SetBranchAddress("ZZPt", &ZZPt);
  tree->SetBranchAddress("ZZEta", &ZZEta);
  tree->SetBranchAddress("ZZPhi", &ZZPhi);
  tree->SetBranchAddress("Z1Mass", &m1);
  tree->SetBranchAddress("Z2Mass", &m2);
  tree->SetBranchAddress("helcosthetaZ1", &h1);
  tree->SetBranchAddress("helcosthetaZ2", &h2);
  tree->SetBranchAddress("helphi", &phi);
  tree->SetBranchAddress("costhetastar", &hs);
  tree->SetBranchAddress("phistarZ1", &phi1);

  TTree* newtree = new TTree("TestTree", "");
  newtree->Branch("pbbh_VAJHU_old_NEW", &pbbh_VAJHU_old_NEW);
  newtree->Branch("ptth_VAJHU_old_NEW", &ptth_VAJHU_old_NEW);
  newtree->Branch("pbbh0minus_VAJHU_old_NEW", &pbbh0minus_VAJHU_old_NEW);
  newtree->Branch("ptth0minus_VAJHU_old_NEW", &ptth0minus_VAJHU_old_NEW);
  newtree->Branch("pbbh_VAJHU_old_NEW_selfD", &pbbh_VAJHU_old_NEW_selfD);
  newtree->Branch("ptth_VAJHU_old_NEW_selfD", &ptth_VAJHU_old_NEW_selfD);
  newtree->Branch("pbbh0minus_VAJHU_old_NEW_selfD", &pbbh0minus_VAJHU_old_NEW_selfD);
  newtree->Branch("ptth0minus_VAJHU_old_NEW_selfD", &ptth0minus_VAJHU_old_NEW_selfD);
  newtree->Branch("ZZMass", &mzz);

  float GenLep1Id=0, GenLep2Id=0, GenLep3Id=0, GenLep4Id=0;
  GenLep1Id=13;
  GenLep2Id=-13;
  GenLep3Id=11;
  GenLep4Id=-11;
  mela.setCandidateDecayMode(TVar::CandidateDecay_ZZ);
  int idOrdered[4] ={ static_cast<int>(GenLep1Id), static_cast<int>(GenLep2Id), static_cast<int>(GenLep3Id), static_cast<int>(GenLep4Id) };

  int nEntries = tree->GetEntries();
  int recorded=0;
  for (int ev = 0; ev < nEntries; ev++){
    if (recorded>=1000) break;
    tree->GetEntry(ev);

    if (JetPt->size()>=2 && NJets30>=2){
      TLorentzVector jet1(0, 0, 1e-3, 1e-3), jet2(0, 0, 1e-3, 1e-3), higgs(0, 0, 0, 0);
      jet1.SetPtEtaPhiM(JetPt->at(0), JetEta->at(0), JetPhi->at(0), JetMass->at(0));
      jet2.SetPtEtaPhiM(JetPt->at(1), JetEta->at(1), JetPhi->at(1), JetMass->at(1));
      higgs.SetPtEtaPhiM(ZZPt, ZZEta, ZZPhi, mzz);
      TVector3 boostH = higgs.BoostVector();

      SimpleParticleCollection_t associated;
      associated.push_back(SimpleParticle_t(0, jet1));
      associated.push_back(SimpleParticle_t(0, jet2));

      TLorentzVector pDaughters[4];
      std::vector<TLorentzVector> daus = mela.calculate4Momentum(mzz, m1, m2, acos(hs), acos(h1), acos(h2), phi1, phi);
      for (int ip=0; ip<min(4, (int)daus.size()); ip++){ pDaughters[ip]=daus.at(ip); pDaughters[ip].Boost(boostH); }
      SimpleParticleCollection_t daughters;
      for (unsigned int idau=0; idau<4; idau++) daughters.push_back(SimpleParticle_t(idOrdered[idau], pDaughters[idau]));
      mela.setInputEvent(&daughters, &associated, (SimpleParticleCollection_t*)0, false);

      // This is some hack to get top daughter four-vectors
      SimpleParticleCollection_t topDaughters;
      topDaughters.push_back(SimpleParticle_t(0, jet1));
      for (unsigned int idau=0; idau<2; idau++) topDaughters.push_back(SimpleParticle_t(0, pDaughters[idau]));
      mela.appendTopCandidate(&topDaughters);
      SimpleParticleCollection_t antitopDaughters;
      antitopDaughters.push_back(SimpleParticle_t(0, jet2));
      for (unsigned int idau=2; idau<4; idau++) antitopDaughters.push_back(SimpleParticle_t(0, pDaughters[idau]));
      mela.appendTopCandidate(&antitopDaughters);

      if (verbosity>=TVar::DEBUG){
        cout << "Mela candidates summary:" << endl;
        for (int ic=0; ic<mela.getNCandidates(); ic++){
          cout << "Candidate " << ic << endl;
          mela.setCurrentCandidateFromIndex(ic);
          TUtil::PrintCandidateSummary(mela.getCurrentCandidate());
        }
        cout << endl;
      }

      mela.setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::ttH);
      mela.computeProdP_ttH(ptth_VAJHU_old_NEW, 2, 0, false);

      mela.setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::bbH);
      mela.computeProdP_ttH(pbbh_VAJHU_old_NEW, 2, 0, false);

      mela.setProcess(TVar::H0minus, TVar::JHUGen, TVar::ttH);
      mela.computeProdP_ttH(ptth0minus_VAJHU_old_NEW, 2, 0, false);

      mela.setProcess(TVar::H0minus, TVar::JHUGen, TVar::bbH);
      mela.computeProdP_ttH(pbbh0minus_VAJHU_old_NEW, 2, 0, false);


      mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ttH);
      mela.selfDHqqcoupl[0][0]=1;
      mela.computeProdP_ttH(ptth_VAJHU_old_NEW_selfD, 2, 0, false);

      mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::bbH);
      mela.selfDHqqcoupl[0][0]=1;
      mela.computeProdP_ttH(pbbh_VAJHU_old_NEW_selfD, 2, 0, false);

      mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ttH);
      mela.selfDHqqcoupl[1][0]=1;
      mela.computeProdP_ttH(ptth0minus_VAJHU_old_NEW_selfD, 2, 0, false);

      mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::bbH);
      mela.selfDHqqcoupl[1][0]=1;
      mela.computeProdP_ttH(pbbh0minus_VAJHU_old_NEW_selfD, 2, 0, false);

      newtree->Fill();
      recorded++;
      mela.resetInputEvent();
    }
  }


  foutput->WriteTObject(newtree);
  delete newtree;
  foutput->Close();
  finput->Close();
}


void testME_Dec_ZZWWComparison_FullSim(){
  int erg_tev=8;
  float mPOLE=125.;
  float wPOLE=4.07e-3;
  TString TREE_NAME = "SelectedTree";

  TVar::VerbosityLevel verbosity = TVar::ERROR;
  Mela mela(erg_tev, mPOLE, verbosity);
  if (verbosity>=TVar::DEBUG) cout << "Mela is initialized" << endl;
  mela.resetMCFM_EWKParameters(1.16639E-05, 1./128., 80.399, 91.1876, 0.23119);

  TString cinput_main = "/scratch0/hep/ianderso/CJLST/140519/PRODFSR_8TeV";
  TFile* finput = new TFile(Form("%s/%s/HZZ4lTree_powheg15jhuGenV3-0PMH125.6.root", cinput_main.Data(), "2mu2e"), "read");
  TFile* foutput = new TFile(Form("HZZ4lTree_powheg15jhuGenV3-0PMH125.6_%s_ZZWWComparison.root", "2l2nu"), "recreate");

  float p0plus_VAJHU_NEW;
  float p0hplus_VAJHU_NEW;
  float p0minus_VAJHU_NEW;
  float p0_g1prime2_VAJHU_NEW;
  float pg1g1prime2_VAJHU_NEW;
  float pg1g1prime2_pi2_VAJHU_NEW;
  float pg1g2_VAJHU_NEW;
  float pg1g2_pi2_VAJHU_NEW;
  float pg1g4_VAJHU_NEW;
  float pg1g4_pi2_VAJHU_NEW;

  float p0hplus_VAMCFM_NEW;
  float p0minus_VAMCFM_NEW;
  float p0_g1prime2_VAMCFM_NEW;
  float pg1g1prime2_VAMCFM_NEW;
  float pg1g1prime2_pi2_VAMCFM_NEW;
  float pg1g2_VAMCFM_NEW;
  float pg1g2_pi2_VAMCFM_NEW;
  float pg1g4_VAMCFM_NEW;
  float pg1g4_pi2_VAMCFM_NEW;

  float p0hplus_VAMCFM_ratio_NEW;
  float p0minus_VAMCFM_ratio_NEW;
  float p0_g1prime2_VAMCFM_ratio_NEW;
  float pg1g1prime2_VAMCFM_ratio_NEW;
  float pg1g1prime2_pi2_VAMCFM_ratio_NEW;
  float pg1g2_VAMCFM_ratio_NEW;
  float pg1g2_pi2_VAMCFM_ratio_NEW;
  float pg1g4_VAMCFM_ratio_NEW;
  float pg1g4_pi2_VAMCFM_ratio_NEW;

  float p0hplus_VAJHU_ratio_NEW;
  float p0minus_VAJHU_ratio_NEW;
  float p0_g1prime2_VAJHU_ratio_NEW;
  float pg1g1prime2_VAJHU_ratio_NEW;
  float pg1g1prime2_pi2_VAJHU_ratio_NEW;
  float pg1g2_VAJHU_ratio_NEW;
  float pg1g2_pi2_VAJHU_ratio_NEW;
  float pg1g4_VAJHU_ratio_NEW;
  float pg1g4_pi2_VAJHU_ratio_NEW;

  float p0plus_VAMCFM_NEW, ggzz_p0plus_VAMCFM_NEW, p0plus_VAMCFM_NEW_BSMOn, ggzz_p0plus_VAMCFM_NEW_BSMOn;
  float bkg_VAMCFM_NEW, ggzz_VAMCFM_NEW;

  float p0plus_WW_VAJHU_NEW;
  float p0hplus_WW_VAJHU_NEW;
  float p0minus_WW_VAJHU_NEW;
  float p0_g1prime2_WW_VAJHU_NEW;
  float pg1g1prime2_WW_VAJHU_NEW;
  float pg1g1prime2_pi2_WW_VAJHU_NEW;
  float pg1g2_WW_VAJHU_NEW;
  float pg1g2_pi2_WW_VAJHU_NEW;
  float pg1g4_WW_VAJHU_NEW;
  float pg1g4_pi2_WW_VAJHU_NEW;

  float p0hplus_WW_VAMCFM_NEW;
  float p0minus_WW_VAMCFM_NEW;
  float p0_g1prime2_WW_VAMCFM_NEW;
  float pg1g1prime2_WW_VAMCFM_NEW;
  float pg1g1prime2_pi2_WW_VAMCFM_NEW;
  float pg1g2_WW_VAMCFM_NEW;
  float pg1g2_pi2_WW_VAMCFM_NEW;
  float pg1g4_WW_VAMCFM_NEW;
  float pg1g4_pi2_WW_VAMCFM_NEW;

  float p0hplus_WW_VAMCFM_ratio_NEW;
  float p0minus_WW_VAMCFM_ratio_NEW;
  float p0_g1prime2_WW_VAMCFM_ratio_NEW;
  float pg1g1prime2_WW_VAMCFM_ratio_NEW;
  float pg1g1prime2_pi2_WW_VAMCFM_ratio_NEW;
  float pg1g2_WW_VAMCFM_ratio_NEW;
  float pg1g2_pi2_WW_VAMCFM_ratio_NEW;
  float pg1g4_WW_VAMCFM_ratio_NEW;
  float pg1g4_pi2_WW_VAMCFM_ratio_NEW;

  float p0hplus_WW_VAJHU_ratio_NEW;
  float p0minus_WW_VAJHU_ratio_NEW;
  float p0_g1prime2_WW_VAJHU_ratio_NEW;
  float pg1g1prime2_WW_VAJHU_ratio_NEW;
  float pg1g1prime2_pi2_WW_VAJHU_ratio_NEW;
  float pg1g2_WW_VAJHU_ratio_NEW;
  float pg1g2_pi2_WW_VAJHU_ratio_NEW;
  float pg1g4_WW_VAJHU_ratio_NEW;
  float pg1g4_pi2_WW_VAJHU_ratio_NEW;

  float p0plus_WW_VAMCFM_NEW, ggzz_p0plus_WW_VAMCFM_NEW, p0plus_WW_VAMCFM_NEW_BSMOn, ggzz_p0plus_WW_VAMCFM_NEW_BSMOn;
  float bkg_WW_VAMCFM_NEW, ggzz_WW_VAMCFM_NEW;

  float mzz = 126.;
  float m1 = 91.471450;
  float m2 = 12.139782;
  float h1 = 0.2682896;
  float h2 = 0.1679779;
  float phi = 1.5969792;
  float hs = -0.727181;
  float phi1 = 1.8828257;

  TTree* tree = (TTree*)finput->Get(TREE_NAME);
  tree->SetBranchAddress("ZZMass", &mzz);
  tree->SetBranchAddress("Z1Mass", &m1);
  tree->SetBranchAddress("Z2Mass", &m2);
  tree->SetBranchAddress("helcosthetaZ1", &h1);
  tree->SetBranchAddress("helcosthetaZ2", &h2);
  tree->SetBranchAddress("helphi", &phi);
  tree->SetBranchAddress("costhetastar", &hs);
  tree->SetBranchAddress("phistarZ1", &phi1);

  TTree* newtree = new TTree("TestTree", "");
  newtree->Branch("ZZMass", &mzz);
  newtree->Branch("Z1Mass", &m1);
  newtree->Branch("Z2Mass", &m2);
  newtree->Branch("helcosthetaZ1", &h1);
  newtree->Branch("helcosthetaZ2", &h2);
  newtree->Branch("helphi", &phi);
  newtree->Branch("costhetastar", &hs);
  newtree->Branch("phistarZ1", &phi1);

  newtree->Branch("p0plus_VAJHU_NEW", &p0plus_VAJHU_NEW);
  newtree->Branch("p0hplus_VAJHU_NEW", &p0hplus_VAJHU_NEW);
  newtree->Branch("p0minus_VAJHU_NEW", &p0minus_VAJHU_NEW);
  newtree->Branch("p0_g1prime2_VAJHU_NEW", &p0_g1prime2_VAJHU_NEW);
  newtree->Branch("pg1g2_VAJHU_NEW", &pg1g2_VAJHU_NEW);
  newtree->Branch("pg1g2_pi2_VAJHU_NEW", &pg1g2_pi2_VAJHU_NEW);
  newtree->Branch("pg1g4_VAJHU_NEW", &pg1g4_VAJHU_NEW);
  newtree->Branch("pg1g4_pi2_VAJHU_NEW", &pg1g4_pi2_VAJHU_NEW);
  newtree->Branch("pg1g1prime2_VAJHU_NEW", &pg1g1prime2_VAJHU_NEW);
  newtree->Branch("pg1g1prime2_pi2_VAJHU_NEW", &pg1g1prime2_pi2_VAJHU_NEW);

  newtree->Branch("bkg_VAMCFM_NEW", &bkg_VAMCFM_NEW);
  newtree->Branch("ggzz_VAMCFM_NEW", &ggzz_VAMCFM_NEW);
  newtree->Branch("p0plus_VAMCFM_NEW", &p0plus_VAMCFM_NEW);
  newtree->Branch("ggzz_p0plus_VAMCFM_NEW", &ggzz_p0plus_VAMCFM_NEW);
  newtree->Branch("p0plus_VAMCFM_NEW_BSMOn", &p0plus_VAMCFM_NEW_BSMOn);
  newtree->Branch("ggzz_p0plus_VAMCFM_NEW_BSMOn", &ggzz_p0plus_VAMCFM_NEW_BSMOn);

  newtree->Branch("p0hplus_VAMCFM_NEW", &p0hplus_VAMCFM_NEW);
  newtree->Branch("p0minus_VAMCFM_NEW", &p0minus_VAMCFM_NEW);
  newtree->Branch("p0_g1prime2_VAMCFM_NEW", &p0_g1prime2_VAMCFM_NEW);
  newtree->Branch("pg1g2_VAMCFM_NEW", &pg1g2_VAMCFM_NEW);
  newtree->Branch("pg1g2_pi2_VAMCFM_NEW", &pg1g2_pi2_VAMCFM_NEW);
  newtree->Branch("pg1g4_VAMCFM_NEW", &pg1g4_VAMCFM_NEW);
  newtree->Branch("pg1g4_pi2_VAMCFM_NEW", &pg1g4_pi2_VAMCFM_NEW);
  newtree->Branch("pg1g1prime2_VAMCFM_NEW", &pg1g1prime2_VAMCFM_NEW);
  newtree->Branch("pg1g1prime2_pi2_VAMCFM_NEW", &pg1g1prime2_pi2_VAMCFM_NEW);

  newtree->Branch("p0hplus_VAMCFM_ratio_NEW", &p0hplus_VAMCFM_ratio_NEW);
  newtree->Branch("p0minus_VAMCFM_ratio_NEW", &p0minus_VAMCFM_ratio_NEW);
  newtree->Branch("p0_g1prime2_VAMCFM_ratio_NEW", &p0_g1prime2_VAMCFM_ratio_NEW);
  newtree->Branch("pg1g2_VAMCFM_ratio_NEW", &pg1g2_VAMCFM_ratio_NEW);
  newtree->Branch("pg1g2_pi2_VAMCFM_ratio_NEW", &pg1g2_pi2_VAMCFM_ratio_NEW);
  newtree->Branch("pg1g4_VAMCFM_ratio_NEW", &pg1g4_VAMCFM_ratio_NEW);
  newtree->Branch("pg1g4_pi2_VAMCFM_ratio_NEW", &pg1g4_pi2_VAMCFM_ratio_NEW);
  newtree->Branch("pg1g1prime2_VAMCFM_ratio_NEW", &pg1g1prime2_VAMCFM_ratio_NEW);
  newtree->Branch("pg1g1prime2_pi2_VAMCFM_ratio_NEW", &pg1g1prime2_pi2_VAMCFM_ratio_NEW);

  newtree->Branch("p0hplus_VAJHU_ratio_NEW", &p0hplus_VAJHU_ratio_NEW);
  newtree->Branch("p0minus_VAJHU_ratio_NEW", &p0minus_VAJHU_ratio_NEW);
  newtree->Branch("p0_g1prime2_VAJHU_ratio_NEW", &p0_g1prime2_VAJHU_ratio_NEW);
  newtree->Branch("pg1g2_VAJHU_ratio_NEW", &pg1g2_VAJHU_ratio_NEW);
  newtree->Branch("pg1g2_pi2_VAJHU_ratio_NEW", &pg1g2_pi2_VAJHU_ratio_NEW);
  newtree->Branch("pg1g4_VAJHU_ratio_NEW", &pg1g4_VAJHU_ratio_NEW);
  newtree->Branch("pg1g4_pi2_VAJHU_ratio_NEW", &pg1g4_pi2_VAJHU_ratio_NEW);
  newtree->Branch("pg1g1prime2_VAJHU_ratio_NEW", &pg1g1prime2_VAJHU_ratio_NEW);
  newtree->Branch("pg1g1prime2_pi2_VAJHU_ratio_NEW", &pg1g1prime2_pi2_VAJHU_ratio_NEW);

  newtree->Branch("p0plus_WW_VAJHU_NEW", &p0plus_WW_VAJHU_NEW);
  newtree->Branch("p0hplus_WW_VAJHU_NEW", &p0hplus_WW_VAJHU_NEW);
  newtree->Branch("p0minus_WW_VAJHU_NEW", &p0minus_WW_VAJHU_NEW);
  newtree->Branch("p0_g1prime2_WW_VAJHU_NEW", &p0_g1prime2_WW_VAJHU_NEW);
  newtree->Branch("pg1g2_WW_VAJHU_NEW", &pg1g2_WW_VAJHU_NEW);
  newtree->Branch("pg1g2_pi2_WW_VAJHU_NEW", &pg1g2_pi2_WW_VAJHU_NEW);
  newtree->Branch("pg1g4_WW_VAJHU_NEW", &pg1g4_WW_VAJHU_NEW);
  newtree->Branch("pg1g4_pi2_WW_VAJHU_NEW", &pg1g4_pi2_WW_VAJHU_NEW);
  newtree->Branch("pg1g1prime2_WW_VAJHU_NEW", &pg1g1prime2_WW_VAJHU_NEW);
  newtree->Branch("pg1g1prime2_pi2_WW_VAJHU_NEW", &pg1g1prime2_pi2_WW_VAJHU_NEW);

  newtree->Branch("bkg_WW_VAMCFM_NEW", &bkg_WW_VAMCFM_NEW);
  newtree->Branch("ggzz_WW_VAMCFM_NEW", &ggzz_WW_VAMCFM_NEW);
  newtree->Branch("p0plus_WW_VAMCFM_NEW", &p0plus_WW_VAMCFM_NEW);
  newtree->Branch("ggzz_p0plus_WW_VAMCFM_NEW", &ggzz_p0plus_WW_VAMCFM_NEW);
  newtree->Branch("p0plus_WW_VAMCFM_NEW_BSMOn", &p0plus_WW_VAMCFM_NEW_BSMOn);
  newtree->Branch("ggzz_p0plus_WW_VAMCFM_NEW_BSMOn", &ggzz_p0plus_WW_VAMCFM_NEW_BSMOn);

  newtree->Branch("p0hplus_WW_VAMCFM_NEW", &p0hplus_WW_VAMCFM_NEW);
  newtree->Branch("p0minus_WW_VAMCFM_NEW", &p0minus_WW_VAMCFM_NEW);
  newtree->Branch("p0_g1prime2_WW_VAMCFM_NEW", &p0_g1prime2_WW_VAMCFM_NEW);
  newtree->Branch("pg1g2_WW_VAMCFM_NEW", &pg1g2_WW_VAMCFM_NEW);
  newtree->Branch("pg1g2_pi2_WW_VAMCFM_NEW", &pg1g2_pi2_WW_VAMCFM_NEW);
  newtree->Branch("pg1g4_WW_VAMCFM_NEW", &pg1g4_WW_VAMCFM_NEW);
  newtree->Branch("pg1g4_pi2_WW_VAMCFM_NEW", &pg1g4_pi2_WW_VAMCFM_NEW);
  newtree->Branch("pg1g1prime2_WW_VAMCFM_NEW", &pg1g1prime2_WW_VAMCFM_NEW);
  newtree->Branch("pg1g1prime2_pi2_WW_VAMCFM_NEW", &pg1g1prime2_pi2_WW_VAMCFM_NEW);

  newtree->Branch("p0hplus_WW_VAMCFM_ratio_NEW", &p0hplus_WW_VAMCFM_ratio_NEW);
  newtree->Branch("p0minus_WW_VAMCFM_ratio_NEW", &p0minus_WW_VAMCFM_ratio_NEW);
  newtree->Branch("p0_g1prime2_WW_VAMCFM_ratio_NEW", &p0_g1prime2_WW_VAMCFM_ratio_NEW);
  newtree->Branch("pg1g2_WW_VAMCFM_ratio_NEW", &pg1g2_WW_VAMCFM_ratio_NEW);
  newtree->Branch("pg1g2_pi2_WW_VAMCFM_ratio_NEW", &pg1g2_pi2_WW_VAMCFM_ratio_NEW);
  newtree->Branch("pg1g4_WW_VAMCFM_ratio_NEW", &pg1g4_WW_VAMCFM_ratio_NEW);
  newtree->Branch("pg1g4_pi2_WW_VAMCFM_ratio_NEW", &pg1g4_pi2_WW_VAMCFM_ratio_NEW);
  newtree->Branch("pg1g1prime2_WW_VAMCFM_ratio_NEW", &pg1g1prime2_WW_VAMCFM_ratio_NEW);
  newtree->Branch("pg1g1prime2_pi2_WW_VAMCFM_ratio_NEW", &pg1g1prime2_pi2_WW_VAMCFM_ratio_NEW);

  newtree->Branch("p0hplus_WW_VAJHU_ratio_NEW", &p0hplus_WW_VAJHU_ratio_NEW);
  newtree->Branch("p0minus_WW_VAJHU_ratio_NEW", &p0minus_WW_VAJHU_ratio_NEW);
  newtree->Branch("p0_g1prime2_WW_VAJHU_ratio_NEW", &p0_g1prime2_WW_VAJHU_ratio_NEW);
  newtree->Branch("pg1g2_WW_VAJHU_ratio_NEW", &pg1g2_WW_VAJHU_ratio_NEW);
  newtree->Branch("pg1g2_pi2_WW_VAJHU_ratio_NEW", &pg1g2_pi2_WW_VAJHU_ratio_NEW);
  newtree->Branch("pg1g4_WW_VAJHU_ratio_NEW", &pg1g4_WW_VAJHU_ratio_NEW);
  newtree->Branch("pg1g4_pi2_WW_VAJHU_ratio_NEW", &pg1g4_pi2_WW_VAJHU_ratio_NEW);
  newtree->Branch("pg1g1prime2_WW_VAJHU_ratio_NEW", &pg1g1prime2_WW_VAJHU_ratio_NEW);
  newtree->Branch("pg1g1prime2_pi2_WW_VAJHU_ratio_NEW", &pg1g1prime2_pi2_WW_VAJHU_ratio_NEW);

  float GenLep1Id=0, GenLep2Id=0, GenLep3Id=0, GenLep4Id=0;
  GenLep1Id=14;
  GenLep2Id=-14;
  GenLep3Id=13;
  GenLep4Id=-13;

  for (int ev = 0; ev < min(1000, (int)tree->GetEntries()); ev++){
    tree->GetEntry(ev);

    int idOrdered[4] ={ static_cast<int>(GenLep1Id), static_cast<int>(GenLep2Id), static_cast<int>(GenLep3Id), static_cast<int>(GenLep4Id) };
    TLorentzVector pOrdered[4];
    std::vector<TLorentzVector> daus = mela.calculate4Momentum(mzz, m1, m2, acos(hs), acos(h1), acos(h2), phi1, phi);
    for (int ip=0; ip<min(4, (int)daus.size()); ip++) pOrdered[ip]=daus.at(ip);
    SimpleParticleCollection_t daughters;
    for (unsigned int idau=0; idau<4; idau++) daughters.push_back(SimpleParticle_t(idOrdered[idau], pOrdered[idau]));
    mela.setCandidateDecayMode(TVar::CandidateDecay_ZZ);
    mela.setInputEvent(&daughters, (SimpleParticleCollection_t*)0, (SimpleParticleCollection_t*)0, false);
    mela.setCandidateDecayMode(TVar::CandidateDecay_WW);
    mela.setInputEvent(&daughters, (SimpleParticleCollection_t*)0, (SimpleParticleCollection_t*)0, false);
    if (verbosity>=TVar::DEBUG){
      cout << "Mela candidate summary:" << endl;
      for (int ic=0; ic<mela.getNCandidates(); ic++){
        cout << "Candidate " << ic << endl;
        mela.setCurrentCandidateFromIndex(ic);
        TUtil::PrintCandidateSummary(mela.getCurrentCandidate());
      }
      cout << endl;
    }

    /***** ZZ *****/
    mela.setCurrentCandidateFromIndex(0);

    if (verbosity>=TVar::DEBUG) cout << "Computing bkg_VAMCFM_NEW" << endl;
    mela.setProcess(TVar::bkgZZ, TVar::MCFM, TVar::ZZQQB);
    mela.computeP(bkg_VAMCFM_NEW, false);

    if (verbosity>=TVar::DEBUG) cout << "Computing ggzz_VAMCFM_NEW" << endl;
    mela.setProcess(TVar::bkgZZ, TVar::MCFM, TVar::ZZGG);
    mela.computeP(ggzz_VAMCFM_NEW, false);

    if (verbosity>=TVar::DEBUG) cout << "Computing p0plus_VAMCFM_NEW" << endl;
    mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.computeP(p0plus_VAMCFM_NEW, false);

    if (verbosity>=TVar::DEBUG) cout << "Computing ggzz_p0plus_VAMCFM_NEW" << endl;
    mela.setProcess(TVar::bkgZZ_SMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.computeP(ggzz_p0plus_VAMCFM_NEW, false);

    if (verbosity>=TVar::DEBUG) cout << "Computing p0plus_VAMCFM_NEW_BSMOn" << endl;
    mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.selfDHzzcoupl[0][0][0]=1;
    mela.computeP(p0plus_VAMCFM_NEW_BSMOn, false);

    if (verbosity>=TVar::DEBUG) cout << "Computing ggzz_p0plus_VAMCFM_NEW_BSMOn" << endl;
    mela.setProcess(TVar::bkgZZ_SMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.selfDHzzcoupl[0][0][0]=1;
    mela.computeP(ggzz_p0plus_VAMCFM_NEW_BSMOn, false);

    mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);

    mela.selfDHggcoupl[0][0]=1;
    mela.selfDHzzcoupl[0][0][0]=1;
    mela.computeP(p0plus_VAJHU_NEW, false);

    mela.selfDHggcoupl[0][0]=1;
    mela.selfDHzzcoupl[0][1][0]=1.638;
    mela.computeP(p0hplus_VAJHU_NEW, false);

    mela.selfDHggcoupl[0][0]=1;
    mela.selfDHzzcoupl[0][3][0]=2.521;
    mela.computeP(p0minus_VAJHU_NEW, false);

    mela.selfDHggcoupl[0][0]=1;
    mela.selfDHzzcoupl[0][11][0]=-12046.01;
    mela.computeP(p0_g1prime2_VAJHU_NEW, false);

    mela.selfDHggcoupl[0][0]=1;
    mela.selfDHzzcoupl[0][0][0]=1;
    mela.selfDHzzcoupl[0][1][0]=1.638;
    mela.computeP(pg1g2_VAJHU_NEW, false);
    pg1g2_VAJHU_NEW -= (p0plus_VAJHU_NEW + p0hplus_VAJHU_NEW);

    mela.selfDHggcoupl[0][0]=1;
    mela.selfDHzzcoupl[0][0][0]=1;
    mela.selfDHzzcoupl[0][1][1]=1.638;
    mela.computeP(pg1g2_pi2_VAJHU_NEW, false);
    pg1g2_pi2_VAJHU_NEW -= (p0plus_VAJHU_NEW + p0hplus_VAJHU_NEW);

    mela.selfDHggcoupl[0][0]=1;
    mela.selfDHzzcoupl[0][0][0]=1;
    mela.selfDHzzcoupl[0][3][0]=2.521;
    mela.computeP(pg1g4_VAJHU_NEW, false);
    pg1g4_VAJHU_NEW -= (p0plus_VAJHU_NEW + p0minus_VAJHU_NEW);

    mela.selfDHggcoupl[0][0]=1;
    mela.selfDHzzcoupl[0][0][0]=1;
    mela.selfDHzzcoupl[0][3][1]=2.521;
    mela.computeP(pg1g4_pi2_VAJHU_NEW, false);
    pg1g4_pi2_VAJHU_NEW -= (p0plus_VAJHU_NEW + p0minus_VAJHU_NEW);

    mela.selfDHggcoupl[0][0]=1;
    mela.selfDHzzcoupl[0][0][0]=1;
    mela.selfDHzzcoupl[0][11][0]=12046.01;
    mela.computeP(pg1g1prime2_VAJHU_NEW, false);
    pg1g1prime2_VAJHU_NEW -= (p0plus_VAJHU_NEW + p0_g1prime2_VAJHU_NEW);

    mela.selfDHggcoupl[0][0]=1;
    mela.selfDHzzcoupl[0][0][0]=1;
    mela.selfDHzzcoupl[0][11][1]=12046.01;
    mela.computeP(pg1g1prime2_pi2_VAJHU_NEW, false);
    pg1g1prime2_pi2_VAJHU_NEW -= (p0plus_VAJHU_NEW + p0_g1prime2_VAJHU_NEW);


    mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);

    mela.selfDHzzcoupl[0][0][0]=1;
    mela.computeP(p0plus_VAMCFM_NEW, false);

    mela.selfDHzzcoupl[0][1][0]=1.638;
    mela.computeP(p0hplus_VAMCFM_NEW, false);

    mela.selfDHzzcoupl[0][3][0]=2.521;
    mela.computeP(p0minus_VAMCFM_NEW, false);

    mela.selfDHzzcoupl[0][11][0]=-12046.01;
    mela.computeP(p0_g1prime2_VAMCFM_NEW, false);

    mela.selfDHzzcoupl[0][0][0]=1;
    mela.selfDHzzcoupl[0][1][0]=1.638;
    mela.computeP(pg1g2_VAMCFM_NEW, false);
    pg1g2_VAMCFM_NEW -= (p0plus_VAMCFM_NEW + p0hplus_VAMCFM_NEW);

    mela.selfDHzzcoupl[0][0][0]=1;
    mela.selfDHzzcoupl[0][1][1]=1.638;
    mela.computeP(pg1g2_pi2_VAMCFM_NEW, false);
    pg1g2_pi2_VAMCFM_NEW -= (p0plus_VAMCFM_NEW + p0hplus_VAMCFM_NEW);

    mela.selfDHzzcoupl[0][0][0]=1;
    mela.selfDHzzcoupl[0][3][0]=2.521;
    mela.computeP(pg1g4_VAMCFM_NEW, false);
    pg1g4_VAMCFM_NEW -= (p0plus_VAMCFM_NEW + p0minus_VAMCFM_NEW);

    mela.selfDHzzcoupl[0][0][0]=1;
    mela.selfDHzzcoupl[0][3][1]=2.521;
    mela.computeP(pg1g4_pi2_VAMCFM_NEW, false);
    pg1g4_pi2_VAMCFM_NEW -= (p0plus_VAMCFM_NEW + p0minus_VAMCFM_NEW);

    mela.selfDHzzcoupl[0][0][0]=1;
    mela.selfDHzzcoupl[0][11][0]=12046.01;
    mela.computeP(pg1g1prime2_VAMCFM_NEW, false);
    pg1g1prime2_VAMCFM_NEW -= (p0plus_VAMCFM_NEW + p0_g1prime2_VAMCFM_NEW);

    mela.selfDHzzcoupl[0][0][0]=1;
    mela.selfDHzzcoupl[0][11][1]=12046.01;
    mela.computeP(pg1g1prime2_pi2_VAMCFM_NEW, false);
    pg1g1prime2_pi2_VAMCFM_NEW -= (p0plus_VAMCFM_NEW + p0_g1prime2_VAMCFM_NEW);


    p0hplus_VAJHU_ratio_NEW = p0hplus_VAJHU_NEW/p0plus_VAJHU_NEW;
    p0minus_VAJHU_ratio_NEW = p0minus_VAJHU_NEW/p0plus_VAJHU_NEW;
    p0_g1prime2_VAJHU_ratio_NEW = p0_g1prime2_VAJHU_NEW/p0plus_VAJHU_NEW;
    pg1g2_VAJHU_ratio_NEW = pg1g2_VAJHU_NEW/p0plus_VAJHU_NEW;
    pg1g4_VAJHU_ratio_NEW = pg1g4_VAJHU_NEW/p0plus_VAJHU_NEW;
    pg1g1prime2_VAJHU_ratio_NEW = pg1g1prime2_VAJHU_NEW/p0plus_VAJHU_NEW;
    pg1g2_pi2_VAJHU_ratio_NEW = pg1g2_pi2_VAJHU_NEW/p0plus_VAJHU_NEW;
    pg1g4_pi2_VAJHU_ratio_NEW = pg1g4_pi2_VAJHU_NEW/p0plus_VAJHU_NEW;
    pg1g1prime2_pi2_VAJHU_ratio_NEW = pg1g1prime2_pi2_VAJHU_NEW/p0plus_VAJHU_NEW;

    p0hplus_VAMCFM_ratio_NEW = p0hplus_VAMCFM_NEW/p0plus_VAMCFM_NEW;
    p0minus_VAMCFM_ratio_NEW = p0minus_VAMCFM_NEW/p0plus_VAMCFM_NEW;
    p0_g1prime2_VAMCFM_ratio_NEW = p0_g1prime2_VAMCFM_NEW/p0plus_VAMCFM_NEW;
    pg1g2_VAMCFM_ratio_NEW = pg1g2_VAMCFM_NEW/p0plus_VAMCFM_NEW;
    pg1g4_VAMCFM_ratio_NEW = pg1g4_VAMCFM_NEW/p0plus_VAMCFM_NEW;
    pg1g1prime2_VAMCFM_ratio_NEW = pg1g1prime2_VAMCFM_NEW/p0plus_VAMCFM_NEW;
    pg1g2_pi2_VAMCFM_ratio_NEW = pg1g2_pi2_VAMCFM_NEW/p0plus_VAMCFM_NEW;
    pg1g4_pi2_VAMCFM_ratio_NEW = pg1g4_pi2_VAMCFM_NEW/p0plus_VAMCFM_NEW;
    pg1g1prime2_pi2_VAMCFM_ratio_NEW = pg1g1prime2_pi2_VAMCFM_NEW/p0plus_VAMCFM_NEW;



    /***** WW *****/
    mela.setCurrentCandidateFromIndex(1);

    if (verbosity>=TVar::DEBUG) cout << "Computing bkg_WW_VAMCFM_NEW" << endl;
    mela.setProcess(TVar::bkgWW, TVar::MCFM, TVar::ZZQQB);
    mela.computeP(bkg_WW_VAMCFM_NEW, false);

    if (verbosity>=TVar::DEBUG) cout << "Computing ggzz_WW_VAMCFM_NEW" << endl;
    mela.setProcess(TVar::bkgWW, TVar::MCFM, TVar::ZZGG);
    mela.computeP(ggzz_WW_VAMCFM_NEW, false);

    if (verbosity>=TVar::DEBUG) cout << "Computing p0plus_WW_VAMCFM_NEW" << endl;
    mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.computeP(p0plus_WW_VAMCFM_NEW, false);

    if (verbosity>=TVar::DEBUG) cout << "Computing ggzz_p0plus_WW_VAMCFM_NEW" << endl;
    mela.setProcess(TVar::bkgWW_SMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.computeP(ggzz_p0plus_WW_VAMCFM_NEW, false);

    if (verbosity>=TVar::DEBUG) cout << "Computing p0plus_WW_VAMCFM_NEW_BSMOn" << endl;
    mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.selfDHzzcoupl[0][0][0]=1;
    mela.computeP(p0plus_WW_VAMCFM_NEW_BSMOn, false);

    if (verbosity>=TVar::DEBUG) cout << "Computing ggzz_p0plus_WW_VAMCFM_NEW_BSMOn" << endl;
    mela.setProcess(TVar::bkgWW_SMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.selfDHzzcoupl[0][0][0]=1;
    mela.computeP(ggzz_p0plus_WW_VAMCFM_NEW_BSMOn, false);

    mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);

    mela.selfDHggcoupl[0][0]=1;
    mela.selfDHzzcoupl[0][0][0]=1;
    mela.computeP(p0plus_WW_VAJHU_NEW, false);

    mela.selfDHggcoupl[0][0]=1;
    mela.selfDHzzcoupl[0][1][0]=1.638;
    mela.computeP(p0hplus_WW_VAJHU_NEW, false);

    mela.selfDHggcoupl[0][0]=1;
    mela.selfDHzzcoupl[0][3][0]=2.521;
    mela.computeP(p0minus_WW_VAJHU_NEW, false);

    mela.selfDHggcoupl[0][0]=1;
    mela.selfDHzzcoupl[0][11][0]=-12046.01;
    mela.computeP(p0_g1prime2_WW_VAJHU_NEW, false);

    mela.selfDHggcoupl[0][0]=1;
    mela.selfDHzzcoupl[0][0][0]=1;
    mela.selfDHzzcoupl[0][1][0]=1.638;
    mela.computeP(pg1g2_WW_VAJHU_NEW, false);
    pg1g2_WW_VAJHU_NEW -= (p0plus_WW_VAJHU_NEW + p0hplus_WW_VAJHU_NEW);

    mela.selfDHggcoupl[0][0]=1;
    mela.selfDHzzcoupl[0][0][0]=1;
    mela.selfDHzzcoupl[0][1][1]=1.638;
    mela.computeP(pg1g2_pi2_WW_VAJHU_NEW, false);
    pg1g2_pi2_WW_VAJHU_NEW -= (p0plus_WW_VAJHU_NEW + p0hplus_WW_VAJHU_NEW);

    mela.selfDHggcoupl[0][0]=1;
    mela.selfDHzzcoupl[0][0][0]=1;
    mela.selfDHzzcoupl[0][3][0]=2.521;
    mela.computeP(pg1g4_WW_VAJHU_NEW, false);
    pg1g4_WW_VAJHU_NEW -= (p0plus_WW_VAJHU_NEW + p0minus_WW_VAJHU_NEW);

    mela.selfDHggcoupl[0][0]=1;
    mela.selfDHzzcoupl[0][0][0]=1;
    mela.selfDHzzcoupl[0][3][1]=2.521;
    mela.computeP(pg1g4_pi2_WW_VAJHU_NEW, false);
    pg1g4_pi2_WW_VAJHU_NEW -= (p0plus_WW_VAJHU_NEW + p0minus_WW_VAJHU_NEW);

    mela.selfDHggcoupl[0][0]=1;
    mela.selfDHzzcoupl[0][0][0]=1;
    mela.selfDHzzcoupl[0][11][0]=12046.01;
    mela.computeP(pg1g1prime2_WW_VAJHU_NEW, false);
    pg1g1prime2_WW_VAJHU_NEW -= (p0plus_WW_VAJHU_NEW + p0_g1prime2_WW_VAJHU_NEW);

    mela.selfDHggcoupl[0][0]=1;
    mela.selfDHzzcoupl[0][0][0]=1;
    mela.selfDHzzcoupl[0][11][1]=12046.01;
    mela.computeP(pg1g1prime2_pi2_WW_VAJHU_NEW, false);
    pg1g1prime2_pi2_WW_VAJHU_NEW -= (p0plus_WW_VAJHU_NEW + p0_g1prime2_WW_VAJHU_NEW);


    mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);

    mela.selfDHzzcoupl[0][0][0]=1;
    mela.computeP(p0plus_WW_VAMCFM_NEW, false);

    mela.selfDHzzcoupl[0][1][0]=1.638;
    mela.computeP(p0hplus_WW_VAMCFM_NEW, false);

    mela.selfDHzzcoupl[0][3][0]=2.521;
    mela.computeP(p0minus_WW_VAMCFM_NEW, false);

    mela.selfDHzzcoupl[0][11][0]=-12046.01;
    mela.computeP(p0_g1prime2_WW_VAMCFM_NEW, false);

    mela.selfDHzzcoupl[0][0][0]=1;
    mela.selfDHzzcoupl[0][1][0]=1.638;
    mela.computeP(pg1g2_WW_VAMCFM_NEW, false);
    pg1g2_WW_VAMCFM_NEW -= (p0plus_WW_VAMCFM_NEW + p0hplus_WW_VAMCFM_NEW);

    mela.selfDHzzcoupl[0][0][0]=1;
    mela.selfDHzzcoupl[0][1][1]=1.638;
    mela.computeP(pg1g2_pi2_WW_VAMCFM_NEW, false);
    pg1g2_pi2_WW_VAMCFM_NEW -= (p0plus_WW_VAMCFM_NEW + p0hplus_WW_VAMCFM_NEW);

    mela.selfDHzzcoupl[0][0][0]=1;
    mela.selfDHzzcoupl[0][3][0]=2.521;
    mela.computeP(pg1g4_WW_VAMCFM_NEW, false);
    pg1g4_WW_VAMCFM_NEW -= (p0plus_WW_VAMCFM_NEW + p0minus_WW_VAMCFM_NEW);

    mela.selfDHzzcoupl[0][0][0]=1;
    mela.selfDHzzcoupl[0][3][1]=2.521;
    mela.computeP(pg1g4_pi2_WW_VAMCFM_NEW, false);
    pg1g4_pi2_WW_VAMCFM_NEW -= (p0plus_WW_VAMCFM_NEW + p0minus_WW_VAMCFM_NEW);

    mela.selfDHzzcoupl[0][0][0]=1;
    mela.selfDHzzcoupl[0][11][0]=12046.01;
    mela.computeP(pg1g1prime2_WW_VAMCFM_NEW, false);
    pg1g1prime2_WW_VAMCFM_NEW -= (p0plus_WW_VAMCFM_NEW + p0_g1prime2_WW_VAMCFM_NEW);

    mela.selfDHzzcoupl[0][0][0]=1;
    mela.selfDHzzcoupl[0][11][1]=12046.01;
    mela.computeP(pg1g1prime2_pi2_WW_VAMCFM_NEW, false);
    pg1g1prime2_pi2_WW_VAMCFM_NEW -= (p0plus_WW_VAMCFM_NEW + p0_g1prime2_WW_VAMCFM_NEW);


    p0hplus_WW_VAJHU_ratio_NEW = p0hplus_WW_VAJHU_NEW/p0plus_WW_VAJHU_NEW;
    p0minus_WW_VAJHU_ratio_NEW = p0minus_WW_VAJHU_NEW/p0plus_WW_VAJHU_NEW;
    p0_g1prime2_WW_VAJHU_ratio_NEW = p0_g1prime2_WW_VAJHU_NEW/p0plus_WW_VAJHU_NEW;
    pg1g2_WW_VAJHU_ratio_NEW = pg1g2_WW_VAJHU_NEW/p0plus_WW_VAJHU_NEW;
    pg1g4_WW_VAJHU_ratio_NEW = pg1g4_WW_VAJHU_NEW/p0plus_WW_VAJHU_NEW;
    pg1g1prime2_WW_VAJHU_ratio_NEW = pg1g1prime2_WW_VAJHU_NEW/p0plus_WW_VAJHU_NEW;
    pg1g2_pi2_WW_VAJHU_ratio_NEW = pg1g2_pi2_WW_VAJHU_NEW/p0plus_WW_VAJHU_NEW;
    pg1g4_pi2_WW_VAJHU_ratio_NEW = pg1g4_pi2_WW_VAJHU_NEW/p0plus_WW_VAJHU_NEW;
    pg1g1prime2_pi2_WW_VAJHU_ratio_NEW = pg1g1prime2_pi2_WW_VAJHU_NEW/p0plus_WW_VAJHU_NEW;

    p0hplus_WW_VAMCFM_ratio_NEW = p0hplus_WW_VAMCFM_NEW/p0plus_WW_VAMCFM_NEW;
    p0minus_WW_VAMCFM_ratio_NEW = p0minus_WW_VAMCFM_NEW/p0plus_WW_VAMCFM_NEW;
    p0_g1prime2_WW_VAMCFM_ratio_NEW = p0_g1prime2_WW_VAMCFM_NEW/p0plus_WW_VAMCFM_NEW;
    pg1g2_WW_VAMCFM_ratio_NEW = pg1g2_WW_VAMCFM_NEW/p0plus_WW_VAMCFM_NEW;
    pg1g4_WW_VAMCFM_ratio_NEW = pg1g4_WW_VAMCFM_NEW/p0plus_WW_VAMCFM_NEW;
    pg1g1prime2_WW_VAMCFM_ratio_NEW = pg1g1prime2_WW_VAMCFM_NEW/p0plus_WW_VAMCFM_NEW;
    pg1g2_pi2_WW_VAMCFM_ratio_NEW = pg1g2_pi2_WW_VAMCFM_NEW/p0plus_WW_VAMCFM_NEW;
    pg1g4_pi2_WW_VAMCFM_ratio_NEW = pg1g4_pi2_WW_VAMCFM_NEW/p0plus_WW_VAMCFM_NEW;
    pg1g1prime2_pi2_WW_VAMCFM_ratio_NEW = pg1g1prime2_pi2_WW_VAMCFM_NEW/p0plus_WW_VAMCFM_NEW;


    mela.resetInputEvent();
    newtree->Fill();
  }


  foutput->WriteTObject(newtree);
  foutput->Close();
  finput->Close();
}

