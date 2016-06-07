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




void testME_Dec_MCFM_Ping(int flavor=0){
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
    else if (flavor == 4){
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
    mela.setInputEvent(&daughters_WWasZZ, (SimpleParticleCollection_t*)0, (SimpleParticleCollection_t*)0, false);
    mela.setCandidateDecayMode(TVar::CandidateDecay_ZZ);
    mela.setInputEvent(&daughters_ZZ, (SimpleParticleCollection_t*)0, (SimpleParticleCollection_t*)0, false);
    mela.setInputEvent(&daughters_ZG, (SimpleParticleCollection_t*)0, (SimpleParticleCollection_t*)0, false);
    mela.setCandidateDecayMode(TVar::CandidateDecay_GG);
    mela.setInputEvent(&daughters_GG, (SimpleParticleCollection_t*)0, (SimpleParticleCollection_t*)0, false);

    int cindex;

    /***** WW *****/
    cindex=0;
    mela.setCurrentCandidateFromIndex(cindex);
    float pVAMCFM_qqWW_bkg;
    mela.setProcess(TVar::bkgWW, TVar::MCFM, TVar::ZZQQB);
    mela.computeP(pVAMCFM_qqWW_bkg, true);
    cout << "pVAMCFM_qqWW_bkg: " << pVAMCFM_qqWW_bkg << '\n' << endl;

    float pVAMCFM_ggVV_total;
    mela.setProcess(TVar::bkgWWZZ_SMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.computeP(pVAMCFM_ggVV_total, true);
    cout << "pVAMCFM_ggVV_total: " << pVAMCFM_ggVV_total << '\n' << endl;
    float pVAMCFM_ggVV_bkg;
    mela.setProcess(TVar::bkgWWZZ, TVar::MCFM, TVar::ZZGG);
    mela.computeP(pVAMCFM_ggVV_bkg, true);
    cout << "pVAMCFM_ggVV_bkg: " << pVAMCFM_ggVV_bkg << '\n' << endl;
    float pVAMCFM_ggVV_sig;
    mela.setProcess(TVar::HSMHiggs_WWZZ, TVar::MCFM, TVar::ZZGG);
    mela.computeP(pVAMCFM_ggVV_sig, true);
    cout << "pVAMCFM_ggVV_sig: " << pVAMCFM_ggVV_sig << '\n' << endl;

    float pVAMCFM_ggWW_total;
    mela.setProcess(TVar::bkgWWZZ_SMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.computeP(pVAMCFM_ggWW_total, true);
    cout << "pVAMCFM_ggWW_total: " << pVAMCFM_ggWW_total << '\n' << endl;
    float pVAMCFM_ggWW_bkg;
    mela.setProcess(TVar::bkgWWZZ, TVar::MCFM, TVar::ZZGG);
    mela.computeP(pVAMCFM_ggWW_bkg, true);
    cout << "pVAMCFM_ggWW_bkg: " << pVAMCFM_ggWW_bkg << '\n' << endl;
    float pVAMCFM_ggWW_sig;
    mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.computeP(pVAMCFM_ggWW_sig, true);
    cout << "pVAMCFM_ggWW_sig: " << pVAMCFM_ggWW_sig << '\n' << endl;

    /***** WW (as ZZ) *****/
    cindex=1;
    mela.setCurrentCandidateFromIndex(cindex);
    float pVAMCFM_ggZZ_total;
    mela.setProcess(TVar::bkgZZ_SMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.computeP(pVAMCFM_ggZZ_total, true);
    cout << "pVAMCFM_ggZZ_total: " << pVAMCFM_ggZZ_total << '\n' << endl;
    float pVAMCFM_ggZZ_bkg;
    mela.setProcess(TVar::bkgZZ, TVar::MCFM, TVar::ZZGG);
    mela.computeP(pVAMCFM_ggZZ_bkg, true);
    cout << "pVAMCFM_ggZZ_bkg: " << pVAMCFM_ggZZ_bkg << '\n' << endl;
    float pVAMCFM_ggWWasZZ_sig;
    mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.computeP(pVAMCFM_ggWWasZZ_sig, true);
    cout << "pVAMCFM_ggWWasZZ_sig: " << pVAMCFM_ggWWasZZ_sig << '\n' << endl;

    /***** ZZ *****/
    cindex=2;
    mela.setCurrentCandidateFromIndex(cindex);
    float pVAMCFM_ggZZ_sig;
    mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
    for (int ii = 0; ii < SIZE_HVV; ii++){ for (int jj = 0; jj < 2; jj++)   mela.selfDHzzcoupl[0][ii][jj] = 0; }
    mela.setMelaHiggsWidth(wPOLE);
    mela.setMelaLeptonInterference(TVar::InterfOn);
    mela.computeP(pVAMCFM_ggZZ_sig, true);
    cout << "pVAMCFM_ggZZ_sig: " << pVAMCFM_ggZZ_sig << '\n' << endl;

    float pVAMCFM_ZZ_sig_selfDg1;
    mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
    for (int ii = 0; ii < SIZE_HVV; ii++){ for (int jj = 0; jj < 2; jj++)   mela.selfDHzzcoupl[0][ii][jj] = 0; }
    mela.selfDHzzcoupl[0][0][0]=1.0;
    mela.setMelaHiggsWidth(wPOLE);
    mela.setMelaLeptonInterference(TVar::InterfOn);
    mela.computeP(pVAMCFM_ZZ_sig_selfDg1, true);
    cout << "pVAMCFM_ggZZ_sig_selfDg1: " << pVAMCFM_ZZ_sig_selfDg1 << '\n' << endl;

    /***** ZG *****/
    cindex=3;
    mela.setCurrentCandidateFromIndex(cindex);
    float pVAMCFM_ZG_bkg;
    mela.setProcess(TVar::bkgZGamma, TVar::MCFM, TVar::ZZQQB);
    mela.computeP(pVAMCFM_ZG_bkg, true);
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
void testME_ProdDec_MCFM_Ping(int flavor=0){
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
  else if (flavor == 4){
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
