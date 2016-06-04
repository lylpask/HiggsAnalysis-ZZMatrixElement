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




void testME_FullMELA_PingWithFourMomenta(int flavor=0){
  int erg_tev=13;
  float mPOLE=125.;
  float wPOLE=4.07e-3;

  TVar::VerbosityLevel verbosity = TVar::DEBUG;
  Mela mela(erg_tev, mPOLE);
  mela.resetMCFM_EWKParameters(1.16639E-05, 1./128., 79.9549392, 91.1876, 0.23119);
  mela.setVerbosity(verbosity);
  mela.setVerbosity(verbosity);
  mela.setCandidateDecayMode(TVar::CandidateDecay_ZZ);

  float GenLep1Id, GenLep2Id, GenLep3Id, GenLep4Id;
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

  for (int ev = 0; ev < nEntries; ev++){
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
    else{
      GenLep1Id=13;
      GenLep2Id=-13;
      GenLep3Id=13;
      GenLep4Id=-13;
    }
    int idOrdered[4] ={ GenLep1Id, GenLep2Id, GenLep3Id, GenLep4Id };
    TLorentzVector pOrdered[4];
    pOrdered[0].SetXYZT(l1_array[ev][1], l1_array[ev][2], l1_array[ev][3], l1_array[ev][0]);
    pOrdered[1].SetXYZT(l2_array[ev][1], l2_array[ev][2], l2_array[ev][3], l2_array[ev][0]);
    pOrdered[2].SetXYZT(l3_array[ev][1], l3_array[ev][2], l3_array[ev][3], l3_array[ev][0]);
    pOrdered[3].SetXYZT(l4_array[ev][1], l4_array[ev][2], l4_array[ev][3], l4_array[ev][0]);
    SimpleParticleCollection_t daughters;
    for (unsigned int idau=0; idau<4; idau++) daughters.push_back(SimpleParticle_t(idOrdered[idau], pOrdered[idau]));
    mela.setInputEvent(&daughters, (SimpleParticleCollection_t*)0, (SimpleParticleCollection_t*)0, false);

    float pVAMCFM_selfDg1;
    mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
    for (int ii = 0; ii < SIZE_HVV; ii++){ for (int jj = 0; jj < 2; jj++)   mela.selfDHzzcoupl[0][ii][jj] = 0; }
    mela.selfDHzzcoupl[0][0][0]=1.0;
    mela.setMelaHiggsWidth(wPOLE);
    mela.setMelaLeptonInterference(TVar::InterfOn);
    mela.computeP_selfDspin0(pVAMCFM_selfDg1, true);

    mela.resetInputEvent();
  }
}
