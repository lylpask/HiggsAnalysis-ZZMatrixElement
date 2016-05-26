// March 28 2011
// S. Jindariani (sergo@fnal.gov)
// Y. Gao (ygao@fnal.gov)
// K. Burkett (burkett@fnal.gov)


#ifndef ZZ_COMMON
#define ZZ_COMMON
#include <TLorentzVector.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH1D.h>
#include <TH2D.h>
#include <string>
#include <vector>
#include <TFile.h>
#include <TF1.h>
#include "TVar.hh"
// Mod_Parameters
#include "TModParameters.hh"
// NNPDF Driver for JHUGen
#include "TNNPDFDriver.hh"
// Mod Kinematics
#include "TModKinematics.hh"
// JHUGenMELA
#include "TModJHUGen.hh"
#include "TModJHUGenMELA.hh"
// Higgs + 0 jet
#include "TModHiggsMatEl.hh"
#include "TModGravitonMatEl.hh"
#include "TModZprimeMatEl.hh"
// Higgs + 1/2 jets
#include "TModHiggsJJMatEl.hh"
#include "TModHiggsJMatEl.hh"
// VH
#include "TModVHiggsMatEl.hh"
// ttH
#include "TModTTBHMatEl.hh"

using namespace std;

// Parameter settings
void SetEwkCouplingParameters();
double InterpretScaleScheme(TVar::Production production, TVar::MatrixElement matrixElement, TVar::EventScaleScheme scheme, TLorentzVector p[mxpart]);
void SetAlphaS(double Q_ren, double Q_fac, double multiplier_ren, double multiplier_fac, int mynloop, int mynflav, string mypartons);
bool MCFM_chooser(TVar::Process process, TVar::Production production, TVar::LeptonInterference leptonInterf, MELACandidate* cand);

// JHUGen-specific wrappers
void InitJHUGenMELA(const char* pathtoPDFSet, int PDFMember);
void SetJHUGenHiggsMassWidth(double MReso, double GaReso);
void SetJHUGenDistinguishWWCouplings(bool doAllow);

// Spin-0 couplings
void SetMCFMSpinZeroVVCouplings(bool useBSM, SpinZeroCouplings* Hcouplings, bool forceZZ);
void SetJHUGenSpinZeroVVCouplings(double Hvvcoupl[SIZE_HVV][2], int Hvvcoupl_cqsq[3], double HvvLambda_qsq[4][3], bool useWWcoupl);
void SetJHUGenSpinZeroGGCouplings(double Hggcoupl[SIZE_HGG][2]);
void SetJHUGenSpinZeroQQCouplings(double Hqqcoupl[SIZE_HQQ][2]);
// Spin-1 couplings
void SetJHUGenSpinOneCouplings(double Zqqcoupl[SIZE_ZQQ][2], double Zvvcoupl[SIZE_ZVV][2]);
// Spin-2 couplings
void SetJHUGenSpinTwoCouplings(double Gacoupl[SIZE_GGG][2], double Gbcoupl[SIZE_GVV][2], double qLeftRightcoupl[SIZE_GQQ][2]);

// PS cuts, unused
bool MCFM_masscuts(double s[][mxpart], TVar::Process process);
bool MCFM_smalls(double s[][mxpart], int npart);

// ME computations
double SumMatrixElementPDF(
  TVar::Process process, TVar::Production production, TVar::MatrixElement matrixElement,
  event_scales_type* event_scales, MelaIO* RcdME,
  double* flux, double EBEAM, double coupling[SIZE_HVV_FREENORM]
  );
double JHUGenMatEl(
  TVar::Process process, TVar::Production production, TVar::MatrixElement matrixElement,
  event_scales_type* event_scales, MelaIO* RcdME,
  double EBEAM
  );
double HJJMatEl(
  TVar::Process process, TVar::Production production, TVar::MatrixElement matrixElement,
  event_scales_type* event_scales,
  MelaIO* RcdME,
  TVar::VerbosityLevel verbosity,
  double EBEAM
  );
double VHiggsMatEl(
  TVar::Process process, TVar::Production production, TVar::MatrixElement matrixElement,
  event_scales_type* event_scales, MelaIO* RcdME,
  TVar::VerbosityLevel verbosity,
  double EBEAM
  );
double TTHiggsMatEl(
  TVar::Process process, TVar::Production production, TVar::MatrixElement matrixElement,
  event_scales_type* event_scales, MelaIO* RcdME,
  TVar::VerbosityLevel verbosity,
  double EBEAM,
  int topDecay, int topProcess
  );
double BBHiggsMatEl(
  TVar::Process process, TVar::Production production, TVar::MatrixElement matrixElement,
  event_scales_type* event_scales, MelaIO* RcdME,
  TVar::VerbosityLevel verbosity,
  double EBEAM,
  int topProcess
  );

void Swap_Momenta(double(&p)[4], double(&q)[4]);
bool CheckPartonMomFraction(const TLorentzVector p0, const TLorentzVector p1, double xx[2], TVar::VerbosityLevel verbosity, double EBEAM);
void ComputePDF(const TLorentzVector p0, const TLorentzVector p1, double fx1[nmsq], double fx2[nmsq], TVar::VerbosityLevel verbosity, double EBEAM);
double SumMEPDF(const TLorentzVector p0, const TLorentzVector p1, double msq[nmsq][nmsq], MelaIO* RcdME, TVar::VerbosityLevel verbosity, double EBEAM);

// Boost the particles with or without associated ones to pT=0 frame and return std::vectors filled with (id, momentum) pairs
void GetBoostedParticleVectors(
  MELACandidate* melaCand,
  simple_event_record& mela_event
  );


#endif
