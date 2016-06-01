// Mod_Parameters
#include <ZZMatrixElement/MELA/interface/TModParameters.hh>
// NNPDF Driver for JHUGen
#include <ZZMatrixElement/MELA/interface/TNNPDFDriver.hh>
// Mod Kinematics
#include <ZZMatrixElement/MELA/interface/TModKinematics.hh>
// JHUGenMELA
#include <ZZMatrixElement/MELA/interface/TModJHUGen.hh>
#include <ZZMatrixElement/MELA/interface/TModJHUGenMELA.hh>
// Higgs + 0 jet
#include <ZZMatrixElement/MELA/interface/TModHiggsMatEl.hh>
#include <ZZMatrixElement/MELA/interface/TModGravitonMatEl.hh>
#include <ZZMatrixElement/MELA/interface/TModZprimeMatEl.hh>
// Higgs + 1/2 jets
#include <ZZMatrixElement/MELA/interface/TModHiggsJJMatEl.hh>
#include <ZZMatrixElement/MELA/interface/TModHiggsJMatEl.hh>
// VH
#include <ZZMatrixElement/MELA/interface/TModVHiggsMatEl.hh>
// ttH
#include <ZZMatrixElement/MELA/interface/TModTTBHMatEl.hh>
#include <ZZMatrixElement/MELA/interface/TUtil.hh>
#include <ZZMatrixElement/MELA/interface/TEvtProb.hh>
#include <ZZMatrixElement/MELA/interface/MELACandidate.h>

#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ nestedclasses;
#pragma link C++ nestedtypedefs;

#pragma link C++ class MELAParticle;
#pragma link C++ class std::vector<MELAParticle*>;
#pragma link C++ class MELATopCandidate;
#pragma link C++ class MELACandidate;
#pragma link C++ class MelaIO;
#pragma link C++ class TVar;

#pragma link C++ namespace TUtil;
#pragma link C++ function TUtil::ConvertVectorFormat;
#pragma link C++ function TUtil::VHiggsMatEl;

#pragma link C++ class TEvtProb;
#pragma link C++ class SpinZeroCouplings;
#pragma link C++ class SpinOneCouplings;
#pragma link C++ class SpinTwoCouplings;



#endif
