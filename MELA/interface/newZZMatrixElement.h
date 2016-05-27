#ifndef newZZMatrixElement_newZZMatrixElement_h
#define newZZMatrixElement_newZZMatrixElement_h

#include <vector>
#include <TLorentzVector.h>
#include "ZZMatrixElement/MELA/interface/TVar.hh"
#include "ZZMatrixElement/MELA/interface/TEvtProb.hh"


class  newZZMatrixElement{
public:
  //pathtoHiggsCSandWidth: path to the textfiles of the HiggsCSandWidth package
  newZZMatrixElement(
    const char* pathtoPDFSet,
    int PDFMember,
    const char* pathtoHiggsCSandWidth,
    double ebeam,
    TVar::VerbosityLevel verbosity
    );

  ~newZZMatrixElement(){ /*std::cout << "End of newZZME" << std::endl;*/ };
  /// Compute KD from masses and angles.
  /// The user must ensure that the order of m1/m2 matches the order of theta1/theta2.
  // flavor 1 for 4e, 2 for 4m, 3 for 2e2mu
  void computeXS(
    TVar::Process process_,
    TVar::MatrixElement me_,
    TVar::Production production_,
    float &mevalue
    );

  void computeProdXS_VVHVV(
    TVar::Process process_,
    TVar::MatrixElement me_,
    TVar::Production production_,
    float& mevalue
    );

  void computeProdXS_JJH(TLorentzVector jet1,
    TVar::Process process_,
    TVar::MatrixElement me_,
    TVar::Production production_,
    float &mevalue
    );

  void computeProdXS_JH(TLorentzVector singleJet,
    TVar::Process process_,
    TVar::MatrixElement me_,
    TVar::Production production_,
    float &mevalue
    );

  void computeProdXS_VH(
    TVar::Process process_,
    TVar::MatrixElement me_,
    TVar::Production production_,
    float &mevalue,
    bool includeHiggsDecay=false
    );

  void computeProdXS_ttH(
    TVar::Process process_,
    TVar::MatrixElement me_,
    TVar::Production production_,
    float &mevalue
    int topProcess,
    int topDecay=0
    );

  void set_Process(TVar::Process process_, TVar::MatrixElement me_, TVar::Production production_);
  void set_mHiggs(float mh_, int index);
  void set_wHiggs(float gah_, int index);
  void set_LeptonInterference(TVar::LeptonInterference myLepInterf);
  void set_LHAgrid(const char* path, int pdfmember=0);
  void set_RenFacScaleMode(TVar::EventScaleScheme renormalizationSch, TVar::EventScaleScheme factorizationSch, double ren_sf, double fac_sf);
  void reset_MCFM_EWKParameters(double ext_Gf, double ext_aemmz, double ext_mW, double ext_mZ, double ext_xW, int ext_ewscheme=3);
  void resetPerEvent();

  void set_SpinZeroCouplings(
    double selfDHvvcoupl_freenorm[SIZE_HVV_FREENORM],
    double selfDHqqcoupl[SIZE_HQQ][2],
    double selfDHggcoupl[SIZE_HGG][2],
    double selfDHzzcoupl[SIZE_HVV][2],
    double selfDHwwcoupl[SIZE_HVV][2],
    double selfDHzzLambda_qsq[4][3],
    double selfDHwwLambda_qsq[4][3],
    int selfDHzzCLambda_qsq[3],
    int selfDHwwCLambda_qsq[3],
    bool diffHWW = false
    );
  void set_SpinOneCouplings(
    double selfDZqqcoupl[SIZE_ZQQ][2],
    double selfDZvvcoupl[SIZE_ZVV][2]
    );
  void set_SpinTwoCouplings(
    double selfDGqqcoupl[SIZE_GQQ][2],
    double selfDGggcoupl[SIZE_GGG][2],
    double selfDGvvcoupl[SIZE_GVV][2]
    );

  // Compute four-momenta from angles - not cos(theta) - only 
  std::vector<TLorentzVector> Calculate4Momentum(double Mx, double M1, double M2, double theta, double theta1, double theta2, double Phi1, double Phi);

  // Get-functions
  MelaIO* get_IORecord();

private:
  
  enum{
    nSupportedHiggses=2
  };

  TVar::VerbosityLevel processVerbosity;
  TVar::LeptonInterference processLeptonInterference;
  TVar::Process processModel;
  TVar::MatrixElement processME;
  TVar::Production processProduction;

  double EBEAM;
  double mHiggs[nSupportedHiggses];
  double wHiggs[nSupportedHiggses];
  TEvtProb Xcal2;

  SpinZeroCouplings* selfD_SpinZeroCouplings;
  SpinOneCouplings* selfD_SpinOneCouplings;
  SpinTwoCouplings* selfD_SpinTwoCouplings;

};
#endif
