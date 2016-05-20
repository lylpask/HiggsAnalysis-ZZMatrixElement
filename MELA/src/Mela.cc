/*

************* CMS MELA interface to MCFM/JHUGen-MELA *************
Last update: 19.05.2016 by U. Sarica

Notes:
1) Each specific type of computeP* function comes with its wrapper for common use.
   Removing these wrappers from Mela will not introduce any bugs, but
   it might affect packages that depend on it (eg. MEMCalculators).
2) ...

Please adhere to the following coding conventions:
1) Never put return statements in the middle of the computeP* functions unless it is a wrapper.
   Functions calling the Mela::ZZME member have to reset the couplings, so an abrupt termination
   does not reset the couplings properly.
2) ...

*/

#include <ZZMatrixElement/MELA/interface/Mela.h>
#include <ZZMatrixElement/MELA/interface/newZZMatrixElement.h>
#include <DataFormats/GeometryVector/interface/Pi.h>
#include <FWCore/ParameterSet/interface/FileInPath.h>

#include <ZZMatrixElement/MELA/interface/computeAngles.h>
#include <ZZMatrixElement/MELA/interface/VectorPdfFactory.h>
#include <ZZMatrixElement/MELA/interface/TensorPdfFactory.h>
#include <ZZMatrixElement/MELA/interface/RooqqZZ_JHU_ZgammaZZ_fast.h>
#include <ZZMatrixElement/MELA/interface/RooqqZZ_JHU.h>
#include <ZZMatrixElement/MELA/interface/SuperMELA.h>

#include <RooMsgService.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TGraph.h>
#include <TSpline.h>
#include <vector>

#include <string>
#include <cstdio>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>


using namespace RooFit;

Mela::Mela(int LHCsqrts, float mh) 
{
  const int maxSqrts = 8;
  setRemoveLeptonMasses(false); // Use Run I scheme for not removing lepton masses

  // Create symlinks to the required files, if these are not already present (do nothing otherwse)
  edm::FileInPath mcfmInput1("ZZMatrixElement/MELA/data/input.DAT");
  edm::FileInPath mcfmInput2("ZZMatrixElement/MELA/data/process.DAT");
  edm::FileInPath mcfmInput3("ZZMatrixElement/MELA/data/Pdfdata/cteq6l1.tbl");  
  edm::FileInPath mcfmInput4("ZZMatrixElement/MELA/data/Pdfdata/cteq6l.tbl");  
  edm::FileInPath mcfmWarning("ZZMatrixElement/MELA/data/ffwarn.dat");
  edm::FileInPath mcfm_brsm_o("ZZMatrixElement/MELA/data/br.sm1");
  edm::FileInPath mcfm_brsm_t("ZZMatrixElement/MELA/data/br.sm2");
  symlink(mcfmWarning.fullPath().c_str(), "ffwarn.dat");
  symlink(mcfm_brsm_o.fullPath().c_str(), "br.sm1");
  symlink(mcfm_brsm_t.fullPath().c_str(), "br.sm2");
  symlink(mcfmInput1.fullPath().c_str(), "input.DAT");
  symlink(mcfmInput2.fullPath().c_str(), "process.DAT");
  mkdir("Pdfdata",S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  symlink(mcfmInput3.fullPath().c_str(), "Pdfdata/cteq6l1.tbl");
  symlink(mcfmInput4.fullPath().c_str(), "Pdfdata/cteq6l.tbl");

  mzz_rrv = new RooRealVar("mzz","m_{ZZ}",(double)mh, 0.,1000.);
  z1mass_rrv = new RooRealVar("z1mass","m_{Z1}",0.,160.);
  z2mass_rrv = new RooRealVar("z2mass","m_{Z2}",0.,200.); 
  costhetastar_rrv = new RooRealVar("costhetastar","cos#theta^{*}",-1.,1.);  
  costheta1_rrv = new RooRealVar("costheta1","cos#theta_{1}",-1.,1.);  
  costheta2_rrv = new RooRealVar("costheta2","cos#theta_{2}",-1.,1.);
  phi_rrv= new RooRealVar("phi", "#Phi", -TMath::Pi(), TMath::Pi());
  phi1_rrv= new RooRealVar("phi1", "#Phi_{1}", -TMath::Pi(), TMath::Pi());
  Y_rrv = new RooRealVar("Yzz", "#Y_{ZZ}", 0, -4, 4);

  upFrac_rrv = new RooRealVar("upFrac", "fraction up-quarks", .5, 0., 1.);

  RooSpinZero::modelMeasurables measurables_;
  measurables_.h1 = costheta1_rrv;
  measurables_.h2 = costheta2_rrv;
  measurables_.Phi = phi_rrv;
  measurables_.m1 = z1mass_rrv;
  measurables_.m2 = z2mass_rrv;
  measurables_.m12 = mzz_rrv;
  measurables_.hs = costhetastar_rrv;
  measurables_.Phi1 = phi1_rrv;
  measurables_.Y = Y_rrv;

  ggSpin0Model = new ScalarPdfFactory_ggH(measurables_, false, 1, 1); // 1,1==ZZ
  spin1Model = new VectorPdfFactory(z1mass_rrv,z2mass_rrv,costhetastar_rrv,costheta1_rrv,costheta2_rrv,phi_rrv,phi1_rrv,mzz_rrv);
  spin2Model = new TensorPdfFactory_HVV(measurables_, RooSpin::kVdecayType_Zll, RooSpin::kVdecayType_Zll);
  qqZZmodel = new RooqqZZ_JHU_ZgammaZZ_fast("qqZZmodel","qqZZmodel",*z1mass_rrv,*z2mass_rrv,*costheta1_rrv,*costheta2_rrv,*phi_rrv,*costhetastar_rrv,*phi1_rrv,*mzz_rrv,*upFrac_rrv);

 //edm::FileInPath HiggsWidthFile("Higgs/Higgs_CS_and_Width/txtFiles/8TeV-ggH.txt");
// edm::FileInPath HiggsWidthFile("Higgs/Higgs_CS_and_Width/txtFiles/YR3/8TeV-ggH.txt");
  edm::FileInPath HiggsWidthFile("ZZMatrixElement/MELA/data/HiggsTotalWidth.txt");
  std::string path = HiggsWidthFile.fullPath();
  //std::cout << path.substr(0,path.length()-12) << std::endl;
 // ZZME = new  newZZMatrixElement(path.substr(0,path.length()-12 ).c_str(),1000.*LHCsqrts/2.,TVar::INFO);
  edm::FileInPath path_nnpdf("ZZMatrixElement/MELA/data/Pdfdata/NNPDF30_lo_as_0130.LHgrid");
  char path_nnpdf_c[] = "Pdfdata/NNPDF30_lo_as_0130.LHgrid";
  int pdfmember = 0;
  symlink(path_nnpdf.fullPath().c_str(), path_nnpdf_c);
  ZZME = new  newZZMatrixElement(path_nnpdf_c, pdfmember, path.substr(0, path.length()-19).c_str(), 1000.*LHCsqrts/2., TVar::INFO);
  ZZME->set_mHiggs(mh);
  ZZME->set_wHiggs(-1);
  ZZME->set_LeptonInterference(TVar::DefaultLeptonInterf);

  // 
  // configure the JHUGEn and MCFM calculations 
  // 
  // load TGraphs for VAMCFM scale factors
  edm::FileInPath ScaleFactorFile("ZZMatrixElement/MELA/data/scaleFactors.root");
  TFile* sf = TFile::Open(ScaleFactorFile.fullPath().c_str(),"r");
  vaScale_4e    = (TGraph*)sf->Get("scaleFactors_4e");
  vaScale_4mu   = (TGraph*)sf->Get("scaleFactors_4mu");
  vaScale_2e2mu = (TGraph*)sf->Get("scaleFactors_2e2mu");
  sf->Close(); 
  edm::FileInPath DggScaleFactorFile("ZZMatrixElement/MELA/data/scalefactor_DggZZ.root");
  //cout << DggScaleFactorFile.fullPath().c_str() << endl;
  TFile* af = new TFile(DggScaleFactorFile.fullPath().c_str(),"r");
  DggZZ_scalefactor = (TGraph*) af->Get("scalefactor");
  af->Close();
  assert(vaScale_4e);
  assert(vaScale_4mu);
  assert(vaScale_2e2mu);
  assert(DggZZ_scalefactor);
  
  //
  // setup supermela
  //

  //deactivate generation messages
  RooMsgService::instance().getStream(1).removeTopic(NumIntegration);
  RooMsgService::instance().setStreamStatus(1,kFALSE);
  RooMsgService::instance().setStreamStatus(0,kFALSE);// silence also the error messages, but should really be looked at.

  myR=new TRandom3(35797);
  //  std::cout << "before supermela" << std::endl;
  
  int superMELA_LHCsqrts = LHCsqrts;
  if (superMELA_LHCsqrts > maxSqrts) superMELA_LHCsqrts = maxSqrts;
  super = new SuperMELA(mh, "4mu", superMELA_LHCsqrts); // preliminary intialization, we adjust the flavor later
  char cardpath[500];
  sprintf(cardpath, "ZZMatrixElement/MELA/data/CombinationInputs/SM_inputs_%dTeV/inputs_4mu.txt", superMELA_LHCsqrts);
  //std::cout << "before supermela, pathToCards: " <<cardpath<< std::endl;
  edm::FileInPath cardfile(cardpath);
  std::string cpath=cardfile.fullPath();
  //std::cout << cpath.substr(0,cpath.length()-14).c_str()  <<std::endl;
  super->SetPathToCards(cpath.substr(0,cpath.length()-14).c_str() );
  super->SetVerbosity(false);
  // std::cout << "starting superMELA initialization" << std::endl;
  super->init();
  //std::cout << "after supermela" << std::endl;
 	edm::FileInPath CTotBkgFile("ZZMatrixElement/MELA/data/ZZ4l-C_TotalBkgM4lGraph.root");
	TFile* finput_ctotbkg = TFile::Open(CTotBkgFile.fullPath().c_str(),"read");
	for (int i=0;i<3;i++){
		tgtotalbkg[i] = 0;
	} 
	setCTotalBkgGraphs(finput_ctotbkg, tgtotalbkg);
	finput_ctotbkg->Close(); 
  for (int i=0; i<3; i++){
    assert(tgtotalbkg[i]);
  }

  reset_SelfDCouplings();
}

Mela::~Mela(){
//  std::cout << "begin destructor" << std::endl;  
  setRemoveLeptonMasses(false); // Use Run I scheme for not removing lepton masses. Notice the switch itself is defined as an extern, so it has to be set to default value at the destructor!

  for (int i=0; i<3; i++){
    if (tgtotalbkg[i] != 0) delete tgtotalbkg[i];
  }
  if (DggZZ_scalefactor!=0) delete DggZZ_scalefactor;

  // Let the derived RooFit objects come first...
  delete ggSpin0Model;
  delete spin1Model;
  delete spin2Model;
  delete qqZZmodel;
  // ...then delete the observables.
  delete mzz_rrv;
  delete z1mass_rrv; 
  delete z2mass_rrv; 
  delete costhetastar_rrv;
  delete costheta1_rrv;
  delete costheta2_rrv;
  delete phi_rrv;
  delete phi1_rrv;
  delete Y_rrv;
  delete upFrac_rrv;

  delete ZZME;
  delete super;
  delete myR;

//  std::cout << "end destructor" << std::endl;
}

void Mela::setProcess(TVar::Process myModel, TVar::MatrixElement myME, TVar::Production myProduction)
{
  myModel_ = myModel;
  myME_ = myME;
  myProduction_ = myProduction;
}

void Mela::setMelaHiggsWidth(float myHiggsWidth){ // Should be called per-event
	ZZME->set_wHiggs(myHiggsWidth);
}

void Mela::setMelaLeptonInterference(TVar::LeptonInterference myLepInterf){
	ZZME->set_LeptonInterference(myLepInterf);
}

void Mela::resetMCFM_EWKParameters(double ext_Gf, double ext_aemmz, double ext_mW, double ext_mZ, double ext_xW, int ext_ewscheme){
  ZZME->reset_MCFM_EWKParameters(ext_Gf, ext_aemmz, ext_mW, ext_mZ, ext_xW, ext_ewscheme);
}

void Mela::setRemoveLeptonMasses(bool MasslessLeptonSwitch){ mela::applyLeptonMassCorrection(MasslessLeptonSwitch); }

void Mela::setRenFacScaleMode(TVar::EventScaleScheme renormalizationSch, TVar::EventScaleScheme factorizationSch, double ren_sf, double fac_sf){
  ZZME->set_RenFacScaleMode(renormalizationSch, factorizationSch, ren_sf, fac_sf);
}

// Re-order masses and angles as needed by likelihoodDiscriminant. 
// This follows a different convention than the usual Z1/Z2 definition!
void Mela::checkZorder(
  float& z1mass, float& z2mass,
  float& costhetastar, float& costheta1,
  float& costheta2, float& phi,
  float& phistar1
  ){

  float tempZ1mass=z1mass;
  float tempZ2mass=z2mass;
  float tempH1=costheta1;
  float tempH2=costheta2;
  float tempHs=costhetastar;
  float tempPhi1=phistar1;
  float tempPhi=phi;

  if (z2mass>z1mass){
    z1mass=tempZ2mass;
    z2mass=tempZ1mass;
    costhetastar=-tempHs;
    costheta1=tempH2;
    costheta2=tempH1;
    phi=tempPhi;
    phistar1=-tempPhi1-tempPhi;
    if (phistar1>3.141592) phistar1=phistar1-2*Geom::pi();
    else if (phistar1<-3.141592) phistar1=phistar1+2*Geom::pi();
  }
  else return;
}

std::vector<TLorentzVector> Mela::calculate4Momentum(double Mx, double M1, double M2, double theta, double theta1, double theta2, double Phi1, double Phi){
	return ZZME->Calculate4Momentum(Mx, M1, M2, theta, theta1, theta2, Phi1, Phi);
}


void Mela::computeD_CP(
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
  ){
  reset_PAux();

  double coupl_mix[SIZE_HVV][2] ={ { 0 } };
  double coupl_1[SIZE_HVV][2] ={ { 0 } };
  double coupl_2[SIZE_HVV][2] ={ { 0 } };

  switch (myType){
  case TVar::D_g1g4:
    coupl_mix[0][0] =1.;
    coupl_mix[3][0] =2.521;
    coupl_1[0][0] =1.;
    coupl_2[3][0] =2.521;
    break;
  case TVar::D_g1g4_pi_2:
    coupl_mix[0][0] =1.;
    coupl_mix[3][1] =2.521;
    coupl_1[0][0] =1.;
    coupl_2[3][1] =2.521;
    break;
  case TVar::D_g1g2:
    coupl_mix[0][0] =1.;
    coupl_mix[1][0] = 1.638;
    coupl_1[0][0] =1.;
    coupl_2[1][0] = 1.638;
    break;
  case TVar::D_g1g2_pi_2:
    coupl_mix[0][0] =1.;
    coupl_mix[1][1] = 1.638;
    coupl_1[0][0] =1.;
    coupl_2[1][1] = 1.638;
    break;
  case TVar::D_g1g1prime2:
    coupl_mix[0][0] =1.;
    coupl_mix[11][0] = 12046.01;
    coupl_1[0][0] =1.;
    coupl_2[11][0] = 12046.01;
    break;
  case TVar::D_zzzg:
    coupl_mix[0][0] =1.;
    coupl_mix[4][0] = 0.0688;
    coupl_1[0][0] =1.;
    coupl_2[4][0] = 0.0688;
    break;
  case TVar::D_zzgg:
    coupl_mix[0][0] =1.;
    coupl_mix[7][0] = -0.0898;
    coupl_1[0][0] =1.;
    coupl_2[7][0] = -0.0898;
    break;
  case TVar::D_zzzg_PS:
    coupl_mix[0][0] =1.;
    coupl_mix[6][0] = 0.0855;
    coupl_1[0][0] =1.;
    coupl_2[6][0] = 0.0855;
    break;
  case TVar::D_zzgg_PS:
    coupl_mix[0][0] =1.;
    coupl_mix[9][0] = -0.0907;
    coupl_1[0][0] =1.;
    coupl_2[9][0] = -0.0907;
    break;
  case TVar::D_zzzg_g1prime2:
    coupl_mix[0][0] =1.;
    coupl_mix[30][0] = -7591.914;
    coupl_1[0][0] =1.;
    coupl_2[30][0] = -7591.914;
    break;
  case TVar::D_zzzg_g1prime2_pi_2:
    coupl_mix[0][0] =1.;
    coupl_mix[30][1] = -7591.914;
    coupl_1[0][0] =1.;
    coupl_2[30][1] = -7591.914;
    break;
  default:
    cout <<"Error: Not supported!"<<endl;
  }
  float pMix, p1, p2;
  setProcess(TVar::SelfDefine_spin0, myME, TVar::ZZGG);
  computeP(mZZ, mZ1, mZ2,
    costhetastar, costheta1, costheta2, phi, phi1, flavor, coupl_mix, pMix);

  setProcess(TVar::SelfDefine_spin0, myME, TVar::ZZGG);
  computeP(mZZ, mZ1, mZ2,
    costhetastar, costheta1, costheta2, phi, phi1, flavor, coupl_1, p1);

  setProcess(TVar::SelfDefine_spin0, myME, TVar::ZZGG);
  computeP(mZZ, mZ1, mZ2,
    costhetastar, costheta1, costheta2, phi, phi1, flavor, coupl_2, p2);

  prob = pMix- p1- p2;
}

 
void Mela::computeP(
  float mZZ, float mZ1, float mZ2, // input kinematics
  float costhetastar,
  float costheta1,
  float costheta2,
  float phi,
  float phi1,
  int flavor,
  float& prob, // output probability
  bool useConstant
  ){
  reset_PAux();
  float constant = 1;


  if (myME_ == TVar::ANALYTICAL) {
    costhetastar_rrv->setVal(costhetastar);
    costheta1_rrv->setVal(costheta1);
    costheta2_rrv->setVal(costheta2);
    phi_rrv->setVal(phi);
    phi1_rrv->setVal(phi1);
    z1mass_rrv->setVal(mZ1);
    z2mass_rrv->setVal(mZ2);
    mzz_rrv->setVal(mZZ);

    configureAnalyticalPDFs();

    if (myProduction_==TVar::ZZINDEPENDENT){
      RooAbsPdf* integral = (RooAbsPdf*)pdf->createIntegral(RooArgSet(*costhetastar_rrv, *phi1_rrv));
      prob = integral->getVal();
      delete integral;
    }
    else prob = pdf->getVal();
  }
  else if (myME_ == TVar::JHUGen || myME_ == TVar::MCFM) {

    //initialize variables
    checkZorder(mZ1, mZ2, costhetastar, costheta1, costheta2, phi, phi1);
    ZZME->computeXS(
      mZZ, mZ1, mZ2,
      costhetastar, costheta1, costheta2,
      phi, phi1, flavor,
      myModel_, myME_, myProduction_,
      prob
      );

    if (myME_ == TVar::MCFM  && myProduction_ == TVar::ZZINDEPENDENT &&  myModel_ == TVar::bkgZZ){
      prob = 0;
      int gridsize_hs = 5;
      double hs_min = 0; //-1.;
      double hs_max = 1;
      double hs_step = (hs_max - hs_min) / double(gridsize_hs);

      int gridsize_phi1 = 5;
      double phi1_min = 0; //-TMath::Pi();
      double phi1_max = TMath::Pi();
      double phi1_step = (phi1_max - phi1_min) / double(gridsize_phi1);

      for (int i_hs = 0; i_hs < gridsize_hs + 1; i_hs++) {
        double hs_val = hs_min + i_hs * hs_step;
        for (int i_phi1 = 0; i_phi1 < gridsize_phi1 +1; i_phi1++) {
          double phi1_val = phi1_min + i_phi1 * phi1_step;
          float temp_prob(0.);
          // calculate the ZZ using MCFM
          ZZME->computeXS(
            mZZ, mZ1, mZ2,
            hs_val, costheta1, costheta2,
            phi, phi1_val, flavor,
            myModel_, myME_, myProduction_,
            temp_prob
            );
          prob += temp_prob;
        }
      }
      prob =  prob / float((gridsize_hs + 1) * (gridsize_phi1 +1));
    }
  }

  if (useConstant){
    constant = getConstant(flavor, mZZ, false);
    prob *= constant;
  }

  reset_SelfDCouplings();
}

float Mela::getConstant(int flavor, float mZZ, bool useOldggZZConstants){
  float constant = 1;

  if (myME_ == TVar::ANALYTICAL) {

    // gg productions 
    if (myProduction_ == TVar::ZZGG) {
      if (flavor == 3) {
        //cout << "ANALYTICAL - GG - flavor=3" << endl;
        if (myModel_ == TVar::H0minus)  constant = 6.4;
        if (myModel_ == TVar::H0hplus)  constant = 2.2;
        if (myModel_ == TVar::H2_g1g5)  constant = 9.5;
        if (myModel_ == TVar::H2_g1)  constant = 5.5;
        if (myModel_ == TVar::H2_g4)  constant = 7.3e7;
        if (myModel_ == TVar::H2_g8)  constant = 1.1e8;
        if (myModel_ == TVar::H2_g5)  constant = 16.3;
        if (myModel_ == TVar::H2_g2) constant = 552385;
        if (myModel_ == TVar::H2_g3) constant = 1.08147e+08;
        if (myModel_ == TVar::H2_g6) constant = 748630;
        if (myModel_ == TVar::H2_g7) constant = 1.2522e+07;
        if (myModel_ == TVar::H2_g9) constant = 3.01467e+06;
        if (myModel_ == TVar::H2_g10) constant = 1.483e+11;
      }
      else{
        //cout << "ANALYTICAL - GG - flavor!=3" << endl;
        if (myModel_ == TVar::H0minus)  constant = 6.5;
        if (myModel_ == TVar::H0hplus)  constant = 2.2;
        if (myModel_ == TVar::H2_g1g5)  constant = 9.3;
        if (myModel_ == TVar::H2_g1)  constant = 5.5;
        if (myModel_ == TVar::H2_g4)  constant = 1.1e8;
        if (myModel_ == TVar::H2_g8)  constant = 1.9e8;
        if (myModel_ == TVar::H2_g5)  constant = 15.6;
        if (myModel_ == TVar::H2_g2) constant = 552385;
        if (myModel_ == TVar::H2_g3) constant = 1.24897e+08;
        if (myModel_ == TVar::H2_g6) constant = 1.10793e+06;
        if (myModel_ == TVar::H2_g7) constant = 2.21423e+07;
        if (myModel_ == TVar::H2_g9) constant = 3.18193e+06;
        if (myModel_ == TVar::H2_g10) constant = 2.63811e+11;
      }
    }
    // qqb productions 
    if (myProduction_ == TVar::ZZQQB) {
      if (flavor == 3) {
        //cout << "ANALYTICAL - QQB - flavor=3" << endl;
        if (myModel_ == TVar::H1minus)  constant = 4.6e5;
        if (myModel_ == TVar::H1plus)  constant = 4.0e5;
        if (myModel_ == TVar::H2_g1g5)  constant = 7.9;
        if (myModel_ == TVar::H2_g1)  constant = 4.5;
        if (myModel_ == TVar::H2_g5) constant = 13.7977;
        if (myModel_ == TVar::H2_g4) constant = 5.12897e+07;
        if (myModel_ == TVar::H2_g2) constant = 477586;
        if (myModel_ == TVar::H2_g3) constant = 1.30907e+08;
        if (myModel_ == TVar::H2_g6) constant = 847461;
        if (myModel_ == TVar::H2_g7) constant = 1.39014e+07;
        if (myModel_ == TVar::H2_g8) constant = 7.08446e+07;
        if (myModel_ == TVar::H2_g9) constant = 2.93583e+06;
        if (myModel_ == TVar::H2_g10) constant = 1.47118e+11;
      }
      else{
        //cout << "ANALYTICAL - QQB - flavor!=3" << endl;
        if (myModel_ == TVar::H1minus)  constant = 4.6e5;
        if (myModel_ == TVar::H1plus)  constant = 4.0e5;
        if (myModel_ == TVar::H2_g1g5)  constant = 7.9;
        if (myModel_ == TVar::H2_g1)  constant = 4.5;
        if (myModel_ == TVar::H2_g5) constant = 13.7289;
        if (myModel_ == TVar::H2_g4) constant = 7.57539e+07;
        if (myModel_ == TVar::H2_g2) constant = 476156;
        if (myModel_ == TVar::H2_g3) constant = 1.44675e+08;
        if (myModel_ == TVar::H2_g6) constant = 1.07303e+06;
        if (myModel_ == TVar::H2_g7) constant = 2.37359e+07;
        if (myModel_ == TVar::H2_g8) constant = 1.35435e+08;
        if (myModel_ == TVar::H2_g9) constant = 2.99514e+06;
        if (myModel_ == TVar::H2_g10) constant = 2.13201e+11;
      }
    }
    // production independent calculations
    if (myProduction_ == TVar::ZZINDEPENDENT) {
      if (flavor == 3) {
        //cout << "ANALYTICAL - INDEP - flavor=3" << endl;
        if (myModel_ == TVar::H1minus)  constant = 3.4e4;
        if (myModel_ == TVar::H1plus)  constant = 3.4e4;
        if (myModel_ == TVar::H2_g1g5)  constant = 0.66;
        if (myModel_ == TVar::H2_g5) constant = 1.15604;
        if (myModel_ == TVar::H2_g4) constant = 4.36662e+06;
        if (myModel_ == TVar::H2_g2) constant = 39994.6;
        if (myModel_ == TVar::H2_g3) constant = 1.0897e+07;
        if (myModel_ == TVar::H2_g6) constant = 61420.6;
        if (myModel_ == TVar::H2_g7) constant = 1.20742e+06;
        if (myModel_ == TVar::H2_g8) constant = 6.07991e+06;
        if (myModel_ == TVar::H2_g9) constant = 239187;
        if (myModel_ == TVar::H2_g10) constant = 1.1843e+10;
      }
      else{
        //cout << "ANALYTICAL - INDEP - flavor!=3" << endl;
        if (myModel_ == TVar::H1minus)  constant = 3.4e4;
        if (myModel_ == TVar::H1plus)  constant = 3.4e4;
        if (myModel_ == TVar::H2_g1g5)  constant = .66;
        if (myModel_ == TVar::H2_g5) constant = 1.15604;
        if (myModel_ == TVar::H2_g4) constant = 5.6237e+06;
        if (myModel_ == TVar::H2_g2) constant = 39715.6;
        if (myModel_ == TVar::H2_g3) constant = 1.16172e+07;
        if (myModel_ == TVar::H2_g6) constant = 77613.8;
        if (myModel_ == TVar::H2_g7) constant = 1.58485e+06;
        if (myModel_ == TVar::H2_g8) constant = 8.71451e+06;
        if (myModel_ == TVar::H2_g9) constant = 241591;
        if (myModel_ == TVar::H2_g10) constant = 1.55139e+10;
      }
    }
  }
  else if (myME_ == TVar::JHUGen) {

    if (myProduction_ == TVar::ZZGG) {
      if (flavor == 3) {
        if (myModel_ == TVar::H0minus)  constant = 6.0;
        if (myModel_ == TVar::H0hplus)  constant = 2.1;
        if (myModel_ == TVar::H2_g1g5)  constant = 0.6;
        if (myModel_ == TVar::H2_g4)  constant = 2.7e10;
        if (myModel_ == TVar::H2_g8)  constant = 4.1e10;
        if (myModel_ == TVar::H2_g5)  constant = .97;
        if (myModel_ == TVar::H2_g2) constant = 4.74608e+08;
        if (myModel_ == TVar::H2_g3) constant = 1.63324e+11;
        if (myModel_ == TVar::H2_g6) constant = 46816.9;
        if (myModel_ == TVar::H2_g7) constant = 783088;
        if (myModel_ == TVar::H2_g9) constant = 1.13853e+09;
        if (myModel_ == TVar::H2_g10) constant = 5.58394e+13;
      }
      else{
        if (myModel_ == TVar::H0minus)  constant = 7.0;
        if (myModel_ == TVar::H0hplus)  constant = 2.3;
        if (myModel_ == TVar::H2_g1g5)  constant = 0.7;
        if (myModel_ == TVar::H2_g4)  constant = 2.6e10;
        if (myModel_ == TVar::H2_g8)  constant = 3.7e10;
        if (myModel_ == TVar::H2_g5)  constant = 1.26;
        if (myModel_ == TVar::H2_g2) constant = 5.96148e+08;
        if (myModel_ == TVar::H2_g3) constant = 1.95534e+11;
        if (myModel_ == TVar::H2_g6) constant = 55160.4;
        if (myModel_ == TVar::H2_g7) constant = 658026;
        if (myModel_ == TVar::H2_g9) constant = 1.23089e+09;
        if (myModel_ == TVar::H2_g10) constant = 5.78862e+13;
      }
    }
    // qqb productions 
    if (myProduction_ == TVar::ZZQQB) {
      if (flavor == 3) {
        if (myModel_ == TVar::H1minus)  constant = 16.;
        if (myModel_ == TVar::H1plus)  constant = 13.;
        if (myModel_ == TVar::H2_g1g5)  constant = 13.;
        if (myModel_ == TVar::H2_g5) constant = 22.8625;
        if (myModel_ == TVar::H2_g4) constant = 8.49013e+07;
        if (myModel_ == TVar::H2_g2) constant = 792938;
        if (myModel_ == TVar::H2_g3) constant = 2.17563e+08;
        if (myModel_ == TVar::H2_g6) constant = 1.40845e+06;
        if (myModel_ == TVar::H2_g7) constant = 2.31499e+07;
        if (myModel_ == TVar::H2_g8) constant = 1.17271e+08;
        if (myModel_ == TVar::H2_g9) constant = 4.86462e+06;
        if (myModel_ == TVar::H2_g10) constant = 2.43529e+11;
      }
      else{
        if (myModel_ == TVar::H1minus)  constant = 19.;
        if (myModel_ == TVar::H1plus)  constant = 14.;
        if (myModel_ == TVar::H2_g1g5)  constant = 15.;
        if (myModel_ == TVar::H2_g5) constant = 26.167;
        if (myModel_ == TVar::H2_g4) constant = 8.9612e+07;
        if (myModel_ == TVar::H2_g2) constant = 984117;
        if (myModel_ == TVar::H2_g3) constant = 2.46285e+08;
        if (myModel_ == TVar::H2_g6) constant = 1.6201e+06;
        if (myModel_ == TVar::H2_g7) constant = 1.95307e+07;
        if (myModel_ == TVar::H2_g8) constant = 1.29087e+08;
        if (myModel_ == TVar::H2_g9) constant = 5.11404e+06;
        if (myModel_ == TVar::H2_g10) constant = 2.4183e+11;
      }
    }
    // production independent calculations
    if (myProduction_ == TVar::ZZINDEPENDENT) {
      if (flavor == 3) {
        if (myModel_ == TVar::H1minus)  constant = 1.3e+10;
        if (myModel_ == TVar::H1plus)  constant = 1.3e+10;
        if (myModel_ == TVar::H2_g1g5)  constant = 1.6e+9;
        if (myModel_ == TVar::H2_g5) constant = 2.71213e+09;
        if (myModel_ == TVar::H2_g4) constant = 1.01932e+16;
        if (myModel_ == TVar::H2_g2) constant = 9.29888e+13;
        if (myModel_ == TVar::H2_g3) constant = 2.53613e+16;
        if (myModel_ == TVar::H2_g6) constant = 1.43664e+14;
        if (myModel_ == TVar::H2_g7) constant = 2.81011e+15;
        if (myModel_ == TVar::H2_g8) constant = 1.40936e+16;
        if (myModel_ == TVar::H2_g9) constant = 5.57788e+14;
        if (myModel_ == TVar::H2_g10) constant = 2.73432e+19;
      }
      else{
        if (myModel_ == TVar::H1minus)  constant = 1.6e+10;
        if (myModel_ == TVar::H1plus)  constant = 1.4e+10;
        if (myModel_ == TVar::H2_g1g5)  constant = 2.0e+9;
        if (myModel_ == TVar::H2_g5) constant = 3.30598e+09;
        if (myModel_ == TVar::H2_g4) constant = 1.01932e+16;
        if (myModel_ == TVar::H2_g2) constant = 1.16336e+14;
        if (myModel_ == TVar::H2_g3) constant = 2.76942e+16;
        if (myModel_ == TVar::H2_g6) constant = 1.68929e+14;
        if (myModel_ == TVar::H2_g7) constant = 2.22827e+15;
        if (myModel_ == TVar::H2_g8) constant = 1.39674e+16;
        if (myModel_ == TVar::H2_g9) constant = 5.63394e+14;
        if (myModel_ == TVar::H2_g10) constant = 2.54946e+19;
      }
    }
    // vbfMELA constants
    if (myProduction_ == TVar::JJGG){
      constant = 1.8e-5;
      if (myModel_ == TVar::H0minus) constant *= 1.0017;
    }
    if (myProduction_ == TVar::JJVBF){
      if (myModel_ == TVar::H0minus) constant = 0.067;
    }
    if (myProduction_ == TVar::ttH || myProduction_ == TVar::bbH){
      if (myModel_ == TVar::H0minus) constant = pow(1.593, 2);
    }

  }
  else if (myME_ == TVar::MCFM){

    if (flavor==1){
      // 4e scale factors
      // for qqZZ
      if (myProduction_ == TVar::ZZQQB || (myProduction_ == TVar::ZZINDEPENDENT &&  myModel_ == TVar::bkgZZ)){
        if (mZZ > 900) constant = vaScale_4e->Eval(900.);
        else if (mZZ <  70) constant = vaScale_4e->Eval(70.);
        else constant = vaScale_4e->Eval(mZZ);
      }
    }
    else if (flavor==2){
      // 4mu scale factors
      // for qqZZ
      if (myProduction_ == TVar::ZZQQB || (myProduction_ == TVar::ZZINDEPENDENT &&  myModel_ == TVar::bkgZZ)){
        if (mZZ > 900) constant = vaScale_4mu->Eval(900.);
        else if (mZZ <  70) constant = vaScale_4mu->Eval(70.);
        else constant = vaScale_4mu->Eval(mZZ);
      }
    }
    else if (flavor==3){
      // 2e2mu scale factors
      // for qqZZ
      if (myProduction_ == TVar::ZZQQB || (myProduction_ == TVar::ZZINDEPENDENT &&  myModel_ == TVar::bkgZZ)){
        if (mZZ > 900) constant = vaScale_2e2mu->Eval(900.);
        else if (mZZ <  70) constant = vaScale_2e2mu->Eval(70.);
        else constant = vaScale_2e2mu->Eval(mZZ);
      }
    }
  }

  if (useOldggZZConstants && (myME_ == TVar::MCFM || myME_ == TVar::JHUGen) && myProduction_ == TVar::ZZGG && myModel_==TVar::bkgZZ_SMHiggs){
    // note: constants are being added to ggZZ ME calculation 
    // for the purpose of building DggZZ to separate qqZZ from
    // ggZZ.  These constants have been tune on the unadulterated
    // qqZZ ME calculation, so the qqZZ scale factors should be 
    // included inorder to cancel those in the qqZZ ME calc

    constant=1.;
    if (flavor==1){
      // 4e scale factors
      if (mZZ > 900) constant = vaScale_4e->Eval(900.)/DggZZ_scalefactor->Eval(900.);
      else if (mZZ <  110) constant = vaScale_4e->Eval(110.)/DggZZ_scalefactor->Eval(110.);
      else constant = vaScale_4e->Eval(mZZ)/DggZZ_scalefactor->Eval(mZZ);
    }
    else if (flavor==2){
      // 4mu scale factors
      if (mZZ > 900) constant = vaScale_4mu->Eval(900.)/DggZZ_scalefactor->Eval(900.);
      else if (mZZ <  110) constant = vaScale_4mu->Eval(110.)/DggZZ_scalefactor->Eval(110.);
      else constant = vaScale_4mu->Eval(mZZ)/DggZZ_scalefactor->Eval(mZZ);
    }
    else if (flavor==3){
      // 2e2mu scale factors
      if (mZZ > 900) constant = vaScale_2e2mu->Eval(900.)/DggZZ_scalefactor->Eval(900.);
      else if (mZZ <  110) constant = vaScale_2e2mu->Eval(110.)/DggZZ_scalefactor->Eval(110.);
      else constant = vaScale_2e2mu->Eval(mZZ)/DggZZ_scalefactor->Eval(mZZ);
    }
  }

  return constant;
}

// Full parton-by-parton ME record
MelaIO* Mela::getIORecord(){ return ZZME->get_IORecord(); }

void Mela::computeP(
  float mZZ, float mZ1, float mZ2, // input kinematics
  float costhetastar,
  float costheta1,
  float costheta2,
  float phi,
  float phi1,
  int flavor,
  double selfDHvvcoupl_input[SIZE_HVV][2],
  float& prob
  ){
  selfDHqqcoupl[0][0] = 1.0; // Don't set these, and you will get 0.
  selfDHggcoupl[0][0] = 1.0;

  for (int im=0; im<2; im++){
    for (int ic=0; ic<SIZE_HVV; ic++){
      selfDHzzcoupl[ic][im] = selfDHvvcoupl_input[ic][im];
      selfDHwwcoupl[ic][im] = selfDHvvcoupl_input[ic][im]; // Just for extra protection since differentiate_HWW_HZZ is set to false.
    }
  }
  computeP_selfDspin0(
    mZZ, mZ1, mZ2,
    costhetastar, costheta1, costheta2, phi, phi1,
    flavor,
    prob
    );
}
// You must absolutely know what you are doing if you are using the function below!
void Mela::computeP_selfDspin0(
  float mZZ, float mZ1, float mZ2, // input kinematics
  float costhetastar,
  float costheta1,
  float costheta2,
  float phi,
  float phi1,
  int flavor,
  float& prob
  ){
  reset_PAux();

  if (!(myModel_ == TVar::SelfDefine_spin0 || myME_ == TVar::MCFM)) cout << " Error: This method only applies to spin0, set Process to SelfDefine_spin0 or ME to MCFM!" << endl;
  else{
    if (myME_ == TVar::JHUGen || myME_ == TVar::MCFM){
      ZZME->set_SpinZeroCouplings(
        selfDHvvcoupl_freenorm,
        selfDHqqcoupl,
        selfDHggcoupl,
        selfDHzzcoupl,
        selfDHwwcoupl,
        selfDHzzLambda_qsq,
        selfDHwwLambda_qsq,
        selfDHzzCLambda_qsq,
        selfDHwwCLambda_qsq,
        differentiate_HWW_HZZ
        );

      //initialize variables
      checkZorder(mZ1, mZ2, costhetastar, costheta1, costheta2, phi, phi1);
      ZZME->computeXS(
        mZZ, mZ1, mZ2,
        costhetastar, costheta1, costheta2,
        phi, phi1, flavor,
        myModel_, myME_, myProduction_, prob
        );
    }
    else if (myME_ == TVar::ANALYTICAL){ // !Needs further generalization!
      costhetastar_rrv->setVal(costhetastar);
      costheta1_rrv->setVal(costheta1);
      costheta2_rrv->setVal(costheta2);
      phi_rrv->setVal(phi);
      phi1_rrv->setVal(phi1);
      z1mass_rrv->setVal(mZ1);
      z2mass_rrv->setVal(mZ2);
      mzz_rrv->setVal(mZZ);
      Y_rrv->setConstant(true); // Just to avoid integrating over this variable unnecessarily

      configureAnalyticalPDFs();

      if (myProduction_==TVar::ZZINDEPENDENT){
        RooAbsPdf* integral = (RooAbsPdf*)pdf->createIntegral(RooArgSet(*costhetastar_rrv, *phi1_rrv));
        prob = integral->getVal();
        delete integral;
      }
      else prob = pdf->getVal();

      Y_rrv->setConstant(false);
    }
    else cout << "ERROR: this method only works for JHUGen or MCFM or ANALYTICAL" << endl;
  }

  reset_SelfDCouplings();
}
void Mela::computeP_selfDspin2(
  float mZZ, float mZ1, float mZ2, // input kinematics
  float costhetastar,
  float costheta1,
  float costheta2,
  float phi,
  float phi1,
  int flavor,
  double selfDGggcoupl_input[SIZE_GGG][2],
  double selfDGvvcoupl_input[SIZE_GVV][2],
  float& prob
  ){
  selfDGqqcoupl[0][0]=1.0; // Set these incorrectly, and you see left-right asymmetries in qqG (or nothing)
  selfDGqqcoupl[1][0]=1.0;
  for (int im=0; im<2; im++){
    for (int ic=0; ic<SIZE_GGG; ic++) selfDGggcoupl[ic][im] = selfDGggcoupl_input[ic][im];
    for (int ic=0; ic<SIZE_GVV; ic++) selfDGvvcoupl[ic][im] = selfDGvvcoupl_input[ic][im];
  }

  computeP_selfDspin2(
    mZZ, mZ1, mZ2,
    costhetastar, costheta1, costheta2, phi, phi1,
    flavor,
    prob
    );
  //reset_SelfDCouplings();
}
void Mela::computeP_selfDspin2(
  float mZZ, float mZ1, float mZ2, // input kinematics
  float costhetastar,
  float costheta1,
  float costheta2,
  float phi,
  float phi1,
  int flavor,
  float& prob
  ){
  reset_PAux();

  if (myModel_ != TVar::SelfDefine_spin2) cout << " Error: This method only applies to spin2, set Process to SelfDefine_spin2!"<<endl;
  else{
    if (myME_ == TVar::JHUGen){
      ZZME->set_SpinTwoCouplings(selfDGqqcoupl, selfDGggcoupl, selfDGvvcoupl);

      //initialize variables
      checkZorder(mZ1, mZ2, costhetastar, costheta1, costheta2, phi, phi1);
      ZZME->computeXS(
        mZZ, mZ1, mZ2,
        costhetastar, costheta1, costheta2,
        phi, phi1, flavor,
        myModel_, myME_, myProduction_,
        prob
        );
    }
    else if (myME_ == TVar::ANALYTICAL){

      costhetastar_rrv->setVal(costhetastar);
      costheta1_rrv->setVal(costheta1);
      costheta2_rrv->setVal(costheta2);
      phi_rrv->setVal(phi);
      phi1_rrv->setVal(phi1);
      z1mass_rrv->setVal(mZ1);
      z2mass_rrv->setVal(mZ2);
      mzz_rrv->setVal(mZZ);

      bool checkImCoupl=false;
      for (int i =0; i<SIZE_GVV; i++){
        if (selfDGvvcoupl[i][1]!=0){
          cout << "Error: MELA does not support complex coupling for the moment! " << endl;
          checkImCoupl = true;
          break;
        }
      }
      if (!checkImCoupl){
        configureAnalyticalPDFs();

        if (myProduction_==TVar::ZZINDEPENDENT){
          RooAbsPdf* integral = (RooAbsPdf*)pdf->createIntegral(RooArgSet(*costhetastar_rrv, *phi1_rrv));
          prob = integral->getVal();
          delete integral;
        }
        else prob = pdf->getVal();
      }
    }
    else cout << "ERROR: this method only works for JHUGen and ANALYTICAL" << endl;
  }

  reset_SelfDCouplings();
}
void Mela::computeP_selfDspin1(
  float mZZ, float mZ1, float mZ2, // input kinematics
  float costhetastar,
  float costheta1,
  float costheta2,
  float phi,
  float phi1,
  int flavor,
  double selfDZvvcoupl_input[SIZE_ZVV][2],
  float& prob
  ){
  // Initialize the quark_left_right couplings
  selfDZqqcoupl[0][0]=1.0;
  selfDZqqcoupl[1][0]=1.0;
  for (int im=0; im<2; im++){
    for (int ic=0; ic<SIZE_ZVV; ic++) selfDZvvcoupl[ic][im] = selfDZvvcoupl_input[ic][im];
  }

  computeP_selfDspin1(
    mZZ, mZ1, mZ2,
    costhetastar, costheta1, costheta2, phi, phi1,
    flavor,
    prob
    );
  //reset_SelfDCouplings();
}
void Mela::computeP_selfDspin1(
  float mZZ, float mZ1, float mZ2, // input kinematics
  float costhetastar,
  float costheta1,
  float costheta2,
  float phi,
  float phi1,
  int flavor,
  float& prob
  ){
  reset_PAux();

  if (myModel_ != TVar::SelfDefine_spin1) cout << " Error: This method only applies to spin1, set Process to SelfDefine_spin1!" << endl;
  else{
    if (myME_ == TVar::JHUGen){
      ZZME->set_SpinOneCouplings(selfDZqqcoupl, selfDZvvcoupl);

      //initialize variables
      checkZorder(mZ1, mZ2, costhetastar, costheta1, costheta2, phi, phi1);
      ZZME->computeXS(
        mZZ, mZ1, mZ2,
        costhetastar, costheta1, costheta2,
        phi, phi1, flavor,
        myModel_, myME_, myProduction_,
        prob
        );
    }
    else if (myME_ == TVar::ANALYTICAL){ // !Needs generalization!
      costhetastar_rrv->setVal(costhetastar);
      costheta1_rrv->setVal(costheta1);
      costheta2_rrv->setVal(costheta2);
      phi_rrv->setVal(phi);
      phi1_rrv->setVal(phi1);
      z1mass_rrv->setVal(mZ1);
      z2mass_rrv->setVal(mZ2);
      mzz_rrv->setVal(mZZ);

      bool checkImCoupl = false;
      for (int i =0; i<SIZE_ZVV; i++){
        if (selfDZvvcoupl[i][1]!=0){
          cout << "Error: MELA does not support complex coupling for the moment! " << endl;
          checkImCoupl = true;
          break;
        }
      }
      if (!checkImCoupl){
        configureAnalyticalPDFs();
        spin1Model->g1Val->setVal(selfDZvvcoupl[0][0]);
        spin1Model->g2Val->setVal(selfDZvvcoupl[1][0]);
        if (myProduction_==TVar::ZZINDEPENDENT){
          RooAbsPdf* integral = (RooAbsPdf*)pdf->createIntegral(RooArgSet(*costhetastar_rrv, *phi1_rrv));
          prob = integral->getVal();
          delete integral;
        }
        else prob = pdf->getVal();
      }
    }
    else cout << "ERROR: this method only works for JHUGen and ANALYTICAL" << endl;
  }

  reset_SelfDCouplings();
}
void Mela::computeP(
  float mZZ, float mZ1, float mZ2, // input kinematics
  float costhetastar,
  float costheta1,
  float costheta2,
  float phi,
  float phi1,
  int flavor,
  double selfDHvvcoupl_freenorm_input[SIZE_HVV_FREENORM],
  float& prob
  ){
  reset_PAux();

  selfDHqqcoupl[0][0] = 1.0;
  selfDHggcoupl[0][0] = 1.0;
  selfDHzzcoupl[0][0] = 1.0;
  selfDHwwcoupl[0][0] = 1.0;

  if ((myME_==TVar::JHUGen || myME_==TVar::MCFM) && myModel_==TVar::bkgZZ_SMHiggs){
    ZZME->set_SpinZeroCouplings(
      selfDHvvcoupl_freenorm,
      selfDHqqcoupl,
      selfDHggcoupl,
      selfDHzzcoupl,
      selfDHwwcoupl,
      selfDHzzLambda_qsq,
      selfDHwwLambda_qsq,
      selfDHzzCLambda_qsq,
      selfDHwwCLambda_qsq,
      differentiate_HWW_HZZ
      );

    float constant=1.;
    //initialize variables
    checkZorder(mZ1, mZ2, costhetastar, costheta1, costheta2, phi, phi1);
    ZZME->computeXS(
      mZZ, mZ1, mZ2,
      costhetastar, costheta1, costheta2,
      phi, phi1, flavor,
      myModel_, myME_, myProduction_,
      prob
      );

    constant = getConstant(flavor, mZZ, true);
    prob *= constant;

    reset_SelfDCouplings();
  }
  else computeP(mZZ, mZ1, mZ2, costhetastar, costheta1, costheta2, phi, phi1, flavor, prob);
}

void Mela::computeP(
  TLorentzVector Z1_lept1, int Z1_lept1Id, // input 4-vectors
  TLorentzVector Z1_lept2, int Z1_lept2Id,
  TLorentzVector Z2_lept1, int Z2_lept1Id,
  TLorentzVector Z2_lept2, int Z2_lept2Id,
  float& prob // output probability
  ){

  // get flavor type
  // NEED TO INCLUDE SOME PROTECTION SO THAT USER CANT                  
  // PASS FOUR-VECTORS IN WRONG ORDER.  FOR NOW ASSUMING                
  // THEY ARE PASSED AS e-,e+,mu-,mu+                                   
  // ------------------ channel ------------------------                

  if (mela::forbidMassiveLeptons){
    mela::constrainedRemoveLeptonMass(Z1_lept1, Z1_lept2);
    mela::constrainedRemoveLeptonMass(Z2_lept1, Z2_lept2);
  }

  int flavor;

  if (abs(Z1_lept1Id)==abs(Z1_lept2Id)&&
    abs(Z1_lept1Id)==abs(Z2_lept1Id)&&
    abs(Z1_lept1Id)==abs(Z2_lept2Id)){
    if (abs(Z1_lept1Id)==11) flavor=1;
    else flavor=2;
  }
  else flavor=3;

  //compute angles  
  float m1=(Z1_lept1 + Z1_lept2).M();
  float m2=(Z2_lept1 + Z2_lept2).M();

  TLorentzVector ZZ = (Z1_lept1 + Z1_lept2 + Z2_lept1 + Z2_lept2);
  float mzz = ZZ.M();

  float costhetastar, costheta1, costheta2, phi, phi1;

  mela::computeAngles(
    Z1_lept1, Z1_lept1Id, Z1_lept2, Z1_lept2Id,
    Z2_lept1, Z2_lept1Id, Z2_lept2, Z2_lept2Id,
    costhetastar, costheta1, costheta2, phi, phi1
    );

  computeP(
    mzz, m1, m2,
    costhetastar,
    costheta1,
    costheta2,
    phi, phi1,
    flavor,
    prob
    );
}

void Mela::computeP( // This is actually a wrapper!
  TLorentzVector Z1_lept1, int Z1_lept1Id, // input 4-vectors
  TLorentzVector Z1_lept2, int Z1_lept2Id, 
  TLorentzVector Z2_lept1, int Z2_lept1Id,
  TLorentzVector Z2_lept2, int Z2_lept2Id,
  double selfDHvvcoupl_freenorm_input[SIZE_HVV_FREENORM],
  float& prob // output probability
  ){

  if (mela::forbidMassiveLeptons){
    mela::constrainedRemoveLeptonMass(Z1_lept1, Z1_lept2);
    mela::constrainedRemoveLeptonMass(Z2_lept1, Z2_lept2);
  }

  if (myME_==TVar::JHUGen && myModel_==TVar::bkgZZ_SMHiggs){
    int flavor;

    if (abs(Z1_lept1Id)==abs(Z1_lept2Id)&&
      abs(Z1_lept1Id)==abs(Z2_lept1Id)&&
      abs(Z1_lept1Id)==abs(Z2_lept2Id)){
      if (abs(Z1_lept1Id)==11) flavor=1;
      else flavor=2;
    }
    else flavor=3;

    //compute angles  
    TLorentzVector ZZ = (Z1_lept1 + Z1_lept2 + Z2_lept1 + Z2_lept2);
    float m1=(Z1_lept1 + Z1_lept2).M();
    float m2=(Z2_lept1 + Z2_lept2).M();
    float mzz = ZZ.M();
    float costhetastar, costheta1, costheta2, phi, phi1;
    mela::computeAngles(Z1_lept1, Z1_lept1Id, Z1_lept2, Z1_lept2Id,
      Z2_lept1, Z2_lept1Id, Z2_lept2, Z2_lept2Id,
      costhetastar, costheta1, costheta2, phi, phi1);

    computeP(
      mzz, m1, m2,
      costhetastar,
      costheta1,
      costheta2,
      phi, phi1,
      flavor, selfDHvvcoupl_freenorm,
      prob
      );
  }
}


void Mela::computeProdDecP(
  TLorentzVector jet[2],
  TLorentzVector Higgs_daughter[4],
  int jet_pdgid[2],
  int Higgs_daughter_pdgid[4],
  double selfDHvvcoupl_input[SIZE_HVV][2],
  double selfDHwwcoupl_input[SIZE_HVV][2],
  float& prob
  ){
  for (int im=0; im<2; im++){
    for (int ic=0; ic<SIZE_HVV_VBF; ic++){
      selfDHzzcoupl[ic][im] = selfDHvvcoupl_input[ic][im];
      selfDHwwcoupl[ic][im] = selfDHwwcoupl_input[ic][im]; // Just for extra protection since differentiate_HWW_HZZ is set to false.
    }
  }
  computeProdDecP(
    jet, Higgs_daughter,
    jet_pdgid, Higgs_daughter_pdgid,
    prob);
}
void Mela::computeProdDecP(
  TLorentzVector jet[2],
  TLorentzVector Higgs_daughter[4],
  int jet_pdgid[2],
  int Higgs_daughter_pdgid[4],
  float& prob
  ){
  reset_PAux();
  float constant=1.;

  bool hasFailed = true;
  if (myME_ != TVar::MCFM){
    cout << "Mela::computeProdDecP ME is not supported for ME " << myME_ << endl;
    hasFailed = true;
  }
  if (myProduction_ != TVar::JJVBF){
    cout << "Mela::computeProdDecP production mode is not supported for production " << myProduction_ << endl;
    hasFailed = true;
  }
  if (hasFailed) prob=0;
  else{
    TLorentzVector jet1massless(0, 0, 0, 0);
    TLorentzVector jet2massless(0, 0, 0, 0);
    TLorentzVector higgs(0, 0, 0, 0);

    if (abs(Higgs_daughter_pdgid[0]) <= 16 && abs(Higgs_daughter_pdgid[0]) >= 11 && abs(Higgs_daughter_pdgid[1]) <= 16 && abs(Higgs_daughter_pdgid[1]) >= 11){
      if (mela::forbidMassiveLeptons) mela::constrainedRemoveLeptonMass(Higgs_daughter[0], Higgs_daughter[1]);
    }
    if (abs(Higgs_daughter_pdgid[2]) <= 16 && abs(Higgs_daughter_pdgid[2]) >= 11 && abs(Higgs_daughter_pdgid[3]) <= 16 && abs(Higgs_daughter_pdgid[3]) >= 11){
      if (mela::forbidMassiveLeptons) mela::constrainedRemoveLeptonMass(Higgs_daughter[2], Higgs_daughter[3]);
    }
    mela::computeJetMassless(jet[0], jet1massless);
    mela::computeJetMassless(jet[1], jet2massless);
    jet[0].SetXYZT(jet1massless.X(), jet1massless.Y(), jet1massless.Z(), jet1massless.T());
    jet[1].SetXYZT(jet2massless.X(), jet2massless.Y(), jet2massless.Z(), jet2massless.T());

    for (int dd=0; dd<4; dd++) higgs = higgs + Higgs_daughter[dd];
    TLorentzVector total=jet[0]+jet[1]+higgs;
    // Preserve Z direction
    higgs.Boost(-total.BoostVector().X(), -total.BoostVector().Y(), 0);
    for (int dd=0; dd<4; dd++) Higgs_daughter[dd].Boost(-total.BoostVector().X(), -total.BoostVector().Y(), 0);
    for (int dd=0; dd<2; dd++) jet[dd].Boost(-total.BoostVector().X(), -total.BoostVector().Y(), 0);

    ZZME->computeProdXS_VVHVV(
      jet, Higgs_daughter,
      jet_pdgid, Higgs_daughter_pdgid,
      myModel_, myME_, myProduction_,
      prob
      );

    constant = getConstant(3, (float)higgs.M());
    prob *= constant;
  }

  reset_SelfDCouplings();
}


void Mela::computeProdP(
  TLorentzVector Jet1, int Jet1_Id,
  TLorentzVector Jet2, int Jet2_Id,
  TLorentzVector Decay1, int Decay1_Id,
  TLorentzVector Decay2, int Decay2_Id,
  double selfDHggcoupl_input[SIZE_HGG][2],
  double selfDHvvcoupl_input[SIZE_HVV_VBF][2],
  double selfDHwwcoupl_input[SIZE_HVV_VBF][2],
  float& prob){
  for (int im=0; im<2; im++){
    for (int ic=0; ic<SIZE_HGG; ic++) selfDHggcoupl[ic][im] = selfDHggcoupl_input[ic][im];
    for (int ic=0; ic<SIZE_HVV_VBF; ic++){
      selfDHzzcoupl_NoGamma[ic][im] = selfDHvvcoupl_input[ic][im];
      selfDHwwcoupl_NoGamma[ic][im] = selfDHwwcoupl_input[ic][im]; // Just for extra protection since differentiate_HWW_HZZ is set to false.
    }
  }
  computeProdP(
    Jet1, Jet1_Id,
    Jet2, Jet2_Id,
    Decay1, Decay1_Id,
    Decay2, Decay2_Id,
    prob
    );
}
void Mela::computeProdP(
  TLorentzVector Jet1, int Jet1_Id,
  TLorentzVector Jet2, int Jet2_Id,
  TLorentzVector Decay1, int Decay1_Id,
  TLorentzVector Decay2, int Decay2_Id,
  float& prob
  ){
  reset_PAux();
  bool isJet2Fake = false;

  float constant=1.;
  TLorentzVector nullFourVector(0, 0, 0, 0);
  TLorentzVector jet1massless(0, 0, 0, 0);
  TLorentzVector jet2massless(0, 0, 0, 0);
  TLorentzVector higgs;
  if (Decay2==nullFourVector || Decay2_Id==0){
    if (Decay1_Id==25) higgs=Decay1;
    if (Decay1_Id!=25){
      cout << "No Higgs event passed. Returning prob=-99." << endl;
      prob=-99;
    }
  }
  else higgs=Decay1+Decay2;

  if (Jet1==nullFourVector && !(Jet2==nullFourVector)){
    TLorentzVector JetTmp = Jet2;
    Jet2=Jet1;
    Jet1=JetTmp;
  }
  mela::computeJetMassless(Jet1, jet1massless);
  if (!(Jet2==nullFourVector) && myProduction_ != TVar::JH){
    mela::computeJetMassless(Jet2, jet2massless);
  }
  else if (myProduction_ == TVar::JJGG || myProduction_ == TVar::JJVBF) {
    mela::computeFakeJet(jet1massless, higgs, jet2massless);
    isJet2Fake=true;
  }
  /*
  cout << "MELA:"
    << " Higgs Pz: " <<  higgs.Pz()
    << " Higgs P: " <<  higgs.P()
    << " Higgs E: " <<  higgs.T()
    << " Jet 1 Pz: " <<  jet1massless.Pz()
    << " Jet 1 P: " <<  jet1massless.P()
    << " Jet 1 E: " <<  jet1massless.T()
    << " Jet 2 Pz: " <<  jet2massless.Pz()
    << " Jet 2 P: " <<  jet2massless.P()
    << " Jet 2 E: " <<  jet2massless.T() << endl;
  */

  //if (myProduction_ == TVar::JH) cout << "Massless jet 2 E: " << jet2massless.T() << ", p: " << jet2massless.P() << endl;
  if (isJet2Fake){
    //cout << "Jet is fake" << endl;
    //int n_MEcomp=0;

    int nGrid=11;
    std::vector<double> etaArray;
    std::vector<double> pArray;
    double eta_max = 15;
    if (jet2massless.Pt()>0) eta_max = max(eta_max, 1.2*fabs(jet2massless.Eta()));
    double eta_min = -eta_max;

    for (int iter=0; iter<nGrid; iter++){
      float prob_temp=-1;
      TLorentzVector higgs_temp;
      TLorentzVector jet1massless_temp;
      TLorentzVector jet2massless_temp;
      higgs_temp.SetXYZT(higgs.X(), higgs.Y(), higgs.Z(), higgs.T());
      jet1massless_temp.SetXYZT(jet1massless.X(), jet1massless.Y(), jet1massless.Z(), jet1massless.T());

      double jet2temp_eta = ((double)iter)*(eta_max-eta_min) / (((double)nGrid) - 1.) + eta_min;
      etaArray.push_back(jet2temp_eta);
      double jet2temp_sinh_eta = TMath::SinH(jet2temp_eta);
      double jet2temp_pz = jet2massless.Pt()*jet2temp_sinh_eta;
      jet2massless_temp.SetZ(jet2temp_pz);
      jet2massless_temp.SetX(jet2massless.X()); jet2massless_temp.SetY(jet2massless.Y()); jet2massless_temp.SetT(jet2massless_temp.P());
      /*
      cout << "MELA:"
        << " Higgs Pz: " <<  higgs_temp.Pz()
        << " Higgs P: " <<  higgs_temp.P()
        << " Higgs E: " <<  higgs_temp.T()
        << " Jet 1 Pz: " <<  jet1massless_temp.Pz()
        << " Jet 1 P: " <<  jet1massless_temp.P()
        << " Jet 1 E: " <<  jet1massless_temp.T()
        << " Jet 2 Pz: " <<  jet2massless_temp.Pz()
        << " Jet 2 P: " <<  jet2massless_temp.P()
        << " Jet 2 Pt: " <<  jet2massless_temp.Pt()
        << " Jet 2 E: " <<  jet2massless_temp.T() << endl;
      */
      TLorentzVector total_temp=jet1massless_temp+jet2massless_temp+higgs_temp;
      jet1massless_temp.Boost(-total_temp.BoostVector().x(), -total_temp.BoostVector().y(), 0);
      jet2massless_temp.Boost(-total_temp.BoostVector().x(), -total_temp.BoostVector().y(), 0);
      higgs_temp.Boost(-total_temp.BoostVector().x(), -total_temp.BoostVector().y(), 0);
      if (myModel_ == TVar::SelfDefine_spin0) ZZME->set_SpinZeroCouplings_NoGamma(
        selfDHvvcoupl_freenorm,
        selfDHqqcoupl,
        selfDHggcoupl,
        selfDHzzcoupl_NoGamma,
        selfDHwwcoupl_NoGamma,
        selfDHzzLambda_qsq,
        selfDHwwLambda_qsq,
        selfDHzzCLambda_qsq,
        selfDHwwCLambda_qsq,
        differentiate_HWW_HZZ
        );
      ZZME->computeProdXS_JJH(
        jet1massless_temp, jet2massless_temp, higgs_temp,
        myModel_, myME_, myProduction_,
        prob_temp
        );
      pArray.push_back((double)prob_temp);
    }

    double* xGrid;
    double* yGrid;
    const double grid_precision = 0.05;
    int ctr_iter=0;
    for (int iG=0; iG<nGrid-1; iG++){ // For each spacing, first compare the average of end points to spline value
      if (pArray[iG]==pArray[iG+1]) continue;
      if (etaArray[iG]==etaArray[iG+1]) continue;

      ctr_iter++;

      xGrid = new double[nGrid];
      yGrid = new double[nGrid];
      for (int iter=0; iter<nGrid; iter++){ // Fill the arrays
        xGrid[iter] = (double)etaArray[iter];
        yGrid[iter] = (double)pArray[iter];
      }

      TGraph* interpolator = new TGraph(nGrid, xGrid, yGrid);
      double derivative_first = (yGrid[1]-yGrid[0])/(xGrid[1]-xGrid[0]);
      double derivative_last = (yGrid[nGrid-1]-yGrid[nGrid-2])/(xGrid[nGrid-1]-xGrid[nGrid-2]);
      TSpline3* spline = new TSpline3("spline", interpolator, "b1e1", derivative_first, derivative_last);
      double x_middle = (xGrid[iG]+xGrid[iG+1])*0.5;
      double y_middle = (yGrid[iG]+yGrid[iG+1])*0.5;
      double y_sp = spline->Eval(x_middle);
      if (y_sp<0) y_sp = 0;

      std::vector<double>::iterator gridIt;

      if (fabs(y_sp-y_middle)<grid_precision*fabs(y_middle) || fabs(xGrid[iG+1]-xGrid[iG])<1e-3){
        gridIt = pArray.begin()+iG+1;
        pArray.insert(gridIt, y_sp);
        gridIt = etaArray.begin()+iG+1;
        etaArray.insert(gridIt, x_middle);
        iG++; // Pass to next bin
      }
      else{
        //cout << "Spline: " << y_sp << "\tmean: " << y_middle << endl;
        //cout << "Current precision: " << fabs(y_sp-y_middle)/fabs(y_middle) << endl;

        float prob_temp=-1;
        TLorentzVector higgs_temp;
        TLorentzVector jet1massless_temp;
        TLorentzVector jet2massless_temp;
        higgs_temp.SetXYZT(higgs.X(), higgs.Y(), higgs.Z(), higgs.T());
        jet1massless_temp.SetXYZT(jet1massless.X(), jet1massless.Y(), jet1massless.Z(), jet1massless.T());

        double jet2temp_eta = x_middle;
        gridIt = etaArray.begin()+iG+1;
        etaArray.insert(gridIt, x_middle);
        double jet2temp_sinh_eta = TMath::SinH(jet2temp_eta);
        double jet2temp_pz = jet2massless.Pt()*jet2temp_sinh_eta;
        jet2massless_temp.SetZ(jet2temp_pz);
        jet2massless_temp.SetX(jet2massless.X()); jet2massless_temp.SetY(jet2massless.Y()); jet2massless_temp.SetT(jet2massless_temp.P());

        TLorentzVector total_temp=jet1massless_temp+jet2massless_temp+higgs_temp;
        jet1massless_temp.Boost(-total_temp.BoostVector().x(), -total_temp.BoostVector().y(), 0);
        jet2massless_temp.Boost(-total_temp.BoostVector().x(), -total_temp.BoostVector().y(), 0);
        higgs_temp.Boost(-total_temp.BoostVector().x(), -total_temp.BoostVector().y(), 0);
        if (myModel_ == TVar::SelfDefine_spin0) ZZME->set_SpinZeroCouplings_NoGamma(
          selfDHvvcoupl_freenorm,
          selfDHqqcoupl,
          selfDHggcoupl,
          selfDHzzcoupl_NoGamma,
          selfDHwwcoupl_NoGamma,
          selfDHzzLambda_qsq,
          selfDHwwLambda_qsq,
          selfDHzzCLambda_qsq,
          selfDHwwCLambda_qsq,
          differentiate_HWW_HZZ
          );
        ZZME->computeProdXS_JJH(jet1massless_temp, jet2massless_temp, higgs_temp,
          myModel_, myME_, myProduction_,
          prob_temp
          );
        //if (ctr_iter>140 && ctr_iter<165){
        //  cout << iG << '\t' << jet2temp_eta << '\t' << prob_temp << endl;
        //}
        //n_MEcomp++;
        gridIt = pArray.begin()+iG+1;
        pArray.insert(gridIt, (double)prob_temp);
        iG--; // Do not pass to next bin, repeat until precision is achieved.
      }
      nGrid++;

      delete spline;
      delete interpolator;
      delete xGrid;
      delete yGrid;
    }

    auxiliaryProb = 0;
    int iGFirst=0, iGLast=nGrid-1;
    for (int iG=1; iG<nGrid; iG++){
      if (pArray[iG]>0 && pArray[iG-1]==0){
        iGFirst = iG-1;
        break;
      }
    }
    for (int iG=nGrid-2; iG>=0; iG--){
      if (pArray[iG]>0 && pArray[iG+1]==0){
        iGLast = iG+1;
        break;
      }
    }

    double dEtaGrid = etaArray[iGLast] - etaArray[iGFirst];
    for (int iG=iGFirst; iG<iGLast-1; iG++){
      double dEta = etaArray[iG+1] - etaArray[iG];
      double sumProb = pArray[iG]+pArray[iG+1];
      sumProb *= 0.5;
      dEta = dEta/dEtaGrid;
      //if (dEtaGrid<0) cout << "dEtaGrid: " << etaArray[nGrid-1] << " - " << etaArray[0] << " = " << dEtaGrid << endl;
      //if (dEta<0) cout << "dEta: " << etaArray[iG+1] << " - " << etaArray[iG] << " / " << dEtaGrid << " = " << dEta << endl;
      //if (sumProb<0) cout << "sumProb: " << pArray[iG] << " + " << pArray[iG+1] << " = " << sumProb << endl;
      double addProb = sumProb*dEta;
      auxiliaryProb += (float)addProb;
    }


    //if (n_MEcomp>=1000){
    //  cout << "N>=1000!:\n";
    //  for (int iG=0; iG<nGrid; iG++){
    //    cout << etaArray[iG] << '\t';
    //  }
    //  cout << endl;
    //  for (int iG=0; iG<nGrid; iG++){
    //    cout << pArray[iG] << '\t';
    //  }
    //  cout << '\n' << endl;
    //}
    //if (auxiliaryProb<=0){
    //  cout << "aux: " << auxiliaryProb << endl;
    //  for (int iG=0; iG<nGrid; iG++){
    //    cout << etaArray[iG] << '\t';
    //  }
    //  cout << endl;
    //  for (int iG=0; iG<nGrid; iG++){
    //    cout << pArray[iG] << '\t';
    //  }
    //  cout << '\n' << endl;
    //}
  }



  TLorentzVector total=jet1massless+jet2massless+higgs;
  jet1massless.Boost(-total.BoostVector().x(), -total.BoostVector().y(), 0);
  jet2massless.Boost(-total.BoostVector().x(), -total.BoostVector().y(), 0);
  higgs.Boost(-total.BoostVector().x(), -total.BoostVector().y(), 0);
  if (myProduction_ == TVar::JJGG || myProduction_ == TVar::JJVBF){
    if (myModel_ == TVar::SelfDefine_spin0) ZZME->set_SpinZeroCouplings_NoGamma(
      selfDHvvcoupl_freenorm,
      selfDHqqcoupl,
      selfDHggcoupl,
      selfDHzzcoupl_NoGamma,
      selfDHwwcoupl_NoGamma,
      selfDHzzLambda_qsq,
      selfDHwwLambda_qsq,
      selfDHzzCLambda_qsq,
      selfDHwwCLambda_qsq,
      differentiate_HWW_HZZ
      );
    ZZME->computeProdXS_JJH(jet1massless, jet2massless, higgs,
      myModel_, myME_, myProduction_,
      prob
      ); // Higgs + 2 jets: SBF or WBF
    if (fabs(prob)>0 && isJet2Fake) auxiliaryProb /= prob;
  }
  else if (myProduction_ == TVar::JH){
    // No anomalous couplings are implemented in HJ
    ZZME->computeProdXS_JH(
      jet1massless,
      higgs,
      myModel_, myME_, myProduction_,
      prob
      ); // Higgs + 1 jet; only SM is supported for now.
  }

  constant = getConstant(3, (float)higgs.M());
  prob *= constant;

  reset_SelfDCouplings();
}


void Mela::computeProdP(
  TLorentzVector V_daughter[2],
  TLorentzVector Higgs_daughter[4],
  int V_daughter_pdgid[2],
  int Higgs_daughter_pdgid[4],
  bool includeHiggsDecay,
  double selfDHvvcoupl_input[SIZE_HVV_VBF][2],
  float& prob
  ){
  // Dedicated function for VH ME
  selfDHggcoupl[0][0] = 1.0; // Don't set this, and you might get 0 in the future fr ggVH.

  for (int im=0; im<2; im++){
    for (int ic=0; ic<SIZE_HVV_VBF; ic++){
      selfDHzzcoupl_NoGamma[ic][im] = selfDHvvcoupl_input[ic][im];
      selfDHwwcoupl_NoGamma[ic][im] = selfDHvvcoupl_input[ic][im]; // Just for extra protection since differentiate_HWW_HZZ is set to false.
    }
  }
  computeProdP(
    V_daughter,
    Higgs_daughter,
    V_daughter_pdgid,
    Higgs_daughter_pdgid,
    includeHiggsDecay,
    prob
    );
}
void Mela::computeProdP(
  TLorentzVector V_daughter[2],
  TLorentzVector Higgs_daughter[4],
  int V_daughter_pdgid[2],
  int Higgs_daughter_pdgid[4],
  bool includeHiggsDecay,
  float& prob
  ){
  // Dedicated function for VH ME
  reset_PAux();

  if (myProduction_ == TVar::ZH || myProduction_ == TVar::WH){
    float constant=1;
    if (abs(V_daughter_pdgid[0]) <= 6){
      TLorentzVector jettmp(0, 0, 0, 0);
      mela::computeJetMassless(V_daughter[0], jettmp);
      V_daughter[0].SetXYZT(jettmp.X(), jettmp.Y(), jettmp.Z(), jettmp.T());
    }
    if (abs(V_daughter_pdgid[1]) <= 6){
      TLorentzVector jettmp(0, 0, 0, 0);
      mela::computeJetMassless(V_daughter[1], jettmp);
      V_daughter[1].SetXYZT(jettmp.X(), jettmp.Y(), jettmp.Z(), jettmp.T());
    }

    if (abs(V_daughter_pdgid[0]) <= 16 && abs(V_daughter_pdgid[0]) >= 11 && abs(V_daughter_pdgid[1]) <= 16 && abs(V_daughter_pdgid[1]) >= 11){
      if (mela::forbidMassiveLeptons) mela::constrainedRemoveLeptonMass(V_daughter[0], V_daughter[1]);
    }
    if (abs(Higgs_daughter_pdgid[0]) <= 16 && abs(Higgs_daughter_pdgid[0]) >= 11 && abs(Higgs_daughter_pdgid[1]) <= 16 && abs(Higgs_daughter_pdgid[1]) >= 11){
      if (mela::forbidMassiveLeptons) mela::constrainedRemoveLeptonMass(Higgs_daughter[0], Higgs_daughter[1]);
    }
    if (abs(Higgs_daughter_pdgid[2]) <= 16 && abs(Higgs_daughter_pdgid[2]) >= 11 && abs(Higgs_daughter_pdgid[3]) <= 16 && abs(Higgs_daughter_pdgid[3]) >= 11){
      if (mela::forbidMassiveLeptons) mela::constrainedRemoveLeptonMass(Higgs_daughter[2], Higgs_daughter[3]);
    }

    if (myModel_ == TVar::SelfDefine_spin0) ZZME->set_SpinZeroCouplings_NoGamma(
      selfDHvvcoupl_freenorm,
      selfDHqqcoupl,
      selfDHggcoupl,
      selfDHzzcoupl_NoGamma,
      selfDHwwcoupl_NoGamma,
      selfDHzzLambda_qsq,
      selfDHwwLambda_qsq,
      selfDHzzCLambda_qsq,
      selfDHwwCLambda_qsq,
      differentiate_HWW_HZZ
      );
    ZZME->computeProdXS_VH(
      V_daughter,
      Higgs_daughter,
      V_daughter_pdgid,
      Higgs_daughter_pdgid,
      includeHiggsDecay,
      myModel_,
      myME_,
      myProduction_,
      prob
      ); // VH

    float mzz = (Higgs_daughter[0]+Higgs_daughter[1]+Higgs_daughter[2]+Higgs_daughter[3]).M();
    constant = getConstant(3, mzz);
    prob *= constant;
  }

  reset_SelfDCouplings();
}


void Mela::computeProdP(
  TLorentzVector vTTH[6],
  TLorentzVector Higgs,
  int ttbar_daughters_pdgid[6],
  double selfDHqqcoupl_input[SIZE_HQQ][2],
  float& prob,
  int topDecay,
  int topProcess
  ){
  for (int im=0; im<2; im++){
    for (int ic=0; ic<SIZE_HQQ; ic++) selfDHqqcoupl[ic][im]=selfDHqqcoupl_input[ic][im];
  }
  computeProdP(vTTH, Higgs, ttbar_daughters_pdgid, prob, topDecay, topProcess);
}
void Mela::computeProdP(
  TLorentzVector vTTH[6],
  TLorentzVector Higgs,
  int ttbar_daughters_pdgid[6],
  float& prob,
  int topDecay,
  int topProcess
  ){
  reset_PAux();
  float constant = 1;

  selfDHggcoupl[0][0] = 1.0;

  ZZME->set_SpinZeroCouplings_NoGamma(
    selfDHvvcoupl_freenorm,
    selfDHqqcoupl,
    selfDHggcoupl,
    selfDHzzcoupl_NoGamma,
    selfDHwwcoupl_NoGamma,
    selfDHzzLambda_qsq,
    selfDHwwLambda_qsq,
    selfDHzzCLambda_qsq,
    selfDHwwCLambda_qsq,
    differentiate_HWW_HZZ
    );
  ZZME->computeProdXS_ttH(vTTH, Higgs, ttbar_daughters_pdgid, topDecay, topProcess, myModel_, myME_, myProduction_, prob);

  constant = getConstant(3, (float)Higgs.M());
  prob *= constant;

  reset_SelfDCouplings();
}

void Mela::computeProdP(
  TLorentzVector p_first, int id_first,
  TLorentzVector p_second, int id_second,
  TLorentzVector Higgs,
  float& prob
  ){
  TLorentzVector nullFourVector(0, 0, 0, 0);

  if (myProduction_ == TVar::ttH || myProduction_ == TVar::bbH){
    int ttbar_daughters_pdgid[6] ={ 0 };
    TLorentzVector vTTH[6];
    ttbar_daughters_pdgid[0] = id_first;
    ttbar_daughters_pdgid[3] = id_second;
    vTTH[0] = p_first;
    vTTH[3] = p_second;
    for (int vv=1; vv<3; vv++) vTTH[vv].SetXYZT(0, 0, 0, 0);
    for (int vv=4; vv<6; vv++) vTTH[vv].SetXYZT(0, 0, 0, 0);
    computeProdP(vTTH, Higgs, ttbar_daughters_pdgid, prob, 0, 2);
  }
  else if (myProduction_ == TVar::ZH || myProduction_ == TVar::WH){
    bool includeHiggsDecay = false;
    int Higgs_daughter_pdgid[4] ={ 25, 0, 0, 0 }; // Pass dummy
    TLorentzVector Higgs_daughter[4] ={ Higgs, nullFourVector, nullFourVector, nullFourVector }; // These are summed as long as includeHiggsDecay is false.
    int V_daughter_pdgid[2] ={ id_first, id_second };
    TLorentzVector V_daughter[2] ={ p_first, p_second };
    computeProdP(V_daughter, Higgs_daughter, V_daughter_pdgid, Higgs_daughter_pdgid, includeHiggsDecay, prob);
  }
  else if (myProduction_ == TVar::JH || myProduction_ == TVar::JJGG || myProduction_ == TVar::JJVBF) computeProdP(p_first, id_first, p_second, id_second, Higgs, 25, nullFourVector, 0, prob);
}



void Mela::computePM4l(
  TLorentzVector Z1_lept1, int Z1_lept1Id,  // input 4-vectors
  TLorentzVector Z1_lept2, int Z1_lept2Id,  // 
  TLorentzVector Z2_lept1, int Z2_lept1Id,
  TLorentzVector Z2_lept2, int Z2_lept2Id,
  TVar::SuperMelaSyst syst,
  float& prob){
  reset_PAux();

  // Notice: No need to correct lepton masses here since m4l is the only information entered. m4l is not supposed to change after mass corrections.
  TLorentzVector ZZ = (Z1_lept1 + Z1_lept2 + Z2_lept1 + Z2_lept2);
  float mzz = ZZ.M();
  TVar::LeptonFlavor flavor = TVar::Flavor_Dummy;

  if (abs(Z1_lept1Id)==11 &&  abs(Z1_lept2Id)==11 &&
    abs(Z2_lept1Id)==11 &&  abs(Z2_lept2Id)==11)
    flavor = TVar::Flavor_4e;

  if (abs(Z1_lept1Id)==13 &&  abs(Z1_lept2Id)==13 &&
    abs(Z2_lept1Id)==13 &&  abs(Z2_lept2Id)==13)
    flavor = TVar::Flavor_4mu;

  if (abs(Z1_lept1Id)==11 &&  abs(Z1_lept2Id)==11 &&
    abs(Z2_lept1Id)==13 &&  abs(Z2_lept2Id)==13)
    flavor = TVar::Flavor_2e2mu;

  if (abs(Z1_lept1Id)==13 &&  abs(Z1_lept2Id)==13 &&
    abs(Z2_lept1Id)==11 &&  abs(Z2_lept2Id)==11)
    flavor = TVar::Flavor_2e2mu;


  computePM4l(mzz, flavor, syst, prob);
}

void Mela::computePM4l(
  float mZZ, TVar::LeptonFlavor flavor, TVar::SuperMelaSyst syst, float& prob
  ){
  reset_PAux();
  prob=-99;//default dummy.

  if (flavor == TVar::Flavor_Dummy) // only compute things if flavor determination succeded
    return;

  switch (flavor){
  case 1: super->SetDecayChannel("4e"); break;
  case 2: super->SetDecayChannel("4mu"); break;
  case 3: super->SetDecayChannel("2e2mu"); break;
  default: std::cout << " unknown flavor: " << flavor << std::endl; exit(0);
  }


  if (syst == TVar::SMSyst_None){
    std::pair<double, double> m4lP = super->M4lProb(mZZ);
    if (myModel_ == TVar::HSMHiggs) // currently only supported signal is H(0+)
      prob = m4lP.first;
    if (myModel_ == TVar::bkgZZ) // currently only supported background is summed paramterization
      prob = m4lP.second;
  }
  else{
    //systematics for p(m4l)
    float mZZtmp=mZZ;
    float meanErr=float(super->GetSigShapeSystematic("meanCB"));
    if (syst == TVar::SMSyst_ScaleUp){
      mZZtmp = mZZ*(1.0+meanErr);
      if (mZZtmp>180.0 || mZZtmp<100)mZZtmp=mZZ;
      std::pair<double, double> m4lPScaleUp = super->M4lProb(mZZtmp);
      if (myModel_ == TVar::HSMHiggs)
        prob = m4lPScaleUp.first;
      if (myModel_ == TVar::bkgZZ)
        prob = m4lPScaleUp.second;
    }

    if (syst == TVar::SMSyst_ScaleDown){
      mZZtmp = mZZ*(1.0-meanErr);
      if (mZZtmp>180 || mZZtmp<100) mZZtmp=mZZ;
      std::pair<double, double> m4lPScaleDown = super->M4lProb(mZZtmp);
      if (myModel_ == TVar::HSMHiggs) prob = m4lPScaleDown.first;
      if (myModel_ == TVar::bkgZZ) prob = m4lPScaleDown.second;
    }

    float sigmaErr=float(super->GetSigShapeSystematic("sigmaCB"));
    float sigmaCB=float(super->GetSigShapeParameter("sigmaCB"));
    if (syst == TVar::SMSyst_ResUp || syst ==  TVar::SMSyst_ResDown){
      mZZtmp= myR->Gaus(mZZ, sigmaErr*sigmaCB);
      if (mZZtmp>180 || mZZtmp<100) mZZtmp=mZZ;
      std::pair<double, double> m4lPResUp = super->M4lProb(mZZtmp);
      if (myModel_ == TVar::HSMHiggs) prob = m4lPResUp.first;
      if (myModel_ == TVar::bkgZZ) prob = m4lPResUp.second;
    }
  }
}

void Mela::computeWeight(
  float mZZ, float mZ1, float mZ2,
  float costhetastar,
  float costheta1,
  float costheta2,
  float phi,
  float phi1,
  float& w
  ){ // Lepton interference using JHUGen
  reset_PAux();

  float dXsec_HZZ_JHU, dXsec_HZZ_JHU_interf; // temporary prob

  // calculate dXsec for 4e/4mu
  setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::ZZGG);
  computeP(mZZ, mZ1, mZ2,
    costhetastar,
    costheta1,
    costheta2,
    phi, phi1,
    1, dXsec_HZZ_JHU_interf);

  // calculate dXsec for 2e2mu
  setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::ZZGG);
  computeP(mZZ, mZ1, mZ2,
    costhetastar,
    costheta1,
    costheta2,
    phi, phi1,
    3, dXsec_HZZ_JHU);

  w = dXsec_HZZ_JHU_interf / dXsec_HZZ_JHU;

  // protect against anomalously large weights
  if (w>5.) w=25./w;
}

void Mela::computeWeight(
  float mZZ, float mZ1, float mZ2,
  float costhetastar,
  float costheta1,
  float costheta2,
  float phi,
  float phi1,
  double selfDHvvcoupl_freenorm_input[SIZE_HVV_FREENORM],
  float& w
  ){
  reset_PAux();

  float dXsec_HZZ_JHU, dXsec_HZZ_JHU_interf; // temporary prob

  // calculate dXsec for 4e/4mu
  setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::ZZGG);
  computeP(mZZ, mZ1, mZ2,
    costhetastar,
    costheta1,
    costheta2,
    phi, phi1,
    1, selfDHvvcoupl_freenorm_input, dXsec_HZZ_JHU_interf);

  // calculate dXsec for 2e2mu
  setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::ZZGG);
  computeP(mZZ, mZ1, mZ2,
    costhetastar,
    costheta1,
    costheta2,
    phi, phi1,
    3, selfDHvvcoupl_freenorm_input, dXsec_HZZ_JHU);

  w = dXsec_HZZ_JHU_interf / dXsec_HZZ_JHU;

  // protect against anomalously large weights
  if (w>5.) w=25./w;
}
void Mela::setCTotalBkgGraphs(TFile* fcontainer, TGraph* tgC[]){ // Hope it has only 3 members in the array
  string tgname = "C_TotalBkgM4l_";

  float rValues[6]={ 1, 5, 10, 15, 20, 25 }; // Possible r Values

  for (int flavor=0; flavor<3; flavor++){
    float myWidth = 1;

    char crValue[20];

    int rCode = 2; // r=10
    myWidth = rValues[rCode];
    sprintf(crValue, "D_Gamma_gg_r%.0f", myWidth);

    string ctgM4L = tgname;
    string strChannel;
    if (flavor==0) strChannel = "4e"; // Check this
    else if (flavor==1) strChannel = "4mu"; // and this
    else strChannel = "2mu2e";
    ctgM4L = ctgM4L + strChannel + "_";
    ctgM4L = ctgM4L + crValue;
    tgC[flavor] = (TGraph*)fcontainer->Get(ctgM4L.c_str());
  }
}
void Mela::constructDggr(float mzz, int flavor, float bkg_VAMCFM_noscale, float ggzz_VAMCFM_noscale, float ggHZZ_prob_pure_noscale, float ggHZZ_prob_int_noscale, float& myDggr){
  float ctotal_bkg = tgtotalbkg[flavor-1]->Eval(mzz);

  float rValues[6]={ 1, 5, 10, 15, 20, 25 };
  float total_sig_ME;
  float total_bkg_ME;
  float myWidth = 1;
  int rCode = 2;
  myWidth = rValues[rCode];

  total_sig_ME = (myWidth * ggHZZ_prob_pure_noscale + sqrt(myWidth) * ggHZZ_prob_int_noscale + ggzz_VAMCFM_noscale);
  total_bkg_ME = bkg_VAMCFM_noscale*ctotal_bkg;
  float kd_denominator = (total_sig_ME+total_bkg_ME);
  float kd = total_sig_ME/kd_denominator;
  myDggr = kd;
}
void Mela::computeD_gg(float mZZ, float mZ1, float mZ2, // input kinematics
  float costhetastar,
  float costheta1,
  float costheta2,
  float phi,
  float phi1,
  int flavor,
  TVar::MatrixElement myME,
  TVar::Process myType,
  float& prob){

  if (myME != TVar::MCFM || myType != TVar::D_gg10){
    cout << "Only support MCFM and D_gg10"<<endl;
    return;
  }

  setMelaLeptonInterference(TVar::DefaultLeptonInterf); // Override lepton interference setting

  float bkg_VAMCFM_noscale, ggzz_VAMCFM_noscale, ggHZZ_prob_pure_noscale, ggHZZ_prob_int_noscale, bkgHZZ_prob_noscale;
  setProcess(TVar::bkgZZ, myME, TVar::ZZGG);
  computeP(mZZ, mZ1, mZ2,
    costhetastar, costheta1, costheta2, phi, phi1, flavor, ggzz_VAMCFM_noscale);
  setProcess(TVar::HSMHiggs, myME, TVar::ZZGG);
  computeP(mZZ, mZ1, mZ2,
    costhetastar, costheta1, costheta2, phi, phi1, flavor, ggHZZ_prob_pure_noscale);
  setProcess(TVar::bkgZZ_SMHiggs, myME, TVar::ZZGG);
  computeP(mZZ, mZ1, mZ2,
    costhetastar, costheta1, costheta2, phi, phi1, flavor, bkgHZZ_prob_noscale);
  setProcess(TVar::bkgZZ, myME, TVar::ZZQQB);
  computeP(mZZ, mZ1, mZ2,
    costhetastar, costheta1, costheta2, phi, phi1, flavor, bkg_VAMCFM_noscale, 0);
  ggHZZ_prob_int_noscale = bkgHZZ_prob_noscale - ggHZZ_prob_pure_noscale -  ggzz_VAMCFM_noscale;
  float myDggr;
  constructDggr(mZZ, flavor, bkg_VAMCFM_noscale, ggzz_VAMCFM_noscale, ggHZZ_prob_pure_noscale, ggHZZ_prob_int_noscale, myDggr);
  prob=myDggr;
}

// Notice that these only set the members of MELA, not TEvtProb. TEvtProb resets itself.
void Mela::reset_SelfDCouplings(){
  //****Spin-0****//
  differentiate_HWW_HZZ=false;
  for (int ic=0; ic<SIZE_HVV_FREENORM; ic++) selfDHvvcoupl_freenorm[ic]=0;
  for (int im=0; im<2; im++){
    for (int ic=0; ic<SIZE_HQQ; ic++) selfDHqqcoupl[ic][im]=0;
    for (int ic=0; ic<SIZE_HGG; ic++) selfDHggcoupl[ic][im]=0;
    for (int ic=0; ic<SIZE_HVV; ic++){
      if (ic<SIZE_HVV_VBF){
        selfDHzzcoupl_NoGamma[ic][im] = 0;
        selfDHwwcoupl_NoGamma[ic][im] = 0;
      }
      selfDHzzcoupl[ic][im] = 0;
      selfDHwwcoupl[ic][im] = 0;
    }
  }
  for (int ik=0; ik<3; ik++){
    selfDHzzCLambda_qsq[ik]=0;
    selfDHwwCLambda_qsq[ik]=0;
    for (int ic=0; ic<4; ic++){ // These default values do not matter as long as the c's are 0.
      selfDHzzLambda_qsq[ic][ik] = 100.;
      selfDHwwLambda_qsq[ic][ik] = 100.;
    }
  }
  //****Spin-1****//
  for (int im=0; im<2; im++){
    for (int ic=0; ic<SIZE_ZVV; ic++) selfDZvvcoupl[ic][im] = 0;
    for (int ic=0; ic<SIZE_ZQQ; ic++) selfDZqqcoupl[ic][im] = 0;
  }
  //****Spin-2****//
  for (int im=0; im<2; im++){
    for (int ic=0; ic<SIZE_GVV; ic++) selfDGvvcoupl[ic][im] = 0;
    for (int ic=0; ic<SIZE_GGG; ic++) selfDGggcoupl[ic][im] = 0;
    for (int ic=0; ic<SIZE_GQQ; ic++) selfDGqqcoupl[ic][im] = 0;
  }
  // Did I tell you that is a lot of them?
}

void Mela::configureAnalyticalPDFs(){
  // 
  // configure the analytical calculations 
  // 
  if (myModel_==TVar::bkgZZ)  pdf = qqZZmodel;
  else if (myProduction_ == TVar::JJGG || myProduction_ == TVar::JJVBF);
  else if (myProduction_ == TVar::ZH || myProduction_ == TVar::WH);
  else if (
    myModel_ == TVar::HSMHiggs
    || myModel_ == TVar::H0minus || myModel_ == TVar::D_g1g4 || myModel_ == TVar::D_g1g4_pi_2
    || myModel_ == TVar::H0hplus || myModel_ == TVar::D_g1g2 || myModel_ == TVar::D_g1g2_pi_2
    || myModel_ == TVar::H0_g1prime2 || myModel_ == TVar::D_g1g1prime2
    || myModel_ == TVar::SelfDefine_spin0
    ){
    pdf = (RooAbsPdf*)ggSpin0Model->getPDF();
    ggSpin0Model->makeParamsConst(false);
    ggSpin0Model->resetHypotheses();

    // Add the hypotheses with best-guess coefficients
    // ZZ/WW
    if (
      myModel_ == TVar::HSMHiggs
      || myModel_ == TVar::D_g1g1prime2
      || myModel_ == TVar::D_g1g2 || myModel_ == TVar::D_g1g2_pi_2
      || myModel_ == TVar::D_g1g4 || myModel_ == TVar::D_g1g4_pi_2
      || myModel_ == TVar::D_zzzg_g1prime2 || myModel_ == TVar::D_zzzg_g1prime2_pi_2 || myModel_ == TVar::D_zzzg || myModel_ == TVar::D_zzzg_PS
      || myModel_ == TVar::D_zzgg || myModel_ == TVar::D_zzgg_PS
      ) ggSpin0Model->addHypothesis(0, 0);
    if (myModel_ == TVar::H0_g1prime2 || myModel_ == TVar::D_g1g1prime2) ggSpin0Model->addHypothesis(0, 2);
    if (myModel_ == TVar::H0hplus || myModel_ == TVar::D_g1g2 || myModel_ == TVar::D_g1g2_pi_2) ggSpin0Model->addHypothesis(1, 0, (myModel_ == TVar::D_g1g2_pi_2 ? TMath::Pi() : 0.));
    if (myModel_ == TVar::H0minus || myModel_ == TVar::D_g1g4 || myModel_ == TVar::D_g1g4_pi_2) ggSpin0Model->addHypothesis(3, 0, (myModel_ == TVar::D_g1g4_pi_2 ? TMath::Pi() : 0.));
    // ZG/ZGs
    if (myModel_ == TVar::H0_Zgsg1prime2 || myModel_ == TVar::D_zzzg_g1prime2 || myModel_ == TVar::D_zzzg_g1prime2_pi_2) ggSpin0Model->addHypothesis(4, 2, (myModel_ == TVar::D_zzzg_g1prime2_pi_2 ? TMath::Pi() : 0.));
    if (myModel_ == TVar::H0_Zgs || myModel_ == TVar::D_zzzg) ggSpin0Model->addHypothesis(5, 0);
    if (myModel_ == TVar::H0_Zgs_PS || myModel_ == TVar::D_zzzg_PS) ggSpin0Model->addHypothesis(7, 0);
    // GG/GGs/GsGs
    if (myModel_ == TVar::H0_gsgs || myModel_ == TVar::D_zzgg) ggSpin0Model->addHypothesis(8, 0);
    if (myModel_ == TVar::H0_gsgs_PS || myModel_ == TVar::D_zzgg_PS) ggSpin0Model->addHypothesis(10, 0);
    // Self-defined spin-0
    if (myModel_ == TVar::SelfDefine_spin0){
      for (int im=0; im<2; im++){
        ((RooRealVar*)ggSpin0Model->parameters.g1List[0][im])->setVal(selfDHzzcoupl[0][im]);
        ((RooRealVar*)ggSpin0Model->parameters.g2List[0][im])->setVal(selfDHzzcoupl[1][im]);
        ((RooRealVar*)ggSpin0Model->parameters.g3List[0][im])->setVal(selfDHzzcoupl[2][im]);
        ((RooRealVar*)ggSpin0Model->parameters.g4List[0][im])->setVal(selfDHzzcoupl[3][im]);

        ((RooRealVar*)ggSpin0Model->parameters.gzgs2List[0][im])->setVal(selfDHzzcoupl[4][im]);
        ((RooRealVar*)ggSpin0Model->parameters.gzgs3List[0][im])->setVal(selfDHzzcoupl[5][im]);
        ((RooRealVar*)ggSpin0Model->parameters.gzgs4List[0][im])->setVal(selfDHzzcoupl[6][im]);
        ((RooRealVar*)ggSpin0Model->parameters.ggsgs2List[0][im])->setVal(selfDHzzcoupl[7][im]);
        ((RooRealVar*)ggSpin0Model->parameters.ggsgs3List[0][im])->setVal(selfDHzzcoupl[8][im]);
        ((RooRealVar*)ggSpin0Model->parameters.ggsgs4List[0][im])->setVal(selfDHzzcoupl[9][im]);

        ((RooRealVar*)ggSpin0Model->parameters.g1List[1][im])->setVal(selfDHzzcoupl[10][im]);
        ((RooRealVar*)ggSpin0Model->parameters.g1List[2][im])->setVal(selfDHzzcoupl[11][im]);
        ((RooRealVar*)ggSpin0Model->parameters.g1List[3][im])->setVal(selfDHzzcoupl[12][im]);
        ((RooRealVar*)ggSpin0Model->parameters.g1List[4][im])->setVal(selfDHzzcoupl[13][im]);
        ((RooRealVar*)ggSpin0Model->parameters.g1List[5][im])->setVal(selfDHzzcoupl[14][im]);

        ((RooRealVar*)ggSpin0Model->parameters.g2List[1][im])->setVal(selfDHzzcoupl[15][im]);
        ((RooRealVar*)ggSpin0Model->parameters.g2List[2][im])->setVal(selfDHzzcoupl[16][im]);
        ((RooRealVar*)ggSpin0Model->parameters.g2List[3][im])->setVal(selfDHzzcoupl[17][im]);
        ((RooRealVar*)ggSpin0Model->parameters.g2List[4][im])->setVal(selfDHzzcoupl[18][im]);
        ((RooRealVar*)ggSpin0Model->parameters.g2List[5][im])->setVal(selfDHzzcoupl[19][im]);

        ((RooRealVar*)ggSpin0Model->parameters.g3List[1][im])->setVal(selfDHzzcoupl[20][im]);
        ((RooRealVar*)ggSpin0Model->parameters.g3List[2][im])->setVal(selfDHzzcoupl[21][im]);
        ((RooRealVar*)ggSpin0Model->parameters.g3List[3][im])->setVal(selfDHzzcoupl[22][im]);
        ((RooRealVar*)ggSpin0Model->parameters.g3List[4][im])->setVal(selfDHzzcoupl[23][im]);
        ((RooRealVar*)ggSpin0Model->parameters.g3List[5][im])->setVal(selfDHzzcoupl[24][im]);

        ((RooRealVar*)ggSpin0Model->parameters.g4List[1][im])->setVal(selfDHzzcoupl[25][im]);
        ((RooRealVar*)ggSpin0Model->parameters.g4List[2][im])->setVal(selfDHzzcoupl[26][im]);
        ((RooRealVar*)ggSpin0Model->parameters.g4List[3][im])->setVal(selfDHzzcoupl[27][im]);
        ((RooRealVar*)ggSpin0Model->parameters.g4List[4][im])->setVal(selfDHzzcoupl[28][im]);
        ((RooRealVar*)ggSpin0Model->parameters.g4List[5][im])->setVal(selfDHzzcoupl[29][im]);

        ((RooRealVar*)ggSpin0Model->parameters.gzgs1List[0][im])->setVal(selfDHzzcoupl[30][im]); // Zgs1_prime2

        ((RooRealVar*)ggSpin0Model->parameters.g1List[6][im])->setVal(selfDHzzcoupl[31][im]);
        ((RooRealVar*)ggSpin0Model->parameters.g1List[7][im])->setVal(selfDHzzcoupl[32][im]);
        ((RooRealVar*)ggSpin0Model->parameters.g2List[6][im])->setVal(selfDHzzcoupl[33][im]);
        ((RooRealVar*)ggSpin0Model->parameters.g2List[7][im])->setVal(selfDHzzcoupl[34][im]);
        ((RooRealVar*)ggSpin0Model->parameters.g3List[6][im])->setVal(selfDHzzcoupl[35][im]);
        ((RooRealVar*)ggSpin0Model->parameters.g3List[7][im])->setVal(selfDHzzcoupl[36][im]);
        ((RooRealVar*)ggSpin0Model->parameters.g4List[6][im])->setVal(selfDHzzcoupl[37][im]);
        ((RooRealVar*)ggSpin0Model->parameters.g4List[7][im])->setVal(selfDHzzcoupl[38][im]);
      }
      for (int qoqtqz=0; qoqtqz<3; qoqtqz++){ // 0==q1, 1==q2, 2==q12
        ((RooRealVar*)ggSpin0Model->parameters.Lambda_z1qsq[qoqtqz])->setVal(selfDHzzLambda_qsq[0][qoqtqz]);
        ((RooRealVar*)ggSpin0Model->parameters.Lambda_z2qsq[qoqtqz])->setVal(selfDHzzLambda_qsq[1][qoqtqz]);
        ((RooRealVar*)ggSpin0Model->parameters.Lambda_z3qsq[qoqtqz])->setVal(selfDHzzLambda_qsq[2][qoqtqz]);
        ((RooRealVar*)ggSpin0Model->parameters.Lambda_z4qsq[qoqtqz])->setVal(selfDHzzLambda_qsq[3][qoqtqz]);
        ((RooRealVar*)ggSpin0Model->parameters.cLambda_qsq[qoqtqz])->setVal(selfDHzzCLambda_qsq[qoqtqz]);
      }
    }
    ggSpin0Model->makeParamsConst(true);
  }
  else if (!spin1Model->configure(myModel_)) pdf = spin1Model->PDF;
  else if (
    myModel_ == TVar::H2_g1
    || myModel_ == TVar::H2_g1g5
    || myModel_ == TVar::H2_g2
    || myModel_ == TVar::H2_g3
    || myModel_ == TVar::H2_g4
    || myModel_ == TVar::H2_g5
    || myModel_ == TVar::H2_g6
    || myModel_ == TVar::H2_g7
    || myModel_ == TVar::H2_g8
    || myModel_ == TVar::H2_g9
    || myModel_ == TVar::H2_g10
    || myModel_ == TVar::SelfDefine_spin2
    ){
    pdf = (RooAbsPdf*)spin2Model->getPDF();
    spin2Model->makeParamsConst(false);
    spin2Model->resetHypotheses();
    // Add the hypotheses with best-guess coefficients
    // ZZ/WW
    if (
      myModel_ == TVar::H2_g1
      || myModel_ == TVar::H2_g1g5
      ) spin2Model->addHypothesis(0, 1.);
    if (
      myModel_ == TVar::H2_g1g5
      || myModel_ == TVar::H2_g5
      ) spin2Model->addHypothesis(4, 1.);
    if (myModel_ == TVar::H2_g2) spin2Model->addHypothesis(1, 1.);
    if (myModel_ == TVar::H2_g3) spin2Model->addHypothesis(2, 1.);
    if (myModel_ == TVar::H2_g4) spin2Model->addHypothesis(3, 1.);
    if (myModel_ == TVar::H2_g5) spin2Model->addHypothesis(4, 1.);
    if (myModel_ == TVar::H2_g6) spin2Model->addHypothesis(5, 1.);
    if (myModel_ == TVar::H2_g7) spin2Model->addHypothesis(6, 1.);
    if (myModel_ == TVar::H2_g8) spin2Model->addHypothesis(7, 1.);
    if (myModel_ == TVar::H2_g9) spin2Model->addHypothesis(8, 1.);
    if (myModel_ == TVar::H2_g10) spin2Model->addHypothesis(10, 1.);
    // Self-defined spin-2
    if (myModel_ == TVar::SelfDefine_spin2){
      for (int ig=0; ig<SIZE_GVV; ig++){
        for (int im=0; im<2; im++) ((RooRealVar*)spin2Model->parameters.bList[ig][im])->setVal(selfDGvvcoupl[ig][im]);
      }
    }
    if (myProduction_ == TVar::ZZQQB){
      spin2Model->setTensorPolarization(1, 1.);
      spin2Model->setTensorPolarization(2, 0.);
    }
    else{
      if (myModel_ == TVar::SelfDefine_spin2){
        double c1 = 2*selfDGggcoupl[0][0] + 2.*selfDGggcoupl[1][0];
        double c2 = -0.5*selfDGggcoupl[0][0] + selfDGggcoupl[2][0] + 2.*selfDGggcoupl[3][0];
        double c5 = 4*selfDGggcoupl[7][0];
        Double_t fppReal = 1./sqrt(6.) * (c1/4.*2. + 2.*c2);
        Double_t fppImag = 1./sqrt(6.) * c5;
        Double_t fmmReal = 1./sqrt(6.) * (c1/4.*2. + 2.*c2);
        Double_t fmmImag = 1./sqrt(6.)* c5;
        Double_t fmpReal = 1./4.*c1*2.;
        Double_t fmpImag = 0;
        Double_t fpp = fppImag*fppImag + fppReal*fppReal;
        Double_t fmm = fmmImag*fmmImag + fmmReal*fmmReal;
        Double_t fmp = fmpImag*fmpImag + fmpReal*fmpReal;
        spin2Model->setTensorPolarization(1, 0.); // This is wrong in the strict sense of what "SelfDefine_spin2" is.
        spin2Model->setTensorPolarization(2, 2.*fmp/(fmm+fpp+2.*fmp));
      }
      else{
        spin2Model->setTensorPolarization(1, 0.);
        spin2Model->setTensorPolarization(2, 1.);
      }
    }
    spin2Model->makeParamsConst(true);
  }
  else if (myME_ == TVar::ANALYTICAL) cout << "Mela::configureAnalyticalPDFs -> ERROR TVar::Process not applicable!!! ME: " << myME_ << ", model: " << myModel_ << endl;
}


