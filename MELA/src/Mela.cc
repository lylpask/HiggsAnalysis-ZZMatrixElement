#include <ZZMatrixElement/MELA/interface/Mela.h>
#include <ZZMatrixElement/MELA/interface/newZZMatrixElement.h>
#include <DataFormats/GeometryVector/interface/Pi.h>
#include <FWCore/ParameterSet/interface/FileInPath.h>

#include "computeAngles.h"
#include "AngularPdfFactory.h"
#include "VectorPdfFactory.h"
#include "TensorPdfFactory.h"
#include "RooqqZZ_JHU_ZgammaZZ_fast.h"
#include "RooqqZZ_JHU.h"
#include "ZZMatrixElement/MELA/interface/SuperMELA.h"

#include <RooMsgService.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TGraph.h>
#include <vector>

#include <string>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>


using namespace RooFit;

Mela::Mela(int LHCsqrts, float mh) 
{
  int maxSqrts = 8;
//  if(LHCsqrts!=8 && LHCsqrts!=7) assert(0);
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

  mzz_rrv = new RooRealVar("mzz","m_{ZZ}",0.,1000.);
  z1mass_rrv = new RooRealVar("z1mass","m_{Z1}",0.,180.);
  z2mass_rrv = new RooRealVar("z2mass","m_{Z2}",0.,120.); 
  costhetastar_rrv = new RooRealVar("costhetastar","cos#theta^{*}",-1.,1.);  
  costheta1_rrv = new RooRealVar("costheta1","cos#theta_{1}",-1.,1.);  
  costheta2_rrv = new RooRealVar("costheta2","cos#theta_{2}",-1.,1.);
  phi_rrv= new RooRealVar("phi", "#Phi", -TMath::Pi(), TMath::Pi());
  phi1_rrv= new RooRealVar("phi1", "#Phi_{1}", -TMath::Pi(), TMath::Pi());
  
  upFrac_rrv = new RooRealVar("upFrac","fraction up-quarks",.5,0.,1.);

  spin0Model = new AngularPdfFactory(z1mass_rrv,z2mass_rrv,costheta1_rrv,costheta2_rrv,phi_rrv,mzz_rrv);
  spin1Model = new VectorPdfFactory(z1mass_rrv,z2mass_rrv,costhetastar_rrv,costheta1_rrv,costheta2_rrv,phi_rrv,phi1_rrv,mzz_rrv);
  spin2Model = new TensorPdfFactory(z1mass_rrv,z2mass_rrv,costhetastar_rrv,costheta1_rrv,costheta2_rrv,phi_rrv,phi1_rrv,mzz_rrv);
  qqZZmodel = new RooqqZZ_JHU_ZgammaZZ_fast("qqZZmodel","qqZZmodel",*z1mass_rrv,*z2mass_rrv,*costheta1_rrv,*costheta2_rrv,*phi_rrv,*costhetastar_rrv,*phi1_rrv,*mzz_rrv,*upFrac_rrv);

 //edm::FileInPath HiggsWidthFile("Higgs/Higgs_CS_and_Width/txtFiles/8TeV-ggH.txt");
// edm::FileInPath HiggsWidthFile("Higgs/Higgs_CS_and_Width/txtFiles/YR3/8TeV-ggH.txt");
  edm::FileInPath HiggsWidthFile("ZZMatrixElement/MELA/data/HiggsTotalWidth.txt");

  std::string path = HiggsWidthFile.fullPath();
  //std::cout << path.substr(0,path.length()-12) << std::endl;
 // ZZME = new  newZZMatrixElement(path.substr(0,path.length()-12 ).c_str(),1000.*LHCsqrts/2.,TVar::INFO);
  ZZME = new  newZZMatrixElement(path.substr(0,path.length()-19 ).c_str(),1000.*LHCsqrts/2.,TVar::INFO);
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

  char jvbf_path[500];
  sprintf(jvbf_path, "ZZMatrixElement/MELA/data/jvbfMELA_fits_wide_%dTeV.root", superMELA_LHCsqrts); // superMELA_LHCsqrts is correct here.
  edm::FileInPath jvbf_file(jvbf_path);
  TFile* finput_jvbfPint = TFile::Open(jvbf_file.fullPath().c_str(), "read");
  jvbf_Pint_par = (TGraph*)finput_jvbfPint->Get("njets1_pjvbf_pars");
  finput_jvbfPint->Close();
  assert(jvbf_Pint_par);
}

Mela::~Mela(){ 
  //std::cout << "begin destructor" << std::endl;  
  setRemoveLeptonMasses(false); // Use Run I scheme for not removing lepton masses. Notice the switch itself is defined as an extern, so it has to be set to default value at the destructor!

  if (jvbf_Pint_par!=0) delete jvbf_Pint_par;
  for (int i=0; i<3; i++){
    if (tgtotalbkg[i] != 0) delete tgtotalbkg[i];
  }
  if (DggZZ_scalefactor!=0) delete DggZZ_scalefactor;

  delete mzz_rrv;
  delete z1mass_rrv; 
  delete z2mass_rrv; 
  delete costhetastar_rrv;
  delete costheta1_rrv;
  delete costheta2_rrv;
  delete phi_rrv;
  delete phi1_rrv;
  delete upFrac_rrv;

  delete spin0Model;
  delete spin1Model;
  delete spin2Model;
  delete qqZZmodel;
  delete ZZME;
  delete super;
  delete myR;
}

void Mela::setProcess(TVar::Process myModel, TVar::MatrixElement myME, TVar::Production myProduction)
{
  myModel_ = myModel;
  myME_ = myME;
  myProduction_ = myProduction;

  // 
  // configure the analytical calculations 
  // 

  if(myModel_==TVar::bkgZZ)  pdf = qqZZmodel;
  else if(myProduction == TVar::JJGG || myProduction == TVar::JJVBF) ;
  else if(myProduction == TVar::ZH || myProduction == TVar::WH ) ;
  else if(!spin0Model->configure(myModel_)) pdf = spin0Model->PDF;
  else if(!spin1Model->configure(myModel_)) pdf = spin1Model->PDF;
  else if(!spin2Model->configure(myModel_,myProduction_)) pdf = spin2Model->PDF;
  else if(myME_ == TVar::ANALYTICAL)
    cout << "Mela::setProcess -> ERROR TVar::Process not found!!! " << myME_ << " "<< myModel_<<endl; 

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


// Re-order masses and angles as needed by likelihoodDiscriminant. 
// This follows a different convention than the usual Z1/Z2 definition!
void Mela::checkZorder(float& z1mass, float& z2mass,
		       float& costhetastar, float& costheta1,
		       float& costheta2, float& phi, 
		       float& phistar1){

  float tempZ1mass=z1mass;
  float tempZ2mass=z2mass;
  float tempH1=costheta1;
  float tempH2=costheta2;
  float tempHs=costhetastar;
  float tempPhi1=phistar1;
  float tempPhi=phi;

  if(z2mass>z1mass){
    //cout<<"inverted"<<endl;
    z1mass=tempZ2mass;
    z2mass=tempZ1mass;
    costhetastar=-tempHs;
    costheta1=tempH2;
    costheta2=tempH1;
    phi=tempPhi;
    phistar1=-tempPhi1-tempPhi;
    if(phistar1>3.1415)
      phistar1=phistar1-2*Geom::pi();
    if(phistar1<-3.1415)
      phistar1=phistar1+2*Geom::pi();

  }else
    return;

}

std::vector<TLorentzVector> Mela::calculate4Momentum(double Mx, double M1, double M2, double theta, double theta1, double theta2, double Phi1, double Phi){
	return ZZME->Calculate4Momentum(Mx, M1, M2, theta, theta1, theta2, Phi1, Phi);
}


void Mela::computeD_CP(float mZZ, float mZ1, float mZ2, // input kinematics
           float costhetastar,
           float costheta1,
           float costheta2,
           float phi,
           float phi1,
           int flavor,
           TVar::MatrixElement myME,
           TVar::Process myType ,
           float& prob){
  reset_PAux();


/******** No analytical for D_CP_T and D_Int_T, ME has no imaginary part now! Only work for JHUGen *******/
// float pMix, p0plus, p_star;
// TVar::Process mixProcess , starProcess;
//
//if(myType == TVar::D_g1g4) { mixProcess = TVar::CPMixHZZ_4l; starProcess = TVar::PSHZZ_g4star;}
//else if(myType == TVar::D_g1g4_pi_2) { mixProcess = TVar::CPMixHZZ_4l_pi_2; starProcess = TVar::PSHZZ_g4star;}
//else if(myType == TVar::D_g1g2) { mixProcess = TVar::HDMixHZZ_4l; starProcess = TVar::HDHZZ_4l_g2star;}
//else if(myType == TVar::D_g1g2_pi_2) { mixProcess = TVar::HDMixHZZ_4l_pi_2; starProcess = TVar::HDHZZ_4l_g2star;}
//else{
// cout<<"Interaction type not supported!"<<endl;
// return;
//}
// setProcess( mixProcess, myME, TVar::ZZGG);
// computeP(mZZ, mZ1, mZ2,
//        costhetastar,costheta1,costheta2,phi,phi1,flavor,
//        pMix);
//
// setProcess(TVar::HSMHiggs, myME, TVar::ZZGG);
// computeP(mZZ, mZ1, mZ2,
//        costhetastar,costheta1,costheta2,phi,phi1,flavor,
//        p0plus);
//
// setProcess(starProcess,myME, TVar::ZZGG);
// computeP(mZZ, mZ1, mZ2,
//        costhetastar,costheta1,costheta2,phi,phi1,flavor,
//        p_star);
//  prob = pMix- p0plus- p_star;
//

// setProcess(TVar::H0minus, myME, TVar::ZZGG);
// computeP(mZZ, mZ1, mZ2,
//        costhetastar,costheta1,costheta2,phi,phi1,flavor,
//        p0minus);                                        

//prob = (p0minus_plus - p0minus_g4star - p0plus ) / (p0minus + p0plus);

double coupl_mix[SIZE_HVV][2] = {{ 0 }};
double coupl_1[SIZE_HVV][2] = {{ 0 }};
double coupl_2[SIZE_HVV][2] = {{ 0 }};

switch (myType){
	case TVar::D_g1g4 :
		coupl_mix[0][0] =1.;
		coupl_mix[3][0] =2.521;
		coupl_1[0][0] =1.;
		coupl_2[3][0] =2.521;
	break; 
	case TVar::D_g1g4_pi_2 :
		coupl_mix[0][0] =1.;
		coupl_mix[3][1] =2.521;
		coupl_1[0][0] =1.;
		coupl_2[3][1] =2.521;
	break; 
	case TVar::D_g1g2 :
		coupl_mix[0][0] =1.;
		coupl_mix[1][0] = 1.638;
		coupl_1[0][0] =1.;
		coupl_2[1][0] = 1.638;
	break; 
	case TVar::D_g1g2_pi_2 :
		coupl_mix[0][0] =1.;
		coupl_mix[1][1] = 1.638 ;
		coupl_1[0][0] =1.;
		coupl_2[1][1] = 1.638;
	break; 
	case TVar::D_g1g1prime2 :
		coupl_mix[0][0] =1.;
		coupl_mix[11][0] = 12046.01;
		coupl_1[0][0] =1.;
		coupl_2[11][0] = 12046.01;
	break; 
	case TVar::D_zzzg :
		coupl_mix[0][0] =1.;
		coupl_mix[4][0] = 0.0688;
		coupl_1[0][0] =1.;
		coupl_2[4][0] = 0.0688;
	break; 
	case TVar::D_zzgg :
		coupl_mix[0][0] =1.;
		coupl_mix[7][0] = -0.0898;
		coupl_1[0][0] =1.;
		coupl_2[7][0] = -0.0898;
	break; 
	case TVar::D_zzzg_PS :
		coupl_mix[0][0] =1.;
		coupl_mix[6][0] = 0.0855;
		coupl_1[0][0] =1.;
		coupl_2[6][0] = 0.0855;
	break; 
	case TVar::D_zzgg_PS :
		coupl_mix[0][0] =1.;
		coupl_mix[9][0] = -0.0907; 
		coupl_1[0][0] =1.;
		coupl_2[9][0] = -0.0907; 
	break; 
	case TVar::D_zzzg_g1prime2 :
		coupl_mix[0][0] =1.;
		coupl_mix[30][0] = -7591.914; 
		coupl_1[0][0] =1.;
		coupl_2[30][0] = -7591.914; 
	break; 
	case TVar::D_zzzg_g1prime2_pi_2 :
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
        costhetastar,costheta1,costheta2,phi,phi1,flavor, coupl_mix,pMix);

 setProcess(TVar::SelfDefine_spin0, myME, TVar::ZZGG);
 computeP(mZZ, mZ1, mZ2,
        costhetastar,costheta1,costheta2,phi,phi1,flavor, coupl_1,p1);

 setProcess(TVar::SelfDefine_spin0, myME, TVar::ZZGG);
 computeP(mZZ, mZ1, mZ2,
        costhetastar,costheta1,costheta2,phi,phi1,flavor, coupl_2,p2);

  prob = pMix- p1- p2;
}

 
 void Mela::computeP(float mZZ, float mZ1, float mZ2, // input kinematics
 		    float costhetastar,
 		    float costheta1, 
 		    float costheta2,
 		    float phi,
 		    float phi1,
 		    int flavor,
 		    float& prob, bool useConstant){                   // output probability
   reset_PAux();

   //cout << "Mela::computeP - begin" << endl;
   //cout << "calculator: " << myME_ << " model: " << myModel_ << " production: " << myProduction_ << endl;

  costhetastar_rrv->setVal(costhetastar);
  costheta1_rrv->setVal(costheta1);
  costheta2_rrv->setVal(costheta2);
  phi_rrv->setVal(phi);
  phi1_rrv->setVal(phi1);
     
   z1mass_rrv->setVal(mZ1);
   z2mass_rrv->setVal(mZ2);
   mzz_rrv->setVal(mZZ);
 
   //cout << "Mela::computeP() - set RooRealVars" << endl;
   
   float constant = 1;
   double couplingvals_NOTggZZ[SIZE_HVV_FREENORM] = { 0 };
   double selfDHvvcoupl[SIZE_HVV][2] = { { 0 } };
   double selfDGqqcoupl[SIZE_GQQ][2] = { { 0 } }; 
   double selfDGggcoupl[SIZE_GGG][2] = { { 0 } };
   double selfDGvvcoupl[SIZE_GVV][2] = { { 0 } };
   double selfDZqqcoupl[SIZE_ZQQ][2] = { { 0 } };
   double selfDZvvcoupl[SIZE_ZVV][2] = { { 0 } };
   
   //
   // analytical calculations
   // 
   if ( myME_ == TVar::ANALYTICAL ) {
     
     //cout << "Mela::computeP() - analytical calc " << endl;
    
     if(mZZ>100.){

       if(myProduction_==TVar::ZZINDEPENDENT){
	RooAbsPdf* integral = (RooAbsPdf*) pdf->createIntegral(RooArgSet(*costhetastar_rrv,*phi1_rrv));
	prob = integral->getVal();
	delete integral;
      }else{
	prob = pdf->getVal();
      }

    }else{
      prob = -99.0;
    }
    
    //cout << "Mela::computeP() - getting analytical c-constants" << endl;
    
    // 
    // define the constants to be used on analytical
    // 
    
    // gg productions 
    if ( myME_ == TVar::ANALYTICAL && myProduction_ == TVar::ZZGG  ) {
      if ( flavor == 3 ) {
	
	//cout << "ANALYTICAL - GG - flavor=3" << endl;
	
				if ( myModel_ == TVar::H0minus)  constant = 6.4;
				if ( myModel_ == TVar::H0hplus)  constant = 2.2;
				if ( myModel_ == TVar::H2_g1g5)  constant = 9.5;
				if ( myModel_ == TVar::H2_g4 )  constant = 7.3e7;
				if ( myModel_ == TVar::H2_g8 )  constant = 1.1e8;
				if ( myModel_ == TVar::H2_g5 )  constant = 16.3;

				if ( myModel_ == TVar::H2_g2 ) constant = 552385;
				if ( myModel_ == TVar::H2_g3 ) constant = 1.08147e+08;
				if ( myModel_ == TVar::H2_g6 ) constant = 748630;
				if ( myModel_ == TVar::H2_g7 ) constant = 1.2522e+07;
				if ( myModel_ == TVar::H2_g9 ) constant = 3.01467e+06;
				if ( myModel_ == TVar::H2_g10 ) constant = 1.483e+11;

      }  else {

	//cout << "ANALYTICAL - GG - flavor!=3" << endl;

				if ( myModel_ == TVar::H0minus )  constant = 6.5;
				if ( myModel_ == TVar::H0hplus )  constant = 2.2;
				if ( myModel_ == TVar::H2_g1g5 )  constant = 9.3;
				if ( myModel_ == TVar::H2_g4 )  constant = 1.1e8;
				if ( myModel_ == TVar::H2_g8 )  constant = 1.9e8;
				if ( myModel_ == TVar::H2_g5 )  constant = 15.6;

				if ( myModel_ == TVar::H2_g2 ) constant = 552385;
				if ( myModel_ == TVar::H2_g3 ) constant = 1.24897e+08;
				if ( myModel_ == TVar::H2_g6 ) constant = 1.10793e+06;
				if ( myModel_ == TVar::H2_g7 ) constant = 2.21423e+07;
				if ( myModel_ == TVar::H2_g9 ) constant = 3.18193e+06;
				if ( myModel_ == TVar::H2_g10 ) constant = 2.63811e+11;

      }

    }
    // qqb productions 
    if ( myME_ == TVar::ANALYTICAL && myProduction_ == TVar::ZZQQB  ) {
      if ( flavor == 3 ) {

	//cout << "ANALYTICAL - QQB - flavor=3" << endl;

				if ( myModel_ == TVar::H1minus )  constant = 4.6e5;
				if ( myModel_ == TVar::H1plus )  constant = 4.0e5;
				if ( myModel_ == TVar::H2_g1g5 )  constant = 7.9;

				if ( myModel_ == TVar::H2_g5 ) constant = 13.7977;
				if ( myModel_ == TVar::H2_g4 ) constant = 5.12897e+07;
				if ( myModel_ == TVar::H2_g2 ) constant = 477586;
				if ( myModel_ == TVar::H2_g3 ) constant = 1.30907e+08;
				if ( myModel_ == TVar::H2_g6 ) constant = 847461;
				if ( myModel_ == TVar::H2_g7 ) constant = 1.39014e+07;
				if ( myModel_ == TVar::H2_g8 ) constant = 7.08446e+07;
				if ( myModel_ == TVar::H2_g9 ) constant = 2.93583e+06;
				if ( myModel_ == TVar::H2_g10 ) constant = 1.47118e+11;
      } else {

	//cout << "ANALYTICAL - QQB - flavor!=3" << endl;

				if ( myModel_ == TVar::H1minus )  constant = 4.6e5;
				if ( myModel_ == TVar::H1plus )  constant = 4.0e5;
				if ( myModel_ == TVar::H2_g1g5 )  constant = 7.9;

				if ( myModel_ == TVar::H2_g5 ) constant = 13.7289;
				if ( myModel_ == TVar::H2_g4 ) constant = 7.57539e+07;
				if ( myModel_ == TVar::H2_g2 ) constant = 476156;
				if ( myModel_ == TVar::H2_g3 ) constant = 1.44675e+08;
				if ( myModel_ == TVar::H2_g6 ) constant = 1.07303e+06;
				if ( myModel_ == TVar::H2_g7 ) constant = 2.37359e+07;
				if ( myModel_ == TVar::H2_g8 ) constant = 1.35435e+08;
				if ( myModel_ == TVar::H2_g9 ) constant = 2.99514e+06;
				if ( myModel_ == TVar::H2_g10 ) constant = 2.13201e+11;
      }
    }
    // production independent calculations
    if ( myME_ == TVar::ANALYTICAL && myProduction_ == TVar::ZZINDEPENDENT  ) {
      if ( flavor == 3) {

	//cout << "ANALYTICAL - INDEP - flavor=3" << endl;

				if ( myModel_ == TVar::H1minus )  constant = 3.4e4;
				if ( myModel_ == TVar::H1plus )  constant = 3.4e4;
				if ( myModel_ == TVar::H2_g1g5 )  constant = 0.66;

				if ( myModel_ == TVar::H2_g5 ) constant = 1.15604;
				if ( myModel_ == TVar::H2_g4 ) constant = 4.36662e+06;
				if ( myModel_ == TVar::H2_g2 ) constant = 39994.6;
				if ( myModel_ == TVar::H2_g3 ) constant = 1.0897e+07;
				if ( myModel_ == TVar::H2_g6 ) constant = 61420.6;
				if ( myModel_ == TVar::H2_g7 ) constant = 1.20742e+06;
				if ( myModel_ == TVar::H2_g8 ) constant = 6.07991e+06;
				if ( myModel_ == TVar::H2_g9 ) constant = 239187;
				if ( myModel_ == TVar::H2_g10 ) constant = 1.1843e+10;
      } else {

	//cout << "ANALYTICAL - INDEP - flavor!=3" << endl;

				if ( myModel_ == TVar::H1minus )  constant = 3.4e4;
				if ( myModel_ == TVar::H1plus )  constant = 3.4e4;
				if ( myModel_ == TVar::H2_g1g5 )  constant = .66;

			if ( myModel_ == TVar::H2_g5 ) constant = 1.15604;
			if ( myModel_ == TVar::H2_g4 ) constant = 5.6237e+06;
				if ( myModel_ == TVar::H2_g2 ) constant = 39715.6;
				if ( myModel_ == TVar::H2_g3 ) constant = 1.16172e+07;
				if ( myModel_ == TVar::H2_g6 ) constant = 77613.8;
				if ( myModel_ == TVar::H2_g7 ) constant = 1.58485e+06;
				if ( myModel_ == TVar::H2_g8 ) constant = 8.71451e+06;
				if ( myModel_ == TVar::H2_g9 ) constant = 241591;
				if ( myModel_ == TVar::H2_g10 ) constant = 1.55139e+10;

      }
    } 


  }

  //cout << "Mela::computeP() - I am out!!" << endl;

  //
  // JHUGen or MCFM 
  //
  if ( myME_ == TVar::JHUGen || myME_ == TVar::MCFM ) {

    //cout << "Mela::computeP() - JHUGen/MCFM calc " << endl;

    //initialize variables
    checkZorder(mZ1,mZ2,costhetastar,costheta1,costheta2,phi,phi1);
    ZZME->computeXS(mZZ,mZ1,mZ2,
		    costhetastar,costheta1,costheta2, 
		    phi, phi1, flavor,
		    myModel_, myME_,  myProduction_, couplingvals_NOTggZZ, selfDHvvcoupl,
           selfDZqqcoupl,
           selfDZvvcoupl,
           selfDGqqcoupl,
           selfDGggcoupl,
           selfDGvvcoupl,prob);

    //cout << "Mela::computeP() - mZZ: " << mZZ << endl;
    //cout << "Mela::computeP() - mZ1: " << mZ1 << endl;
    //cout << "Mela::computeP() - mZ2: " << mZ2 << endl;
    //cout << "Mela::computeP() - costheta1: " << costheta1 << endl;
    //cout << "Mela::computeP() - costheta2: " << costheta2 << endl;
    //cout << "Mela::computeP() - costhetastar: " << costhetastar << endl;
    //cout << "Mela::computeP() - phi: " << phi << endl;
    //cout << "Mela::computeP() - phi1: " << phi1 << endl;    
    //cout << "Mela::computeP() - prob: " << prob << endl;

    // adding scale factors for MCMF calculation
    // -- taken from old code --

    // note: constants are being added to ggZZ ME calculation 
    // for the purpose of building DggZZ to separate qqZZ from
    // ggZZ.  These constants have been tune on the unadulterated
    // qqZZ ME calculation, so the qqZZ scale factors should be 
    // included inorder to cancel those in the qqZZ ME calc

    // 4e scale factors
if (useConstant){
    if(flavor==1 && myME_ == TVar::MCFM){

      // for ggZZ 
 //     if(myProduction_ == TVar::ZZGG){

 // if(mZZ > 900)
 //   prob *=vaScale_4e->Eval(900.)/DggZZ_scalefactor->Eval(900.);
 // else if (mZZ <  110 )
 //   prob *=vaScale_4e->Eval(110.)/DggZZ_scalefactor->Eval(110.);
 // else
 //   prob *=vaScale_4e->Eval(mZZ)/DggZZ_scalefactor->Eval(mZZ);
 // 
 //     }// end GG

      // for qqZZ
      if(myProduction_ == TVar::ZZQQB){

	if(mZZ > 900)
	  prob *= vaScale_4e->Eval(900.);
	else if (mZZ <  70 )
	  prob *= vaScale_4e->Eval(70.);
	else
	  prob *= vaScale_4e->Eval(mZZ);

      }// end QQB 

    }// end 4e scale factors

    // 4mu scale factors
    if(flavor==2 && myME_ == TVar::MCFM){
      
      // for ggZZ 
//      if(myProduction_ == TVar::ZZGG){
//	if(mZZ > 900)                   
//	  prob *=vaScale_4mu->Eval(900.)/DggZZ_scalefactor->Eval(900.);
//	else if (mZZ <  110 )
//	  prob *=vaScale_4mu->Eval(110.)/DggZZ_scalefactor->Eval(110.);
//	else
//	  prob *=vaScale_4mu->Eval(mZZ)/DggZZ_scalefactor->Eval(mZZ);
//      }// end GG
            
      // for qqZZ
      if(myProduction_ == TVar::ZZQQB){
	if(mZZ > 900)                   
	  prob *= vaScale_4mu->Eval(900.);
	else if (mZZ <  70 )
	  prob *= vaScale_4mu->Eval(70.);
	else
	  prob *= vaScale_4mu->Eval(mZZ);
      }// end qqZZ

    }// end 4mu scale factors

    // 2e2mu scale factors
    if(flavor==3 && myME_ == TVar::MCFM){
  //    cout<<"BEFORE scale: "<<prob<<endl;

      // for ggZZ 
//      if(myProduction_ == TVar::ZZGG){
//	if(mZZ > 900) 
//	  prob *=vaScale_2e2mu->Eval(900.)/DggZZ_scalefactor->Eval(900.);
//      else if (mZZ <  110 )
//	  prob *=vaScale_2e2mu->Eval(110.)/DggZZ_scalefactor->Eval(110.);
//      else
//	  prob *=vaScale_2e2mu->Eval(mZZ)/DggZZ_scalefactor->Eval(mZZ);
//
//      }// end GG

      // for qqZZ
      if(myProduction_ == TVar::ZZQQB){
	if(mZZ > 900)                   
	  prob *= vaScale_2e2mu->Eval(900.);
	else if (mZZ <  70 )
	  prob *= vaScale_2e2mu->Eval(70.);
	else
	  prob *= vaScale_2e2mu->Eval(mZZ);
      }// end qqZZ

    }// end 2e2mu scale factors
 }
    //cout << "Mela::computeP() - getting JHUGen c-constants" << endl;

    // 
    // define the constants to be used on JHUGen
    // 
    
    // gg productions 
    if ( myME_ == TVar::JHUGen && myProduction_ == TVar::ZZGG  ) {
      if ( flavor == 3 ) {
				if ( myModel_ == TVar::H0minus )  constant = 6.0;
				if ( myModel_ == TVar::H0hplus )  constant = 2.1;
				if ( myModel_ == TVar::H2_g1g5 )  constant = 0.6;
				if ( myModel_ == TVar::H2_g4 )  constant = 2.7e10;
				if ( myModel_ == TVar::H2_g8 )  constant = 4.1e10;
				if ( myModel_ == TVar::H2_g5 )  constant = .97;

				if ( myModel_ == TVar::H2_g2 ) constant = 4.74608e+08;
				if ( myModel_ == TVar::H2_g3 ) constant = 1.63324e+11;
				if ( myModel_ == TVar::H2_g6 ) constant = 46816.9;
				if ( myModel_ == TVar::H2_g7 ) constant = 783088;
				if ( myModel_ == TVar::H2_g9 ) constant = 1.13853e+09;
				if ( myModel_ == TVar::H2_g10 ) constant = 5.58394e+13;

      }  else {
				if ( myModel_ == TVar::H0minus )  constant = 7.0;
				if ( myModel_ == TVar::H0hplus )  constant = 2.3;
				if ( myModel_ == TVar::H2_g1g5 )  constant = 1.4/2.;
				if ( myModel_ == TVar::H2_g4 )  constant = 2.6e10;
				if ( myModel_ == TVar::H2_g8 )  constant = 3.7e10;
				if ( myModel_ == TVar::H2_g5 )  constant = 1.26;

				if ( myModel_ == TVar::H2_g2 ) constant = 5.96148e+08;
				if ( myModel_ == TVar::H2_g3 ) constant = 1.95534e+11;
				if ( myModel_ == TVar::H2_g6 ) constant = 55160.4;
				if ( myModel_ == TVar::H2_g7 ) constant = 658026;
				if ( myModel_ == TVar::H2_g9 ) constant = 1.23089e+09;
				if ( myModel_ == TVar::H2_g10 ) constant = 5.78862e+13;

			}
    } 
    // qqb productions 
    if ( myME_ == TVar::JHUGen && myProduction_ == TVar::ZZQQB  ) {
      if ( flavor == 3 ) {
				if ( myModel_ == TVar::H1minus )  constant = 16.;
				if ( myModel_ == TVar::H1plus )  constant = 13.;
				if ( myModel_ == TVar::H2_g1g5 )  constant = 13.;

				if ( myModel_ == TVar::H2_g5 ) constant = 22.8625;
				if ( myModel_ == TVar::H2_g4 ) constant = 8.49013e+07;
				if ( myModel_ == TVar::H2_g2 ) constant = 792938;
				if ( myModel_ == TVar::H2_g3 ) constant = 2.17563e+08;
				if ( myModel_ == TVar::H2_g6 ) constant = 1.40845e+06;
				if ( myModel_ == TVar::H2_g7 ) constant = 2.31499e+07;
				if ( myModel_ == TVar::H2_g8 ) constant = 1.17271e+08;
				if ( myModel_ == TVar::H2_g9 ) constant = 4.86462e+06;
				if ( myModel_ == TVar::H2_g10 ) constant = 2.43529e+11;
      } else {
				if ( myModel_ == TVar::H1minus )  constant = 38/2.;
				if ( myModel_ == TVar::H1plus )  constant = 28/2.;
				if ( myModel_ == TVar::H2_g1g5 )  constant = 30/2.;

				if ( myModel_ == TVar::H2_g5 ) constant = 26.167;
				if ( myModel_ == TVar::H2_g4 ) constant = 8.9612e+07;
				if ( myModel_ == TVar::H2_g2 ) constant = 984117;
				if ( myModel_ == TVar::H2_g3 ) constant = 2.46285e+08;
				if ( myModel_ == TVar::H2_g6 ) constant = 1.6201e+06;
				if ( myModel_ == TVar::H2_g7 ) constant = 1.95307e+07;
				if ( myModel_ == TVar::H2_g8 ) constant = 1.29087e+08;
				if ( myModel_ == TVar::H2_g9 ) constant = 5.11404e+06;
				if ( myModel_ == TVar::H2_g10 ) constant = 2.4183e+11;
      }
    }
    // production independent calculations
    if ( myME_ == TVar::JHUGen && myProduction_ == TVar::ZZINDEPENDENT  ) {
      if ( flavor == 3) {
				if ( myModel_ == TVar::H1minus )  constant = 1.3e+10;
				if ( myModel_ == TVar::H1plus )  constant = 1.3e+10;
				if ( myModel_ == TVar::H2_g1g5 )  constant = 1.6e+9;

				if ( myModel_ == TVar::H2_g5 ) constant = 2.71213e+09;
				if ( myModel_ == TVar::H2_g4 ) constant = 1.01932e+16;
				if ( myModel_ == TVar::H2_g2 ) constant = 9.29888e+13;
				if ( myModel_ == TVar::H2_g3 ) constant = 2.53613e+16;
				if ( myModel_ == TVar::H2_g6 ) constant = 1.43664e+14;
				if ( myModel_ == TVar::H2_g7 ) constant = 2.81011e+15;
				if ( myModel_ == TVar::H2_g8 ) constant = 1.40936e+16;
				if ( myModel_ == TVar::H2_g9 ) constant = 5.57788e+14;
				if ( myModel_ == TVar::H2_g10 ) constant = 2.73432e+19;
      } else {
				if ( myModel_ == TVar::H1minus )  constant = 1.6e+10;
				if ( myModel_ == TVar::H1plus )  constant = 1.4e+10;
				if ( myModel_ == TVar::H2_g1g5 )  constant = 2.0e+9;

				if ( myModel_ == TVar::H2_g5 ) constant = 3.30598e+09;
				if ( myModel_ == TVar::H2_g4 ) constant = 1.01932e+16;
				if ( myModel_ == TVar::H2_g2 ) constant = 1.16336e+14;
				if ( myModel_ == TVar::H2_g3 ) constant = 2.76942e+16;
				if ( myModel_ == TVar::H2_g6 ) constant = 1.68929e+14;
				if ( myModel_ == TVar::H2_g7 ) constant = 2.22827e+15;
				if ( myModel_ == TVar::H2_g8 ) constant = 1.39674e+16;
				if ( myModel_ == TVar::H2_g9 ) constant = 5.63394e+14;
				if ( myModel_ == TVar::H2_g10 ) constant = 2.54946e+19;
      }
    } 
    
    // ***
    // experimental for the ZZ decay 
    // ****
    
    //cout << "Mela::computeP() - production indep MCFM" << endl;
    
    if ( myME_ == TVar::MCFM 
	 && myProduction_ == TVar::ZZINDEPENDENT 
	 &&  myModel_ == TVar::bkgZZ
	 )
      {
	prob = 0;
	int gridsize_hs = 5; 
	double hs_min = 0; //-1.;
	double hs_max = 1;
	double hs_step = ( hs_max - hs_min ) / double (gridsize_hs); 
	
	int gridsize_phi1 = 5; 
	double phi1_min = 0; //-TMath::Pi();
	double phi1_max = TMath::Pi();
	double phi1_step = ( phi1_max - phi1_min ) / double (gridsize_phi1); 
	
	for ( int i_hs = 0; i_hs < gridsize_hs + 1; i_hs ++ ) {
	  double hs_val = hs_min + i_hs * hs_step; 
	  for ( int i_phi1 = 0; i_phi1 < gridsize_phi1 +1 ; i_phi1 ++ ) {
	    double phi1_val = phi1_min + i_phi1 * phi1_step; 
	    float temp_prob(0.); 
	    // calculate the ZZ using MCFM
	    ZZME->computeXS(mZZ,mZ1,mZ2,
			    hs_val,costheta1,costheta2, 
			    phi, phi1_val, flavor,
			    myModel_, myME_,  myProduction_, couplingvals_NOTggZZ,selfDHvvcoupl, 
           selfDZqqcoupl,
           selfDZvvcoupl,
           selfDGqqcoupl,
           selfDGggcoupl,
           selfDGvvcoupl,temp_prob);
	    prob += temp_prob;
	  }
	}
	prob =  prob / float ( (gridsize_hs + 1) * (gridsize_phi1 +1 )); 
      
	// adding scale factors for MCMF calculation
	// -- taken from old code --
if(useConstant){	
	if(flavor==1){
	  if(mZZ > 900)                   
	    prob *= vaScale_4e->Eval(900.);
	  else if (mZZ <  70 )
	    prob *= vaScale_4e->Eval(70.);
	  else
	    prob *= vaScale_4e->Eval(mZZ);
	}
	
	if(flavor==2){
	  if(mZZ > 900)                   
	    prob *= vaScale_4mu->Eval(900.);
	  else if (mZZ <  70 )
	    prob *= vaScale_4mu->Eval(70.);
	  else
	    prob *= vaScale_4mu->Eval(mZZ);
	}
	
	if(flavor==3){
	  if(mZZ > 900)                   
	    prob *= vaScale_2e2mu->Eval(900.);
	  else if (mZZ <  70 )
	    prob *= vaScale_2e2mu->Eval(70.);
	  else
	    prob *= vaScale_2e2mu->Eval(mZZ);
	}
}	
      }
  }

if(useConstant)
  prob *= constant; 

}

void Mela::computeP(float mZZ, float mZ1, float mZ2, // input kinematics
		    float costhetastar,
		    float costheta1, 
		    float costheta2,
		    float phi,
		    float phi1,
		    int flavor,
		    double selfDHvvcoupl[SIZE_HVV][2],
		    float& prob){ 
   reset_PAux();

   double couplingvals_NOTggZZ[SIZE_HVV_FREENORM] = {0,0};
   double selfDGqqcoupl[SIZE_GQQ][2]= {{0}};
   double selfDGggcoupl[SIZE_GGG][2]= {{0}};
   double selfDGvvcoupl[SIZE_GVV][2]= {{0}};
   double selfDZqqcoupl[SIZE_ZQQ][2]= {{0}};
   double selfDZvvcoupl[SIZE_ZVV][2]= {{0}};

 if ( !( myModel_ == TVar::SelfDefine_spin0 || myME_ == TVar::MCFM ) ){
	cout << " Error: This method only applies to spin0, set Process to SelfDefine or ME to MCFM!"<<endl;
  return;
}
  if ( myME_ == TVar::JHUGen || myME_ == TVar::MCFM){ 
    //initialize variables
    checkZorder(mZ1,mZ2,costhetastar,costheta1,costheta2,phi,phi1);
    ZZME->computeXS(mZZ,mZ1,mZ2,
        costhetastar,costheta1,costheta2,
        phi, phi1, flavor,
        myModel_, myME_,  myProduction_, couplingvals_NOTggZZ,     
	 	selfDHvvcoupl,
           selfDZqqcoupl,
           selfDZvvcoupl,
           selfDGqqcoupl,
           selfDGggcoupl,
           selfDGvvcoupl, prob);
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

 for (int i =0 ;i<SIZE_HVV;i++){
		if(selfDHvvcoupl[i][1]!=0){
			cout << "Error: MELA does not support complex coupling for the moment! "<<endl;
			return;
		}
	}
	spin0Model->useGTerm->setVal(1);
	spin0Model->modelIndex =-1;

	spin0Model->g1Val->setVal(selfDHvvcoupl[0][0]); 
	spin0Model->g2Val->setVal(selfDHvvcoupl[1][0]); 
	spin0Model->g3Val->setVal(selfDHvvcoupl[2][0]); 
	spin0Model->g4Val->setVal(selfDHvvcoupl[3][0]); 

	spin0Model->g1_primeVal->setVal(selfDHvvcoupl[10][0]); 
	spin0Model->g1_prime2Val->setVal(selfDHvvcoupl[11][0]); 
	spin0Model->g1_prime3Val->setVal(selfDHvvcoupl[12][0]); 
	spin0Model->g1_prime4Val->setVal(selfDHvvcoupl[13][0]); 
	spin0Model->g1_prime5Val->setVal(selfDHvvcoupl[14][0]); 

	spin0Model->g2_primeVal->setVal(selfDHvvcoupl[15][0]); 
	spin0Model->g2_prime2Val->setVal(selfDHvvcoupl[16][0]); 
	spin0Model->g2_prime3Val->setVal(selfDHvvcoupl[17][0]); 
	spin0Model->g2_prime4Val->setVal(selfDHvvcoupl[18][0]); 
	spin0Model->g2_prime5Val->setVal(selfDHvvcoupl[19][0]); 

	spin0Model->g3_primeVal->setVal(selfDHvvcoupl[20][0]); 
	spin0Model->g3_prime2Val->setVal(selfDHvvcoupl[21][0]); 
	spin0Model->g3_prime3Val->setVal(selfDHvvcoupl[22][0]); 
	spin0Model->g3_prime4Val->setVal(selfDHvvcoupl[23][0]); 
	spin0Model->g3_prime5Val->setVal(selfDHvvcoupl[24][0]); 

	spin0Model->g4_primeVal->setVal(selfDHvvcoupl[25][0]); 
	spin0Model->g4_prime2Val->setVal(selfDHvvcoupl[26][0]); 
	spin0Model->g4_prime3Val->setVal(selfDHvvcoupl[27][0]); 
	spin0Model->g4_prime4Val->setVal(selfDHvvcoupl[28][0]); 
	spin0Model->g4_prime5Val->setVal(selfDHvvcoupl[29][0]); 

	spin0Model->g1_prime6Val->setVal(selfDHvvcoupl[31][0]); 
	spin0Model->g1_prime7Val->setVal(selfDHvvcoupl[32][0]); 
	spin0Model->g2_prime6Val->setVal(selfDHvvcoupl[33][0]); 
	spin0Model->g2_prime7Val->setVal(selfDHvvcoupl[34][0]); 
	spin0Model->g3_prime6Val->setVal(selfDHvvcoupl[35][0]); 
	spin0Model->g3_prime7Val->setVal(selfDHvvcoupl[36][0]); 
	spin0Model->g4_prime6Val->setVal(selfDHvvcoupl[37][0]); 
	spin0Model->g4_prime7Val->setVal(selfDHvvcoupl[38][0]); 

 if(myProduction_==TVar::ZZINDEPENDENT){
  RooAbsPdf* integral = (RooAbsPdf*) pdf->createIntegral(RooArgSet(*costhetastar_rrv,*phi1_rrv));
  prob = integral->getVal();
  delete integral;
      }else{
  prob = pdf->getVal();
      }
}

else cout<<"ERROR: this method only works for JHUGen or MCFM or ANALYTICAL";

}
void Mela::computeP_selfDspin2(float mZZ, float mZ1, float mZ2, // input kinematics
        float costhetastar,
        float costheta1,
        float costheta2,
        float phi,
        float phi1,
        int flavor,
        double selfDGggcoupl[SIZE_GGG][2],
        double selfDGvvcoupl[SIZE_GVV][2], 
        float& prob){
   reset_PAux();

   double couplingvals_NOTggZZ[SIZE_HVV_FREENORM] = { 0 };
   double selfDHvvcoupl[SIZE_HVV][2] = { { 0 } };
   double selfDGqqcoupl[SIZE_GQQ][2] = { { 0 } };
   double selfDZqqcoupl[SIZE_ZQQ][2] = { { 0 } };
   double selfDZvvcoupl[SIZE_ZVV][2] = { { 0 } };

	 if(myModel_ != TVar::SelfDefine_spin2){
	  cout << " Error: This method only applies to spin2, set Process to SelfDefine_spin2!"<<endl;
	  return;
	 }
  if ( myME_ == TVar::JHUGen){


    //initialize variables
    checkZorder(mZ1,mZ2,costhetastar,costheta1,costheta2,phi,phi1);
    ZZME->computeXS(mZZ,mZ1,mZ2,
        costhetastar,costheta1,costheta2,
        phi, phi1, flavor,
        myModel_, myME_,  myProduction_, couplingvals_NOTggZZ,
           selfDHvvcoupl,
           selfDZqqcoupl,
           selfDZvvcoupl,
           selfDGqqcoupl,
           selfDGggcoupl,
           selfDGvvcoupl, prob);
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

 for (int i =0 ;i<SIZE_GVV;i++){
    if(selfDGvvcoupl[i][1]!=0){
      cout << "Error: MELA does not support complex coupling for the moment! "<<endl;
      return;
    }
  }
	
if(myProduction_ == TVar::ZZGG || myProduction_==TVar::ZZINDEPENDENT){
  spin2Model->fz1Val->setVal(0.);
  spin2Model->fz2Val->setVal(1.);
}
if(myProduction_ == TVar::ZZQQB ){
  spin2Model->fz1Val->setVal(1.);
  spin2Model->fz2Val->setVal(0.);
}
  spin2Model->g1Val->setVal(selfDGvvcoupl[0][0]);
  spin2Model->g2Val->setVal(selfDGvvcoupl[1][0]);
  spin2Model->g3Val->setVal(selfDGvvcoupl[2][0]);
  spin2Model->g4Val->setVal(selfDGvvcoupl[3][0]);
  spin2Model->g5Val->setVal(selfDGvvcoupl[4][0]);
  spin2Model->g6Val->setVal(selfDGvvcoupl[5][0]);
  spin2Model->g7Val->setVal(selfDGvvcoupl[6][0]);
  spin2Model->g8Val->setVal(selfDGvvcoupl[7][0]);
  spin2Model->g9Val->setVal(selfDGvvcoupl[8][0]);
  spin2Model->g10Val->setVal(selfDGvvcoupl[9][0]);
	
	spin2Model->calculatefz2();

 if(myProduction_==TVar::ZZINDEPENDENT){
  RooAbsPdf* integral = (RooAbsPdf*) pdf->createIntegral(RooArgSet(*costhetastar_rrv,*phi1_rrv));
  prob = integral->getVal();
  delete integral;
      }else{
  prob = pdf->getVal();
      }
}

else cout<<"ERROR: this method only works for JHUGen and ANALYTICAL";

}
void Mela::computeP_selfDspin1(float mZZ, float mZ1, float mZ2, // input kinematics
        float costhetastar,
        float costheta1,
        float costheta2,
        float phi,
        float phi1,
        int flavor,
        double selfDZvvcoupl[SIZE_ZVV][2],
        float& prob){
  reset_PAux();

	double couplingvals_NOTggZZ[SIZE_HVV_FREENORM] = { 0 };
	double selfDHvvcoupl[SIZE_HVV][2] = { { 0 } };
	double selfDGqqcoupl[SIZE_GQQ][2] = { { 0 } }; 
	double selfDGggcoupl[SIZE_GGG][2] = { { 0 } };
	double selfDGvvcoupl[SIZE_GVV][2] = { { 0 } };
	double selfDZqqcoupl[SIZE_ZQQ][2] = { { 0 } };

 if(myModel_ != TVar::SelfDefine_spin1){
  cout << " Error: This method only applies to spin1, set Process to SelfDefine_spin1!"<<endl;
  return;
 }
  if ( myME_ == TVar::JHUGen){


    //initialize variables
    checkZorder(mZ1,mZ2,costhetastar,costheta1,costheta2,phi,phi1);
    ZZME->computeXS(mZZ,mZ1,mZ2,
        costhetastar,costheta1,costheta2,
        phi, phi1, flavor,
        myModel_, myME_,  myProduction_, couplingvals_NOTggZZ, 
		selfDHvvcoupl,
        selfDZqqcoupl,
        selfDZvvcoupl,
        selfDGqqcoupl,
        selfDGggcoupl,
        selfDGvvcoupl,
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

 for (int i =0 ;i<SIZE_ZVV;i++){
    if(selfDZvvcoupl[i][1]!=0){
      cout << "Error: MELA does not support complex coupling for the moment! "<<endl;
      return;
    }
  }
  spin1Model->g1Val->setVal(selfDZvvcoupl[0][0]);
  spin1Model->g2Val->setVal(selfDZvvcoupl[1][0]);
 if(myProduction_==TVar::ZZINDEPENDENT){
  RooAbsPdf* integral = (RooAbsPdf*) pdf->createIntegral(RooArgSet(*costhetastar_rrv,*phi1_rrv));
  prob = integral->getVal();
  delete integral;
      }else{
  prob = pdf->getVal();
      }
}

else cout<<"ERROR: this method only works for JHUGen and ANALYTICAL";

}
void Mela::computeP(float mZZ, float mZ1, float mZ2, // input kinematics
		    float costhetastar,
		    float costheta1, 
		    float costheta2,
		    float phi,
		    float phi1,
		    int flavor,
		    double couplingvals[SIZE_HVV_FREENORM],
		    float& prob){
  reset_PAux();

	double selfDHvvcoupl[SIZE_HVV][2] = { { 0 } };
	double selfDGqqcoupl[SIZE_GQQ][2] = { { 0 } }; 
	double selfDGggcoupl[SIZE_GGG][2] = { { 0 } };
	double selfDGvvcoupl[SIZE_GVV][2] = { { 0 } };
	double selfDZqqcoupl[SIZE_ZQQ][2] = { { 0 } };
	double selfDZvvcoupl[SIZE_ZVV][2] = { { 0 } };

  if ( (myME_==TVar::JHUGen || myME_==TVar::MCFM)&& myModel_==TVar::bkgZZ_SMHiggs ){
    //initialize variables
    checkZorder(mZ1,mZ2,costhetastar,costheta1,costheta2,phi,phi1);
    ZZME->computeXS(mZZ,mZ1,mZ2,
		    costhetastar,costheta1,costheta2, 
		    phi, phi1, flavor,
		    myModel_, myME_,  myProduction_, couplingvals,selfDHvvcoupl, 
           selfDZqqcoupl,
           selfDZvvcoupl,
           selfDGqqcoupl,
           selfDGggcoupl,
           selfDGvvcoupl,prob);

    // note: constants are being added to ggZZ ME calculation 
    // for the purpose of building DggZZ to separate qqZZ from
    // ggZZ.  These constants have been tune on the unadulterated
    // qqZZ ME calculation, so the qqZZ scale factors should be 
    // included inorder to cancel those in the qqZZ ME calc

    // 4e scale factors
    if(flavor==1){
      // for ggZZ 
      if(myProduction_ == TVar::ZZGG){

	if(mZZ > 900)
	  prob *=vaScale_4e->Eval(900.)/DggZZ_scalefactor->Eval(900.);
	else if (mZZ <  110 )
	  prob *=vaScale_4e->Eval(110.)/DggZZ_scalefactor->Eval(110.);
	else
	  prob *=vaScale_4e->Eval(mZZ)/DggZZ_scalefactor->Eval(mZZ);
	
      }// end GG
    }// end 4e scale factors

    // 4mu scale factors
    if(flavor==2){
      // for ggZZ 
      if(myProduction_ == TVar::ZZGG){
	if(mZZ > 900)                   
	  prob *=vaScale_4mu->Eval(900.)/DggZZ_scalefactor->Eval(900.);
	else if (mZZ <  110 )
	  prob *=vaScale_4mu->Eval(110.)/DggZZ_scalefactor->Eval(110.);
	else
	  prob *=vaScale_4mu->Eval(mZZ)/DggZZ_scalefactor->Eval(mZZ);
      }// end GG
    }// end 4mu scale factors

    // 2e2mu scale factors
    if(flavor==3){
//      cout<<"BEFORE scale: "<<prob<<endl;

      // for ggZZ 
      if(myProduction_ == TVar::ZZGG){
	if(mZZ > 900) 
	  prob *=vaScale_2e2mu->Eval(900.)/DggZZ_scalefactor->Eval(900.);
      else if (mZZ <  110 )
	  prob *=vaScale_2e2mu->Eval(110.)/DggZZ_scalefactor->Eval(110.);
      else
	  prob *=vaScale_2e2mu->Eval(mZZ)/DggZZ_scalefactor->Eval(mZZ);

      }// end GG
    }// end 2e2mu scale factors
  }
//	else if (myME_==TVar::MCFM ){ //new MCFM calculation
//
//  }
  else{
    computeP(mZZ, mZ1, mZ2, costhetastar,costheta1, costheta2,phi,phi1,flavor,prob);
  }
}

void Mela::computeP(TLorentzVector Z1_lept1, int Z1_lept1Id,  // input 4-vectors
		    TLorentzVector Z1_lept2, int Z1_lept2Id,  // 
		    TLorentzVector Z2_lept1, int Z2_lept1Id,
		    TLorentzVector Z2_lept2, int Z2_lept2Id,
		    float& prob){                             // output probability
    
  // get flavor type
  // NEED TO INCLUDE SOME PROTECTION SO THAT USER CANT                  
  // PASS FOUR-VECTORS IN WRONG ORDER.  FOR NOW ASSUMING                
  // THEY ARE PASSED AS e-,e+,mu-,mu+                                   
  // ------------------ channel ------------------------                

  if (mela::forbidMassiveLeptons){
    mela::constrainedRemoveLeptonMass(Z1_lept1,Z1_lept2);
    mela::constrainedRemoveLeptonMass(Z2_lept1,Z2_lept2);
  }

  int flavor;

  if(abs(Z1_lept1Id)==abs(Z1_lept2Id)&&
     abs(Z1_lept1Id)==abs(Z2_lept1Id)&&
     abs(Z1_lept1Id)==abs(Z2_lept2Id)){

    if(abs(Z1_lept1Id)==11) flavor=1;
    else flavor=2;

  }else flavor=3;

  //compute angles  
  float m1=(Z1_lept1 + Z1_lept2).M();
  float m2=(Z2_lept1 + Z2_lept2).M();
    
  TLorentzVector ZZ = (Z1_lept1 + Z1_lept2 + Z2_lept1 + Z2_lept2);
  float mzz = ZZ.M();
    
  // Skip candidates where KD is irrelevant.
//  if (mzz<100.){
//    prob = -99.0;
//    return;
//  }

  float costhetastar, costheta1, costheta2, phi, phi1;

  mela::computeAngles(Z1_lept1, Z1_lept1Id, Z1_lept2, Z1_lept2Id, 
		      Z2_lept1, Z2_lept1Id, Z2_lept2, Z2_lept2Id,
		      costhetastar,costheta1,costheta2,phi,phi1);

  computeP(mzz, m1, m2,
	   costhetastar,
	   costheta1,
	   costheta2,
	   phi,phi1,
	   flavor,
	   prob);

}

void Mela::computeP(TLorentzVector Z1_lept1, int Z1_lept1Id,  // input 4-vectors
		    TLorentzVector Z1_lept2, int Z1_lept2Id,  // 
		    TLorentzVector Z2_lept1, int Z2_lept1Id,
		    TLorentzVector Z2_lept2, int Z2_lept2Id,
		    double couplingvals[SIZE_HVV_FREENORM],
		    float& prob){                             // output probability

  if (mela::forbidMassiveLeptons){
    mela::constrainedRemoveLeptonMass(Z1_lept1,Z1_lept2);
    mela::constrainedRemoveLeptonMass(Z2_lept1,Z2_lept2);
  }

  if(myME_==TVar::JHUGen && myModel_==TVar::bkgZZ_SMHiggs){
    int flavor;
    
    if(abs(Z1_lept1Id)==abs(Z1_lept2Id)&&
       abs(Z1_lept1Id)==abs(Z2_lept1Id)&&
       abs(Z1_lept1Id)==abs(Z2_lept2Id)){
      
      if(abs(Z1_lept1Id)==11) flavor=1;
      else flavor=2;
      
    }else flavor=3;
    
    //compute angles  
    float m1=(Z1_lept1 + Z1_lept2).M();
    float m2=(Z2_lept1 + Z2_lept2).M();
    
    TLorentzVector ZZ = (Z1_lept1 + Z1_lept2 + Z2_lept1 + Z2_lept2);
    float mzz = ZZ.M();
    
    // Skip candidates where KD is irrelevant.
   // if (mzz<100.){
   //   prob = -99.0;
   //   return;
   // }
    
    float costhetastar, costheta1, costheta2, phi, phi1;
    
    mela::computeAngles(Z1_lept1, Z1_lept1Id, Z1_lept2, Z1_lept2Id, 
			Z2_lept1, Z2_lept1Id, Z2_lept2, Z2_lept2Id,
			costhetastar,costheta1,costheta2,phi,phi1);
    
    computeP(mzz, m1, m2,
	     costhetastar,
	     costheta1,
	     costheta2,
	     phi,phi1,
	     flavor,couplingvals,
	     prob);
  }
  else{
  }
  
}

void Mela::computeProdP(TLorentzVector Jet1, int Jet1_Id,
        TLorentzVector Jet2, int Jet2_Id,
        TLorentzVector Decay1, int Decay1_Id,
        TLorentzVector Decay2, int Decay2_Id,
        double selfDHggcoupl[SIZE_HGG][2],
        double selfDHvvcoupl[SIZE_HVV_VBF][2],
        double selfDHwwcoupl[SIZE_HWW_VBF][2],
        float& prob){
  reset_PAux();
  float constant=1.;
  TLorentzVector nullFourVector(0, 0, 0, 0);
  TLorentzVector jet1massless(0, 0, 0, 0);
  TLorentzVector jet2massless(0, 0, 0, 0);
  TLorentzVector higgs;
  if (Decay2==nullFourVector || Decay2_Id==0){
    if (Decay1_Id==25) higgs=Decay1;
    if (Decay1_Id!=25){
      cout<<"No Higgs event passed. Returning prob=-99."<<endl;
    }
  }
  else{
    higgs=Decay1+Decay2;
  }

  if (Jet1==nullFourVector && !(Jet2==nullFourVector)){
    TLorentzVector JetTmp = Jet2;
    Jet2=Jet1;
    Jet1=JetTmp;
  }
  mela::computeJetMassless(Jet1, jet1massless);
  if (!(Jet2==nullFourVector) && myProduction_ != TVar::JH){
    mela::computeJetMassless(Jet2, jet2massless);
  }
  else if(myProduction_ == TVar::JJGG || myProduction_ == TVar::JJVBF) {
    mela::computeFakeJet(jet1massless, higgs, jet2massless);
    double* jvbfpar = jvbf_Pint_par->GetY();
    auxiliaryProb = 1. / ( jvbfpar[0] * TMath::Gaus(fabs(jet2massless.Eta()), jvbfpar[1], jvbfpar[2], false) * ( 1.+jvbfpar[3]*TMath::Gaus(fabs(jet2massless.Eta()), jvbfpar[4], jvbfpar[5], false) ) );
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

//  if (myProduction_ == TVar::JH) cout << "Massless jet 2 E: " << jet2massless.T() << ", p: " << jet2massless.P() << endl;
  TLorentzVector total=jet1massless+jet2massless+higgs;
  jet1massless.Boost(-total.BoostVector().x(), -total.BoostVector().y(), 0);
  jet2massless.Boost(-total.BoostVector().x(), -total.BoostVector().y(), 0);
  higgs.Boost(-total.BoostVector().x(), -total.BoostVector().y(), 0);
  if (myProduction_ == TVar::JJGG || myProduction_ == TVar::JJVBF) ZZME->computeProdXS_JJH(jet1massless, jet2massless, higgs,
    myModel_, myProduction_,
    selfDHggcoupl,
    selfDHvvcoupl,
    selfDHwwcoupl,
    prob
    ); // Higgs + 2 jets: SBF or WBF
  else if (myProduction_ == TVar::JH) ZZME->computeProdXS_JH(
    jet1massless,
    higgs,
    myModel_,
    myProduction_,
    prob
    ); // Higgs + 1 jet; only SM is supported for now.

  if (myME_==TVar::JHUGen){
    if (myProduction_ == TVar::JJGG){
      constant = 1.8e-5;
      if (myModel_ == TVar::H0minus) constant *= 1.0017;
    }
    if (myProduction_ == TVar::JJVBF){
      if (myModel_ == TVar::H0minus) constant = 0.067;
    }
  }
  if (myME_==TVar::ANALYTICAL){
    //To be added later
  }

  prob*=constant;
}

void Mela::computeProdP(TLorentzVector Jet1, int Jet1_Id,
			TLorentzVector Jet2, int Jet2_Id,
			TLorentzVector Decay1, int Decay1_Id,
			TLorentzVector Decay2, int Decay2_Id,
			float& prob){

	double selfDHggcoupl[SIZE_HGG][2] = { { 0 } };
	double selfDHvvcoupl[SIZE_HVV_VBF][2] = { { 0 } };
	double selfDHwwcoupl[SIZE_HWW_VBF][2] = { { 0 } };

	computeProdP(Jet1, Jet1_Id,
				Jet2, Jet2_Id,
				Decay1, Decay1_Id,
				Decay2, Decay2_Id,
				selfDHggcoupl,
				selfDHvvcoupl,
				selfDHwwcoupl,
				prob);
}


void Mela::computeProdP(
			TLorentzVector V_daughter[2],
			TLorentzVector Higgs_daughter[4],
			int V_daughter_pdgid[2],
			int Higgs_daughter_pdgid[4],

			bool includeHiggsDecay,

			double selfDHvvcoupl[SIZE_HVV_VBF][2],
			float& prob){
    // Dedicated function for VH ME
    reset_PAux();

    float constant=1;
    double energy,p3sq,ratio;
    if (abs(V_daughter_pdgid[0]) <= 6){
      energy = V_daughter[0].T();
      p3sq = V_daughter[0].P();
      ratio = (p3sq > 0 ? (energy / p3sq) : 1);
      V_daughter[0].SetVect(V_daughter[0].Vect() * ratio);
    }
    if (abs(V_daughter_pdgid[1]) <= 6){
      energy = V_daughter[1].T();
      p3sq = V_daughter[1].P();
      ratio = (p3sq > 0 ? (energy / p3sq) : 1);
      V_daughter[1].SetVect(V_daughter[1].Vect() * ratio);
    }

    if ( abs(V_daughter_pdgid[0]) <= 16 && abs(V_daughter_pdgid[0]) >= 11 && abs(V_daughter_pdgid[1]) <= 16 && abs(V_daughter_pdgid[1]) >= 11){
      if (mela::forbidMassiveLeptons){
        mela::constrainedRemoveLeptonMass(V_daughter[0],V_daughter[1]);
      }
    }
    if ( abs(Higgs_daughter_pdgid[0]) <= 16 && abs(Higgs_daughter_pdgid[0]) >= 11 && abs(Higgs_daughter_pdgid[1]) <= 16 && abs(Higgs_daughter_pdgid[1]) >= 11){
      if (mela::forbidMassiveLeptons){
        mela::constrainedRemoveLeptonMass(Higgs_daughter[0],Higgs_daughter[1]);
      }
    }
    if ( abs(Higgs_daughter_pdgid[2]) <= 16 && abs(Higgs_daughter_pdgid[2]) >= 11 && abs(Higgs_daughter_pdgid[3]) <= 16 && abs(Higgs_daughter_pdgid[3]) >= 11){
      if (mela::forbidMassiveLeptons){
        mela::constrainedRemoveLeptonMass(Higgs_daughter[2],Higgs_daughter[3]);
      }
    }

    if (myProduction_ == TVar::ZH || myProduction_ == TVar::WH) ZZME->computeProdXS_VH(
		V_daughter,
		Higgs_daughter,
		V_daughter_pdgid,
		Higgs_daughter_pdgid,
		includeHiggsDecay,
		myModel_,
		myME_,
		myProduction_,
		selfDHvvcoupl,
		prob
		); // VH

	if(myME_==TVar::JHUGen){
		constant = 1; // No c-constant yet
	}
	else if(myME_==TVar::ANALYTICAL){
	//To be added later
	}

	prob*=constant;
}



void Mela::computePM4l(TLorentzVector Z1_lept1, int Z1_lept1Id,  // input 4-vectors
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

  if( abs(Z1_lept1Id)==11 &&  abs(Z1_lept2Id)==11 &&
      abs(Z2_lept1Id)==11 &&  abs(Z2_lept2Id)==11 )
    flavor = TVar::Flavor_4e;
  
  if( abs(Z1_lept1Id)==13 &&  abs(Z1_lept2Id)==13 &&
      abs(Z2_lept1Id)==13 &&  abs(Z2_lept2Id)==13 )
    flavor = TVar::Flavor_4mu;

  if( abs(Z1_lept1Id)==11 &&  abs(Z1_lept2Id)==11 &&
      abs(Z2_lept1Id)==13 &&  abs(Z2_lept2Id)==13 )
    flavor = TVar::Flavor_2e2mu;

  if( abs(Z1_lept1Id)==13 &&  abs(Z1_lept2Id)==13 &&
      abs(Z2_lept1Id)==11 &&  abs(Z2_lept2Id)==11 )
    flavor = TVar::Flavor_2e2mu;
  
  
  computePM4l(mzz,flavor,syst,prob);
}

void Mela::computePM4l(float mZZ, TVar::LeptonFlavor flavor, TVar::SuperMelaSyst syst, float& prob){
  reset_PAux();
  prob=-99;//default dummy.
  
  if(flavor == TVar::Flavor_Dummy) // only compute things if flavor determination succeded
    return;
  
  switch(flavor){
	  case 1: super->SetDecayChannel("4e")   ;break;
	  case 2: super->SetDecayChannel("4mu")  ;break;
	  case 3: super->SetDecayChannel("2e2mu");break;
	  default: std::cout << " unknown flavor: " << flavor << std::endl; exit(0);
  }


  if(syst == TVar::SMSyst_None){
    std::pair<double,double> m4lP = super->M4lProb(mZZ);
    if(myModel_ == TVar::HSMHiggs) // currently only supported signal is H(0+)
      prob = m4lP.first;
    if(myModel_ == TVar::bkgZZ) // currently only supported background is summed paramterization
      prob = m4lP.second;
  }
  else{
    //systematics for p(m4l)
    float mZZtmp=mZZ;
    float meanErr=float(super->GetSigShapeSystematic("meanCB") );
    if( syst == TVar::SMSyst_ScaleUp ){
      mZZtmp = mZZ*(1.0+meanErr);
      if(mZZtmp>180.0 || mZZtmp<100)mZZtmp=mZZ;      
      std::pair<double,double> m4lPScaleUp = super->M4lProb(mZZtmp);
      if(myModel_ == TVar::HSMHiggs)
	prob = m4lPScaleUp.first; 
      if(myModel_ == TVar::bkgZZ)
	prob = m4lPScaleUp.second;
    }

    if( syst == TVar::SMSyst_ScaleDown ){    
      mZZtmp = mZZ*(1.0-meanErr);
      if(mZZtmp>180 || mZZtmp<100) mZZtmp=mZZ;
      std::pair<double,double> m4lPScaleDown = super->M4lProb(mZZtmp);
      if(myModel_ == TVar::HSMHiggs) prob = m4lPScaleDown.first; 
      if(myModel_ == TVar::bkgZZ) prob = m4lPScaleDown.second;
    }
    
    float sigmaErr=float( super->GetSigShapeSystematic("sigmaCB") );
    float sigmaCB=float( super->GetSigShapeParameter("sigmaCB") );    
    if( syst == TVar::SMSyst_ResUp || syst ==  TVar::SMSyst_ResDown ){
      mZZtmp= myR->Gaus(mZZ,sigmaErr*sigmaCB);
      if(mZZtmp>180 || mZZtmp<100) mZZtmp=mZZ;
      std::pair<double,double> m4lPResUp = super->M4lProb(mZZtmp);
      if(myModel_ == TVar::HSMHiggs) prob = m4lPResUp.first; 
      if(myModel_ == TVar::bkgZZ) prob = m4lPResUp.second;      
    }
  }
}

void Mela::computeWeight(float mZZ, float mZ1, float mZ2, 
			 float costhetastar,
			 float costheta1, 
			 float costheta2,
			 float phi,
			 float phi1,
			 // return variables:
			 float& w
			 ){ // Lepton interference using JHUGen
  reset_PAux();

  float dXsec_HZZ_JHU,dXsec_HZZ_JHU_interf; // temporary prob
  
  // calculate dXsec for 4e/4mu
  setProcess(TVar::HSMHiggs,TVar::JHUGen,TVar::ZZGG);
  computeP(mZZ,mZ1,mZ2,
	   costhetastar,
	   costheta1,
	   costheta2,
	   phi,phi1,
	   1,dXsec_HZZ_JHU_interf);

  // calculate dXsec for 2e2mu
  setProcess(TVar::HSMHiggs,TVar::JHUGen,TVar::ZZGG);
  computeP(mZZ,mZ1,mZ2,
	   costhetastar,
	   costheta1,
	   costheta2,
	   phi,phi1,
	   3,dXsec_HZZ_JHU);
  
  w = dXsec_HZZ_JHU_interf / dXsec_HZZ_JHU;

  // protect against anomalously large weights
  if (w>5.) w=25./w;

}

void Mela::computeWeight(float mZZ, float mZ1, float mZ2, 
			 float costhetastar,
			 float costheta1, 
			 float costheta2,
			 float phi,
			 float phi1,
			 double couplingvals[SIZE_HVV_FREENORM],
			 // return variables:
			 float& w
			 ){
  reset_PAux();

  float dXsec_HZZ_JHU,dXsec_HZZ_JHU_interf; // temporary prob
  
  // calculate dXsec for 4e/4mu
  setProcess(TVar::HSMHiggs,TVar::JHUGen,TVar::ZZGG);
  computeP(mZZ,mZ1,mZ2,
	   costhetastar,
	   costheta1,
	   costheta2,
	   phi,phi1,
	   1,couplingvals,dXsec_HZZ_JHU_interf);

  // calculate dXsec for 2e2mu
  setProcess(TVar::HSMHiggs,TVar::JHUGen,TVar::ZZGG);
  computeP(mZZ,mZ1,mZ2,
	   costhetastar,
	   costheta1,
	   costheta2,
	   phi,phi1,
	   3,couplingvals,dXsec_HZZ_JHU);
  
  w = dXsec_HZZ_JHU_interf / dXsec_HZZ_JHU;

  // protect against anomalously large weights
  if (w>5.) w=25./w;

}
void Mela::setCTotalBkgGraphs(TFile* fcontainer, TGraph* tgC[]){ // Hope it has only 3 members in the array
	string tgname = "C_TotalBkgM4l_";

	float rValues[6]={1,5,10,15,20,25}; // Possible r Values

	for(int flavor=0;flavor<3;flavor++){
		float myWidth = 1;

		char crValue[20];

		int rCode = 2; // r=10
		myWidth = rValues[rCode];
		sprintf(crValue,"D_Gamma_gg_r%.0f",myWidth);

		string ctgM4L = tgname;
		string strChannel;
		if(flavor==0) strChannel = "4e"; // Check this
		else if(flavor==1) strChannel = "4mu"; // and this
		else strChannel = "2mu2e";
		ctgM4L = ctgM4L + strChannel + "_";
		ctgM4L = ctgM4L + crValue;
		tgC[flavor] = (TGraph*) fcontainer->Get(ctgM4L.c_str());
	}
}
void Mela::constructDggr(float mzz, int flavor, float bkg_VAMCFM_noscale, float ggzz_VAMCFM_noscale, float ggHZZ_prob_pure_noscale, float ggHZZ_prob_int_noscale, float& myDggr){
	float ctotal_bkg = tgtotalbkg[flavor-1]->Eval(mzz);

	float rValues[6]={1,5,10,15,20,25};
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
  reset_PAux();

	if(myME != TVar::MCFM || myType != TVar::D_gg10){
		cout << "Only support MCFM and D_gg10"<<endl;
		return;
	}

  setMelaLeptonInterference(TVar::DefaultLeptonInterf); // Override lepton interference setting

	float bkg_VAMCFM_noscale, ggzz_VAMCFM_noscale, ggHZZ_prob_pure_noscale, ggHZZ_prob_int_noscale, bkgHZZ_prob_noscale;
	setProcess(TVar::bkgZZ, myME, TVar::ZZGG);
	computeP(mZZ, mZ1, mZ2,
			costhetastar,costheta1,costheta2,phi,phi1,flavor, ggzz_VAMCFM_noscale);
	setProcess(TVar::HSMHiggs, myME, TVar::ZZGG);
	computeP(mZZ, mZ1, mZ2,
			costhetastar,costheta1,costheta2,phi,phi1,flavor, ggHZZ_prob_pure_noscale);
	setProcess(TVar::bkgZZ_SMHiggs, myME, TVar::ZZGG);
	computeP(mZZ, mZ1, mZ2,
			costhetastar,costheta1,costheta2,phi,phi1,flavor, bkgHZZ_prob_noscale);
	setProcess(TVar::bkgZZ, myME, TVar::ZZQQB);
	computeP(mZZ, mZ1, mZ2,
			costhetastar,costheta1,costheta2,phi,phi1,flavor, bkg_VAMCFM_noscale,0);
	ggHZZ_prob_int_noscale = bkgHZZ_prob_noscale - ggHZZ_prob_pure_noscale -  ggzz_VAMCFM_noscale;
	float myDggr;
	constructDggr(mZZ, flavor, bkg_VAMCFM_noscale, ggzz_VAMCFM_noscale, ggHZZ_prob_pure_noscale, ggHZZ_prob_int_noscale, myDggr);
	prob=myDggr;
}
