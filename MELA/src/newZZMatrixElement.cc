#include "ZZMatrixElement/MELA/interface/newZZMatrixElement.h"
#include "TLorentzRotation.h"

newZZMatrixElement::newZZMatrixElement(const char* pathtoHiggsCSandWidth,
				       double ebeam,
				       TVar::VerbosityLevel verbosity):

  verb(verbosity),
  Xcal2(pathtoHiggsCSandWidth,ebeam),
  hzz4l_event(),
  vh_event(),
  EBEAM(ebeam)
				       {
// Set default parameters explicitly
  mHiggs = 125.6;
  wHiggs = -1;
  myLeptonInterference = TVar::DefaultLeptonInterf;

  //std::cout << "TEST" << std::endl;
}

std::vector<TLorentzVector> newZZMatrixElement::Calculate4Momentum(double Mx,double M1,double M2,double theta,double theta1,double theta2,double Phi1,double Phi){
    double phi1,phi2;
    phi1=TMath::Pi()-Phi1;
    phi2=Phi1+Phi;

    double gamma1,gamma2,beta1,beta2;

    gamma1=(Mx*Mx+M1*M1-M2*M2)/(2*Mx*M1);
    gamma2=(Mx*Mx-M1*M1+M2*M2)/(2*Mx*M2);
    beta1=sqrt(1-1/(gamma1*gamma1));
    beta2=sqrt(1-1/(gamma2*gamma2));

    //gluon 4 vectors
    TLorentzVector p1CM(0,0,Mx/2,Mx/2);
    TLorentzVector p2CM(0,0,-Mx/2,Mx/2);

    //vector boson 4 vectors
    TLorentzVector kZ1(gamma1*M1*sin(theta)*beta1,0, gamma1*M1*cos(theta)*beta1,gamma1*M1*1);
    TLorentzVector kZ2(-gamma2*M2*sin(theta)*beta2,0, -gamma2*M2*cos(theta)*beta2,gamma2*M2*1);

    //Rotation and Boost matrices. Note gamma1*beta1*M1=gamma2*beta2*M2.

    TLorentzRotation Z1ToZ,Z2ToZ;

    Z1ToZ.Boost(0,0,beta1);
    Z2ToZ.Boost(0,0,beta2);
    Z1ToZ.RotateY(theta);
    Z2ToZ.RotateY(TMath::Pi()+theta);


    //fermion 4 vectors in vector boson rest frame

    TLorentzVector p3Z1((M1/2)*sin(theta1)*cos(phi1),(M1/2)*sin(theta1)*sin(phi1),(M1/2)*cos(theta1),(M1/2)*1);
    TLorentzVector p4Z1(-(M1/2)*sin(theta1)*cos(phi1),-(M1/2)*sin(theta1)*sin(phi1),-(M1/2)*cos(theta1),(M1/2)*1);
    TLorentzVector p5Z2((M2/2)*sin(theta2)*cos(phi2),(M2/2)*sin(theta2)*sin(phi2),(M2/2)*cos(theta2),(M2/2)*1);
    TLorentzVector p6Z2(-(M2/2)*sin(theta2)*cos(phi2),-(M2/2)*sin(theta2)*sin(phi2),-(M2/2)*cos(theta2),(M2/2)*1);

    // fermions 4 vectors in CM frame

    TLorentzVector p3CM,p4CM,p5CM,p6CM;

    p3CM=Z1ToZ*p3Z1;
    p4CM=Z1ToZ*p4Z1;
    p5CM=Z2ToZ*p5Z2;
    p6CM=Z2ToZ*p6Z2;

    vector<TLorentzVector> p;

    p.push_back(p3CM);
    p.push_back(p4CM);
    p.push_back(p5CM);
    p.push_back(p6CM);

    return p;
}

void newZZMatrixElement::set_mHiggs(float myPoleMass){
	mHiggs = myPoleMass;
}
void newZZMatrixElement::set_wHiggs(float myPoleWidth){
	wHiggs = myPoleWidth;
}
void newZZMatrixElement::set_LeptonInterference(TVar::LeptonInterference myLepInterf){
	myLeptonInterference = myLepInterf;
}
void newZZMatrixElement::reset_MCFM_EWKParameters(double ext_Gf, double ext_aemmz, double ext_mW, double ext_mZ, double ext_xW, int ext_ewscheme){
  Xcal2.ResetMCFM_EWKParameters(ext_Gf, ext_aemmz, ext_mW, ext_mZ, ext_xW, ext_ewscheme);
}


void newZZMatrixElement::computeXS(float mZZ, float mZ1, float mZ2,
				 float costhetastar,
				 float costheta1,
				 float costheta2,
				 float phi,
				 float phistar1,
				 int flavor,
				 TVar::Process myModel,
				 TVar::MatrixElement myME,
				 TVar::Production myProduction,
				 double couplingvals[SIZE_HVV_FREENORM],
				 double selfDHvvcoupl[SIZE_HVV][2],
				 double selfDZqqcoupl[SIZE_ZQQ][2],
				 double selfDZvvcoupl[SIZE_ZVV][2],
				 double selfDGqqcoupl[SIZE_GQQ][2],
				 double selfDGggcoupl[SIZE_GGG][2],
				 double selfDGvvcoupl[SIZE_GVV][2],
				 float &mevalue
			){
// Higgs + 2 jets: SBF or WBF

  std::vector<TLorentzVector> p;
  p=Calculate4Momentum(mZZ,mZ1,mZ2,acos(costhetastar),acos(costheta1),acos(costheta2),phistar1,phi);
  TLorentzVector Z1_lept1 = p[0];
  TLorentzVector Z1_lept2 = p[1];
  TLorentzVector Z2_lept1 = p[2];
  TLorentzVector Z2_lept2 = p[3];
  
  int Z1_lept1Id, Z1_lept2Id, Z2_lept1Id, Z2_lept2Id; 
  assert(flavor==1 || flavor==2 || flavor==3 );
  if ( flavor == 1 ) {
    Z1_lept1Id = 11;
    Z1_lept2Id = -11;
    Z2_lept1Id = 11;
    Z2_lept2Id = -11;
  }
  if ( flavor == 2 ) {
    Z1_lept1Id = 13;
    Z1_lept2Id = -13;
    Z2_lept1Id = 13;
    Z2_lept2Id = -13;
  }
  if ( flavor == 3 ) {
    Z1_lept1Id = 11;
    Z1_lept2Id = -11;
    Z2_lept1Id = 13;
    Z2_lept2Id = -13;
  }
  
  hzz4l_event.p[0].SetXYZM(Z1_lept1.Px(), Z1_lept1.Py(), Z1_lept1.Pz(), 0.);
  hzz4l_event.p[1].SetXYZM(Z1_lept2.Px(), Z1_lept2.Py(), Z1_lept2.Pz(), 0.);
  hzz4l_event.p[2].SetXYZM(Z2_lept1.Px(), Z2_lept1.Py(), Z2_lept1.Pz(), 0.);
  hzz4l_event.p[3].SetXYZM(Z2_lept2.Px(), Z2_lept2.Py(), Z2_lept2.Pz(), 0.);
  
  hzz4l_event.PdgCode[0] =  Z1_lept1Id;
  hzz4l_event.PdgCode[1] =  Z1_lept2Id;
  hzz4l_event.PdgCode[2] =  Z2_lept1Id;
  hzz4l_event.PdgCode[3] =  Z2_lept2Id;
  
  double z1mass = (hzz4l_event.p[0]+hzz4l_event.p[1]).M();
  double z2mass = (hzz4l_event.p[2]+hzz4l_event.p[3]).M();
  double zzmass = (hzz4l_event.p[0]+hzz4l_event.p[1]+hzz4l_event.p[2]+hzz4l_event.p[3]).M();
  
  if (verb >= TVar::DEBUG) {
    std::cout << "four lepton p4 for ME calculations: ===========================" <<endl;
    std::cout << "PDG code \n";
    std::cout << "hzz4l_event.PdgCode[0] = " << hzz4l_event.PdgCode[0] << "\n";
    std::cout << "hzz4l_event.PdgCode[1] = " << hzz4l_event.PdgCode[1] << "\n";
    std::cout << "hzz4l_event.PdgCode[2] = " << hzz4l_event.PdgCode[2] << "\n";
    std::cout << "hzz4l_event.PdgCode[3] = " << hzz4l_event.PdgCode[3] << "\n";
    printf("lep1 p3 = (%4.4f, %4.4f, %4.4f)  lep2 p3 = (%4.4f, %4.4f, %4.4f)\n",
	   Z1_lept1.Px(), Z1_lept1.Py(), Z1_lept1.Pz(), Z1_lept2.Px(), Z1_lept2.Py(), Z1_lept2.Pz());
    printf("lep3 p3 = (%4.4f, %4.4f, %4.4f)  lep4 p3 = (%4.4f, %4.4f, %4.4f)\n",
	   Z2_lept1.Px(), Z2_lept1.Py(), Z2_lept1.Pz(), Z2_lept2.Px(), Z2_lept2.Py(), Z2_lept2.Pz());
    std::cout << "ZZ system (pX, pY, pZ, E, mass) = ( "
	      << (hzz4l_event.p[0]+hzz4l_event.p[1]+hzz4l_event.p[2]+hzz4l_event.p[3]).Px() << ", "
	      << (hzz4l_event.p[0]+hzz4l_event.p[1]+hzz4l_event.p[2]+hzz4l_event.p[3]).Py() << ", "
	      << (hzz4l_event.p[0]+hzz4l_event.p[1]+hzz4l_event.p[2]+hzz4l_event.p[3]).Pz() << ", "
	      << (hzz4l_event.p[0]+hzz4l_event.p[1]+hzz4l_event.p[2]+hzz4l_event.p[3]).Energy()  << ", "
	      << zzmass << ")\n";
    std::cout << "Z1 mass = " << z1mass << "\tz2mass = " << z2mass << "\n";
    std::cout << "=========================================================\n";
  }
  
  // ==== Begin the differential cross-section calculation
  if(myME==TVar::MCFM || myModel==TVar::bkgZZ_SMHiggs) Xcal2.SetHiggsMass(mHiggs,wHiggs);
  else Xcal2.SetHiggsMass(zzmass,wHiggs);
  Xcal2.SetMatrixElement(myME);
  Xcal2.SetProduction(myProduction);
  Xcal2.SetProcess(myModel);
  Xcal2.SetLeptonInterf(myLeptonInterference);
  mevalue = Xcal2.XsecCalc(myModel,myProduction,hzz4l_event,verb,couplingvals,selfDHvvcoupl,selfDZqqcoupl,selfDZvvcoupl,selfDGqqcoupl,selfDGggcoupl,selfDGvvcoupl);
  if (wHiggs>=0){
	  set_wHiggs(-1); // Protection against forgetfulness; custom width has to be set per-event
	  Xcal2.SetHiggsMass(mHiggs,-1);
  }
  if( myLeptonInterference != TVar::DefaultLeptonInterf ){
	  set_LeptonInterference(TVar::DefaultLeptonInterf); // Return back to default lepton interference settings after each calculation
	  Xcal2.SetLeptonInterf(TVar::DefaultLeptonInterf);
  }
  return;
}


void newZZMatrixElement::computeProdXS_JJH(TLorentzVector jet1,
				       TLorentzVector jet2,
				       TLorentzVector higgs,
				       TVar::Process myModel,
				       TVar::Production myProduction,
					   double selfDHggcoupl[SIZE_HGG][2],
					   double selfDHvvcoupl[SIZE_HVV_VBF][2],
					   double selfDHwwcoupl[SIZE_HWW_VBF][2],
				       float &mevalue){

  TLorentzVector p4[3];
  p4[0]=jet1;
  p4[1]=jet2;
  p4[2]=higgs;
  mevalue  = Xcal2.XsecCalcXJJ(myModel,myProduction, p4,
	  verb,
	  selfDHggcoupl,
	  selfDHvvcoupl,
	  selfDHwwcoupl
	  );
  
  return;

}

void newZZMatrixElement::computeProdXS_JH(TLorentzVector singleJet,
				       TLorentzVector higgs,
				       TVar::Process myModel,
				       TVar::Production myProduction,
				       float &mevalue){

// Higgs + 1 jet: Only SM is supported for now.

  TLorentzVector p4[2];
  p4[0]=higgs;
  p4[1]=singleJet;
  mevalue  = Xcal2.XsecCalcXJ(myModel,myProduction, p4,
	  verb
	  );
  
  return;
}

void newZZMatrixElement::computeProdXS_VH(
					   TLorentzVector V_daughter[2],
				       TLorentzVector Higgs_daughter[4],
					   int V_daughter_pdgid[2],
				       int Higgs_daughter_pdgid[4],
					   bool includeHiggsDecay,

					   TVar::Process myModel,
					   TVar::MatrixElement myME,
					   TVar::Production myProduction,

					   double selfDHvvcoupl[SIZE_HVV_VBF][2],
				       float &mevalue){

// Dedicated VH function

  Xcal2.SetMatrixElement(myME);
  Xcal2.SetProduction(myProduction);
  Xcal2.SetProcess(myModel);

  vh_event.p[1] = V_daughter[0];
  vh_event.p[2] = V_daughter[1];
  vh_event.PdgCode[1] = V_daughter_pdgid[0];
  vh_event.PdgCode[2] = V_daughter_pdgid[1];

  TLorentzVector nullVector(0,0,0,0);
  TLorentzVector higgs(0,0,0,0);
  for (int hd = 0; hd < 4; hd++){
	  higgs = higgs + Higgs_daughter[hd];
	  if (includeHiggsDecay){
		  vh_event.pHdecay[hd] = Higgs_daughter[hd];
		  vh_event.PdgCode_Hdecay[hd] = Higgs_daughter_pdgid[hd];
	  }
	  else{
		  vh_event.pHdecay[hd] = nullVector;
		  vh_event.PdgCode_Hdecay[hd] = 0;
	  }
  }
  vh_event.p[0] = higgs;
  vh_event.PdgCode[0] = 25;

//  Xcal2.SetHiggsMass(higgs.M(),wHiggs);
  Xcal2.SetHiggsMass(mHiggs,wHiggs);
  mevalue  = Xcal2.XsecCalc_VX(myModel,myProduction,
	  vh_event,
	  verb,
	  selfDHvvcoupl
	  );
  if (wHiggs>=0){
	  set_wHiggs(-1); // Protection against forgetfulness; custom width has to be set per-event
	  Xcal2.SetHiggsMass(mHiggs,-1);
  }
  
  return;
}
