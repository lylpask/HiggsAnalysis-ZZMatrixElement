#include "ZZMatrixElement/PythonWrapper/interface/MEMCalculatorsWrapper.h"

#include "ZZMatrixElement/MEMCalculators/interface/MEMCalculators.h"
#include "ZZMatrixElement/MELA/src/computeAngles.h"


MEMCalculatorsWrapper::MEMCalculatorsWrapper(double collisionEnergy, double sKD_mass) {
    mem_ = new MEMs(collisionEnergy,sKD_mass);
}

MEMCalculatorsWrapper::~MEMCalculatorsWrapper() {
    if(mem_ !=0) delete mem_;
}

MEMCalculatorsWrapper::Angles 
MEMCalculatorsWrapper::computeAngles(TLorentzVector Z1_lept1, int Z1_lept1Id,
		   TLorentzVector Z1_lept2, int Z1_lept2Id,
		   TLorentzVector Z2_lept1, int Z2_lept1Id,
		   TLorentzVector Z2_lept2, int Z2_lept2Id) {
    Angles ret;
    mela::computeAngles(Z1_lept1,Z1_lept1Id,
                        Z1_lept2,Z1_lept2Id,
                        Z2_lept1,Z2_lept1Id,
                        Z2_lept2,Z2_lept2Id,
                    ret.costhetastar,ret.costheta1,ret.costheta2,ret.phi,ret.phistar1);
    return ret;
  }


void  
MEMCalculatorsWrapper::computeAll(TLorentzVector Z1_lept1, int Z1_lept1Id,
		   TLorentzVector Z1_lept2, int Z1_lept2Id,
		   TLorentzVector Z2_lept1, int Z2_lept1Id,
		   TLorentzVector Z2_lept2, int Z2_lept2Id) {

    std::vector<TLorentzVector> ps;
    ps.push_back(Z1_lept1);
    ps.push_back(Z1_lept2);
    ps.push_back(Z2_lept1);
    ps.push_back(Z2_lept2);

    if ( Z2_lept1Id == Z2_lept2Id)
      Z2_lept2Id=-Z2_lept2Id;


    std::vector<int> id;
    id.push_back(Z1_lept1Id);
    id.push_back(Z1_lept2Id);
    id.push_back(Z2_lept1Id);
    id.push_back(Z2_lept2Id);

    mem_->computeMEs(ps,id);

    //Now the SuperKD part
    pm4l_sig_=0.0;
    pm4l_bkg_=0.0;
    mem_->computePm4l(ps,id,kNone,pm4l_sig_,pm4l_bkg_);
}


float
MEMCalculatorsWrapper:: getKD() {
    using namespace MEMNames;
    double KD,ME_ggHiggs,ME_qqZZ;
    mem_->computeKD(kSMHiggs, kJHUGen, kqqZZ, kMCFM, &MEMs::probRatio, KD, ME_ggHiggs, ME_qqZZ);
    return KD;
}

float
MEMCalculatorsWrapper:: getSuperKD() {
    using namespace MEMNames;
    double KD,ME_ggHiggs,ME_qqZZ;
    mem_->computeKD(kSMHiggs, kJHUGen, kqqZZ, kMCFM, &MEMs::probRatio, KD, ME_ggHiggs, ME_qqZZ);
    return pm4l_sig_/(pm4l_sig_+pm4l_bkg_*(1./KD-1));
}

float
MEMCalculatorsWrapper:: getGG0KD() {
    using namespace MEMNames;
    double KD,ME_ggHiggs,ME_gg0Minus;
    mem_->computeKD(kSMHiggs, kJHUGen, k0minus, kJHUGen, &MEMs::probRatio, KD, ME_ggHiggs, ME_gg0Minus);
    return KD;
}

float
MEMCalculatorsWrapper:: getGG0HKD() {
    using namespace MEMNames;
    double KD,ME_ggHiggs,ME_gg0hPlus;
    mem_->computeKD(kSMHiggs, kJHUGen, k0hplus, kJHUGen, &MEMs::probRatio, KD, ME_ggHiggs, ME_gg0hPlus);
    return KD;
}

float
MEMCalculatorsWrapper:: getQQ1MinusKD() {
    using namespace MEMNames;
    double KD,ME_ggHiggs,ME_qq1Minus;
    mem_->computeKD(kSMHiggs, kJHUGen, k1minus, kJHUGen, &MEMs::probRatio, KD, ME_ggHiggs, ME_qq1Minus);
    return KD;
}

float
MEMCalculatorsWrapper:: getQQ1PlusKD() {
    using namespace MEMNames;
    double KD,ME_ggHiggs,ME_qq2Plus;
    mem_->computeKD(kSMHiggs, kJHUGen, k1plus, kJHUGen, &MEMs::probRatio, KD, ME_ggHiggs, ME_qq2Plus);
    return KD;
}

float
MEMCalculatorsWrapper:: getGG2PlusKD() {
    using namespace MEMNames;
    double KD,ME_ggHiggs,ME_gg2Plus;
    mem_->computeKD(kSMHiggs, kJHUGen, k2mplus_gg, kJHUGen, &MEMs::probRatio, KD, ME_ggHiggs, ME_gg2Plus);
    return KD;
}

float
MEMCalculatorsWrapper:: getQQ2PlusKD() {
    using namespace MEMNames;
    double KD,ME_ggHiggs,ME_qq2Plus;
    mem_->computeKD(kSMHiggs, kJHUGen, k2mplus_qqbar, kJHUGen, &MEMs::probRatio, KD, ME_ggHiggs, ME_qq2Plus);
    return KD;
}



float
MEMCalculatorsWrapper:: getInterferenceWeight() {
    return mem_->getMELAWeight();
}




