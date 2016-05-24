#include <ZZMatrixElement/MELA/interface/MELACandidate.h>

void MELACandidate::sortDaughters(){
  if (debugVars::debugFlag) std::cout << "Starting MELACandidate::sortDaughtersInitial" << std::endl;
  sortDaughtersInitial();
  if (debugVars::debugFlag) std::cout << "Starting MELACandidate::sortDaughtersByBestZ1" << std::endl;
  sortDaughtersByBestZ1();
  if (debugVars::debugFlag) std::cout << "Starting MELACandidate::createSortedVs" << std::endl;
  createSortedVs();
}

std::vector<int> MELACandidate::getDaughterIds()const{
  std::vector<int> result;
  for (unsigned int idau=0; idau<sortedDaughters.size(); idau++){
    if (sortedDaughters.at(idau)!=0) result.push_back(sortedDaughters.at(idau)->id);
  }
  return result;
}
std::vector<int> MELACandidate::getAssociatedParticleIds()const{
  std::vector<int> result;
  for (unsigned int ip=0; ip<associatedLeptons.size(); ip++){
    if (associatedLeptons.at(ip)!=0) result.push_back(associatedLeptons.at(ip)->id);
  }
  for (unsigned int ip=0; ip<associatedPhotons.size(); ip++){
    if (associatedPhotons.at(ip)!=0) result.push_back(associatedPhotons.at(ip)->id);
  }
  for (unsigned int ip=0; ip<associatedJets.size(); ip++){
    if (associatedJets.at(ip)!=0) result.push_back(associatedJets.at(ip)->id);
  }
  return result;
}
MELAParticle* MELACandidate::getSortedDaughter(int index) const{
  if ((int)sortedDaughters.size()>index) return sortedDaughters.at(index);
  else return 0;
}
MELAParticle* MELACandidate::getSortedV(int index) const{
  if ((int)sortedVs.size()>index) return sortedVs.at(index);
  else return 0;
}
MELAParticle* MELACandidate::getSortedTop(int index) const{
  if ((int)sortedTops.size()>index) return sortedTops.at(index);
  else return 0;
}
MELAParticle* MELACandidate::getAssociatedLepton(int index)const{
  if ((int)associatedLeptons.size()>index) return associatedLeptons.at(index);
  else return 0;
}
MELAParticle* MELACandidate::getAssociatedNeutrino(int index)const{
  if ((int)associatedNeutrinos.size()>index) return associatedNeutrinos.at(index);
  else return 0;
}
MELAParticle* MELACandidate::getAssociatedPhoton(int index)const{
  if ((int)associatedPhotons.size()>index) return associatedPhotons.at(index);
  else return 0;
}
MELAParticle* MELACandidate::getAssociatedJet(int index)const{
  if ((int)associatedJets.size()>index) return associatedJets.at(index);
  else return 0;
}
void MELACandidate::sortDaughtersInitial(){
  int tmpDindex[2]={ 0 };
  MELAParticle* df[2] = { getDaughter(0), 0 };
  for (int j=1; j<getNDaughters(); j++){
    MELAParticle* dtmp = getDaughter(j);
    if (
      (dtmp->charge()+df[0]->charge()==0 && (PDGHelpers::HVVmass==PDGHelpers::Zmass || PDGHelpers::HVVmass==PDGHelpers::Zeromass))
      ||
      (std::abs(dtmp->charge()+df[0]->charge())==1 && PDGHelpers::HVVmass==PDGHelpers::Wmass)
      ){
      df[1] = dtmp;
      tmpDindex[1] = j;
      break;
    }
  }
  MELAParticle* ds[2] ={ 0 };
  int sindex=0;
  for (int j=1; j<getNDaughters(); j++){
    if (j==tmpDindex[1]) continue;
    MELAParticle* dtmp = getDaughter(j);
    ds[sindex] = dtmp;
    sindex++;
  }
  if (
    (df[0]!=0 && df[1]!=0)
    &&
    (
    // Order by ubar(0)v(1)
    (df[0]->id<df[1]->id && (PDGHelpers::HVVmass==PDGHelpers::Zmass || PDGHelpers::HVVmass==PDGHelpers::Zeromass))
    ||
    (df[0]->id<df[1]->id && df[0]->id<0 && PDGHelpers::HVVmass==PDGHelpers::Wmass)
    )
    ){
    MELAParticle* dtmp = df[0];
    df[0] = df[1];
    df[1] = dtmp;
  }
  if (
    (ds[0]!=0 && ds[1]!=0)
    &&
    (
    // Order by ubar(0)v(1)
    (ds[0]->id<ds[1]->id && (PDGHelpers::HVVmass==PDGHelpers::Zmass || PDGHelpers::HVVmass==PDGHelpers::Zeromass))
    ||
    (ds[0]->id<ds[1]->id && ds[0]->id<0 && PDGHelpers::HVVmass==PDGHelpers::Wmass)
    )
    ){
    MELAParticle* dtmp = ds[0];
    ds[0] = ds[1];
    ds[1] = dtmp;
  }
  if (df[1]==0 && df[0]!=0 && ds[0]!=0 && ds[1]!=0){
    for (int ip=0; ip<2; ip++){
      MELAParticle* dtmp = ds[ip];
      ds[ip] = df[ip];
      df[ip] = dtmp;
    }
  }
  for (int i=0; i<2; i++){
    if (df[i]!=0) sortedDaughters.push_back(df[i]);
  }
  for (int i=0; i<2; i++){
    if (ds[i]!=0) sortedDaughters.push_back(ds[i]);
  }
}
void MELACandidate::sortDaughtersByBestZ1(){
  MELAParticle* orderedDs[2][2]={ { 0 } };

  TLorentzVector pZ1(0, 0, 0, 0);
  TLorentzVector pZ2(0, 0, 0, 0);
  if (sortedDaughters.size()>2){ // WW, ZZ, ZG
    bool dauDiffType = true;
    if (debugVars::debugFlag) std::cout << "Ndaughters>2" << std::endl;

    for (int d=0; d<std::min(2, (int)sortedDaughters.size()); d++){
      if (sortedDaughters.at(d)!=0) pZ1 = pZ1 + sortedDaughters.at(d)->p4;
    }
    for (int d=std::min(2, (int)sortedDaughters.size()); d<std::min(4, (int)sortedDaughters.size()); d++){
      if (sortedDaughters.at(d)!=0) pZ2 = pZ2 + sortedDaughters.at(d)->p4;
    }

    if (debugVars::debugFlag) std::cout << "Preliminary pZ1 and pZ2 calculated!" << std::endl;

    if (sortedDaughters.size()>=4){
      if (
        (PDGHelpers::HVVmass==PDGHelpers::Zmass || PDGHelpers::HVVmass==PDGHelpers::Zeromass)
        &&
        (
        (isALepton(sortedDaughters.at(0)->id) && isALepton(sortedDaughters.at(1)->id) && isALepton(sortedDaughters.at(2)->id) && isALepton(sortedDaughters.at(3)->id))
        ||
        (isANeutrino(sortedDaughters.at(0)->id) && isANeutrino(sortedDaughters.at(1)->id) && isANeutrino(sortedDaughters.at(2)->id) && isANeutrino(sortedDaughters.at(3)->id))
        ||
        (isAPhoton(sortedDaughters.at(0)->id) && isAPhoton(sortedDaughters.at(1)->id) && isAPhoton(sortedDaughters.at(2)->id) && isAPhoton(sortedDaughters.at(3)->id))
        ||
        (isAJet(sortedDaughters.at(0)->id) && isAJet(sortedDaughters.at(1)->id) && isAJet(sortedDaughters.at(2)->id) && isAJet(sortedDaughters.at(3)->id))
        )
         ) dauDiffType=false;
    }

    if (
      (dauDiffType && (PDGHelpers::HVVmass==PDGHelpers::Zmass || PDGHelpers::HVVmass==PDGHelpers::Zeromass) && sortedDaughters.size()<4)
      ||
      (
      dauDiffType && (PDGHelpers::HVVmass==PDGHelpers::Zmass || PDGHelpers::HVVmass==PDGHelpers::Zeromass) && sortedDaughters.size()>=4 && (
      isALepton(sortedDaughters.at(0)->id) ||
      (isANeutrino(sortedDaughters.at(0)->id) && !isALepton(sortedDaughters.at(2)->id)) ||
      (isAPhoton(sortedDaughters.at(0)->id) && !isALepton(sortedDaughters.at(2)->id) && !isANeutrino(sortedDaughters.at(2)->id)) ||
      (isAJet(sortedDaughters.at(0)->id) && !isALepton(sortedDaughters.at(2)->id) && !isANeutrino(sortedDaughters.at(2)->id) && !isAPhoton(sortedDaughters.at(2)->id))
      )
      )
      ||
      (std::abs(pZ1.M() - PDGHelpers::HVVmass)<std::abs(pZ2.M() - PDGHelpers::HVVmass) && !dauDiffType && (PDGHelpers::HVVmass==PDGHelpers::Zmass || PDGHelpers::HVVmass==PDGHelpers::Zeromass)) // Z1 / Z2
      ||
      ((sortedDaughters.at(0)!=0 && sortedDaughters.at(1)!=0 && PDGHelpers::HVVmass==PDGHelpers::Wmass) && sortedDaughters.at(0)->charge()+sortedDaughters.at(1)->charge()>0) // W+ / W-
      ){
      if (debugVars::debugFlag) std::cout << "pZ1 is closer to HVVmass " << PDGHelpers::HVVmass << std::endl;
      orderedDs[0][0]=sortedDaughters.at(0);
      orderedDs[0][1]=sortedDaughters.at(1);
      orderedDs[1][0]=((int)sortedDaughters.size()>2 ? sortedDaughters.at(2) : 0);
      orderedDs[1][1]=((int)sortedDaughters.size()>3 ? sortedDaughters.at(3) : 0);
    }
    else{
      if (debugVars::debugFlag) std::cout << "pZ2 is closer to HVVmass " << PDGHelpers::HVVmass << std::endl;
      orderedDs[0][0]=((int)sortedDaughters.size()>2 ? sortedDaughters.at(2) : 0);
      orderedDs[0][1]=((int)sortedDaughters.size()>3 ? sortedDaughters.at(3) : 0);
      orderedDs[1][0]=sortedDaughters.at(0);
      orderedDs[1][1]=sortedDaughters.at(1);
      TLorentzVector ptmp = pZ1;
      pZ1 = pZ2;
      pZ2 = ptmp;
    }
  }
  else if (sortedDaughters.size()==2){ // GG, ffbar
    if (debugVars::debugFlag) std::cout << "Ndaughters==2" << std::endl;

    if (sortedDaughters.at(0)!=0) pZ1 = pZ1 + sortedDaughters.at(0)->p4;
    if (sortedDaughters.at(1)!=0) pZ2 = pZ2 + sortedDaughters.at(1)->p4;
    orderedDs[0][0]=sortedDaughters.at(0);
    orderedDs[0][1]=0;
    orderedDs[1][0]=sortedDaughters.at(1);
    orderedDs[1][1]=0;
  }
  if (
    (
    (orderedDs[0][1]!=0 && orderedDs[1][1]!=0)
    &&
    (orderedDs[0][1]->id == orderedDs[1][1]->id)
    )
    ||
    (
    (orderedDs[1][0]!=0 && orderedDs[0][0]!=0)
    &&
    (orderedDs[1][0]->id == orderedDs[0][0]->id)
    )
    ){
    if (debugVars::debugFlag) std::cout << "Checking alternative pairings." << std::endl;

    MELAParticle* orderedDps[2][2]={ { 0 } };

    TLorentzVector pZ1p(0, 0, 0, 0);
    TLorentzVector pZ2p(0, 0, 0, 0);
    for (int d=0; d<2; d++){
      if (orderedDs[d][d]!=0) pZ1p = pZ1p + orderedDs[d][d]->p4;
    }
    for (int d=0; d<2; d++){
      if (orderedDs[1-d][d]!=0) pZ2p = pZ2p + orderedDs[1-d][d]->p4;
    }

    if (std::abs(pZ1p.M() - PDGHelpers::HVVmass)<std::abs(pZ2p.M() - PDGHelpers::HVVmass)){
      orderedDps[0][0]=orderedDs[0][0];
      orderedDps[0][1]=orderedDs[1][1];
      orderedDps[1][0]=orderedDs[1][0];
      orderedDps[1][1]=orderedDs[0][1];
    }
    else{
      orderedDps[0][0]=orderedDs[1][0];
      orderedDps[0][1]=orderedDs[0][1];
      orderedDps[1][0]=orderedDs[0][0];
      orderedDps[1][1]=orderedDs[1][1];
      TLorentzVector ptmp = pZ1p;
      pZ1p = pZ2p;
      pZ2p = ptmp;
    }
    if (std::abs(pZ1p.M() - PDGHelpers::HVVmass)<std::abs(pZ1.M() - PDGHelpers::HVVmass) || (std::abs(pZ1p.M() - PDGHelpers::HVVmass)==std::abs(pZ1.M() - PDGHelpers::HVVmass) && pZ2p.Pt()>pZ2.Pt()) ){
      for (int i=0; i<2; i++){
        for (int j=0; j<2; j++) orderedDs[i][j] = orderedDps[i][j];
      }
    }
  }
  sortedDaughters.clear();
  for (int i=0; i<2; i++){
    for (int j=0; j<2; j++){
      if (orderedDs[i][j]!=0) sortedDaughters.push_back(orderedDs[i][j]);
    }
  }
  if (debugVars::debugFlag) std::cout << "Final number of daughters in sortedDaughters: " << sortedDaughters.size() << std::endl;
}
void MELACandidate::createSortedVs(){
  int VID = 23;
  if (PDGHelpers::HVVmass==PDGHelpers::Wmass) VID = 24;
  else if (PDGHelpers::HVVmass==PDGHelpers::Zeromass) VID = 21;

  TLorentzVector pZ1(0, 0, 0, 0);
  TLorentzVector pZ2(0, 0, 0, 0);
  int icutoff = (sortedDaughters.size()>2 ? 2 : (sortedDaughters.size()>0 ? 1 : 0));
  int imax = (sortedDaughters.size()>4 ? 4 : (int)sortedDaughters.size());
  int V1id=0, V2id=0;
  if (icutoff==0){
    pZ1=this->p4;
    V1id=25;
  }
  else{
    for (int d=0; d<icutoff; d++){
      if (sortedDaughters.at(d)!=0){
        pZ1 = pZ1 + sortedDaughters.at(d)->p4;
        if (icutoff==2) V1id=VID;
        else if (icutoff==1) V1id=sortedDaughters.at(d)->id;
      }
    }
    for (int d=icutoff; d<imax; d++){
      if (sortedDaughters.at(d)!=0){
        pZ2 = pZ2 + sortedDaughters.at(d)->p4;
        if ((imax-icutoff)==2) V2id=VID;
        else if ((imax-icutoff)==1) V2id=sortedDaughters.at(d)->id;
      }
    }
  }

  // If the number of Zs is less than 2, should still create empty particles
  MELAParticle* Z1 = new MELAParticle(V1id, pZ1);
  Z1->addMother(this);
  for (int d=0; d<icutoff; d++){
    if (sortedDaughters.at(d)!=0) Z1->addDaughter(sortedDaughters.at(d));
  }
  addSortedV(Z1);

  MELAParticle* Z2 = new MELAParticle(V2id, pZ2);
  Z2->addMother(this);
  for (int d=icutoff; d<imax; d++){
    if (sortedDaughters.at(d)!=0) Z2->addDaughter(sortedDaughters.at(d));
  }
  addSortedV(Z2);
}
TLorentzVector MELACandidate::getAlternativeVMomentum(int index)const{
  TLorentzVector nullFourVector(0, 0, 0, 0);
  if (sortedDaughters.size()>=3){
    TLorentzVector pZ1 = sortedDaughters.at(0)->p4;
    if (sortedDaughters.size()>3) pZ1 = pZ1 + sortedDaughters.at(3)->p4;
    TLorentzVector pZ2 = sortedDaughters.at(2)->p4+sortedDaughters.at(1)->p4;
    if (std::abs(pZ1.M() - PDGHelpers::HVVmass)>std::abs(pZ2.M() - PDGHelpers::HVVmass)){
      TLorentzVector pZtmp = pZ1;
      pZ1 = pZ2;
      pZ2 = pZtmp;
    }
    return (index==0 ? pZ1 : pZ2);
  }
  else return nullFourVector;
}
bool MELACandidate::daughtersInterfere()const{
  bool doInterfere = false;
  if (sortedDaughters.size()>3) doInterfere = (
    (sortedDaughters.at(0))->id == (sortedDaughters.at(2))->id
    &&
    (sortedDaughters.at(1))->id == (sortedDaughters.at(3))->id
    &&
    !isAnUnknownJet(sortedDaughters.at(0)) && !isAnUnknownJet(sortedDaughters.at(2))
    );
  return doInterfere;
}


bool MELACandidate::checkDaughtership(MELAParticle* myParticle)const{
  for (int dd=0; dd<getNDaughters(); dd++){
    if (myParticle==getDaughter(dd)) return true;
  }
  return false;
}

void MELACandidate::addAssociatedLeptons(MELAParticle* myParticle){
  if (!checkDaughtership(myParticle)) addByHighestPt(myParticle, associatedLeptons);
}
void MELACandidate::addAssociatedNeutrinos(MELAParticle* myParticle){
  if (!checkDaughtership(myParticle)){
    addByHighestPt(myParticle, associatedLeptons); // Neutrinos are leptons at the ZZ candidate level
    addByHighestPt(myParticle, associatedNeutrinos);
  }
}
void MELACandidate::addAssociatedPhotons(MELAParticle* myParticle){
  if (!checkDaughtership(myParticle)) addByHighestPt(myParticle, associatedPhotons);
}
void MELACandidate::addAssociatedJets(MELAParticle* myParticle){
  if (!checkDaughtership(myParticle)) addByHighestPt(myParticle, associatedJets);
}
void MELACandidate::addByHighestPt(MELAParticle* myParticle, std::vector<MELAParticle*>& particleArray){
  bool inserted = checkParticleExists(myParticle, particleArray); // Test if the particle is already in the vector
  if (!inserted){
    if (!associatedByHighestPt){ particleArray.push_back(myParticle); return; }

    for (std::vector<MELAParticle*>::iterator it = particleArray.begin(); it<particleArray.end(); it++){
      if ((*it)->pt()<myParticle->pt()){
        inserted=true;
        particleArray.insert(it, myParticle);
        break;
      }
    }
    if (!inserted) particleArray.push_back(myParticle);
  }
}
void MELACandidate::addAssociatedVs(){
  createAssociatedVs(associatedJets);
  createAssociatedVs(associatedLeptons);
}
void MELACandidate::addAssociatedTops(){ createAssociatedTops(); }
void MELACandidate::createAssociatedVs(std::vector<MELAParticle*>& particleArray){
  for (unsigned int i = 0; i<particleArray.size(); i++){
    double Qi = particleArray.at(i)->charge();
    int id_i = particleArray.at(i)->id;

    for (unsigned int j = 1; j<particleArray.size(); j++){
      if (j<=i) continue;
      double Qj = particleArray.at(j)->charge();
      int id_j = particleArray.at(j)->id;

      int bosonId=-1;
      if ((Qi+Qj)==0 && (id_i+id_j)==0) bosonId = (id_i==0 ? 0 : 23);
      else if (
        abs(Qi+Qj)==1 // W boson
        &&
        (
        ((PDGHelpers::isALepton(id_i) || PDGHelpers::isALepton(id_j)) && std::abs(id_i+id_j)==1) // Require SF lnu particle-antiparticle pairs
        ||
        (PDGHelpers::isUpTypeQuark(id_i) && PDGHelpers::isDownTypeQuark(id_j)) // Require ud- or du-type pairs, qqbar requirement is satisfied with charge.
        ||
        (PDGHelpers::isDownTypeQuark(id_i) && PDGHelpers::isUpTypeQuark(id_j))
        )
        ) bosonId=24*(Qi+Qj);

      if (bosonId!=-1){
        TLorentzVector pV = particleArray.at(i)->p4+particleArray.at(j)->p4;
        MELAParticle* boson = new MELAParticle(bosonId, pV);
        int firstdaughter = i, seconddaughter = j;
        if (
          (particleArray.at(firstdaughter)->id<particleArray.at(seconddaughter)->id && !PDGHelpers::isAWBoson(bosonId))
          ||
          //(std::abs(particleArray.at(firstdaughter)->id)<std::abs(particleArray.at(seconddaughter)->id) && PDGHelpers::isAWBoson(bosonId)) // Order by nu-l / u-d
          (particleArray.at(firstdaughter)->id<particleArray.at(seconddaughter)->id && particleArray.at(firstdaughter)->id<0 && PDGHelpers::isAWBoson(bosonId)) // Order by f-f'b
          ){
          firstdaughter = j; seconddaughter = i;
        }
        boson->addDaughter(particleArray.at(firstdaughter));
        boson->addDaughter(particleArray.at(seconddaughter));
        addSortedV(boson);
      }
    }
  }
}
void MELACandidate::createAssociatedTops(){
  for (unsigned int i = 2; i<sortedVs.size(); i++){
    double Qi = sortedVs.at(i)->charge();
    int id_i = sortedVs.at(i)->id;
    if (id_i==23 || id_i==21) continue; // Skip definite Z or gamma
    if (sortedVs.at(i)->getNDaughters()<2) continue; // Skip in case the intermediate V has no daughters
    MELAParticle* Vdau[2]={ 0 };
    for (int idau=0; idau<2; idau++) Vdau[idau] = sortedVs.at(i)->getDaughter(idau);
    if (Vdau[0]==0 || Vdau[1]==0) continue; // Pass empty particles
    if (Vdau[0]->topTag==0 || Vdau[1]->topTag==0 || Vdau[0]->topTag!=Vdau[1]->topTag) continue; // Ensure to use the same top-tagged fermions, leptons
    int Vtoptag = Vdau[0]->topTag;
    if ((Vtoptag==6 && Qi<0) || (Vtoptag==-6 && Qi>0)) continue; // Ensure that the top tag is consistent with the W+- charge
    if (Vdau[0]->id==6 || Vdau[1]->id==6 || Vdau[0]->id==-6 || Vdau[1]->id==-6) continue; // Veto Ws with tops as daughters, such tops should be standalone

    for (unsigned int j = 0; j<associatedJets.size(); j++){
      MELAParticle* theBot = associatedJets.at(j);
      double Qj = theBot->charge();
      int id_j = theBot->id;
      if (PDGHelpers::isUpTypeQuark(id_j)) continue; // Up-type quarks cannot come from top decays
      // Veto untagged ones, inconsistent charges, and those with incompatible tags with the tag of the "W"
      if (theBot->topTag==0 || (theBot->topTag==6 && Qj>0) || (theBot->topTag==-6 && Qj<0) || !(theBot->topTag==Vtoptag || theBot->topTag==1 || Vtoptag==1)) continue;
      bool isUnique = !(Vdau[0]==theBot || Vdau[1]==theBot);
      if (!isUnique) continue;

      int topId = -1;
      if ((Qi+Qj)==0 || id_i==0 || id_j==0){
         topId = 0;
         if (topTag==6) topId=6;
         else if (topTag==-6) topId=-6;
      }
      else if ((Qi+Qj)<0) topId=-6;
      else if ((Qi+Qj)>0) topId=6;


      if (topId!=-1){
        TLorentzVector pV = sortedVs.at(i)->p4+associatedJets.at(j)->p4;
        MELAParticle* theTop = new MELAParticle(topId, pV);
        theTop->addDaughter(Vdau[0]);
        theTop->addDaughter(Vdau[1]);
        theTop->addDaughter(theBot);
        int theTopTag;
        if (Vtoptag==1 || theBot->topTag==1) theTopTag=1;
        else theTopTag = Vtoptag;
        theTop->tagTop(theTopTag);
        bool success = addSortedTop(theTop);
        if (!success) delete theTop; // Avoid memory leak if adding to the array is not successful.
      }
    }
  }
}
void MELACandidate::addSortedTop(MELAParticle* myParticle){
  bool foundEquivalent=false;
  for (unsigned int itop=0; itop<sortedTops.size(); itop++){
    MELAParticle* existingTop = sortedTops.at(itop);
    if (existingTop->getNDaughters()!=myParticle->getNDaughters()) continue;
    int nidentical=0;
    for (int idau=0; idau<existingTop->getNDaughters(); idau++){
      bool dauIdentical=false;
      for (int jdau=0; jdau<myParticle->getNDaughters(); jdau++){
        if (myParticle->getDaughter(jdau)==existingTop->getDaughter(idau)){ dauIdentical=true; break; }
      }
      if (dauIdentical) nidentical++;
    }
    if (nidentical==existingTop->getNDaughters()){ foundEquivalent=true; break; } // If all daughters are identical, myParticle is just a different permutation.
  }
  if (!foundEquivalent) sortedTops.push_back(myParticle);
  return (!foundEquivalent);
}
void MELACandidate::tagAsTopDaughter(int collection, int position, int tag){
// collection=0 for jet, =1 for lepton, =2 for neutrino;  position specifies which particle in the array
// Notice that lepton and neutrino are handled separately, differing from how neutrinos are also leptons here.
// The reasoning is that the user should not have to deal with how leptons are defined.
// Neutrinos are defined as leptons here for convenience in iteration.
  if (collection<0 || collection>2) return;

  std::vector<MELAParticle*>* particleArray;
  if (collection==0) particleArray = &associatedJets;
  else particleArray = &associatedLeptons;

  int postracker=-1;
  for (int ip=0; ip<particleArray->size(); ip++){
    int idtmp = particleArray->at(ip)->id;
    if (
      PDGHelpers::isAJet(idtmp)
      ||
      (collection==1 && PDGHelpers::isALepton(idtmp))
      ||
      (collection==2 && PDGHelpers::isANeutrino(idtmp))
      ) postracker++;
    if (postracker==position){ particleArray->at(ip)->tagTop(tag); break; }
  }
}


void MELACandidate::testPreSelectedDaughters(){
  for (int i = 0; i<getNDaughters(); i++){
    if (!(daughters.at(i)->passSelection)){
      passSelection=false;
      break;
    }
  }
}


