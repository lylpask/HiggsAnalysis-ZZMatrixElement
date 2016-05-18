#ifndef MELACANDIDATE_H
#define MELACANDIDATE_H

#include <ZZMatrixElement/MELA/interface/MELAParticle.h>

class MELACandidate : public MELAParticle{
public:
  MELACandidate(int id_, TLorentzVector p4_) : MELAParticle(id_, p4_) {};

  ~MELACandidate(){ for (unsigned int i=0; i<sortedVs.size(); i++) delete sortedVs.at(i); sortedVs.clear(); sortedDaughters.clear(); associatedJets.clear(); associatedLeptons.clear(); associatedNeutrinos.clear(); associatedPhotons.clear(); };


  // Member functions

  MELAParticle* getSortedDaughter(int index)const;
  MELAParticle* getSortedV(int index)const;
  MELAParticle* getAssociatedLepton(int index)const;
  MELAParticle* getAssociatedNeutrino(int index)const;
  MELAParticle* getAssociatedPhoton(int index)const;
  MELAParticle* getAssociatedJet(int index)const;
  TLorentzVector getAlternativeVMomentum(int index)const;

  int getNAssociatedLeptons()const{ return associatedLeptons.size(); };
  int getNAssociatedNeutrinos()const{ return associatedNeutrinos.size(); };
  int getNAssociatedPhotons()const{ return associatedPhotons.size(); };
  int getNAssociatedJets()const{ return associatedJets.size(); };
  int getNSortedVs()const{ return sortedVs.size(); };

  void addAssociatedLeptons(MELAParticle* myMELAParticle);
  void addAssociatedNeutrinos(MELAParticle* myMELAParticle);
  void addAssociatedPhotons(MELAParticle* myMELAParticle);
  void addAssociatedJets(MELAParticle* myMELAParticle);
  void addSortedV(MELAParticle* myMELAParticle){ sortedVs.push_back(myMELAParticle); };
  void addAssociatedVs();

  void sortDaughters();
  void testPreSelectedDaughters();

private:
  std::vector<MELAParticle*> associatedLeptons;
  std::vector<MELAParticle*> associatedNeutrinos;
  std::vector<MELAParticle*> associatedPhotons;
  std::vector<MELAParticle*> associatedJets;

  std::vector<MELAParticle*> sortedDaughters;
  std::vector<MELAParticle*> sortedVs;


  void sortDaughtersInitial();
  void sortDaughtersByBestZ1();
  void createSortedVs();
  bool checkDaughtership(MELAParticle* myMELAParticle)const;
  void createAssociatedVs(std::vector<MELAParticle*>& particleArray);
  void addByHighestPt(MELAParticle* myMELAParticle, std::vector<MELAParticle*>& particleArray);
};



#endif
