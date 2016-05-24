#ifndef MELACANDIDATE_H
#define MELACANDIDATE_H

#include <ZZMatrixElement/MELA/interface/MELAParticle.h>

class MELACandidate : public MELAParticle{
public:
  MELACandidate(int id_, TLorentzVector p4_, associatedByHighestPt_=false) : MELAParticle(id_, p4_),associatedByHighestPt(associatedByHighestPt_) {};

  ~MELACandidate(){ for (unsigned int i=0; i<sortedTops.size(); i++) delete sortedTops.at(i); sortedTops.clear(); for (unsigned int i=0; i<sortedVs.size(); i++) delete sortedVs.at(i); sortedVs.clear(); sortedDaughters.clear(); associatedJets.clear(); associatedLeptons.clear(); associatedNeutrinos.clear(); associatedPhotons.clear(); };


  // Member functions

  MELAParticle* getSortedDaughter(int index)const;
  MELAParticle* getSortedV(int index)const;
  MELAParticle* getSortedTop(int index)const;
  MELAParticle* getAssociatedLepton(int index)const;
  MELAParticle* getAssociatedNeutrino(int index)const;
  MELAParticle* getAssociatedPhoton(int index)const;
  MELAParticle* getAssociatedJet(int index)const;
  TLorentzVector getAlternativeVMomentum(int index)const;

  virtual std::vector<int> getDaughterIds()const;
  std::vector<int> getAssociatedParticleIds()const;

  int getNAssociatedLeptons()const{ return associatedLeptons.size(); }
  int getNAssociatedNeutrinos()const{ return associatedNeutrinos.size(); }
  int getNAssociatedPhotons()const{ return associatedPhotons.size(); }
  int getNAssociatedJets()const{ return associatedJets.size(); }
  int getNSortedVs()const{ return sortedVs.size(); }
  int getNSortedTops()const{ return sortedTops.size(); }

  void addAssociatedLeptons(MELAParticle* myParticle);
  void addAssociatedNeutrinos(MELAParticle* myParticle);
  void addAssociatedPhotons(MELAParticle* myParticle);
  void addAssociatedJets(MELAParticle* myParticle);
  void addSortedV(MELAParticle* myParticle){ sortedVs.push_back(myParticle); }
  void addAssociatedVs();
  void addAssociatedTops();

  void sortDaughters();
  void testPreSelectedDaughters();

  bool daughtersInterfere()const;

  void setAddAssociatedByHighestPt(bool associatedByHighestPt_){ associatedByHighestPt=associatedByHighestPt_; }
  void tagAsTopDaughter(int collection, int position, int tag); // collection=0 for jet, 1 for lepton, 2 for neutrino;  position specifies which particle in the array

private:
  bool associatedByHighestPt;

  std::vector<MELAParticle*> associatedLeptons;
  std::vector<MELAParticle*> associatedNeutrinos;
  std::vector<MELAParticle*> associatedPhotons;
  std::vector<MELAParticle*> associatedJets;

  std::vector<MELAParticle*> sortedDaughters;
  std::vector<MELAParticle*> sortedVs;
  std::vector<MELAParticle*> sortedTops;


  void sortDaughtersInitial();
  void sortDaughtersByBestZ1();
  void createSortedVs();
  bool checkDaughtership(MELAParticle* myParticle)const;
  void createAssociatedVs(std::vector<MELAParticle*>& particleArray);
  void createAssociatedTops();
  void addByHighestPt(MELAParticle* myParticle, std::vector<MELAParticle*>& particleArray);
  bool addSortedTop(MELAParticle* myParticle);
};



#endif
