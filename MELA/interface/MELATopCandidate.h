#ifndef MELATOPCANDIDATE_H
#define MELATOPCANDIDATE_H

#include <ZZMatrixElement/MELA/interface/MELAParticle.h>

class MELATopCandidate : public MELAParticle{
public:
  MELATopCandidate(int id_, TLorentzVector p4_) : MELAParticle(id_, p4_), lightQuark(0), Wferm(0), Wfermbar(0) {}
  MELATopCandidate(MELAParticle* lightQuark_, MELAParticle* Wferm_, MELAParticle* Wfermbar_);
  ~MELATopCandidate(){}

  MELAParticle* setLightQuark(MELAParticle* myParticle);
  MELAParticle* setWFermion(MELAParticle* myParticle);
  MELAParticle* setWAntifermion(MELAParticle* myParticle);

  MELAParticle* getLightQuark(){ return lightQuark; }
  MELAParticle* getWFermion(){ return Wferm; }
  MELAParticle* getWAntifermion(){ return Wfermbar; }

protected:
  MELAParticle* lightQuark;
  MELAParticle* Wferm;
  MELAParticle* Wfermbar;

};


#endif
