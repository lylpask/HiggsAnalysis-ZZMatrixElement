#ifndef MELAINPUTOUTPUT_H
#define MELAINPUTOUTPUT_H

#include "TMCFM.hh"
#include "MELACandidate.h"
#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

class MelaIO{
private:

  double partonWeight[2][nmsq];
  double MEsq[nmsq][nmsq];
  double weightedMEsq[nmsq][nmsq];
  double sumME;

public:

  MELACandidate* melaCand; // Persistent container of the four-vectors, not owned by MelaIO

  void reset(){
    sumME=0;
    for (int ix=0; ix<nmsq; ix++){
      for (int pp=0; pp<2; pp++) partonWeight[pp][ix]=0;
      for (int iy=0; iy<nmsq; iy++){
        MEsq[ix][iy]=0;
        weightedMEsq[ix][iy]=0;
      }
    }
  }
  MelaIO* getRef(){ return this; }
  void setPartonWeights(
    double partonOneWeight_[nmsq],
    double partonTwoWeight_[nmsq]
    ){
    for (int ix=0; ix<nmsq; ix++){
      partonWeight[0][ix]=partonOneWeight_[ix];
      partonWeight[1][ix]=partonTwoWeight_[ix];
    }
  }
  void setMEArray(double MEsq_[nmsq][nmsq], bool transpose=false){
    for (int ix=0; ix<nmsq; ix++){
      for (int iy=0; iy<nmsq; iy++){
        int jx=(transpose ? iy : ix);
        int jy=(transpose ? ix : iy);
        MEsq[ix][iy]=MEsq_[jx][jy];
      }
    }
  }
  void addMEArray(double MEsq_[nmsq][nmsq], double factor=1., bool transpose=false){
    for (int ix=0; ix<nmsq; ix++){
      for (int iy=0; iy<nmsq; iy++){
        int jx=(transpose ? iy : ix);
        int jy=(transpose ? ix : iy);
        MEsq[ix][iy]+=MEsq_[jx][jy]*factor;
      }
    }
    if (sumME!=0) computeWeightedMEArray();
  }
  void addMERecord(MelaIO* rcd, double factor=1., bool overwrite=false){
    double MEsq_[nmsq][nmsq]={ { 0 } };
    double partonOneWeight_[nmsq]={ 0 };
    double partonTwoWeight_[nmsq]={ 0 };
    rcd->getUnweightedMEArray(MEsq_);
    rcd->getPartonWeights(partonOneWeight_, partonTwoWeight_);

    if (overwrite) reset();
    setPartonWeights(partonOneWeight_, partonTwoWeight_);
    addMEArray(MEsq_, factor);
  }

  void computeWeightedMEArray(){
    sumME=0;
    for (int ix=0; ix<nmsq; ix++){
      for (int iy=0; iy<nmsq; iy++){
        weightedMEsq[ix][iy]=partonWeight[0][ix]*MEsq[ix][iy]*partonWeight[1][iy];
        sumME += weightedMEsq[ix][iy];
      }
    }
  }

  MelaIO(){ melaCand=0; reset(); }

  double getSumME(){ return sumME; }
  void getWeightedMEArray(double MEsq_[nmsq][nmsq]){
    for (int ix=0; ix<nmsq; ix++){
      for (int iy=0; iy<nmsq; iy++) MEsq_[ix][iy] = weightedMEsq[ix][iy];
    }
  }
  void getUnweightedMEArray(double MEsq_[nmsq][nmsq]){
    for (int ix=0; ix<nmsq; ix++){
      for (int iy=0; iy<nmsq; iy++) MEsq_[ix][iy] = MEsq[ix][iy];
    }
  }
  void getPartonWeights(
    double partonOneWeight_[nmsq],
    double partonTwoWeight_[nmsq]
    ){
    for (int ix=0; ix<nmsq; ix++){
      partonOneWeight_[ix] = partonWeight[0][ix];
      partonTwoWeight_[ix] = partonWeight[1][ix];
    }
  }

  ClassDef(MelaIO, 0)
};

ClassImp(MelaIO)

#endif
