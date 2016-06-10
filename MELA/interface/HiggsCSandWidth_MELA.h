#ifndef HIGGSCSANDWIDTH_H
#define HIGGSCSANDWIDTH_H

#include <string>
#include "TGraph.h"
#include "TSpline.h"


class HiggsCSandWidth_MELA{
public:
  HiggsCSandWidth_MELA(std::string fileLoc = "../txtFiles", std::string strAppend="YR3");
  ~HiggsCSandWidth_MELA();
  double HiggsWidth(double mH);
protected:
  std::string fileName;
  std::vector<double> mass_BR;
  std::vector<double> BR;
  double* xmhW;
  double* sigW;
  TGraph* graphW;
  TSpline3* gsW;
};

#endif

