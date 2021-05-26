//
//  ScistErrRateInf.hpp
//
//
//  Created by Yufeng Wu on 6/20/18.
//
//

#ifndef ScistErrRateInf_hpp
#define ScistErrRateInf_hpp

#include "ScistGenotype.hpp"
#include "UtilsNumerical.h"

// *************************************************************************************
// Inf error rate

class ScistErrRateInf {
public:
  ScistErrRateInf(ScistGenGenotypeMat &genos);
  void Infer();
  void SetVerbose(bool f) { fVerbose = f; }

private:
  double CalcMaxProbFor(
      double rateFN, double rateFP,
      std::set<std::pair<std::pair<int, int>, int> > &listChangedPlaces);
  double CalcMaxProbForMat(
      ScistGenGenotypeMat &genosTest,
      std::set<std::pair<std::pair<int, int>, int> > &listChangedPlaces);
  void UpdateEstimates(
      const std::set<std::pair<std::pair<int, int>, int> > &listChangedPlaces);

  ScistGenGenotypeMat &genosInput;
  double rateFNMin;
  double rateFNMax;
  double rateFPMin;
  double rateFPMax;
  double rateFNOpt;
  double rateFPOpt;
  bool fVerbose;
};

#endif /* ScistErrRateInf_hpp */
