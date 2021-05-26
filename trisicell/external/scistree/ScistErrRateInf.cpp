//
//  ScistErrRateInf.cpp
//
//
//  Created by Yufeng Wu on 6/20/18.
//
//

#include "ScistErrRateInf.hpp"
#include "ScistPerfPhyImp.hpp"
#include "UtilsNumerical.h"

// *************************************************************************************
// Inf error rate

const double DEF_RATE_FN_MIN = 0.000001;
const double DEF_RATE_FN_MAX = 0.5;
const double DEF_RATE_FP_MIN = 0.0000001;
const double DEF_RATE_FP_MAX = 0.05;

ScistErrRateInf ::ScistErrRateInf(ScistGenGenotypeMat &genos)
    : genosInput(genos), rateFNMin(DEF_RATE_FN_MIN), rateFNMax(DEF_RATE_FN_MAX),
      rateFPMin(DEF_RATE_FP_MIN), rateFPMax(DEF_RATE_FP_MAX), fVerbose(false) {
  //
  rateFNOpt = 0.5 * (rateFNMin + rateFNMax);
  rateFPOpt = 0.5 * (rateFPMin + rateFPMax);
}

void ScistErrRateInf ::Infer() {
  // EM algorithm.
  const double THRES_LARGER_RATIO = 1.05;
  double likeliMaxAll = -1.0 * HAP_MAX_INT;
  while (true) {
    // now search for rateFP then we are done
    std::set<std::pair<std::pair<int, int>, int> > listChangedPlaces;
    double likeliMax2 =
        CalcMaxProbFor(this->rateFNOpt, this->rateFPOpt, listChangedPlaces);

    if (fVerbose) {
      cout << "Current likelihood for optimizing false positive rate is "
           << likeliMax2 << ", FN estimate: " << this->rateFNOpt
           << ", FP estimate: " << this->rateFPOpt << endl;
    }
    if (NumericalAlgoUtils::IsLikeliSignificantlyLargeThresNum(
            likeliMax2, likeliMaxAll, 1, THRES_LARGER_RATIO) == false) {
      break;
    }
    likeliMaxAll = likeliMax2;
    UpdateEstimates(listChangedPlaces);
  }

  cout << "Optimal false negative rate is " << this->rateFNOpt
       << ", and optimal false positive rate is " << this->rateFPOpt << endl;
}

double ScistErrRateInf ::CalcMaxProbFor(
    double rateFN, double rateFP,
    std::set<std::pair<std::pair<int, int>, int> > &listChangedPlaces) {
  // cout << "rateFN: " << rateFN << ", rateFP: " << rateFP << endl;
  //
  ScistGenGenotypeMat *pGenosMatTest = genosInput.Copy();

  // setup prob based on the rate
  for (int s = 0; s < genosInput.GetNumSites(); ++s) {
    for (int c = 0; c < genosInput.GetNumHaps(); ++c) {
      int allele = genosInput.GetGenotypeAt(c, s);
      double prob0 = 1.0 - rateFN;
      if (allele == 1) {
        prob0 = rateFP;
      }
      // cout << "Setting cell " << c << ", site " << s << ", prob0: " << prob0
      // << endl;

      pGenosMatTest->SetGenotypeProbAt(c, s, prob0);
    }
  }
  // cout << "Genotype matrix to test: " << endl;
  // pGenosMatTest->Dump();

  double probMax = CalcMaxProbForMat(*pGenosMatTest, listChangedPlaces);
  // cout << "For rateFN: " << rateFN << ", rateFP: " << rateFP << "
  // CalcMaxProbFor: " << probMax << endl;

  delete pGenosMatTest;

  return probMax;
}

double ScistErrRateInf ::CalcMaxProbForMat(
    ScistGenGenotypeMat &genosTest,
    std::set<std::pair<std::pair<int, int>, int> > &listChangedPlaces) {
  //
  ScistPerfPhyMLE phInf1(genosTest);
  phInf1.SetVerbose(false);
  phInf1.SetOutput(false);
  double res = phInf1.Infer(&listChangedPlaces);
  // cout << "In CalcMaxProbForMat: prob=" << res << ", matrix: \n";
  // genosTest.Dump();
  return res;
}

void ScistErrRateInf ::UpdateEstimates(
    const std::set<std::pair<std::pair<int, int>, int> > &listChangedPlaces) {
  //
  int num0to1 = 0, num1to0 = 0;
  for (set<pair<pair<int, int>, int> >::const_iterator it =
           listChangedPlaces.begin();
       it != listChangedPlaces.end(); ++it) {
    //
    if (it->second == 0) {
      ++num1to0;
    } else {
      ++num0to1;
    }
  }
  int num0Tot = this->genosInput.GetGenotypeNumOf(0);
  int num1Tot = this->genosInput.GetGenotypeNumOf(1);
  // cout << "In UpdateEsimate: num0to1: " << num0to1 << ", num1to0: " <<
  // num1to0 << ", num0Tot: " << num0Tot << ", num1Tot: " << num1Tot << endl;
  this->rateFNOpt = ((double)(num0to1 + 1)) / (num0to1 + num1Tot + 2);
  this->rateFPOpt = ((double)(num1to0 + 1)) / (num1to0 + num0Tot + 2);
}
