/*
 * DreamCF.h
 *
 *  Created on: Aug 22, 2018
 *      Author: hohlweger
 */

#ifndef DREAMCF_H_
#define DREAMCF_H_
#include <vector>
#include "DreamPair.h"
#include "TRatioPlot.h"
#include "TH1F.h"

class DreamCF {
 public:
  DreamCF();
  virtual ~DreamCF();
  void SetPairs(DreamPair* pairOne, DreamPair* pairTwo, const bool& direct_sum=false) {
    fPairOne = pairOne;
    fPairTwo = pairTwo;
    fDirectSum = direct_sum;
  }
  ;
  void SetPairsBBar(DreamPair* pairOne) {
    fPairOne = pairOne;
  }
  ;
  void GetCorrelations(const char* pairName = "");
  void LoopCorrelations(std::vector<DreamDist*> PairOne,
                        std::vector<DreamDist*> PairTwo, const char* name);
  void LoopCorrelations(std::vector<DreamDist*> Pair, const char* name);
  void WriteOutput(const char* name);
  void WriteOutput(TFile* output, bool closeFile);
  //required of rebinned data!
  TH1F* AddCF(DreamDist* DD1, DreamDist* DD2, const char* name);
  TH1F* AddCF(TH1F* CF1, TH1F* CF2, const char* name);
  static TGraphAsymmErrors* AddCF(TH1F* histSum, std::vector<DreamPair*> pairs, const char* name);
  static TH1F* ConvertToOtherUnit(TH1F* HistCF, int Scale, const char* name);
  static TGraphAsymmErrors* ConvertToOtherUnit(TGraphAsymmErrors* HistCF, int Scale, const char* name);
  std::vector<TH1F*> GetCorrelationFunctions() {
    return fCF;
  }
  ;
  std::vector<TGraphAsymmErrors*> GetCorrelationFunctionGraphs() {
    return fGrCF;
  }
  unsigned int GetFemtoPairs(float kMin, float kMax) {
    return
        (fPairOne && fPairTwo) ?
            (fPairOne->GetFemtoPairs(kMin, kMax)
                + fPairTwo->GetFemtoPairs(kMin, kMax)) :
            0;
  }
  unsigned int GetFemtoPairsBBar(float kMin, float kMax) {
    return (fPairOne) ? (fPairOne->GetFemtoPairs(kMin, kMax)) : 0;
  }
  DreamPair* GetPairOne() {
    return fPairOne;
  }
  DreamPair* GetPairTwo() {
    return fPairTwo;
  }
  TH1F* FindCorrelationFunction(TString name);
 private:
  std::vector<TH1F*> fCF;
  std::vector<TGraphAsymmErrors*> fGrCF;
  std::vector<TH1F*> fRatio;
  DreamPair* fPairOne;
  DreamPair* fPairTwo;
  //if false, the two correlations will be combined using the
  //weighted mean based on the errors of the correlations (WRONG)
  //the above option has to be kept the default one for compatibility reasons
  //in any future analysis this HAS TO BE TRUE!!!
  bool fDirectSum;
};

#endif /* DREAMCF_H_ */
