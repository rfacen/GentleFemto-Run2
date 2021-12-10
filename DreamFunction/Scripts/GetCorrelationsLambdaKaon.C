#include "TROOT.h"
#include "ReadDreamFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TSystem.h"
#include <iostream>

void GetCorrelationsLambdaKaon(const char *filename,
                              const char *prefix, const char *addon = "", double_t norm1 = 0.24, double_t norm2 = 0.34)
{
  //gStyle->SetOptStat(0);
  ReadDreamFile* DreamFile = new ReadDreamFile(4, 4);
  DreamFile->SetAnalysisFile(filename, prefix, addon);

  DreamCF* CF_LK = new DreamCF();/// part-part: ΛK+ + antiΛK-
  DreamCF *CF_ALK = new DreamCF();/// part-antipart: ΛK- + antiΛK+
  DreamCF *CF_LKPlus = new DreamCF();
  DreamCF *CF_ALKMin = new DreamCF();
  DreamCF *CF_LKMin = new DreamCF();
  DreamCF *CF_ALKPlus = new DreamCF();

  DreamPair *LKPlus = new DreamPair("Part", norm1, norm2);
  DreamPair *ALKMin = new DreamPair("AntiPart", norm1, norm2);
  DreamPair *LKMin = new DreamPair("PartAntiPart", norm1, norm2);
  DreamPair *ALKPlus = new DreamPair("AntiPartPart", norm1, norm2);

  DreamCF* CF_LAL = new DreamCF();
  DreamPair* LAL = new DreamPair("PartAntiPart", norm1, norm2);

  std::cout << "=========================" << std::endl;
  std::cout << "========Pair Set=========" << std::endl;
  std::cout << "=========================" << std::endl;
  LKPlus->SetPair(DreamFile->GetPairDistributions(0, 2, ""));
  ALKMin->SetPair(DreamFile->GetPairDistributions(1, 3, ""));

  LKMin->SetPair(DreamFile->GetPairDistributions(1, 2, ""));
  ALKPlus->SetPair(DreamFile->GetPairDistributions(0, 3, ""));

  LAL->SetPair(DreamFile->GetPairDistributions(2, 3, ""));

  std::cout << "=========================" << std::endl;
  std::cout << "======Pair Shifted=======" << std::endl;
  std::cout << "=========================" << std::endl;
  LKPlus->ShiftForEmpty(LKPlus->GetPair());
  ALKMin->ShiftForEmpty(ALKMin->GetPair());

  LKMin->ShiftForEmpty(LKMin->GetPair());
  ALKPlus->ShiftForEmpty(ALKPlus->GetPair());

  LAL->ShiftForEmpty(LAL->GetPair());

  std::cout << "=========================" << std::endl;
  std::cout << "====Pair Fix Shifted=====" << std::endl;
  std::cout << "=========================" << std::endl;
  LKPlus->FixShift(LKPlus->GetPairShiftedEmpty(0), ALKMin->GetPairShiftedEmpty(0),
                   ALKMin->GetFirstBin());
  ALKMin->FixShift(ALKMin->GetPairShiftedEmpty(0), LKPlus->GetPairShiftedEmpty(0),
                   LKPlus->GetFirstBin());

  LKMin->FixShift(LKMin->GetPairShiftedEmpty(0), ALKPlus->GetPairShiftedEmpty(0),
                  ALKPlus->GetFirstBin());
  ALKPlus->FixShift(ALKPlus->GetPairShiftedEmpty(0), LKMin->GetPairShiftedEmpty(0),
                    LKMin->GetFirstBin());

  LAL->FixShift(LAL->GetPairShiftedEmpty(0), LAL->GetPairShiftedEmpty(0),
                LAL->GetFirstBin());

  std::cout << "=========================" << std::endl;
  std::cout << "==Rebinning & Weighting==" << std::endl;
  std::cout << "=========================" << std::endl;

  std::vector<int> rebinVec = {{4, 5}};

  for (int iReb = 4; iReb < 6; ++iReb) {
    std::cout << "==Rebinning LK+==" << std::endl;
    LKPlus->Rebin(LKPlus->GetPairFixShifted(0), rebinVec[iReb]);
    std::cout << "==Weighting LK+==" << std::endl;
    LKPlus->ReweightMixedEvent(LKPlus->GetPairRebinned(iReb-4), 0.2, 0.9);
    std::cout << "==Rebinning ALK-==" << std::endl;
    ALKMin->Rebin(ALKMin->GetPairFixShifted(0), rebinVec[iReb]);
    std::cout << "==Weighting ALK-==" << std::endl;
    ALKMin->ReweightMixedEvent(ALKMin->GetPairRebinned(iReb-4), 0.2, 0.9);
    std::cout << "==Rebinning LK-==" << std::endl;
    LKMin->Rebin(LKMin->GetPairFixShifted(0), rebinVec[iReb]);
    std::cout << "==Weighting LK-==" << std::endl;
    LKMin->ReweightMixedEvent(LKMin->GetPairRebinned(iReb-4), 0.2, 0.9);
    std::cout << "==Rebinning ALK+==" << std::endl;
    ALKPlus->Rebin(ALKPlus->GetPairFixShifted(0), rebinVec[iReb]);
    std::cout << "==Weighting ALK+==" << std::endl;
    ALKPlus->ReweightMixedEvent(ALKPlus->GetPairRebinned(iReb-4), 0.2, 0.9);
  }


  std::cout << "=========================" << std::endl;
  std::cout << "=========CFs=============" << std::endl;
  std::cout << "=========================" << std::endl;

  TString foldername = filename;
  foldername.ReplaceAll("AnalysisResults.root", "");

  std::cout << "==========CF LK===============" << std::endl;
  CF_LK->SetPairs(LKPlus, ALKMin);
  CF_LK->GetCorrelations();
  double nfpairsLK = CF_LK->GetFemtoPairs(0., 0.2);
  CF_LK->WriteOutput(Form("%s/CFOutput_LK_%s.root", foldername.Data(), addon));

  CF_LKPlus->SetPairs(LKPlus,nullptr);
  CF_LKPlus->GetCorrelationsSingle();
  CF_LKPlus->WriteOutput(Form("%s/CFOutput_LKPlus_%s.root", foldername.Data(), addon));

  CF_ALKMin->SetPairs(ALKMin, nullptr);
  CF_ALKMin->GetCorrelationsSingle();
  CF_ALKMin->WriteOutput(Form("%s/CFOutput_CF_ALKMin_%s.root", foldername.Data(), addon));

  CF_ALK->SetPairs(LKMin, ALKPlus);
  CF_ALK->GetCorrelations();
  double nfpairsALK = CF_ALK->GetFemtoPairs(0., 0.2);
  CF_ALK->WriteOutput(Form("%s/CFOutput_ALK_%s.root", foldername.Data(), addon));

  CF_LKMin->SetPairs(LKMin, nullptr);
  CF_LKMin->GetCorrelationsSingle();
  CF_LKMin->WriteOutput(Form("%s/CFOutput_LKMin_%s.root", foldername.Data(), addon));

  CF_ALKPlus->SetPairs(ALKPlus, nullptr);
  CF_ALKPlus->GetCorrelationsSingle();
  CF_ALKPlus->WriteOutput(Form("%s/CFOutput_CF_ALKPlus_%s.root", foldername.Data(), addon));

  std::cout << "=========================" << std::endl;
  std::cout << "Femto Pairs ΛK is " << nfpairsLK << std::endl;
  std::cout << "Femto Pairs ΛK+ is " << CF_LKPlus->GetFemtoPairsBBar(0.,0.2) << std::endl;
  std::cout << "Femto Pairs antiΛK- is " << CF_ALKMin->GetFemtoPairsBBar(0., 0.2) << std::endl;
  std::cout << "=========================" << std::endl;
  std::cout << "Femto Pairs antiΛK is " << nfpairsALK << std::endl;
  std::cout << "Femto Pairs ΛK- is " << CF_LKMin->GetFemtoPairsBBar(0., 0.2) << std::endl;
  std::cout << "Femto Pairs antiΛK+ is " << CF_ALKPlus->GetFemtoPairsBBar(0., 0.2) << std::endl;
  std::cout << "=========================" << std::endl;

  std::cout << "==========CF LAL===============" << std::endl;
  CF_LAL->SetPairs(LAL, nullptr);
  CF_LAL->GetCorrelations();
  CF_LAL->WriteOutput(Form("%s/CFOutput_LAL_%s.root", foldername.Data(), addon));
}
