#include "TROOT.h"
#include "ReadDreamFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TSystem.h"
#include <iostream>

void GetCorrelationsLambdaPhi(const char *filename,
                              const char *prefix, const char *addon = "", double_t norm1 = 0.24, double_t norm2 = 0.34)
{
  //gStyle->SetOptStat(0);
  ReadDreamFile* DreamFile = new ReadDreamFile(5, 5);
  DreamFile->SetAnalysisFile(filename, prefix, addon);

  DreamCF* CF_LPhi_ALPhi = new DreamCF();
  DreamCF *CF_LPhi = new DreamCF();
  DreamCF *CF_ALPhi = new DreamCF();
  DreamPair *LPhi = new DreamPair("PartAntiPart", norm1, norm2);
  DreamPair *ALPhi = new DreamPair("AntiPartPart", norm1, norm2);

  DreamCF* CF_LAL = new DreamCF();
  DreamPair* LAL = new DreamPair("PartAntiPart", norm1, norm2);

  std::cout << "=========================" << std::endl;
  std::cout << "========Pair Set=========" << std::endl;
  std::cout << "=========================" << std::endl;
  LPhi->SetPair(DreamFile->GetPairDistributions(0, 2, ""));
  ALPhi->SetPair(DreamFile->GetPairDistributions(1, 2, ""));

  LAL->SetPair(DreamFile->GetPairDistributions(0, 1, ""));

  std::cout << "=========================" << std::endl;
  std::cout << "======Pair Shifted=======" << std::endl;
  std::cout << "=========================" << std::endl;
  LPhi->ShiftForEmpty(LPhi->GetPair());
  ALPhi->ShiftForEmpty(ALPhi->GetPair());

  LAL->ShiftForEmpty(LAL->GetPair());

  std::cout << "=========================" << std::endl;
  std::cout << "====Pair Fix Shifted=====" << std::endl;
  std::cout << "=========================" << std::endl;

  LPhi->FixShift(LPhi->GetPairShiftedEmpty(0), LPhi->GetPairShiftedEmpty(0),
                 LPhi->GetFirstBin());
  ALPhi->FixShift(ALPhi->GetPairShiftedEmpty(0), LPhi->GetPairShiftedEmpty(0),
                  LPhi->GetFirstBin());
  LAL->FixShift(LAL->GetPairShiftedEmpty(0), LAL->GetPairShiftedEmpty(0),
                LAL->GetFirstBin());

  std::cout << "=========================" << std::endl;
  std::cout << "==Rebinning & Weighting==" << std::endl;
  std::cout << "=========================" << std::endl;

  for (int iReb = 4; iReb < 6; ++iReb) {
    std::cout << "==Rebinning==" << std::endl;
    LPhi->Rebin(LPhi->GetPairFixShifted(0), iReb);
    LAL->Rebin(LAL->GetPairFixShifted(0), iReb);
    std::cout << "==Weighting==" << std::endl;
    LPhi->ReweightMixedEvent(LPhi->GetPairRebinned(iReb - 4), 0.2, 0.9);
    LAL->ReweightMixedEvent(LAL->GetPairRebinned(iReb - 4), 0.2, 0.9);
    std::cout << "==Rebinning==" << std::endl;
    ALPhi->Rebin(ALPhi->GetPairFixShifted(0), iReb);
    std::cout << "==Weighting==" << std::endl;
    ALPhi->ReweightMixedEvent(ALPhi->GetPairRebinned(iReb - 4), 0.2, 0.9);
  }


  std::cout << "=========================" << std::endl;
  std::cout << "=========CFs=============" << std::endl;
  std::cout << "=========================" << std::endl;

  TString foldername = filename;
  foldername.ReplaceAll("AnalysisResults.root", "");

  std::cout << "==========CF LPhi===============" << std::endl;
  CF_LPhi_ALPhi->SetPairs(LPhi, ALPhi);
  CF_LPhi_ALPhi->GetCorrelations();
  int nfpairs = CF_LPhi_ALPhi->GetFemtoPairs(0.,0.2);
  CF_LPhi_ALPhi->WriteOutput(Form("%s/CFOutput_LPhi_%s.root", foldername.Data(), addon));

  CF_LPhi->SetPairs(LPhi,nullptr);
  CF_LPhi->GetCorrelations();
  int nfpairs_Lphi = CF_LPhi->GetFemtoPairsBBar(0., 0.2);
  CF_LPhi->WriteOutput(Form("%s/CFOutput_LPhiOnly_%s.root", foldername.Data(), addon));

  CF_ALPhi->SetPairs(ALPhi, nullptr);
  CF_ALPhi->GetCorrelations();
  int nfpairs_ALphi = CF_ALPhi->GetFemtoPairsBBar(0., 0.2);
  CF_ALPhi->WriteOutput(Form("%s/CFOutput_ALPhiOnly_%s.root", foldername.Data(), addon));

  std::cout << "=========================" << std::endl;
  std::cout << "Femto Pairs Λφ is " << nfpairs_Lphi << std::endl;
  std::cout << "Femto Pairs antiΛφ is " << nfpairs_ALphi << std::endl;
  std::cout << "Femto Pairs inclusive is " << nfpairs << std::endl;
  std::cout << "=========================" << std::endl;

  std::cout << "==========CF LAL===============" << std::endl;
  CF_LAL->SetPairs(LAL, nullptr);
  CF_LAL->GetCorrelations();
  CF_LAL->WriteOutput(Form("%s/CFOutput_LAL_%s.root", foldername.Data(), addon));
}
