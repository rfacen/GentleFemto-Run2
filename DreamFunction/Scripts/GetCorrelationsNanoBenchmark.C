#include "TROOT.h"
#include "ReadDreamFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TSystem.h"
#include <iostream>

void GetCorrelationsNanoBenchmark(const char* filename,
                     const char* prefix, const char* addon = "", double norm1 = 0.24, double norm2 = 0.34) {
  //gStyle->SetOptStat(0);
  ReadDreamFile* DreamFile = new ReadDreamFile(6, 6);
  DreamFile->SetAnalysisFile(filename, prefix, addon);

  // Double_t norm1=0.24;//default 0.2-0.4, 0.18-0.28 where CFs look flatter
  // Double_t norm2=0.34;// 0.24-0.34 from BB analysis

  DreamCF* CF_pp = new DreamCF();
  DreamCF *CF_ppOnly = new DreamCF();
  DreamCF *CF_apapOnly = new DreamCF();

  DreamPair* pp = new DreamPair("Part", norm1, norm2);
  DreamPair* ApAp = new DreamPair("AntiPart", norm1, norm2);



  std::cout << "=========================" << std::endl;
  std::cout << "========Pair Set=========" << std::endl;
  std::cout << "=========================" << std::endl;

  pp->SetPair(DreamFile->GetPairDistributions(0, 0, ""));
  ApAp->SetPair(DreamFile->GetPairDistributions(1, 1, ""));

  std::cout << "=========================" << std::endl;
  std::cout << "======Pair Shifted=======" << std::endl;
  std::cout << "=========================" << std::endl;

  std::cout << "==pp==" << std::endl;
  pp->ShiftForEmpty(pp->GetPair());
  ApAp->ShiftForEmpty(ApAp->GetPair());


  std::cout << "=========================" << std::endl;
  std::cout << "====Pair Fix Shifted=====" << std::endl;
  std::cout << "=========================" << std::endl;
  pp->FixShift(pp->GetPairShiftedEmpty(0), ApAp->GetPairShiftedEmpty(0),
               ApAp->GetFirstBin());
  ApAp->FixShift(ApAp->GetPairShiftedEmpty(0), pp->GetPairShiftedEmpty(0),
                 pp->GetFirstBin());


  std::cout << "=========================" << std::endl;
  std::cout << "==Rebinning & Weighting==" << std::endl;
  std::cout << "=========================" << std::endl;
// To be applied only to p-antiL and L-antiL
  pp->ReweightMixedEvent(pp->GetPairFixShifted(0), 0.2, 0.9);
  ApAp->ReweightMixedEvent(ApAp->GetPairFixShifted(0), 0.2, 0.9);


  std::cout << "=========================" << std::endl;
  std::cout << "=========CFs=============" << std::endl;
  std::cout << "=========================" << std::endl;

  TString foldername = filename;
  foldername.ReplaceAll("AnalysisResults.root", "");

  std::cout << "==========CF pAp===============" << std::endl;
  CF_pp->SetPairs(pp, ApAp);
  CF_pp->GetCorrelations();
  CF_pp->WriteOutput(Form("%s/CFOutput_pp_%s.root",foldername.Data(),addon));

  CF_ppOnly->SetPairs(pp, nullptr);
  CF_apapOnly->SetPairs(ApAp, nullptr);
  CF_ppOnly->GetCorrelationsSingle();
  CF_apapOnly->GetCorrelationsSingle();
  CF_ppOnly->WriteOutput(Form("%s/CFOutput_ppOnly_%s.root", foldername.Data(), addon));
  CF_apapOnly->WriteOutput(Form("%s/CFOutput_ApApOnly_%s.root", foldername.Data(), addon));
}
