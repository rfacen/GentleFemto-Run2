#include "TROOT.h"
#include "ReadDreamFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TSystem.h"
#include <iostream>

void GetCorrelationsLKCommonAncestors(const char* filename,
                     const char* prefix, const char* addon = "",
                     double_t norm1 = 0.24, double_t norm2 = 0.34) {
  //gStyle->SetOptStat(0);
  std::cout << "bug1" << std::endl;
  ReadDreamFile* DreamFile = new ReadDreamFile(4, 4);
  char* addon2 = Form("(null)%s", addon);
  DreamFile->SetAnalysisFileAncestors(filename, prefix, addon2);

  
 
  //create the CFs
  DreamCF* CF_LKPlus_common = new DreamCF();
  DreamCF* CF_ALKMin_common = new DreamCF();
  DreamCF* CF_LKMin_common = new DreamCF();
  DreamCF* CF_ALKPlus_common = new DreamCF();
  DreamCF* CF_attractive_common = new DreamCF(); //sum of the 2
  DreamCF* CF_repulsive_common = new DreamCF(); //sum of the 2

  DreamCF* CF_LKPlus_uncommon = new DreamCF();
  DreamCF* CF_ALKMin_uncommon = new DreamCF();
  DreamCF* CF_LKMin_uncommon = new DreamCF();
  DreamCF* CF_ALKPlus_uncommon = new DreamCF();
  DreamCF* CF_attractive_uncommon = new DreamCF();
  DreamCF* CF_repulsive_uncommon = new DreamCF();

  //create the pairs
  DreamPair* LKPlus_common = new DreamPair("Part", norm1, norm2);
  DreamPair* LKPlus_uncommon = new DreamPair("Part", norm1, norm2);
  DreamPair* ALKMin_common = new DreamPair("AntiPart", norm1, norm2);
  DreamPair* ALKMin_uncommon = new DreamPair("AntiPart", norm1, norm2);
  DreamPair* LKMin_common = new DreamPair("PartAntiPart", norm1, norm2);
  DreamPair* LKMin_uncommon = new DreamPair("PartAntiPart", norm1, norm2);  
  DreamPair* ALKPlus_common = new DreamPair("AntiPartAntiPart", norm1, norm2);
  DreamPair* ALKPlus_uncommon = new DreamPair("AntiPartAntiPart", norm1, norm2);

  DreamCF* CF_LAL_common = new DreamCF(); //check because we already studied Lambda-Antilambda
  DreamCF* CF_LAL_uncommon = new DreamCF(); //check because we already studied Lambda-Antilambda
  
  DreamPair* LAL_common = new DreamPair("PartAntiPart", norm1, norm2);
  DreamPair* LAL_uncommon = new DreamPair("PartAntiPart", norm1, norm2);

  std::cout << "=========================" << std::endl;
  std::cout << "========Pair Set=========" << std::endl;
  std::cout << "=========================" << std::endl;

  //get pair distribution, for common and uncommon ancestors
  LKPlus_common->SetPair(DreamFile->GetPairDistributionsCommon(0, 2, ""));
  LKPlus_uncommon->SetPair(DreamFile->GetPairDistributionsNonCommon(0, 2, ""));
  
  ALKMin_common->SetPair(DreamFile->GetPairDistributionsCommon(1, 3, ""));
  ALKMin_uncommon->SetPair(DreamFile->GetPairDistributionsNonCommon(1, 3, ""));
  
  LKMin_common->SetPair(DreamFile->GetPairDistributionsCommon(1, 2, ""));
  LKMin_uncommon->SetPair(DreamFile->GetPairDistributionsNonCommon(1, 2, ""));
  
  ALKPlus_common->SetPair(DreamFile->GetPairDistributionsCommon(0, 3, ""));
  ALKPlus_uncommon->SetPair(DreamFile->GetPairDistributionsNonCommon(0, 3, ""));

  LAL_common->SetPair(DreamFile->GetPairDistributionsCommon(2, 3, ""));
  LAL_uncommon->SetPair(DreamFile->GetPairDistributionsNonCommon(2, 3, ""));


  std::cout << "=========================" << std::endl;
  std::cout << "======Pair Shifted=======" << std::endl;
  std::cout << "=========================" << std::endl;

  LKPlus_common->ShiftForEmpty(LKPlus_common->GetPair());
  LKPlus_uncommon->ShiftForEmpty(LKPlus_uncommon->GetPair());
  
  ALKMin_common->ShiftForEmpty(ALKMin_common->GetPair());
  ALKMin_uncommon->ShiftForEmpty(ALKMin_uncommon->GetPair());

  LKMin_common->ShiftForEmpty(LKMin_common->GetPair());
  LKMin_uncommon->ShiftForEmpty(LKMin_uncommon->GetPair());
  
  ALKPlus_common->ShiftForEmpty(ALKPlus_common->GetPair());
  ALKPlus_uncommon->ShiftForEmpty(ALKPlus_uncommon->GetPair());

  LAL_common->ShiftForEmpty(LAL_common->GetPair());
  LAL_uncommon->ShiftForEmpty(LAL_uncommon->GetPair());

  
  std::cout << "=========================" << std::endl;
  std::cout << "====Pair Fix Shifted=====" << std::endl;
  std::cout << "=========================" << std::endl;



  LKPlus_common->FixShift(LKPlus_common->GetPairShiftedEmpty(0), ALKMin_common->GetPairShiftedEmpty(0),
                   ALKMin_common->GetFirstBin());
  ALKMin_common->FixShift(ALKMin_common->GetPairShiftedEmpty(0), LKPlus_common->GetPairShiftedEmpty(0),
                   LKPlus_common->GetFirstBin());

  ALKPlus_common->FixShift(ALKPlus_common->GetPairShiftedEmpty(0), LKMin_common->GetPairShiftedEmpty(0),
                   LKMin_common->GetFirstBin());
  LKMin_common->FixShift(LKMin_common->GetPairShiftedEmpty(0), ALKPlus_common->GetPairShiftedEmpty(0),
                   ALKPlus_common->GetFirstBin());


  /* LKMin_common->FixShift(LKMin_common->GetPairShiftedEmpty(0), ALKPlus_common->GetPairShiftedEmpty(0),
                  ALKPlus_common->GetFirstBin());
  ALKPlus_common->FixShift(ALKPlus_common->GetPairShiftedEmpty(0), LKMin_common->GetPairShiftedEmpty(0),
                    LKMin_common->GetFirstBin()); */

  LAL_common->FixShift(LAL_common->GetPairShiftedEmpty(0), LAL_common->GetPairShiftedEmpty(0),
                LAL_common->GetFirstBin());


  LKPlus_uncommon->FixShift(LKPlus_uncommon->GetPairShiftedEmpty(0), ALKMin_uncommon->GetPairShiftedEmpty(0),
                   ALKMin_uncommon->GetFirstBin());
  ALKMin_uncommon->FixShift(ALKMin_uncommon->GetPairShiftedEmpty(0), LKPlus_uncommon->GetPairShiftedEmpty(0),
                   LKPlus_uncommon->GetFirstBin());

  LKMin_uncommon->FixShift(LKMin_uncommon->GetPairShiftedEmpty(0), ALKPlus_uncommon->GetPairShiftedEmpty(0),
                  ALKPlus_uncommon->GetFirstBin());
  ALKPlus_uncommon->FixShift(ALKPlus_uncommon->GetPairShiftedEmpty(0), LKMin_uncommon->GetPairShiftedEmpty(0),
                    LKMin_uncommon->GetFirstBin());

  LAL_uncommon->FixShift(LAL_uncommon->GetPairShiftedEmpty(0), LAL_uncommon->GetPairShiftedEmpty(0),
                LAL_uncommon->GetFirstBin());


  std::cout << "=========================" << std::endl;
  std::cout << "==Rebinning & Weighting==" << std::endl;
  std::cout << "=========================" << std::endl;
  
  std::vector<int> rebinVec = {{1,4, 5}};

  // To be applied only to p-antiL and L-antiL
  for (size_t iReb = 0; iReb < rebinVec.size(); ++iReb) {
    std::cout << "==Rebinning LK+==" << std::endl;
    LKPlus_common->Rebin(LKPlus_common->GetPairFixShifted(0), rebinVec[iReb]);
    std::cout << "==Weighting LK+==" << std::endl;
    LKPlus_common->ReweightMixedEvent(LKPlus_common->GetPairRebinned(iReb), 0.2, 0.9);
    std::cout << "==Rebinning LK+==" << std::endl;
 
    LKPlus_uncommon->Rebin(LKPlus_uncommon->GetPairFixShifted(0), rebinVec[iReb]);
    LKPlus_uncommon->ReweightMixedEvent(LKPlus_uncommon->GetPairRebinned(iReb), 0.2, 0.9);


    ALKMin_common->Rebin(ALKMin_common->GetPairFixShifted(0), rebinVec[iReb]);
    std::cout << "==Weighting ALK-==" << std::endl;
    ALKMin_common->ReweightMixedEvent(ALKMin_common->GetPairRebinned(iReb), 0.2, 0.9);
    std::cout << "==Rebinning ALK-==" << std::endl;
    
    ALKMin_uncommon->Rebin(ALKMin_uncommon->GetPairFixShifted(0), rebinVec[iReb]);
    ALKMin_uncommon->ReweightMixedEvent(ALKMin_uncommon->GetPairRebinned(iReb), 0.2, 0.9);
    

    LKMin_common->Rebin(LKMin_common->GetPairFixShifted(0), rebinVec[iReb]);
    std::cout << "==Weighting LK-==" << std::endl;
    LKMin_common->ReweightMixedEvent(LKMin_common->GetPairRebinned(iReb), 0.2, 0.9);
    std::cout << "==Rebinning LK-==" << std::endl;
    
    LKMin_uncommon->Rebin(LKMin_uncommon->GetPairFixShifted(0), rebinVec[iReb]);
    LKMin_uncommon->ReweightMixedEvent(LKMin_uncommon->GetPairRebinned(iReb), 0.2, 0.9);
    

    /* ALKPlus_common->Rebin(ALKPlus_common->GetPairFixShifted(0), rebinVec[iReb]);
    ALKPlus_common->ReweightMixedEvent(ALKPlus_common->GetPairRebinned(iReb), 0.2, 0.9);
     */
    ALKPlus_common->Rebin(ALKPlus_common->GetPairFixShifted(0), rebinVec[iReb]);
    ALKPlus_common->ReweightMixedEvent(ALKPlus_common->GetPairRebinned(iReb), 0.2, 0.9);


    ALKPlus_uncommon->Rebin(ALKPlus_uncommon->GetPairFixShifted(0), rebinVec[iReb]);
    ALKPlus_uncommon->ReweightMixedEvent(ALKPlus_uncommon->GetPairRebinned(iReb), 0.2, 0.9);
  }
  
  std::cout << "=========================" << std::endl;
  std::cout << "=========CFs=============" << std::endl;
  std::cout << "=========================" << std::endl;
  

 TString foldername = filename;
  foldername.ReplaceAll("AnalysisResults.root", "");
  //  TString foldername = gSystem->pwd();

  std::cout << "==========CF LK===============" << std::endl;

  CF_LKPlus_common->SetPairs(LKPlus_common,nullptr);
  CF_LKPlus_common->GetCorrelationsSingle();
  CF_LKPlus_common->WriteOutput(Form("%s/CFOutput_LKPlus_%s_norm_%.2f-%.2f_Common.root", foldername.Data(), addon, norm1, norm2));
  
  CF_LKPlus_uncommon->SetPairs(LKPlus_uncommon,nullptr);
  CF_LKPlus_uncommon->GetCorrelationsSingle();
  CF_LKPlus_uncommon->WriteOutput(Form("%s/CFOutput_LKPlus_%s_norm_%.2f-%.2f_Uncommon.root", foldername.Data(), addon, norm1, norm2));

  CF_ALKMin_common->SetPairs(ALKMin_common, nullptr);
  CF_ALKMin_common->GetCorrelationsSingle();
  CF_ALKMin_common->WriteOutput(Form("%s/CFOutput_ALKMin_%s_norm_%.2f-%.2f_Common.root", foldername.Data(), addon, norm1, norm2));

  CF_ALKMin_uncommon->SetPairs(ALKMin_uncommon, nullptr);
  CF_ALKMin_uncommon->GetCorrelationsSingle();
  CF_ALKMin_uncommon->WriteOutput(Form("%s/CFOutput_ALKMin_%s_norm_%.2f-%.2f_Uncommon.root", foldername.Data(), addon, norm1, norm2));


  CF_repulsive_common->SetPairs(LKPlus_common, ALKMin_common);
  CF_repulsive_common->GetCorrelations();
  CF_repulsive_common->WriteOutput(Form("%s/CFOutput_repulsive_%s_norm_%.2f-%.2f_Common.root", foldername.Data(), addon, norm1, norm2));

  CF_repulsive_uncommon->SetPairs(LKPlus_uncommon, ALKMin_uncommon);
  CF_repulsive_uncommon->GetCorrelations();
  CF_repulsive_uncommon->WriteOutput(Form("%s/CFOutput_repulsive_%s_norm_%.2f-%.2f_Uncommon.root", foldername.Data(), addon, norm1, norm2));

 
  CF_LKMin_common->SetPairs(LKMin_common, nullptr);
  CF_LKMin_common->GetCorrelationsSingle();
  CF_LKMin_common->WriteOutput(Form("%s/CFOutput_LKMin_%s_norm_%.2f-%.2f_Common.root", foldername.Data(), addon, norm1, norm2));

  CF_LKMin_uncommon->SetPairs(LKMin_uncommon, nullptr);
  CF_LKMin_uncommon->GetCorrelationsSingle();
  CF_LKMin_uncommon->WriteOutput(Form("%s/CFOutput_LKMin_%s_norm_%.2f-%.2f_Uncommon.root", foldername.Data(), addon, norm1, norm2));


  CF_ALKPlus_common->SetPairs(ALKPlus_common, nullptr);
  CF_ALKPlus_common->GetCorrelationsSingle();
  CF_ALKPlus_common->WriteOutput(Form("%s/CFOutput_ALKPlus_%s_norm_%.2f-%.2f_Common.root", foldername.Data(), addon, norm1, norm2));

  CF_ALKPlus_uncommon->SetPairs(ALKPlus_uncommon, nullptr);
  CF_ALKPlus_uncommon->GetCorrelationsSingle();
  CF_ALKPlus_uncommon->WriteOutput(Form("%s/CFOutput_ALKPlus_%s_norm_%.2f-%.2f_Uncommon.root", foldername.Data(), addon, norm1, norm2));
  std::cout << "bug3" << std::endl;

  CF_attractive_common->SetPairs(ALKPlus_common, LKMin_common);
  CF_attractive_common->GetCorrelations();
  CF_attractive_common->WriteOutput(Form("%s/CFOutput_attractive_%s_norm_%.2f-%.2f_Common.root", foldername.Data(), addon, norm1, norm2));

  std::cout << "bug4" << std::endl;


  CF_attractive_uncommon->SetPairs(LKMin_uncommon, ALKPlus_uncommon);
  CF_attractive_uncommon->GetCorrelations();
  CF_attractive_uncommon->WriteOutput(Form("%s/CFOutput_attractive_%s_norm_%.2f-%.2f_Uncommon.root", foldername.Data(), addon, norm1, norm2));
  std::cout << "bug5" << std::endl;


  double nfpairsLK = CF_LKMin_common->GetFemtoPairs(0., 0.2);
  double nfpairsALK = CF_ALKMin_common->GetFemtoPairs(0., 0.2);

  std::cout << "=========================" << std::endl;
  std::cout << "Femto Pairs ΛK is " << nfpairsLK << std::endl;
  std::cout << "Femto Pairs ΛK+ is " << CF_LKPlus_common->GetFemtoPairsBBar(0.,0.2) << std::endl;
  std::cout << "Femto Pairs antiΛK- is " << CF_ALKMin_common->GetFemtoPairsBBar(0., 0.2) << std::endl;
  std::cout << "=========================" << std::endl;
  std::cout << "Femto Pairs antiΛK is " << nfpairsALK << std::endl;
  std::cout << "Femto Pairs ΛK- is " << CF_LKMin_common->GetFemtoPairsBBar(0., 0.2) << std::endl;
  std::cout << "Femto Pairs antiΛK+ is " << CF_ALKPlus_common->GetFemtoPairsBBar(0., 0.2) << std::endl;
  std::cout << "=========================" << std::endl;

  std::cout << "==========CF LAL===============" << std::endl;
  CF_LAL_common->SetPairs(LAL_common, nullptr);
  CF_LAL_common->GetCorrelations();
  CF_LAL_common->WriteOutput(Form("%s/CFOutput_LAL_%s_norm_%.2f-%.2f_Ancestors.root", foldername.Data(), addon, norm1, norm2));
}
