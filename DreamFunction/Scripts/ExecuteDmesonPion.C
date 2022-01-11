#include "ReadDreamFile.h"
#include "TLegend.h"
#include "METoSEReweighting.C"

int main(int argc, char* argv[]) {
  const char* filename = argv[1];
  TString appendix = TString::Format("%s", argv[2]);
  TString suffix = TString::Format("%s", argv[3]);
  TString filenamesuffix = TString::Format("%s", argv[4]);
  ReadDreamFile* DreamFile = new ReadDreamFile(4, 4);
  DreamFile->SetAnalysisFile(filename, appendix.Data(), suffix.Data());

  const double normLower = 1.0;
  const double normUpper = 1.5;

  DreamCF* CF_PIplusDminus = new DreamCF();
  DreamPair* PIplusDminus = new DreamPair("Part_opp", normLower, normUpper);
  DreamPair* PIminusDplus = new DreamPair("AntiPart_opp", normLower, normUpper);

  DreamCF* CF_PIplusDplus = new DreamCF();
  DreamPair* PIplusDplus = new DreamPair("Part_same", normLower, normUpper);
  DreamPair* PIminusDminus = new DreamPair("AntiPart_same", normLower, normUpper);


  PIplusDminus->SetPair(DreamFile->GetPairDistributions(0, 3, ""));
  PIminusDplus->SetPair(DreamFile->GetPairDistributions(1, 2, ""));

 /* PIplusDminus->ShiftForEmpty(PIplusDminus->GetPair());
  PIminusDplus->ShiftForEmpty(PIminusDplus->GetPair());*/

  PIplusDplus->SetPair(DreamFile->GetPairDistributions(0, 2, ""));
  PIminusDminus->SetPair(DreamFile->GetPairDistributions(1, 3, ""));

/*  PIplusDplus->ShiftForEmpty(PIplusDplus->GetPair());
  PIminusDminus->ShiftForEmpty(PIminusDminus->GetPair());*/

/*  PIplusDminus->FixShift(PIplusDminus->GetPairShiftedEmpty(0),
                    PIminusDplus->GetPairShiftedEmpty(0), PIminusDplus->GetFirstBin());
  PIminusDplus->FixShift(PIminusDplus->GetPairShiftedEmpty(0),
                    PIplusDminus->GetPairShiftedEmpty(0), PIplusDminus->GetFirstBin());*/

  PIplusDminus->FixShift(PIplusDminus->GetPair(), PIminusDplus->GetPair(), 0.0, true);
  PIminusDplus->FixShift(PIminusDplus->GetPair(), PIplusDminus->GetPair(), 0.0, true);

  PIplusDplus->FixShift(PIplusDplus->GetPair(), PIminusDminus->GetPair(), 0.0, true);
  PIminusDminus->FixShift(PIminusDminus->GetPair(), PIplusDplus->GetPair(), 0.0, true);


/*  PIplusDplus->FixShift(PIplusDplus->GetPairShiftedEmpty(0),
                    PIminusDminus->GetPairShiftedEmpty(0), PIminusDminus->GetFirstBin());
  PIminusDminus->FixShift(PIminusDminus->GetPairShiftedEmpty(0),
                    PIplusDplus->GetPairShiftedEmpty(0), PIplusDplus->GetFirstBin());*/

  std::vector<int> rebin = { { 1, 4, 5, 10, 20} };

  for (size_t iReb = 0; iReb < rebin.size(); ++iReb) {
    PIplusDminus->Rebin(PIplusDminus->GetPairFixShifted(0), rebin[iReb], true);
    PIminusDplus->Rebin(PIminusDplus->GetPairFixShifted(0), rebin[iReb], true);

    PIplusDminus->ReweightMixedEvent(PIplusDminus->GetPairRebinned(iReb), 0.2, 0.9,
                                PIplusDminus->GetPair());
    PIminusDplus->ReweightMixedEvent(PIminusDplus->GetPairRebinned(iReb), 0.2, 0.9,
                                PIminusDplus->GetPair());

    PIplusDplus->Rebin(PIplusDplus->GetPairFixShifted(0), rebin[iReb], true);
    PIminusDminus->Rebin(PIminusDminus->GetPairFixShifted(0), rebin[iReb], true);

    PIplusDplus->ReweightMixedEvent(PIplusDplus->GetPairRebinned(iReb), 0.2, 0.9,
                                PIplusDplus->GetPair());
    PIminusDminus->ReweightMixedEvent(PIminusDminus->GetPairRebinned(iReb), 0.2, 0.9,
                                PIminusDminus->GetPair());
  }

  TString foldername = filename;
  foldername.ReplaceAll("AnalysisResults_all.root", "");

  TString fileAppendix =
      (suffix == "0") ? "" : TString::Format("_%s", suffix.Data());

  CF_PIplusDminus->SetPairs(PIplusDminus, PIminusDplus, true);
  CF_PIplusDminus->GetCorrelations();
  CF_PIplusDminus->WriteOutput(
      Form("%s/CFOutput_PIplusDminus%s_%s.root", foldername.Data(),
           fileAppendix.Data(), filenamesuffix.Data()));


  CF_PIplusDplus->SetPairs(PIplusDplus, PIminusDminus, true);
  CF_PIplusDplus->GetCorrelations();
  CF_PIplusDplus->WriteOutput(
      Form("%s/CFOutput_PIplusDplus%s_%s.root", foldername.Data(),
           fileAppendix.Data(), filenamesuffix.Data()));

/*  TString FileName = Form("%s/CFOutput_PIplusDminus%s.root", foldername.Data(),
                          fileAppendix.Data());
  TFile* file = TFile::Open(FileName, "update");
  TList* PairDist = (TList*) file->Get("PairDist");

  if (PairDist) {
    ReweightingQA(PairDist);
  } else {
    file->ls();
  }
  TList* AntiPairDist = (TList*) file->Get("AntiPairDist");
  if (AntiPairDist) {
    ReweightingQA(AntiPairDist);
  }
  return 1;*/
}
