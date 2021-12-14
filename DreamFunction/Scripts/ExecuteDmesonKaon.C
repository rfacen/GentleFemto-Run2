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

  const double normLower = 1.5;
  const double normUpper = 2.;

  DreamCF* CF_KplusDminus = new DreamCF();
  DreamPair* KplusDminus = new DreamPair("Part_opp", normLower, normUpper);
  DreamPair* KminusDplus = new DreamPair("AntiPart_opp", normLower, normUpper);

  DreamCF* CF_KplusDplus = new DreamCF();
  DreamPair* KplusDplus = new DreamPair("Part_same", normLower, normUpper);
  DreamPair* KminusDminus = new DreamPair("AntiPart_same", normLower, normUpper);


  KplusDminus->SetPair(DreamFile->GetPairDistributions(0, 3, ""));
  KminusDplus->SetPair(DreamFile->GetPairDistributions(1, 2, ""));

 /* KplusDminus->ShiftForEmpty(KplusDminus->GetPair());
  KminusDplus->ShiftForEmpty(KminusDplus->GetPair());*/

  KplusDplus->SetPair(DreamFile->GetPairDistributions(0, 2, ""));
  KminusDminus->SetPair(DreamFile->GetPairDistributions(1, 3, ""));

/*  KplusDplus->ShiftForEmpty(KplusDplus->GetPair());
  KminusDminus->ShiftForEmpty(KminusDminus->GetPair());*/

/*  KplusDminus->FixShift(KplusDminus->GetPairShiftedEmpty(0),
                    KminusDplus->GetPairShiftedEmpty(0), KminusDplus->GetFirstBin());
  KminusDplus->FixShift(KminusDplus->GetPairShiftedEmpty(0),
                    KplusDminus->GetPairShiftedEmpty(0), KplusDminus->GetFirstBin());*/

  KplusDminus->FixShift(KplusDminus->GetPair(), KminusDplus->GetPair(), 0.0, true);
  KminusDplus->FixShift(KminusDplus->GetPair(), KplusDminus->GetPair(), 0.0, true);

  KplusDplus->FixShift(KplusDplus->GetPair(), KminusDminus->GetPair(), 0.0, true);
  KminusDminus->FixShift(KminusDminus->GetPair(), KplusDplus->GetPair(), 0.0, true);


/*  KplusDplus->FixShift(KplusDplus->GetPairShiftedEmpty(0),
                    KminusDminus->GetPairShiftedEmpty(0), KminusDminus->GetFirstBin());
  KminusDminus->FixShift(KminusDminus->GetPairShiftedEmpty(0),
                    KplusDplus->GetPairShiftedEmpty(0), KplusDplus->GetFirstBin());*/

  std::vector<int> rebin = { { 1, 4, 5, 10, 20} };

  for (size_t iReb = 0; iReb < rebin.size(); ++iReb) {
    KplusDminus->Rebin(KplusDminus->GetPairFixShifted(0), rebin[iReb], true);
    KminusDplus->Rebin(KminusDplus->GetPairFixShifted(0), rebin[iReb], true);

    KplusDminus->ReweightMixedEvent(KplusDminus->GetPairRebinned(iReb), 0.2, 0.9,
                                KplusDminus->GetPair());
    KminusDplus->ReweightMixedEvent(KminusDplus->GetPairRebinned(iReb), 0.2, 0.9,
                                KminusDplus->GetPair());

    KplusDplus->Rebin(KplusDplus->GetPairFixShifted(0), rebin[iReb], true);
    KminusDminus->Rebin(KminusDminus->GetPairFixShifted(0), rebin[iReb], true);

    KplusDplus->ReweightMixedEvent(KplusDplus->GetPairRebinned(iReb), 0.2, 0.9,
                                KplusDplus->GetPair());
    KminusDminus->ReweightMixedEvent(KminusDminus->GetPairRebinned(iReb), 0.2, 0.9,
                                KminusDminus->GetPair());
  }

  TString foldername = filename;
  foldername.ReplaceAll("AnalysisResults_all.root", "");

  TString fileAppendix =
      (suffix == "0") ? "" : TString::Format("_%s", suffix.Data());

  CF_KplusDminus->SetPairs(KplusDminus, KminusDplus, true);
  CF_KplusDminus->GetCorrelations();
  CF_KplusDminus->WriteOutput(
      Form("%s/CFOutput_KplusDminus%s_%s.root", foldername.Data(),
           fileAppendix.Data(), filenamesuffix.Data()));


  CF_KplusDplus->SetPairs(KplusDplus, KminusDminus, true);
  CF_KplusDplus->GetCorrelations();
  CF_KplusDplus->WriteOutput(
      Form("%s/CFOutput_KplusDplus%s_%s.root", foldername.Data(),
           fileAppendix.Data(), filenamesuffix.Data()));

/*  TString FileName = Form("%s/CFOutput_KplusDminus%s.root", foldername.Data(),
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
