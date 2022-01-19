#include "ForgivingReader.h"
#include "MakeHistosGreat.h"
#include "EventQA.h"
#include "TrackQA.h"
#include "DecayQA.h"

int main(int argc, char* argv[]) {
  const char* filename = argv[1];
  const char *prefix = argv[2];
  const char* addon = (argv[3]) ? argv[3] : "";
  MakeHistosGreat::SetStyle(false);
  ForgivingReader* reader = new ForgivingReader(filename, prefix, addon);
  auto file = reader->GetFile();
  auto dirRes = file->GetDirectory(Form("%sResults%s", prefix, addon));
  auto dirQA = file->GetDirectory(Form("%sResultsQA%s", prefix, addon));
  TList *listQA;
  dirQA->GetObject(Form("%sResultsQA%s", prefix, addon), listQA);

  EventQA* evtQA = new EventQA();
  evtQA->SetLooseMargin();
  evtQA->SetQAList(reader->GetQA());
  evtQA->SetEventCuts(reader->GetEventCuts());
  int nevts = evtQA->GetNumberOfEvents();
  std::cout << "#######################" << std::endl;
  std::cout << "Number of events is " << nevts << std::endl;
  std::cout << "#######################" << std::endl;
  evtQA->PlotEventCounter(prefix,addon);
  evtQA->PlotEventProperties(200);
  evtQA->PlotPileUpRejection();
  evtQA->SetTightMargin();
  evtQA->SetQAList(listQA);
  evtQA->PlotStatsTrackCleaner({}, {"#Lambda-#varphi", "#bar{#Lambda}-#varphi"},7);

  //Kaons
  auto dirTracks = file->GetDirectory(Form("%sKaonPlusCuts%s", prefix, addon));
  TList *listTracks;
  dirTracks->GetObject(Form("%sKaonPlusCuts%s", prefix, addon), listTracks);
  TrackQA *KaonPlusQA = new TrackQA();
  KaonPlusQA->SetTrackCuts(listTracks);
  KaonPlusQA->PlotKinematicTracks("KaonPlus");
  KaonPlusQA->PlotPIDTracks("KaonPlus");

  auto dirAntiTracks = file->GetDirectory(Form("%sKaonMinusCuts%s", prefix, addon));
  TList *listAntiTracks;
  dirAntiTracks->GetObject(Form("%sKaonMinusCuts%s", prefix, addon), listAntiTracks);
  TrackQA *KaonMinusQA = new TrackQA();
  KaonMinusQA->SetAntiTrackCuts(listAntiTracks);
  KaonMinusQA->PlotKinematicAntiTracks("KaonMinus");
  KaonMinusQA->PlotPIDAntiTracks("KaonMinus");


  //Lambda
  auto dirLambda = file->GetDirectory(Form("%sLambdaCuts%s", prefix, addon));
  TList *listLambda;
  dirLambda->GetObject(Form("%sLambdaCuts%s", prefix, addon), listLambda);
  DecayQA *v0QA = new DecayQA("#Lambda", "p#pi");
  v0QA->SetCanvasDivisions(4, 2);
  v0QA->SetDecayCuts(listLambda);
  v0QA->SetIMHistoScale(1.75, 0.8, 0.35);
  v0QA->SetRangesFitting(1.1075, 1.1235, 1.09, 1.15);
  v0QA->InvariantMassPartLambda(1.112, 1.120, false, 0.48, 0.515);
  v0QA->PlotQATopologyPartLambda();

  auto dirAntiLambda = file->GetDirectory(Form("%sAntiLambdaCuts%s", prefix, addon));
  TList *listAntiLambda;
  dirAntiLambda->GetObject(Form("%sAntiLambdaCuts%s", prefix, addon), listAntiLambda);
  DecayQA *antiv0QA = new DecayQA("#bar{#Lambda}", "#bar{p}#pi");
  antiv0QA->SetCanvasDivisions(4, 2);
  antiv0QA->SetIMHistoScale(1.75, 0.8, 0.35);
  antiv0QA->SetAntiDecayCuts(listAntiLambda);
  antiv0QA->SetRangesFitting(1.1075, 1.1235, 1.09, 1.15);
  antiv0QA->InvariantMassAntiPartLambda(1.112, 1.120, false, 0.48, 0.515);
  antiv0QA->PlotQATopologyAntiPartLambda();

  //Phi
  auto dirPhi = file->GetDirectory(Form("%sPhiCuts%s", prefix, addon));
  TList *listPhi;
  dirPhi->GetObject(Form("%sPhiCuts%s", prefix, addon), listPhi);
  float MphiPDG = 1.019461;
  DecayQA *PhiQA = new DecayQA("#phi", "K^{-}K^{+}");
  PhiQA->SetCanvasDivisions(4, 2);
  PhiQA->SetDecayCuts(listPhi);
  PhiQA->SetRangesFitting(0.99, 1.08, 0.99, 1.08);
  PhiQA->InvariantMassPartPhi(0.008);
  PhiQA->PlotQATopologyLambda(listPhi, "Phi");
}
