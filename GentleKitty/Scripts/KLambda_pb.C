#include "TROOT.h"
#include "TSystem.h"
#include "DLM_Source.h"
#include "DLM_Potentials.h"
#include "DLM_CkModels.h"
#include "DLM_WfModel.h"
#include "CATS.h"
#include "DLM_Ck.h"
#include "DLM_CkDecomposition.h"
#include "DLM_Fitters.h"
#include "TidyCats.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TGraph.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TF1.h"
#include "TCanvas.h"
#include "CATSLambdaParam.h"
#include "TMinuit.h"
#include "TMath.h"
#include "TPaveText.h"
#include "DreamPlot.h"
#include "TNtuple.h"
#include "DreamSystematics.h"
#include "TVirtualFitter.h"
#include "TSpline.h"
#include "TDatabasePDG.h"
#include "TString.h"

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <Riostream.h>
#include "TLegendEntry.h"

#include "Math/IFunction.h"
#include <cmath>
#include "TSystem.h"
#include "TLatex.h"

using std::string;
using std::vector;

/// =====================================================================================
void FillCkGraph(DLM_Ck *ck, TGraph *gr)
{
  for (unsigned int i = 0; i < ck->GetNbins(); ++i)
  {
    const float mom = ck->GetBinCenter(0, i);
    gr->SetPoint(i, mom, ck->Eval(mom));
  }
}
/// =====================================================================================
void FillCkGraph(DLM_Ck *ck, DLM_CkDecomposition &ckGraph, TGraph *gr)
{
  for (unsigned int i = 0; i < ck->GetNbins(); ++i)
  {
    const float mom = ck->GetBinCenter(0, i);
    gr->SetPoint(i, mom, ckGraph.EvalCk(mom));
  }
}
/// =====================================================================================
void FillCkGraph(CATS &kitty, TGraph *gr)
{
  for (unsigned int i = 0; i < kitty.GetNumMomBins(); ++i)
  {
    const float mom = kitty.GetMomentum(i);
    gr->SetPoint(i, mom, kitty.GetCorrFun(i));
  }
}

TH2F* GetMatix(const char* rootfile, const char* filename) {
  auto file = TFile::Open(rootfile);
  file->ls();
  TH2F* matrix = (TH2F*)file->FindObjectAny(Form("%s", filename));
  return matrix;

}


TH2F *TransformToMeV(const TH2F *input)
{
  auto histMeV = new TH2F(
      TString::Format("%s_MeV", input->GetName()).Data(), input->GetTitle(),
      input->GetNbinsX(), 1000. * input->GetXaxis()->GetBinLowEdge(1),
      1000. * input->GetXaxis()->GetBinUpEdge(input->GetNbinsX()),
      input->GetNbinsY(), 1000. * input->GetYaxis()->GetBinLowEdge(1),
      1000. * input->GetYaxis()->GetBinUpEdge(input->GetNbinsY()));
  for (int iBinX = 1; iBinX <= input->GetNbinsX(); iBinX++)
  {
    for (int iBinY = 1; iBinY <= input->GetNbinsY(); iBinY++)
    {
      histMeV->SetBinContent(iBinX, iBinY, input->GetBinContent(iBinX, iBinY));
    }
  }
  return histMeV;
}

TH2F *CutMeV(const TH2F *input)
{
  double final_binx = input->GetXaxis()->FindFixBin(1000); //last bin we are interested in, corrisponding to p = 1000 MeV
  double final_biny = input->GetYaxis()->FindFixBin(1000);

  auto histCut = new TH2F(TString::Format("%s_cut", input->GetName()).Data(),
  input->GetTitle(), final_binx, input->GetXaxis()->GetBinLowEdge(1),
  input->GetXaxis()->GetBinUpEdge(input->GetXaxis()->FindFixBin(1000)), final_biny, input->GetYaxis()->GetBinLowEdge(1),
  input->GetYaxis()->GetBinUpEdge(input->GetYaxis()->FindFixBin(1000)));
  for (int iBinX = 1; iBinX <= input->GetNbinsX(); iBinX++)
  {
    for (int iBinY = 1; iBinY <= input->GetNbinsY(); iBinY++)
    {
      histCut->SetBinContent(iBinX, iBinY, input->GetBinContent(iBinX, iBinY));
    }
  }
  return histCut;

}


TF1* fTemplate_MC_NonCommon;
TF1* fTemplate_MC_Common;


double fit_Common_mTThreeGauss(Double_t *x,Double_t *par) {
  double t = *x;
  double fitval =  par[0]*(1.+par[1]*exp(-pow((t-par[2])/par[3],2.)))*(1.+par[4]*exp(-pow((t-par[5])/par[6],2.)))*(1.+par[7]*exp(-pow((t-par[8])/par[9],2.)))*(1.+par[10]*t+par[11]*t*t);
  return fitval;
}


double fit_NonCommon_mT(Double_t *x,Double_t *par) {
  double t = *x;
  double fitval =
  par[0]*(1.+par[1]*exp(-pow((t-par[2])/par[3],2.)))*(1.+par[4]*exp(-pow((t-par[5])/par[6],2.)))*(1.+par[7]*t+par[8]*t*t);
  return fitval;
}


double PreFit_pol1(Double_t *x,Double_t *par) {
  double t = *x;
  double fitval = 
  par[0]*(par[1]*fTemplate_MC_Common->Eval(t) + (1. - par[1])*fTemplate_MC_NonCommon->Eval(t) + (par[2]+par[3]*t));
  return fitval;
}


double PreFit_pol2(Double_t *x,Double_t *par) {
  double t = *x;
  double fitval = 
  par[0]*(par[1]*fTemplate_MC_Common->Eval(t) + (1. - par[1])*fTemplate_MC_NonCommon->Eval(t) + (par[2]+par[3]*t+par[4]*t*t));
  return fitval;
}

DLM_CkDecomposition* FITTER_DECOMLedn;


//total fit
double fit_functionLedn(double* x, double* par) {
    double t = *x;

  FITTER_DECOMLedn->GetCk()->SetPotPar(0, par[0]);//scatlen RE
  FITTER_DECOMLedn->GetCk()->SetPotPar(1, par[1]);//scattlen IM
  FITTER_DECOMLedn->GetCk()->SetPotPar(2, par[2]);//effrange

  //cout<<FITTER_DECOMLedn->GetCk()->GetPotPar(2)<<endl;

  FITTER_DECOMLedn->GetContribution("CkDec_Sig0")->GetCk()->SetPotPar(0, par[0]);//scatlen RE
  FITTER_DECOMLedn->GetContribution("CkDec_Sig0")->GetCk()->SetPotPar(1, par[1]);//scattlen IM
  FITTER_DECOMLedn->GetContribution("CkDec_Sig0")->GetCk()->SetPotPar(2, par[2]);//effrange

  //cout<<FITTER_DECOMLedn->GetContribution("CkDec_Sig0")->GetCk()->GetPotPar(2)<<endl;

   FITTER_DECOMLedn->GetContribution("CkDec_Xi0")->GetCk()->SetPotPar(0, par[0]);//scatlen RE
  FITTER_DECOMLedn->GetContribution("CkDec_Xi0")->GetCk()->SetPotPar(1, par[1]);//scattlen IM
  FITTER_DECOMLedn->GetContribution("CkDec_Xi0")->GetCk()->SetPotPar(2, par[2]);//effrange
  
  //cout<<FITTER_DECOMLedn->GetContribution("CkDec_Xi0")->GetCk()->GetPotPar(2)<<endl;

  FITTER_DECOMLedn->Update(true, true);

  double yval= par[3]*FITTER_DECOMLedn->EvalCk(t)*(par[4]*fTemplate_MC_Common->Eval(t) + (1. - par[4])*fTemplate_MC_NonCommon->Eval(t) + (par[5]+par[6]*t+par[7]*t*t));
  // cout<<*x<<" :"<<FITTER_DECOMLedn->EvalCk(*x)<<endl;
  return yval;
}

int main(int argc, char* argv[]) {
 //include the folders etc that we need
  const char* path_ancestors = "/home/rossanafacen/Analysis/LambdaKaon/TestTask_nosph/Common_Ancestors/Trains_MC";
  const char* path_CF = "/home/rossanafacen/Analysis/LambdaKaon/TestTask_nosph/Classical/Trains/HM/Norm 0.24-0.34";

  double sourceRadius = 1.1;
  double scattlenRE = 0.0;
  double scattlenIM = 0.0;
  double effrange = 0.0;


  const auto pdgDatabase = TDatabasePDG::Instance();
  const double massLambda = pdgDatabase->GetParticle(3122)->Mass() * 1000; 
  const double massKaon = pdgDatabase->GetParticle(321)->Mass() * 1000; 
  const double massXi = pdgDatabase->GetParticle(3312)->Mass() * 1000; 
  const double masspiminus = pdgDatabase->GetParticle(211)->Mass() * 1000; 

  cout<<"mass lambda pdg: "<< massLambda <<endl;
  cout<<"mass kaon pdg: "<< massKaon <<endl;
  cout<<"mass xi- pdg: "<< massXi <<endl;
  cout<<"mass pi- pdg: "<< masspiminus <<endl;

  //const char *pair[4] = {"LKPlus", "CF_ALKMin", "LKMin", "CF_ALKPlus"};
  const char *pair[2] = {"LKPlus+ALKMin", "LKMin+ALKPlus"};
   int i=0;

  double norm1[4] = {0.24, 0.50, 0.34, 0.80}; 
  double chargecombi[4] = {-1, -1, 1, 1}; 
  const char* finalfolder[2] = {"Pb-Pb method/Prefit", "Pb-Pb method/Everything open"};
  const char* path_Results = TString::Format("/home/rossanafacen/Analysis/LambdaKaon/TestTask_nosph/Fit/%s/%s", pair[i],finalfolder[0]);

 
  auto outfile = new TFile(TString::Format("/home/rossanafacen/Analysis/LambdaKaon/TestTask_nosph/Fit/pb_model_%s.root", pair[0]), "RECREATE");   

  //CF of HM
  auto *file_CF = TFile::Open(TString::Format("%s/CFOutput_%s_00_norm_%.2f-%.2f.root", path_CF, pair[i], norm1[i], norm1[i+2])); //it opens the analysis results
  auto *h_HM = (TH1F *)file_CF->FindObjectAny("hCk_ReweightedMeV_1"); 
  

  //common and uncommon CF
  auto *file_common = TFile::Open(TString::Format("%s/CFOutput_%s_0_norm_%.2f-%.2f_Common.root", path_ancestors, pair[0], norm1[0], norm1[0+2]));
  auto *file_uncommon = TFile::Open(TString::Format("%s/CFOutput_%s_0_norm_%.2f-%.2f_Uncommon.root", path_ancestors, pair[0], norm1[0], norm1[0+2])); //it opens the analysis results
  //take a histogram from the root file of CF
  auto *h_MC_Common = (TH1F *)file_common->FindObjectAny("hCk_ReweightedMeV_1"); 
  auto *h_MC_NonCommon = (TH1F *)file_uncommon->FindObjectAny("hCk_ReweightedMeV_1"); 
  
  
  const double LambdaPurity = 0.940459;;
  const double Sigma_feeding = 0.192; //feeddown from the sigma
  const double Csicharged_feeding = 0.232/2; //feeddown from the csi
  const double Csi0_feeding = 0.232/2; //feeddown from the csi

  const double LambdaPrimary = 1.f - Sigma_feeding - Csi0_feeding-Csicharged_feeding;
  const Particle Lambda(LambdaPurity, LambdaPrimary,
                        { { Csicharged_feeding, Csi0_feeding, Sigma_feeding} });


 /// Lambda parameters Kaons
  const double KaonPurity = 0.995;
  const double KaonPrimary = 0.998*0.94;
  const double KaonPhi = 0.998*0.06;

  const double KaonSec = 1.f-KaonPrimary-KaonPhi;

  const Particle Kaon(
      KaonPurity,
      KaonPrimary,
      { { KaonSec, KaonPhi} });


  //const CATSLambdaParam lambdaParam(DBuddy, dmeson);
  const CATSLambdaParam lambdaParam(Kaon, Lambda);

  lambdaParam.PrintLambdaParams();

  //here we have to see what contributions we will combine... this is now nonsense...

  float lPrim = lambdaParam.GetLambdaParam(CATSLambdaParam::Primary);
  float lSigma = lambdaParam.GetLambdaParam(CATSLambdaParam::Primary,
                                                   CATSLambdaParam::FeedDown, 0,
                                                   2);
  float lCsicharged = lambdaParam.GetLambdaParam(CATSLambdaParam::Primary,
                                                   CATSLambdaParam::FeedDown, 0,
                                                   0);  
  float lCsi0 = lambdaParam.GetLambdaParam(CATSLambdaParam::Primary,
                                                   CATSLambdaParam::FeedDown, 0,
                                                   1);                                                                                                
  
  float lFlat = 1.f - lPrim - lCsicharged - lSigma - lCsi0;
  std::cout << "---------------------------------------------------\n";
  std::cout << "---------------------------------------------------\n";

  std::cout << "Lambda parameters\n";

  std::cout << "Genuine K-Lambda " << lPrim << "\n";
  std::cout << "K-Sigma0 -> K-Lambda " << lSigma << "\n";
  std::cout << "K-Xi- -> K-Lambda " << lCsicharged << "\n";
  std::cout << "K-Xi0 -> K-Lambda " << lCsi0 << "\n";
  std::cout << "Flat in K-D+ " << lFlat << "\n";
  std::cout << "---------------------------------------------------\n";
  std::cout << "---------------------------------------------------\n";


  //I build a template of common and common correlation, to model the 2 contributions. I model them with a 3 gaussian
  fTemplate_MC_Common = new TF1(TString::Format("fTemplate_MC_Common"),fit_Common_mTThreeGauss,0,2500,12);
  
  /* fTemplate_MC_Common->SetName(TString::Format("fTemplate_MC_Common_var%i",
  uSyst));*/

  fTemplate_MC_Common->SetParameter(0,0.4);
  fTemplate_MC_Common->SetParLimits(0,0.1,2.);

  fTemplate_MC_Common->SetParameter(1,0.87);
  fTemplate_MC_Common->SetParLimits(1,-2.,2.);

  fTemplate_MC_Common->SetParameter(2,-90.);
  fTemplate_MC_Common->SetParLimits(2,-1000.,100);

  fTemplate_MC_Common->SetParameter(3,300.);
  fTemplate_MC_Common->SetParLimits(3,100.,2000);

  fTemplate_MC_Common->SetParameter(4,0.5);
  fTemplate_MC_Common->SetParLimits(4,-50,+50);

  fTemplate_MC_Common->SetParameter(5,100.);
  fTemplate_MC_Common->SetParLimits(5,-200., +100.);

  fTemplate_MC_Common->SetParameter(6,340);
  fTemplate_MC_Common->SetParLimits(6,0.,2000.);

  fTemplate_MC_Common->SetParameter(7,0.44);
  fTemplate_MC_Common->SetParLimits(7,0.,10.);

  fTemplate_MC_Common->SetParameter(8,400);
  fTemplate_MC_Common->SetParLimits(8,100.,1500.);

  fTemplate_MC_Common->SetParameter(9,285);
  fTemplate_MC_Common->SetParLimits(9,100.,800.);

  fTemplate_MC_Common->SetParameter(10,-0.0000123);
  fTemplate_MC_Common->SetParameter(11,9.9e-8);
  
  h_MC_Common->Fit(fTemplate_MC_Common," S, N, R, M", 0, 2500);
  

  outfile->cd();
  fTemplate_MC_Common->Write("Fit Common");
  fTemplate_MC_Common->SetLineColor(kRed + 2);

  fTemplate_MC_Common->SetName("Fit Common");    
  h_MC_Common->Write("Data Common");
  h_MC_Common->SetLineColor(kOrange + 2); 
 
  double chi_common = fTemplate_MC_Common->GetChisquare()/fTemplate_MC_Common->GetNDF();
  cout<<"Chi2 reduced common: "<< chi_common <<endl;


  fTemplate_MC_NonCommon = new TF1(TString::Format("fTemplate_MC_NonCommon"),fit_NonCommon_mT, 0, 2500, 9);

  fTemplate_MC_NonCommon->SetParameter(0,0.9);
  fTemplate_MC_NonCommon->SetParLimits(0,0.1,2.);

  fTemplate_MC_NonCommon->SetParameter(1,-0.02);
  fTemplate_MC_NonCommon->SetParLimits(1,-1.,+1.);

  fTemplate_MC_NonCommon->SetParameter(2,390);
  fTemplate_MC_NonCommon->SetParLimits(2,-1000.,500.);

  fTemplate_MC_NonCommon->SetParameter(3,185);
  fTemplate_MC_NonCommon->SetParLimits(3,100.,3000.);

  fTemplate_MC_NonCommon->SetParameter(4,0.1);
  fTemplate_MC_NonCommon->SetParLimits(4,1.e-3,3.);

  fTemplate_MC_NonCommon->SetParameter(5,-168);
  fTemplate_MC_NonCommon->SetParLimits(5,-200.,2000.);

  fTemplate_MC_NonCommon->SetParameter(6,216);
  fTemplate_MC_NonCommon->SetParLimits(6,100.,2000.);

  fTemplate_MC_NonCommon->SetParameter(7,0.0001);
  fTemplate_MC_NonCommon->SetParameter(8,3.5e-9);


  h_MC_NonCommon->Fit(fTemplate_MC_NonCommon,"S, N, R, M", 0 , 2500);
  
  fTemplate_MC_NonCommon->SetLineColor(kGreen + 2); 
  fTemplate_MC_NonCommon->SetTitle("Fit NonCommon");


  h_MC_NonCommon->SetLineColor(kViolet + 2); 
  h_MC_NonCommon->SetTitle("Data NonCommon");

  h_MC_NonCommon->Write("Data NonCommon");
  fTemplate_MC_NonCommon->Write("Fit NonCommon");


  double chi_noncommon = fTemplate_MC_NonCommon->GetChisquare()/fTemplate_MC_NonCommon->GetNDF();
  cout<<"Chi2 reduced noncommon: "<<chi_noncommon<<endl;


  TCanvas *c_common_data = new TCanvas();
  h_MC_NonCommon->Draw("");
  h_MC_Common->Draw("same");
  auto legend = new TLegend(0.65, 0.7, 0.85, 0.85, NULL, "brNDC");
  legend->AddEntry(h_MC_Common, "Common Ancestors", "l");
  legend->AddEntry(h_MC_NonCommon, "Uncommon Ancestors", "l");
  legend->Draw();
  h_MC_NonCommon->GetXaxis()->SetRangeUser(0, 1800);
  h_MC_NonCommon->GetYaxis()->SetRangeUser(0.4, 1.3);
  h_MC_NonCommon->GetXaxis()->SetTitle("k* (MeV/c)");
  h_MC_NonCommon->GetYaxis()->SetTitle("C(k*)");
  gStyle->SetOptStat(0);
  c_common_data->Print("/home/rossanafacen/Analysis/LambdaKaon/TestTask_nosph/Fit/Common_data.pdf");

  TCanvas *c_common_fit = new TCanvas();
  fTemplate_MC_NonCommon->Draw("");
  fTemplate_MC_Common->Draw("same"); 
  h_MC_NonCommon->Draw("same");
  h_MC_Common->Draw("same");
  fTemplate_MC_NonCommon->GetXaxis()->SetRangeUser(0, 1800);
  fTemplate_MC_NonCommon->GetYaxis()->SetRangeUser(0.4, 1.3);
  
  legend->Draw();
  gStyle->SetOptStat(0);
  c_common_fit->Print("/home/rossanafacen/Analysis/LambdaKaon/TestTask_nosph/Fit/Common_fit.pdf");

  //prefit of the background (with pol1 and pol2): i have 4 parameters free, and I take the template from common and uncommom
  TF1* fPreFit_pol1 = new TF1("fPreFit_pol1", PreFit_pol1, 0, 2500, 4);
  fPreFit_pol1->SetParNames("N_d", "w_c", "a", "b");

  cout << "\n---------------------------------\n Prefit with pol1 \n--------------------------------\n" << endl;

  h_HM->Fit(fPreFit_pol1, "S, N, R, M", "", 400, 2500);
  double chi_prefit_pol1 = fPreFit_pol1->GetChisquare()/fPreFit_pol1->GetNDF();
  cout<<"Chi2 reduced prefit pol1: "<< chi_prefit_pol1 <<endl;

  fPreFit_pol1->SetLineColor(kMagenta);
  fPreFit_pol1->SetParLimits (1, 0, 1);

  //same thing, but for pol2  
  TF1* fPreFit_pol2 = new TF1("fPreFit_pol2", PreFit_pol2, 0, 2500, 5);
  fPreFit_pol2->SetParameter(0, fPreFit_pol1->GetParameter(0));
  fPreFit_pol2->SetParameter(1, fPreFit_pol1->GetParameter(1));
  fPreFit_pol2->SetParameter(2, fPreFit_pol1->GetParameter(2));
  fPreFit_pol2->SetParameter(3, fPreFit_pol1->GetParameter(3));
  
  fPreFit_pol2->SetLineColor(kGreen + 2); 
  fPreFit_pol2->SetParLimits(1, 5e-2, 1);
  fPreFit_pol2->SetParNames("N_d", "w_c", "a", "b", "c");

  cout << "\n---------------------------------\n Prefit with pol2 \n--------------------------------\n" << endl;

  h_HM->Fit(fPreFit_pol2, "S, N, R, M", "", 400, 2500);
  double chi_prefit_pol2 = fPreFit_pol2->GetChisquare()/fPreFit_pol2->GetNDF();
  cout<<"Chi2 reduced prefit pol2: "<< chi_prefit_pol2 <<endl;

 
  //to study the behaviour of this fit also at low k*, I build a graph, and evaluate the function
  TGraph* gPreFitfull_pol1 = new TGraph();
  gPreFitfull_pol1->SetLineColor(kMagenta +2);
  gPreFitfull_pol1->SetLineWidth(2);

  TGraph* gPreFitfull_pol2 = new TGraph();
  gPreFitfull_pol2->SetLineColor(kGreen);


  for (int x = 0; x < 2500; x++){
    double y_pol1 = fPreFit_pol1->Eval(x);
    gPreFitfull_pol1->SetPoint(x,x,y_pol1);
    double y_pol2 = fPreFit_pol2->Eval(x);
    gPreFitfull_pol2->SetPoint(x,x,y_pol2);
  } 

  gPreFitfull_pol2->Write("PrefitFull_pol2");
  gPreFitfull_pol1->Write("PrefitFull_pol1");
  h_HM->Write("Raw Data");
  fPreFit_pol1->Write("Prefit_pol1");
  fPreFit_pol2->Write("Prefit_pol2");


  TCanvas *c_prefit = new TCanvas();
  h_HM->Draw("");
  gPreFitfull_pol2->Draw("same");
  gPreFitfull_pol1->Draw("same");
  h_HM->GetXaxis()->SetRangeUser(0, 2500);
  h_HM->GetYaxis()->SetRangeUser(0.79, 1.06);
  h_HM->GetXaxis()->SetTitle("k* (MeV/c)");
  h_HM->GetYaxis()->SetTitle("C(k*)");
  auto legend2 = new TLegend(0.65, 0.65, 0.86, 0.86, NULL, "brNDC");
  legend2->AddEntry(h_HM, "Data", "LEP");
  legend2->AddEntry(gPreFitfull_pol1, "Full prefit with Pol1");
  legend2->AddEntry(gPreFitfull_pol2, "Full prefit with Pol2");
  legend2->SetTextSize(gStyle->GetTextSize()*0.6); 
  legend2->Draw();
  gPreFitfull_pol2->SetLineColor(kGreen + 2);
  gPreFitfull_pol1->SetLineColor(kMagenta +2);
  gPreFitfull_pol1->SetLineWidth(2);

  c_prefit->Show();
  c_prefit->Print("/home/rossanafacen/Analysis/LambdaKaon/TestTask_nosph/Fit/Prefit_total.pdf");
 
  int NumMomBins=50;
  int kMin=0;
  int kMax=800;

//model the coulomb
CATS catsXi; //ask dimi if ok for K Xi- coulomb
  catsXi.SetMomBins(NumMomBins, kMin, kMax);
  catsXi.SetUseAnalyticSource(true);
  CATSparameters source_func(CATSparameters::tSource, 1, true);
  catsXi.SetAnaSource(GaussSource, source_func);
  catsXi.SetMomentumDependentSource(false);
  catsXi.SetThetaDependentSource(false);
  catsXi.SetNumChannels(1);
  catsXi.SetQuantumStatistics(0);
  // catsXi.SetNumPW(0, 1);
  // catsXi.SetSpin(0, 0);
  catsXi.SetChannelWeight(0, 1.);
  catsXi.SetQ1Q2(chargecombi[i]);
  catsXi.SetPdgId(321, (chargecombi[i]*3312)); 
  catsXi.SetRedMass((massKaon * massXi) / (massKaon + massXi));
  catsXi.SetAnaSource(0, sourceRadius);
  catsXi.KillTheCat();  

  auto grXiCoulomb = new TGraph();
  grXiCoulomb->SetTitle(";#it{k}* (MeV/#it{c}); #it{C}(#it{k}*)");

  auto grXiSmeared = new TGraph();
  grXiSmeared->SetTitle(";#it{k}* (MeV/#it{c}); #it{C}(#it{k}*)");


  FillCkGraph(catsXi, grXiCoulomb);
  grXiCoulomb->Write("Xi Coulomb before");


  //from feeddown.cpp in downloads TODO a
  
  TString CalibBaseDir = "/home/rossanafacen/Analysis/LambdaKaon/TestTask_nosph/Fit";
      auto calibFile = TFile::Open(TString::Format("%s/histDecayKinematics7.root", CalibBaseDir.Data()));
      auto histDecayKinematicsXi = (TH2F*) calibFile->Get("KXi_KLambda"); //matrix for the feeddown  
      auto histDecayKinematicsSigma = (TH2F*) calibFile->Get("KSigma_KLambda"); //matrix for the feeddown  
      auto histDecayKinematicsXi0 = (TH2F*) calibFile->Get("KXi0_KLambda"); //matrix for the feeddown  

  auto DLM_Xi = new DLM_Ck(1, 0, catsXi);
  DLM_Xi->Update(); 

  const char *path_MC = "/home/rossanafacen/Analysis/LambdaKaon/TestTask_nosph/Classical/Trains/MC/";
  auto *file_MC = TFile::Open(TString::Format("%sAnalysisResults.root", path_MC));

  TList *dirResults;
  
  //find and use the smear matrix, to compute the momentum resolution
  TDirectoryFile *dir = (TDirectoryFile*)(file_MC->FindObjectAny("HMResultsQA00"));
  dir->GetObject("HMResultsQA00", dirResults);
  TList* PairQAList =  (TList*)dirResults->FindObject("PairQA");
  TList* ParticlesList =  (TList*)PairQAList->FindObject("QA_Particle0_Particle2");
  auto Smear_Matrix = (TH2F*)ParticlesList->FindObject("MomentumResolutionME_Particle0_Particle2");   

  //transform to MeV the matrix
  auto Smear_transf = (TH2F*)TransformToMeV(Smear_Matrix);
  //cut the matrix up to 1000 MeV
  auto Smear_transf_cut = (TH2F*)CutMeV(Smear_transf);  

  DLM_CkDecomposition CkDec_Xi_Smeared("KXi_smear", 2,
                                              *DLM_Xi,
                                              Smear_transf_cut);  //TODO mom smearing matrix Ana results
  CkDec_Xi_Smeared.Update();  

  
  FillCkGraph(DLM_Xi, CkDec_Xi_Smeared, grXiSmeared);
  outfile->cd();  
  grXiSmeared->Write("Xi Smeared");

 DLM_Ck *LKModel = new DLM_Ck(1, 3, NumMomBins, kMin, kMax,
                                 ComplexLednicky_Singlet); //only scatt length?
 LKModel->SetSourcePar(0, sourceRadius);
 LKModel->SetPotPar(0, scattlenRE);
 LKModel->SetPotPar(1, scattlenIM);
 LKModel->SetPotPar(2, effrange);
 LKModel->Update();
 //grXiSmeared  

  DLM_Ck *Sig0KModel = new DLM_Ck(1, 3, NumMomBins, kMin, kMax,
                                 ComplexLednicky_Singlet); //only scatt length?
 Sig0KModel->SetSourcePar(0, sourceRadius);
 Sig0KModel->SetPotPar(0, scattlenRE);
 Sig0KModel->SetPotPar(1, scattlenIM);
 Sig0KModel->SetPotPar(2, effrange);
 Sig0KModel->Update();

  DLM_CkDecomposition CkDec_Sig0("CkDec_Sig0", 0,
                                              *Sig0KModel,
                                              nullptr);  //TODO mom smearing matrix Ana results
  CkDec_Sig0.Update();  

  DLM_Ck *Xi0KModel = new DLM_Ck(1, 3, NumMomBins, kMin, kMax,
                                 ComplexLednicky_Singlet); //only scatt length?
 Xi0KModel->SetSourcePar(0, sourceRadius);
 Xi0KModel->SetPotPar(0, scattlenRE);
 Xi0KModel->SetPotPar(1, scattlenIM);
 Xi0KModel->SetPotPar(2, effrange);
 Xi0KModel->Update();

  DLM_CkDecomposition CkDec_Xi0("CkDec_Xi0", 0,
                                              *Xi0KModel,
                                              nullptr);  //TODO mom smearing matrix Ana results
  CkDec_Xi0.Update();  

 DLM_CkDecomposition LKFullCF("LK", 4, *LKModel, Smear_transf_cut);//TODO get mom smearing matrix


 LKFullCF.AddContribution(0, lCsicharged, DLM_CkDecomposition::cFeedDown,
                             &CkDec_Xi_Smeared,histDecayKinematicsXi);
 LKFullCF.AddContribution(1, lSigma, DLM_CkDecomposition::cFeedDown, &CkDec_Sig0, histDecayKinematicsSigma );
 LKFullCF.AddContribution(2, lCsi0, DLM_CkDecomposition::cFeedDown, &CkDec_Xi0, histDecayKinematicsXi0 );
 LKFullCF.AddContribution(3, lFlat, DLM_CkDecomposition::cFeedDown);

 LKFullCF.Update();

 FITTER_DECOMLedn = &LKFullCF;

  cout << "debug0" << endl;

  //total function, using model*background (we already did a prefit before)
  unsigned NumberOfFitPars = 8;            
  TF1* fitter = new TF1("fitter", fit_functionLedn, 0, 500, NumberOfFitPars);            
  fitter->SetParameter(0, scattlenRE);
  fitter->SetParameter(1, scattlenIM); 
  fitter->SetParLimits(1, 0, 100);           
  fitter->SetParameter(2, effrange);
  fitter->SetParLimits(2, 0, 100);                  
  fitter->SetParameter(3, fPreFit_pol2->GetParameter(0)); //constant of normalization 
  //fitter->SetParLimits(3, fPreFit_pol2->GetParameter(0) - fPreFit_pol2->GetParError(0), fPreFit_pol2->GetParameter(0) + fPreFit_pol2->GetParError(0)); //constant of normalization 
  fitter->SetParameter(4, fPreFit_pol2->GetParameter(1)); //weight of common 
  //fitter->SetParLimits(4, fPreFit_pol2->GetParameter(1) - fPreFit_pol2->GetParError(1), fPreFit_pol2->GetParameter(1) + fPreFit_pol2->GetParError(1)); //constant of normalization 
  fitter->SetParLimits(4, 0, 1);                  
  fitter->FixParameter(5, fPreFit_pol2->GetParameter(2)); //pol2  coefficients, to be kept fixed from the perfit    
  fitter->FixParameter(6, fPreFit_pol2->GetParameter(3));      
  fitter->FixParameter(7, fPreFit_pol2->GetParameter(4)); 
  fitter->SetParNames("Re(f0)","Im(f0)","d0", "N_d", "w_c", "a", "b", "c");


  cout << "\n---------------------------------\n Total fit \n--------------------------------\n" << endl;
  int status= h_HM->Fit(fitter, "NR", "", 0, 500);

  double chi_fitter = fitter->GetChisquare()/fitter->GetNDF();
  cout<<"Chi2 reduced: "<< chi_fitter <<endl;

  DLM_Ck *LKModelResult = new DLM_Ck(1, 3, NumMomBins, kMin, kMax,
                                  ComplexLednicky_Singlet); //only scatt length?
  LKModelResult->SetSourcePar(0, sourceRadius);
  LKModelResult->SetPotPar(0, fitter->GetParameter(0));
  LKModelResult->SetPotPar(1, fitter->GetParameter(1));
  LKModelResult->SetPotPar(2, fitter->GetParameter(1));
  LKModelResult->Update();
  
  auto  grGenuine_new = new TGraph();
  FillCkGraph(LKModelResult, grGenuine_new);

 
  TF1* fitter_bg = new TF1("fitter_bg", PreFit_pol2, 0, 2500, 5);
  fitter_bg->FixParameter(0, fitter->GetParameter(3));
  fitter_bg->FixParameter(1, fitter->GetParameter(4));
  fitter_bg->FixParameter(2, fitter->GetParameter(5));
  fitter_bg->FixParameter(3, fitter->GetParameter(6));
  fitter_bg->FixParameter(4, fitter->GetParameter(7));
  fitter_bg->SetParNames("N_d", "w_c", "a", "b", "c");
  
  fitter_bg->SetLineColor(fitter->GetLineColor());
  fitter_bg->SetLineStyle(6);

 
  TCanvas *c_fitter = new TCanvas();
  h_HM->Draw("");
  fitter->Draw("same");
  h_HM->GetXaxis()->SetRangeUser(0, 2500);
  h_HM->GetYaxis()->SetRangeUser(0.79, 1.06);
  h_HM->GetXaxis()->SetTitle("k* (MeV/c)");
  h_HM->GetYaxis()->SetTitle("C(k*)");
  auto legend_fit = new TLegend(0.65, 0.65, 0.86, 0.86, NULL, "brNDC");
  legend_fit->AddEntry(h_HM, "Data", "LEP");
  legend_fit->AddEntry(fitter, "Fit with Lednicky");
  legend_fit->SetTextSize(gStyle->GetTextSize()*0.6); 
  legend_fit->Draw();
  c_fitter->Show();
  c_fitter->Print("/home/rossanafacen/Analysis/LambdaKaon/TestTask_nosph/Fit/Lednic_fit.pdf");

  Smear_Matrix->Write("Original Matrix");
  Smear_transf->Write("Matrix MeV transformed");
  Smear_transf_cut->Write("Matrix cut");

  fitter->Write("fitter");
  fitter_bg->Write("fitter full");
  grXiCoulomb->Write();
  grXiSmeared->Write();

  double scattlenRE_pb[2] = {-0.6, 0.27}; double scattlenRE_pb_err[2] = {0.12, 0.12};
  double scattlenIM_pb[2] = {0.51, 0.40}; double  scattlenIM_pb_err[2] = {0.15, 0.11}; 
  double effrange_pb[2] = {0.83, -5.23}; double effrange_pb_err[2] = {0.47, 2.13}; 

  double FitMin = 0.;
  double FitMax = 500.;
  //fit with Pb parameters
  TF1* fitter_pb = new TF1("fitter_pb", fit_functionLedn, FitMin, FitMax, NumberOfFitPars);            
  fitter_pb->FixParameter(0, scattlenRE_pb[i]);
  //fitter_pb->SetParLimits(0, scattlenRE_pb[i] - scattlenRE_pb_err[i], scattlenRE_pb [i] + scattlenRE_pb_err[i]);
  fitter_pb->FixParameter(1, scattlenIM_pb[i]);            
  //fitter_pb->SetParLimits(1, scattlenIM_pb[i] -scattlenIM_pb_err[i], scattlenIM_pb[i] + scattlenIM_pb_err[i]);
  fitter_pb->FixParameter(2, effrange_pb[i]);            
  //fitter_pb->SetParLimits(2, effrange_pb[i] - effrange_pb_err[i], effrange_pb[i] + effrange_pb_err[i]);
  fitter_pb->SetParameter(3, fPreFit_pol2->GetParameter(0)); //constant of normalization           
  fitter_pb->SetParameter(4, fPreFit_pol2->GetParameter(1)); //weight of common 
  fitter_pb->SetParLimits(4, 0, 1); //weight of common 
  fitter_pb->FixParameter(5, fPreFit_pol2->GetParameter(2)); //pol2  coefficients, to be kept fixed from the perfit    
  fitter_pb->FixParameter(6, fPreFit_pol2->GetParameter(3));      
  fitter_pb->FixParameter(7, fPreFit_pol2->GetParameter(4));
  fitter_pb->SetParNames("Re(f0)_pb","Im(f0)_pb","d0_pb", "N_d", "w_c", "a", "b", "c");

  fitter_pb->SetLineColor(kAzure+7);
  
  h_HM->Fit(fitter_pb, "MR", "", FitMin, FitMax);

  DLM_Ck *LKModelResultPb = new DLM_Ck(1, 3, NumMomBins, kMin, kMax,
                                  ComplexLednicky_Singlet); //only scatt length?
  LKModelResultPb->SetSourcePar(0, sourceRadius);
  LKModelResultPb->SetPotPar(0, fitter_pb->GetParameter(0));
  LKModelResultPb->SetPotPar(1, fitter_pb->GetParameter(1));
  LKModelResultPb->SetPotPar(2, fitter_pb->GetParameter(1));
  LKModelResultPb->Update();
  
  auto  grGenuine_new_Pb = new TGraph();
  FillCkGraph(LKModelResultPb, grGenuine_new_Pb);

  outfile->cd();
  grGenuine_new_Pb->Write("grGenuine_new_Pb");
  
  TGraph* gfitterfull_pb = new TGraph();
  gfitterfull_pb->SetLineColor(fitter_pb->GetLineColor());
  gfitterfull_pb->SetLineWidth(2);
  gfitterfull_pb->SetLineStyle(6);

  for (int x = 0; x < 2500; x++){
    double y = fitter_pb->Eval(x);
    gfitterfull_pb->SetPoint(x,x,y);
    }

  
  TCanvas *c_fitter_pb = new TCanvas();
  h_HM->Draw("");
  fitter_pb->Draw("same");
  h_HM->GetXaxis()->SetRangeUser(0, 2500);
  h_HM->GetYaxis()->SetRangeUser(0.87, 1.06);
  h_HM->GetXaxis()->SetTitle("k* (MeV/c)");
  h_HM->GetYaxis()->SetTitle("C(k*)");
  auto legend_fit_pb = new TLegend(0.65, 0.65, 0.86, 0.86, NULL, "brNDC");
  legend_fit_pb->AddEntry(h_HM, "Data", "LEP");
  legend_fit_pb->AddEntry(fitter_pb, "Fit with Pb parms");
  legend_fit_pb->SetTextSize(gStyle->GetTextSize()*0.6); 
  legend_fit_pb->Draw();
  c_fitter_pb->Show();
  c_fitter_pb->Print("/home/rossanafacen/Analysis/LambdaKaon/TestTask_nosph/Fit/Pb_fit.pdf");

  double chi_fitter_pb = fitter_pb->GetChisquare()/fitter_pb->GetNDF();
  cout<<"Chi2 reduced pb: "<< chi_fitter_pb <<endl;

  std::fstream output5;
  output5.open(Form("/home/rossanafacen/Analysis/LambdaKaon/TestTask_nosph/Fit/outputPAR_%s_Pb_everything_open.dat", pair[i]), std::fstream::in | std::fstream::out | std::fstream::app);
  for(int i=0;i<NumberOfFitPars;i++){
    output5 << fitter->GetParName(i)<<" "<<fitter->GetParameter(i)<<" "<<fitter->GetParError(i)<< "\n";
  }

  output5 << "--------------------------------\n";
  output5 << "With Pb-Pb parameters \n";
  for(int i=0;i<NumberOfFitPars;i++){
    output5 << fitter->GetParName(i)<<" "<<fitter_pb->GetParameter(i)<<" "<<fitter_pb->GetParError(i)<< "\n";
  }
  output5.close();
  fitter_pb->Write("fitter pb");
  gfitterfull_pb->Write("fitter full pb");
  outfile->Close();
  return 0;

}

