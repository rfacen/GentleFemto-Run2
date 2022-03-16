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

TH1F* TransformToMev1D(TH1F *input){
  auto histMeV1D = new TH1F(TString::Format("%s_MeV", input->GetName()).Data(), input->GetTitle(),
  input->GetNbinsX(), input->GetXaxis()->GetBinLowEdge(1) * 1000, input->GetXaxis()->GetBinUpEdge(input->GetNbinsX())*1000);
  
  for(int i = 1; i <= input->GetNbinsX(); i++){
    histMeV1D->SetBinContent(i, input->GetBinContent(i));
  }
  return histMeV1D;
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

double fit_Common_mTThreeGauss_gau(Double_t *x,Double_t *par) {
  double t = *x;
  double fitval =  par[0]*(1.+par[1]*exp(-pow((t-par[2])/par[3],2.)))*(1.+par[4]*exp(-pow((t-par[5])/par[6],2.)))
  *(1.+par[7]*exp(-pow((t-par[8])/par[9],2.)))*(1.+par[10]*t+par[11]*t*t) + par[12]*exp(-pow((t-par[13])/par[14],2.))
  +par[15]*exp(-pow((t-par[16])/par[17],2.));
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
  FITTER_DECOMLedn->Update(true, true);

  double yval= par[3]*FITTER_DECOMLedn->EvalCk(t)*(par[4]*fTemplate_MC_Common->Eval(t) + (1. - par[4])*fTemplate_MC_NonCommon->Eval(t) + (par[5]+par[6]*t+par[7]*t*t));
  // cout<<*x<<" :"<<FITTER_DECOMLedn->EvalCk(*x)<<endl;
  return yval;
}


  TH1F *getBootstrapHisto(TH1F *histo){
  auto histoOut = (TH1F *)histo->Clone(Form("bootstrap_%s_%i", histo->GetName(),int(gRandom->Uniform() * 10000.f)));
  static double xVal, yVal;
  
  for (int i = 1; i <= histo->GetNbinsX(); ++i)
  {
    //histo->GetPoint(i, xVal, yVal);
    yVal = histo->GetBinContent(i);
    histoOut->SetBinContent(i, gRandom->Gaus(yVal, histo->GetBinError(i)));
  }
  return histoOut;
}

int main(int argc, char* argv[]) {

/*       //flag
     char* fit   
     if (fit == "pol2") {
      cout << "Pol2" << endl;
      systvar = false;

  } else if (fit == "pol1") {
      cout << "Pol1" << endl;
      systvar = true;

  } else if (fit == "pb ") {
      cout << "Pb comparision" << endl;
      evalchi2 = true;

  } else {
      cout<<"error"<<endl;
  }
 */

  //include the folders etc that we need
  const char* path_ancestors = "/home/rossanafacen/Analysis/LambdaKaon/TestTask_nosph/Common_Ancestors/Trains_MC";
  const char* path_CF = "/home/rossanafacen/Analysis/LambdaKaon/TestTask_nosph/Classical/Trains/HM";

  double sourceRadius = 1.11251; //
  //L+K+(1.354):   low 1.07421 mean 1.11251 up 1.14903
  //L+K-(1.355):   low 1.07366 mean 1.11193 up 1.14843

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
  const char* norm[2] = {"Norm 0.24-0.34", "Norm 0.50-0.80"};
  
  double norm1[4] = {0.24, 0.50, 0.34, 0.80}; 
  //double chargecombi[4] = {-1, -1, 1, 1}; 
  double chargecombi[2] = {-1, 1}; 

  int i=0;

  auto outfile = new TFile(TString::Format("/home/rossanafacen/Analysis/LambdaKaon/TestTask_nosph/Fit/%s/%s_Fit.root", pair[i], pair[i]), "RECREATE");   
  const double LambdaPurity = 0.940459;;
  const double Sigma_feeding = 0.192; //feedwon from the sigma
  const double Csicharged_feeding = 0.232/2; //feedown from the csi
  const double Csi0_feeding = 0.232/2; //feedown from the csi

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
  
  float lFlat = 1.f - lPrim - lCsicharged;
  std::cout << "---------------------------------------------------\n";
  std::cout << "---------------------------------------------------\n";

  std::cout << "Lambda parameters\n";

  std::cout << "Genuine K-Lambda " << lPrim << "\n";
  //std::cout << "K-Sigma0 -> K-Lambda " << lSigma << "\n";
  std::cout << "K-Xi- -> K-Lambda " << lCsicharged << "\n";
  //std::cout << "K-Xi0 -> K-Lambda " << lCsi0 << "\n";
  std::cout << "Flat in K-D+ " << lFlat << "\n";
  std::cout << "---------------------------------------------------\n";
  std::cout << "---------------------------------------------------\n";


  auto *file_CF = TFile::Open(TString::Format("%s/%s/CFOutput_%s_00_norm_%.2f-%.2f.root", path_CF, norm[i], pair[i], norm1[i], norm1[i+2])); //it opens the analysis results
  auto *h_HM_orig = (TH1F *)file_CF->FindObjectAny("hCk_ReweightedMeV_1"); 
  //common and uncommon CF
  auto *file_common = TFile::Open(TString::Format("%s/CFOutput_%s_0_norm_%.2f-%.2f_Common.root", path_ancestors, pair[i], norm1[i], norm1[i+2]));
  auto *file_uncommon = TFile::Open(TString::Format("%s/CFOutput_%s_0_norm_%.2f-%.2f_Uncommon.root", path_ancestors, pair[i], norm1[i], norm1[i+2])); //it opens the analysis results
  //take a histogram from the root file of CF
  auto *h_MC_Common_orig = (TH1F *)file_common->FindObjectAny("hCk_ReweightedMeV_1"); 
  auto *h_MC_NonCommon_orig = (TH1F *)file_uncommon->FindObjectAny("hCk_ReweightedMeV_1");   

  int nbins = h_HM_orig->GetNbinsX();
  
  //arrays which save the results of the fit, with the dimension of our histogram: we ll save the y value here after the fit
  unsigned* entries = new unsigned[nbins];
  double* Fit_mean = new double[nbins];
  double* Fit_StDev = new double[nbins];  
  
  // initialize to 0
  for (unsigned uBin = 0; uBin < nbins; uBin++) {
    Fit_mean[uBin] = 0;
    Fit_StDev[uBin] = 0;
    entries[uBin] = 0;
  }
  
  int boot = 1000;

  //histograms which save the scattering parameters for every bootstrap
  auto* Re_scatlen = new TH1F("Re_scatlen", "Re_scatlen", 30, -0.5, 0.5);
  Re_scatlen->GetXaxis()->SetTitle("Re(f_{0})");
  Re_scatlen->GetYaxis()->SetTitle("Entries");  
  
  auto* Im_scatlen = new TH1F("Im_scatlen", "Im_scatlen", 30, 0, 0.3);
  Im_scatlen->GetXaxis()->SetTitle("Im(f_{0})");
  Im_scatlen->GetYaxis()->SetTitle("Entries");  
  
  auto* effran = new TH1F("effran", "effran", 30, -0, 1);
  effran->GetXaxis()->SetTitle("d_{0}");
  effran->GetYaxis()->SetTitle("Entries"); 
  
  auto* chi2ndf = new TH1F("chi2ndf", "chi2ndf", 10, 0.0, 40.0);
  chi2ndf->GetXaxis()->SetTitle("#Chi^{2}/ndf");
  chi2ndf->GetYaxis()->SetTitle("Entries");



  for(int j = 0; j < boot; j++){
      auto *h_HM = (TH1F*)getBootstrapHisto(h_HM_orig);      
      auto *h_MC_Common = (TH1F*)getBootstrapHisto(h_MC_Common_orig);
      auto *h_MC_NonCommon = (TH1F*)getBootstrapHisto(h_MC_NonCommon_orig);

      //I build a template of common and common correlation, to model the 2 contributions. I model them with a 3 gaussian
      if(i == 0) fTemplate_MC_Common = new TF1("fTemplate_MC_Common",fit_Common_mTThreeGauss, 0, 2500, 12);
      else if (i == 1) fTemplate_MC_Common = new TF1("fTemplate_MC_Common",fit_Common_mTThreeGauss_gau, 0, 2500, 18);

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
      fTemplate_MC_Common->SetParameter(12,0.2); 
      fTemplate_MC_Common->SetParLimits(12, 0, 0.3);
      fTemplate_MC_Common->SetParameter(13, 220); //position of the peak of the omega resonance
      fTemplate_MC_Common->SetParLimits(13, 210, 230);      
      fTemplate_MC_Common->SetParameter(14, 10); //width of the peak 
      fTemplate_MC_Common->SetParLimits(14, 2, 20);      
      fTemplate_MC_Common->SetParameter(15,0.5); 
      fTemplate_MC_Common->SetParLimits(15, 0, 1);
      fTemplate_MC_Common->SetParameter(16, 220); //position of the peak of the omega resonance
      fTemplate_MC_Common->SetParLimits(16, 210, 230);     
      fTemplate_MC_Common->SetParameter(17, 10); //width of the peak 
      fTemplate_MC_Common->SetParLimits(17, 2, 20);
      
      cout << "center : " << h_MC_Common->GetXaxis()->GetBinCenter(2) << endl;
      h_MC_Common->Fit(fTemplate_MC_Common, "S, N, R, M", "", h_MC_Common->GetXaxis()->GetBinCenter(2), 2500);
      
      outfile->cd();
      fTemplate_MC_Common->Write("Fit Common");
      fTemplate_MC_Common->SetLineColor(kRed + 2);   
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

      h_MC_NonCommon->Fit(fTemplate_MC_NonCommon,"S, N, R, M", "", h_MC_NonCommon->GetXaxis()->GetBinCenter(2), 2500);
      
      fTemplate_MC_NonCommon->SetLineColor(kGreen + 2); 
      fTemplate_MC_NonCommon->SetTitle("Fit NonCommon");

      h_MC_NonCommon->SetLineColor(kViolet + 2); 
      h_MC_NonCommon->SetTitle("Data NonCommon");

      h_MC_NonCommon->Write("Data NonCommon");
      fTemplate_MC_NonCommon->Write("Fit NonCommon");

      double chi_noncommon = fTemplate_MC_NonCommon->GetChisquare()/fTemplate_MC_NonCommon->GetNDF();
      cout<<"Chi2 reduced noncommon: "<<chi_noncommon<<endl;

      /* TCanvas *c_common_data = new TCanvas();
      h_MC_NonCommon->Draw("");
      h_MC_Common->Draw("same");
      auto legend = new TLegend(0.65, 0.7, 0.88, 0.88, NULL, "brNDC");
      legend->AddEntry(h_MC_Common, "Common Ancestors", "l");
      legend->AddEntry(h_MC_NonCommon, "Uncommon Ancestors", "l");
      legend->SetTextSize(gStyle->GetTextSize()*0.55);
      legend->Draw();
      h_MC_NonCommon->GetXaxis()->SetRangeUser(0, 1800);
      h_MC_NonCommon->GetYaxis()->SetRangeUser(0.4, 2);
      h_MC_NonCommon->SetTitle("; k* (MeV/c); C(k*)");
      gStyle->SetOptStat(0);
      c_common_data->Print(TString::Format("/home/rossanafacen/Analysis/LambdaKaon/TestTask_nosph/Fit/%s/Common_data_%s.pdf", pair[i], pair[i]));
 */
     /*  TCanvas *c_common_fit = new TCanvas();
      fTemplate_MC_NonCommon->Draw("");
      fTemplate_MC_Common->Draw("same"); 
      h_MC_NonCommon->Draw("same");
      h_MC_Common->Draw("same");
      fTemplate_MC_NonCommon->GetXaxis()->SetRangeUser(0, 1800);
      fTemplate_MC_NonCommon->GetYaxis()->SetRangeUser(0.4, 1.3);
      fTemplate_MC_NonCommon->SetTitle("; k* (MeV/c); C(k*)");
      
      legend->Draw();
      gStyle->SetOptStat(0);
      c_common_fit->Print(TString::Format("/home/rossanafacen/Analysis/LambdaKaon/TestTask_nosph/Fit/%s/Common_fit_%s.pdf", pair[i], pair[i]));
 */

      TGraph* gPreFitfull_pol1 = new TGraph();
      gPreFitfull_pol1->SetLineColor(kMagenta +2);
      gPreFitfull_pol1->SetLineWidth(2);

      TGraph* gPreFitfull_pol2 = new TGraph();
      gPreFitfull_pol2->SetLineColor(kGreen);

      //prefit of the background (with pol1 and pol2): i have 4 parameters free, and I take the template from common and uncommom
      TF1* fPreFit_pol1 = new TF1("fPreFit_pol1", PreFit_pol1, 350, 2000, 4);
      
      fPreFit_pol1->SetParNames("N_d", "w_c", "a", "b");

      cout << "\n---------------------------------\n Prefit with pol1 \n--------------------------------\n" << endl;

      h_HM->Fit(fPreFit_pol1, "S, N, R, M", "", 350, 2000); 
      
      double chi_prefit_pol1 = fPreFit_pol1->GetChisquare()/fPreFit_pol1->GetNDF();
      cout<<"Chi2 reduced prefit pol1: "<< chi_prefit_pol1 <<endl;

      fPreFit_pol1->SetLineColor(kMagenta);
      fPreFit_pol1->SetParLimits (1, 0, 1);

      //same thing, but for pol2  
      TF1* fPreFit_pol2 = new TF1("fPreFit_pol2", PreFit_pol2, 350, 2000, 5);
      fPreFit_pol2->SetParameter(0, fPreFit_pol1->GetParameter(0));
      fPreFit_pol2->SetParameter(1, fPreFit_pol1->GetParameter(1));
      fPreFit_pol2->SetParameter(2, fPreFit_pol1->GetParameter(2));
      fPreFit_pol2->SetParameter(3, fPreFit_pol1->GetParameter(3));


      fPreFit_pol2->SetLineColor(kGreen + 2); 
      fPreFit_pol2->SetParLimits(1, 5e-2, 1);
      fPreFit_pol2->SetParNames("N_d", "w_c", "a", "b", "c");

      cout << "\n---------------------------------\n Prefit with pol2 \n--------------------------------\n" << endl;

      h_HM->Fit(fPreFit_pol2, "S, N, R, M", "", 350, 2000);
      double chi_prefit_pol2 = fPreFit_pol2->GetChisquare()/fPreFit_pol2->GetNDF();
      cout<<"Chi2 reduced prefit pol2: "<< chi_prefit_pol2 <<endl;

      //to study the behaviour of this fit also at low k*, I build a graph, and evaluate the function
      
      for (int x = 0; x < 2000; x++){
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
      h_HM->SetTitle("; k* (MeV/c); C(k*)");
      auto legend2 = new TLegend(0.65, 0.65, 0.86, 0.86, NULL, "brNDC");
      legend2->AddEntry(h_HM, "Data", "LEP");
      legend2->AddEntry(gPreFitfull_pol1, "Full prefit with Pol1");
      legend2->AddEntry(gPreFitfull_pol2, "Full prefit with Pol2");
      legend2->SetTextSize(gStyle->GetTextSize()*0.6); 
      legend2->Draw();
      gPreFitfull_pol2->SetLineColor(kGreen + 2);
      gPreFitfull_pol1->SetLineColor(kMagenta +2);
      gPreFitfull_pol1->SetLineWidth(2);

      /* c_prefit->Show();
      c_prefit->Print(TString::Format("/home/rossanafacen/Analysis/LambdaKaon/TestTask_nosph/Fit/%s/Prefit_total_%s.pdf", pair[i], pair[i]));
     */
      int NumMomBins=50;
      int kMin=0;
      int kMax=800;


      //model the coulomb, which consists on K+ Xi-
      CATS catsXi; //catsXi is an object belonging to CATS class
      catsXi.SetMomBins(NumMomBins, kMin, kMax); //number of bins, and limits
      catsXi.SetUseAnalyticSource(true);
      CATSparameters source_func(CATSparameters::tSource, 1, true); //compute the source function here, and then use it later
      catsXi.SetAnaSource(GaussSource, source_func);
      catsXi.SetMomentumDependentSource(false);
      catsXi.SetThetaDependentSource(false);
      catsXi.SetNumChannels(1);
      catsXi.SetQuantumStatistics(0);
      // catsXi.SetNumPW(0, 1);
      // catsXi.SetSpin(0, 0);
      catsXi.SetChannelWeight(0, 1.);
      catsXi.SetQ1Q2(chargecombi[i]); //product between the two charged particles
      catsXi.SetPdgId(321, (chargecombi[i]*3312)); //pdg of the particles: kaon and csi-
      catsXi.SetRedMass((massKaon * massXi) / (massKaon + massXi)); //reduced mass
      catsXi.SetAnaSource(0, sourceRadius);
      catsXi.KillTheCat();  

      auto grXiCoulomb = new TGraph();
      grXiCoulomb->SetTitle(";#it{k}* (MeV/#it{c}); #it{C}(#it{k}*)");

      auto grXiSmeared = new TGraph();
      grXiSmeared->SetTitle(";#it{k}* (MeV/#it{c}); #it{C}(#it{k}*)");


      FillCkGraph(catsXi, grXiCoulomb); //fill the XiCoulomb graph with catsXi values (in the form of a graph)
      grXiCoulomb->Write("Xi Coulomb before");

      TString CalibBaseDir = "/home/rossanafacen/Analysis/LambdaKaon/TestTask_nosph/Fit";
      auto calibFile = TFile::Open(TString::Format("%s/KXi_KLmabda7.root", CalibBaseDir.Data()));
      auto histDecayKinematicsXi = (TH2F*) calibFile->Get("KXi_KLambda"); //matrix for the feeddown  


      auto DLM_Xi = new DLM_Ck(1, 0, catsXi); //1: number of source patricles; 0: number of pot particles; cats object
      DLM_Xi->Update(); 
      const char *path_MC = "/home/rossanafacen/Analysis/LambdaKaon/TestTask_nosph/Classical/Trains/MC/";
      auto *file_MC = TFile::Open(TString::Format("%sAnalysisResults.root", path_MC));

      TList *dirResultsQA;
      
      //find and use the smear matrix, to compute the momentum resolution
      TDirectoryFile *dirQA = (TDirectoryFile*)(file_MC->FindObjectAny("HMResultsQA00"));
      dirQA->GetObject("HMResultsQA00", dirResultsQA);
      TList* PairQAList =  (TList*)dirResultsQA->FindObject("PairQA");
      TList* ParticlesQAList =  (TList*)PairQAList->FindObject(TString::Format("QA_Particle%d_Particle%d", i, i+2));
      auto Smear_Matrix = (TH2F*)ParticlesQAList->FindObject(TString::Format("MomentumResolutionME_Particle%d_Particle%d", i, i+2));   

      //transform to MeV the matrix of momentum smearing, taken from analysis results
      auto Smear_transf = (TH2F*)TransformToMeV(Smear_Matrix);
      //cut the matrix of smearing up to 1000 MeV
      auto Smear_transf_cut = (TH2F*)CutMeV(Smear_transf);  

      
      DLM_CkDecomposition CkDec_Xi("KXi_presmear", 0,
                                                  *DLM_Xi,
                                                  nullptr);  //class which allows you to do a decomposition. 
                                                  //title, numberofchilder, original Ck, matrix for the smearing
      
      TList *dirResults;
      TDirectoryFile *dir = (TDirectoryFile*)(file_MC->FindObjectAny("HMResults00"));
      dir->GetObject("HMResults00", dirResults);
      TList* ParticlesList =  (TList*)dirResults->FindObject(TString::Format("Particle%d_Particle%d", i, i+2));
      auto hME_GeV = (TH1F*)ParticlesList->FindObject(TString::Format("MEDist_Particle%d_Particle%d", i, i+2)); //mixed event distribution     
      auto hME = (TH1F*)TransformToMev1D(hME_GeV);
       
      outfile->cd();
      hME->Write("hME_MeV");
      hME_GeV->Write("hME_GeV");
      

      FillCkGraph(DLM_Xi, CkDec_Xi, grXiSmeared);
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

      DLM_CkDecomposition LKFullCF("LK", 2, *LKModel, Smear_transf_cut);//TODO enough to put mom smear here? or additionally also for contributions

      LKFullCF.AddContribution(0, lCsicharged, DLM_CkDecomposition::cFeedDown,
                                &CkDec_Xi, histDecayKinematicsXi); //first contribution has the lambda parameter from csi
     LKFullCF.AddContribution(1, lFlat, DLM_CkDecomposition::cFeedDown); //second contribution has the lambda parameter for flat contribution

      LKFullCF.AddPhaseSpace(hME); //
      LKFullCF.AddPhaseSpace(0,hME);

      LKFullCF.Update();

      FITTER_DECOMLedn = &LKFullCF;

      histDecayKinematicsXi->Write("decay Matrix");
      //total function, using model*background (we already did a prefit before)
      unsigned NumberOfFitPars = 8;            
      TF1* fitter = new TF1("fitter", fit_functionLedn, 0, 500, NumberOfFitPars);            
      fitter->SetParameter(0, scattlenRE);
      fitter->SetParameter(1, scattlenIM); 
      fitter->SetParLimits(1, 0, 100);           
      fitter->SetParameter(2, effrange);
      fitter->SetParLimits(2, 0, 100);                  
      fitter->SetParameter(3, fPreFit_pol2->GetParameter(0)); //constant of normalization           
      fitter->SetParameter(4, fPreFit_pol2->GetParameter(1)); //weight of common 
      fitter->SetParLimits(4, 0, 1);                  
      fitter->FixParameter(5, fPreFit_pol2->GetParameter(2)); //pol2  coefficients, to be kept fixed from the perfit    
      fitter->FixParameter(6, fPreFit_pol2->GetParameter(3));      
      fitter->FixParameter(7, fPreFit_pol2->GetParameter(4)); 
      fitter->SetParNames("Re(f0)","Im(f0)","d0", "N_d", "w_c", "a", "b", "c");
      
      cout << "\n---------------------------------\n Total fit pol2 \n--------------------------------\n" << endl;
      int status= h_HM->Fit(fitter, "MR", "", 0, 500);
      
      double chi_fitter = fitter->GetChisquare()/fitter->GetNDF();
      cout << "Chi reduced total fit " << chi_fitter << endl;
      
      Re_scatlen->Fill(fitter->GetParameter(0));
      Im_scatlen->Fill(fitter->GetParameter(1));
      effran->Fill(fitter->GetParameter(2));
    

      //we save y results, that s different for every bootstrap
      for (int i_save = 0; i_save < h_HM->GetNbinsX(); i_save++){
        Fit_mean[i_save] += fitter->Eval(h_HM->GetXaxis()->GetBinCenter(i_save));
        Fit_StDev[i_save] += (fitter->Eval(h_HM->GetXaxis()->GetBinCenter(i_save)) * fitter->Eval(h_HM->GetXaxis()->GetBinCenter(i_save)));  // N*rms^2
        entries[i_save]++; //counting
      }


      TGraph grGenuine;
      grGenuine.SetName("GenuineCF");
      
      TGraph grFeedCF;
      grFeedCF.SetName("GenuineCFSmearedLambda");
      grGenuine.SetLineColor(kRed);
      for (int i = 0; i < NumMomBins; ++i) {
        const float mom = LKModel->GetBinCenter(0, i);
        grGenuine.SetPoint(i, mom, LKFullCF.EvalMain(mom));  // genuine p-Phi CF with the parameters obtained in the fit
        grFeedCF.SetPoint(i, mom, LKFullCF.EvalMainFeed(mom));  // same as above, scaled by lambda params and momentum smearing
        //   grSidebandCF.SetPoint(i, mom, sidebandFullCF.EvalMain(mom));
      }

      outfile->cd();
      grGenuine.Write();
      grFeedCF.Write();

      TF1* fitter_bg = new TF1("fitter_bg", PreFit_pol2, 0, 2500, 5);
      fitter_bg->SetParameter(0, fitter->GetParameter(3));
      fitter_bg->SetParameter(1, fitter->GetParameter(4));
      fitter_bg->FixParameter(2, fitter->GetParameter(5));
      fitter_bg->FixParameter(3, fitter->GetParameter(6));
      fitter_bg->FixParameter(4, fitter->GetParameter(7));
      fitter_bg->SetParNames("N_d", "w_c", "a", "b", "c");

      fitter_bg->SetLineColor(fitter->GetLineColor());
      fitter_bg->SetLineStyle(6);

     /*  TCanvas *c_fitter = new TCanvas();
      h_HM->Draw("");
      fitter->Draw("same");
      fitter_bg->Draw("same");
      h_HM->GetXaxis()->SetRangeUser(0, 2500);
      h_HM->GetYaxis()->SetRangeUser(0.79, 1.06);
      h_HM->SetTitle("; k* (MeV/c); C(k*)");
      auto legend_fit = new TLegend(0.65, 0.65, 0.88, 0.88, NULL, "brNDC");
      legend_fit->AddEntry(h_HM, "Data", "LEP");
      legend_fit->AddEntry(fitter, "Fit Lednicky - pol2");
      legend_fit->AddEntry(fitter_bg, "Background");
      legend_fit->SetTextSize(gStyle->GetTextSize()*0.6); 
      legend_fit->Draw();
      c_fitter->Show();
      c_fitter->Print(TString::Format("/home/rossanafacen/Analysis/LambdaKaon/TestTask_nosph/Fit/%s/Lednic_fit_bg_%s.pdf", pair[i], pair[i]));
     */

      TF1* fitter_pol1 = new TF1("fitter_pol1", fit_functionLedn, 0, 500, NumberOfFitPars);            
      fitter_pol1->SetParameter(0, scattlenRE);
      fitter_pol1->SetParameter(1, scattlenIM);            
      fitter_pol1->SetParameter(2, effrange);        
      fitter_pol1->SetParameter(3, fPreFit_pol1->GetParameter(0)); //constant of normalization           
      fitter_pol1->SetParameter(4, fPreFit_pol1->GetParameter(1)); //weight of common             
      fitter_pol1->FixParameter(5, fPreFit_pol1->GetParameter(2)); //pol2  coefficients, to be kept fixed from the perfit    
      fitter_pol1->FixParameter(6, fPreFit_pol1->GetParameter(3));
      fitter_pol1->FixParameter(7, 0);      
      fitter_pol1->SetParNames("Re(f0)","Im(f0)","d0", "N_d", "w_c", "a", "b");
      fitter_pol1->SetLineColor(kGreen);  

      cout << "\n---------------------------------\n Total fit pol1 \n--------------------------------\n" << endl;
      h_HM->Fit(fitter_pol1, "MR", "", 0, 500);
      
      TGraph grGenuine_pol1;
      grGenuine_pol1.SetName("GenuineCF pol1");
      grGenuine_pol1.SetLineColor(kBlue);
      
      TGraph grFeedCF_pol1;
      grFeedCF_pol1.SetName("GenuineCFSmearedLambda pol1");
      grFeedCF_pol1.SetLineColor(kBlue + 5);
      
      for (int i = 0; i < NumMomBins; ++i) {
        const float mom = LKModel->GetBinCenter(0, i);
        grGenuine_pol1.SetPoint(i, mom, LKFullCF.EvalMain(mom));  // genuine p-Phi CF with the parameters obtained in the fit
        grFeedCF_pol1.SetPoint(i, mom, LKFullCF.EvalMainFeed(mom));  // same as above, scaled by lambda params and momentum smearing
      }

      grGenuine_pol1.Write();
      grFeedCF_pol1.Write();
      

      double chi_fitter_pol1 = fitter_pol1->GetChisquare()/fitter_pol1->GetNDF();
      cout << "Chi reduced total fit pol1: " << chi_fitter_pol1 << endl;
      
      outfile->cd();

      TF1* fitter_bg_pol1 = new TF1("fitter_bg_pol1", PreFit_pol1, 0, 2500, 5);
      fitter_bg_pol1->FixParameter(0, fitter_pol1->GetParameter(3));
      fitter_bg_pol1->FixParameter(1, fitter_pol1->GetParameter(4));
      fitter_bg_pol1->FixParameter(2, fitter_pol1->GetParameter(5));
      fitter_bg_pol1->FixParameter(3, fitter_pol1->GetParameter(6));
      //fitter_bg_pol1->FixParameter(4, fitter_pol1->GetParameter(7));
      
      cout << "\n---------------------------------\n Background of pol1 \n--------------------------------\n" << endl;
      fitter_bg_pol1->SetLineColor(fitter_pol1->GetLineColor());
      fitter_bg_pol1->SetLineStyle(6);

     /*  TCanvas *c_fitter_pol1 = new TCanvas();
      h_HM->Draw("");
      fitter_pol1->Draw("same");
      fitter_bg_pol1->Draw("same");
      //gfitterfull->Draw("same");
      h_HM->GetXaxis()->SetRangeUser(0, 2500);
      h_HM->GetYaxis()->SetRangeUser(0.79, 1.06);
      h_HM->SetTitle("; k* (MeV/c); C(k*)");
      auto legend_fit_pol1 = new TLegend(0.65, 0.65, 0.88, 0.88, NULL, "brNDC");
      legend_fit_pol1->AddEntry(h_HM, "Data", "LEP");
      legend_fit_pol1->AddEntry(fitter_pol1, "Fit Lednicky - pol1");
      legend_fit_pol1->AddEntry(fitter_bg_pol1, "Background");
      legend_fit_pol1->SetTextSize(gStyle->GetTextSize()*0.6); 
      legend_fit_pol1->Draw();
      c_fitter_pol1->Show();
      c_fitter_pol1->Print(TString::Format("/home/rossanafacen/Analysis/LambdaKaon/TestTask_nosph/Fit/%s/Lednic_fit_pol1_%s.pdf", pair[i], pair[i]));
     */

      Smear_Matrix->Write("Original Matrix");
      Smear_transf->Write("Matrix MeV transformed");
      Smear_transf_cut->Write("Matrix cut");

      /* TCanvas *c_matrix = new TCanvas();
      Smear_transf_cut->Draw("");
      c_matrix->Print(TString::Format("/home/rossanafacen/Analysis/LambdaKaon/TestTask_nosph/Fit/%s/Smear matrix.pdf", pair[i]));
 */

      fitter->Write("fitter");
      fitter_pol1->Write("fitter pol1");
      fitter_bg->Write("bacgkround");
      fitter_bg_pol1->Write("background pol1");
      grXiCoulomb->Write();
      grXiSmeared->Write();

      double scattlenRE_pb[2] = {-0.6, 0.27}; double scattlenRE_pb_err[] = {sqrt(pow(0.12,2) + pow(0.11,2)), sqrt(pow(0.12,2) +pow(0.07,2))};
      double scattlenIM_pb[] = {0.51, 0.40}; double  scattlenIM_pb_err[] = {sqrt(pow(0.15,2)+pow(0.12,2)), sqrt(pow(0.11,2)+pow(0.07,2))}; 
      double effrange_pb[] = {0.83, -5.23}; double effrange_pb_err[] = {sqrt(pow(0.47,2)+pow(1.23,2)), sqrt(pow(2.13,2)+pow(4.8,2))}; 


      //fit with Pb parameters
      TF1* fitter_pb = new TF1("fitter_pb", fit_functionLedn, 0, 500, NumberOfFitPars);            
      fitter_pb->FixParameter(0, scattlenRE_pb[i]);
      //fitter_pb->SetParLimits(0, scattlenRE_pb[i] - scattlenRE_pb_err[i]*2, scattlenRE_pb[i] + scattlenRE_pb_err[i]*2);
      fitter_pb->FixParameter(1, scattlenIM_pb[i]);            
      //fitter_pb->SetParLimits(1, scattlenIM_pb[i] - scattlenIM_pb_err[i]*2, scattlenIM_pb[i] + scattlenIM_pb_err[i]*2);
      fitter_pb->FixParameter(2, effrange_pb[i]);            
      //fitter_pb->SetParLimits(2, effrange_pb[i] - effrange_pb_err[i]*2, effrange_pb[i] + effrange_pb_err[i]*2);
      fitter_pb->SetParameter(3, fPreFit_pol2->GetParameter(0)); //constant of normalization           
      fitter_pb->SetParameter(4, fPreFit_pol2->GetParameter(1)); //weight of common             
      fitter_pb->FixParameter(5, fPreFit_pol2->GetParameter(2)); //pol2  coefficients, to be kept fixed from the perfit    
      fitter_pb->FixParameter(6, fPreFit_pol2->GetParameter(3));      
      fitter_pb->FixParameter(7, fPreFit_pol2->GetParameter(4));
      fitter_pb->SetParNames("Re(f0)_pb","Im(f0)_pb","d0_pb", "N_d_pb", "w_c_pb", "a", "b", "c");

      fitter_pb->SetLineColor(kAzure+7);

      h_HM->Fit(fitter_pb, "MR", "", 0, 500);
      double chi_fitter_pb = fitter_pb->GetChisquare()/fitter_pb->GetNDF();
      cout << "chi square for pb parameters: " << chi_fitter_pb << endl;

      TGraph grGenuine_pb;
      grGenuine_pb.SetName("GenuineCF pb");
      grGenuine_pb.SetLineColor(kGreen);
      
      TGraph grFeedCF_pb;
      grFeedCF_pb.SetName("GenuineCFSmearedLambda pb");
      grFeedCF_pb.SetLineColor(kViolet);
      
      for (int i = 0; i < NumMomBins; ++i) {
        const float mom = LKModel->GetBinCenter(0, i);
        grGenuine_pb.SetPoint(i, mom, LKFullCF.EvalMain(mom));  // genuine p-Phi CF with the parameters obtained in the fit
        grFeedCF_pb.SetPoint(i, mom, LKFullCF.EvalMainFeed(mom));  // same as above, scaled by lambda params and momentum smearing
      }

      grGenuine_pb.Write();
      grFeedCF_pb.Write();
      
      TGraph* gfitterfull_pb = new TGraph();
      gfitterfull_pb->SetLineColor(fitter_pb->GetLineColor());
      gfitterfull_pb->SetLineWidth(2);
      gfitterfull_pb->SetLineStyle(6);

      for (int x = 0; x < 2500; x++){
        double y = fitter_pb->Eval(x);
        gfitterfull_pb->SetPoint(x,x,y);
        }


      /* TCanvas *c_fitter_pb = new TCanvas();
      h_HM->Draw("");
      fitter_pb->Draw("same");
      h_HM->GetXaxis()->SetRangeUser(0, 2500);
      h_HM->GetYaxis()->SetRangeUser(0.87, 1.06);
      h_HM->SetTitle("; k* (MeV/c); C(k*)");
      auto legend_fit_pb = new TLegend(0.65, 0.65, 0.86, 0.86, NULL, "brNDC");
      legend_fit_pb->AddEntry(h_HM, "Data", "LEP");
      legend_fit_pb->AddEntry(fitter_pb, "Fit with Pb parms");
      legend_fit_pb->SetTextSize(gStyle->GetTextSize()*0.6); 
      legend_fit_pb->Draw();
      c_fitter_pb->Show();
      c_fitter_pb->Print(TString::Format("/home/rossanafacen/Analysis/LambdaKaon/TestTask_nosph/Fit/%s/Pb_fit_%s.pdf", pair[i], pair[i]));

 */
      fitter_pb->Write("fitter pb");
      gfitterfull_pb->Write("fitter full pb");
      
    }
        
  outfile->cd();      
  Re_scatlen->Write("Re Scatt length %d");
  Im_scatlen->Write("Im Scatt length");
  effran->Write("Eff range");  

  for (unsigned uBin = 0; uBin < h_HM_orig->GetNbinsX(); uBin++) {
    Fit_mean[uBin] /= double(entries[uBin]);   // mean of the band
    Fit_StDev[uBin] /= double(entries[uBin]);  // this is the rms^2
    Fit_StDev[uBin] = sqrt(Fit_StDev[uBin] - Fit_mean[uBin] * Fit_mean[uBin]);  // standard deviation
  }

  auto Cf_tot = new TH1F("Cf_tot", "Cf_tot", nbins, h_HM_orig->GetXaxis()->GetBinLowEdge(1), h_HM_orig->GetXaxis()->GetBinUpEdge(nbins));
  //auto Cf_tot = (TH1F *)h_HM_orig->Clone("Cf_tot");

  for (unsigned uBin = 0; uBin < nbins; uBin++) {    
    Cf_tot->SetBinContent(uBin + 1, Fit_mean[uBin]);
    Cf_tot->SetBinError(uBin + 1, Fit_StDev[uBin]);

    
  } 

  Cf_tot->Write("Lednicky fit");  
  outfile->Close();
  return 0;

}

