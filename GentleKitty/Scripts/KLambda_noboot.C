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
/// =====================================================================================
TH2F* GetMatix(const char* rootfile, const char* filename) {
  auto file = TFile::Open(rootfile);
  file->ls();
  TH2F* matrix = (TH2F*)file->FindObjectAny(Form("%s", filename));
  return matrix;
}
/// =====================================================================================
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

double weights(Double_t *x, Double_t *par){
  double t = *x;
  double fitval = par[0]*(par[1]*fTemplate_MC_Common->Eval(t) + (1. - par[1])*fTemplate_MC_NonCommon->Eval(t));
  return fitval;
}


double fit_NonCommon_mT(Double_t *x,Double_t *par) {
  double t = *x;
  double fitval =
  par[0]*(1.+par[1]*exp(-pow((t-par[2])/par[3],2.)))*(1.+par[4]*exp(-pow((t-par[5])/par[6],2.)))*(1.+par[7]*t+par[8]*t*t);
  return fitval;
}
/// =====================================================================================
double PreFit_pol1(Double_t *x,Double_t *par) {
  double t = *x;
  double fitval = 
  par[0]*(par[1]*fTemplate_MC_Common->Eval(t) + (1. - par[1])*fTemplate_MC_NonCommon->Eval(t) + (par[2]+par[3]*t));
  return fitval;
}
/// =====================================================================================
double PreFit_pol2(Double_t *x,Double_t *par) {
  double t = *x;
  double fitval = 
  par[0]*(par[1]*fTemplate_MC_Common->Eval(t) + (1. - par[1])*fTemplate_MC_NonCommon->Eval(t) + (par[2]+par[3]*t+par[4]*t*t));
  return fitval;
}
/// =====================================================================================
double PreFit_pol3(Double_t *x,Double_t *par) {
  double t = *x;
  double fitval = 
  par[0]*(par[1]*fTemplate_MC_Common->Eval(t) + (1. - par[1])*fTemplate_MC_NonCommon->Eval(t) + (par[2]+par[3]*t*t+par[4]*t*t*t));
  return fitval;
}
/// =====================================================================================
DLM_CkDecomposition* FITTER_DECOMLedn;


//total fit
double fit_functionLedn(double* x, double* par) {
    double t = *x;

  FITTER_DECOMLedn->GetCk()->SetPotPar(0, par[0]);//scatlen RE
  FITTER_DECOMLedn->GetCk()->SetPotPar(1, par[1]);//scattlen IM
  FITTER_DECOMLedn->GetCk()->SetPotPar(2, par[2]);//effrange
  FITTER_DECOMLedn->Update(true, true);

  double yval= par[3]*FITTER_DECOMLedn->EvalCk(t)*(par[4]*fTemplate_MC_Common->Eval(t) + (1. - par[4])*fTemplate_MC_NonCommon->Eval(t) + (par[5]+par[6]*t+par[7]*t*t));
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

void Scale(TH1F* histo, float scale) {
    std::cout << histo->GetNbinsX() << endl;
  for (int i = 0; i < histo->GetNbinsX(); i++) {
    float bincont = histo->GetBinContent(i);
    std::cout << i << " " << bincont << " scaled: " << (1 + (bincont - 1) * scale) << endl;
    histo->SetBinContent(i, (1 + (bincont - 1) * scale));
  }
}


int main(int argc, char* argv[]) {
  ///char* action1(argv[1]);
  const char* action1 = "pol2"; 
  std::cout << action1 << std::endl;
  
  //decide which bl we want
  bool pol1BL = false;
  bool pol3BL = false;
  bool pol2BL = false;
  bool polpb = false;

  if(action1=="pol3"){
    pol3BL=true;
  }
  else if(action1=="pol2"){
    pol2BL=true;
    std::cout<<"good"<<endl;
  }
  else if(action1 =="pb"){
    polpb=true;
  }
  else{
    std::cout<<"error"<<endl;
  }

  const char *pair[2] = {"LKPlus+ALKMin", "LKMin+ALKPlus"};
  const char* norm[2] = {"Norm 0.24-0.34", "Norm 0.50-0.80"};
  double norm1[4] = {0.24, 0.50, 0.34, 0.80}; 
  int i=0; //either LK+ or LK-

  //include the paths etc that we need
  const char* path_ancestors = "/home/rossanafacen/Analysis/LambdaKaon/TestTask_nosph/Common_Ancestors/Trains_MC";
  const char* path_CF = "/home/rossanafacen/Analysis/LambdaKaon/TestTask_nosph/Classical/Trains/HM";
  const char* path_results = "/home/rossanafacen/Analysis/LambdaKaon/TestTask_nosph/Fit/LKPlus+ALKMin/P-p method/Prefit";
  const char* finalfolder[2] = {"P-p method/Prefit", "P-p method/Everything open"};
  
  auto *file_CF = TFile::Open(TString::Format("%s/%s/CFOutput_%s_00_norm_%.2f-%.2f.root", path_CF, norm[i], pair[i], norm1[i], norm1[i+2])); //it opens the analysis results
  auto *h_HM = (TH1F *)file_CF->FindObjectAny("hCk_ReweightedMeV_1"); 
  
  //common and uncommon CF
  auto *file_common = TFile::Open(TString::Format("%s/CFOutput_%s_0_norm_%.2f-%.2f_Common.root", path_ancestors, pair[i], norm1[i], norm1[i+2]));
  auto *file_uncommon = TFile::Open(TString::Format("%s/CFOutput_%s_0_norm_%.2f-%.2f_Uncommon.root", path_ancestors, pair[i], norm1[i], norm1[i+2])); //it opens the analysis results
  
  //take a histogram from the root file of CF
  auto *h_MC_Common = (TH1F *)file_common->FindObjectAny("hCk_ReweightedMeV_1"); 
  auto *h_MC_NonCommon = (TH1F *)file_uncommon->FindObjectAny("hCk_ReweightedMeV_1");   

  //root file output
  auto outfile = new TFile(TString::Format("%s/Fit_%s_%s.root", path_results, action1, pair[i]), "RECREATE");   
  h_HM->Write("Raw Data");
  
  //take the matrix for the kinematic feeddown, considering KXi->KLambda
  TString CalibBaseDir = "/home/rossanafacen/Analysis/LambdaKaon/TestTask_nosph/Fit";
  auto calibFile = TFile::Open(TString::Format("%s/KXi_KLmabda7.root", CalibBaseDir.Data()));
  auto histDecayKinematicsXi = (TH2F*) calibFile->Get("KXi_KLambda"); 

  //take the smear matrix, to compute the momentum resolution
  const char *path_MC = "/home/rossanafacen/Analysis/LambdaKaon/TestTask_nosph/Classical/Trains/MC/";
  auto *file_MC = TFile::Open(TString::Format("%sAnalysisResults.root", path_MC));
  TList *dirResultsQA;  
  TDirectoryFile *dirQA = (TDirectoryFile*)(file_MC->FindObjectAny("HMResultsQA00"));
  dirQA->GetObject("HMResultsQA00", dirResultsQA);
  TList* PairQAList =  (TList*)dirResultsQA->FindObject("PairQA");
  TList* ParticlesQAList =  (TList*)PairQAList->FindObject(TString::Format("QA_Particle%d_Particle%d", i, i+2));
  auto Smear_Matrix = (TH2F*)ParticlesQAList->FindObject(TString::Format("MomentumResolutionME_Particle%d_Particle%d", i, i+2));   

  //transform to MeV the matrix of momentum smearing, taken from analysis results
  auto Smear_transf = (TH2F*)TransformToMeV(Smear_Matrix);
  //cut the matrix of smearing up to 1000 MeV
  auto Smear_transf_cut = (TH2F*)CutMeV(Smear_transf); 
  outfile->cd(); 
  Smear_transf->Write("Matrix MeV transformed");
  Smear_transf_cut->Write("Matrix cut");

  //mixed event distribution, for phase space contribution
  TList *dirResults;
  TDirectoryFile *dir = (TDirectoryFile*)(file_MC->FindObjectAny("HMResults00"));
  dir->GetObject("HMResults00", dirResults);
  TList* ParticlesList =  (TList*)dirResults->FindObject(TString::Format("Particle%d_Particle%d", i, i+2));
  auto hME_GeV = (TH1F*)ParticlesList->FindObject(TString::Format("MEDist_Particle%d_Particle%d", i, i+2));      
  auto hME = (TH1F*)TransformToMev1D(hME_GeV);

  //output file for scat parameters
  std::fstream output5;
  output5.open(Form("%s/Scatterin_pars_%s.dat", path_results, pair[i]), std::fstream::in | std::fstream::out | std::fstream::app);
  

  
  //mass of particles from dpg
  const auto pdgDatabase = TDatabasePDG::Instance();
  const double massLambda = pdgDatabase->GetParticle(3122)->Mass() * 1000; 
  const double massKaon = pdgDatabase->GetParticle(321)->Mass() * 1000; 
  const double massXi = pdgDatabase->GetParticle(3312)->Mass() * 1000; 
  const double masspiminus = pdgDatabase->GetParticle(211)->Mass() * 1000; 

  std::cout<<"mass lambda pdg: "<< massLambda <<endl;
  std::cout<<"mass kaon pdg: "<< massKaon <<endl;
  std::cout<<"mass xi- pdg: "<< massXi <<endl;
  std::cout<<"mass pi- pdg: "<< masspiminus <<endl;

  double chargecombi[2] = {-1, 1}; 
  
  //fractions and purity for lambda
  const double LambdaPurity = 0.940459;;
  const double Sigma_feeding = 0.192; //feeddown from the sigma
  const double Csicharged_feeding = 0.232/2; //feedown from the csi
  const double Csi0_feeding = 0.232/2; //feedown from the csi
  const double LambdaPrimary = 1.f - Sigma_feeding - Csi0_feeding-Csicharged_feeding; //primary fraction of the lambda

  const Particle Lambda(LambdaPurity, LambdaPrimary,
                        { { Csicharged_feeding, Csi0_feeding, Sigma_feeding} });


  //fractions and purity for kaon
  const double KaonPurity = 0.995;
  const double KaonPrimary = 0.998*0.94; // from Emma and Daniel
  const double KaonPhi = 0.998*0.06; //feeddown from phi 
  const double KaonSec = 1.f-KaonPrimary-KaonPhi;

  const Particle Kaon(KaonPurity, KaonPrimary,
      { { KaonSec, KaonPhi} });


  const CATSLambdaParam lambdaParam(Kaon, Lambda); //calculate the lambda parameters of kaon and lambda through cats

  lambdaParam.PrintLambdaParams();

  float lPrim = lambdaParam.GetLambdaParam(CATSLambdaParam::Primary); //take the genuine lambda
  float lSigma = lambdaParam.GetLambdaParam(CATSLambdaParam::Primary, 
                                                   CATSLambdaParam::FeedDown, 0, 2); //contribution from feeddown. 0 from k (secondary) and 2 from L (sigma feeding) 
  float lCsicharged = lambdaParam.GetLambdaParam(CATSLambdaParam::Primary, 
                                                   CATSLambdaParam::FeedDown, 0, 0);  
  float lCsi0 = lambdaParam.GetLambdaParam(CATSLambdaParam::Primary,
                                                   CATSLambdaParam::FeedDown, 0, 1);                                                                                                
  
  //print lambda parameters
  float lFlat = 1.f - lPrim - lCsicharged;
  std::cout << "---------------------------------------------------\n";
  std::cout << "---------------------------------------------------\n";
  std::cout << "Lambda parameters\n";
  std::cout << "Genuine K-Lambda " << lPrim << "\n";
  std::cout << "K-Xi- -> K-Lambda " << lCsicharged << "\n";
  std::cout << "Flat in K-Lambda " << lFlat << "\n";
  std::cout << "---------------------------------------------------\n";
  std::cout << "---------------------------------------------------\n";



  //Build a template of common and common correlation, to model the 2 contributions. I model them with a 3 gaussian
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
  
  fTemplate_MC_Common->SetLineColor(kRed + 2);  
  h_MC_Common->SetLineColor(kOrange + 2); 

  //perform the fit, excludig the first two bins
  h_MC_Common->Fit(fTemplate_MC_Common, "S, N, R, M", "", (h_MC_Common->GetXaxis()->GetBinCenter(2))+8, 2500);
    
  double chi_common = fTemplate_MC_Common->GetChisquare()/fTemplate_MC_Common->GetNDF();
  std::cout<<"Chi2 reduced common: "<< chi_common <<endl;

  outfile->cd();
  fTemplate_MC_Common->Write("Fit Common");
  h_MC_Common->Write("Data Common");
  

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
  
  fTemplate_MC_NonCommon->SetLineColor(kGreen + 2); 
  h_MC_NonCommon->SetLineColor(kViolet + 2); 
    
  h_MC_NonCommon->Fit(fTemplate_MC_NonCommon,"S, N, R, M", "", (h_MC_NonCommon->GetXaxis()->GetBinCenter(2))+8, 2500);
  
  double chi_noncommon = fTemplate_MC_NonCommon->GetChisquare()/fTemplate_MC_NonCommon->GetNDF();
  std::cout<<"Chi2 reduced noncommon: "<<chi_noncommon<<endl;
  
  h_MC_NonCommon->Write("Data NonCommon");
  fTemplate_MC_NonCommon->Write("Fit NonCommon");


  //prefit of the background (with pol1 and pol2): i have 4 parameters free, and I take the template from common and uncommom
  TF1* fPreFit_pol1 = new TF1("fPreFit_pol1", PreFit_pol1, 0, 2000, 4);
  fPreFit_pol1->SetParLimits(1, 0, 1);
  fPreFit_pol1->SetParameter(2, 0.93);
  fPreFit_pol1->SetParLimits(2, 0.5, 1);
  fPreFit_pol1->SetParameter(3, 2);
  fPreFit_pol1->SetParLimits(3, 0, 4);

  fPreFit_pol1->SetParNames("N_d", "w_c", "a", "b");
  fPreFit_pol1->SetLineColor(kMagenta);
  std::cout << "\n---------------------------------\n Prefit with pol1 \n--------------------------------\n" << endl;

  h_HM->Fit(fPreFit_pol1, "S, N, R, M", "", 400, 2000); 
  
  double chi_prefit_pol1 = fPreFit_pol1->GetChisquare()/fPreFit_pol1->GetNDF();
  std::cout<<"Chi2 reduced prefit pol1: "<< chi_prefit_pol1 <<endl;
 
  //same thing, but for pol2  
  TF1* fPreFit_pol2 = new TF1("fPreFit_pol2", PreFit_pol2, 0, 2000, 5);
  fPreFit_pol2->SetParameter(0, fPreFit_pol1->GetParameter(0));
  fPreFit_pol2->SetParameter(1, fPreFit_pol1->GetParameter(1));
  fPreFit_pol2->SetParLimits(1, 5e-2, 1);
  fPreFit_pol2->SetParameter(2, fPreFit_pol1->GetParameter(2));
  fPreFit_pol2->SetParameter(3, fPreFit_pol1->GetParameter(3));
  fPreFit_pol2->SetParNames("N_d", "w_c", "a", "b", "c");
  fPreFit_pol2->SetLineColor(kGreen + 2); 
  std::cout << "\n---------------------------------\n Prefit with pol2 \n--------------------------------\n" << endl;

  h_HM->Fit(fPreFit_pol2, "S, N, R, M", "", 400, 2000);
  
  double chi_prefit_pol2 = fPreFit_pol2->GetChisquare()/fPreFit_pol2->GetNDF();
  std::cout<<"Chi2 reduced prefit pol2: "<< chi_prefit_pol2 <<endl;

  outfile->cd();
  fPreFit_pol2->Write("prefit pol2");
  fPreFit_pol1->Write("prefit pol1");
  

  //starting using cats
  int NumMomBins=50;
  int kMin=0;
  int kMax=800;

  //initial scattering parameters
  double sourceRadius = 1.11251; // //L+K+(1.354):   low 1.07421 mean 1.11251 up 1.14903 //L+K-(1.355):   low 1.07366 mean 1.11193 up 1.14843
  double scattlenRE = 0.0;
  double scattlenIM = 0.0;
  double effrange = 0.0;

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
  catsXi.SetChannelWeight(0, 1.);
  catsXi.SetQ1Q2(chargecombi[i]); //product between the two charged particles
  catsXi.SetPdgId(321, (chargecombi[i]*3312)); //pdg of the particles: kaon and csi-
  catsXi.SetRedMass((massKaon * massXi) / (massKaon + massXi)); //reduced mass
  catsXi.SetAnaSource(0, sourceRadius);
  catsXi.KillTheCat();  


  auto grXiCoulomb = new TGraph();
  grXiCoulomb->SetTitle(";#it{k}* (MeV/#it{c}); #it{C}(#it{k}*)");
  FillCkGraph(catsXi, grXiCoulomb); //fill the XiCoulomb graph with catsXi values (in the form of a graph)
  grXiCoulomb->Write("Xi Coulomb before smearing");

  histDecayKinematicsXi->RebinX(2);
  histDecayKinematicsXi->RebinY(2);

  auto DLM_Xi = new DLM_Ck(1, 0, catsXi); //1: number of source patricles; 0: number of pot particles; cats object
  DLM_Xi->Update(); 
  
  DLM_CkDecomposition CkDec_Xi("KXi_presmear", 0,
                                              *DLM_Xi,
                                              nullptr);  //class which allows you to do a decomposition. 
                                              //title, numberofchilder, original Ck, matrix for the smearing
  
  DLM_CkDecomposition CkDec_Xi_Smeared("KXi_smear", 0,
                                              *DLM_Xi,
                                              histDecayKinematicsXi);  //smearing KXi for the kinematics feeddown

  CkDec_Xi_Smeared.AddPhaseSpace(hME); //
  
  auto grXiSmeared = new TGraph();
  grXiSmeared->SetTitle(";#it{k}* (MeV/#it{c}); #it{C}(#it{k}*)");
  FillCkGraph(DLM_Xi, CkDec_Xi_Smeared, grXiSmeared);
  
  auto grXi_before = new TGraph();
  grXi_before->SetTitle(";#it{k}* (MeV/#it{c}); #it{C}(#it{k}*)");
  FillCkGraph(DLM_Xi, grXi_before);
  
  outfile->cd();  
  histDecayKinematicsXi->Write("Decay matrix");
  grXiSmeared->Write("Xi Smeared");
  grXi_before->Write("Xi before");
  

  
  //LKModel
  DLM_Ck *LKModel = new DLM_Ck(1, 3, NumMomBins, kMin, kMax, ComplexLednicky_Singlet); 
  LKModel->SetSourcePar(0, sourceRadius);
  LKModel->SetPotPar(0, scattlenRE);
  LKModel->SetPotPar(1, scattlenIM);
  LKModel->SetPotPar(2, effrange);
  LKModel->Update(); 

  //momentum smearing
  DLM_CkDecomposition LKFullCF("LK", 2, *LKModel, Smear_transf_cut);
  //add contribution
  LKFullCF.AddContribution(0, lCsicharged, DLM_CkDecomposition::cFeedDown, &CkDec_Xi, histDecayKinematicsXi); //first contribution has the lambda parameter from csi
  LKFullCF.AddContribution(1, lFlat, DLM_CkDecomposition::cFeedDown); //second contribution has the lambda parameter for flat contribution

  LKFullCF.AddPhaseSpace(hME); 
  LKFullCF.AddPhaseSpace(0,hME);
  LKFullCF.Update();

  FITTER_DECOMLedn = &LKFullCF;
  unsigned NumberOfFitPars = 8;  
  
  if(pol2BL == true){  
    //total function, using model*background (we already did a prefit before)       
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
    
    std::cout << "\n---------------------------------\n Total fit pol2 \n--------------------------------\n" << endl;
    int status= h_HM->Fit(fitter, "MR", "", 0, 500);
    
    double chi_fitter = fitter->GetChisquare()/fitter->GetNDF();
    std::cout << "Chi reduced total fit " << chi_fitter << endl;
    std::cout<<"Re scat len after pol 2 fit"  <<LKFullCF.GetCk()->GetPotPar(0)<<endl;
    std::cout<<"Re scat len after pol 2 fit"  <<LKModel->GetPotPar(0)<<endl;

    //background
    TF1* fitter_bg = new TF1("fitter_bg", PreFit_pol2, 0, 2500, 5);
    fitter_bg->FixParameter(0, fitter->GetParameter(3));
    fitter_bg->FixParameter(1, fitter->GetParameter(4));
    fitter_bg->FixParameter(2, fitter->GetParameter(5));
    fitter_bg->FixParameter(3, fitter->GetParameter(6));
    fitter_bg->FixParameter(4, fitter->GetParameter(7));
    fitter_bg->SetParNames("N_d", "w_c", "a", "b", "c");
    fitter_bg->SetLineColor(fitter->GetLineColor());
    fitter_bg->SetLineStyle(6);
    fitter->Write("fitter");
    fitter_bg->Write("bacgkround");
    
    TF1* fitter_bg_pol = new TF1("fitter_bg_pol", "pol2", 0, 2500);
    fitter_bg_pol->FixParameter(0, fitter->GetParameter(5));
    fitter_bg_pol->FixParameter(1, fitter->GetParameter(6));
    fitter_bg_pol->FixParameter(2, fitter->GetParameter(7));
    fitter_bg_pol->Write("pol part");

    TF1* fitter_bg_comm = new TF1("fitter_bg_comm", weights, 0, 2500, 2);
    fitter_bg_comm->FixParameter(0, fitter->GetParameter(3));
    fitter_bg_comm->FixParameter(1, fitter->GetParameter(4));
    fitter_bg_comm->Write("weight part");

    //TF1* fitter_back = ;

    //fill the output file
    output5 << "With Pol2 \n";
    for(int i=0;i<NumberOfFitPars;i++){
      output5 << fitter->GetParName(i)<<" "<<fitter->GetParameter(i)<<" "<<fitter->GetParError(i)<< "\n";
    }

    //get the genuine correlation
    TGraph grGenuine;
    grGenuine.SetName("Genuine_CF_pol2");
    TGraph grFeedCF;
    grFeedCF.SetName("Genuine_CF_Smeared_pol2");
    grGenuine.SetLineColor(kRed);
    for (int i = 0; i < NumMomBins; ++i) {
      const float mom = LKModel->GetBinCenter(0, i);
      grGenuine.SetPoint(i, mom, LKFullCF.EvalMain(mom));  // genuine p-Phi CF with the parameters obtained in the fit
      grFeedCF.SetPoint(i, mom, LKFullCF.EvalMainFeed(mom));  // same as above, scaled by lambda params and momentum smearing
    }

    outfile->cd();
    grGenuine.Write("Genuine_CF_pol2");
    grFeedCF.Write("Genuine_CF_Smeared_pol2");

    
  }
  
  if(pol3BL == true){
      TF1* fitter_pol3 = new TF1("fitter_pol3", fit_functionLedn, 0, 500, NumberOfFitPars);            
      fitter_pol3->SetParameter(0, scattlenRE);
      fitter_pol3->SetParameter(1, scattlenIM);            
      fitter_pol3->SetParameter(2, effrange);        
      fitter_pol3->SetParameter(3, fPreFit_pol1->GetParameter(0)); //constant of normalization           
      fitter_pol3->SetParameter(4, fPreFit_pol1->GetParameter(1)); //weight of common             
      fitter_pol3->FixParameter(5, fPreFit_pol1->GetParameter(2)); //pol2  coefficients, to be kept fixed from the perfit    
      fitter_pol3->FixParameter(6, fPreFit_pol1->GetParameter(3));
      fitter_pol3->FixParameter(7, 0);      
      fitter_pol3->SetParNames("Re(f0)","Im(f0)","d0", "N_d", "w_c", "a", "b");
      fitter_pol3->SetLineColor(kGreen);  

      std::cout << "\n---------------------------------\n Total fit pol3 \n--------------------------------\n" << endl;
      h_HM->Fit(fitter_pol3, "MR", "", 0, 500);
      double chi_fitter_pol1 = fitter_pol3->GetChisquare()/fitter_pol3->GetNDF();
      std::cout << "Chi reduced total fit pol3: " << chi_fitter_pol1 << endl;
      
      TF1* fitter_bg_pol3 = new TF1("fitter_bg_pol3", PreFit_pol1, 0, 2500, 5);
      fitter_bg_pol3->FixParameter(0, fitter_pol3->GetParameter(3));
      fitter_bg_pol3->FixParameter(1, fitter_pol3->GetParameter(4));
      fitter_bg_pol3->FixParameter(2, fitter_pol3->GetParameter(5));
      fitter_bg_pol3->FixParameter(3, fitter_pol3->GetParameter(6));
      fitter_bg_pol3->SetLineColor(fitter_pol3->GetLineColor());
      fitter_bg_pol3->SetLineStyle(6);

      fitter_pol3->Write("pol1");
      fitter_bg_pol3->Write("background pol1");

      //write scattering parameters
      output5 << "--------------------------------\n";
      output5 << "With Pol1 \n";
      for(int i=0;i<NumberOfFitPars;i++){
        output5 << fitter_pol3->GetParName(i)<<" "<<fitter_pol3->GetParameter(i)<<" "<<fitter_pol3->GetParError(i)<< "\n";
      }
      
    }  
      
  
  if(pol1BL == true){
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

    std::cout << "\n---------------------------------\n Total fit pol1 \n--------------------------------\n" << endl;
    h_HM->Fit(fitter_pol1, "MR", "", 0, 500);
    double chi_fitter_pol1 = fitter_pol1->GetChisquare()/fitter_pol1->GetNDF();
    std::cout << "Chi reduced total fit pol1: " << chi_fitter_pol1 << endl;
    
    TF1* fitter_bg_pol1 = new TF1("fitter_bg_pol1", PreFit_pol1, 0, 2500, 5);
    fitter_bg_pol1->FixParameter(0, fitter_pol1->GetParameter(3));
    fitter_bg_pol1->FixParameter(1, fitter_pol1->GetParameter(4));
    fitter_bg_pol1->FixParameter(2, fitter_pol1->GetParameter(5));
    fitter_bg_pol1->FixParameter(3, fitter_pol1->GetParameter(6));
    fitter_bg_pol1->SetLineColor(fitter_pol1->GetLineColor());
    fitter_bg_pol1->SetLineStyle(6);

    fitter_pol1->Write("pol1");
    fitter_bg_pol1->Write("background pol1");

    //write scattering parameters
    output5 << "--------------------------------\n";
    output5 << "With Pol1 \n";
    for(int i=0;i<NumberOfFitPars;i++){
      output5 << fitter_pol1->GetParName(i)<<" "<<fitter_pol1->GetParameter(i)<<" "<<fitter_pol1->GetParError(i)<< "\n";
    }
    
  }  
  
  if(polpb == true){
    double scattlenRE_pb[2] = {-0.6, 0.27}; double scattlenRE_pb_err[] = {sqrt(pow(0.12,2) + pow(0.11,2)), sqrt(pow(0.12,2) +pow(0.07,2))};
    double scattlenIM_pb[] = {0.51, 0.40}; double  scattlenIM_pb_err[] = {sqrt(pow(0.15,2)+pow(0.12,2)), sqrt(pow(0.11,2)+pow(0.07,2))}; 
    double effrange_pb[] = {0.83, -5.23}; double effrange_pb_err[] = {sqrt(pow(0.47,2)+pow(1.23,2)), sqrt(pow(2.13,2)+pow(4.8,2))}; 

    std::cout << "\n---------------------------------\n Total fit pb \n--------------------------------\n" << endl;
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
    std::cout << "chi square for pb parameters: " << chi_fitter_pb << endl;

    //background  
    TF1* fitter_bg_pb = new TF1("fitter_bg_pb", PreFit_pol2, 0, 2500, 5);
    fitter_bg_pb->FixParameter(0, fitter_pb->GetParameter(3));
    fitter_bg_pb->FixParameter(1, fitter_pb->GetParameter(4));
    fitter_bg_pb->FixParameter(2, fitter_pb->GetParameter(5));
    fitter_bg_pb->FixParameter(3, fitter_pb->GetParameter(6));
    fitter_bg_pb->SetLineColor(fitter_pb->GetLineColor());
    fitter_bg_pb->SetLineStyle(6);

    outfile->cd();
    fitter_pb->Write("fitter pb");
    fitter_bg_pb->Write("fitter background pb");
    
    TGraph grGenuine_pb;
    TGraph grFeedCF_pb;
    grGenuine_pb.SetLineColor(kGreen);

    for (int i = 0; i < NumMomBins; ++i) {
      const float mom = LKModel->GetBinCenter(0, i);
      grGenuine_pb.SetPoint(i, mom, LKFullCF.EvalMain(mom));  // genuine p-Phi CF with the parameters obtained in the fit
      grFeedCF_pb.SetPoint(i, mom, LKFullCF.EvalMainFeed(mom));  // same as above, scaled by lambda params and momentum smearing
    }

    grGenuine_pb.Write("Genuine pb before");
    grFeedCF_pb.Write("Genuine pb after");
    
    output5 << "--------------------------------\n";
    output5 << "With Pb-Pb parameters \n";
    for(int i=0;i<NumberOfFitPars;i++){
      output5 << fitter_pb->GetParName(i)<<" "<<fitter_pb->GetParameter(i)<<" "<<fitter_pb->GetParError(i)<< "\n";
    }
  }
  

  output5.close();
  outfile->Close();
  return 0;

}


















/* DLM_Ck *LKModelResult = new DLM_Ck(1, 3, NumMomBins, kMin, kMax,
                                    ComplexLednicky_Singlet); //only scatt length?
    LKModelResult->SetSourcePar(0, sourceRadius);
    LKModelResult->SetPotPar(0, fitter->GetParameter(0));
    LKModelResult->SetPotPar(1, fitter->GetParameter(1));
    LKModelResult->SetPotPar(2, fitter->GetParameter(1));
    LKModelResult->Update();
    
    auto  grGenuine_new = new TGraph();
    FillCkGraph(LKModelResult, grGenuine_new);
    outfile->cd();
    grGenuine_new->Write("Genuine new");
 */