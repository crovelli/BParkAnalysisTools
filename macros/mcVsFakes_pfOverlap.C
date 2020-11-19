#define mcVsFakes_cxx

// ROOT includes 
#include <TROOT.h>
#include <TStyle.h>
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <iostream>

using namespace std;

void mcVsFakes_pfOverlap()
{
  // ----------------------------------------------------     
  // Files: 
  TFile fileFakesData("filesNew/pfOverlapStudy/myFileFakesInData__noWeights.root");     // from fakes selection applied to data (prepareInputsFromFakes, options isLowPt=-1 and studyOverlap=1)
  TFile fileFakesMc("filesNew/pfOverlapStudy/myFileFakesInMc__dataMcWeightsForOverlap.root");     // from fakes selection applied to (B->KJPsi->Kmm) MC (prepareInputsFromFakes, options isLowPt=-1 and studyOverlap=1)

  // ----------------------------------------------------

  // Data histos, fake selection
  TH1F *mvaFakesNoOverlapData = (TH1F*)fileFakesData.Get("mvaFakes_LptNotPfOverlap");
  TH1F *ptFakesNoOverlapData  = (TH1F*)fileFakesData.Get("ptFakes_LptNotPfOverlap");
  TH1F *etaFakesNoOverlapData = (TH1F*)fileFakesData.Get("etaFakes_LptNotPfOverlap");
  mvaFakesNoOverlapData->Sumw2();   
  ptFakesNoOverlapData->Sumw2();   
  etaFakesNoOverlapData->Sumw2();   
  //
  TH1F *mvaFakesWithOverlapData = (TH1F*)fileFakesData.Get("mvaFakes_LptPfOverlap");
  TH1F *ptFakesWithOverlapData  = (TH1F*)fileFakesData.Get("ptFakes_LptPfOverlap");
  TH1F *etaFakesWithOverlapData = (TH1F*)fileFakesData.Get("etaFakes_LptPfOverlap");
  mvaFakesWithOverlapData->Sumw2();   
  ptFakesWithOverlapData->Sumw2();   
  etaFakesWithOverlapData->Sumw2();   

  // MC histos, fake selection
  TH1F *mvaFakesNoOverlapMc = (TH1F*)fileFakesMc.Get("mvaFakes_LptNotPfOverlap");
  TH1F *ptFakesNoOverlapMc  = (TH1F*)fileFakesMc.Get("ptFakes_LptNotPfOverlap");
  TH1F *etaFakesNoOverlapMc = (TH1F*)fileFakesMc.Get("etaFakes_LptNotPfOverlap");
  mvaFakesNoOverlapMc->Sumw2();   
  ptFakesNoOverlapMc->Sumw2();   
  etaFakesNoOverlapMc->Sumw2();   
  //
  TH1F *mvaFakesWithOverlapMc = (TH1F*)fileFakesMc.Get("mvaFakes_LptPfOverlap");
  TH1F *ptFakesWithOverlapMc  = (TH1F*)fileFakesMc.Get("ptFakes_LptPfOverlap");
  TH1F *etaFakesWithOverlapMc = (TH1F*)fileFakesMc.Get("etaFakes_LptPfOverlap");
  mvaFakesWithOverlapMc->Sumw2();   
  ptFakesWithOverlapMc->Sumw2();   
  etaFakesWithOverlapMc->Sumw2();   

  // ----------------------------------------------------

  // Cosmetics
  mvaFakesWithOverlapData -> SetLineWidth(2);  
  mvaFakesWithOverlapData -> SetLineColor(3);  
  ptFakesWithOverlapData  -> SetLineWidth(2);  
  ptFakesWithOverlapData  -> SetLineColor(3);  
  etaFakesWithOverlapData -> SetLineWidth(2);  
  etaFakesWithOverlapData -> SetLineColor(3);  
  //
  mvaFakesNoOverlapData -> SetLineWidth(2);  
  mvaFakesNoOverlapData -> SetLineColor(3);  
  ptFakesNoOverlapData  -> SetLineWidth(2);  
  ptFakesNoOverlapData  -> SetLineColor(3);  
  etaFakesNoOverlapData -> SetLineWidth(2);  
  etaFakesNoOverlapData -> SetLineColor(3);  
  //
  mvaFakesNoOverlapMc -> SetLineWidth(2);  
  mvaFakesNoOverlapMc -> SetLineColor(6);  
  ptFakesNoOverlapMc  -> SetLineWidth(2);  
  ptFakesNoOverlapMc  -> SetLineColor(6);  
  etaFakesNoOverlapMc -> SetLineWidth(2);  
  etaFakesNoOverlapMc -> SetLineColor(6);  
  //
  mvaFakesWithOverlapMc -> SetLineWidth(2);  
  mvaFakesWithOverlapMc -> SetLineColor(6);  
  ptFakesWithOverlapMc  -> SetLineWidth(2);  
  ptFakesWithOverlapMc  -> SetLineColor(6);  
  etaFakesWithOverlapMc -> SetLineWidth(2);  
  etaFakesWithOverlapMc -> SetLineColor(6);  

  // ----------------------------------------------------  
  // Plots: with fake selection: data vs Mc 
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);

  TLegend *legA;
  legA = new TLegend(0.55,0.65,0.85,0.90);
  legA->SetFillStyle(0);
  legA->SetBorderSize(0);
  legA->SetTextSize(0.05);
  legA->SetFillColor(0);
  legA->AddEntry(mvaFakesNoOverlapData, "Fake sel., data", "lp");
  legA->AddEntry(mvaFakesNoOverlapMc,   "Fake sel., MC", "lp");
  //
  TLegend *legA2;
  legA2 = new TLegend(0.30,0.10,0.70,0.40);
  legA2->SetFillStyle(0);
  legA2->SetBorderSize(0);
  legA2->SetTextSize(0.05);
  legA2->SetFillColor(0);
  legA2->AddEntry(mvaFakesNoOverlapData, "Fake sel., data", "lp");
  legA2->AddEntry(mvaFakesNoOverlapMc,   "Fake sel., MC", "lp");

  TCanvas cmvaA("cmvaA","cmvaA",1);
  mvaFakesNoOverlapData->SetTitle("");
  mvaFakesNoOverlapMc->SetTitle("");
  mvaFakesNoOverlapData->GetXaxis()->SetTitle("Id BDT");
  mvaFakesNoOverlapMc->GetXaxis()->SetTitle("Id BDT");
  mvaFakesNoOverlapMc->DrawNormalized("histE");
  mvaFakesNoOverlapData->DrawNormalized("samehistE");
  legA->Draw();
  cmvaA.SaveAs("outputBDT_noOverlap_dataVsMc.png");
  //
  TCanvas cmvaA1("cmvaA1","cmvaA1",1); 
  mvaFakesWithOverlapData->Rebin(4);
  mvaFakesWithOverlapMc->Rebin(4);
 mvaFakesWithOverlapData->SetTitle("");
  mvaFakesWithOverlapMc->SetTitle("");
  mvaFakesWithOverlapData->GetXaxis()->SetTitle("Id BDT");
  mvaFakesWithOverlapMc->GetXaxis()->SetTitle("Id BDT");
  mvaFakesWithOverlapMc->DrawNormalized("histE");
  mvaFakesWithOverlapData->DrawNormalized("samehistE");
  legA->Draw();
  cmvaA1.SaveAs("outputBDT_withOverlap_dataVsMc.png");

  TCanvas cptA("cptA","cptA",1);
  ptFakesNoOverlapData->SetTitle("");
  ptFakesNoOverlapMc->SetTitle("");
  ptFakesNoOverlapData->GetXaxis()->SetTitle("p_{T} [GeV]");
  ptFakesNoOverlapMc->GetXaxis()->SetTitle("p_{T} [GeV]");
  ptFakesNoOverlapMc->DrawNormalized("histE");
  ptFakesNoOverlapData->DrawNormalized("samehistE");
  legA->Draw();
  cptA.SaveAs("pt_noOverlap_dataVsMc.png");
  //
  TCanvas cptA1("cptA1","cptA1",1);
  ptFakesWithOverlapData->Rebin(4);
  ptFakesWithOverlapMc->Rebin(4);
  ptFakesWithOverlapData->SetTitle("");
  ptFakesWithOverlapMc->SetTitle("");
  ptFakesWithOverlapData->GetXaxis()->SetTitle("p_{T} [GeV]");
  ptFakesWithOverlapMc->GetXaxis()->SetTitle("p_{T} [GeV]");
  ptFakesWithOverlapData->DrawNormalized("histE");
  ptFakesWithOverlapMc->DrawNormalized("samehistE");
  legA->Draw();
  cptA1.SaveAs("pt_withOverlap_dataVsMc.png");

  TCanvas cetaA("cetaA","cetaA",1);
  etaFakesNoOverlapData->SetTitle("");
  etaFakesNoOverlapMc->SetTitle("");
  etaFakesNoOverlapData->GetXaxis()->SetTitle("#eta");
  etaFakesNoOverlapMc->GetXaxis()->SetTitle("#eta");
  etaFakesNoOverlapMc->DrawNormalized("histE");
  etaFakesNoOverlapData->DrawNormalized("samehistE");
  legA->Draw();
  cetaA.SaveAs("eta_noOverlap_dataVsMc.png");
  //
  TCanvas cetaA1("cetaA1","cetaA1",1);
  etaFakesWithOverlapData->Rebin(4);
  etaFakesWithOverlapMc->Rebin(4);
  etaFakesWithOverlapData->SetTitle("");
  etaFakesWithOverlapMc->SetTitle("");
  etaFakesWithOverlapData->GetXaxis()->SetTitle("#eta");
  etaFakesWithOverlapMc->GetXaxis()->SetTitle("#eta");
  etaFakesWithOverlapMc->DrawNormalized("histE");
  etaFakesWithOverlapData->DrawNormalized("samehistE");
  legA->Draw();
  cetaA1.SaveAs("eta_withOverlap_dataVsMc.png");

}
