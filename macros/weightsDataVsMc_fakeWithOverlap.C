#define weightsDataVsMc_fakeWithOverlap_cxx

// ROOT includes 
#include <TROOT.h>
#include <TStyle.h>
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <iostream>

using namespace std;

// To compute weights mc->data for low-pt electrons matching PF after the fake selection

void weightsDataVsMc_fakeWithOverlap()
{
  // Files - from fakes selection (prepareInputsFromFakes) applied to data or MC with options isLowPt=-1 and studyOverlap=1
  TFile fileData("filesNew/pfOverlapStudy/myFileFakesInData__noWeights.root");  
  TFile fileMc("filesNew/pfOverlapStudy/myFileFakesInMc__dataMcWeightsForOverlap.root");  

  // Histos: pT
  TH1F *ptFakesWithOverlapMc   = (TH1F*)fileMc.Get("ptFakes_LptPfOverlap");
  TH1F *ptFakesWithOverlapMcWW = (TH1F*)fileMc.Get("ptFakes_LptPfOverlapWW");
  TH1F *ptFakesWithOverlapData = (TH1F*)fileData.Get("ptFakes_LptPfOverlap");
  ptFakesWithOverlapMc   -> Sumw2();   
  ptFakesWithOverlapMcWW -> Sumw2();   
  ptFakesWithOverlapData -> Sumw2(); 
  ptFakesWithOverlapMc   -> Scale(1./ptFakesWithOverlapMc->Integral());   
  ptFakesWithOverlapMcWW -> Scale(1./ptFakesWithOverlapMcWW->Integral());   
  ptFakesWithOverlapData -> Scale(1./ptFakesWithOverlapData->Integral());   
  //
  // Histos: eta
  TH1F *etaFakesWithOverlapMc   = (TH1F*)fileMc.Get("etaFakes_LptPfOverlap");
  TH1F *etaFakesWithOverlapMcWW = (TH1F*)fileMc.Get("etaFakes_LptPfOverlapWW");
  TH1F *etaFakesWithOverlapData = (TH1F*)fileData.Get("etaFakes_LptPfOverlap");
  etaFakesWithOverlapMc   -> Sumw2();   
  etaFakesWithOverlapMcWW -> Sumw2();   
  etaFakesWithOverlapData -> Sumw2(); 
  etaFakesWithOverlapMc   -> Scale(1./etaFakesWithOverlapMc->Integral());   
  etaFakesWithOverlapMcWW -> Scale(1./etaFakesWithOverlapMcWW->Integral());   
  etaFakesWithOverlapData -> Scale(1./etaFakesWithOverlapData->Integral());   
  //
  // Histos: mva
  TH1F *mvaFakesWithOverlapMc   = (TH1F*)fileMc.Get("mvaFakes_LptPfOverlap");
  TH1F *mvaFakesWithOverlapMcWW = (TH1F*)fileMc.Get("mvaFakes_LptPfOverlapWW");
  TH1F *mvaFakesWithOverlapData = (TH1F*)fileData.Get("mvaFakes_LptPfOverlap");
  mvaFakesWithOverlapMc   -> Sumw2();   
  mvaFakesWithOverlapMcWW -> Sumw2();   
  mvaFakesWithOverlapData -> Sumw2(); 
  mvaFakesWithOverlapMc   -> Scale(1./mvaFakesWithOverlapMc->Integral());   
  mvaFakesWithOverlapMcWW -> Scale(1./mvaFakesWithOverlapMcWW->Integral());   
  mvaFakesWithOverlapData -> Scale(1./mvaFakesWithOverlapData->Integral());   

  // 2Dim histos: pT vs eta 
  TH2F *ptVsEtaFakesWithOverlapMc   = (TH2F*)fileMc.Get("ptVsEtaFakes_LptPfOverlap");
  TH2F *ptVsEtaFakesWithOverlapMcWW = (TH2F*)fileMc.Get("ptVsEtaFakes_LptPfOverlapWW");
  TH2F *ptVsEtaFakesWithOverlapData = (TH2F*)fileData.Get("ptVsEtaFakes_LptPfOverlap");
  cout << " D " << endl;

  // --------------------------------------------------------
  // Compute weights
  TH2F *ptVsEtaWeights = (TH2F*)ptVsEtaFakesWithOverlapData->Clone("ptVsEtaWeights");
  ptVsEtaWeights->Sumw2(); 
  ptVsEtaWeights->Divide(ptVsEtaFakesWithOverlapMc);
  ptVsEtaWeights->SetTitle("ptVsEtaWeights");
  ptVsEtaWeights->SetName("ptVsEtaWeights");
  cout << "signal pre: min = " << ptVsEtaWeights->GetMinimum() << ", max = " << ptVsEtaWeights->GetMaximum() << endl;
  for (int iBinEta=0; iBinEta<(ptVsEtaWeights->GetNbinsX()); iBinEta++) {
    for (int iBinPt=0; iBinPt<(ptVsEtaWeights->GetNbinsY()); iBinPt++) {
      int iBinEtaP = iBinEta+1;
      int iBinPtP  = iBinPt+1;
      // if ((ptVsEtaWeights->GetBinContent(iBinEtaP,iBinPtP))>2) ptVsEtaWeights->SetBinContent(iBinEtaP,iBinPtP,2);
    }}
  cout << "signal post: min = " << ptVsEtaWeights->GetMinimum() << ", max = " << ptVsEtaWeights->GetMaximum() << endl;
  ptVsEtaWeights->Scale(1./ptVsEtaWeights->Integral()); 
  cout << "signal scaled: min = " << ptVsEtaWeights->GetMinimum() << ", max = " << ptVsEtaWeights->GetMaximum() << endl;

  // --------------------------------------------------------
  // Cosmetics
  ptFakesWithOverlapMc   -> SetLineWidth(2);  
  ptFakesWithOverlapMcWW -> SetLineWidth(2);  
  ptFakesWithOverlapData -> SetLineWidth(2);  
  ptFakesWithOverlapMc   -> SetLineColor(2);  
  ptFakesWithOverlapMcWW -> SetLineColor(3);  
  ptFakesWithOverlapData -> SetLineColor(4);  
  ptFakesWithOverlapMc   -> SetTitle("");
  ptFakesWithOverlapMcWW -> SetTitle("");
  ptFakesWithOverlapData -> SetTitle("");
  ptFakesWithOverlapMc   -> GetXaxis()->SetTitle("pT [GeV]");
  ptFakesWithOverlapMcWW -> GetXaxis()->SetTitle("pT [GeV]");
  ptFakesWithOverlapData -> GetXaxis()->SetTitle("pT [GeV]");
  //
  etaFakesWithOverlapMc   -> SetLineWidth(2);  
  etaFakesWithOverlapMcWW -> SetLineWidth(2);  
  etaFakesWithOverlapData -> SetLineWidth(2);  
  etaFakesWithOverlapMc   -> SetLineColor(2);  
  etaFakesWithOverlapMcWW -> SetLineColor(3);  
  etaFakesWithOverlapData -> SetLineColor(4);  
  etaFakesWithOverlapMc   -> SetTitle("");
  etaFakesWithOverlapMcWW -> SetTitle("");
  etaFakesWithOverlapData -> SetTitle("");
  etaFakesWithOverlapMc   -> GetXaxis()->SetTitle("#eta");
  etaFakesWithOverlapMcWW -> GetXaxis()->SetTitle("#eta");
  etaFakesWithOverlapData -> GetXaxis()->SetTitle("#eta");
  //
  mvaFakesWithOverlapMc   -> SetLineWidth(2);  
  mvaFakesWithOverlapMcWW -> SetLineWidth(2);  
  mvaFakesWithOverlapData -> SetLineWidth(2);  
  mvaFakesWithOverlapMc   -> SetLineColor(2);  
  mvaFakesWithOverlapMcWW -> SetLineColor(3);  
  mvaFakesWithOverlapData -> SetLineColor(4);  
  mvaFakesWithOverlapMc   -> SetTitle("");
  mvaFakesWithOverlapMcWW -> SetTitle("");
  mvaFakesWithOverlapData -> SetTitle("");
  mvaFakesWithOverlapMc   -> GetXaxis()->SetTitle("ID BDT output");
  mvaFakesWithOverlapMcWW -> GetXaxis()->SetTitle("ID BDT output");
  mvaFakesWithOverlapData -> GetXaxis()->SetTitle("ID BDT output");

  // -------------------------------------------------------------
  // Plots
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  // 
  TLegend *leg;
  leg = new TLegend(0.65,0.65,0.90,0.90);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.05);
  leg->SetFillColor(0);
  leg->AddEntry(mvaFakesWithOverlapMc,   "MC",   "lp");
  leg->AddEntry(mvaFakesWithOverlapData, "Data", "lp");
  //
  TLegend *legww;
  legww = new TLegend(0.60,0.65,0.80,0.90);
  legww->SetFillStyle(0);
  legww->SetBorderSize(0);
  legww->SetTextSize(0.05);
  legww->SetFillColor(0);
  legww->AddEntry(mvaFakesWithOverlapMcWW, "MC, weighted", "lp");
  legww->AddEntry(mvaFakesWithOverlapData, "Data", "lp");

  TLegend *leg2;
  leg2 = new TLegend(0.15,0.65,0.40,0.90);
  leg2->SetFillStyle(0);
  leg2->SetBorderSize(0);
  leg2->SetTextSize(0.05);
  leg2->SetFillColor(0);
  leg2->AddEntry(mvaFakesWithOverlapMc,   "MC", "lp");
  leg2->AddEntry(mvaFakesWithOverlapData, "Data", "lp");
  //
  TLegend *leg2ww;
  leg2ww = new TLegend(0.15,0.65,0.40,0.90);
  leg2ww->SetFillStyle(0);
  leg2ww->SetBorderSize(0);
  leg2ww->SetTextSize(0.05);
  leg2ww->SetFillColor(0);
  leg2ww->AddEntry(mvaFakesWithOverlapMcWW, "MC, weighted", "lp");
  leg2ww->AddEntry(mvaFakesWithOverlapData, "Data", "lp");

  TLegend *leg3;
  leg3 = new TLegend(0.70,0.70,0.90,0.90);
  leg3->SetFillStyle(0);
  leg3->SetBorderSize(0);
  leg3->SetTextSize(0.05);
  leg3->SetFillColor(0);
  leg3->AddEntry(mvaFakesWithOverlapMc,   "MC", "lp");
  leg3->AddEntry(mvaFakesWithOverlapData, "Data", "lp");
  //
  TLegend *leg3ww;
  leg3ww = new TLegend(0.55,0.70,0.85,0.90);
  leg3ww->SetFillStyle(0);
  leg3ww->SetBorderSize(0);
  leg3ww->SetTextSize(0.05);
  leg3ww->SetFillColor(0);
  leg3ww->AddEntry(mvaFakesWithOverlapMcWW, "MC, weighted", "lp");
  leg3ww->AddEntry(mvaFakesWithOverlapData, "Data", "lp");

  TCanvas cs1("cs1","cs1",1);
  ptFakesWithOverlapData->Rebin(4); // 6
  ptFakesWithOverlapMc->Rebin(4);
  ptFakesWithOverlapData->DrawNormalized("histe");
  ptFakesWithOverlapMc->DrawNormalized("samehiste");
  leg->Draw();
  cs1.SaveAs("pt_withOverlap_dataVsMc.png");

  TCanvas cs2("cs2","cs2",1);
  ptFakesWithOverlapMcWW->Rebin(4);
  ptFakesWithOverlapMcWW->DrawNormalized("histe");
  ptFakesWithOverlapData->DrawNormalized("samehiste");
  legww->Draw();
  cs2.SaveAs("pt_withOverlap_dataVsMc_withWeights.png");
  
  // -----------

  TCanvas cs3("cs3","cs3",1);
  etaFakesWithOverlapData->Rebin(4);
  etaFakesWithOverlapMc->Rebin(4);
  etaFakesWithOverlapData->DrawNormalized("histe");
  etaFakesWithOverlapMc->DrawNormalized("samehiste");
  leg3->Draw();
  cs3.SaveAs("eta_withOverlap_dataVsMc.png");

  TCanvas cs4("cs4","cs4",1);
  etaFakesWithOverlapMcWW->Rebin(4);
  etaFakesWithOverlapMcWW->DrawNormalized("histe");
  etaFakesWithOverlapData->DrawNormalized("samehiste");
  leg3ww->Draw();
  cs4.SaveAs("eta_withOverlap_dataVsMc_withWeights.png");

  // -----------

  TCanvas cs5("cs5","cs5",1);
  mvaFakesWithOverlapData->Rebin(4); // 6
  mvaFakesWithOverlapMc->Rebin(4);
  mvaFakesWithOverlapData->DrawNormalized("histe");
  mvaFakesWithOverlapMc->DrawNormalized("samehiste");
  leg2->Draw();
  cs5.SaveAs("mva_withOverlap_dataVsMc.png");

  TCanvas cs6("cs6","cs6",1);
  mvaFakesWithOverlapMcWW->Rebin(4);
  mvaFakesWithOverlapMcWW->DrawNormalized("histe");
  mvaFakesWithOverlapData->DrawNormalized("samehiste");
  leg2ww->Draw();
  cs6.SaveAs("mva_withOverlap_dataVsMc_withWeights.png");

  // -----------

  TCanvas csw("csw","cs",1);
  ptVsEtaWeights->Draw("colz");
  csw.SaveAs("ptVsEtaWeights.png"); 

  // -----------

  // Check weights 2dim
  TH2F *check2dSignal = (TH2F*)ptVsEtaFakesWithOverlapMcWW->Clone("check2dSignal");  
  check2dSignal->Sumw2();
  check2dSignal->Divide(ptVsEtaFakesWithOverlapData);

  // ----------------------------------------------  
  TFile myFile("weightFile_withPFoverlap_datiVsMc.root","RECREATE");
  ptVsEtaWeights->Write();
  check2dSignal->Write();
  myFile.Close();
}
