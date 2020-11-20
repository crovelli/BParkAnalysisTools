#define mcVsSplots_cxx

// Data vs MC signal distributions

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

void mcVsSplots_pfOverlap()
{
  // Input files (after TnP selection)
  TFile *fileMC = new TFile("filesNew/pfOverlapStudy/myFileMcAfterTnp.root");                                                   // from prepareInputsFromMcWithTnP 
  TFile *fileSPlotsWithOverlap = new TFile("filesNew/pfOverlapStudy/myFileSPlots_lowPtWithPfOverlap_inclusive.root");           // from jpsi_splot.C 
  TFile *fileSPlotsNoOverlap   = new TFile("filesNew/pfOverlapStudy/myFileSPlots_lowPtWithNoPfOverlap_inclusive.root");         // from jpsi_splot.C 

  // MC histos
  TH1F *probeMvaSignalMc_noOverlap = (TH1F*)fileMC->Get("probeMvaSignalMc_LptNotPfOverlap"); 
  TH1F *probePtSignalMc_noOverlap  = (TH1F*)fileMC->Get("probePtSignalMc_LptNotPfOverlap");
  TH1F *probeEtaSignalMc_noOverlap = (TH1F*)fileMC->Get("probeEtaSignalMc_LptNotPfOverlap");
  probeMvaSignalMc_noOverlap->Sumw2();   
  probePtSignalMc_noOverlap->Sumw2(); 
  probeEtaSignalMc_noOverlap->Sumw2(); 

  TH1F *probeMvaSignalMc_withOverlap = (TH1F*)fileMC->Get("probeMvaSignalMc_LptPfOverlap"); 
  TH1F *probePtSignalMc_withOverlap  = (TH1F*)fileMC->Get("probePtSignalMc_LptPfOverlap");
  TH1F *probeEtaSignalMc_withOverlap = (TH1F*)fileMC->Get("probeEtaSignalMc_LptPfOverlap");
  probeMvaSignalMc_withOverlap->Sumw2();   
  probePtSignalMc_withOverlap->Sumw2(); 
  probeEtaSignalMc_withOverlap->Sumw2(); 

  // s-Plots output
  TH1F *probeMvaSignalData_noOverlap = (TH1F*)fileSPlotsNoOverlap->Get("h1_probeMvaId_jpsi__probeMvaId");
  TH1F *probePtSignalData_noOverlap  = (TH1F*)fileSPlotsNoOverlap->Get("h1_probePt_jpsi__probePt");
  TH1F *probeEtaSignalData_noOverlap = (TH1F*)fileSPlotsNoOverlap->Get("h1_probeEta_jpsi__probeEta");
  probeMvaSignalData_noOverlap->Sumw2();
  probePtSignalData_noOverlap->Sumw2();
  probeEtaSignalData_noOverlap->Sumw2();

  TH1F *probeMvaSignalData_withOverlap = (TH1F*)fileSPlotsWithOverlap->Get("h1_probeMvaId_jpsi__probeMvaId");
  TH1F *probePtSignalData_withOverlap  = (TH1F*)fileSPlotsWithOverlap->Get("h1_probePt_jpsi__probePt");
  TH1F *probeEtaSignalData_withOverlap = (TH1F*)fileSPlotsWithOverlap->Get("h1_probeEta_jpsi__probeEta");
  probeMvaSignalData_withOverlap->Sumw2();
  probePtSignalData_withOverlap->Sumw2();
  probeEtaSignalData_withOverlap->Sumw2();

  // Cosmetics
  probeMvaSignalMc_withOverlap -> SetLineWidth(2);  
  probeMvaSignalMc_withOverlap -> SetLineColor(6);  
  probePtSignalMc_withOverlap  -> SetLineWidth(2);  
  probePtSignalMc_withOverlap  -> SetLineColor(6);  
  probeEtaSignalMc_withOverlap -> SetLineWidth(2);  
  probeEtaSignalMc_withOverlap -> SetLineColor(6);  
  // 
  probeMvaSignalMc_noOverlap -> SetLineWidth(2);  
  probeMvaSignalMc_noOverlap -> SetLineColor(4);  
  probePtSignalMc_noOverlap  -> SetLineWidth(2);  
  probePtSignalMc_noOverlap  -> SetLineColor(4);  
  probeEtaSignalMc_noOverlap -> SetLineWidth(2);  
  probeEtaSignalMc_noOverlap -> SetLineColor(4);  
  //
  probeMvaSignalData_withOverlap -> SetLineWidth(2);  
  probeMvaSignalData_withOverlap -> SetLineColor(6);  
  probePtSignalData_withOverlap  -> SetLineWidth(2);  
  probePtSignalData_withOverlap  -> SetLineColor(6);  
  probeEtaSignalData_withOverlap -> SetLineWidth(2);  
  probeEtaSignalData_withOverlap -> SetLineColor(6);  
  // 
  probeMvaSignalData_noOverlap -> SetLineWidth(2);  
  probeMvaSignalData_noOverlap -> SetLineColor(4);  
  probePtSignalData_noOverlap  -> SetLineWidth(2);  
  probePtSignalData_noOverlap  -> SetLineColor(4);  
  probeEtaSignalData_noOverlap -> SetLineWidth(2);  
  probeEtaSignalData_noOverlap -> SetLineColor(4);  


  // Rebin
  probeEtaSignalMc_noOverlap->Rebin(2);
  probeEtaSignalData_noOverlap->Rebin(2);
  probeMvaSignalMc_noOverlap->Rebin(2);
  probeMvaSignalData_noOverlap->Rebin(2);
  probePtSignalMc_noOverlap->Rebin(2);
  probePtSignalData_noOverlap->Rebin(2);
  probeEtaSignalMc_withOverlap->Rebin(2);
  probeEtaSignalData_withOverlap->Rebin(2);
  probeMvaSignalMc_withOverlap->Rebin(2);
  probeMvaSignalData_withOverlap->Rebin(2);
  probePtSignalMc_withOverlap->Rebin(2);
  probePtSignalData_withOverlap->Rebin(2);


  // -------------------------------------------------------------
  // Plots: with vs without overlap
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);

  TLegend *legA1;
  legA1 = new TLegend(0.10,0.65,0.45,0.90);
  legA1->SetFillStyle(0);
  legA1->SetBorderSize(0);
  legA1->SetTextSize(0.05);
  legA1->SetFillColor(0);
  legA1->AddEntry(probeMvaSignalData_noOverlap,   "No PF overlap", "lp");
  legA1->AddEntry(probeMvaSignalData_withOverlap, "PF overlap", "lp");
  //
  TLegend *legB1;
  legB1 = new TLegend(0.40,0.65,0.75,0.90);
  legB1->SetFillStyle(0);
  legB1->SetBorderSize(0);
  legB1->SetTextSize(0.05);
  legB1->SetFillColor(0);
  legB1->AddEntry(probeMvaSignalData_noOverlap,   "No PF overlap", "lp");
  legB1->AddEntry(probeMvaSignalData_withOverlap, "PF overlap", "lp");
  //
  TLegend *legC1;
  legC1 = new TLegend(0.25,0.15,0.65,0.30);
  legC1->SetFillStyle(0);
  legC1->SetBorderSize(0);
  legC1->SetTextSize(0.05);
  legC1->SetFillColor(0);
  legC1->AddEntry(probeMvaSignalData_noOverlap,   "No PF overlap", "lp");
  legC1->AddEntry(probeMvaSignalData_withOverlap, "PF overlap", "lp");
  
  TCanvas ca1("ca1","ca1",1);
  probeMvaSignalMc_withOverlap -> SetTitle(""); 
  probeMvaSignalMc_noOverlap   -> SetTitle(""); 
  probeMvaSignalMc_withOverlap -> GetXaxis()->SetTitle("Id BDT");
  probeMvaSignalMc_noOverlap   -> GetXaxis()->SetTitle("Id BDT");
  probeMvaSignalMc_withOverlap -> DrawNormalized("histe");
  probeMvaSignalMc_noOverlap   -> DrawNormalized("samehiste");
  legA1->Draw();
  ca1.SaveAs("outputBDTMc_w-wo-overlap.png");
  //
  TCanvas ca2("ca2","ca2",1);
  probePtSignalMc_withOverlap -> SetTitle(""); 
  probePtSignalMc_noOverlap   -> SetTitle(""); 
  probePtSignalMc_withOverlap -> GetXaxis()->SetTitle("p_{T} [GeV]");
  probePtSignalMc_noOverlap   -> GetXaxis()->SetTitle("p_{T} [GeV]");
  probePtSignalMc_noOverlap -> DrawNormalized("histe");
  probePtSignalMc_withOverlap   -> DrawNormalized("samehiste");
  legB1->Draw();
  ca2.SaveAs("ptMc_w-wo-overlap.png");
  //
  TCanvas ca3("ca3","ca3",1);
  probeEtaSignalMc_withOverlap -> SetTitle(""); 
  probeEtaSignalMc_noOverlap   -> SetTitle(""); 
  probeEtaSignalMc_withOverlap -> GetXaxis()->SetTitle("#eta");
  probeEtaSignalMc_noOverlap   -> GetXaxis()->SetTitle("#eta");
  probeEtaSignalMc_withOverlap -> DrawNormalized("histe");
  probeEtaSignalMc_noOverlap   -> DrawNormalized("samehiste");
  legC1->Draw();
  ca3.SaveAs("etaMc_w-wo-overlap.png");
  //
  //
  //
  TCanvas cb1("cb1","cb1",1);
  probeMvaSignalData_withOverlap -> SetTitle(""); 
  probeMvaSignalData_noOverlap   -> SetTitle(""); 
  probeMvaSignalData_withOverlap -> GetXaxis()->SetTitle("Id BDT");
  probeMvaSignalData_noOverlap   -> GetXaxis()->SetTitle("Id BDT");
  probeMvaSignalData_withOverlap -> DrawNormalized("histe");
  probeMvaSignalData_noOverlap   -> DrawNormalized("samehiste");
  legA1->Draw();
  cb1.SaveAs("outputBDTData_w-wo-overlap.png");
  //
  TCanvas cb2("cb2","cb2",1);
  probePtSignalData_withOverlap -> SetTitle(""); 
  probePtSignalData_noOverlap   -> SetTitle(""); 
  probePtSignalData_withOverlap -> GetXaxis()->SetTitle("p_{T} [GeV]");
  probePtSignalData_noOverlap   -> GetXaxis()->SetTitle("p_{T} [GeV]");
  probePtSignalData_noOverlap -> DrawNormalized("histe");
  probePtSignalData_withOverlap   -> DrawNormalized("samehiste");
  legB1->Draw();
  cb2.SaveAs("ptData_w-wo-overlap.png");
  //
  TCanvas cb3("cb3","cb3",1);
  probeEtaSignalData_withOverlap -> SetTitle(""); 
  probeEtaSignalData_noOverlap   -> SetTitle(""); 
  probeEtaSignalData_withOverlap -> GetXaxis()->SetTitle("#eta");
  probeEtaSignalData_noOverlap   -> GetXaxis()->SetTitle("#eta");
  probeEtaSignalData_noOverlap -> DrawNormalized("histe");
  probeEtaSignalData_withOverlap   -> DrawNormalized("samehiste");
  legC1->Draw();
  cb3.SaveAs("etaData_w-wo-overlap.png");



  // -------------------------------------------------------------
  // Plots: data vs mc

  // Cosmetics
  probeMvaSignalMc_withOverlap -> SetLineColor(2);  
  probePtSignalMc_withOverlap  -> SetLineColor(2);  
  probeEtaSignalMc_withOverlap -> SetLineColor(2);  
  probeMvaSignalMc_noOverlap   -> SetLineColor(2);  
  probePtSignalMc_noOverlap    -> SetLineColor(2);  
  probeEtaSignalMc_noOverlap   -> SetLineColor(2);  
  //
  probeMvaSignalData_withOverlap -> SetLineColor(4);  
  probePtSignalData_withOverlap  -> SetLineColor(4);  
  probeEtaSignalData_withOverlap -> SetLineColor(4);  
  probeMvaSignalData_noOverlap   -> SetLineColor(4);  
  probePtSignalData_noOverlap    -> SetLineColor(4);  
  probeEtaSignalData_noOverlap   -> SetLineColor(4);  

  TLegend *leg;
  leg = new TLegend(0.65,0.65,0.90,0.90);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.05);
  leg->SetFillColor(0);
  leg->AddEntry(probeMvaSignalMc_withOverlap,   "MC",   "lp");
  leg->AddEntry(probeMvaSignalData_withOverlap, "Data", "lp");
  //
  TLegend *leg2;
  leg2 = new TLegend(0.15,0.65,0.40,0.90);
  leg2->SetFillStyle(0);
  leg2->SetBorderSize(0);
  leg2->SetTextSize(0.05);
  leg2->SetFillColor(0);
  leg2->AddEntry(probeMvaSignalMc_withOverlap,   "MC", "lp");
  leg2->AddEntry(probeMvaSignalData_withOverlap, "Data", "lp");
  //
  TLegend *leg3;
  leg3 = new TLegend(0.70,0.70,0.90,0.90);
  leg3->SetFillStyle(0);
  leg3->SetBorderSize(0);
  leg3->SetTextSize(0.05);
  leg3->SetFillColor(0);
  leg3->AddEntry(probeMvaSignalMc_withOverlap,   "MC", "lp");
  leg3->AddEntry(probeMvaSignalData_withOverlap, "Data", "lp");
  //
  //
  TCanvas cc1("cc1","cc1",1);
  probeEtaSignalData_noOverlap->DrawNormalized("histe");
  probeEtaSignalMc_noOverlap->DrawNormalized("samehiste");
  leg3->Draw();
  cc1.SaveAs("eta_noOverlap_dataVsMc.png");

  TCanvas cc2("cc2","cc2",1);
  probeMvaSignalMc_noOverlap->DrawNormalized("histe");
  probeMvaSignalData_noOverlap->DrawNormalized("samehiste");
  leg2->Draw();
  cc2.SaveAs("mva_noOverlap_dataVsMc.png");

  TCanvas cc3("cc3","cc3",1);
  probePtSignalData_noOverlap->DrawNormalized("histe");
  probePtSignalMc_noOverlap->DrawNormalized("samehiste");
  leg->Draw();
  cc3.SaveAs("pt_noOverlap_dataVsMc.png");

  TCanvas cd1("cd1","cd1",1);
  probeEtaSignalData_withOverlap->DrawNormalized("histe");
  probeEtaSignalMc_withOverlap->DrawNormalized("samehiste");
  leg3->Draw();
  cd1.SaveAs("eta_withOverlap_dataVsMc.png");

  TCanvas cd2("cd2","cd2",1);
  probeMvaSignalMc_withOverlap->DrawNormalized("histe");
  probeMvaSignalData_withOverlap->DrawNormalized("samehiste");
  leg2->Draw();
  cd2.SaveAs("mva_withOverlap_dataVsMc.png");

  TCanvas cd3("cd3","cd3",1);
  probePtSignalMc_withOverlap->DrawNormalized("histe");
  probePtSignalData_withOverlap->DrawNormalized("samehiste");
  leg->Draw();
  cd3.SaveAs("pt_withOverlap_dataVsMc.png");


}
