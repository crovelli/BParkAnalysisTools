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

void mcVsSplots()
{
  // Input files (after TnP selection)
  TFile fileMC("filesNew/myFileMcAfterTnp.root");                               // from prepareInputsFromMcWithTnP
  TFile fileSPlotsInclusive("filesNew/fileSPlots_inclusive.root");              // from jpsi_splot.C 
  TFile fileSPlotsEB0("filesNew/fileSPlots_EB_0d5-to-1d0.root");       
  TFile fileSPlotsEB1("filesNew/fileSPlots_EB_1d0-to-1d5.root");        
  TFile fileSPlotsEB2("filesNew/fileSPlots_EB_1d5-to-2d0.root");        
  TFile fileSPlotsEB3("filesNew/fileSPlots_EB_2d0-to-5d0.root");        
  TFile fileSPlotsEB4("filesNew/fileSPlots_EB_gt5d0.root");             
  TFile fileSPlotsEE0("filesNew/fileSPlots_EE_0d5-to-1d0.root");        
  TFile fileSPlotsEE1("filesNew/fileSPlots_EE_1d0-to-1d5.root");        
  TFile fileSPlotsEE2("filesNew/fileSPlots_EE_1d5-to-2d0.root");        
  TFile fileSPlotsEE3("filesNew/fileSPlots_EE_2d0-to-5d0.root");        
  TFile fileSPlotsEE4("filesNew/fileSPlots_EE_gt5d0.root");            

  // MC histos
  TH1F *probeMvaSignalMc    = (TH1F*)fileMC.Get("probeMvaSignalMc");
  TH1F *probeMvaFakeMc      = (TH1F*)fileMC.Get("probeMvaFakeMc");
  TH1F *probePtSignalMc     = (TH1F*)fileMC.Get("probePtSignalMc");
  TH1F *probePtFakeMc       = (TH1F*)fileMC.Get("probePtFakeMc");
  TH1F *probeEtaSignalMc    = (TH1F*)fileMC.Get("probeEtaSignalMc");
  TH1F *probeEtaFakeMc      = (TH1F*)fileMC.Get("probeEtaFakeMc");
  TH1F *probeFBremSignalMc  = (TH1F*)fileMC.Get("probeFBremSignalMc");
  TH1F *probeFBremFakeMc    = (TH1F*)fileMC.Get("probeFBremFakeMc");
  TH1F *probeDxySigSignalMc = (TH1F*)fileMC.Get("probeDxySigSignalMc");
  TH1F *probeDxySigFakeMc   = (TH1F*)fileMC.Get("probeDxySigFakeMc");
  TH1F *probeDzSigSignalMc  = (TH1F*)fileMC.Get("probeDzSigSignalMc");
  TH1F *probeDzSigFakeMc    = (TH1F*)fileMC.Get("probeDzSigFakeMc");
  probeMvaSignalMc->Sumw2();   
  probeMvaFakeMc->Sumw2();   
  probePtSignalMc->Sumw2(); 
  probePtFakeMc->Sumw2(); 
  probeEtaSignalMc->Sumw2(); 
  probeEtaFakeMc->Sumw2(); 
  probeFBremSignalMc->Sumw2(); 
  probeFBremFakeMc->Sumw2(); 
  probeDxySigSignalMc->Sumw2(); 
  probeDxySigFakeMc->Sumw2(); 
  probeDzSigSignalMc->Sumw2(); 
  probeDzSigFakeMc->Sumw2(); 

  // MC histos in eta/pT bins
  TH1F *probeMvaSignalMcEB0 = (TH1F*)fileMC.Get("mvaSignalEBMc0");
  TH1F *probeMvaSignalMcEB1 = (TH1F*)fileMC.Get("mvaSignalEBMc1");
  TH1F *probeMvaSignalMcEB2 = (TH1F*)fileMC.Get("mvaSignalEBMc2");
  TH1F *probeMvaSignalMcEB3 = (TH1F*)fileMC.Get("mvaSignalEBMc3");
  TH1F *probeMvaSignalMcEB4 = (TH1F*)fileMC.Get("mvaSignalEBMc4");
  TH1F *probeMvaSignalMcEE0 = (TH1F*)fileMC.Get("mvaSignalEEMc0");
  TH1F *probeMvaSignalMcEE1 = (TH1F*)fileMC.Get("mvaSignalEEMc1");
  TH1F *probeMvaSignalMcEE2 = (TH1F*)fileMC.Get("mvaSignalEEMc2");
  TH1F *probeMvaSignalMcEE3 = (TH1F*)fileMC.Get("mvaSignalEEMc3");
  TH1F *probeMvaSignalMcEE4 = (TH1F*)fileMC.Get("mvaSignalEEMc4");
  probeMvaSignalMcEB0->Sumw2();  
  probeMvaSignalMcEB1->Sumw2();  
  probeMvaSignalMcEB2->Sumw2();  
  probeMvaSignalMcEB3->Sumw2();  
  probeMvaSignalMcEB4->Sumw2();  
  probeMvaSignalMcEE0->Sumw2();  
  probeMvaSignalMcEE1->Sumw2();  
  probeMvaSignalMcEE2->Sumw2();  
  probeMvaSignalMcEE3->Sumw2();  
  probeMvaSignalMcEE4->Sumw2();  

  // s-Plots output
  TH1F *probeMvaSignalData    = (TH1F*)fileSPlotsInclusive.Get("h1_probeMvaId_jpsi__probeMvaId");
  TH1F *probeMvaBkgData       = (TH1F*)fileSPlotsInclusive.Get("h1_probeMvaId_bkg__probeMvaId");
  TH1F *probePtSignalData     = (TH1F*)fileSPlotsInclusive.Get("h1_probePt_jpsi__probePt");
  TH1F *probePtBkgData        = (TH1F*)fileSPlotsInclusive.Get("h1_probePt_bkg__probePt");
  TH1F *probeEtaSignalData    = (TH1F*)fileSPlotsInclusive.Get("h1_probeEta_jpsi__probeEta");
  TH1F *probeEtaBkgData       = (TH1F*)fileSPlotsInclusive.Get("h1_probeEta_bkg__probeEta");
  TH1F *probeFBremSignalData  = (TH1F*)fileSPlotsInclusive.Get("h1_probeFBrem_jpsi__probeFBrem");
  TH1F *probeFBremBkgData     = (TH1F*)fileSPlotsInclusive.Get("h1_probeFBrem_bkg__probeFBrem");
  TH1F *probeDxySigSignalData = (TH1F*)fileSPlotsInclusive.Get("h1_probeDxySig_jpsi__probeDxySig");
  TH1F *probeDxySigBkgData    = (TH1F*)fileSPlotsInclusive.Get("h1_probeDxySig_bkg__probeDxySig");
  TH1F *probeDzSigSignalData  = (TH1F*)fileSPlotsInclusive.Get("h1_probeDzSig_jpsi__probeDzSig");
  TH1F *probeDzSigBkgData     = (TH1F*)fileSPlotsInclusive.Get("h1_probeDzSig_bkg__probeDzSig");
  probeMvaSignalData->Sumw2();   
  probeMvaBkgData->Sumw2();   
  probePtSignalData->Sumw2(); 
  probePtBkgData->Sumw2(); 
  probeEtaSignalData->Sumw2(); 
  probeEtaBkgData->Sumw2(); 
  probeFBremSignalData->Sumw2(); 
  probeFBremBkgData->Sumw2(); 
  probeDxySigSignalData->Sumw2(); 
  probeDxySigBkgData->Sumw2(); 
  probeDzSigSignalData->Sumw2(); 
  probeDzSigBkgData->Sumw2(); 

  // s-Plots output in eta/pT bins
  TH1F *probeMvaSignalDataEB0 = (TH1F*)fileSPlotsEB0.Get("h1_probeMvaId_jpsi__probeMvaId");
  TH1F *probeMvaSignalDataEB1 = (TH1F*)fileSPlotsEB1.Get("h1_probeMvaId_jpsi__probeMvaId");
  TH1F *probeMvaSignalDataEB2 = (TH1F*)fileSPlotsEB2.Get("h1_probeMvaId_jpsi__probeMvaId");
  TH1F *probeMvaSignalDataEB3 = (TH1F*)fileSPlotsEB3.Get("h1_probeMvaId_jpsi__probeMvaId");
  TH1F *probeMvaSignalDataEB4 = (TH1F*)fileSPlotsEB4.Get("h1_probeMvaId_jpsi__probeMvaId");
  TH1F *probeMvaSignalDataEE0 = (TH1F*)fileSPlotsEE0.Get("h1_probeMvaId_jpsi__probeMvaId");
  TH1F *probeMvaSignalDataEE1 = (TH1F*)fileSPlotsEE1.Get("h1_probeMvaId_jpsi__probeMvaId");
  TH1F *probeMvaSignalDataEE2 = (TH1F*)fileSPlotsEE2.Get("h1_probeMvaId_jpsi__probeMvaId");
  TH1F *probeMvaSignalDataEE3 = (TH1F*)fileSPlotsEE3.Get("h1_probeMvaId_jpsi__probeMvaId");
  TH1F *probeMvaSignalDataEE4 = (TH1F*)fileSPlotsEE4.Get("h1_probeMvaId_jpsi__probeMvaId");
  probeMvaSignalDataEB0->Sumw2();   
  probeMvaSignalDataEB1->Sumw2();   
  probeMvaSignalDataEB2->Sumw2();   
  probeMvaSignalDataEB3->Sumw2();   
  probeMvaSignalDataEB4->Sumw2();   
  probeMvaSignalDataEE0->Sumw2();   
  probeMvaSignalDataEE1->Sumw2();   
  probeMvaSignalDataEE2->Sumw2();   
  probeMvaSignalDataEE3->Sumw2();   
  probeMvaSignalDataEE4->Sumw2();   

  // Cosmetics
  probeMvaSignalMc    -> SetLineWidth(2);  
  probeMvaFakeMc      -> SetLineWidth(2);  
  probeMvaSignalMc    -> SetLineColor(6);  
  probeMvaFakeMc      -> SetLineColor(4);  
  probePtSignalMc     -> SetLineWidth(2);  
  probePtFakeMc       -> SetLineWidth(2);  
  probePtSignalMc     -> SetLineColor(6);  
  probePtFakeMc       -> SetLineColor(4);  
  probeEtaSignalMc    -> SetLineWidth(2);  
  probeEtaFakeMc      -> SetLineWidth(2);  
  probeEtaSignalMc    -> SetLineColor(6);  
  probeEtaFakeMc      -> SetLineColor(4);  
  probeFBremSignalMc  -> SetLineWidth(2);  
  probeFBremFakeMc    -> SetLineWidth(2);  
  probeFBremSignalMc  -> SetLineColor(6);  
  probeFBremFakeMc    -> SetLineColor(4);  
  probeDxySigSignalMc -> SetLineWidth(2);
  probeDxySigFakeMc   -> SetLineWidth(2);
  probeDxySigSignalMc -> SetLineColor(6);
  probeDxySigFakeMc   -> SetLineColor(4);
  probeDzSigSignalMc  -> SetLineWidth(2);
  probeDzSigFakeMc    -> SetLineWidth(2);
  probeDzSigSignalMc  -> SetLineColor(6);
  probeDzSigFakeMc    -> SetLineColor(4);
  //
  probeMvaSignalMcEB0->SetLineWidth(2);  
  probeMvaSignalMcEB0->SetLineColor(6);  
  probeMvaSignalMcEB1->SetLineWidth(2);  
  probeMvaSignalMcEB1->SetLineColor(6);  
  probeMvaSignalMcEB2->SetLineWidth(2);  
  probeMvaSignalMcEB2->SetLineColor(6);  
  probeMvaSignalMcEB3->SetLineWidth(2);  
  probeMvaSignalMcEB3->SetLineColor(6);  
  probeMvaSignalMcEB4->SetLineWidth(2);  
  probeMvaSignalMcEB4->SetLineColor(6);  
  probeMvaSignalMcEE0->SetLineWidth(2);  
  probeMvaSignalMcEE0->SetLineColor(6);  
  probeMvaSignalMcEE1->SetLineWidth(2);  
  probeMvaSignalMcEE1->SetLineColor(6);  
  probeMvaSignalMcEE2->SetLineWidth(2);  
  probeMvaSignalMcEE2->SetLineColor(6);  
  probeMvaSignalMcEE3->SetLineWidth(2);  
  probeMvaSignalMcEE3->SetLineColor(6);  
  probeMvaSignalMcEE4->SetLineWidth(2);  
  probeMvaSignalMcEE4->SetLineColor(6);  
  //
  probeMvaSignalData    -> SetLineWidth(2);  
  probeMvaBkgData       -> SetLineWidth(2);  
  probeMvaSignalData    -> SetLineColor(3);  
  probeMvaBkgData       -> SetLineColor(2);  
  probePtSignalData     -> SetLineWidth(2);  
  probePtBkgData        -> SetLineWidth(2);  
  probePtSignalData     -> SetLineColor(3);  
  probePtBkgData        -> SetLineColor(2);  
  probeEtaSignalData    -> SetLineWidth(2);  
  probeEtaBkgData       -> SetLineWidth(2);  
  probeEtaSignalData    -> SetLineColor(3);  
  probeEtaBkgData       -> SetLineColor(2);  
  probeFBremSignalData  -> SetLineWidth(2);  
  probeFBremBkgData     -> SetLineWidth(2);  
  probeFBremSignalData  -> SetLineColor(3);  
  probeFBremBkgData     -> SetLineColor(2);  
  probeDxySigSignalData -> SetLineWidth(2);
  probeDxySigBkgData    -> SetLineWidth(2);
  probeDxySigSignalData -> SetLineColor(3);
  probeDxySigBkgData    -> SetLineColor(2);
  probeDzSigSignalData  -> SetLineWidth(2);
  probeDzSigBkgData     -> SetLineWidth(2);
  probeDzSigSignalData  -> SetLineColor(3);
  probeDzSigBkgData     -> SetLineColor(2);
  //
  probeMvaSignalDataEB0 -> SetLineWidth(2);
  probeMvaSignalDataEB0 -> SetLineColor(3);
  probeMvaSignalDataEB1 -> SetLineWidth(2);
  probeMvaSignalDataEB1 -> SetLineColor(3);
  probeMvaSignalDataEB2 -> SetLineWidth(2);
  probeMvaSignalDataEB2 -> SetLineColor(3);
  probeMvaSignalDataEB3 -> SetLineWidth(2);
  probeMvaSignalDataEB3 -> SetLineColor(3);
  probeMvaSignalDataEB4 -> SetLineWidth(2);
  probeMvaSignalDataEB4 -> SetLineColor(3);
  probeMvaSignalDataEE0 -> SetLineWidth(2);
  probeMvaSignalDataEE0 -> SetLineColor(3);
  probeMvaSignalDataEE1 -> SetLineWidth(2);
  probeMvaSignalDataEE1 -> SetLineColor(3);
  probeMvaSignalDataEE2 -> SetLineWidth(2);
  probeMvaSignalDataEE2 -> SetLineColor(3);
  probeMvaSignalDataEE3 -> SetLineWidth(2);
  probeMvaSignalDataEE3 -> SetLineColor(3);
  probeMvaSignalDataEE4 -> SetLineWidth(2);
  probeMvaSignalDataEE4 -> SetLineColor(3);

  
  // Plots
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);

  // Outputs
  // --------------------------------------------------
  TLegend *leg;
  leg = new TLegend(0.15,0.55,0.50,0.80);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.05);
  leg->SetFillColor(0);
  leg->AddEntry(probeMvaSignalMc,   "Matched, from MC", "lp");
  leg->AddEntry(probeMvaFakeMc,     "Fake, from MC", "lp");
  leg->AddEntry(probeMvaSignalData, "Signal, from sPlot", "lp");
  leg->AddEntry(probeMvaBkgData,    "Bkg, from sPlot", "lp");
  //
  TLegend *legS;
  legS = new TLegend(0.15,0.55,0.50,0.80);
  legS->SetFillStyle(0);
  legS->SetBorderSize(0);
  legS->SetTextSize(0.05);
  legS->SetFillColor(0);
  legS->AddEntry(probeMvaSignalMc,   "Matched, from MC", "lp");
  legS->AddEntry(probeMvaSignalData, "Signal, from sPlot", "lp");
  //
  TLegend *legB;
  legB = new TLegend(0.50,0.55,0.75,0.80);
  legB->SetFillStyle(0);
  legB->SetBorderSize(0);
  legB->SetTextSize(0.05);
  legB->SetFillColor(0);
  legB->AddEntry(probeMvaFakeMc,  "Fake, from MC", "lp");
  legB->AddEntry(probeMvaBkgData, "Bkg, from sPlot", "lp");

  TCanvas cmva("cmva","cmva",1);
  probeMvaSignalMc->SetTitle("");
  probeMvaFakeMc->SetTitle("");
  probeMvaSignalData->SetTitle("");
  probeMvaBkgData->SetTitle("");
  probeMvaSignalMc->GetXaxis()->SetTitle("Id BDT");
  probeMvaFakeMc->GetXaxis()->SetTitle("Id BDT");
  probeMvaSignalData->GetXaxis()->SetTitle("Id BDT");
  probeMvaBkgData->GetXaxis()->SetTitle("Id BDT");
  probeMvaSignalMc->DrawNormalized("hist");
  probeMvaFakeMc->DrawNormalized("samehist");
  probeMvaSignalData->DrawNormalized("samehist");
  probeMvaBkgData->DrawNormalized("samehist");
  leg->Draw();
  cmva.SaveAs("outputBDT.png");
  
  TCanvas cmvas("cmvas","cmvas",1);
  probeMvaSignalMc->DrawNormalized("hist");
  probeMvaSignalData->DrawNormalized("samehist");
  legS->Draw();
  cmvas.SaveAs("outputBDT_signal.png");

  TCanvas cmvab("cmvab","cmvab",1);
  probeMvaFakeMc->DrawNormalized("hist");
  probeMvaBkgData->DrawNormalized("samehist");
  legB->Draw();
  cmvab.SaveAs("outputBDT_back.png");


  // -----------------------------------------------
  TLegend *leg1;
  leg1 = new TLegend(0.50,0.55,0.75,0.80);
  leg1->SetFillStyle(0);
  leg1->SetBorderSize(0);
  leg1->SetTextSize(0.05);
  leg1->SetFillColor(0);
  leg1->AddEntry(probeMvaSignalMc,   "Matched, from MC", "lp");
  leg1->AddEntry(probeMvaFakeMc,     "Fake, from MC", "lp");
  leg1->AddEntry(probeMvaSignalData, "Signal, from sPlot", "lp");
  leg1->AddEntry(probeMvaBkgData,    "Bkg, from sPlot", "lp");
  //
  TLegend *leg1S;
  leg1S = new TLegend(0.50,0.55,0.75,0.80);
  leg1S->SetFillStyle(0);
  leg1S->SetBorderSize(0);
  leg1S->SetTextSize(0.05);
  leg1S->SetFillColor(0);
  leg1S->AddEntry(probeMvaSignalMc,   "Matched, from MC", "lp");
  leg1S->AddEntry(probeMvaSignalData, "Signal, from sPlot", "lp");
  //
  TLegend *leg1B;
  leg1B = new TLegend(0.50,0.55,0.75,0.80);
  leg1B->SetFillStyle(0);
  leg1B->SetBorderSize(0);
  leg1B->SetTextSize(0.05);
  leg1B->SetFillColor(0);
  leg1B->AddEntry(probeMvaFakeMc,  "Fake, from MC", "lp");
  leg1B->AddEntry(probeMvaBkgData, "Bkg, from sPlot", "lp");

  TCanvas cprobePts("cprobePts","cprobePts",1);
  probePtSignalData->DrawNormalized("hist");
  probePtSignalMc->DrawNormalized("samehist");
  leg1S->Draw();
  cprobePts.SaveAs("probePt_signal.png");

  TCanvas cprobePtb("cprobePtb","cprobePtb",1);
  probePtBkgData->DrawNormalized("hist");
  probePtFakeMc->DrawNormalized("samehist");
  leg1B->Draw();
  cprobePtb.SaveAs("probePt_back.png");


  // -----------------------------------------------
  TLegend *leg2;
  leg2 = new TLegend(0.30,0.10,0.65,0.40);
  leg2->SetFillStyle(0);
  leg2->SetBorderSize(0);
  leg2->SetTextSize(0.05);
  leg2->SetFillColor(0);
  leg2->AddEntry(probeMvaSignalMc,   "Matched, from MC", "lp");
  leg2->AddEntry(probeMvaFakeMc,     "Fake, from MC", "lp");
  leg2->AddEntry(probeMvaSignalData, "Signal, from sPlot", "lp");
  leg2->AddEntry(probeMvaBkgData,    "Bkg, from sPlot", "lp");
  //
  TLegend *leg2S;
  leg2S = new TLegend(0.30,0.10,0.65,0.40);
  leg2S->SetFillStyle(0);
  leg2S->SetBorderSize(0);
  leg2S->SetTextSize(0.05);
  leg2S->SetFillColor(0);
  leg2S->AddEntry(probeMvaSignalMc,   "Matched, from MC", "lp");
  leg2S->AddEntry(probeMvaSignalData, "Signal, from sPlot", "lp");
  //
  TLegend *leg2B;
  leg2B = new TLegend(0.30,0.10,0.65,0.40);
  leg2B->SetFillStyle(0);
  leg2B->SetBorderSize(0);
  leg2B->SetTextSize(0.05);
  leg2B->SetFillColor(0);
  leg2B->AddEntry(probeMvaFakeMc,  "Fake, from MC", "lp");
  leg2B->AddEntry(probeMvaBkgData, "Bkg, from sPlot", "lp");

  TCanvas cprobeEtas("cprobeEtas","cprobeEtas",1);
  probeEtaSignalData->DrawNormalized("hist");
  probeEtaSignalMc->DrawNormalized("samehist");
  leg2S->Draw();
  cprobeEtas.SaveAs("probeEta_signal.png");

  TCanvas cprobeEtab("cprobeEtab","cprobeEtab",1);
  probeEtaBkgData->DrawNormalized("hist");
  probeEtaFakeMc->DrawNormalized("samehist");
  leg2B->Draw();
  cprobeEtab.SaveAs("probeEta_back.png");


  // -----------------------------------------------
  TLegend *leg3;
  leg3 = new TLegend(0.40,0.50,0.75,0.80);
  leg3->SetFillStyle(0);
  leg3->SetBorderSize(0);
  leg3->SetTextSize(0.05);
  leg3->SetFillColor(0);
  leg3->AddEntry(probeMvaSignalMc,   "Matched, from MC", "lp");
  leg3->AddEntry(probeMvaFakeMc,     "Fake, from MC", "lp");
  leg3->AddEntry(probeMvaSignalData, "Signal, from sPlot", "lp");
  leg3->AddEntry(probeMvaBkgData,    "Bkg, from sPlot", "lp");
  //
  TLegend *leg3S;
  leg3S = new TLegend(0.40,0.50,0.75,0.80);
  leg3S->SetFillStyle(0);
  leg3S->SetBorderSize(0);
  leg3S->SetTextSize(0.05);
  leg3S->SetFillColor(0);
  leg3S->AddEntry(probeMvaSignalMc,   "Matched, from MC", "lp");
  leg3S->AddEntry(probeMvaSignalData, "Signal, from sPlot", "lp");
  //
  TLegend *leg3B;
  leg3B = new TLegend(0.40,0.50,0.75,0.80);
  leg3B->SetFillStyle(0);
  leg3B->SetBorderSize(0);
  leg3B->SetTextSize(0.05);
  leg3B->SetFillColor(0);
  leg3B->AddEntry(probeMvaFakeMc,  "Fake, from MC", "lp");
  leg3B->AddEntry(probeMvaBkgData, "Bkg, from sPlot", "lp");

  TCanvas cprobeFBrems("cprobeFBrems","cprobeFBrems",1);
  probeFBremSignalData->DrawNormalized("hist");
  probeFBremSignalMc->DrawNormalized("samehist");
  leg3S->Draw();
  cprobeFBrems.SaveAs("probeFBrem_signal.png");

  TCanvas cprobeFBremb("cprobeFBremb","cprobeFBremb",1);
  probeFBremBkgData->DrawNormalized("hist");
  probeFBremFakeMc->DrawNormalized("samehist");
  leg3B->Draw();
  cprobeFBremb.SaveAs("probeFBrem_back.png");


  // -----------------------------------------------
  TLegend *leg4;
  leg4 = new TLegend(0.50,0.50,0.80,0.80);
  leg4->SetFillStyle(0);
  leg4->SetBorderSize(0);
  leg4->SetTextSize(0.05);
  leg4->SetFillColor(0);
  leg4->AddEntry(probeMvaSignalMc,   "Matched, from MC", "lp");
  leg4->AddEntry(probeMvaFakeMc,     "Fake, from MC", "lp");
  leg4->AddEntry(probeMvaSignalData, "Signal, from sPlot", "lp");
  leg4->AddEntry(probeMvaBkgData,    "Bkg, from sPlot", "lp");
  //
  TLegend *leg4S;
  leg4S = new TLegend(0.50,0.50,0.80,0.80);
  leg4S->SetFillStyle(0);
  leg4S->SetBorderSize(0);
  leg4S->SetTextSize(0.05);
  leg4S->SetFillColor(0);
  leg4S->AddEntry(probeMvaSignalMc,   "Matched, from MC", "lp");
  leg4S->AddEntry(probeMvaSignalData, "Signal, from sPlot", "lp");
  //
  TLegend *leg4B;
  leg4B = new TLegend(0.50,0.50,0.80,0.80);
  leg4B->SetFillStyle(0);
  leg4B->SetBorderSize(0);
  leg4B->SetTextSize(0.05);
  leg4B->SetFillColor(0);
  leg4B->AddEntry(probeMvaFakeMc,  "Fake, from MC", "lp");
  leg4B->AddEntry(probeMvaBkgData, "Bkg, from sPlot", "lp");

  TCanvas cprobeDxySigs("cprobeDxySigs","cprobeDxySigs",1);
  probeDxySigSignalData->DrawNormalized("hist");
  probeDxySigSignalMc->DrawNormalized("samehist");
  leg4S->Draw();
  cprobeDxySigs.SaveAs("probeDxySig_signal.png");

  TCanvas cprobeDxySigb("cprobeDxySigb","cprobeDxySigb",1);
  probeDxySigBkgData->DrawNormalized("hist");
  probeDxySigFakeMc->DrawNormalized("samehist");
  leg4B->Draw();
  cprobeDxySigb.SaveAs("probeDxySig_back.png");


  // -----------------------------------------------
  TCanvas cprobeDzSigs("cprobeDzSigs","cprobeDzSigs",1);
  probeDzSigSignalMc->DrawNormalized("hist");
  probeDzSigSignalData->DrawNormalized("samehist");
  leg4S->Draw();
  cprobeDzSigs.SaveAs("probeDzSig_signal.png");

  TCanvas cprobeDzSigb("cprobeDzSigb","cprobeDzSigb",1);
  probeDzSigBkgData->DrawNormalized("hist");
  probeDzSigFakeMc->DrawNormalized("samehist");
  leg4B->Draw();
  cprobeDzSigb.SaveAs("probeDzSig_back.png");


  // -----------------------------------------------
  TCanvas cmvaseb0("cmvaseb0","cmvaseb0",1);
  probeMvaSignalDataEB0->Rebin();
  probeMvaSignalMcEB0->Rebin();
  probeMvaSignalDataEB0->SetTitle("EB, 0.5 < pT < 1");
  probeMvaSignalMcEB0->SetTitle("EB, 0.5 < pT < 1");
  probeMvaSignalDataEB0->DrawNormalized("hist");
  probeMvaSignalMcEB0->DrawNormalized("samehist");
  legS->Draw();
  cmvaseb0.SaveAs("outputBDT_EB0_signal.png");
  //
  TCanvas cmvaseb1("cmvaseb1","cmvaseb1",1);
  probeMvaSignalMcEB1->Rebin();
  probeMvaSignalDataEB1->Rebin();
  probeMvaSignalDataEB1->SetTitle("EB, 1.0 < pT < 1.5");
  probeMvaSignalMcEB1->SetTitle("EB, 1.0 < pT < 1.5");
  probeMvaSignalDataEB1->DrawNormalized("hist");
  probeMvaSignalMcEB1->DrawNormalized("samehist");
  legS->Draw();
  cmvaseb1.SaveAs("outputBDT_EB1_signal.png");
  //
  TCanvas cmvaseb2("cmvaseb2","cmvaseb2",1);
  probeMvaSignalMcEB2->Rebin();
  probeMvaSignalDataEB2->Rebin();
  probeMvaSignalDataEB2->SetTitle("EB, 1.5 < pT < 2.0");
  probeMvaSignalMcEB2->SetTitle("EB, 1.5 < pT < 2.0");
  probeMvaSignalMcEB2->DrawNormalized("hist");
  probeMvaSignalDataEB2->DrawNormalized("samehist");
  legS->Draw();
  cmvaseb2.SaveAs("outputBDT_EB2_signal.png");
  //
  TCanvas cmvaseb3("cmvaseb3","cmvaseb3",1);
  probeMvaSignalMcEB3->Rebin();
  probeMvaSignalDataEB3->Rebin();
  probeMvaSignalDataEB3->SetTitle("EB, 2.0 < pT < 5.0");
  probeMvaSignalMcEB3->SetTitle("EB, 2.0 < pT < 5.0");
  probeMvaSignalMcEB3->DrawNormalized("hist");
  probeMvaSignalDataEB3->DrawNormalized("samehist");
  legS->Draw();
  cmvaseb3.SaveAs("outputBDT_EB3_signal.png");
  //
  TCanvas cmvaseb4("cmvaseb4","cmvaseb4",1);
  probeMvaSignalMcEB4->Rebin();
  probeMvaSignalDataEB4->Rebin();
  probeMvaSignalDataEB4->SetTitle("EB, pT > 5.0");
  probeMvaSignalMcEB4->SetTitle("EB, pT > 5.0");
  probeMvaSignalMcEB4->DrawNormalized("hist");
  probeMvaSignalDataEB4->DrawNormalized("samehist");
  legS->Draw();
  cmvaseb4.SaveAs("outputBDT_EB4_signal.png");
  //
  TCanvas cmvasee0("cmvasee0","cmvasee0",1);
  probeMvaSignalMcEE0->Rebin();
  probeMvaSignalDataEE0->Rebin();
  probeMvaSignalDataEE0->SetTitle("EE, 0.5 < pT < 1.0");
  probeMvaSignalMcEE0->SetTitle("EE, 0.5 < pT < 1.0");
  probeMvaSignalDataEE0->DrawNormalized("hist");
  probeMvaSignalMcEE0->DrawNormalized("samehist");
  legS->Draw();
  cmvasee0.SaveAs("outputBDT_EE0_signal.png");
  //
  TCanvas cmvasee1("cmvasee1","cmvasee1",1);
  probeMvaSignalMcEE1->Rebin();
  probeMvaSignalDataEE1->Rebin();
  probeMvaSignalDataEE1->SetTitle("EE, 1.0 < pT < 1.5");
  probeMvaSignalMcEE1->SetTitle("EE, 1.0 < pT < 1.5");
  probeMvaSignalDataEE1->DrawNormalized("hist");
  probeMvaSignalMcEE1->DrawNormalized("samehist");
  legS->Draw();
  cmvasee1.SaveAs("outputBDT_EE1_signal.png");
  //
  TCanvas cmvasee2("cmvasee2","cmvasee2",1);
  probeMvaSignalMcEE2->Rebin();
  probeMvaSignalDataEE2->Rebin();
  probeMvaSignalDataEE2->SetTitle("EE, 1.5 < pT < 2.0");
  probeMvaSignalMcEE2->SetTitle("EE, 1.5 < pT < 2.0");
  probeMvaSignalDataEE2->DrawNormalized("hist");
  probeMvaSignalMcEE2->DrawNormalized("samehist");
  legS->Draw();
  cmvasee2.SaveAs("outputBDT_EE2_signal.png");
  //
  TCanvas cmvasee3("cmvasee3","cmvasee3",1);
  probeMvaSignalMcEE3->Rebin();
  probeMvaSignalDataEE3->Rebin();
  probeMvaSignalDataEE3->SetTitle("EE, 2.0 < pT < 5.0");
  probeMvaSignalMcEE3->SetTitle("EE, 2.0 < pT < 5.0");
  probeMvaSignalMcEE3->DrawNormalized("hist");
  probeMvaSignalDataEE3->DrawNormalized("samehist");
  legS->Draw();
  cmvasee3.SaveAs("outputBDT_EE3_signal.png");
  //
  TCanvas cmvasee4("cmvasee4","cmvasee4",1);
  probeMvaSignalMcEE4->Rebin();
  probeMvaSignalDataEE4->Rebin();
  probeMvaSignalDataEE3->SetTitle("EE, pT > 5.0");
  probeMvaSignalMcEE3->SetTitle("EE, pT > 5.0");
  probeMvaSignalMcEE4->DrawNormalized("hist");
  probeMvaSignalDataEE4->DrawNormalized("samehist");
  legS->Draw();
  cmvasee4.SaveAs("outputBDT_EE4_signal.png");
}
