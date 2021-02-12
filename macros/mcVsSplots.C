#define mcVsSplots_cxx

// Data vs MC signal distributions

// ROOT includes 
#include <TROOT.h>
#include <TStyle.h>
#include <TF1.h>
#include <TH2.h>
#include <TFile.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <iostream>

using namespace std;

float drawTH1pair(TH1* h1, TH1* h2, 
		  const string& xAxisNameTmp = "", const string& yAxisName = "Events", 
		  float lumi=-1, const string& canvasName = "default", 
		  const string& outputDIR = "./", int mycolor=2,
		  const string& legEntry1 = "data", const string& legEntry2 = "MC", const string& ratioPadYaxisName = "data/MC") 
{
  string xAxisName = "";
  string separator = "::";
  Bool_t setXAxisRangeFromUser = false;
  Double_t xmin = 0;
  Double_t xmax = 0;

  size_t pos = xAxisNameTmp.find(separator);
  if (pos != string::npos) {
    string xrange = "";
    setXAxisRangeFromUser = true;
    xAxisName.assign(xAxisNameTmp, 0, pos); 
    xrange.assign(xAxisNameTmp, pos + separator.size(), string::npos);
    separator = ",";
    pos = xrange.find(separator);
    string numString = ""; 
    numString.assign(xrange,0,pos);
    xmin = std::stod(numString);
    numString.assign(xrange,pos + separator.size(), string::npos);
    xmax = std::stod(numString);
  } else {
    xAxisName = xAxisNameTmp;
  }

  if (yAxisName == "a.u.") {
    h1->Scale(1./h1->Integral());
    h2->Scale(1./h2->Integral());
  }
  else if (lumi>-1) {
    h1->Scale(lumi/h1->Integral());
    h2->Scale(lumi/h2->Integral());
  }

  // To deal with splots
  for (int ii=0; ii<h1->GetNbinsX(); ii++){
    int iip1 = ii+1;
    if (h1->GetBinContent(iip1)<=0) { 
      h1->SetBinContent(iip1,0);
      h1->SetBinError(iip1,0);
    } 
    if (h2->GetBinContent(iip1)<=0) {
      h2->SetBinContent(iip1,0);
      h2->SetBinError(iip1,0);
    } 
  }
  
  h1->SetStats(0);
  h2->SetStats(0);

  TCanvas* canvas = new TCanvas("canvas","",600,700);
  canvas->cd();
  canvas->SetTickx(1);
  canvas->SetTicky(1);
  canvas->cd();
  canvas->SetBottomMargin(0.3);
  canvas->SetRightMargin(0.06);

  TPad *pad2 = new TPad("pad2","pad2",0,0.,1,0.9);
  pad2->SetTopMargin(0.7);
  pad2->SetRightMargin(0.06);
  pad2->SetFillColor(0);
  pad2->SetGridy(1);
  pad2->SetFillStyle(0);

  TH1* frame =  (TH1*) h1->Clone("frame");
  frame->SetTitle("");
  frame->GetXaxis()->SetLabelSize(0.04);
  frame->SetStats(0);

  h1->SetLineColor(kBlack);
  h1->SetMarkerColor(kBlack);
  h1->SetMarkerStyle(20);
  h1->SetMarkerSize(1);

  h1->SetTitle("");
  h1->GetXaxis()->SetLabelSize(0);
  h1->GetYaxis()->SetTitle(yAxisName.c_str());
  h1->GetYaxis()->SetTitleOffset(1.1);
  h1->GetYaxis()->SetTitleSize(0.05);
  h1->GetYaxis()->SetRangeUser(0.0, max(h1->GetMaximum(),h2->GetMaximum()) * 1.2);
  if (setXAxisRangeFromUser) h1->GetXaxis()->SetRangeUser(xmin,xmax);
  h1->Draw("EP");

  h2->SetTitle("");
  h2->SetLineColor(mycolor);
  h2->SetLineWidth(2);
  h2->Draw("hist E same");

  TLegend leg (0.15,0.65,0.45,0.85);
  leg.SetFillColor(0);
  leg.SetFillStyle(0);
  leg.SetBorderSize(0);
  leg.AddEntry(h1,legEntry1.c_str(),"PLE");
  leg.AddEntry(h2,legEntry2.c_str(),"L");
  leg.Draw("same");
  canvas->RedrawAxis("sameaxis");

  pad2->Draw();
  pad2->cd();

  frame->Reset("ICES");
  frame->GetYaxis()->SetRangeUser(0.,4.);
  frame->GetYaxis()->SetNdivisions(5);
  frame->GetYaxis()->SetTitle(ratioPadYaxisName.c_str());
  frame->GetYaxis()->SetTitleOffset(1.2);
  frame->GetYaxis()->CenterTitle();
  frame->GetXaxis()->SetTitle(xAxisName.c_str());
  if (setXAxisRangeFromUser) frame->GetXaxis()->SetRangeUser(xmin,xmax);
  frame->GetXaxis()->SetTitleSize(0.05);
  
  TH1D* ratio = (TH1D*) h1->Clone("ratio");
  ratio->Sumw2();
  TH1D* den_noerr = (TH1D*) h2->Clone("den_noerr");
  den_noerr->Sumw2();
  TH1D* den = (TH1D*) h2->Clone("den");
  den->Sumw2();
  for(int iBin = 1; iBin < den->GetNbinsX()+1; iBin++) den_noerr->SetBinError(iBin,0.);

  ratio->Divide(den_noerr);
  den->Divide(den_noerr);
  den->SetFillColor(kGray);
  frame->Draw();
  ratio->SetMarkerSize(0.85);
  ratio->Draw("EPsame");
  den->Draw("E2same");

  TF1* line = new TF1("horiz_line","1",ratio->GetXaxis()->GetBinLowEdge(1),ratio->GetXaxis()->GetBinLowEdge(ratio->GetNbinsX()+1));
  line->SetLineColor(mycolor);
  line->SetLineWidth(2);
  line->Draw("Lsame");
  ratio->Draw("EPsame");
  pad2->RedrawAxis("sameaxis");

  canvas->SaveAs((outputDIR + canvasName + ".png").c_str());

  if (yAxisName == "a.u.") h1->GetYaxis()->SetRangeUser(max(0.0001,min(h1->GetMinimum(),h2->GetMinimum())*0.8),max(h1->GetMaximum(),h2->GetMaximum())*100);
  else h1->GetYaxis()->SetRangeUser(max(0.001,min(h1->GetMinimum(),h2->GetMinimum())*0.8),max(h1->GetMaximum(),h2->GetMaximum())*100);
  canvas->SetLogy(0);

  delete canvas;
  frame->Reset("ICES");
}

void mcVsSplots(int doLowPt)
{
  // Input files (after TnP selection)
  TFile *fileMC;
  TFile *fileSPlotsInclusive;
  if (doLowPt==1) {
    fileMC = new TFile("files_v2/probeLowPt/myFileMcAfterTnp__withPUweight.root");                  // from prepareInputsFromMcWithTnP
    fileSPlotsInclusive = new TFile("files_v2/probeLowPt/myFileSPlots_fullRange.root");             // from jpsi_splot.C 
  }
  if (doLowPt==0) {
    fileMC = new TFile("files_v2/probePF/myFileMcAfterTnp.root");                               // from prepareInputsFromMcWithTnP
    fileSPlotsInclusive = new TFile("files_v2/probePF/fileSPlots_inclusive.root");              // from jpsi_splot.C 
  }

  TFile *fileSPlotsEB0;
  TFile *fileSPlotsEB1;
  TFile *fileSPlotsEB2;
  TFile *fileSPlotsEB3;
  TFile *fileSPlotsEE0;
  TFile *fileSPlotsEE1;
  TFile *fileSPlotsEE2;
  if (doLowPt==1) {
    fileSPlotsEB0 = new TFile("files_v2/probeLowPt/myFileSPlots_EB_0d5-to-1d5.root");       
    fileSPlotsEB1 = new TFile("files_v2/probeLowPt/myFileSPlots_EB_1d5-to-2d0.root");        
    fileSPlotsEB2 = new TFile("files_v2/probeLowPt/myFileSPlots_EB_2d0-to-5d0.root");        
    fileSPlotsEB3 = new TFile("files_v2/probeLowPt/myFileSPlots_EB_5d0-up.root");                 
    fileSPlotsEE0 = new TFile("files_v2/probeLowPt/myFileSPlots_EE_0d5-to-2d0.root");       
    fileSPlotsEE1 = new TFile("files_v2/probeLowPt/myFileSPlots_EE_2d0-to-5d0.root");        
    fileSPlotsEE2 = new TFile("files_v2/probeLowPt/myFileSPlots_EE_5d0-up.root");             
  }
  
  // MC histos
  TH1F *probeFBremSignalMc = (TH1F*)fileMC->Get("probeFBremSignalMc");
  TH1F *probeFBremFakeMc   = (TH1F*)fileMC->Get("probeFBremFakeMc");
  if (doLowPt==1) {
    probeFBremSignalMc = (TH1F*)fileMC->Get("probeFBremSignalMc");
    probeFBremFakeMc   = (TH1F*)fileMC->Get("probeFBremFakeMc");
  }    
  TH1F *probeMvaSignalMc    = (TH1F*)fileMC->Get("probeMvaSignalMc");
  TH1F *probeMvaFakeMc      = (TH1F*)fileMC->Get("probeMvaFakeMc");
  TH1F *probePtSignalMc     = (TH1F*)fileMC->Get("probePtSignalMc");
  TH1F *probePtFakeMc       = (TH1F*)fileMC->Get("probePtFakeMc");
  TH1F *probeEtaSignalMc    = (TH1F*)fileMC->Get("probeEtaSignalMc");
  TH1F *probeEtaFakeMc      = (TH1F*)fileMC->Get("probeEtaFakeMc");

  // MC histos in eta/pT bins
  TH1F *probeMvaSignalMcEB0;
  TH1F *probeMvaSignalMcEB1;
  TH1F *probeMvaSignalMcEB2;
  TH1F *probeMvaSignalMcEB3;
  TH1F *probeMvaSignalMcEE0;
  TH1F *probeMvaSignalMcEE1;
  TH1F *probeMvaSignalMcEE2;
  if (doLowPt==1) {
    // MC histos in eta/pT bins
    probeMvaSignalMcEB0 = (TH1F*)fileMC->Get("mvaSignalEBMc0");
    probeMvaSignalMcEB1 = (TH1F*)fileMC->Get("mvaSignalEBMc1");
    probeMvaSignalMcEB2 = (TH1F*)fileMC->Get("mvaSignalEBMc2");
    probeMvaSignalMcEB3 = (TH1F*)fileMC->Get("mvaSignalEBMc3");
    probeMvaSignalMcEE0 = (TH1F*)fileMC->Get("mvaSignalEEMc0");
    probeMvaSignalMcEE1 = (TH1F*)fileMC->Get("mvaSignalEEMc1");
    probeMvaSignalMcEE2 = (TH1F*)fileMC->Get("mvaSignalEEMc2");
  }

  // s-Plots output
  TH1F *probeFBremSignalData;
  TH1F *probeFBremBkgData;
  TH1F *probeMvaSignalData;
  TH1F *probeMvaBkgData;
  TH1F *probePtSignalData     = (TH1F*)fileSPlotsInclusive->Get("h1_probePt_jpsi__probePt");
  TH1F *probePtBkgData        = (TH1F*)fileSPlotsInclusive->Get("h1_probePt_bkg__probePt");
  TH1F *probeEtaSignalData    = (TH1F*)fileSPlotsInclusive->Get("h1_probeEta_jpsi__probeEta");
  TH1F *probeEtaBkgData       = (TH1F*)fileSPlotsInclusive->Get("h1_probeEta_bkg__probeEta");
  if (doLowPt==1) {
    probeMvaSignalData   = (TH1F*)fileSPlotsInclusive->Get("h1_probeMvaId_jpsi__probeMvaId");
    probeMvaBkgData      = (TH1F*)fileSPlotsInclusive->Get("h1_probeMvaId_bkg__probeMvaId");
    probeFBremSignalData = (TH1F*)fileSPlotsInclusive->Get("h1_probeFBrem_jpsi__probeFBrem");
    probeFBremBkgData    = (TH1F*)fileSPlotsInclusive->Get("h1_probeFBrem_bkg__probeFBrem");
  } else {
    probeMvaSignalData = (TH1F*)fileSPlotsInclusive->Get("h1_probePfmvaId_jpsi__probePfmvaId");
    probeMvaBkgData    = (TH1F*)fileSPlotsInclusive->Get("h1_probePfmvaId_bkg__probePfmvaId");
  }

  // s-Plots output in eta/pT bins
  TH1F *probeMvaSignalDataEB0;
  TH1F *probeMvaSignalDataEB1;
  TH1F *probeMvaSignalDataEB2;
  TH1F *probeMvaSignalDataEB3;
  TH1F *probeMvaSignalDataEE0;
  TH1F *probeMvaSignalDataEE1;
  TH1F *probeMvaSignalDataEE2;
  if (doLowPt==1) {
    // s-Plots output in eta/pT bins
    probeMvaSignalDataEB0 = (TH1F*)fileSPlotsEB0->Get("h1_probeMvaId_jpsi__probeMvaId");
    probeMvaSignalDataEB1 = (TH1F*)fileSPlotsEB1->Get("h1_probeMvaId_jpsi__probeMvaId");
    probeMvaSignalDataEB2 = (TH1F*)fileSPlotsEB2->Get("h1_probeMvaId_jpsi__probeMvaId");
    probeMvaSignalDataEB3 = (TH1F*)fileSPlotsEB3->Get("h1_probeMvaId_jpsi__probeMvaId");
    probeMvaSignalDataEE0 = (TH1F*)fileSPlotsEE0->Get("h1_probeMvaId_jpsi__probeMvaId");
    probeMvaSignalDataEE1 = (TH1F*)fileSPlotsEE1->Get("h1_probeMvaId_jpsi__probeMvaId");
    probeMvaSignalDataEE2 = (TH1F*)fileSPlotsEE2->Get("h1_probeMvaId_jpsi__probeMvaId");
  }

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
  if (doLowPt==1) {
    probeFBremSignalMc  -> SetLineWidth(2);  
    probeFBremFakeMc    -> SetLineWidth(2);  
    probeFBremSignalMc  -> SetLineColor(6);  
    probeFBremFakeMc    -> SetLineColor(4);  
  }
  //
  if (doLowPt==1) {
    probeMvaSignalMcEB0->SetLineWidth(2);  
    probeMvaSignalMcEB0->SetLineColor(6);  
    probeMvaSignalMcEB1->SetLineWidth(2);  
    probeMvaSignalMcEB1->SetLineColor(6);  
    probeMvaSignalMcEB2->SetLineWidth(2);  
    probeMvaSignalMcEB2->SetLineColor(6);  
    probeMvaSignalMcEB3->SetLineWidth(2);  
    probeMvaSignalMcEB3->SetLineColor(6);  
    probeMvaSignalMcEE0->SetLineWidth(2);  
    probeMvaSignalMcEE0->SetLineColor(6);  
    probeMvaSignalMcEE1->SetLineWidth(2);  
    probeMvaSignalMcEE1->SetLineColor(6);  
    probeMvaSignalMcEE2->SetLineWidth(2);  
    probeMvaSignalMcEE2->SetLineColor(6);  
  } 
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
  if (doLowPt==1) {
    probeFBremSignalData  -> SetLineWidth(2);  
    probeFBremBkgData     -> SetLineWidth(2);  
    probeFBremSignalData  -> SetLineColor(3);  
    probeFBremBkgData     -> SetLineColor(2);  
  }
  //
  if (doLowPt==1) {
    probeMvaSignalDataEB0 -> SetLineWidth(2);
    probeMvaSignalDataEB0 -> SetLineColor(3);
    probeMvaSignalDataEB1 -> SetLineWidth(2);
    probeMvaSignalDataEB1 -> SetLineColor(3);
    probeMvaSignalDataEB2 -> SetLineWidth(2);
    probeMvaSignalDataEB2 -> SetLineColor(3);
    probeMvaSignalDataEB3 -> SetLineWidth(2);
    probeMvaSignalDataEB3 -> SetLineColor(3);
    probeMvaSignalDataEE0 -> SetLineWidth(2);
    probeMvaSignalDataEE0 -> SetLineColor(3);
    probeMvaSignalDataEE1 -> SetLineWidth(2);
    probeMvaSignalDataEE1 -> SetLineColor(3);
    probeMvaSignalDataEE2 -> SetLineWidth(2);
    probeMvaSignalDataEE2 -> SetLineColor(3);
  }

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

  if (doLowPt==1) {
    TCanvas cprobeFBrems("cprobeFBrems","cprobeFBrems",1);
    probeFBremSignalData->DrawNormalized("hist");
    probeFBremSignalMc->DrawNormalized("samehist");
    leg3S->Draw();
    cprobeFBrems.SaveAs("probeFBrem_signal.png");
    // 
    TCanvas cprobeFBremb("cprobeFBremb","cprobeFBremb",1);
    probeFBremBkgData->DrawNormalized("hist");
    probeFBremFakeMc->DrawNormalized("samehist");
    leg3B->Draw();
    cprobeFBremb.SaveAs("probeFBrem_back.png");
  }


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

  // -----------------------------------------------
  if (doLowPt==1) {
    probeMvaSignalDataEB0->Rebin(2);
    probeMvaSignalDataEB1->Rebin(2);
    probeMvaSignalDataEB2->Rebin(2);
    probeMvaSignalDataEB3->Rebin(2);
    probeMvaSignalDataEE0->Rebin(2);
    probeMvaSignalDataEE1->Rebin(2);
    probeMvaSignalDataEE2->Rebin(2);
    //
    probeMvaSignalMcEB0->Rebin(2);
    probeMvaSignalMcEB1->Rebin(2);
    probeMvaSignalMcEB2->Rebin(2);
    probeMvaSignalMcEB3->Rebin(2);
    probeMvaSignalMcEE0->Rebin(2);
    probeMvaSignalMcEE1->Rebin(2);
    probeMvaSignalMcEE2->Rebin(2);
    // 
    drawTH1pair(probeMvaSignalDataEB0, probeMvaSignalMcEB0,  "ele ID BDT","a.u.",1.,"EB-0d5-1d5","./",2,"Data","MC");
    drawTH1pair(probeMvaSignalDataEB1, probeMvaSignalMcEB1,  "ele ID BDT","a.u.",1.,"EB-1d5-2d0","./",2,"Data","MC");
    drawTH1pair(probeMvaSignalDataEB2, probeMvaSignalMcEB2,  "ele ID BDT","a.u.",1.,"EB-2d0-5d0","./",2,"Data","MC");
    drawTH1pair(probeMvaSignalDataEB3, probeMvaSignalMcEB3,  "ele ID BDT","a.u.",1.,"EB-gt5","./",2,"Data","MC");
    //
    drawTH1pair(probeMvaSignalDataEE0, probeMvaSignalMcEE0,  "ele ID BDT","a.u.",1.,"EE-0g5-2d0","./",4,"Data","MC");
    drawTH1pair(probeMvaSignalDataEE1, probeMvaSignalMcEE1,  "ele ID BDT","a.u.",1.,"EE-2d0-5d0","./",4,"Data","MC");
    drawTH1pair(probeMvaSignalDataEE2, probeMvaSignalMcEE2,  "ele ID BDT","a.u.",1.,"EE-gt5","./",4,"Data","MC");
  }

}
