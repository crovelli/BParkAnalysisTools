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
		  const string& legEntry1 = "data", const string& legEntry2 = "MC", bool isleg1 = true, const string& ratioPadYaxisName = "data/MC", const string& outputFILE = "outFile.root") 
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

  TLegend leg2 (0.55,0.65,0.85,0.85);
  leg2.SetFillColor(0);
  leg2.SetFillStyle(0);
  leg2.SetBorderSize(0);
  leg2.AddEntry(h1,legEntry1.c_str(),"PLE");
  leg2.AddEntry(h2,legEntry2.c_str(),"L");

  if (isleg1)  leg.Draw("same");
  if (!isleg1) leg2.Draw("same");
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
  TH1D* den_noerr = (TH1D*) h2->Clone("den_noerr");
  TH1D* den = (TH1D*) h2->Clone("den");
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

  TFile fileOut(outputFILE.c_str(), "UPDATE");
  fileOut.cd();
  ratio->Write( ("ratio_" + canvasName).c_str());
  h2->Write( ("MC_" + canvasName).c_str());
  fileOut.Close();
}

void mcVsSplots(int doLowPt)
{
  // Input files (after TnP selection)
  TFile *fileMC;
  if (doLowPt==1) {
    fileMC = new TFile("files_marchWithReg2/probeLowPt/TightSel/myFileMcAfterTnp_TightSelection_noMcMatch_noKineWeights.root");         // from prepareInputsFromMcWithTnP
    //fileMC = new TFile("files_marchWithReg2/probeLowPt/TightSel/myFileMcAfterTnp_TightSelection_withMcMatch_noKineWeights.root");       // from prepareInputsFromMcWithTnP
  }
  if (doLowPt==0) {
    fileMC = new TFile("files_marchWithReg2/probePF/TightSel/myFileMcAfterTnp_noMcMatch_noKineWeights.root");         // from prepareInputsFromMcWithTnP
    //fileMC = new TFile("files_marchWithReg2/probePF/TightSel/myFileMcAfterTnp_withMcMatch_noKineWeights.root");         // from prepareInputsFromMcWithTnP
  }
  
  TFile *fileSPlotsEB0;
  TFile *fileSPlotsEB1;
  TFile *fileSPlotsEB2;
  TFile *fileSPlotsEB3;
  TFile *fileSPlotsEE0;
  TFile *fileSPlotsEE1;
  if (doLowPt==1) {
    fileSPlotsEB0 = new TFile("files_marchWithReg2/probeLowPt/TightSel/fileSPlots_EB_1d0-to-1d5.root");       
    fileSPlotsEB1 = new TFile("files_marchWithReg2/probeLowPt/TightSel/fileSPlots_EB_1d5-to-2d0.root");        
    fileSPlotsEB2 = new TFile("files_marchWithReg2/probeLowPt/TightSel/fileSPlots_EB_2d0-to-5d0.root");        
    fileSPlotsEB3 = new TFile("files_marchWithReg2/probeLowPt/TightSel/fileSPlots_EB_gt5d0.root");                 
    fileSPlotsEE0 = new TFile("files_marchWithReg2/probeLowPt/TightSel/fileSPlots_EE_2d0-to-5d0.root");       
    fileSPlotsEE1 = new TFile("files_marchWithReg2/probeLowPt/TightSel/fileSPlots_EE_gt5d0.root");        
  } else {
    fileSPlotsEB0 = new TFile("files_marchWithReg2/probePF/TightSel/fileSPlots_EB_2d0-to-5d0.root");        
    fileSPlotsEB1 = new TFile("files_marchWithReg2/probePF/TightSel/fileSPlots_EB_gt5d0.root");                 
    fileSPlotsEE0 = new TFile("files_marchWithReg2/probePF/TightSel/fileSPlots_EE_2d0-to-5d0.root");       
    fileSPlotsEE1 = new TFile("files_marchWithReg2/probePF/TightSel/fileSPlots_EE_gt5d0.root");        
  }
  
  // MC histos in eta/pT bins
  TH1F *probeMvaSignalMcEB0 = (TH1F*)fileMC->Get("mvaSignalEBMc0");
  TH1F *probeMvaSignalMcEB1 = (TH1F*)fileMC->Get("mvaSignalEBMc1");
  TH1F *probeMvaSignalMcEE0 = (TH1F*)fileMC->Get("mvaSignalEEMc0");
  TH1F *probeMvaSignalMcEE1 = (TH1F*)fileMC->Get("mvaSignalEEMc1");
  TH1F *probeMvaSignalMcEB2, *probeMvaSignalMcEB3;
  if (doLowPt==1) {
    probeMvaSignalMcEB2 = (TH1F*)fileMC->Get("mvaSignalEBMc2");
    probeMvaSignalMcEB3 = (TH1F*)fileMC->Get("mvaSignalEBMc3");
  }
  //
  TH1F *probeDxysigSignalMcEB0 = (TH1F*)fileMC->Get("dxysigSignalEBMc0");
  TH1F *probeDxysigSignalMcEB1 = (TH1F*)fileMC->Get("dxysigSignalEBMc1");
  TH1F *probeDxysigSignalMcEE0 = (TH1F*)fileMC->Get("dxysigSignalEEMc0");
  TH1F *probeDxysigSignalMcEE1 = (TH1F*)fileMC->Get("dxysigSignalEEMc1");
  TH1F *probeDxysigSignalMcEB2, *probeDxysigSignalMcEB3;
  if (doLowPt==1) {
    probeDxysigSignalMcEB2 = (TH1F*)fileMC->Get("dxysigSignalEBMc2");
    probeDxysigSignalMcEB3 = (TH1F*)fileMC->Get("dxysigSignalEBMc3");
  }  
  //
  TH1F *probeDztrgSignalMcEB0 = (TH1F*)fileMC->Get("dztrgSignalEBMc0");
  TH1F *probeDztrgSignalMcEB1 = (TH1F*)fileMC->Get("dztrgSignalEBMc1");
  TH1F *probeDztrgSignalMcEE0 = (TH1F*)fileMC->Get("dztrgSignalEEMc0");
  TH1F *probeDztrgSignalMcEE1 = (TH1F*)fileMC->Get("dztrgSignalEEMc1");
  TH1F *probeDztrgSignalMcEB2, *probeDztrgSignalMcEB3;
  if (doLowPt==1) {
    probeDztrgSignalMcEB2 = (TH1F*)fileMC->Get("dztrgSignalEBMc2");
    probeDztrgSignalMcEB3 = (TH1F*)fileMC->Get("dztrgSignalEBMc3");
  }
  //
  TH1F *probeIso04relSignalMcEB0 = (TH1F*)fileMC->Get("iso04relSignalEBMc0");
  TH1F *probeIso04relSignalMcEB1 = (TH1F*)fileMC->Get("iso04relSignalEBMc1");
  TH1F *probeIso04relSignalMcEE0 = (TH1F*)fileMC->Get("iso04relSignalEEMc0");
  TH1F *probeIso04relSignalMcEE1 = (TH1F*)fileMC->Get("iso04relSignalEEMc1");
  TH1F *probeIso04relSignalMcEB2, *probeIso04relSignalMcEB3;
  if (doLowPt==1) {
    probeIso04relSignalMcEB2 = (TH1F*)fileMC->Get("iso04relSignalEBMc2");
    probeIso04relSignalMcEB3 = (TH1F*)fileMC->Get("iso04relSignalEBMc3");
  }

  // s-Plots output in eta/pT bins
  TH1F *probeMvaSignalDataEB0;
  TH1F *probeMvaSignalDataEB1;
  TH1F *probeMvaSignalDataEB2;
  TH1F *probeMvaSignalDataEB3;
  TH1F *probeMvaSignalDataEE0;
  TH1F *probeMvaSignalDataEE1;
  if (doLowPt==1) {
    probeMvaSignalDataEB0 = (TH1F*)fileSPlotsEB0->Get("h1_probeMvaId_jpsi__probeMvaId");
    probeMvaSignalDataEB1 = (TH1F*)fileSPlotsEB1->Get("h1_probeMvaId_jpsi__probeMvaId");
    probeMvaSignalDataEB2 = (TH1F*)fileSPlotsEB2->Get("h1_probeMvaId_jpsi__probeMvaId");
    probeMvaSignalDataEB3 = (TH1F*)fileSPlotsEB3->Get("h1_probeMvaId_jpsi__probeMvaId");
    probeMvaSignalDataEE0 = (TH1F*)fileSPlotsEE0->Get("h1_probeMvaId_jpsi__probeMvaId");
    probeMvaSignalDataEE1 = (TH1F*)fileSPlotsEE1->Get("h1_probeMvaId_jpsi__probeMvaId");
  } else {
    probeMvaSignalDataEB0 = (TH1F*)fileSPlotsEB0->Get("h1_probePfmvaId_jpsi__probePfmvaId");
    probeMvaSignalDataEB1 = (TH1F*)fileSPlotsEB1->Get("h1_probePfmvaId_jpsi__probePfmvaId");
    probeMvaSignalDataEE0 = (TH1F*)fileSPlotsEE0->Get("h1_probePfmvaId_jpsi__probePfmvaId");
    probeMvaSignalDataEE1 = (TH1F*)fileSPlotsEE1->Get("h1_probePfmvaId_jpsi__probePfmvaId");
  }
  //
  TH1F *probeDxySigSignalDataEB0 = (TH1F*)fileSPlotsEB0->Get("h1_probeDxySig_jpsi__probeDxySig");
  TH1F *probeDxySigSignalDataEB1 = (TH1F*)fileSPlotsEB1->Get("h1_probeDxySig_jpsi__probeDxySig");
  TH1F *probeDxySigSignalDataEE0 = (TH1F*)fileSPlotsEE0->Get("h1_probeDxySig_jpsi__probeDxySig");
  TH1F *probeDxySigSignalDataEE1 = (TH1F*)fileSPlotsEE1->Get("h1_probeDxySig_jpsi__probeDxySig");
  TH1F *probeDxySigSignalDataEB2, *probeDxySigSignalDataEB3;
  if (doLowPt==1) {
    probeDxySigSignalDataEB2 = (TH1F*)fileSPlotsEB2->Get("h1_probeDxySig_jpsi__probeDxySig"); 
    probeDxySigSignalDataEB3 = (TH1F*)fileSPlotsEB3->Get("h1_probeDxySig_jpsi__probeDxySig"); 
  }
  //
  TH1F *probeDzTrgSignalDataEB0 = (TH1F*)fileSPlotsEB0->Get("h1_probeDzTrg_jpsi__probeDzTrg");
  TH1F *probeDzTrgSignalDataEB1 = (TH1F*)fileSPlotsEB1->Get("h1_probeDzTrg_jpsi__probeDzTrg");
  TH1F *probeDzTrgSignalDataEE0 = (TH1F*)fileSPlotsEE0->Get("h1_probeDzTrg_jpsi__probeDzTrg");
  TH1F *probeDzTrgSignalDataEE1 = (TH1F*)fileSPlotsEE1->Get("h1_probeDzTrg_jpsi__probeDzTrg");
  TH1F *probeDzTrgSignalDataEB2, *probeDzTrgSignalDataEB3;
  if (doLowPt==1) { 
    probeDzTrgSignalDataEB2 = (TH1F*)fileSPlotsEB2->Get("h1_probeDzTrg_jpsi__probeDzTrg");
    probeDzTrgSignalDataEB3 = (TH1F*)fileSPlotsEB3->Get("h1_probeDzTrg_jpsi__probeDzTrg");
  }
  //
  TH1F *probeIso04RelSignalDataEB0 = (TH1F*)fileSPlotsEB0->Get("h1_probeIso04Rel_jpsi__probeIso04Rel");
  TH1F *probeIso04RelSignalDataEB1 = (TH1F*)fileSPlotsEB1->Get("h1_probeIso04Rel_jpsi__probeIso04Rel");
  TH1F *probeIso04RelSignalDataEE0 = (TH1F*)fileSPlotsEE0->Get("h1_probeIso04Rel_jpsi__probeIso04Rel");
  TH1F *probeIso04RelSignalDataEE1 = (TH1F*)fileSPlotsEE1->Get("h1_probeIso04Rel_jpsi__probeIso04Rel");
  TH1F *probeIso04RelSignalDataEB2, *probeIso04RelSignalDataEB3;
  if (doLowPt==1) {
    probeIso04RelSignalDataEB2 = (TH1F*)fileSPlotsEB2->Get("h1_probeIso04Rel_jpsi__probeIso04Rel");
    probeIso04RelSignalDataEB3 = (TH1F*)fileSPlotsEB3->Get("h1_probeIso04Rel_jpsi__probeIso04Rel");
  }

  // Plots
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);


  // Id
  probeMvaSignalDataEB0->Rebin(2);
  probeMvaSignalMcEB0->Rebin(2);
  if (doLowPt==1) drawTH1pair(probeMvaSignalDataEB0, probeMvaSignalMcEB0,  "ele ID BDT","a.u.",1.,"Id_EB-1d0-1d5","./",2,"Data","MC",0);
  else            drawTH1pair(probeMvaSignalDataEB0, probeMvaSignalMcEB0,  "ele ID BDT","a.u.",1.,"Id_EB-2d0-5d0","./",2,"Data","MC",1);

  probeMvaSignalDataEB1->Rebin(2);
  probeMvaSignalMcEB1->Rebin(2);
  if (doLowPt==1) drawTH1pair(probeMvaSignalDataEB1, probeMvaSignalMcEB1,  "ele ID BDT","a.u.",1.,"Id_EB-1d5-2d0","./",2,"Data","MC",0);
  else            drawTH1pair(probeMvaSignalDataEB1, probeMvaSignalMcEB1,  "ele ID BDT","a.u.",1.,"Id_EB-gt5d0","./",2,"Data","MC",1);

  if (doLowPt==1) {
    probeMvaSignalDataEB2->Rebin(2);
    probeMvaSignalMcEB2->Rebin(2);
    drawTH1pair(probeMvaSignalDataEB2, probeMvaSignalMcEB2,  "ele ID BDT","a.u.",1.,"Id_EB-2d0-5d0","./",2,"Data","MC",0);

    probeMvaSignalDataEB3->Rebin(2);
    probeMvaSignalMcEB3->Rebin(2);
    drawTH1pair(probeMvaSignalDataEB3, probeMvaSignalMcEB3,  "ele ID BDT","a.u.",1.,"Id_EB-gt5d0","./",2,"Data","MC",0);
  }

  probeMvaSignalDataEE0->Rebin(4);
  probeMvaSignalMcEE0->Rebin(4);
  drawTH1pair(probeMvaSignalDataEE0, probeMvaSignalMcEE0,  "ele ID BDT","a.u.",1.,"Id_EE-2g0-5d0","./",4,"Data","MC");

  probeMvaSignalDataEE1->Rebin(4);
  probeMvaSignalMcEE1->Rebin(4);
  drawTH1pair(probeMvaSignalDataEE1, probeMvaSignalMcEE1,  "ele ID BDT","a.u.",1.,"Id_EE-gt5d0","./",4,"Data","MC");
}
