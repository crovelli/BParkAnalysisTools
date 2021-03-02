#define finalRocs_cxx

// ROOT includes 
#include <TROOT.h>
#include <TStyle.h>
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <iostream>

// Utilities
#include "./FiguresOfMeritEvaluator.cc"

using namespace std;

void finalRocs(bool isLowPt)
{
  // Files
  TFile fileSPlotsData("files_v2/probeLowPt/myFileSPlots_fullRange.root");                           // TnP selection, sPlots             (jpsi_splot.C)
  TFile fileTnpMc("files_v2/probeLowPt/myFileMcAfterTnp__noPUweight__withTnpVsNaniWeight.root");     // TnP selection, match to MC-truth  (prepareInputsFromMcWithTnP.C)
  TFile fileFakeSelData("files_v2/probeLowPt/myFileFakesDataFull2018.root");                         // Fakes selection applied to data   (prepareInputsFromFakes.C)
  TFile fileFakeSelMc("files_v2/probeLowPt/myFileFakesMC_noPuWeight__withFakeVsNaniWeight.root");    // Fakes selection applied to MC     (prepareInputsFromFakes.C)   
  TFile fileNani("files_v2/probeLowPt/myFileFromNani__noPUweight.root");                             // Match to MC-truth                 (prepareInputsFromNaniInMc.C)

  // Histos
  TH1F *mvaSignalData;
  if (isLowPt==1) mvaSignalData = (TH1F*)fileSPlotsData.Get("h1_probeMvaId_jpsi__probeMvaId");
  if (isLowPt==0) mvaSignalData = (TH1F*)fileSPlotsData.Get("h1_probePfmvaId_jpsi__probePfmvaId");
  TH1F *mvaBkgData    = (TH1F*)fileFakeSelData.Get("mvaFakes"); 
  TH1F *mvaSignalMc   = (TH1F*)fileTnpMc.Get("probeMvaSignalMc");
  TH1F *mvaBkgMc      = (TH1F*)fileFakeSelMc.Get("mvaFakes"); 
  TH1F *mvaSignalMcWW = (TH1F*)fileTnpMc.Get("probeMvaSignalMcWW");
  TH1F *mvaBkgMcWW    = (TH1F*)fileFakeSelMc.Get("mvaFakesWW"); 
  TH1F *mvaSignalNani = (TH1F*)fileNani.Get("mvaSignalMc");
  TH1F *mvaBkgNani    = (TH1F*)fileNani.Get("mvaFakeMc");
  mvaSignalData->Sumw2();   
  mvaBkgData->Sumw2();   
  mvaSignalMc->Sumw2();   
  mvaBkgMc->Sumw2();   
  mvaSignalMcWW->Sumw2();   
  mvaBkgMcWW->Sumw2();   
  mvaSignalNani->Sumw2();
  mvaBkgNani->Sumw2();

  // Histos in bin
  TH1F *mvaBkgDataEB0;
  TH1F *mvaBkgMcEB0;
  TH1F *mvaSignalMcEB0;
  TH1F *mvaBkgDataEB1;
  TH1F *mvaBkgMcEB1;
  TH1F *mvaSignalMcEB1;
  TH1F *mvaBkgDataEB2;
  TH1F *mvaBkgMcEB2;
  TH1F *mvaSignalMcEB2;
  TH1F *mvaBkgDataEB3;
  TH1F *mvaBkgMcEB3;
  TH1F *mvaSignalMcEB3;
  TH1F *mvaBkgDataEE0;
  TH1F *mvaBkgMcEE0;
  TH1F *mvaSignalMcEE0;
  TH1F *mvaBkgDataEE1;
  TH1F *mvaBkgMcEE1;
  TH1F *mvaSignalMcEE1;
  TH1F *mvaBkgDataEE2;
  TH1F *mvaBkgMcEE2;
  TH1F *mvaSignalMcEE2;
  if (isLowPt==1) {
    mvaBkgDataEB0  = (TH1F*)fileFakeSelData.Get("mvaFakeEB0"); 
    mvaBkgMcEB0    = (TH1F*)fileFakeSelMc.Get("mvaFakeEB0"); 
    mvaSignalMcEB0 = (TH1F*)fileTnpMc.Get("mvaSignalEBMc0");
    mvaBkgDataEB1  = (TH1F*)fileFakeSelData.Get("mvaFakeEB1"); 
    mvaBkgMcEB1    = (TH1F*)fileFakeSelMc.Get("mvaFakeEB1"); 
    mvaSignalMcEB1 = (TH1F*)fileTnpMc.Get("mvaSignalEBMc1");
    mvaBkgDataEB2  = (TH1F*)fileFakeSelData.Get("mvaFakeEB2"); 
    mvaBkgMcEB2    = (TH1F*)fileFakeSelMc.Get("mvaFakeEB2"); 
    mvaSignalMcEB2 = (TH1F*)fileTnpMc.Get("mvaSignalEBMc2");
    mvaBkgDataEB3  = (TH1F*)fileFakeSelData.Get("mvaFakeEB3"); 
    mvaBkgMcEB3    = (TH1F*)fileFakeSelMc.Get("mvaFakeEB3"); 
    mvaSignalMcEB3 = (TH1F*)fileTnpMc.Get("mvaSignalEBMc3");
    mvaBkgDataEE0  = (TH1F*)fileFakeSelData.Get("mvaFakeEE0"); 
    mvaBkgMcEE0    = (TH1F*)fileFakeSelMc.Get("mvaFakeEE0"); 
    mvaSignalMcEE0 = (TH1F*)fileTnpMc.Get("mvaSignalEEMc0");
    mvaBkgDataEE1  = (TH1F*)fileFakeSelData.Get("mvaFakeEE1"); 
    mvaBkgMcEE1    = (TH1F*)fileFakeSelMc.Get("mvaFakeEE1"); 
    mvaSignalMcEE1 = (TH1F*)fileTnpMc.Get("mvaSignalEEMc1");
    mvaBkgDataEE2  = (TH1F*)fileFakeSelData.Get("mvaFakeEE2"); 
    mvaBkgMcEE2    = (TH1F*)fileFakeSelMc.Get("mvaFakeEE2"); 
    mvaSignalMcEE2 = (TH1F*)fileTnpMc.Get("mvaSignalEEMc2");
  }

  // Cosmetics
  mvaSignalData -> SetLineWidth(2);  
  mvaSignalData -> SetLineColor(1);  
  mvaBkgData    -> SetLineWidth(2);  
  mvaBkgData    -> SetLineColor(1);  
  mvaSignalMc   -> SetLineWidth(2);  
  mvaSignalMc   -> SetLineColor(4);  
  mvaBkgMc      -> SetLineWidth(2);  
  mvaBkgMc      -> SetLineColor(4);  

  // Plots
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);

  TLegend *lega;
  lega = new TLegend(0.20,0.65,0.50,0.90);
  lega->SetFillStyle(0);
  lega->SetBorderSize(0);
  lega->SetTextSize(0.05);
  lega->SetFillColor(0);
  lega->AddEntry(mvaBkgMc,   "MC", "lp");
  lega->AddEntry(mvaBkgData, "Data", "lp");

  TLegend *legb;
  legb = new TLegend(0.65,0.65,0.90,0.90);
  legb->SetFillStyle(0);
  legb->SetBorderSize(0);
  legb->SetTextSize(0.05);
  legb->SetFillColor(0);
  legb->AddEntry(mvaBkgMc,   "MC", "lp");
  legb->AddEntry(mvaBkgData, "Data", "lp");

  TCanvas cmvas("cmvas","cmvas",1);
  mvaSignalMc->DrawNormalized("hist");
  mvaSignalData->DrawNormalized("samehist");
  lega->Draw();
  cmvas.SaveAs("signalCheck.png");

  TCanvas cmvab("cmvab","cmvab",1);
  mvaBkgMc->DrawNormalized("hist");
  mvaBkgData->DrawNormalized("samehist");
  legb->Draw();
  cmvab.SaveAs("bkgCheck.png");


  // ----------------------------------------------------  
  // Output file with ROCs
  TFile fileROC("outputROCs.root","RECREATE");  
  fileROC.cd();

  // ----------------------------------------------------
  // FOMs
  FiguresOfMeritEvaluator myRocBDTmc;
  myRocBDTmc.addSignal("BDT", mvaSignalMc);
  myRocBDTmc.addBackgrounds(mvaBkgMc);
  myRocBDTmc.setCutDirection(">");
  TGraph *myGraphBDTmc= myRocBDTmc.getFOM("BDT",2);
  myGraphBDTmc->SetTitle("BDT");
  myGraphBDTmc->SetName("BDT");
  //
  FiguresOfMeritEvaluator myRocBDTmcww;
  myRocBDTmcww.addSignal("BDT", mvaSignalMcWW);
  myRocBDTmcww.addBackgrounds(mvaBkgMcWW);
  myRocBDTmcww.setCutDirection(">");
  TGraph *myGraphBDTmcww= myRocBDTmcww.getFOM("BDT",2);
  myGraphBDTmcww->SetTitle("BDT");
  myGraphBDTmcww->SetName("BDT");
  //
  FiguresOfMeritEvaluator myRocBDTdata;
  myRocBDTdata.addSignal("BDT", mvaSignalMc);      // use mc for signal
  myRocBDTdata.addBackgrounds(mvaBkgData);
  myRocBDTdata.setCutDirection(">");
  TGraph *myGraphBDTdata= myRocBDTdata.getFOM("BDT",2);
  myGraphBDTdata->SetTitle("BDT");
  myGraphBDTdata->SetName("BDT");
  //
  FiguresOfMeritEvaluator myRocBDTnani;
  myRocBDTnani.addSignal("BDT", mvaSignalNani);
  myRocBDTnani.addBackgrounds(mvaBkgNani);
  myRocBDTnani.setCutDirection(">");
  TGraph *myGraphBDTnani= myRocBDTnani.getFOM("BDT",2);
  myGraphBDTnani->SetTitle("BDT");
  myGraphBDTnani->SetName("BDT");

  // ----------------------------------------------------
  // Cosmetics
  myGraphBDTdata->SetMarkerColor(1);
  myGraphBDTdata->SetMarkerStyle(20);
  myGraphBDTdata->SetMarkerSize(1);
  myGraphBDTmc->SetMarkerColor(4);
  myGraphBDTmc->SetMarkerStyle(20);
  myGraphBDTmc->SetMarkerSize(1);
  myGraphBDTmcww->SetMarkerColor(2);
  myGraphBDTmcww->SetMarkerStyle(20);
  myGraphBDTmcww->SetMarkerSize(1);
  myGraphBDTnani->SetMarkerColor(3);
  myGraphBDTnani->SetMarkerStyle(20);
  myGraphBDTnani->SetMarkerSize(1);

  // ----------------------------------------------------
  // Plots
  TH2F *myH = new TH2F("myH","",100, 0., 1., 100, 0.,1.);
  myH->GetXaxis()->SetTitle("Mistag Rate");
  myH->GetYaxis()->SetTitle("Efficiency");
  TH2F *myHz = new TH2F("myHz","",100, 0., 1., 100, 0.7,1.);
  myHz->GetXaxis()->SetTitle("Mistag Rate");
  myHz->GetYaxis()->SetTitle("Efficiency");
  TH2F *myH2 = new TH2F("myH2","",100, 0.0001, 1., 100, 0.,1.);
  myH2->GetXaxis()->SetTitle("Mistag Rate");
  myH2->GetYaxis()->SetTitle("Efficiency");

  // Data vs Mc with same selection (tnp and fakes)
  TLegend *leg2;
  leg2 = new TLegend(0.50,0.25,0.85,0.60);
  leg2->SetFillStyle(0);
  leg2->SetBorderSize(0);
  leg2->SetTextSize(0.05);
  leg2->SetFillColor(0);
  leg2->AddEntry(myGraphBDTdata, "Data B, MC S", "lp");
  leg2->AddEntry(myGraphBDTmc,   "From MC", "lp");
  //
  TCanvas croc("roc","",1);
  croc.SetGrid();
  myHz->Draw();
  myGraphBDTdata->Draw("sameP");
  myGraphBDTmc->Draw("sameP");
  leg2->Draw();
  croc.SaveAs("roc_dataVsMc.png");  
  myH2->Draw();
  myGraphBDTdata->Draw("sameP");
  myGraphBDTmc->Draw("sameP");
  leg2->Draw();
  croc.SetLogx();
  croc.SaveAs("roc_dataVsMc_log.png");  


  fileROC.cd();
  myGraphBDTdata->Write("myGraphBDTdata");

  // -----------------------------------
  // MC after selection vs Mc-from-nani
  TLegend *leg3;
  leg3 = new TLegend(0.50,0.25,0.85,0.60);
  leg3->SetFillStyle(0);
  leg3->SetBorderSize(0);
  leg3->SetTextSize(0.05);
  leg3->SetFillColor(0);
  leg3->AddEntry(myGraphBDTnani, "MC-truth", "lp");
  leg3->AddEntry(myGraphBDTmc,   "MC after sel", "lp");
  //
  TCanvas croc3("roc","",1);
  croc3.SetGrid();
  myHz->Draw();
  myGraphBDTmc->Draw("sameP");
  myGraphBDTnani->Draw("sameP");
  leg3->Draw();
  croc3.SaveAs("roc_mcVsNani_noweight.png");  
  myH2->Draw();
  myGraphBDTmc->Draw("sameP");
  myGraphBDTnani->Draw("sameP");
  leg3->Draw();
  croc3.SetLogx();
  croc3.SaveAs("roc_mcVsNani_noweight_log.png");  

  // -----------------------------------
  // MC after selection vs Mc-from-nani with weights
  TLegend *leg4;
  leg4 = new TLegend(0.50,0.25,0.85,0.60);
  leg4->SetFillStyle(0);
  leg4->SetBorderSize(0);
  leg4->SetTextSize(0.05);
  leg4->SetFillColor(0);
  leg4->AddEntry(myGraphBDTnani, "MC-truth", "lp");
  leg4->AddEntry(myGraphBDTmcww, "MC after sel", "lp");
  //
  TCanvas croc4("roc","",1);
  croc4.SetGrid();
  myHz->Draw();
  myGraphBDTmcww->Draw("sameP");
  myGraphBDTnani->Draw("sameP");
  leg4->Draw();
  croc4.SaveAs("roc_mcVsNani_withweight.png");  

  myH2->Draw();
  myGraphBDTmcww->Draw("sameP");
  myGraphBDTnani->Draw("sameP");
  leg4->Draw();
  croc4.SetLogx();
  croc4.SaveAs("roc_mcVsNani_withweight_log.png");  


  // ---------------------------------------------------
  if (isLowPt==1) {
    FiguresOfMeritEvaluator myRocBDTdataEB0;
    myRocBDTdataEB0.addSignal("BDT", mvaSignalMcEB0);      
    myRocBDTdataEB0.addBackgrounds(mvaBkgDataEB0);
    myRocBDTdataEB0.setCutDirection(">");
    TGraph *myGraphBDTdataEB0 = myRocBDTdataEB0.getFOM("BDT",2);
    myGraphBDTdataEB0->SetTitle("BDT");
    myGraphBDTdataEB0->SetName("BDT");
    myGraphBDTdataEB0->SetMarkerColor(1);
    myGraphBDTdataEB0->SetMarkerStyle(20);
    myGraphBDTdataEB0->SetMarkerSize(1);
    //
    FiguresOfMeritEvaluator myRocBDTmcEB0;
    myRocBDTmcEB0.addSignal("BDT", mvaSignalMcEB0);      
    myRocBDTmcEB0.addBackgrounds(mvaBkgMcEB0);
    myRocBDTmcEB0.setCutDirection(">");
    TGraph *myGraphBDTmcEB0 = myRocBDTmcEB0.getFOM("BDT",2);
    myGraphBDTmcEB0->SetTitle("BDT");
    myGraphBDTmcEB0->SetName("BDT");
    myGraphBDTmcEB0->SetMarkerColor(4);
    myGraphBDTmcEB0->SetMarkerStyle(20);
    myGraphBDTmcEB0->SetMarkerSize(1);
    //
    TCanvas croceb0("roc","",1);
    croceb0.SetGrid();
    myH->Draw();
    myGraphBDTdataEB0->Draw("sameP");
    myGraphBDTmcEB0->Draw("sameP");
    leg2->Draw();
    croceb0.SaveAs("roc_EB0_dataVsMc.png");  
    myH2->Draw();
    myGraphBDTdataEB0->Draw("sameP");
    myGraphBDTmcEB0->Draw("sameP");
    leg2->Draw();
    croceb0.SetLogx();
    croceb0.SaveAs("roc_EB0_dataVsMc_log.png");  

    fileROC.cd();
    myGraphBDTdataEB0->Write("myGraphBDTdataEB0");
    
    // ---------------------------------------------------
    FiguresOfMeritEvaluator myRocBDTdataEB1;
    myRocBDTdataEB1.addSignal("BDT", mvaSignalMcEB1);      
    myRocBDTdataEB1.addBackgrounds(mvaBkgDataEB1);
    myRocBDTdataEB1.setCutDirection(">");
    TGraph *myGraphBDTdataEB1 = myRocBDTdataEB1.getFOM("BDT",2);
    myGraphBDTdataEB1->SetTitle("BDT");
    myGraphBDTdataEB1->SetName("BDT");
    myGraphBDTdataEB1->SetMarkerColor(1);
    myGraphBDTdataEB1->SetMarkerStyle(20);
    myGraphBDTdataEB1->SetMarkerSize(1);
    //
    FiguresOfMeritEvaluator myRocBDTmcEB1;
    myRocBDTmcEB1.addSignal("BDT", mvaSignalMcEB1);      
    myRocBDTmcEB1.addBackgrounds(mvaBkgMcEB1);
    myRocBDTmcEB1.setCutDirection(">");
    TGraph *myGraphBDTmcEB1 = myRocBDTmcEB1.getFOM("BDT",2);
    myGraphBDTmcEB1->SetTitle("BDT");
    myGraphBDTmcEB1->SetName("BDT");
    myGraphBDTmcEB1->SetMarkerColor(4);
    myGraphBDTmcEB1->SetMarkerStyle(20);
    myGraphBDTmcEB1->SetMarkerSize(1);
    //
    TCanvas croceb1("roc","",1);
    croceb1.SetGrid();
    myH->Draw();
    myGraphBDTdataEB1->Draw("sameP");
    myGraphBDTmcEB1->Draw("sameP");
    leg2->Draw();
    croceb1.SaveAs("roc_EB1_dataVsMc.png");  
    myH2->Draw();
    myGraphBDTdataEB1->Draw("sameP");
    myGraphBDTmcEB1->Draw("sameP");
    leg2->Draw();
    croceb1.SetLogx();
    croceb1.SaveAs("roc_EB1_dataVsMc_log.png");  

    fileROC.cd();
    myGraphBDTdataEB1->Write("myGraphBDTdataEB1");
    
    // ---------------------------------------------------
    FiguresOfMeritEvaluator myRocBDTdataEB2;
    myRocBDTdataEB2.addSignal("BDT", mvaSignalMcEB1);      
    myRocBDTdataEB2.addBackgrounds(mvaBkgDataEB2);
    myRocBDTdataEB2.setCutDirection(">");
    TGraph *myGraphBDTdataEB2 = myRocBDTdataEB2.getFOM("BDT",2);
    myGraphBDTdataEB2->SetTitle("BDT");
    myGraphBDTdataEB2->SetName("BDT");
    myGraphBDTdataEB2->SetMarkerColor(1);
    myGraphBDTdataEB2->SetMarkerStyle(20);
    myGraphBDTdataEB2->SetMarkerSize(1);
    //
    FiguresOfMeritEvaluator myRocBDTmcEB2;
    myRocBDTmcEB2.addSignal("BDT", mvaSignalMcEB2);      
    myRocBDTmcEB2.addBackgrounds(mvaBkgMcEB2);
    myRocBDTmcEB2.setCutDirection(">");
    TGraph *myGraphBDTmcEB2 = myRocBDTmcEB2.getFOM("BDT",2);
    myGraphBDTmcEB2->SetTitle("BDT");
    myGraphBDTmcEB2->SetName("BDT");
    myGraphBDTmcEB2->SetMarkerColor(4);
    myGraphBDTmcEB2->SetMarkerStyle(20);
    myGraphBDTmcEB2->SetMarkerSize(1);
    //
    TCanvas croceb2("roc","",1);
    croceb2.SetGrid();
    myH->Draw();
    myGraphBDTdataEB2->Draw("sameP");
    myGraphBDTmcEB2->Draw("sameP");
    leg2->Draw();
    croceb2.SaveAs("roc_EB2_dataVsMc.png");  
    myH2->Draw();
    myGraphBDTdataEB2->Draw("sameP");
    myGraphBDTmcEB2->Draw("sameP");
    leg2->Draw();
    croceb2.SetLogx();
    croceb2.SaveAs("roc_EB2_dataVsMc_log.png");  

    fileROC.cd();
    myGraphBDTdataEB2->Write("myGraphBDTdataEB2");
    
    // ---------------------------------------------------
    FiguresOfMeritEvaluator myRocBDTdataEB3;
    myRocBDTdataEB3.addSignal("BDT", mvaSignalMcEB3);      
    myRocBDTdataEB3.addBackgrounds(mvaBkgDataEB3);
    myRocBDTdataEB3.setCutDirection(">");
    TGraph *myGraphBDTdataEB3 = myRocBDTdataEB3.getFOM("BDT",2);
    myGraphBDTdataEB3->SetTitle("BDT");
    myGraphBDTdataEB3->SetName("BDT");
    myGraphBDTdataEB3->SetMarkerColor(1);
    myGraphBDTdataEB3->SetMarkerStyle(20);
    myGraphBDTdataEB3->SetMarkerSize(1);
    //
    FiguresOfMeritEvaluator myRocBDTmcEB3;
    myRocBDTmcEB3.addSignal("BDT", mvaSignalMcEB3);      
    myRocBDTmcEB3.addBackgrounds(mvaBkgMcEB3);
    myRocBDTmcEB3.setCutDirection(">");
    TGraph *myGraphBDTmcEB3 = myRocBDTmcEB3.getFOM("BDT",2);
    myGraphBDTmcEB3->SetTitle("BDT");
    myGraphBDTmcEB3->SetName("BDT");
    myGraphBDTmcEB3->SetMarkerColor(4);
    myGraphBDTmcEB3->SetMarkerStyle(20);
    myGraphBDTmcEB3->SetMarkerSize(1);
    //
    TCanvas croceb3("roc","",1);
    croceb3.SetGrid();
    myH->Draw();
    myGraphBDTdataEB3->Draw("sameP");
    myGraphBDTmcEB3->Draw("sameP");
    leg2->Draw();
    croceb3.SaveAs("roc_EB3_dataVsMc.png");  
    myH2->Draw();
    myGraphBDTdataEB3->Draw("sameP");
    myGraphBDTmcEB3->Draw("sameP");
    leg2->Draw();
    croceb3.SetLogx();
    croceb3.SaveAs("roc_EB3_dataVsMc_log.png");  

    fileROC.cd();
    myGraphBDTdataEB3->Write("myGraphBDTdataEB3");
    
    // ---------------------------------------------------
    FiguresOfMeritEvaluator myRocBDTdataEE0;
    myRocBDTdataEE0.addSignal("BDT", mvaSignalMcEE0);      
    myRocBDTdataEE0.addBackgrounds(mvaBkgDataEE0);
    myRocBDTdataEE0.setCutDirection(">");
    TGraph *myGraphBDTdataEE0 = myRocBDTdataEE0.getFOM("BDT",2);
    myGraphBDTdataEE0->SetTitle("BDT");
    myGraphBDTdataEE0->SetName("BDT");
    myGraphBDTdataEE0->SetMarkerColor(1);
    myGraphBDTdataEE0->SetMarkerStyle(20);
    myGraphBDTdataEE0->SetMarkerSize(1);
    //
    FiguresOfMeritEvaluator myRocBDTmcEE0;
    myRocBDTmcEE0.addSignal("BDT", mvaSignalMcEE0);      
    myRocBDTmcEE0.addBackgrounds(mvaBkgMcEE0);
    myRocBDTmcEE0.setCutDirection(">");
    TGraph *myGraphBDTmcEE0 = myRocBDTmcEE0.getFOM("BDT",2);
    myGraphBDTmcEE0->SetTitle("BDT");
    myGraphBDTmcEE0->SetName("BDT");
    myGraphBDTmcEE0->SetMarkerColor(4);
    myGraphBDTmcEE0->SetMarkerStyle(20);
    myGraphBDTmcEE0->SetMarkerSize(1);
    //
    TCanvas crocee0("roc","",1);
    crocee0.SetGrid();
    myH->Draw();
    myGraphBDTdataEE0->Draw("sameP");
    myGraphBDTmcEE0->Draw("sameP");
    leg2->Draw();
    crocee0.SaveAs("roc_EE0_dataVsMc.png");  
    myH2->Draw();
    myGraphBDTdataEE0->Draw("sameP");
    myGraphBDTmcEE0->Draw("sameP");
    leg2->Draw();
    crocee0.SetLogx();
    crocee0.SaveAs("roc_EE0_dataVsMc_log.png");  

    fileROC.cd();
    myGraphBDTdataEE0->Write("myGraphBDTdataEE0");
    
    // ---------------------------------------------------
    FiguresOfMeritEvaluator myRocBDTdataEE1;
    myRocBDTdataEE1.addSignal("BDT", mvaSignalMcEE1);      
    myRocBDTdataEE1.addBackgrounds(mvaBkgDataEE1);
    myRocBDTdataEE1.setCutDirection(">");
    TGraph *myGraphBDTdataEE1 = myRocBDTdataEE1.getFOM("BDT",2);
    myGraphBDTdataEE1->SetTitle("BDT");
    myGraphBDTdataEE1->SetName("BDT");
    myGraphBDTdataEE1->SetMarkerColor(1);
    myGraphBDTdataEE1->SetMarkerStyle(20);
    myGraphBDTdataEE1->SetMarkerSize(1);
    //
    FiguresOfMeritEvaluator myRocBDTmcEE1;
    myRocBDTmcEE1.addSignal("BDT", mvaSignalMcEE1);      
    myRocBDTmcEE1.addBackgrounds(mvaBkgMcEE1);
    myRocBDTmcEE1.setCutDirection(">");
    TGraph *myGraphBDTmcEE1 = myRocBDTmcEE1.getFOM("BDT",2);
    myGraphBDTmcEE1->SetTitle("BDT");
    myGraphBDTmcEE1->SetName("BDT");
    myGraphBDTmcEE1->SetMarkerColor(4);
    myGraphBDTmcEE1->SetMarkerStyle(20);
    myGraphBDTmcEE1->SetMarkerSize(1);
    //
    TCanvas crocee1("roc","",1);
    crocee1.SetGrid();
    myH->Draw();
    myGraphBDTdataEE1->Draw("sameP");
    myGraphBDTmcEE1->Draw("sameP");
    leg2->Draw();
    crocee1.SaveAs("roc_EE1_dataVsMc.png");  
    myH2->Draw();
    myGraphBDTdataEE1->Draw("sameP");
    myGraphBDTmcEE1->Draw("sameP");
    leg2->Draw();
    crocee1.SetLogx();
    crocee1.SaveAs("roc_EE1_dataVsMc_log.png");  

    fileROC.cd();
    myGraphBDTdataEE1->Write("myGraphBDTdataEE1");
    
    // ---------------------------------------------------
    FiguresOfMeritEvaluator myRocBDTdataEE2;
    myRocBDTdataEE2.addSignal("BDT", mvaSignalMcEE1);      
    myRocBDTdataEE2.addBackgrounds(mvaBkgDataEE2);
    myRocBDTdataEE2.setCutDirection(">");
    TGraph *myGraphBDTdataEE2 = myRocBDTdataEE2.getFOM("BDT",2);
    myGraphBDTdataEE2->SetTitle("BDT");
    myGraphBDTdataEE2->SetName("BDT");
    myGraphBDTdataEE2->SetMarkerColor(1);
    myGraphBDTdataEE2->SetMarkerStyle(20);
    myGraphBDTdataEE2->SetMarkerSize(1);
    //
    FiguresOfMeritEvaluator myRocBDTmcEE2;
    myRocBDTmcEE2.addSignal("BDT", mvaSignalMcEE2);      
    myRocBDTmcEE2.addBackgrounds(mvaBkgMcEE2);
    myRocBDTmcEE2.setCutDirection(">");
    TGraph *myGraphBDTmcEE2 = myRocBDTmcEE2.getFOM("BDT",2);
    myGraphBDTmcEE2->SetTitle("BDT");
    myGraphBDTmcEE2->SetName("BDT");
    myGraphBDTmcEE2->SetMarkerColor(4);
    myGraphBDTmcEE2->SetMarkerStyle(20);
    myGraphBDTmcEE2->SetMarkerSize(1);
    //
    TCanvas crocee2("roc","",1);
    crocee2.SetGrid();
    myH->Draw();
    myGraphBDTdataEE2->Draw("sameP");
    myGraphBDTmcEE2->Draw("sameP");
    leg2->Draw();
    crocee2.SaveAs("roc_EE2_dataVsMc.png");  
    myH2->Draw();
    myGraphBDTdataEE2->Draw("sameP");
    myGraphBDTmcEE2->Draw("sameP");
    leg2->Draw();
    crocee2.SetLogx();
    crocee2.SaveAs("roc_EE2_dataVsMc_log.png");  

    fileROC.cd();
    myGraphBDTdataEE2->Write("myGraphBDTdataEE2");
  }
}
