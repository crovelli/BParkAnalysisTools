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

void mcVsFakes()
{
  // Files: 
  TFile fileFakesData("filesNew/myFileFakesInData__weightsBasedOnMc.root");     // from fakes selection applied to data (prepareInputsFromFakes)
  TFile fileFakesMc("filesNew/myFileFakesInMc__weightsBasedOnMc.root");         // from fakes selection applied to (B->KJPsi->Kmm) MC (prepareInputsFromFakes)

  // ----------------------------------------------------

  // Data histos, fake selection
  TH1F *mvaFakesData = (TH1F*)fileFakesData.Get("mvaFakes");
  TH1F *ptFakesData  = (TH1F*)fileFakesData.Get("ptFakes");
  TH1F *etaFakesData = (TH1F*)fileFakesData.Get("etaFakes");
  mvaFakesData->Sumw2();   
  ptFakesData->Sumw2();   
  etaFakesData->Sumw2();   

  // MC histos, fake selection
  TH1F *mvaFakesMc = (TH1F*)fileFakesMc.Get("mvaFakes");
  TH1F *ptFakesMc  = (TH1F*)fileFakesMc.Get("ptFakes");
  TH1F *etaFakesMc = (TH1F*)fileFakesMc.Get("etaFakes");
  mvaFakesMc->Sumw2();   
  ptFakesMc->Sumw2();   
  etaFakesMc->Sumw2();   

  // ----------------------------------------------------

  // Small eta/pt bins: data histos, fake selection 
  TH1F *mvaFakesData_EB0 = (TH1F*)fileFakesData.Get("mvaFakeEB0");
  TH1F *mvaFakesData_EB1 = (TH1F*)fileFakesData.Get("mvaFakeEB1");
  TH1F *mvaFakesData_EB2 = (TH1F*)fileFakesData.Get("mvaFakeEB2");
  TH1F *mvaFakesData_EB3 = (TH1F*)fileFakesData.Get("mvaFakeEB3");
  TH1F *mvaFakesData_EB4 = (TH1F*)fileFakesData.Get("mvaFakeEB4");
  TH1F *mvaFakesData_EE0 = (TH1F*)fileFakesData.Get("mvaFakeEE0");
  TH1F *mvaFakesData_EE1 = (TH1F*)fileFakesData.Get("mvaFakeEE1");
  TH1F *mvaFakesData_EE2 = (TH1F*)fileFakesData.Get("mvaFakeEE2");
  TH1F *mvaFakesData_EE3 = (TH1F*)fileFakesData.Get("mvaFakeEE3");
  TH1F *mvaFakesData_EE4 = (TH1F*)fileFakesData.Get("mvaFakeEE4");

  // Small eta/pt bins: MC histos, fake selection 
  TH1F *mvaFakesMc_EB0 = (TH1F*)fileFakesMc.Get("mvaFakeEB0");
  TH1F *mvaFakesMc_EB1 = (TH1F*)fileFakesMc.Get("mvaFakeEB1");
  TH1F *mvaFakesMc_EB2 = (TH1F*)fileFakesMc.Get("mvaFakeEB2");
  TH1F *mvaFakesMc_EB3 = (TH1F*)fileFakesMc.Get("mvaFakeEB3");
  TH1F *mvaFakesMc_EB4 = (TH1F*)fileFakesMc.Get("mvaFakeEB4");
  TH1F *mvaFakesMc_EE0 = (TH1F*)fileFakesMc.Get("mvaFakeEE0");
  TH1F *mvaFakesMc_EE1 = (TH1F*)fileFakesMc.Get("mvaFakeEE1");
  TH1F *mvaFakesMc_EE2 = (TH1F*)fileFakesMc.Get("mvaFakeEE2");
  TH1F *mvaFakesMc_EE3 = (TH1F*)fileFakesMc.Get("mvaFakeEE3");
  TH1F *mvaFakesMc_EE4 = (TH1F*)fileFakesMc.Get("mvaFakeEE4");

  // ----------------------------------------------------

  // Cosmetics
  mvaFakesData -> SetLineWidth(2);  
  mvaFakesData -> SetLineColor(3);  
  ptFakesData  -> SetLineWidth(2);  
  ptFakesData  -> SetLineColor(3);  
  etaFakesData -> SetLineWidth(2);  
  etaFakesData -> SetLineColor(3);  
  //
  mvaFakesMc -> SetLineWidth(2);  
  mvaFakesMc -> SetLineColor(6);  
  ptFakesMc  -> SetLineWidth(2);  
  ptFakesMc  -> SetLineColor(6);  
  etaFakesMc -> SetLineWidth(2);  
  etaFakesMc -> SetLineColor(6);  
  //
  mvaFakesData_EB0 -> SetLineWidth(2);
  mvaFakesData_EB0 -> SetLineColor(3);
  mvaFakesData_EB1 -> SetLineWidth(2);
  mvaFakesData_EB1 -> SetLineColor(3);
  mvaFakesData_EB2 -> SetLineWidth(2);
  mvaFakesData_EB2 -> SetLineColor(3);
  mvaFakesData_EB3 -> SetLineWidth(2);
  mvaFakesData_EB3 -> SetLineColor(3);
  mvaFakesData_EB4 -> SetLineWidth(2);
  mvaFakesData_EB4 -> SetLineColor(3);
  //
  mvaFakesData_EE0 -> SetLineWidth(2);
  mvaFakesData_EE0 -> SetLineColor(3);
  mvaFakesData_EE1 -> SetLineWidth(2);
  mvaFakesData_EE1 -> SetLineColor(3);
  mvaFakesData_EE2 -> SetLineWidth(2);
  mvaFakesData_EE2 -> SetLineColor(3);
  mvaFakesData_EE3 -> SetLineWidth(2);
  mvaFakesData_EE3 -> SetLineColor(3);
  mvaFakesData_EE4 -> SetLineWidth(2);
  mvaFakesData_EE4 -> SetLineColor(3);
  //
  mvaFakesMc_EB0 -> SetLineWidth(2);
  mvaFakesMc_EB0 -> SetLineColor(6);
  mvaFakesMc_EB1 -> SetLineWidth(2);
  mvaFakesMc_EB1 -> SetLineColor(6);
  mvaFakesMc_EB2 -> SetLineWidth(2);
  mvaFakesMc_EB2 -> SetLineColor(6);
  mvaFakesMc_EB3 -> SetLineWidth(2);
  mvaFakesMc_EB3 -> SetLineColor(6);
  mvaFakesMc_EB4 -> SetLineWidth(2);
  mvaFakesMc_EB4 -> SetLineColor(6);
  //
  mvaFakesMc_EE0 -> SetLineWidth(2);
  mvaFakesMc_EE0 -> SetLineColor(6);
  mvaFakesMc_EE1 -> SetLineWidth(2);
  mvaFakesMc_EE1 -> SetLineColor(6);
  mvaFakesMc_EE2 -> SetLineWidth(2);
  mvaFakesMc_EE2 -> SetLineColor(6);
  mvaFakesMc_EE3 -> SetLineWidth(2);
  mvaFakesMc_EE3 -> SetLineColor(6);
  mvaFakesMc_EE4 -> SetLineWidth(2);
  mvaFakesMc_EE4 -> SetLineColor(6);

  // Rebinning
  mvaFakesData_EB0 -> Rebin();
  mvaFakesData_EB1 -> Rebin();
  mvaFakesData_EB2 -> Rebin();
  mvaFakesData_EB3 -> Rebin();
  mvaFakesData_EB4 -> Rebin();
  mvaFakesData_EE0 -> Rebin();
  mvaFakesData_EE1 -> Rebin();
  mvaFakesData_EE2 -> Rebin();
  mvaFakesData_EE3 -> Rebin();
  mvaFakesData_EE4 -> Rebin();
  mvaFakesMc_EB0 -> Rebin();
  mvaFakesMc_EB1 -> Rebin();
  mvaFakesMc_EB2 -> Rebin();
  mvaFakesMc_EB3 -> Rebin();
  mvaFakesMc_EB4 -> Rebin();
  mvaFakesMc_EE0 -> Rebin();
  mvaFakesMc_EE1 -> Rebin();
  mvaFakesMc_EE2 -> Rebin();
  mvaFakesMc_EE3 -> Rebin();
  mvaFakesMc_EE4 -> Rebin();


  // --------------------------------------------------  
  // Plots: with fake selection: data vs Mc 
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);

  TLegend *legA;
  legA = new TLegend(0.55,0.65,0.85,0.90);
  legA->SetFillStyle(0);
  legA->SetBorderSize(0);
  legA->SetTextSize(0.05);
  legA->SetFillColor(0);
  legA->AddEntry(mvaFakesData, "Fake sel., data", "lp");
  legA->AddEntry(mvaFakesMc, "Fake sel., MC", "lp");
  //
  TLegend *legA2;
  legA2 = new TLegend(0.30,0.10,0.70,0.40);
  legA2->SetFillStyle(0);
  legA2->SetBorderSize(0);
  legA2->SetTextSize(0.05);
  legA2->SetFillColor(0);
  legA2->AddEntry(mvaFakesData, "Fake sel., data", "lp");
  legA2->AddEntry(mvaFakesMc, "Fake sel., MC", "lp");

  TCanvas cmvaA("cmvaA","cmvaA",1);
  mvaFakesData->SetTitle("");
  mvaFakesMc->SetTitle("");
  mvaFakesData->GetXaxis()->SetTitle("Id BDT");
  mvaFakesMc->GetXaxis()->SetTitle("Id BDT");
  mvaFakesMc->DrawNormalized("hist");
  mvaFakesData->DrawNormalized("samehist");
  legA->Draw();
  cmvaA.SaveAs("outputBDT_dataVsMc.png");

  TCanvas cptA("cptA","cptA",1);
  ptFakesData->SetTitle("");
  ptFakesMc->SetTitle("");
  ptFakesData->GetXaxis()->SetTitle("p_{T} [GeV]");
  ptFakesMc->GetXaxis()->SetTitle("p_{T} [GeV]");
  ptFakesMc->DrawNormalized("hist");
  ptFakesData->DrawNormalized("samehist");
  legA->Draw();
  cptA.SaveAs("pT_dataVsMc.png.png");

  TCanvas cetaA("cetaA","cetaA",1);
  etaFakesData->SetTitle("");
  etaFakesMc->SetTitle("");
  etaFakesData->GetXaxis()->SetTitle("#eta");
  etaFakesMc->GetXaxis()->SetTitle("#eta");
  etaFakesData->DrawNormalized("hist");
  etaFakesMc->DrawNormalized("samehist");
  legA2->Draw();
  cetaA.SaveAs("Eta_dataVsMc.png");

  TCanvas cmvaeb0("cmvaeb0","cmvaeb0",1);
  mvaFakesData_EB0->SetTitle("EB, 0.5 < pT < 1.0 GeV");
  mvaFakesMc_EB0->SetTitle("EB, 0.5 < pT < 1.0 GeV");
  mvaFakesData_EB0->GetXaxis()->SetTitle("Id BDT");
  mvaFakesMc_EB0->GetXaxis()->SetTitle("Id BDT");
  mvaFakesMc_EB0->DrawNormalized("hist");
  mvaFakesData_EB0->DrawNormalized("samehist");
  legA->Draw();
  cmvaeb0.SaveAs("outputBDT_dataVsMc_EB0.png");

  TCanvas cmvaeb1("cmvaeb1","cmvaeb1",1);
  mvaFakesData_EB1->SetTitle("EB, 1.0 < pT < 1.5 GeV");
  mvaFakesMc_EB1->SetTitle("EB, 1.0 < pT < 1.5 GeV");
  mvaFakesData_EB1->GetXaxis()->SetTitle("Id BDT");
  mvaFakesMc_EB1->GetXaxis()->SetTitle("Id BDT");
  mvaFakesMc_EB1->DrawNormalized("hist");
  mvaFakesData_EB1->DrawNormalized("samehist");
  legA->Draw();
  cmvaeb1.SaveAs("outputBDT_dataVsMc_EB1.png");

  TCanvas cmvaeb2("cmvaeb2","cmvaeb2",1);
  mvaFakesData_EB2->SetTitle("EB, 1.5 < pT < 2.0 GeV");
  mvaFakesMc_EB2->SetTitle("EB, 1.5 < pT < 2.0 GeV");
  mvaFakesData_EB2->GetXaxis()->SetTitle("Id BDT");
  mvaFakesMc_EB2->GetXaxis()->SetTitle("Id BDT");
  mvaFakesMc_EB2->DrawNormalized("hist");
  mvaFakesData_EB2->DrawNormalized("samehist");
  legA->Draw();
  cmvaeb2.SaveAs("outputBDT_dataVsMc_EB2.png");

  TCanvas cmvaeb3("cmvaeb3","cmvaeb3",1);
  mvaFakesData_EB3->SetTitle("EB, 2.0 < pT < 5.0 GeV");
  mvaFakesMc_EB3->SetTitle("EB, 2.0 < pT < 5.0 GeV");
  mvaFakesData_EB3->GetXaxis()->SetTitle("Id BDT");
  mvaFakesMc_EB3->GetXaxis()->SetTitle("Id BDT");
  mvaFakesMc_EB3->DrawNormalized("hist");
  mvaFakesData_EB3->DrawNormalized("samehist");
  legA->Draw();
  cmvaeb3.SaveAs("outputBDT_dataVsMc_EB3.png");

  TCanvas cmvaeb4("cmvaeb4","cmvaeb4",1);
  mvaFakesData_EB4->SetTitle("EB, pT > 5 GeV");
  mvaFakesMc_EB4->SetTitle("EB, pT > 5 GeV");
  mvaFakesData_EB4->GetXaxis()->SetTitle("Id BDT");
  mvaFakesMc_EB4->GetXaxis()->SetTitle("Id BDT");
  mvaFakesMc_EB4->DrawNormalized("hist");
  mvaFakesData_EB4->DrawNormalized("samehist");
  legA->Draw();
  cmvaeb4.SaveAs("outputBDT_dataVsMc_EB4.png");

  TCanvas cmvaee0("cmvaee0","cmvaee0",1);
  mvaFakesData_EE0->SetTitle("EE, 0.5 < pT < 1.0 GeV");
  mvaFakesMc_EE0->SetTitle("EE, 0.5 < pT < 1.0 GeV");
  mvaFakesData_EE0->GetXaxis()->SetTitle("Id BDT");
  mvaFakesMc_EE0->GetXaxis()->SetTitle("Id BDT");
  mvaFakesMc_EE0->DrawNormalized("hist");
  mvaFakesData_EE0->DrawNormalized("samehist");
  legA->Draw();
  cmvaee0.SaveAs("outputBDT_dataVsMc_EE0.png");

  TCanvas cmvaee1("cmvaee1","cmvaee1",1);
  mvaFakesData_EE1->SetTitle("EE, 1.0 < pT < 1.5 GeV");
  mvaFakesMc_EE1->SetTitle("EE, 1.0 < pT < 1.5 GeV");
  mvaFakesData_EE1->GetXaxis()->SetTitle("Id BDT");
  mvaFakesMc_EE1->GetXaxis()->SetTitle("Id BDT");
  mvaFakesMc_EE1->DrawNormalized("hist");
  mvaFakesData_EE1->DrawNormalized("samehist");
  legA->Draw();
  cmvaee1.SaveAs("outputBDT_dataVsMc_EE1.png");

  TCanvas cmvaee2("cmvaee2","cmvaee2",1);
  mvaFakesData_EE2->SetTitle("EE, 1.5 < pT < 2.0 GeV");
  mvaFakesMc_EE2->SetTitle("EE, 1.5 < pT < 2.0 GeV");
  mvaFakesData_EE2->GetXaxis()->SetTitle("Id BDT");
  mvaFakesMc_EE2->GetXaxis()->SetTitle("Id BDT");
  mvaFakesMc_EE2->DrawNormalized("hist");
  mvaFakesData_EE2->DrawNormalized("samehist");
  legA->Draw();
  cmvaee2.SaveAs("outputBDT_dataVsMc_EE2.png");

  TCanvas cmvaee3("cmvaee3","cmvaee3",1);
  mvaFakesData_EE3->SetTitle("EE, 2.0 < pT < 5.0 GeV");
  mvaFakesMc_EE3->SetTitle("EE, 2.0 < pT < 5.0 GeV");
  mvaFakesData_EE3->GetXaxis()->SetTitle("Id BDT");
  mvaFakesMc_EE3->GetXaxis()->SetTitle("Id BDT");
  mvaFakesMc_EE3->DrawNormalized("hist");
  mvaFakesData_EE3->DrawNormalized("samehist");
  legA->Draw();
  cmvaee3.SaveAs("outputBDT_dataVsMc_EE3.png");

  TCanvas cmvaee4("cmvaee4","cmvaee4",1);
  mvaFakesData_EE4->SetTitle("EE, pT > 5 GeV");
  mvaFakesMc_EE4->SetTitle("EE, pT > 5 GeV");
  mvaFakesData_EE4->GetXaxis()->SetTitle("Id BDT");
  mvaFakesMc_EE4->GetXaxis()->SetTitle("Id BDT");
  mvaFakesMc_EE4->DrawNormalized("hist");
  mvaFakesData_EE4->DrawNormalized("samehist");
  legA->Draw();
  cmvaee4.SaveAs("outputBDT_dataVsMc_EE4.png");
}
