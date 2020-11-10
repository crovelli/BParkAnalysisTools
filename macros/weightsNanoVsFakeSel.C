#define weightsNanoVsFakeSel_cxx

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

// To compute weights nani->fake selection for fakes distributions
// and to make control distributions
// before and after applying weights

void weightsNanoVsFakeSel(bool isLowPt)
{
  // Files
  TFile fileFakes("myFileFakesInMc.root");           // using fake selection applied to MC (prepareInputsFromFakes)
  TFile fileNoFakes("myFileFromNani.root");          // without fake selection, directy applied to nanpAODs (prepareInputsFromNaniInMc)

  // Histos: pT
  TH1F *ptFromFakesSel     = (TH1F*)fileFakes.Get("ptFakes");
  TH1F *ptFromFakesSelWW   = (TH1F*)fileFakes.Get("ptFakesWW");
  TH1F *ptFakeMcNoFakesSel = (TH1F*)fileNoFakes.Get("ptFakeMc");
  ptFromFakesSel     -> Sumw2();   
  ptFromFakesSelWW   -> Sumw2();   
  ptFakeMcNoFakesSel -> Sumw2(); 
  ptFromFakesSel     -> Scale(1./ptFromFakesSel->Integral());   
  ptFromFakesSelWW   -> Scale(1./ptFromFakesSelWW->Integral());   
  ptFakeMcNoFakesSel -> Scale(1./ptFakeMcNoFakesSel->Integral());   
  //
  // Histos: eta
  TH1F *etaFromFakesSel     = (TH1F*)fileFakes.Get("etaFakes");
  TH1F *etaFromFakesSelWW   = (TH1F*)fileFakes.Get("etaFakesWW");
  TH1F *etaFakeMcNoFakesSel = (TH1F*)fileNoFakes.Get("etaFakeMc");
  etaFromFakesSel     -> Sumw2();   
  etaFromFakesSelWW   -> Sumw2();   
  etaFakeMcNoFakesSel -> Sumw2(); 
  etaFromFakesSel     -> Scale(1./etaFromFakesSel->Integral());   
  etaFromFakesSelWW   -> Scale(1./etaFromFakesSelWW->Integral());   
  etaFakeMcNoFakesSel -> Scale(1./etaFakeMcNoFakesSel->Integral());   
  //
  // Histos: mva
  TH1F *mvaFromFakesSel     = (TH1F*)fileFakes.Get("mvaFakes");
  TH1F *mvaFromFakesSelWW   = (TH1F*)fileFakes.Get("mvaFakesWW");
  TH1F *mvaFakeMcNoFakesSel = (TH1F*)fileNoFakes.Get("mvaFakeMc");
  mvaFromFakesSel     -> Sumw2();   
  mvaFromFakesSelWW   -> Sumw2();   
  mvaFakeMcNoFakesSel -> Sumw2(); 
  mvaFromFakesSel     -> Scale(1./mvaFromFakesSel->Integral());   
  mvaFromFakesSelWW   -> Scale(1./mvaFromFakesSelWW->Integral());   
  mvaFakeMcNoFakesSel -> Scale(1./mvaFakeMcNoFakesSel->Integral());   

  // 2Dim histos: pT vs eta
  TH2F *ptVsEtaFromFakesSel     = (TH2F*)fileFakes.Get("probePtVsEtaFakeMc");
  TH2F *ptVsEtaFromFakesSelWW   = (TH2F*)fileFakes.Get("probePtVsEtaFakeMcWW");
  TH2F *ptVsEtaFakeMcNoFakesSel = (TH2F*)fileNoFakes.Get("ptVsEtaFakeMc");

  // --------------------------------------------------------
  // Compute weights
  TH2F *ptVsEtaFakeWeights = (TH2F*)ptVsEtaFakeMcNoFakesSel->Clone("ptVsEtaFakesWeights");
  ptVsEtaFakeWeights->Sumw2(); 
  ptVsEtaFakeWeights->Divide(ptVsEtaFromFakesSel);
  ptVsEtaFakeWeights->SetTitle("ptVsEtaFakeWeights");
  ptVsEtaFakeWeights->SetName("ptVsEtaFakeWeights");
  cout << "signal pre: min = " << ptVsEtaFakeWeights->GetMinimum() << ", max = " << ptVsEtaFakeWeights->GetMaximum() << endl;
  for (int iBinEta=0; iBinEta<(ptVsEtaFakeWeights->GetNbinsX()); iBinEta++) {
    for (int iBinPt=0; iBinPt<(ptVsEtaFakeWeights->GetNbinsY()); iBinPt++) {
      int iBinEtaP = iBinEta+1;
      int iBinPtP  = iBinPt+1;
      // if ((ptVsEtaFakeWeights->GetBinContent(iBinEtaP,iBinPtP))>2) ptVsEtaFakeWeights->SetBinContent(iBinEtaP,iBinPtP,2);
    }}
  cout << "signal post: min = " << ptVsEtaFakeWeights->GetMinimum() << ", max = " << ptVsEtaFakeWeights->GetMaximum() << endl;
  ptVsEtaFakeWeights->Scale(1./ptVsEtaFakeWeights->Integral()); 
  cout << "signal scaled: min = " << ptVsEtaFakeWeights->GetMinimum() << ", max = " << ptVsEtaFakeWeights->GetMaximum() << endl;

  // --------------------------------------------------------
  // Cosmetics
  ptFromFakesSel     -> SetLineWidth(2);  
  ptFromFakesSelWW   -> SetLineWidth(2);  
  ptFakeMcNoFakesSel -> SetLineWidth(2);  
  ptFromFakesSel     -> SetLineColor(2);  
  ptFromFakesSelWW   -> SetLineColor(3);  
  ptFakeMcNoFakesSel -> SetLineColor(4);  
  ptFromFakesSel     -> SetTitle("");
  ptFromFakesSelWW   -> SetTitle("");
  ptFakeMcNoFakesSel -> SetTitle("");
  ptFromFakesSel     -> GetXaxis()->SetTitle("pT [GeV]");
  ptFromFakesSelWW   -> GetXaxis()->SetTitle("pT [GeV]");
  ptFakeMcNoFakesSel -> GetXaxis()->SetTitle("pT [GeV]");
  //
  etaFromFakesSel     -> SetLineWidth(2);  
  etaFromFakesSelWW   -> SetLineWidth(2);  
  etaFakeMcNoFakesSel -> SetLineWidth(2);  
  etaFromFakesSel     -> SetLineColor(2);  
  etaFromFakesSelWW   -> SetLineColor(3);  
  etaFakeMcNoFakesSel -> SetLineColor(4);  
  etaFromFakesSel     -> SetTitle("");
  etaFromFakesSelWW   -> SetTitle("");
  etaFakeMcNoFakesSel -> SetTitle("");
  etaFromFakesSel     -> GetXaxis()->SetTitle("#eta");
  etaFromFakesSelWW   -> GetXaxis()->SetTitle("#eta");
  etaFakeMcNoFakesSel -> GetXaxis()->SetTitle("#eta");
  //
  mvaFromFakesSel     -> SetLineWidth(2);  
  mvaFromFakesSelWW   -> SetLineWidth(2);  
  mvaFakeMcNoFakesSel -> SetLineWidth(2);  
  mvaFromFakesSel     -> SetLineColor(2);  
  mvaFromFakesSelWW   -> SetLineColor(3);  
  mvaFakeMcNoFakesSel -> SetLineColor(4);  
  mvaFromFakesSel     -> SetTitle("");
  mvaFromFakesSelWW   -> SetTitle("");
  mvaFakeMcNoFakesSel -> SetTitle("");
  mvaFromFakesSel     -> GetXaxis()->SetTitle("ID BDT output");
  mvaFromFakesSelWW   -> GetXaxis()->SetTitle("ID BDT output");
  mvaFakeMcNoFakesSel -> GetXaxis()->SetTitle("ID BDT output");

  // -------------------------------------------------------------
  // Plots
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  // 
  TLegend *leg;
  leg = new TLegend(0.45,0.55,0.75,0.80);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.05);
  leg->SetFillColor(0);
  leg->AddEntry(ptFromFakesSel,     "With fakes sel", "lp");
  leg->AddEntry(ptFakeMcNoFakesSel, "No fakes sel", "lp");
  //
  TLegend *legww;
  legww = new TLegend(0.45,0.55,0.75,0.80);
  legww->SetFillStyle(0);
  legww->SetBorderSize(0);
  legww->SetTextSize(0.05);
  legww->SetFillColor(0);
  legww->AddEntry(ptFromFakesSelWW,   "With fakes sel, weighted", "lp");
  legww->AddEntry(ptFakeMcNoFakesSel, "No tnp", "lp");

  TLegend *leg2;
  leg2 = new TLegend(0.10,0.65,0.40,0.90);
  leg2->SetFillStyle(0);
  leg2->SetBorderSize(0);
  leg2->SetTextSize(0.05);
  leg2->SetFillColor(0);
  leg2->AddEntry(ptFromFakesSel,     "With fakes sel", "lp");
  leg2->AddEntry(ptFakeMcNoFakesSel, "No fakes sel", "lp");
  //
  TLegend *leg2ww;
  leg2ww = new TLegend(0.10,0.65,0.40,0.90);
  leg2ww->SetFillStyle(0);
  leg2ww->SetBorderSize(0);
  leg2ww->SetTextSize(0.05);
  leg2ww->SetFillColor(0);
  leg2ww->AddEntry(ptFromFakesSelWW,   "With fakes sel, weighted", "lp");
  leg2ww->AddEntry(ptFakeMcNoFakesSel, "No fakes sel", "lp");

  TLegend *leg3;
  leg3 = new TLegend(0.35,0.10,0.65,0.35);
  leg3->SetFillStyle(0);
  leg3->SetBorderSize(0);
  leg3->SetTextSize(0.05);
  leg3->SetFillColor(0);
  leg3->AddEntry(ptFromFakesSel,     "With fakes sel", "lp");
  leg3->AddEntry(ptFakeMcNoFakesSel, "No fakes sel", "lp");
  //
  TLegend *leg3ww;
  leg3ww = new TLegend(0.35,0.10,0.65,0.35);
  leg3ww->SetFillStyle(0);
  leg3ww->SetBorderSize(0);
  leg3ww->SetTextSize(0.05);
  leg3ww->SetFillColor(0);
  leg3ww->AddEntry(ptFromFakesSelWW,   "With fakes sel, weighted", "lp");
  leg3ww->AddEntry(ptFakeMcNoFakesSel, "No fakes sel", "lp");

  TCanvas cs1("cs1","cs1",1);
  if (isLowPt==0) ptFakeMcNoFakesSel->Rebin(4);
  if (isLowPt==0) ptFromFakesSel->Rebin(4);
  ptFakeMcNoFakesSel->DrawNormalized("hist");
  ptFromFakesSel->DrawNormalized("samehist");
  leg->Draw();
  cs1.SaveAs("ptSignal_naniVsFakeSel.png");

  TCanvas cs2("cs2","cs2",1);
  if (isLowPt==0) ptFromFakesSelWW->Rebin(4);
  ptFakeMcNoFakesSel->DrawNormalized("hist");
  ptFromFakesSelWW->DrawNormalized("samehist");
  legww->Draw();
  cs2.SaveAs("ptSignal_naniVsFakeSel_withWeight.png");
  
  // -----------

  TCanvas cs3("cs3","cs3",1);
  if (isLowPt==0) etaFakeMcNoFakesSel->Rebin(4);
  if (isLowPt==0) etaFromFakesSel->Rebin(4);
  etaFromFakesSel->DrawNormalized("hist");
  etaFakeMcNoFakesSel->DrawNormalized("samehist");
  leg3->Draw();
  cs3.SaveAs("etaSignal_naniVsFakeSel.png");

  TCanvas cs4("cs4","cs4",1);
  if (isLowPt==0) etaFromFakesSelWW->Rebin(4);
  etaFromFakesSelWW->DrawNormalized("hist");
  etaFakeMcNoFakesSel->DrawNormalized("samehist");
  leg3ww->Draw();
  cs4.SaveAs("etaSignal_naniVsFakeSel_withWeight.png");

  // -----------

  TCanvas cs5("cs5","cs5",1);
  if (isLowPt==0) mvaFakeMcNoFakesSel->Rebin(4);
  if (isLowPt==0) mvaFromFakesSel->Rebin(4);
  mvaFakeMcNoFakesSel->DrawNormalized("hist");
  mvaFromFakesSel->DrawNormalized("samehist");
  leg2->Draw();
  cs5.SaveAs("mvaSignal_naniVsFakeSel.png");

  TCanvas cs6("cs6","cs6",1);
  if (isLowPt==0) mvaFromFakesSelWW->Rebin(4);
  mvaFakeMcNoFakesSel->DrawNormalized("hist");
  mvaFromFakesSelWW->DrawNormalized("samehist");
  leg2ww->Draw();
  cs6.SaveAs("mvaSignal_naniVsFakeSel_withWeight.png");

  // -----------

  TCanvas csw("csw","cs",1);
  ptVsEtaFakeWeights->Draw("colz");
  csw.SaveAs("ptVsEtaFakeWeights.png"); 

  // -----------

  // Check weights 2dim
  TH2F *check2dSignal = (TH2F*)ptVsEtaFromFakesSelWW->Clone("check2dSignal");  
  check2dSignal->Sumw2();
  check2dSignal->Divide(ptVsEtaFakeMcNoFakesSel);

  // ----------------------------------------------  
  TFile myFile("weightFile_fakeVsNani.root","RECREATE");
  ptVsEtaFakeWeights->Write();
  check2dSignal->Write();
  myFile.Close();
}
