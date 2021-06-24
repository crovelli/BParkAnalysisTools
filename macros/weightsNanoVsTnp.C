#define weightsNanoVsTnp_cxx

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

// To compute weights nani->tnp for signal distributions
// and to make control distributions
// before and after applying weights

void weightsNanoVsTnp()
{
  // Files
  TFile fileTnP("files_marchNoReg/probeLowPt/myFileMcAfterTnp_noKineWeights.root");   // using TnP selection applied to MC (prepareInputsFromMcWithTnP)
  //TFile fileTnP("files_marchNoReg/probeLowPt/myFileMcAfterTnp_withKineWeights.root");   // ""
  TFile fileNoTnP("files_marchNoReg/probeLowPt/myFileFromNani.root");                   // without TnP selection, directy applied to nanpAODs (prepareInputsFromNaniInMc)

  // Histos: pT
  TH1F *ptSignalMcFromTnP   = (TH1F*)fileTnP.Get("probePtSignalMc");
  TH1F *ptSignalMcFromTnPWW = (TH1F*)fileTnP.Get("probePtSignalMcWW");
  TH1F *ptSignalMcNoTnP     = (TH1F*)fileNoTnP.Get("ptSignalMc");
  ptSignalMcFromTnP   -> Scale(1./ptSignalMcFromTnP->Integral());   
  ptSignalMcFromTnPWW -> Scale(1./ptSignalMcFromTnPWW->Integral());   
  ptSignalMcNoTnP     -> Scale(1./ptSignalMcNoTnP->Integral());   

  // Histos: eta
  TH1F *etaSignalMcFromTnP   = (TH1F*)fileTnP.Get("probeEtaSignalMc");
  TH1F *etaSignalMcFromTnPWW = (TH1F*)fileTnP.Get("probeEtaSignalMcWW");
  TH1F *etaSignalMcNoTnP     = (TH1F*)fileNoTnP.Get("etaSignalMc");
  etaSignalMcFromTnP   -> Scale(1./etaSignalMcFromTnP->Integral());   
  etaSignalMcFromTnPWW -> Scale(1./etaSignalMcFromTnPWW->Integral());   
  etaSignalMcNoTnP     -> Scale(1./etaSignalMcNoTnP->Integral());   
  
  // Histos: mva
  TH1F *mvaSignalMcFromTnP   = (TH1F*)fileTnP.Get("probeMvaSignalMc");
  TH1F *mvaSignalMcFromTnPWW = (TH1F*)fileTnP.Get("probeMvaSignalMcWW");
  TH1F *mvaSignalMcNoTnP     = (TH1F*)fileNoTnP.Get("mvaSignalMc");
  mvaSignalMcFromTnP   -> Scale(1./mvaSignalMcFromTnP->Integral());   
  mvaSignalMcFromTnPWW -> Scale(1./mvaSignalMcFromTnPWW->Integral());   
  mvaSignalMcNoTnP     -> Scale(1./mvaSignalMcNoTnP->Integral());   

  // dxy significance
  TH1F *dxysigSignalMcFromTnP   = (TH1F*)fileTnP.Get("probeDxysigSignalMc");
  TH1F *dxysigSignalMcFromTnPWW = (TH1F*)fileTnP.Get("probeDxysigSignalMcWW");
  TH1F *dxysigSignalMcNoTnP     = (TH1F*)fileNoTnP.Get("dxysigSignalMc");
  dxysigSignalMcFromTnP   -> Scale(1./dxysigSignalMcFromTnP->Integral());   
  dxysigSignalMcFromTnPWW -> Scale(1./dxysigSignalMcFromTnPWW->Integral());   
  dxysigSignalMcNoTnP     -> Scale(1./dxysigSignalMcNoTnP->Integral());   

  // dz trg
  TH1F *dztrgSignalMcFromTnP   = (TH1F*)fileTnP.Get("probeDztrgSignalMc");
  TH1F *dztrgSignalMcFromTnPWW = (TH1F*)fileTnP.Get("probeDztrgSignalMcWW");
  TH1F *dztrgSignalMcNoTnP     = (TH1F*)fileNoTnP.Get("dztrgSignalMc");
  dztrgSignalMcFromTnP   -> Scale(1./dztrgSignalMcFromTnP->Integral());   
  dztrgSignalMcFromTnPWW -> Scale(1./dztrgSignalMcFromTnPWW->Integral());   
  dztrgSignalMcNoTnP     -> Scale(1./dztrgSignalMcNoTnP->Integral());   

  // 2Dim histos: pT vs eta
  TH2F *ptVsEtaSignalMcFromTnP = (TH2F*)fileTnP.Get("probePtVsEtaSignalMc");
  TH2F *ptVsEtaSignalMcNoTnP   = (TH2F*)fileNoTnP.Get("ptVsEtaSignalMc");
  //
  TH2F *ptVsEtaSignalMcFromTnPWW = (TH2F*)fileTnP.Get("probePtVsEtaSignalMcWW");

  // --------------------------------------------------------
  // Compute weights
  TH2F *ptVsEtaSignalWeights = (TH2F*)ptVsEtaSignalMcNoTnP->Clone("ptVsEtaSignalWeights");
  ptVsEtaSignalWeights->Divide(ptVsEtaSignalMcFromTnP);
  ptVsEtaSignalWeights->SetTitle("ptVsEtaSignalWeights");
  ptVsEtaSignalWeights->SetName("ptVsEtaSignalWeights");
  cout << "signal pre: min = " << ptVsEtaSignalWeights->GetMinimum() << ", max = " << ptVsEtaSignalWeights->GetMaximum() << endl;
  for (int iBinEta=0; iBinEta<(ptVsEtaSignalWeights->GetNbinsX()); iBinEta++) {
    for (int iBinPt=0; iBinPt<(ptVsEtaSignalWeights->GetNbinsY()); iBinPt++) {
      int iBinEtaP = iBinEta+1;
      int iBinPtP  = iBinPt+1;
      if ((ptVsEtaSignalWeights->GetBinContent(iBinEtaP,iBinPtP))>2) ptVsEtaSignalWeights->SetBinContent(iBinEtaP,iBinPtP,2);
    }}
  cout << "signal post: min = " << ptVsEtaSignalWeights->GetMinimum() << ", max = " << ptVsEtaSignalWeights->GetMaximum() << endl;
  ptVsEtaSignalWeights->Scale(1./ptVsEtaSignalWeights->Integral()); 
  cout << "signal scaled: min = " << ptVsEtaSignalWeights->GetMinimum() << ", max = " << ptVsEtaSignalWeights->GetMaximum() << endl;

  // --------------------------------------------------------
  // Cosmetics
  ptSignalMcFromTnP   -> SetLineWidth(2);  
  ptSignalMcFromTnPWW -> SetLineWidth(2);  
  ptSignalMcNoTnP     -> SetLineWidth(2);  
  ptSignalMcFromTnP   -> SetLineColor(2);  
  ptSignalMcFromTnPWW -> SetLineColor(3);  
  ptSignalMcNoTnP     -> SetLineColor(4);  
  ptSignalMcFromTnP   -> SetTitle("");
  ptSignalMcFromTnPWW -> SetTitle("");
  ptSignalMcNoTnP     -> SetTitle("");
  ptSignalMcFromTnP   -> GetXaxis()->SetTitle("pT [GeV]");
  ptSignalMcFromTnPWW -> GetXaxis()->SetTitle("pT [GeV]");
  ptSignalMcNoTnP     -> GetXaxis()->SetTitle("pT [GeV]");
  //
  etaSignalMcFromTnP   -> SetLineWidth(2);  
  etaSignalMcFromTnPWW -> SetLineWidth(2);  
  etaSignalMcNoTnP     -> SetLineWidth(2);  
  etaSignalMcFromTnP   -> SetLineColor(2);  
  etaSignalMcFromTnPWW -> SetLineColor(3);  
  etaSignalMcNoTnP     -> SetLineColor(4);  
  etaSignalMcFromTnP   -> SetTitle("");
  etaSignalMcFromTnPWW -> SetTitle("");
  etaSignalMcNoTnP     -> SetTitle("");
  etaSignalMcFromTnP   -> GetXaxis()->SetTitle("#eta");
  etaSignalMcFromTnPWW -> GetXaxis()->SetTitle("#eta");
  etaSignalMcNoTnP     -> GetXaxis()->SetTitle("#eta");
  //
  mvaSignalMcFromTnP   -> SetLineWidth(2);  
  mvaSignalMcFromTnPWW -> SetLineWidth(2);  
  mvaSignalMcNoTnP     -> SetLineWidth(2);  
  mvaSignalMcFromTnP   -> SetLineColor(2);  
  mvaSignalMcFromTnPWW -> SetLineColor(3);  
  mvaSignalMcNoTnP     -> SetLineColor(4);  
  mvaSignalMcFromTnP   -> SetTitle("");
  mvaSignalMcFromTnPWW -> SetTitle("");
  mvaSignalMcNoTnP     -> SetTitle("");
  mvaSignalMcFromTnP   -> GetXaxis()->SetTitle("ID BDT output");
  mvaSignalMcFromTnPWW -> GetXaxis()->SetTitle("ID BDT output");
  mvaSignalMcNoTnP     -> GetXaxis()->SetTitle("ID BDT output");
  //
  dxysigSignalMcFromTnP   -> SetLineWidth(2);  
  dxysigSignalMcFromTnPWW -> SetLineWidth(2);  
  dxysigSignalMcNoTnP     -> SetLineWidth(2);  
  dxysigSignalMcFromTnP   -> SetLineColor(2);  
  dxysigSignalMcFromTnPWW -> SetLineColor(3);  
  dxysigSignalMcNoTnP     -> SetLineColor(4);  
  dxysigSignalMcFromTnP   -> SetTitle("");
  dxysigSignalMcFromTnPWW -> SetTitle("");
  dxysigSignalMcNoTnP     -> SetTitle("");
  dxysigSignalMcFromTnP   -> GetXaxis()->SetTitle("dXY signif");
  dxysigSignalMcFromTnPWW -> GetXaxis()->SetTitle("dXY signif");
  dxysigSignalMcNoTnP     -> GetXaxis()->SetTitle("dXY signif");
  //
  dztrgSignalMcFromTnP   -> SetLineWidth(2);  
  dztrgSignalMcFromTnPWW -> SetLineWidth(2);  
  dztrgSignalMcNoTnP     -> SetLineWidth(2);  
  dztrgSignalMcFromTnP   -> SetLineColor(2);  
  dztrgSignalMcFromTnPWW -> SetLineColor(3);  
  dztrgSignalMcNoTnP     -> SetLineColor(4);  
  dztrgSignalMcFromTnP   -> SetTitle("");
  dztrgSignalMcFromTnPWW -> SetTitle("");
  dztrgSignalMcNoTnP     -> SetTitle("");
  dztrgSignalMcFromTnP   -> GetXaxis()->SetTitle("dZ trg");
  dztrgSignalMcFromTnPWW -> GetXaxis()->SetTitle("dZ trg");
  dztrgSignalMcNoTnP     -> GetXaxis()->SetTitle("dZ trg");


  // -------------------------------------------------------------
  // Plots
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);

  TLegend *leg;
  leg = new TLegend(0.45,0.55,0.75,0.80);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.05);
  leg->SetFillColor(0);
  leg->AddEntry(ptSignalMcFromTnP, "With tnp", "lp");
  leg->AddEntry(ptSignalMcNoTnP,   "No tnp", "lp");
  //
  TLegend *legww;
  legww = new TLegend(0.45,0.55,0.75,0.80);
  legww->SetFillStyle(0);
  legww->SetBorderSize(0);
  legww->SetTextSize(0.05);
  legww->SetFillColor(0);
  legww->AddEntry(ptSignalMcFromTnPWW, "With tnp, weighted", "lp");
  legww->AddEntry(ptSignalMcNoTnP,   "No tnp", "lp");

  TLegend *leg2;
  leg2 = new TLegend(0.10,0.65,0.40,0.90);
  leg2->SetFillStyle(0);
  leg2->SetBorderSize(0);
  leg2->SetTextSize(0.05);
  leg2->SetFillColor(0);
  leg2->AddEntry(ptSignalMcFromTnP, "With tnp", "lp");
  leg2->AddEntry(ptSignalMcNoTnP,   "No tnp", "lp");
  //
  TLegend *leg2ww;
  leg2ww = new TLegend(0.10,0.65,0.40,0.90);
  leg2ww->SetFillStyle(0);
  leg2ww->SetBorderSize(0);
  leg2ww->SetTextSize(0.05);
  leg2ww->SetFillColor(0);
  leg2ww->AddEntry(ptSignalMcFromTnPWW, "With tnp, weighted", "lp");
  leg2ww->AddEntry(ptSignalMcNoTnP,   "No tnp", "lp");

  TLegend *leg3;
  leg3 = new TLegend(0.35,0.10,0.65,0.35);
  leg3->SetFillStyle(0);
  leg3->SetBorderSize(0);
  leg3->SetTextSize(0.05);
  leg3->SetFillColor(0);
  leg3->AddEntry(ptSignalMcFromTnP, "With tnp", "lp");
  leg3->AddEntry(ptSignalMcNoTnP,   "No tnp", "lp");
  //
  TLegend *leg3ww;
  leg3ww = new TLegend(0.35,0.10,0.65,0.35);
  leg3ww->SetFillStyle(0);
  leg3ww->SetBorderSize(0);
  leg3ww->SetTextSize(0.05);
  leg3ww->SetFillColor(0);
  leg3ww->AddEntry(ptSignalMcFromTnPWW, "With tnp, weighted", "lp");
  leg3ww->AddEntry(ptSignalMcNoTnP,   "No tnp", "lp");

  TCanvas cs1("cs1","cs1",1);
  ptSignalMcNoTnP->DrawNormalized("hist");
  ptSignalMcFromTnP->DrawNormalized("samehist");
  leg->Draw();
  cs1.SaveAs("ptSignal_naniVsTnp.png");

  TCanvas cs2("cs2","cs2",1);
  ptSignalMcNoTnP->DrawNormalized("hist");
  ptSignalMcFromTnPWW->DrawNormalized("samehist");
  legww->Draw();
  cs2.SaveAs("ptSignal_naniVsTnp_withWeight.png");
  
  // -----------

  TCanvas cs3("cs3","cs3",1);
  etaSignalMcFromTnP->DrawNormalized("hist");
  etaSignalMcNoTnP->DrawNormalized("samehist");
  leg3->Draw();
  cs3.SaveAs("etaSignal_naniVsTnp.png");

  TCanvas cs4("cs4","cs4",1);
  etaSignalMcFromTnPWW->DrawNormalized("hist");
  etaSignalMcNoTnP->DrawNormalized("samehist");
  leg3ww->Draw();
  cs4.SaveAs("etaSignal_naniVsTnp_withWeight.png");

  TCanvas cs5("cs5","cs5",1);
  mvaSignalMcNoTnP->DrawNormalized("hist");
  mvaSignalMcFromTnP->DrawNormalized("samehist");
  leg2->Draw();
  cs5.SaveAs("mvaSignal_naniVsTnp.png");

  TCanvas cs6("cs6","cs6",1);
  mvaSignalMcNoTnP->DrawNormalized("hist");
  mvaSignalMcFromTnPWW->DrawNormalized("samehist");
  leg2ww->Draw();
  cs6.SaveAs("mvaSignal_naniVsTnp_withWeight.png");

  TCanvas cs7("cs7","cs7",1);
  dxysigSignalMcNoTnP->DrawNormalized("hist");
  dxysigSignalMcFromTnP->DrawNormalized("samehist");
  leg2->Draw();
  cs7.SaveAs("dxysigSignal_naniVsTnp.png");

  TCanvas cs8("cs8","cs8",1);
  dxysigSignalMcNoTnP->DrawNormalized("hist");
  dxysigSignalMcFromTnPWW->DrawNormalized("samehist");
  leg2ww->Draw();
  cs8.SaveAs("dxysigSignal_naniVsTnp_withWeight.png");

  TCanvas cs9("cs9","cs9",1);
  dztrgSignalMcNoTnP->DrawNormalized("hist");
  dztrgSignalMcFromTnP->DrawNormalized("samehist");
  leg2->Draw();
  cs9.SaveAs("dztrgSignal_naniVsTnp.png");

  TCanvas cs10("cs10","cs10",1);
  dztrgSignalMcNoTnP->DrawNormalized("hist");
  dztrgSignalMcFromTnPWW->DrawNormalized("samehist");
  leg2ww->Draw();
  cs10.SaveAs("dztrgSignal_naniVsTnp_withWeight.png");

  TCanvas csw("csw","cs",1);
  ptVsEtaSignalWeights->Draw("colz");
  csw.SaveAs("ptVsEtaSignalWeights.png"); 

  // -----------------------------------------------
  // Check weights 2dim
  TH2F *check2dSignal = (TH2F*)ptVsEtaSignalMcFromTnPWW->Clone("check2dSignal");  
  check2dSignal->Divide(ptVsEtaSignalMcNoTnP);

  // ----------------------------------------------  
  TFile myFile("weightFile_tnpVsNani.root","RECREATE");
  ptVsEtaSignalWeights->Write();
  check2dSignal->Write();
  myFile.Close();
}
