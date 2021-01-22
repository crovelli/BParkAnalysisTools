#include "TH1F.h"
#include "TFile.h"
#include "TStyle.h"
#include "TCanvas.h"

// For PU reweighting using number of vertices
// Take in input the outputs of nvtxInput.C
// Give the output to be used with TaPJpsiSelectionNaod (with proper options)

void nvtxWeights() {

  TFile fileData("fileNumberVtx_2018D1_upto37.root");    
  // TFile fileMC("fileNumberMC.root");
  TFile fileMC("fileNumberMC_withWeights.root");
  
  TH1F *mcNvtx = (TH1F*)fileMC.Get("H_nvtx");
  mcNvtx->Sumw2();
  mcNvtx->Scale(1./mcNvtx->Integral());
  mcNvtx->SetTitle("mcNvtx");
  mcNvtx->SetName("mcNvtx");

  TH1F *dataNvtx = (TH1F*)fileData.Get("H_nvtx");
  dataNvtx->Sumw2();
  dataNvtx->Scale(1./dataNvtx->Integral());
  dataNvtx->SetTitle("dataNvtx");
  dataNvtx->SetName("dataNvtx");

  // Cosmetics
  mcNvtx->SetTitle("");
  mcNvtx->GetXaxis()->SetTitle("Number of vertices");
  mcNvtx->SetLineColor(2);
  mcNvtx->SetLineWidth(2);
  dataNvtx->SetTitle("");
  dataNvtx->GetXaxis()->SetTitle("Number of vertices");
  dataNvtx->SetLineColor(4);
  dataNvtx->SetLineWidth(2);
  
  // Plots
  TCanvas c1("c1","",1);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);
  dataNvtx->Draw("hist");
  mcNvtx->Draw("samehist");
  c1.SaveAs("vertices.png");

  TH1F *myClone = (TH1F*)dataNvtx->Clone("myClone");
  myClone->Divide(mcNvtx);
  myClone->SetTitle("weights");
  myClone->SetName("weights");
  myClone->Scale(1./myClone->Integral());
  // for (int ii=0; ii<53; ii++) myClone->SetBinError(ii,0);
  
  TFile *fileOut = new TFile("nvtxWeights.root","RECREATE");
  fileOut->cd();
  dataNvtx->Write();
  mcNvtx->Write();
  myClone->Write();
  fileOut->Close();
  
}
