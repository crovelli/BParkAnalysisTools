#include "TH1F.h"
#include "TFile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"

// For PU reweighting using number of vertices
// Take in input the outputs of nvtxInput.C
// Give the output to be used with TaPJpsiSelectionNaod (with proper options)

void nvtxWeights() {

  TFile fileData("/eos/cms/store/user/crovelli/LowPtEle/Vertices/March21/fileNumberVtx2018ALL.root");    
  // TFile fileMC("/eos/cms/store/user/crovelli/LowPtEle/Vertices/March21/fileNumberVtxMC.root");
  TFile fileMC("/eos/cms/store/user/crovelli/LowPtEle/Vertices/March21/fileNumberVtxWithWeights.root");
  
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

  TLegend* leg = new TLegend(0.15,0.65,0.50,0.90);
  leg->SetFillStyle(0); 
  leg->SetBorderSize(0); 
  leg->SetTextSize(0.03); 
  leg->SetFillColor(0);
  leg->AddEntry(dataNvtx,"DATA","l");
  leg->AddEntry(mcNvtx,"MC","l");

  // Plots
  TCanvas c1("c1","",1);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);
  dataNvtx->Draw("histE");
  mcNvtx->Draw("samehist");
  leg->Draw();
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
