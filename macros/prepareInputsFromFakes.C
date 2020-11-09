#define prepareInputsFromFakes_cxx
#include "prepareInputsFromFakes.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TLorentzVector.h"
#include <iostream>

// To be run on data selected with fakes selection.
// Produce basic distributions to be compared with MC

using namespace std;

void prepareInputsFromFakes::Loop(bool applyWeight)
{
  if (fChain == 0) return;
  
  // B
  TH1F *bMassH      = new TH1F("bMassH",      "bMassH",      50, 4.4, 6.);
  TH1F *bMassAfterH = new TH1F("bMassAfterH", "bMassAfterH", 50, 4.4, 6.);
  bMassH->Sumw2();
  bMassAfterH->Sumw2();
  
  // MuMu
  TH1F *mumuMassH      = new TH1F("mumuMassH",      "mumuMassH",      50, 0., 5.);
  TH1F *mumuMassAfterH = new TH1F("mumuMassAfterH", "mumuMassAfterH", 50, 0., 5.);
  mumuMassH->Sumw2();
  mumuMassAfterH->Sumw2();

  // Fakes selection
  TH1F *dREleK       = new TH1F("dREleK",       "dREleK",       50, 0., 6.);
  TH1F *dREleMu1     = new TH1F("dREleMu1",     "dREleMu1",     50, 0., 6.);
  TH1F *dREleMu2     = new TH1F("dREleMu2",     "dREleMu2",     50, 0., 6.);
  TH1F *dREleMuMu    = new TH1F("dREleMuMu",    "dREleMuMu",    50, 0., 6.);
  TH1F *dREleKzoom   = new TH1F("dREleKzoom",   "dREleKzoom",   50, 0., 0.5);
  TH1F *dREleMu1zoom = new TH1F("dREleMu1zoom", "dREleMu1zoom", 50, 0., 0.5);
  dREleK->Sumw2();
  dREleMu1->Sumw2();
  dREleMu2->Sumw2();
  dREleMuMu->Sumw2();
  dREleKzoom->Sumw2();
  dREleMu1zoom->Sumw2();

  // To compute weights
  TH2F *probePtVsEtaFakeMc   = new TH2F("probePtVsEtaFakeMc","probePtVsEtaFakeMc",     40, -2.4, 2.4, 60, 0., 15.);
  TH2F *probePtVsEtaFakeMcWW = new TH2F("probePtVsEtaFakeMcWW","probePtVsEtaFakeMcWW", 40, -2.4, 2.4, 60, 0., 15.);
  probePtVsEtaFakeMc->Sumw2();
  probePtVsEtaFakeMcWW->Sumw2();
  
  // Full pT range - no weight
  TH1F *mvaFakes = new TH1F("mvaFakes", "mvaFakes",  60, -10., 10.); 
  TH1F *ptFakes  = new TH1F("ptFakes",  "ptFakes",   60,  0.,  15.);
  TH1F *etaFakes = new TH1F("etaFakes", "etaFakes",  40, -2.4, 2.4);
  mvaFakes->Sumw2();  
  ptFakes->Sumw2();
  etaFakes->Sumw2();  

  // Full pT range - with weight
  TH1F *mvaFakesWW = new TH1F("mvaFakesWW", "mvaFakesWW",  60, -10., 10.); 
  TH1F *ptFakesWW  = new TH1F("ptFakesWW",  "ptFakesWW",   60,  0.,  15.);
  TH1F *etaFakesWW = new TH1F("etaFakesWW", "etaFakesWW",  40, -2.4, 2.4);
  mvaFakesWW->Sumw2();  
  ptFakesWW->Sumw2();
  etaFakesWW->Sumw2();  

  // Small eta/pt bins - no weight
  TH1F *mvaFakeEB0 = new TH1F("mvaFakeEB0", "mvaFakeEB0", 60, -10., 10.);  // pT: 0.5-1.0
  TH1F *mvaFakeEB1 = new TH1F("mvaFakeEB1", "mvaFakeEB1", 60, -10., 10.);  // 1.0-1.5
  TH1F *mvaFakeEB2 = new TH1F("mvaFakeEB2", "mvaFakeEB2", 60, -10., 10.);  // 1.5-2.0
  TH1F *mvaFakeEB3 = new TH1F("mvaFakeEB3", "mvaFakeEB3", 60, -10., 10.);  // 2.0-5.0
  TH1F *mvaFakeEB4 = new TH1F("mvaFakeEB4", "mvaFakeEB4", 60, -10., 10.);  // >5
  TH1F *mvaFakeEE0 = new TH1F("mvaFakeEE0", "mvaFakeEE0", 60, -10., 10.);  // pT: 0.5-1.0
  TH1F *mvaFakeEE1 = new TH1F("mvaFakeEE1", "mvaFakeEE1", 60, -10., 10.);  // 1.0-1.5
  TH1F *mvaFakeEE2 = new TH1F("mvaFakeEE2", "mvaFakeEE2", 60, -10., 10.);  // 1.5-2.0
  TH1F *mvaFakeEE3 = new TH1F("mvaFakeEE3", "mvaFakeEE3", 60, -10., 10.);  // 2.0-5.0
  TH1F *mvaFakeEE4 = new TH1F("mvaFakeEE4", "mvaFakeEE4", 60, -10., 10.);  // >5
  mvaFakeEB0->Sumw2();
  mvaFakeEB1->Sumw2();
  mvaFakeEB2->Sumw2();
  mvaFakeEB3->Sumw2();
  mvaFakeEB4->Sumw2();
  mvaFakeEE0->Sumw2();
  mvaFakeEE1->Sumw2();
  mvaFakeEE2->Sumw2();
  mvaFakeEE3->Sumw2();
  mvaFakeEE4->Sumw2();

  // --------------------------
  // weights (pt vs eta)
  float minBPt[100], maxBPt[100], minBEta[100], maxBEta[100];      
  int nBinsWPt  = -999;
  int nBinsWEta = -999;
  for (int iBin=0; iBin<100; iBin++) { 
    minBPt[iBin]=-999.; 
    maxBPt[iBin]=999.; 
    minBEta[iBin]=-999.; 
    maxBEta[iBin]=999.; 
  }

  float weightsFakes[100][100];
  for (int iBinEta=0; iBinEta<100; iBinEta++) { 
    for (int iBinPt=0; iBinPt<100; iBinPt++) { 
      weightsFakes[iBinEta][iBinPt]=1.;
    }
  }

  if (applyWeight) {
    TFile fileWeight("weightFile_fakeVsNani.root"); 
    TH2F *fakeWeights = (TH2F*)fileWeight.Get("ptVsEtaFakeWeights"); 
    nBinsWPt  = fakeWeights->GetNbinsY(); 
    nBinsWEta = fakeWeights->GetNbinsX(); 
    
    for (int iBinEta=0; iBinEta<=nBinsWEta; iBinEta++) {
      for (int iBinPt=0; iBinPt<=nBinsWPt; iBinPt++) {
	minBPt[iBinPt]   = fakeWeights->GetYaxis()->GetBinLowEdge(iBinPt);
	maxBPt[iBinPt]   = fakeWeights->GetYaxis()->GetBinUpEdge(iBinPt);
	minBEta[iBinEta] = fakeWeights->GetXaxis()->GetBinLowEdge(iBinEta);
	maxBEta[iBinEta] = fakeWeights->GetXaxis()->GetBinUpEdge(iBinEta);
	weightsFakes[iBinEta][iBinPt] = fakeWeights->GetBinContent(iBinEta,iBinPt);    // 0=underflow, N=overflow
      }
    }

    // over/under flow ranges
    minBPt[0]  = -999.;
    minBEta[0] = -999.;
    maxBPt[nBinsWPt] = 999.;
    maxBEta[nBinsWEta] = 999.;

  } // apply weight



  // Loop over entries
  Long64_t nentries = fChain->GetEntriesFast();  
  Long64_t nbytes = 0, nb = 0;

  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    // Trigger
    if (hlt_9==0) continue;

    // Acceptance
    if (fabs(eleEta)>2.4 || (elePt<0.5)) cout << "problems with acceptance" << endl;

    // B peak only
    bMassH -> Fill(bMass);
    if (bMass<5 || bMass>5.5) continue;

    // MuMu
    TLorentzVector mu1TLV(0,0,0,0);     
    TLorentzVector mu2TLV(0,0,0,0);     
    mu1TLV.SetPtEtaPhiM(mu1Pt, mu1Eta, mu1Phi, 0);
    mu2TLV.SetPtEtaPhiM(mu2Pt, mu2Eta, mu2Phi, 0);  
    TLorentzVector mumuTLV = mu1TLV + mu2TLV;
    float mumuMass = mumuTLV.M();

    // JPsi peak only
    mumuMassH -> Fill(mumuMass);
    if (mumuMass<2.7 || mumuMass>3.5) continue;

    // Selected events
    bMassAfterH -> Fill(bMass);
    mumuMassAfterH -> Fill(mumuMass);

    // Ele-particle match
    dREleK    -> Fill(dR_ele_k);    
    dREleMu1  -> Fill(dR_ele_mu1);    
    dREleMu2  -> Fill(dR_ele_mu2);
    dREleMuMu -> Fill(dR_ele_mu1mu2);
    dREleKzoom   -> Fill(dR_ele_k);    
    dREleMu1zoom -> Fill(dR_ele_mu1);    

    // Matching ele only
    if (dR_ele_mu1mu2>0.5) continue;

    // pT vs eta weight   
    float ptEtaFakeWeight=1.;
    if (applyWeight) {
      for (int iBinEta=0; iBinEta<=nBinsWEta; iBinEta++) {
	for (int iBinPt=0; iBinPt<=nBinsWPt; iBinPt++) {
	  bool thisBin = false;
	  if (elePt>=minBPt[iBinPt] && elePt<maxBPt[iBinPt] && eleEta>=minBEta[iBinEta] && eleEta<maxBEta[iBinEta]) thisBin = true;
	  if (thisBin) ptEtaFakeWeight = weightsFakes[iBinEta][iBinPt];
	}
      }
    }

    // Fakes - before weight
    mvaFakes -> Fill(eleMvaId);
    ptFakes  -> Fill(elePt);
    etaFakes -> Fill(eleEta);
    probePtVsEtaFakeMc -> Fill(eleEta, elePt); 

    // Fakes - with weight
    mvaFakesWW -> Fill(eleMvaId, ptEtaFakeWeight);
    ptFakesWW  -> Fill(elePt, ptEtaFakeWeight);
    etaFakesWW -> Fill(eleEta, ptEtaFakeWeight);
    probePtVsEtaFakeMcWW -> Fill(eleEta, elePt, ptEtaFakeWeight); 

    // Finer bins
    if (fabs(eleEta)<1.5) {
      if (elePt>0.5 && elePt<1.0) mvaFakeEB0->Fill(eleMvaId);
      if (elePt>1.0 && elePt<1.5) mvaFakeEB1->Fill(eleMvaId);
      if (elePt>1.5 && elePt<2.0) mvaFakeEB2->Fill(eleMvaId);
      if (elePt>2.0 && elePt<5.0) mvaFakeEB3->Fill(eleMvaId);
      if (elePt>5.0) mvaFakeEB4->Fill(eleMvaId);
    } else {
      if (elePt>0.5 && elePt<1.0) mvaFakeEE0->Fill(eleMvaId);
      if (elePt>1.0 && elePt<1.5) mvaFakeEE1->Fill(eleMvaId);
      if (elePt>1.5 && elePt<2.0) mvaFakeEE2->Fill(eleMvaId);
      if (elePt>2.0 && elePt<5.0) mvaFakeEE3->Fill(eleMvaId);
      if (elePt>5.0) mvaFakeEE4->Fill(eleMvaId);
    }

  } // loop over entries


  // Plots
  gStyle->SetOptStat(0);

  TCanvas c0("c0","",1);
  bMassH->SetLineColor(1);
  bMassH->SetLineWidth(2);
  bMassH->SetTitle("");
  bMassH->GetXaxis()->SetTitle("m_{B}");
  bMassH->DrawNormalized("hist");
  c0.SaveAs("bMass.png");

  TCanvas c0a("c0a","",1);
  mumuMassH->SetLineColor(1);
  mumuMassH->SetLineWidth(2);
  mumuMassH->SetTitle("");
  mumuMassH->GetXaxis()->SetTitle("m_{#mu#mu}");
  mumuMassH->DrawNormalized("hist");
  c0a.SaveAs("mumuMass.png");

  TCanvas c0b("c0b","",1);
  bMassAfterH->SetLineColor(1);
  bMassAfterH->SetLineWidth(2);
  bMassAfterH->SetTitle("");
  bMassAfterH->GetXaxis()->SetTitle("m_{B}");
  bMassAfterH->DrawNormalized("hist");
  c0b.SaveAs("bMassAfterCuts.png");

  TCanvas c0c("c0c","",1);
  mumuMassAfterH->SetLineColor(1);
  mumuMassAfterH->SetLineWidth(2);
  mumuMassAfterH->SetTitle("");
  mumuMassAfterH->GetXaxis()->SetTitle("m_{#mu#mu}");
  mumuMassAfterH->DrawNormalized("hist");
  c0c.SaveAs("mumuMassAfterCuts.png");

  TCanvas c1a("c1a","",1);
  dREleK->SetLineColor(1);
  dREleK->SetLineWidth(2);
  dREleK->SetTitle("");
  dREleK->GetXaxis()->SetTitle("#Delta R (ele, K)");
  dREleK->DrawNormalized("hist");
  c1a.SaveAs("deltaR_eleK.png");

  TCanvas c1aa("c1aa","",1);
  dREleKzoom->SetLineColor(1);
  dREleKzoom->SetLineWidth(2);
  dREleKzoom->SetTitle("");
  dREleKzoom->GetXaxis()->SetTitle("#Delta R (ele, K)");
  dREleKzoom->DrawNormalized("hist");
  c1aa.SaveAs("deltaR_eleK_zoom.png");

  TCanvas c1b("c1b","",1);
  dREleMu1->SetLineColor(1);
  dREleMu1->SetLineWidth(2);
  dREleMu1->SetTitle("");
  dREleMu1->GetXaxis()->SetTitle("#Delta R (ele, #mu1)");
  dREleMu1->DrawNormalized("hist");
  c1b.SaveAs("deltaR_eleMu1.png");

  TCanvas c1bb("c1bb","",1);
  dREleMu1zoom->SetLineColor(1);
  dREleMu1zoom->SetLineWidth(2);
  dREleMu1zoom->SetTitle("");
  dREleMu1zoom->GetXaxis()->SetTitle("#Delta R (ele, #mu1)");
  dREleMu1zoom->DrawNormalized("hist");
  c1bb.SaveAs("deltaR_eleMu1zoom.png");

  TCanvas c1c("c1c","",1);
  dREleMu2->SetLineColor(1);
  dREleMu2->SetLineWidth(2);
  dREleMu2->SetTitle("");
  dREleMu2->GetXaxis()->SetTitle("#Delta R (ele, #mu2)");
  dREleMu2->DrawNormalized("hist");
  c1c.SaveAs("deltaR_eleMu2.png");

  TCanvas c1d("c1d","",1);
  dREleMuMu->SetLineColor(1);
  dREleMuMu->SetLineWidth(2);
  dREleMuMu->SetTitle("");
  dREleMuMu->GetXaxis()->SetTitle("#Delta R (ele, #mu#mu)");
  dREleMuMu->DrawNormalized("hist");
  c1d.SaveAs("deltaR_eleMuMu.png");


  // Dave histos
  TFile filefakes("myFileFakes.root","RECREATE");
  // 
  dREleK    -> Write(); 
  dREleMu1  -> Write(); 
  dREleMu2  -> Write(); 
  dREleMuMu -> Write(); 
  dREleKzoom   -> Write(); 
  dREleMu1zoom -> Write(); 
  // 
  mvaFakes -> Write();
  ptFakes  -> Write();
  etaFakes -> Write();
  // 
  mvaFakesWW -> Write();
  ptFakesWW  -> Write();
  etaFakesWW -> Write();
  //
  probePtVsEtaFakeMc -> Write();
  probePtVsEtaFakeMcWW -> Write();
  //
  mvaFakeEB0 -> Write();
  mvaFakeEB1 -> Write();
  mvaFakeEB2 -> Write();
  mvaFakeEB3 -> Write();
  mvaFakeEB4 -> Write();
  mvaFakeEE0 -> Write();
  mvaFakeEE1 -> Write();
  mvaFakeEE2 -> Write();
  mvaFakeEE3 -> Write();
  mvaFakeEE4 -> Write();
}
