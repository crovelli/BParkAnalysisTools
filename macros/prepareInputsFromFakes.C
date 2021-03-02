#define prepareInputsFromFakes_cxx
#include "prepareInputsFromFakes.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include "TLorentzVector.h"
#include <iostream>

// To be run on data selected with fakes selection.
// Produce basic distributions to be compared with MC

using namespace std;

void prepareInputsFromFakes::Loop(bool applyKineWeight, int isLowPt=-1, bool studyOverlap=0, bool applyWeightOverlap=0, bool applyPuWeight=1)
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
  TH1F *dREleMuMu      = new TH1F("dREleMuMu",      "dREleMuMu",       50, 0., 6.);
  TH1F *dREleMuMuAfter = new TH1F("dREleMuMuAfter", "dREleMuMuAfter",  60, 0., 0.6);
  TH1F *dREleK         = new TH1F("dREleK",         "dREleK",         200, 0., 2.);
  TH1F *dREleMu1       = new TH1F("dREleMu1",       "dREleMu1",       200, 0., 2.);
  TH1F *dREleMu2       = new TH1F("dREleMu2",       "dREleMu2",       200, 0., 2.);
  TH1F *dREleMin       = new TH1F("dREleMin",       "dREleMin",       200, 0., 2.);
  TH1F *dREleKzoom     = new TH1F("dREleKzoom",     "dREleKzoom",      50, 0., 0.5);
  TH1F *dREleMu1zoom   = new TH1F("dREleMu1zoom",   "dREleMu1zoom",    50, 0., 0.5);
  TH1F *dREleMu2zoom   = new TH1F("dREleMu2zoom",   "dREleMu2zoom",    50, 0., 0.5);
  TH1F *dREleMinzoom   = new TH1F("dREleMinzoom",   "dREleMinzoom",    50, 0., 0.5);
  dREleMuMu->Sumw2();
  dREleMuMuAfter->Sumw2();
  dREleK->Sumw2();
  dREleMu1->Sumw2();
  dREleMu2->Sumw2();
  dREleMin->Sumw2();
  dREleKzoom->Sumw2();
  dREleMu1zoom->Sumw2();
  dREleMu2zoom->Sumw2();
  dREleMinzoom->Sumw2();

  // Closest particle to
  TH1F *dRMinParticleCutGt002 = new TH1F("dRMinParticleCutGt002", "dRMinParticleCutGt002", 3, 0.5, 3.5);
  TH1F *dRMinParticleCutLt002 = new TH1F("dRMinParticleCutLt002", "dRMinParticleCutLt002", 3, 0.5, 3.5);
  dRMinParticleCutGt002->Sumw2();
  dRMinParticleCutLt002->Sumw2();

  // Isolation cone around muons, for cases where the muon is the closest particle and dR with mu > 0.2 
  TH1F *dRMinEleProbe_Iso0d3mu1 = new TH1F("dRMinEleProbe_Iso0d3mu1", "dRMinEleProbe_Iso0d3mu1", 50, 0., 0.5);
  TH1F *dRMinEleProbe_Iso0d3mu2 = new TH1F("dRMinEleProbe_Iso0d3mu2", "dRMinEleProbe_Iso0d3mu2", 50, 0., 0.5);
  dRMinEleProbe_Iso0d3mu1->Sumw2();
  dRMinEleProbe_Iso0d3mu2->Sumw2();
  
  // dR between ele and muon when the mu is the closest particle dut dR>0.02
  TH1F *dREleMu1_mu1Closest = new TH1F("dREleMu1_mu1Closest", "dREleMu1_mu1Closest", 50, 0., 0.5);
  TH1F *dREleMu2_mu2Closest = new TH1F("dREleMu2_mu2Closest", "dREleMu2_mu2Closest", 50, 0., 0.5);

  // 2dim plots
  TH2F *TwoDimDr_eleMu1_vs_eleProbeIso03 = new TH2F("TwoDimDr_eleMu1_vs_eleProbeIso03", "TwoDimDr_eleMu1_vs_eleProbeIso03", 50, 0., 0.5, 50, 0., 0.5);
  TH2F *TwoDimDr_eleMu2_vs_eleProbeIso03 = new TH2F("TwoDimDr_eleMu2_vs_eleProbeIso03", "TwoDimDr_eleMu2_vs_eleProbeIso03", 50, 0., 0.5, 50, 0., 0.5);
  
  // To compute weights
  TH2F *probePtVsEtaFakeMc   = new TH2F("probePtVsEtaFakeMc","probePtVsEtaFakeMc",     40, -2.4, 2.4, 60, 0., 15.);
  TH2F *probePtVsEtaFakeMcWW = new TH2F("probePtVsEtaFakeMcWW","probePtVsEtaFakeMcWW", 40, -2.4, 2.4, 60, 0., 15.);
  probePtVsEtaFakeMc->Sumw2();
  probePtVsEtaFakeMcWW->Sumw2();

  // To compute weights 
  TH2F *ptVsEtaFakes_LptPfOverlap   = new TH2F("ptVsEtaFakes_LptPfOverlap",  "ptVsEtaFakes_LptPfOverlap",   40, -2.4, 2.4, 60, 0., 15.);   
  TH2F *ptVsEtaFakes_LptPfOverlapWW = new TH2F("ptVsEtaFakes_LptPfOverlapWW","ptVsEtaFakes_LptPfOverlapWW", 40, -2.4, 2.4, 60, 0., 15.);   
  
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
  TH1F *mvaFakeEE0 = new TH1F("mvaFakeEE0", "mvaFakeEE0", 60, -10., 10.);  // pT: 0.5-1.0
  TH1F *mvaFakeEE1 = new TH1F("mvaFakeEE1", "mvaFakeEE1", 60, -10., 10.);  // 1.0-1.5
  TH1F *mvaFakeEE2 = new TH1F("mvaFakeEE2", "mvaFakeEE2", 60, -10., 10.);  // 1.5-2.0
  mvaFakeEB0->Sumw2();
  mvaFakeEB1->Sumw2();
  mvaFakeEB2->Sumw2();
  mvaFakeEB3->Sumw2();
  mvaFakeEE0->Sumw2();
  mvaFakeEE1->Sumw2();
  mvaFakeEE2->Sumw2();

  // Fake rate from muons and from Ks
  TH1F *mvaFakesMu1 = new TH1F("mvaFakesMu1", "mvaFakesMu1",  60, -10., 10.);   
  TH1F *mvaFakesMu2 = new TH1F("mvaFakesMu2", "mvaFakesMu2",  60, -10., 10.);   
  TH1F *mvaFakesK   = new TH1F("mvaFakesK",   "mvaFakesK",    60, -10., 10.);   
  TH1F *ptFakesMu1  = new TH1F("ptFakesMu1",  "ptFakesMu1",   60,  0.,  15.);
  TH1F *ptFakesMu2  = new TH1F("ptFakesMu2",  "ptFakesMu2",   60,  0.,  15.);
  TH1F *ptFakesK    = new TH1F("ptFakesK",    "ptFakesK",     60,  0.,  15.);
  TH1F *etaFakesMu1 = new TH1F("etaFakesMu1", "etaFakesMu1",  40, -2.4, 2.4);
  TH1F *etaFakesMu2 = new TH1F("etaFakesMu2", "etaFakesMu2",  40, -2.4, 2.4);
  TH1F *etaFakesK   = new TH1F("etaFakesK",   "etaFakesK",    40, -2.4, 2.4);
  mvaFakesMu1->Sumw2();
  mvaFakesMu2->Sumw2();
  mvaFakesK->Sumw2();
  ptFakesMu1->Sumw2();
  ptFakesMu2->Sumw2();
  ptFakesK->Sumw2();
  etaFakesMu1->Sumw2();
  etaFakesMu2->Sumw2();
  etaFakesK->Sumw2();

  // Full pT range, PF-LowPt overlap studies
  TH1F *mvaFakes_LptPfOverlap = new TH1F("mvaFakes_LptPfOverlap", "mvaFakes_LptPfOverlap",  60, -10., 10.); 
  TH1F *ptFakes_LptPfOverlap  = new TH1F("ptFakes_LptPfOverlap",  "ptFakes_LptPfOverlap",   60,  0.,  15.);
  TH1F *etaFakes_LptPfOverlap = new TH1F("etaFakes_LptPfOverlap", "etaFakes_LptPfOverlap",  40, -2.4, 2.4);
  mvaFakes_LptPfOverlap->Sumw2();  
  ptFakes_LptPfOverlap->Sumw2();
  etaFakes_LptPfOverlap->Sumw2();  
  //
  TH1F *mvaFakes_LptNotPfOverlap = new TH1F("mvaFakes_LptNotPfOverlap", "mvaFakes_LptNotPfOverlap",  60, -10., 10.); 
  TH1F *ptFakes_LptNotPfOverlap  = new TH1F("ptFakes_LptNotPfOverlap",  "ptFakes_LptNotPfOverlap",   60,  0.,  15.);
  TH1F *etaFakes_LptNotPfOverlap = new TH1F("etaFakes_LptNotPfOverlap", "etaFakes_LptNotPfOverlap",  40, -2.4, 2.4);
  mvaFakes_LptNotPfOverlap->Sumw2();  
  ptFakes_LptNotPfOverlap->Sumw2();
  etaFakes_LptNotPfOverlap->Sumw2();  
  //
  TH1F *mvaFakes_LptPfOverlapWW = new TH1F("mvaFakes_LptPfOverlapWW", "mvaFakes_LptPfOverlapWW",  60, -10., 10.); 
  TH1F *ptFakes_LptPfOverlapWW  = new TH1F("ptFakes_LptPfOverlapWW",  "ptFakes_LptPfOverlapWW",   60,  0.,  15.);
  TH1F *etaFakes_LptPfOverlapWW = new TH1F("etaFakes_LptPfOverlapWW", "etaFakes_LptPfOverlapWW",  40, -2.4, 2.4);
  mvaFakes_LptPfOverlapWW->Sumw2();  
  ptFakes_LptPfOverlapWW->Sumw2();
  etaFakes_LptPfOverlapWW->Sumw2();  


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

  float minBPtOv[100], maxBPtOv[100], minBEtaOv[100], maxBEtaOv[100];      
  int nBinsWPtOv  = -999;
  int nBinsWEtaOv = -999;
  for (int iBin=0; iBin<100; iBin++) { 
    minBPtOv[iBin]=-999.; 
    maxBPtOv[iBin]=999.; 
    minBEtaOv[iBin]=-999.; 
    maxBEtaOv[iBin]=999.; 
  }

  float weightsFakes[100][100];
  for (int iBinEta=0; iBinEta<100; iBinEta++) { 
    for (int iBinPt=0; iBinPt<100; iBinPt++) { 
      weightsFakes[iBinEta][iBinPt]=1.;
    }
  }

  float weightsFakesLptPfOverlap[100][100];
  for (int iBinEta=0; iBinEta<100; iBinEta++) { 
    for (int iBinPt=0; iBinPt<100; iBinPt++) { 
      weightsFakesLptPfOverlap[iBinEta][iBinPt]=1.;
    }
  }

  if (applyKineWeight) {
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
	weightsFakes[iBinEta][iBinPt] = fakeWeights->GetBinContent(iBinEta,iBinPt); 
      }
    }
    minBPt[0]  = -999.;
    minBEta[0] = -999.;
    maxBPt[nBinsWPt] = 999.;
    maxBEta[nBinsWEta] = 999.;
  }   // apply weight
  
  if (applyWeightOverlap==1 && studyOverlap==1) {
    TFile fileWeightOv("weightFile_withPFoverlap_datiVsMc.root"); 
    TH2F *fakeWeightsOv = (TH2F*)fileWeightOv.Get("ptVsEtaWeights"); 
    nBinsWPtOv  = fakeWeightsOv->GetNbinsY(); 
    nBinsWEtaOv = fakeWeightsOv->GetNbinsX(); 
    for (int iBinEta=0; iBinEta<=nBinsWEtaOv; iBinEta++) {
      for (int iBinPt=0; iBinPt<=nBinsWPtOv; iBinPt++) {
	minBPtOv[iBinPt]   = fakeWeightsOv->GetYaxis()->GetBinLowEdge(iBinPt);
	maxBPtOv[iBinPt]   = fakeWeightsOv->GetYaxis()->GetBinUpEdge(iBinPt);
	minBEtaOv[iBinEta] = fakeWeightsOv->GetXaxis()->GetBinLowEdge(iBinEta);
	maxBEtaOv[iBinEta] = fakeWeightsOv->GetXaxis()->GetBinUpEdge(iBinEta);
	weightsFakesLptPfOverlap[iBinEta][iBinPt] = fakeWeightsOv->GetBinContent(iBinEta,iBinPt); 
      }
    }
    minBPtOv[0]  = -999.;
    minBEtaOv[0] = -999.;
    maxBPtOv[nBinsWPtOv] = 999.;
    maxBEtaOv[nBinsWEtaOv] = 999.;
  } // apply weight



  // Loop over entries
  Long64_t nentries = fChain->GetEntriesFast();  
  Long64_t nbytes = 0, nb = 0;

  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    float theWeight = 1.;
    if (applyPuWeight) theWeight = weight;

    // Trigger - already applied to data only when formatting
    if (hlt_9==0) continue;

    // Acceptance - should already be applied
    if (fabs(eleEta)>2.4 || (elePt<0.5)) cout << "problems with acceptance" << endl;
    
    // Should already be applied
    if (isLowPt==1 && eleIsLowPt==0) continue; 
    if (isLowPt==0 && eleIsPF==0) continue;   
    
    // Sanity check
    if (studyOverlap==1 && isLowPt!=-1) { cout << "should not happen!" << endl; continue; }
    if (studyOverlap==0 && applyWeightOverlap==1) { cout << "should not happen!" << endl; continue; }
    
    // B peak only
    bMassH -> Fill(bMass, theWeight);
    if (bMass<5 || bMass>5.5) continue;

    // Ele
    TLorentzVector eleTLV(0,0,0,0); 
    eleTLV.SetPtEtaPhiM(elePt, eleEta, elePhi, 0);   

    // MuMu
    TLorentzVector mu1TLV(0,0,0,0);     
    TLorentzVector mu2TLV(0,0,0,0);     
    mu1TLV.SetPtEtaPhiM(mu1Pt, mu1Eta, mu1Phi, 0);
    mu2TLV.SetPtEtaPhiM(mu2Pt, mu2Eta, mu2Phi, 0);  
    TLorentzVector mumuTLV = mu1TLV + mu2TLV;
    float mumuMass = mumuTLV.M();

    // JPsi peak only
    mumuMassH -> Fill(mumuMass, theWeight);
    if (mumuMass<2.7 || mumuMass>3.5) continue;

    // Selected events
    bMassAfterH -> Fill(bMass, theWeight);
    mumuMassAfterH -> Fill(mumuMass, theWeight);

    // Ele-mm match 
    dREleMuMu -> Fill(dR_ele_mu1mu2, theWeight);

    // Matching ele only
    if (dR_ele_mu1mu2>0.5) continue;

    // Ele-particle match
    dREleMuMuAfter -> Fill(dR_ele_mu1mu2, theWeight);
    if (dR_ele_k>2 || dR_ele_mu1>2 || dR_ele_mu2>2) cout << "deltaR>range" << endl;
    dREleK    -> Fill(dR_ele_k, theWeight);    
    dREleMu1  -> Fill(dR_ele_mu1, theWeight);    
    dREleMu2  -> Fill(dR_ele_mu2, theWeight);
    dREleKzoom   -> Fill(dR_ele_k, theWeight);    
    dREleMu1zoom -> Fill(dR_ele_mu1, theWeight);    
    dREleMu2zoom -> Fill(dR_ele_mu2, theWeight);    

    // closest particle
    float dR_min = 999.;
    int p_drmin  = -999;
    if (dR_ele_k<dR_ele_mu1 && dR_ele_k<dR_ele_mu2)   { dR_min = dR_ele_k;   p_drmin = 1; }
    if (dR_ele_mu1<dR_ele_k && dR_ele_mu1<dR_ele_mu2) { dR_min = dR_ele_mu1; p_drmin = 2; }
    if (dR_ele_mu2<dR_ele_k && dR_ele_mu2<dR_ele_mu1) { dR_min = dR_ele_mu2; p_drmin = 3; }

    dREleMin     -> Fill(dR_min, theWeight);
    dREleMinzoom -> Fill(dR_min, theWeight);
    if (dR_min>0.02) dRMinParticleCutGt002 -> Fill(p_drmin, theWeight);
    if (dR_min<0.02) dRMinParticleCutLt002 -> Fill(p_drmin, theWeight);

    // if closest particle is mu1, but dR>0.02
    if(p_drmin==2 && dR_min>0.02) { // probe closest to electron in the isolation cone around mu1 
      dREleMu1_mu1Closest -> Fill(dR_ele_mu1, theWeight);
      if (probeCloseToEleForMu1Pt>=0) {
	TLorentzVector probeTLV(0,0,0,0);     
	probeTLV.SetPtEtaPhiM(probeCloseToEleForMu1Pt, probeCloseToEleForMu1Eta, probeCloseToEleForMu1Phi, 0);
	float dRmu1probe = probeTLV.DeltaR(mu1TLV);
	float dReleprobe = probeTLV.DeltaR(eleTLV);
	if (dRmu1probe<0.3) { 
	  dRMinEleProbe_Iso0d3mu1->Fill(dReleprobe, theWeight);
	  TwoDimDr_eleMu1_vs_eleProbeIso03 -> Fill(dReleprobe,dR_ele_mu1, theWeight);
	}
	if (dRmu1probe>0.4) cout << "problem with isolation cone around mu1" << endl;
      }
    }

    // if closest particle is mu2, but dR>0.02 
    if(p_drmin==3 && dR_min>0.02) { // probe closest to electron in the isolation cone around mu2
      dREleMu2_mu2Closest -> Fill(dR_ele_mu2, theWeight);
      if (probeCloseToEleForMu2Pt>=0) {
	TLorentzVector probeTLV(0,0,0,0);     
	probeTLV.SetPtEtaPhiM(probeCloseToEleForMu2Pt, probeCloseToEleForMu2Eta, probeCloseToEleForMu2Phi, 0);
	float dRmu2probe = probeTLV.DeltaR(mu2TLV);
	float dReleprobe = probeTLV.DeltaR(eleTLV);
	if (dRmu2probe<0.3) { 
	  dRMinEleProbe_Iso0d3mu2->Fill(dReleprobe, theWeight);
	  TwoDimDr_eleMu2_vs_eleProbeIso03 -> Fill(dReleprobe,dR_ele_mu2, theWeight);
	}
	if (dRmu2probe>0.4) cout << "problem with isolation cone around mu2" << endl;
      }
    }

    // pT vs eta weight   
    float ptEtaFakeWeight=1.;
    if (applyKineWeight) {
      for (int iBinEta=0; iBinEta<=nBinsWEta; iBinEta++) {
	for (int iBinPt=0; iBinPt<=nBinsWPt; iBinPt++) {
	  bool thisBin = false;
	  if (elePt>=minBPt[iBinPt] && elePt<maxBPt[iBinPt] && eleEta>=minBEta[iBinEta] && eleEta<maxBEta[iBinEta]) thisBin = true;
	  if (thisBin) ptEtaFakeWeight = weightsFakes[iBinEta][iBinPt];
	}
      }
    }

    // pT vs eta weight for lowPt electrons with PF overlap
    float ptEtaFakeWeightWithOverlap=1.;    
    if (applyWeightOverlap==1 && studyOverlap==1) {
      for (int iBinEta=0; iBinEta<=nBinsWEtaOv; iBinEta++) {
	for (int iBinPt=0; iBinPt<=nBinsWPtOv; iBinPt++) {
	  bool thisBin = false;
	  if (elePt>=minBPtOv[iBinPt] && elePt<maxBPtOv[iBinPt] && eleEta>=minBEtaOv[iBinEta] && eleEta<maxBEtaOv[iBinEta]) thisBin = true;
	  if (thisBin) ptEtaFakeWeightWithOverlap = weightsFakesLptPfOverlap[iBinEta][iBinPt];
	}
      }
    }
  
    // which id? Low Pt or PF
    float eleId = -999.;
    if (isLowPt==1) eleId = eleMvaId;
    if (isLowPt==0) eleId = elePfmvaId;

    // Fakes - before weight
    mvaFakes -> Fill(eleId, theWeight);
    ptFakes  -> Fill(elePt, theWeight);
    etaFakes -> Fill(eleEta, theWeight);
    probePtVsEtaFakeMc -> Fill(eleEta, elePt, theWeight); 

    // Fakes before weight: split FR from K wrt muons
    if (dR_min<0.02) {   // electron matching k or mu
      if(p_drmin==2) {                 // leading mu
	mvaFakesMu1 -> Fill(eleId, theWeight);
	ptFakesMu1  -> Fill(elePt, theWeight);
	etaFakesMu1 -> Fill(eleEta, theWeight);
      }
      if(p_drmin==3) {                 // subleading mu
	mvaFakesMu2 -> Fill(eleId, theWeight);
	ptFakesMu2  -> Fill(elePt, theWeight);
	etaFakesMu2 -> Fill(eleEta, theWeight);
      }
      if(p_drmin==1) {                 // K
	mvaFakesK -> Fill(eleId, theWeight);
	ptFakesK  -> Fill(elePt, theWeight);
	etaFakesK -> Fill(eleEta, theWeight);
      }
    }
      
    // Fakes - with weight
    mvaFakesWW -> Fill(eleId, ptEtaFakeWeight*theWeight);
    ptFakesWW  -> Fill(elePt, ptEtaFakeWeight*theWeight);
    etaFakesWW -> Fill(eleEta, ptEtaFakeWeight*theWeight);
    probePtVsEtaFakeMcWW -> Fill(eleEta, elePt, ptEtaFakeWeight*theWeight); 

    // Finer bins
    if (fabs(eleEta)<1.5) {
      if (elePt>0.5 && elePt<1.5) mvaFakeEB0->Fill(eleId, theWeight);
      if (elePt>1.5 && elePt<2.0) mvaFakeEB1->Fill(eleId, theWeight);
      if (elePt>2.0 && elePt<5.0) mvaFakeEB2->Fill(eleId, theWeight);
      if (elePt>5.0) mvaFakeEB3->Fill(eleId, theWeight);
    } else {
      if (elePt>0.5 && elePt<2.0) mvaFakeEE0->Fill(eleId, theWeight);
      if (elePt>2.0 && elePt<5.0) mvaFakeEE1->Fill(eleId, theWeight);
      if (elePt>5.0) mvaFakeEE2->Fill(eleId, theWeight);
    }


    // Low pt with/wo overlap with PF studies - if required
    if (studyOverlap==1) {
      if (eleIsLowPt==1 && eleIsPFOverlap==1) {
	mvaFakes_LptPfOverlap -> Fill(eleMvaId, theWeight);
	ptFakes_LptPfOverlap  -> Fill(elePt, theWeight);
	etaFakes_LptPfOverlap -> Fill(eleEta, theWeight);
	mvaFakes_LptPfOverlapWW -> Fill(eleMvaId, ptEtaFakeWeightWithOverlap*theWeight);
	ptFakes_LptPfOverlapWW  -> Fill(elePt, ptEtaFakeWeightWithOverlap*theWeight);
	etaFakes_LptPfOverlapWW -> Fill(eleEta, ptEtaFakeWeightWithOverlap*theWeight);
	ptVsEtaFakes_LptPfOverlap   -> Fill(eleEta, elePt, theWeight);
	ptVsEtaFakes_LptPfOverlapWW -> Fill(eleEta, elePt, ptEtaFakeWeightWithOverlap*theWeight);
      }
      if (eleIsLowPt==1 && eleIsPFOverlap==0) {
	mvaFakes_LptNotPfOverlap -> Fill(eleMvaId, theWeight);
	ptFakes_LptNotPfOverlap  -> Fill(elePt, theWeight);
	etaFakes_LptNotPfOverlap -> Fill(eleEta, theWeight);
      }
    }

  } // loop over entries


  // Summary
  cout << endl;
  cout << "1and2 bin : k = " << dREleKzoom->GetBinContent(1)+dREleKzoom->GetBinContent(2) 
       << ", mu1 = " << dREleMu1zoom->GetBinContent(1)+dREleMu1zoom->GetBinContent(2) 
       << ", mu2 = " << dREleMu2zoom->GetBinContent(1)+dREleMu2zoom->GetBinContent(2) << endl;
  cout << "tot : k = "   << dREleKzoom->GetEntries()     << ", mu1 = " << dREleMu1zoom->GetEntries()     << ", mu2 = " << dREleMu2zoom->GetEntries() << endl;
  cout << "Fraction : k = " << (dREleKzoom->GetBinContent(1)+dREleKzoom->GetBinContent(2))/dREleKzoom->GetEntries() 
       << ", mu1 = " << (dREleMu1zoom->GetBinContent(1)+dREleMu1zoom->GetBinContent(2))/ dREleMu1zoom->GetEntries()
       << ", mu2 = " << (dREleMu2zoom->GetBinContent(1)+dREleMu2zoom->GetBinContent(2))/ dREleMu2zoom->GetEntries() << endl;
  cout << "min : 1and2 bin = " << dREleMinzoom->GetBinContent(1)+dREleMinzoom->GetBinContent(2) << ", tot = " << dREleMinzoom->GetEntries() << endl;
  cout << endl;
  cout << "dRMinParticleCutLt002: 1bin = " << dRMinParticleCutLt002->GetBinContent(1) 
       << ", 2 = " << dRMinParticleCutLt002->GetBinContent(2) << ", 3 = " << dRMinParticleCutLt002->GetBinContent(3) 
       << ", tot = " << dRMinParticleCutLt002->GetEntries() << endl;
  cout << endl;
  cout << "Isolation, mu1, 03: 1and2bin (all) = " << dRMinEleProbe_Iso0d3mu1->GetBinContent(1)+dRMinEleProbe_Iso0d3mu1->GetBinContent(2) 
       << ", tot = " << dRMinEleProbe_Iso0d3mu1->GetEntries() 
       << ", fraction = " << (dRMinEleProbe_Iso0d3mu1->GetBinContent(1)+dRMinEleProbe_Iso0d3mu1->GetBinContent(2))/dRMinEleProbe_Iso0d3mu1->GetEntries() << endl;
  cout << "Isolation, mu2, 03: 1and2bin (all) = " << dRMinEleProbe_Iso0d3mu2->GetBinContent(1)+dRMinEleProbe_Iso0d3mu2->GetBinContent(2) 
       << ", tot = " << dRMinEleProbe_Iso0d3mu2->GetEntries() 
       << ", fraction = " << (dRMinEleProbe_Iso0d3mu2->GetBinContent(1)+dRMinEleProbe_Iso0d3mu2->GetBinContent(2))/dRMinEleProbe_Iso0d3mu2->GetEntries() << endl;

  // ---------------------------------------------------------
  // Save histos
  TFile filefakes("myFileFakes.root","RECREATE");
  // 
  dREleK        -> Write(); 
  dREleMu1      -> Write(); 
  dREleMu2      -> Write(); 
  dREleMuMu     -> Write(); 
  dREleKzoom    -> Write(); 
  dREleMu1zoom  -> Write(); 
  dREleMin      -> Write();
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
  mvaFakeEE0 -> Write();
  mvaFakeEE1 -> Write();
  mvaFakeEE2 -> Write();
  //
  mvaFakes_LptNotPfOverlap -> Write();
  ptFakes_LptNotPfOverlap -> Write();
  etaFakes_LptNotPfOverlap -> Write();
  mvaFakes_LptPfOverlap -> Write();
  ptFakes_LptPfOverlap -> Write();
  etaFakes_LptPfOverlap -> Write();
  //
  mvaFakes_LptPfOverlapWW -> Write();
  ptFakes_LptPfOverlapWW -> Write();
  etaFakes_LptPfOverlapWW -> Write();
  //
  ptVsEtaFakes_LptPfOverlap->Write();
  ptVsEtaFakes_LptPfOverlapWW->Write();
  // ---------------------------------------------------------


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
  dREleK->Scale(1./dREleK->GetEntries());
  dREleK->Draw("hist");
  c1a.SaveAs("deltaR_eleK.png");

  TCanvas c1aa("c1aa","",1);
  dREleKzoom->SetLineColor(1);
  dREleKzoom->SetLineWidth(2);
  dREleKzoom->SetTitle("");
  dREleKzoom->GetXaxis()->SetTitle("#Delta R (ele, K)");
  dREleKzoom->Scale(1./dREleKzoom->GetEntries());
  dREleKzoom->Draw("hist");
  c1aa.SaveAs("deltaR_eleK_zoom.png");

  TCanvas c1b("c1b","",1);
  dREleMu1->SetLineColor(1);
  dREleMu1->SetLineWidth(2);
  dREleMu1->SetTitle("");
  dREleMu1->GetXaxis()->SetTitle("#Delta R (ele, #mu1)");
  dREleMu1->Scale(1./dREleMu1->GetEntries());
  dREleMu1->Draw("hist");
  c1b.SaveAs("deltaR_eleMu1.png");

  TCanvas c1bb("c1bb","",1);
  dREleMu1zoom->SetLineColor(1);
  dREleMu1zoom->SetLineWidth(2);
  dREleMu1zoom->SetTitle("");
  dREleMu1zoom->GetXaxis()->SetTitle("#Delta R (ele, #mu1)");
  dREleMu1zoom->Scale(1./dREleMu1zoom->GetEntries());
  dREleMu1zoom->Draw("hist");
  c1bb.SaveAs("deltaR_eleMu1zoom.png");

  TCanvas c1c("c1c","",1);
  dREleMu2->SetLineColor(1);
  dREleMu2->SetLineWidth(2);
  dREleMu2->SetTitle("");
  dREleMu2->GetXaxis()->SetTitle("#Delta R (ele, #mu2)");
  dREleMu2->Scale(1./dREleMu2->GetEntries());
  dREleMu2->Draw("hist");
  c1c.SaveAs("deltaR_eleMu2.png");

  TCanvas c1cc("c1cc","",1);
  dREleMu2zoom->SetLineColor(1);
  dREleMu2zoom->SetLineWidth(2);
  dREleMu2zoom->SetTitle("");
  dREleMu2zoom->GetXaxis()->SetTitle("#Delta R (ele, #mu2)");
  dREleMu2zoom->Scale(1./dREleMu2zoom->GetEntries());
  dREleMu2zoom->Draw("hist");
  c1cc.SaveAs("deltaR_eleMu2zoom.png");

  TCanvas c1d("c1d","",1);
  dREleMuMu->SetLineColor(1);
  dREleMuMu->SetLineWidth(2);
  dREleMuMu->SetTitle("");
  dREleMuMu->GetXaxis()->SetTitle("#Delta R (ele, #mu#mu)");
  dREleMuMu->DrawNormalized("hist");
  c1d.SaveAs("deltaR_eleMuMu.png");

  TCanvas c1dd("c1dd","",1);
  dREleMuMuAfter->SetLineColor(1);
  dREleMuMuAfter->SetLineWidth(2);
  dREleMuMuAfter->SetTitle("");
  dREleMuMuAfter->GetXaxis()->SetTitle("#Delta R (ele, #mu#mu)");
  dREleMuMuAfter->DrawNormalized("hist");
  c1dd.SaveAs("deltaR_eleMuMu_after.png");

  TCanvas c1e("c1e","",1);
  dREleMin->SetLineColor(1);
  dREleMin->SetLineWidth(2);
  dREleMin->SetTitle("");
  dREleMin->GetXaxis()->SetTitle("min #Delta R");
  dREleMin->Scale(1./dREleMin->GetEntries());
  dREleMin->Draw("hist");
  c1e.SaveAs("deltaR_eleMin.png");

  TCanvas c1ee("c1ee","",1);
  dREleMinzoom->SetLineColor(1);
  dREleMinzoom->SetLineWidth(2);
  dREleMinzoom->SetTitle("");
  dREleMinzoom->GetXaxis()->SetTitle("min #Delta R");
  dREleMinzoom->Scale(1./dREleMinzoom->GetEntries());
  dREleMinzoom->Draw("hist");
  c1ee.SaveAs("deltaR_eleMinzoom.png");

  TCanvas c1f("c1f","",1);
  dRMinParticleCutGt002->SetLineColor(1);
  dRMinParticleCutGt002->SetLineWidth(2);
  dRMinParticleCutGt002->SetTitle("");
  dRMinParticleCutGt002->GetXaxis()->SetTitle("min #Delta R");
  dRMinParticleCutGt002->DrawNormalized("hist");
  c1f.SaveAs("minDrParticle_DrMinGt0d02.png");

  TCanvas c1ff("c1ff","",1);
  dRMinParticleCutLt002->SetLineColor(1);
  dRMinParticleCutLt002->SetLineWidth(2);
  dRMinParticleCutLt002->SetTitle("");
  dRMinParticleCutLt002->GetXaxis()->SetTitle("min #Delta R");
  dRMinParticleCutLt002->DrawNormalized("hist");
  c1ff.SaveAs("minDrParticle_DrMinLt0d02.png");

  TCanvas c1f1("c1f1","",1);
  dRMinEleProbe_Iso0d3mu1->SetLineColor(1);
  dRMinEleProbe_Iso0d3mu1->SetLineWidth(2);
  dRMinEleProbe_Iso0d3mu1->SetTitle("");
  dRMinEleProbe_Iso0d3mu1->GetXaxis()->SetTitle("#DeltaR(ele, probe), #mu1 isol");
  dRMinEleProbe_Iso0d3mu1->Scale(1./dRMinEleProbe_Iso0d3mu1->GetEntries());
  dRMinEleProbe_Iso0d3mu1->Draw("hist");
  c1f1.SaveAs("dRMinEleProbe_Iso0d3mu1.png");

  TCanvas c1f2("c1f2","",1);
  TwoDimDr_eleMu1_vs_eleProbeIso03 ->SetTitle("");
  TwoDimDr_eleMu1_vs_eleProbeIso03 ->GetXaxis()->SetTitle("#DeltaR(ele, probe), #mu1 isol");
  TwoDimDr_eleMu1_vs_eleProbeIso03 ->GetYaxis()->SetTitle("#Delta R (ele, #mu1)");   
  TwoDimDr_eleMu1_vs_eleProbeIso03 ->Draw("colz");
  c1f2.SaveAs("TwoDimDr_eleMu1_vs_eleProbeIso03.png");

  TCanvas c1f3("c1f3","",1);
  dREleMu1_mu1Closest -> SetLineColor(1);
  dREleMu1_mu1Closest -> SetLineWidth(2);
  dREleMu1_mu1Closest ->SetTitle("");
  dREleMu1_mu1Closest ->GetXaxis()->SetTitle("#Delta R (ele, #mu1)");
  dREleMu1_mu1Closest ->DrawNormalized("hist");
  c1f3.SaveAs("dREleMu1_mu1IsClosest.png");

  TCanvas c1ff1("c1ff1","",1);
  dRMinEleProbe_Iso0d3mu2->SetLineColor(1);
  dRMinEleProbe_Iso0d3mu2->SetLineWidth(2);
  dRMinEleProbe_Iso0d3mu2->SetTitle("");
  dRMinEleProbe_Iso0d3mu2->GetXaxis()->SetTitle("#DeltaR(ele, probe), #mu2 isol");
  dRMinEleProbe_Iso0d3mu2->Scale(1./dRMinEleProbe_Iso0d3mu2->GetEntries());
  dRMinEleProbe_Iso0d3mu2->Draw("hist");
  c1ff1.SaveAs("dRMinEleProbe_Iso0d3mu2.png");

  TCanvas c1ff2("c1ff2","",1);
  TwoDimDr_eleMu2_vs_eleProbeIso03 ->SetTitle("");
  TwoDimDr_eleMu2_vs_eleProbeIso03 ->GetXaxis()->SetTitle("#DeltaR(ele, probe), #mu2 isol");
  TwoDimDr_eleMu2_vs_eleProbeIso03 ->GetYaxis()->SetTitle("#Delta R (ele, #mu2)");   
  TwoDimDr_eleMu2_vs_eleProbeIso03 ->Draw("colz");
  c1ff2.SaveAs("TwoDimDr_eleMu2_vs_eleProbeIso03.png");

  TCanvas c1ff3("c1ff3","",1);
  dREleMu2_mu2Closest -> SetLineColor(1);
  dREleMu2_mu2Closest -> SetLineWidth(2);
  dREleMu2_mu2Closest ->SetTitle("");
  dREleMu2_mu2Closest ->GetXaxis()->SetTitle("#Delta R (ele, #mu2)");
  dREleMu2_mu2Closest ->DrawNormalized("hist");
  c1ff3.SaveAs("dREleMu2_mu2IsClosest.png");

  // ----------------------------------
  // Splitting fake rates
  mvaFakesMu1 -> SetLineColor(2);
  mvaFakesMu1 -> SetLineWidth(2);
  mvaFakesMu1 -> SetTitle("");
  mvaFakesMu1 -> GetXaxis()->SetTitle("ID BDT output");
  mvaFakesMu2 -> SetLineColor(3);
  mvaFakesMu2 -> SetLineWidth(2);
  mvaFakesMu2 -> SetTitle("");
  mvaFakesMu2 -> GetXaxis()->SetTitle("ID BDT output");
  mvaFakesK   -> SetLineColor(4);
  mvaFakesK   -> SetLineWidth(2);
  mvaFakesK   -> SetTitle("");
  mvaFakesK   -> GetXaxis()->SetTitle("ID BDT output");
  ptFakesMu1  -> SetLineColor(2);
  ptFakesMu1  -> SetLineWidth(2);
  ptFakesMu1  -> SetTitle("");
  ptFakesMu1  -> GetXaxis()->SetTitle("p_{T} [GeV]");
  ptFakesMu2  -> SetLineColor(3);
  ptFakesMu2  -> SetLineWidth(2);
  ptFakesMu2  -> SetTitle("");
  ptFakesMu2  -> GetXaxis()->SetTitle("p_{T} [GeV]");
  ptFakesK    -> SetLineColor(4);
  ptFakesK    -> SetLineWidth(2);
  ptFakesK    -> SetTitle("");
  ptFakesK    -> GetXaxis()->SetTitle("p_{T} [GeV]");
  etaFakesMu1 -> SetLineColor(2);
  etaFakesMu1 -> SetLineWidth(2);
  etaFakesMu1 -> SetTitle("");
  etaFakesMu1 -> GetXaxis()->SetTitle("#eta");
  etaFakesMu2 -> SetLineColor(3);
  etaFakesMu2 -> SetLineWidth(2);
  etaFakesMu2 -> SetTitle("");
  etaFakesMu2 -> GetXaxis()->SetTitle("#eta");
  etaFakesK   -> SetLineColor(4);
  etaFakesK   -> SetLineWidth(2);
  etaFakesK   -> SetTitle("");
  etaFakesK   -> GetXaxis()->SetTitle("#eta");

  TLegend *leg;
  leg = new TLegend(0.55,0.60,0.85,0.85);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.05);
  leg->SetFillColor(0);
  leg->AddEntry(mvaFakesMu1, "1st #mu", "lp");
  leg->AddEntry(mvaFakesMu2, "2nd #mu", "lp");
  leg->AddEntry(mvaFakesK,   "K", "lp");

  TCanvas c1h1("c1h1","",1);
  mvaFakesMu1 -> Rebin();
  mvaFakesMu2 -> Rebin();
  mvaFakesK   -> Rebin();
  mvaFakesMu2 -> DrawNormalized("hist");
  mvaFakesMu1 -> DrawNormalized("samehist");
  mvaFakesK   -> DrawNormalized("samehist");
  leg->Draw();
  c1h1.SaveAs("mvaId_KvsMu.png");

  TCanvas c1h2("c1h2","",1);
  ptFakesMu1 -> Rebin();
  ptFakesMu2 -> Rebin();
  ptFakesK   -> Rebin();
  ptFakesMu2 -> DrawNormalized("hist");
  ptFakesMu1 -> DrawNormalized("samehist");
  ptFakesK   -> DrawNormalized("samehist");
  leg->Draw();
  c1h2.SaveAs("pt_KvsMu.png");

  TCanvas c1h3("c1h3","",1);
  etaFakesMu1 -> Rebin();
  etaFakesMu2 -> Rebin();
  etaFakesK   -> Rebin();
  etaFakesMu1 -> DrawNormalized("hist");
  etaFakesMu2 -> DrawNormalized("samehist");
  etaFakesK   -> DrawNormalized("samehist");
  leg->Draw();
  c1h3.SaveAs("eta_KvsMu.png");

  // -----------------------------
  // PF overlap or not
  if (studyOverlap==1) {

    mvaFakes_LptPfOverlap    -> SetLineWidth(2);  
    mvaFakes_LptPfOverlap    -> SetLineColor(6);  
    ptFakes_LptPfOverlap     -> SetLineWidth(2);  
    ptFakes_LptPfOverlap     -> SetLineColor(6);  
    etaFakes_LptPfOverlap    -> SetLineWidth(2);  
    etaFakes_LptPfOverlap    -> SetLineColor(6);  
    mvaFakes_LptNotPfOverlap -> SetLineWidth(2);  
    mvaFakes_LptNotPfOverlap -> SetLineColor(4);  
    ptFakes_LptNotPfOverlap  -> SetLineWidth(2);  
    ptFakes_LptNotPfOverlap  -> SetLineColor(4);  
    etaFakes_LptNotPfOverlap -> SetLineWidth(2);  
    etaFakes_LptNotPfOverlap -> SetLineColor(4);      
    
    TLegend *legA1;
    legA1 = new TLegend(0.60,0.65,0.90,0.90);
    legA1->SetFillStyle(0);
    legA1->SetBorderSize(0);
    legA1->SetTextSize(0.05);
    legA1->SetFillColor(0);
    legA1->AddEntry(mvaFakes_LptNotPfOverlap, "No PF overlap", "lp");
    legA1->AddEntry(mvaFakes_LptPfOverlap,    "PF overlap", "lp");
    //
    TLegend *legB1;
    legB1 = new TLegend(0.60,0.65,0.90,0.90);
    legB1->SetFillStyle(0);
    legB1->SetBorderSize(0);
    legB1->SetTextSize(0.05);
    legB1->SetFillColor(0);
    legB1->AddEntry(mvaFakes_LptNotPfOverlap, "No PF overlap", "lp");
    legB1->AddEntry(mvaFakes_LptPfOverlap,    "PF overlap", "lp");
    //
    TLegend *legC1;
    legC1 = new TLegend(0.20,0.10,0.60,0.25);
    legC1->SetFillStyle(0);
    legC1->SetBorderSize(0);
    legC1->SetTextSize(0.05);
    legC1->SetFillColor(0);
    legC1->AddEntry(mvaFakes_LptNotPfOverlap, "No PF overlap", "lp");
    legC1->AddEntry(mvaFakes_LptPfOverlap,    "PF overlap", "lp");

    TCanvas ca1("ca1","ca1",1);
    mvaFakes_LptNotPfOverlap -> Rebin(4);
    mvaFakes_LptPfOverlap    -> Rebin(4);
    mvaFakes_LptNotPfOverlap -> SetTitle(""); 
    mvaFakes_LptPfOverlap    -> SetTitle(""); 
    mvaFakes_LptNotPfOverlap -> GetXaxis()->SetTitle("Id BDT");
    mvaFakes_LptPfOverlap    -> GetXaxis()->SetTitle("Id BDT");
    mvaFakes_LptNotPfOverlap -> DrawNormalized("histe");
    mvaFakes_LptPfOverlap    -> DrawNormalized("samehiste");
    legA1->Draw();
    ca1.SaveAs("outputBDT_fakes_w-wo-overlap.png");
    //
    TCanvas ca2("ca2","ca2",1);
    ptFakes_LptNotPfOverlap -> Rebin(4);
    ptFakes_LptPfOverlap    -> Rebin(4);
    ptFakes_LptNotPfOverlap -> SetTitle(""); 
    ptFakes_LptPfOverlap    -> SetTitle(""); 
    ptFakes_LptNotPfOverlap -> GetXaxis()->SetTitle("pt [GeV]");
    ptFakes_LptPfOverlap    -> GetXaxis()->SetTitle("pt [GeV]");
    ptFakes_LptNotPfOverlap -> DrawNormalized("histe");
    ptFakes_LptPfOverlap    -> DrawNormalized("samehiste");
    legB1->Draw();
    ca2.SaveAs("pt_fakes_w-wo-overlap.png");
    //
    TCanvas ca3("ca3","ca3",1);
    etaFakes_LptNotPfOverlap -> Rebin(4);
    etaFakes_LptPfOverlap    -> Rebin(4);
    etaFakes_LptNotPfOverlap -> SetTitle(""); 
    etaFakes_LptPfOverlap    -> SetTitle(""); 
    etaFakes_LptNotPfOverlap -> GetXaxis()->SetTitle("#eta");
    etaFakes_LptPfOverlap    -> GetXaxis()->SetTitle("#eta");
    etaFakes_LptPfOverlap    -> DrawNormalized("histe");
    etaFakes_LptNotPfOverlap -> DrawNormalized("samehiste");
    legC1->Draw();
    ca3.SaveAs("eta_fakes_w-wo-overlap.png");
  }    
 
}
