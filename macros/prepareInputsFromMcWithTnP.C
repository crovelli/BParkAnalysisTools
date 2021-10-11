#define prepareInputsFromMcWithTnP_cxx
#include "prepareInputsFromMcWithTnP.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <iostream>


// To be run on MC tnp formatted ntuples to get 
// - Id distribution before/after weighting the kinematics to match kine without TnP
// - All electron-related distributions for comparison with sPlots
// Select signal and background based on match with MC-truth after the TnP selection 

using namespace std;

void prepareInputsFromMcWithTnP::Loop(bool applyWeight, bool testLPT)
{
  if (fChain == 0) return;

  // -----------------------------------------------------------------------
  // To compute weights
  TH2F *probePtVsEtaSignalMc   = new TH2F("probePtVsEtaSignalMc","probePtVsEtaSignalMc",     40, -2.4, 2.4, 60, 0., 15.);
  TH2F *probePtVsEtaSignalMcWW = new TH2F("probePtVsEtaSignalMcWW","probePtVsEtaSignalMcWW", 40, -2.4, 2.4, 60, 0., 15.);
  probePtVsEtaSignalMc->Sumw2();
  probePtVsEtaSignalMcWW->Sumw2();
  
  // Full pT range: distributions to compute/test weights and for data/MC
  // 
  // No weight
  TH1F *probePtSignalMc  = new TH1F("probePtSignalMc",  "probePtSignalMc",  90,  0.,  30.);
  TH1F *probePtFakeMc    = new TH1F("probePtFakeMc",    "probePtFakeMc",    90,  0.,  30.);
  TH1F *probeEtaSignalMc = new TH1F("probeEtaSignalMc", "probeEtaSignalMc", 40, -2.4, 2.4);
  TH1F *probeEtaFakeMc   = new TH1F("probeEtaFakeMc",   "probeEtaFakeMc",   40, -2.4, 2.4);
  TH1F *probeMvaSignalMc, *probeMvaFakeMc;
  if (testLPT==1) {
    probeMvaSignalMc = new TH1F("probeMvaSignalMc", "probeMvaSignalMc", 40, -4., 16.);    
    probeMvaFakeMc   = new TH1F("probeMvaFakeMc",   "probeMvaFakeMc",   40, -4., 16.);
  } else {
    probeMvaSignalMc = new TH1F("probeMvaSignalMc", "probeMvaSignalMc", 44, -12., 10.);    
    probeMvaFakeMc   = new TH1F("probeMvaFakeMc",   "probeMvaFakeMc",   44, -12., 10.);
  }
  TH1F *probeDxysigSignalMc   = new TH1F("probeDxysigSignalMc",   "probeDxysigSignalMc",   120, -30., 30.);  
  TH1F *probeDxysigFakeMc     = new TH1F("probeDxysigFakeMc",     "probeDxysigFakeMc",     120, -30., 30.);  
  TH1F *probeDztrgSignalMc    = new TH1F("probeDztrgSignalMc",    "probeDztrgSignalMc",    100, -1., 1.); 
  TH1F *probeDztrgFakeMc      = new TH1F("probeDztrgFakeMc",      "probeDztrgFakeMc",      100, -1., 1.);  
  TH1F *probeIso04relSignalMc = new TH1F("probeIso04relSignalMc", "probeIso04relSignalMc",  50,  0., 50.);  
  TH1F *probeIso04relFakeMc   = new TH1F("probeIso04relFakeMc",   "probeIso04relFakeMc",    50,  0., 50.);  
  probePtSignalMc->Sumw2();
  probePtFakeMc->Sumw2();
  probeEtaSignalMc->Sumw2();
  probeEtaFakeMc->Sumw2();
  probeMvaSignalMc->Sumw2(); 
  probeMvaFakeMc->Sumw2();  
  probeDxysigSignalMc->Sumw2();
  probeDxysigFakeMc->Sumw2();
  probeDztrgSignalMc->Sumw2();
  probeDztrgFakeMc->Sumw2();
  probeIso04relSignalMc->Sumw2();
  probeIso04relFakeMc->Sumw2();
  //
  // With weight (kine only)
  TH1F *probePtSignalMcWW  = new TH1F("probePtSignalMcWW",  "probePtSignalMcWW",  90,  0.,  30.);
  TH1F *probeEtaSignalMcWW = new TH1F("probeEtaSignalMcWW", "probeEtaSignalMcWW", 40, -2.4, 2.4);
  TH1F *probeMvaSignalMcWW;
  if (testLPT==1) {
    probeMvaSignalMcWW = new TH1F("probeMvaSignalMcWW", "probeMvaSignalMcWW", 40, -4., 16.); 
  } else {
    probeMvaSignalMcWW = new TH1F("probeMvaSignalMcWW", "probeMvaSignalMcWW", 44, -12., 10.); 
  }
  TH1F *probeDxysigSignalMcWW   = new TH1F("probeDxysigSignalMcWW",   "probeDxysigSignalMcWW",   120, -30., 30.);  
  TH1F *probeDztrgSignalMcWW    = new TH1F("probeDztrgSignalMcWW",    "probeDztrgSignalMcWW",    100, -1., 1.); 
  TH1F *probeIso04relSignalMcWW = new TH1F("probeIso04relSignalMcWW", "probeIso04relSignalMcWW",  50,  0., 50.);  
  probePtSignalMcWW->Sumw2();
  probeEtaSignalMcWW->Sumw2();
  probeMvaSignalMcWW->Sumw2(); 
  probeDxysigSignalMcWW->Sumw2();
  probeDztrgSignalMcWW->Sumw2();
  probeIso04relSignalMcWW->Sumw2();

  // -----------------------------------------------------------------------
  // Many pT/eta bins (no weight): ID output
  TH1F *mvaSignalEBMc0, *mvaSignalEBMc1, *mvaSignalEBMc2, *mvaSignalEBMc3;
  TH1F *mvaSignalEEMc0, *mvaSignalEEMc1;
  if (testLPT==1) {
    mvaSignalEBMc0 = new TH1F("mvaSignalEBMc0", "mvaSignalEBMc0", 40, -4., 16.);  // pT: 1.0-1.5
    mvaSignalEBMc1 = new TH1F("mvaSignalEBMc1", "mvaSignalEBMc1", 40, -4., 16.);  // 1.5-2.0
    mvaSignalEBMc2 = new TH1F("mvaSignalEBMc2", "mvaSignalEBMc2", 40, -4., 16.);  // 2.0-5.0
    mvaSignalEBMc3 = new TH1F("mvaSignalEBMc3", "mvaSignalEBMc3", 40, -4., 16.);  // >5
    mvaSignalEEMc0 = new TH1F("mvaSignalEEMc0", "mvaSignalEEMc0", 40, -4., 16.);  // pT: 2.0-5.0  
    mvaSignalEEMc1 = new TH1F("mvaSignalEEMc1", "mvaSignalEEMc1", 40, -4., 16.);  // >5
  } else {
    mvaSignalEBMc0 = new TH1F("mvaSignalEBMc0", "mvaSignalEBMc0", 44, -12., 10.);  // 2.0-5.0
    mvaSignalEBMc1 = new TH1F("mvaSignalEBMc1", "mvaSignalEBMc1", 44, -12., 10.);  // >5
    mvaSignalEEMc0 = new TH1F("mvaSignalEEMc0", "mvaSignalEEMc0", 44, -12., 10.);  // pT: 2.0-5.0  
    mvaSignalEEMc1 = new TH1F("mvaSignalEEMc1", "mvaSignalEEMc1", 44, -12., 10.);  // >5
  }

  mvaSignalEBMc0->Sumw2();
  mvaSignalEBMc1->Sumw2();
  mvaSignalEEMc0->Sumw2();
  mvaSignalEEMc1->Sumw2();
  if (testLPT==1) {
    mvaSignalEBMc2->Sumw2();
    mvaSignalEBMc3->Sumw2();
  }
  //
  TH1F *mvaFakeEBMc0, *mvaFakeEBMc1, *mvaFakeEBMc2, *mvaFakeEBMc3;
  TH1F *mvaFakeEEMc0, *mvaFakeEEMc1;
  if (testLPT==1) {
    mvaFakeEBMc0 = new TH1F("mvaFakeEBMc0", "mvaFakeEBMc0", 40, -4., 16.);  // pT: 1.0-1.5
    mvaFakeEBMc1 = new TH1F("mvaFakeEBMc1", "mvaFakeEBMc1", 40, -4., 16.);  // 1.5-2.0
    mvaFakeEBMc2 = new TH1F("mvaFakeEBMc2", "mvaFakeEBMc2", 40, -4., 16.);  // 2.0-5.0
    mvaFakeEBMc3 = new TH1F("mvaFakeEBMc3", "mvaFakeEBMc3", 40, -4., 16.);  // >5
    mvaFakeEEMc0 = new TH1F("mvaFakeEEMc0", "mvaFakeEEMc0", 40, -4., 16.);  // pT: 2.0-5.0
    mvaFakeEEMc1 = new TH1F("mvaFakeEEMc1", "mvaFakeEEMc1", 40, -4., 16.);  // >5
  } else {
    mvaFakeEBMc0 = new TH1F("mvaFakeEBMc0", "mvaFakeEBMc0", 44, -12., 10.);  // 2.0-5.0
    mvaFakeEBMc1 = new TH1F("mvaFakeEBMc1", "mvaFakeEBMc1", 44, -12., 10.);  // >5
    mvaFakeEEMc0 = new TH1F("mvaFakeEEMc0", "mvaFakeEEMc0", 44, -12., 10.);  // pT: 2.0-5.0
    mvaFakeEEMc1 = new TH1F("mvaFakeEEMc1", "mvaFakeEEMc1", 44, -12., 10.);  // >5
  }
  mvaFakeEBMc0->Sumw2();
  mvaFakeEBMc1->Sumw2();
  mvaFakeEEMc0->Sumw2();
  mvaFakeEEMc1->Sumw2();
  if (testLPT==1) {
    mvaFakeEBMc2->Sumw2();
    mvaFakeEBMc3->Sumw2();
  }
  
  // -----------------------------------------------------------------------
  // Many pT/eta bins (with weight): ID output
  TH1F *mvaSignalEBMc0WW, *mvaSignalEBMc1WW, *mvaSignalEBMc2WW, *mvaSignalEBMc3WW;
  TH1F *mvaSignalEEMc0WW, *mvaSignalEEMc1WW;
  if (testLPT==0) {
    mvaSignalEBMc0WW = new TH1F("mvaSignalEBMc0WW", "mvaSignalEBMc0WW", 44, -12., 10.);  // pT: 1.0-1.5
    mvaSignalEBMc1WW = new TH1F("mvaSignalEBMc1WW", "mvaSignalEBMc1WW", 44, -12., 10.);  // 1.5-2.0
    mvaSignalEBMc2WW = new TH1F("mvaSignalEBMc2WW", "mvaSignalEBMc2WW", 44, -12., 10.);  // 2.0-5.0
    mvaSignalEBMc3WW = new TH1F("mvaSignalEBMc3WW", "mvaSignalEBMc3WW", 44, -12., 10.);  // >5
    mvaSignalEEMc0WW = new TH1F("mvaSignalEEMc0WW", "mvaSignalEEMc0WW", 44, -12., 10.);  // pT: 2.0-5.0 
    mvaSignalEEMc1WW = new TH1F("mvaSignalEEMc1WW", "mvaSignalEEMc1WW", 44, -12., 10.);  // >5
  } else {
    mvaSignalEBMc0WW = new TH1F("mvaSignalEBMc0WW", "mvaSignalEBMc0WW", 40, -4., 16.);  // 2.0-5.0
    mvaSignalEBMc1WW = new TH1F("mvaSignalEBMc1WW", "mvaSignalEBMc1WW", 40, -4., 16.);  // >5
    mvaSignalEEMc0WW = new TH1F("mvaSignalEEMc0WW", "mvaSignalEEMc0WW", 40, -4., 16.);  // pT: 2.0-5.0 
    mvaSignalEEMc1WW = new TH1F("mvaSignalEEMc1WW", "mvaSignalEEMc1WW", 40, -4., 16.);  // >5
  }
  mvaSignalEBMc0WW->Sumw2();
  mvaSignalEBMc1WW->Sumw2();
  mvaSignalEEMc0WW->Sumw2();
  mvaSignalEEMc1WW->Sumw2();
  if (testLPT==1) {
    mvaSignalEBMc2WW->Sumw2();
    mvaSignalEBMc3WW->Sumw2();
  }

  // -----------------------------------------------------------------------
  // Many pT/eta bins (no weight): dxy significance
  TH1F *dxysigSignalEBMc0 = new TH1F("dxysigSignalEBMc0", "dxysigSignalEBMc0", 120, -30., 30.); 
  TH1F *dxysigSignalEBMc1 = new TH1F("dxysigSignalEBMc1", "dxysigSignalEBMc1", 120, -30., 30.);
  TH1F *dxysigSignalEEMc0 = new TH1F("dxysigSignalEEMc0", "dxysigSignalEEMc0", 120, -30., 30.); 
  TH1F *dxysigSignalEEMc1 = new TH1F("dxysigSignalEEMc1", "dxysigSignalEEMc1", 120, -30., 30.); 
  dxysigSignalEBMc0->Sumw2();
  dxysigSignalEBMc1->Sumw2();
  dxysigSignalEEMc0->Sumw2();
  dxysigSignalEEMc1->Sumw2();

  TH1F *dxysigSignalEBMc2, *dxysigSignalEBMc3;
  if (testLPT==1) {
    dxysigSignalEBMc2 = new TH1F("dxysigSignalEBMc2", "dxysigSignalEBMc2", 120, -30., 30.); 
    dxysigSignalEBMc3 = new TH1F("dxysigSignalEBMc3", "dxysigSignalEBMc3", 120, -30., 30.); 
    dxysigSignalEBMc2->Sumw2();
    dxysigSignalEBMc3->Sumw2();
  }

  // -----------------------------------------------------------------------
  // Many pT/eta bins (no weight): dzTrg
  TH1F *dztrgSignalEBMc0 = new TH1F("dztrgSignalEBMc0", "dztrgSignalEBMc0", 100, -1., 1.); 
  TH1F *dztrgSignalEBMc1 = new TH1F("dztrgSignalEBMc1", "dztrgSignalEBMc1", 100, -1., 1.); 
  TH1F *dztrgSignalEEMc0 = new TH1F("dztrgSignalEEMc0", "dztrgSignalEEMc0", 100, -1., 1.); 
  TH1F *dztrgSignalEEMc1 = new TH1F("dztrgSignalEEMc1", "dztrgSignalEEMc1", 100, -1., 1.); 
  dztrgSignalEBMc0->Sumw2();
  dztrgSignalEBMc1->Sumw2();
  dztrgSignalEEMc0->Sumw2();
  dztrgSignalEEMc1->Sumw2();

  TH1F *dztrgSignalEBMc2, *dztrgSignalEBMc3;
  if (testLPT==1) {
    dztrgSignalEBMc2 = new TH1F("dztrgSignalEBMc2", "dztrgSignalEBMc2", 100, -1., 1.); 
    dztrgSignalEBMc3 = new TH1F("dztrgSignalEBMc3", "dztrgSignalEBMc3", 100, -1., 1.);
    dztrgSignalEBMc2->Sumw2();
    dztrgSignalEBMc3->Sumw2();
  }

  // -----------------------------------------------------------------------
  // Many pT/eta bins (no weight): iso04rel
  TH1F *iso04relSignalEBMc0 = new TH1F("iso04relSignalEBMc0", "iso04relSignalEBMc0", 50, 0., 50.); 
  TH1F *iso04relSignalEBMc1 = new TH1F("iso04relSignalEBMc1", "iso04relSignalEBMc1", 50, 0., 50.); 
  TH1F *iso04relSignalEEMc0 = new TH1F("iso04relSignalEEMc0", "iso04relSignalEEMc0", 50, 0., 50.); 
  TH1F *iso04relSignalEEMc1 = new TH1F("iso04relSignalEEMc1", "iso04relSignalEEMc1", 50, 0., 50.); 
  iso04relSignalEBMc0->Sumw2();
  iso04relSignalEBMc1->Sumw2();
  iso04relSignalEEMc0->Sumw2();
  iso04relSignalEEMc1->Sumw2();

  TH1F *iso04relSignalEBMc2, *iso04relSignalEBMc3;
  if (testLPT==1) {
    iso04relSignalEBMc2 = new TH1F("iso04relSignalEBMc2", "iso04relSignalEBMc2", 50, 0., 50.); 
    iso04relSignalEBMc3 = new TH1F("iso04relSignalEBMc3", "iso04relSignalEBMc3", 50, 0., 50.); 
    iso04relSignalEBMc2->Sumw2();
    iso04relSignalEBMc3->Sumw2();
  }

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

  float weightsSignal[100][100];
  for (int iBinEta=0; iBinEta<100; iBinEta++) { 
    for (int iBinPt=0; iBinPt<100; iBinPt++) { 
      weightsSignal[iBinEta][iBinPt]=1.;
    }
  }

  if (applyWeight) {
    //TFile fileWeight("files_marchNoReg2/probeLowPt/weightFile_tnpVsNani.root"); 
    TFile fileWeight("files_marchNoReg2_TC/probeLowPt/weightFile_tnpVsNani.root"); 
    TH2F *signalWeights = (TH2F*)fileWeight.Get("ptVsEtaSignalWeights"); 
    nBinsWPt  = signalWeights->GetNbinsY(); 
    nBinsWEta = signalWeights->GetNbinsX(); 
    
    for (int iBinEta=0; iBinEta<=nBinsWEta; iBinEta++) {
      for (int iBinPt=0; iBinPt<=nBinsWPt; iBinPt++) {
	minBPt[iBinPt]   = signalWeights->GetYaxis()->GetBinLowEdge(iBinPt);
	maxBPt[iBinPt]   = signalWeights->GetYaxis()->GetBinUpEdge(iBinPt);
	minBEta[iBinEta] = signalWeights->GetXaxis()->GetBinLowEdge(iBinEta);
	maxBEta[iBinEta] = signalWeights->GetXaxis()->GetBinUpEdge(iBinEta);
	weightsSignal[iBinEta][iBinPt] = signalWeights->GetBinContent(iBinEta,iBinPt);    // 0=underflow, N=overflow
      }
    }

    // over/under flow ranges
    minBPt[0]  = -999.;
    minBEta[0] = -999.;
    maxBPt[nBinsWPt] = 999.;
    maxBEta[nBinsWEta] = 999.;

  } // apply weight


  // --------------------------------------------------------
  // Loop over entries
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    // HLT
    if (hlt_9ip6==0) continue;
    
    // Acceptance
    if (fabs(probeEta)>2.4) continue;
    if (probePt<1.0)        continue;

    // LowPt vs PF selection already applied in formatted-tnp ntuples

    // Fitting range applied in data
    if (pair_mass<2.6 || pair_mass>3.4) continue;

    // Apply same loose cuts as done for data - chiara
    if (testLPT==1) {
      if (probeMvaId<-4 || probeMvaId>16) continue;  
    } else {
      if (probePfmvaId<-12 || probePfmvaId>10) continue;    
    }
    if (probePt>1000) continue;
    if (fabs(probeEta)>2.4) continue;  
    if (fabs(probeDxySig)>200) continue;   
    if (fabs(probeDzTrg)>10) continue;   
    if (probeIso04Rel<-1 || probeIso04Rel>5000) continue; 

    // pT vs eta weight
    float ptEtaSignalWeight=1.;

    if (applyWeight) {
      for (int iBinEta=0; iBinEta<=nBinsWEta; iBinEta++) {
	for (int iBinPt=0; iBinPt<=nBinsWPt; iBinPt++) {
	  bool thisBin = false;
	  if (probePt>=minBPt[iBinPt] && probePt<maxBPt[iBinPt] && probeEta>=minBEta[iBinEta] && probeEta<maxBEta[iBinEta]) thisBin = true;
	  if (probeMatchMcFromJPsi==1 && thisBin) ptEtaSignalWeight = weightsSignal[iBinEta][iBinPt];
	}
      }
    }

    // matching mc-truth in MC: full pT range before weight
    if (probeMatchMcFromJPsi==1) { 
      if (testLPT==1) probeMvaSignalMc->Fill(probeMvaId, weight);
      if (testLPT==0) probeMvaSignalMc->Fill(probePfmvaId, weight);
      probePtSignalMc -> Fill(probePt, weight);
      probeEtaSignalMc->Fill(probeEta, weight);
      probePtVsEtaSignalMc->Fill(probeEta,probePt, weight);
      probeDxysigSignalMc -> Fill(probeDxySig, weight);
      probeDztrgSignalMc -> Fill(probeDzTrg, weight);
      probeIso04relSignalMc -> Fill(probeIso04Rel, weight);
    }
    
    // matching mc-truth in MC: full pT range after weight
    if (probeMatchMcFromJPsi==1) { 
      if (testLPT==1) probeMvaSignalMcWW->Fill(probeMvaId, ptEtaSignalWeight*weight);       
      if (testLPT==0) probeMvaSignalMcWW->Fill(probePfmvaId, ptEtaSignalWeight*weight);       
      probePtSignalMcWW -> Fill(probePt, ptEtaSignalWeight*weight);
      probeEtaSignalMcWW->Fill(probeEta, ptEtaSignalWeight*weight);
      probePtVsEtaSignalMcWW->Fill(probeEta,probePt,ptEtaSignalWeight*weight);
      probeDxysigSignalMcWW -> Fill(probeDxySig, ptEtaSignalWeight*weight);
      probeDztrgSignalMcWW -> Fill(probeDzTrg, ptEtaSignalWeight*weight);
      probeIso04relSignalMcWW -> Fill(probeIso04Rel, ptEtaSignalWeight*weight);
    }

    // not matching mc-truth in MC: full pT range before weight
    if (probeMatchMcFromJPsi==0) { 
      if (testLPT==1) probeMvaFakeMc->Fill(probeMvaId, weight);
      if (testLPT==0) probeMvaFakeMc->Fill(probePfmvaId, weight);
      probePtFakeMc -> Fill(probePt, weight);
      probeEtaFakeMc->Fill(probeEta, weight);
      probeDxysigFakeMc -> Fill(probeDxySig, weight);
      probeDztrgFakeMc -> Fill(probeDzTrg, weight);
      probeIso04relFakeMc -> Fill(probeIso04Rel, weight);
    }

    // matching mc-truth in MC: eta/pT bins
    if (probeMatchMcFromJPsi==1) {  // signal

      float theId=-1;
      if (testLPT==1) theId=probeMvaId;
      if (testLPT==0) theId=probePfmvaId;

      if (fabs(probeEta)<1.5) {   // barrel

	if (testLPT==1) {  // LPT barrel	

	  if (probePt>=1.0 && probePt<1.5) {
	    mvaSignalEBMc0      -> Fill(theId, weight);
	    mvaSignalEBMc0WW    -> Fill(theId,ptEtaSignalWeight*weight);
	    dxysigSignalEBMc0   -> Fill(probeDxySig, weight);
	    dztrgSignalEBMc0    -> Fill(probeDzTrg, weight);
	    iso04relSignalEBMc0 -> Fill(probeIso04Rel, weight);
	  }
	  if (probePt>=1.5 && probePt<2.0) {
	    mvaSignalEBMc1      -> Fill(theId, weight);
	    mvaSignalEBMc1WW    -> Fill(theId,ptEtaSignalWeight*weight);
	    dxysigSignalEBMc1   -> Fill(probeDxySig, weight);
	    dztrgSignalEBMc1    -> Fill(probeDzTrg, weight);
	    iso04relSignalEBMc1 -> Fill(probeIso04Rel, weight);
	  }
	  if (probePt>=2.0 && probePt<5.0) {
	    mvaSignalEBMc2      -> Fill(theId, weight);
	    mvaSignalEBMc2WW    -> Fill(theId,ptEtaSignalWeight*weight);
	    dxysigSignalEBMc2   -> Fill(probeDxySig, weight);
	    dztrgSignalEBMc2    -> Fill(probeDzTrg, weight);
	    iso04relSignalEBMc2 -> Fill(probeIso04Rel, weight);
	  }
	  if (probePt>=5.0) {
	    mvaSignalEBMc3      -> Fill(theId, weight);
	    mvaSignalEBMc3WW    -> Fill(theId,ptEtaSignalWeight*weight);
	    dxysigSignalEBMc3   -> Fill(probeDxySig, weight);
	    dztrgSignalEBMc3    -> Fill(probeDzTrg, weight);
	    iso04relSignalEBMc3 -> Fill(probeIso04Rel, weight);
	  }

	} else {  // PF barrel

	  if (probePt>=2.0 && probePt<5.0) {
	    mvaSignalEBMc0      -> Fill(theId, weight);
	    mvaSignalEBMc0WW    -> Fill(theId,ptEtaSignalWeight*weight);
	    dxysigSignalEBMc0   -> Fill(probeDxySig, weight);
	    dztrgSignalEBMc0    -> Fill(probeDzTrg, weight);
	    iso04relSignalEBMc0 -> Fill(probeIso04Rel, weight);
	  }
	  if (probePt>=5.0) {
	    mvaSignalEBMc1      -> Fill(theId, weight);
	    mvaSignalEBMc1WW    -> Fill(theId,ptEtaSignalWeight*weight);
	    dxysigSignalEBMc1   -> Fill(probeDxySig, weight);
	    dztrgSignalEBMc1    -> Fill(probeDzTrg, weight);
	    iso04relSignalEBMc1 -> Fill(probeIso04Rel, weight);
	  }
	}

      } else {  // endcap (PF and LPT)

	if (probePt>=2.0 && probePt<5.0) {
	  mvaSignalEEMc0      -> Fill(theId, weight);
	  mvaSignalEEMc0WW    -> Fill(theId,ptEtaSignalWeight*weight);
	  dxysigSignalEEMc0   -> Fill(probeDxySig, weight);
	  dztrgSignalEEMc0    -> Fill(probeDzTrg, weight);
	  iso04relSignalEEMc0 -> Fill(probeIso04Rel, weight);
	}
	if (probePt>=5.0) {
	  mvaSignalEEMc1      -> Fill(theId, weight);
	  mvaSignalEEMc1WW    -> Fill(theId,ptEtaSignalWeight*weight);
	  dxysigSignalEEMc1   -> Fill(probeDxySig, weight);
	  dztrgSignalEEMc1    -> Fill(probeDzTrg, weight);
	  iso04relSignalEEMc1 -> Fill(probeIso04Rel, weight);
	}
      }
    }

    // not matching mc-truth in MC: eta/pT bins
    if (probeMatchMcFromJPsi==0) { 

      float theId=-1;
      if (testLPT==1) theId=probeMvaId;
      if (testLPT==0) theId=probePfmvaId;

      if (fabs(probeEta)<1.5) {  // barrel
	if (testLPT==1) { 
	  if (probePt>=1.0 && probePt<1.5) mvaFakeEBMc0->Fill(theId, weight);
	  if (probePt>=1.5 && probePt<2.0) mvaFakeEBMc1->Fill(theId, weight);
	  if (probePt>=2.0 && probePt<5.0) mvaFakeEBMc2->Fill(theId, weight);
	  if (probePt>=5.0) mvaFakeEBMc3->Fill(theId, weight);
	} else {
	  if (probePt>=2.0 && probePt<5.0) mvaFakeEBMc0->Fill(theId, weight);
	  if (probePt>=5.0)                mvaFakeEBMc1->Fill(theId, weight);	
	}
      } else {  // endcap
	if (probePt>=2.0 && probePt<5.0) mvaFakeEEMc0->Fill(theId, weight);
	if (probePt>=5.0)                mvaFakeEEMc1->Fill(theId, weight);
      }
    }

  } // loop over entries


  // -----------------------------------------------------------------------
  // Cosmetics
  probeMvaSignalMc -> SetLineWidth(2);
  probeMvaSignalMc -> SetLineColor(6);
  probeMvaFakeMc   -> SetLineWidth(2);
  probeMvaFakeMc   -> SetLineColor(4);
  probePtSignalMc  -> SetLineWidth(2);
  probePtSignalMc  -> SetLineColor(6);
  probePtFakeMc    -> SetLineWidth(2);
  probePtFakeMc    -> SetLineColor(4);
  probeDxysigSignalMc -> SetLineWidth(2); 
  probeDxysigSignalMc -> SetLineColor(6); 
  probeDxysigFakeMc   -> SetLineWidth(2); 
  probeDxysigFakeMc   -> SetLineColor(4); 
  probeDztrgSignalMc  -> SetLineWidth(2); 
  probeDztrgSignalMc  -> SetLineColor(6); 
  probeDztrgFakeMc    -> SetLineWidth(2); 
  probeDztrgFakeMc    -> SetLineColor(4); 
  probeIso04relSignalMc  -> SetLineWidth(2); 
  probeIso04relSignalMc  -> SetLineColor(6); 
  probeIso04relFakeMc    -> SetLineWidth(2); 
  probeIso04relFakeMc    -> SetLineColor(4); 
  probeEtaSignalMc -> SetLineWidth(2);
  probeEtaSignalMc -> SetLineColor(6);
  probeEtaFakeMc   -> SetLineWidth(2);
  probeEtaFakeMc   -> SetLineColor(4);
  //
  mvaSignalEBMc0 -> SetLineWidth(2); 
  mvaSignalEBMc1 -> SetLineWidth(2); 
  mvaSignalEBMc0 -> SetLineColor(6); 
  mvaSignalEBMc1 -> SetLineColor(6); 
  if (testLPT==1) {
    mvaSignalEBMc2 -> SetLineWidth(2); 
    mvaSignalEBMc3 -> SetLineWidth(2); 
    mvaSignalEBMc2 -> SetLineColor(6); 
    mvaSignalEBMc3 -> SetLineColor(6); 
  }
  //
  mvaSignalEBMc0WW -> SetLineWidth(2); 
  mvaSignalEBMc1WW -> SetLineWidth(2); 
  mvaSignalEBMc0WW -> SetLineColor(6); 
  mvaSignalEBMc1WW -> SetLineColor(6); 
  if (testLPT==1) {
    mvaSignalEBMc2WW -> SetLineWidth(2); 
    mvaSignalEBMc3WW -> SetLineWidth(2); 
    mvaSignalEBMc2WW -> SetLineColor(6); 
    mvaSignalEBMc3WW -> SetLineColor(6); 
  }
  //
  mvaSignalEEMc0 -> SetLineWidth(2); 
  mvaSignalEEMc1 -> SetLineWidth(2); 
  mvaSignalEEMc0 -> SetLineColor(6); 
  mvaSignalEEMc1 -> SetLineColor(6); 
  //
  mvaSignalEEMc0WW -> SetLineWidth(2); 
  mvaSignalEEMc1WW -> SetLineWidth(2); 
  mvaSignalEEMc0WW -> SetLineColor(6); 
  mvaSignalEEMc1WW -> SetLineColor(6); 
  //
  mvaFakeEBMc0 -> SetLineWidth(2); 
  mvaFakeEBMc1 -> SetLineWidth(2); 
  mvaFakeEBMc0 -> SetLineColor(4); 
  mvaFakeEBMc1 -> SetLineColor(4); 
  if (testLPT==1) {
    mvaFakeEBMc2 -> SetLineWidth(2); 
    mvaFakeEBMc3 -> SetLineWidth(2); 
    mvaFakeEBMc2 -> SetLineColor(4); 
    mvaFakeEBMc3 -> SetLineColor(4); 
  }
  //
  mvaFakeEEMc0 -> SetLineWidth(2); 
  mvaFakeEEMc1 -> SetLineWidth(2); 
  mvaFakeEEMc0 -> SetLineColor(4); 
  mvaFakeEEMc1 -> SetLineColor(4); 


  // -----------------------------------------------------------------------
  // Save outputs
  TFile myFile("myFileMcAfterTnp.root","RECREATE");
  myFile.cd();
  //
  probePtSignalMc->Write();
  probePtFakeMc->Write();
  probeEtaSignalMc->Write();
  probeEtaFakeMc->Write();
  probeMvaSignalMc->Write(); 
  probeMvaFakeMc->Write();  
  probeDxysigSignalMc   -> Write(); 
  probeDxysigFakeMc     -> Write();
  probeDztrgSignalMc    -> Write();
  probeDztrgFakeMc      -> Write();
  probeIso04relSignalMc -> Write();
  probeIso04relFakeMc   -> Write();
  //
  probePtSignalMcWW->Write();
  probeEtaSignalMcWW->Write();
  probeMvaSignalMcWW->Write(); 
  probeDxysigSignalMcWW->Write(); 
  probeDztrgSignalMcWW->Write();
  probeIso04relSignalMcWW->Write();
  //
  probePtVsEtaSignalMc->Write();
  probePtVsEtaSignalMcWW->Write();
  //
  mvaSignalEBMc0->Write();
  mvaSignalEBMc1->Write();
  if (testLPT==1) {
    mvaSignalEBMc2->Write();
    mvaSignalEBMc3->Write();
  }
  mvaSignalEEMc0->Write();
  mvaSignalEEMc1->Write();
  //
  mvaSignalEBMc0WW->Write();
  mvaSignalEBMc1WW->Write();
  if (testLPT==1) {
    mvaSignalEBMc2WW->Write();
    mvaSignalEBMc3WW->Write();
  }
  mvaSignalEEMc0WW->Write();
  mvaSignalEEMc1WW->Write();
  //
  mvaFakeEBMc0->Write();
  mvaFakeEBMc1->Write();
  if (testLPT==1) {
    mvaFakeEBMc2->Write();
    mvaFakeEBMc3->Write();
  }
  mvaFakeEEMc0->Write();
  mvaFakeEEMc1->Write();
  //
  dxysigSignalEBMc0->Write();
  dxysigSignalEBMc1->Write();
  if (testLPT==1) {
    dxysigSignalEBMc2->Write();
    dxysigSignalEBMc3->Write();
  }
  dxysigSignalEEMc0->Write();
  dxysigSignalEEMc1->Write();
  //
  dztrgSignalEBMc0->Write();
  dztrgSignalEBMc1->Write();
  if (testLPT==1) {
    dztrgSignalEBMc2->Write();
    dztrgSignalEBMc3->Write();
  }
  dztrgSignalEEMc0->Write();
  dztrgSignalEEMc1->Write();
  //
  iso04relSignalEBMc0->Write();
  iso04relSignalEBMc1->Write();
  if (testLPT==1) {
    iso04relSignalEBMc2->Write();
    iso04relSignalEBMc3->Write();
  }
  iso04relSignalEEMc0->Write();
  iso04relSignalEEMc1->Write();
  //

  myFile.Close();


  // -----------------------------------------------------------------------
  // Rebin
  mvaSignalEBMc0->Rebin();
  mvaSignalEBMc1->Rebin();
  if (testLPT==1) {
    mvaSignalEBMc2->Rebin();
    mvaSignalEBMc3->Rebin();
  }
  mvaSignalEEMc0->Rebin();
  mvaSignalEEMc1->Rebin();

  mvaFakeEBMc0->Rebin();
  mvaFakeEBMc1->Rebin();
  if (testLPT==1) {
    mvaFakeEBMc2->Rebin();
    mvaFakeEBMc3->Rebin();
  }
  mvaFakeEEMc0->Rebin();
  mvaFakeEEMc1->Rebin();


  // -----------------------------------------------------------------------
  // Plots
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);

  TLegend *leg;
  leg = new TLegend(0.10,0.65,0.45,0.90);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.05);
  leg->SetFillColor(0);
  leg->AddEntry(probeMvaSignalMc, "Matching", "lp");
  leg->AddEntry(probeMvaFakeMc,   "Not matching", "lp");
  //
  TLegend *legB;
  legB = new TLegend(0.60,0.65,0.95,0.90);
  legB->SetFillStyle(0);
  legB->SetBorderSize(0);
  legB->SetTextSize(0.05);
  legB->SetFillColor(0);
  legB->AddEntry(probeMvaSignalMc, "Matching", "lp");
  legB->AddEntry(probeMvaFakeMc,   "Not matching", "lp");
  //
  TLegend *legC;
  legC = new TLegend(0.25,0.15,0.65,0.30);
  legC->SetFillStyle(0);
  legC->SetBorderSize(0);
  legC->SetTextSize(0.05);
  legC->SetFillColor(0);
  legC->AddEntry(probeMvaSignalMc, "Matching", "lp");
  legC->AddEntry(probeMvaFakeMc,   "Not matching", "lp");

  TCanvas cmvamc("cmvamc","cmvamc",1);
  probeMvaSignalMc->SetTitle("");
  probeMvaFakeMc->SetTitle("");
  probeMvaSignalMc->GetXaxis()->SetTitle("Id BDT");
  probeMvaFakeMc->GetXaxis()->SetTitle("Id BDT");
  probeMvaFakeMc->DrawNormalized("hist");
  probeMvaSignalMc->DrawNormalized("samehist");
  legB->Draw();
  cmvamc.SaveAs("outputBDT_matchVsFake_withTnP.png");
 
  TCanvas cptmc("cptmc","cptmc",1);
  probePtSignalMc->SetTitle("");
  probePtFakeMc->SetTitle("");
  probePtSignalMc->GetXaxis()->SetTitle("pT [GeV]");
  probePtFakeMc->GetXaxis()->SetTitle("pT [GeV]");
  probePtFakeMc->DrawNormalized("hist");
  probePtSignalMc->DrawNormalized("samehist");
  legB->Draw();
  cptmc.SaveAs("pt_matchVsFake_withTnP.png");

  TCanvas cetamc("cetamc","cetamc",1);
  probeEtaSignalMc->SetTitle("");
  probeEtaFakeMc->SetTitle("");
  probeEtaSignalMc->GetXaxis()->SetTitle("#eta");
  probeEtaFakeMc->GetXaxis()->SetTitle("#eta");
  probeEtaSignalMc->DrawNormalized("hist");
  probeEtaFakeMc->DrawNormalized("samehist");
  legC->Draw();
  cetamc.SaveAs("eta_matchVsFake_withTnP.png");

  TCanvas cdxymc("cdxymc","cdxymc",1);
  probeDxysigSignalMc->SetTitle("");
  probeDxysigFakeMc->SetTitle("");
  probeDxysigSignalMc->GetXaxis()->SetTitle("dXY sig");
  probeDxysigFakeMc->GetXaxis()->SetTitle("dXY sig");
  probeDxysigFakeMc->DrawNormalized("hist");
  probeDxysigSignalMc->DrawNormalized("samehist");
  legB->Draw();
  cdxymc.SaveAs("dxySig_matchVsFake_withTnP.png");

  TCanvas cdztmc("cdztmc","cdztmc",1);
  probeDztrgSignalMc->SetTitle("");
  probeDztrgFakeMc->SetTitle("");
  probeDztrgSignalMc->GetXaxis()->SetTitle("dz trg");
  probeDztrgFakeMc->GetXaxis()->SetTitle("dz trg");
  probeDztrgFakeMc->DrawNormalized("hist");
  probeDztrgSignalMc->DrawNormalized("samehist");
  legB->Draw();
  cdztmc.SaveAs("dzTrg_matchVsFake_withTnP.png");

  TCanvas cis04mc("cis04mc","cis04mc",1);
  probeIso04relSignalMc->SetTitle("");
  probeIso04relFakeMc->SetTitle("");
  probeIso04relSignalMc->GetXaxis()->SetTitle("Rel Iso 04");
  probeIso04relFakeMc->GetXaxis()->SetTitle("Rel Iso 04");
  probeIso04relSignalMc->DrawNormalized("hist");
  probeIso04relFakeMc->DrawNormalized("samehist");
  legB->Draw();
  cis04mc.SaveAs("iso04rel_matchVsFake_withTnP.png");

  TCanvas cmvaeb0("cmvaeb0","cmvaeb0",1);
  if (testLPT==1) {
    mvaSignalEBMc0 -> SetTitle("EB, 1.0 < pT < 1.5");
    mvaFakeEBMc0   -> SetTitle("EB, 1.0 < pT < 1.5");
  } else {
    mvaSignalEBMc0 -> SetTitle("EB, 2.0 < pT < 5.0");
    mvaFakeEBMc0   -> SetTitle("EB, 2.0 < pT < 5.0");
  }
  mvaSignalEBMc0 -> GetXaxis()->SetTitle("Id BDT");
  mvaFakeEBMc0   -> GetXaxis()->SetTitle("Id BDT");
  mvaFakeEBMc0   -> DrawNormalized("hist");
  mvaSignalEBMc0 -> DrawNormalized("samehist");
  legB->Draw();
  cmvaeb0.SaveAs("outputBDT_matchVsFake_withTnP_EB0.png");  
  //
  TCanvas cmvaeb1("cmvaeb1","cmvaeb1",1);
  if (testLPT==1) {
    mvaSignalEBMc1 -> SetTitle("EB, 1.5 < pT < 2.0");
    mvaFakeEBMc1   -> SetTitle("EB, 1.5 < pT < 2.0");
  } else {
    mvaSignalEBMc1 -> SetTitle("EB, pT >= 5");
    mvaFakeEBMc1   -> SetTitle("EB, pT >= 5");
  }
  mvaSignalEBMc1 -> GetXaxis()->SetTitle("Id BDT");
  mvaFakeEBMc1   -> GetXaxis()->SetTitle("Id BDT");
  mvaFakeEBMc1   -> DrawNormalized("hist");
  mvaSignalEBMc1 -> DrawNormalized("samehist");
  legB->Draw();
  cmvaeb1.SaveAs("outputBDT_matchVsFake_withTnP_EB1.png");  
  //
  if (testLPT==1) {
    TCanvas cmvaeb2("cmvaeb2","cmvaeb2",1);
    mvaSignalEBMc2 -> SetTitle("EB, 2.0 < pT < 5.0");
    mvaFakeEBMc2   -> SetTitle("EB, 2.0 < pT < 5.0");
    mvaSignalEBMc2 -> GetXaxis()->SetTitle("Id BDT");
    mvaFakeEBMc2   -> GetXaxis()->SetTitle("Id BDT");
    mvaFakeEBMc2   -> DrawNormalized("hist");
    mvaSignalEBMc2 -> DrawNormalized("samehist");
    legB->Draw();
    cmvaeb2.SaveAs("outputBDT_matchVsFake_withTnP_EB2.png");  
    //
    TCanvas cmvaeb3("cmvaeb3","cmvaeb3",1);
    mvaSignalEBMc3 -> SetTitle("EB, pT >= 5");
    mvaFakeEBMc3   -> SetTitle("EB, pT >= 5");
    mvaSignalEBMc3 -> GetXaxis()->SetTitle("Id BDT");
    mvaFakeEBMc3   -> GetXaxis()->SetTitle("Id BDT");
    mvaSignalEBMc3 -> DrawNormalized("hist");
    mvaFakeEBMc3   -> DrawNormalized("samehist");
    legB->Draw();
    cmvaeb3.SaveAs("outputBDT_matchVsFake_withTnP_EB3.png");  
  }

  TCanvas cmvaee0("cmvaee0","cmvaee0",1);
  mvaSignalEEMc0 -> SetTitle("EE, 2.0 < pT < 5.0");
  mvaFakeEEMc0   -> SetTitle("EE, 2.0 < pT < 5.0");
  mvaSignalEEMc0 -> GetXaxis()->SetTitle("Id BDT");
  mvaFakeEEMc0   -> GetXaxis()->SetTitle("Id BDT");
  mvaFakeEEMc0   -> DrawNormalized("hist");
  mvaSignalEEMc0 -> DrawNormalized("samehist");
  legB->Draw();
  cmvaee0.SaveAs("outputBDT_matchVsFake_withTnP_EE0.png");  
  //
  TCanvas cmvaee1("cmvaee1","cmvaee1",1);
  mvaSignalEEMc1 -> SetTitle("EE, pT >= 5.0");
  mvaFakeEEMc1   -> SetTitle("EE, pT >= 5.0");
  mvaSignalEEMc1 -> GetXaxis()->SetTitle("Id BDT");
  mvaFakeEEMc1   -> GetXaxis()->SetTitle("Id BDT");
  mvaFakeEEMc1->DrawNormalized("hist");
  mvaSignalEEMc1->DrawNormalized("samehist");
  legB->Draw();
  cmvaee1.SaveAs("outputBDT_matchVsFake_withTnP_EE1.png");  
}
