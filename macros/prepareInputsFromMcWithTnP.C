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

void prepareInputsFromMcWithTnP::Loop(bool applyWeight, bool testLPT, bool studyOverlap=0)
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
  TH1F *probePtSignalMc  = new TH1F("probePtSignalMc",  "probePtSignalMc",  60,  0.,  15.);
  TH1F *probePtFakeMc    = new TH1F("probePtFakeMc",    "probePtFakeMc",    60,  0.,  15.);
  TH1F *probeEtaSignalMc = new TH1F("probeEtaSignalMc", "probeEtaSignalMc", 40, -2.4, 2.4);
  TH1F *probeEtaFakeMc   = new TH1F("probeEtaFakeMc",   "probeEtaFakeMc",   40, -2.4, 2.4);
  TH1F *probeMvaSignalMc = new TH1F("probeMvaSignalMc", "probeMvaSignalMc", 54, -3., 15.);           // 60, -10., 10. for january production
  TH1F *probeMvaFakeMc   = new TH1F("probeMvaFakeMc",   "probeMvaFakeMc",   54, -3., 15.);
  TH1F *probeDxysigSignalMc   = new TH1F("probeDxysigSignalMc",   "probeDxysigSignalMc",   50, -10., 10.);  
  TH1F *probeDxysigFakeMc     = new TH1F("probeDxysigFakeMc",     "probeDxysigFakeMc",     50, -10., 10.);  
  TH1F *probeDztrgSignalMc    = new TH1F("probeDztrgSignalMc",    "probeDztrgSignalMc",    50, -0.5, 0.5); 
  TH1F *probeDztrgFakeMc      = new TH1F("probeDztrgFakeMc",      "probeDztrgFakeMc",      50, -0.5, 0.5);  
  TH1F *probeIso04relSignalMc = new TH1F("probeIso04relSignalMc", "probeIso04relSignalMc", 50,  0.,  10.);  
  TH1F *probeIso04relFakeMc   = new TH1F("probeIso04relFakeMc",   "probeIso04relFakeMc",   50,  0.,  10.);  
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
  TH1F *probePtSignalMcWW  = new TH1F("probePtSignalMcWW",  "probePtSignalMcWW",  60,  0.,  15.);
  TH1F *probeEtaSignalMcWW = new TH1F("probeEtaSignalMcWW", "probeEtaSignalMcWW", 40, -2.4, 2.4);
  TH1F *probeMvaSignalMcWW = new TH1F("probeMvaSignalMcWW", "probeMvaSignalMcWW", 54, -3., 15.); 
  TH1F *probeDxysigSignalMcWW   = new TH1F("probeDxysigSignalMcWW",   "probeDxysigSignalMcWW",   50, -10., 10.);  
  TH1F *probeDztrgSignalMcWW    = new TH1F("probeDztrgSignalMcWW",    "probeDztrgSignalMcWW",    50, -0.5, 0.5); 
  TH1F *probeIso04relSignalMcWW = new TH1F("probeIso04relSignalMcWW", "probeIso04relSignalMcWW", 50,  0.,  10.);  
  probePtSignalMcWW->Sumw2();
  probeEtaSignalMcWW->Sumw2();
  probeMvaSignalMcWW->Sumw2(); 
  probeDxysigSignalMcWW->Sumw2();
  probeDztrgSignalMcWW->Sumw2();
  probeIso04relSignalMcWW->Sumw2();

  // -----------------------------------------------------------------------
  // Many pT/eta bins (no weight): ID output
  TH1F *mvaSignalEBMc0 = new TH1F("mvaSignalEBMc0", "mvaSignalEBMc0", 54, -3., 15.);  // pT: 1.0-1.5
  TH1F *mvaSignalEBMc1 = new TH1F("mvaSignalEBMc1", "mvaSignalEBMc1", 54, -3., 15.);  // 1.5-2.0
  TH1F *mvaSignalEBMc2 = new TH1F("mvaSignalEBMc2", "mvaSignalEBMc2", 54, -3., 15.);  // 2.0-5.0
  TH1F *mvaSignalEBMc3 = new TH1F("mvaSignalEBMc3", "mvaSignalEBMc3", 54, -3., 15.);  // >5
  //
  TH1F *mvaSignalEEMc0 = new TH1F("mvaSignalEEMc0", "mvaSignalEEMc0", 54, -3., 15.);  // pT: 1.0-2.0  
  TH1F *mvaSignalEEMc1 = new TH1F("mvaSignalEEMc1", "mvaSignalEEMc1", 54, -3., 15.);  // pT: 2.0-5.0  
  TH1F *mvaSignalEEMc2 = new TH1F("mvaSignalEEMc2", "mvaSignalEEMc2", 54, -3., 15.);  // >5
  mvaSignalEBMc0->Sumw2();
  mvaSignalEBMc1->Sumw2();
  mvaSignalEBMc2->Sumw2();
  mvaSignalEBMc3->Sumw2();
  mvaSignalEEMc0->Sumw2();
  mvaSignalEEMc1->Sumw2();
  mvaSignalEEMc2->Sumw2();
  //
  TH1F *mvaFakeEBMc0 = new TH1F("mvaFakeEBMc0", "mvaFakeEBMc0", 54, -3., 15.);  // pT: 1.0-1.5
  TH1F *mvaFakeEBMc1 = new TH1F("mvaFakeEBMc1", "mvaFakeEBMc1", 54, -3., 15.);  // 1.5-2.0
  TH1F *mvaFakeEBMc2 = new TH1F("mvaFakeEBMc2", "mvaFakeEBMc2", 54, -3., 15.);  // 2.0-5.0
  TH1F *mvaFakeEBMc3 = new TH1F("mvaFakeEBMc3", "mvaFakeEBMc3", 54, -3., 15.);  // >5
  TH1F *mvaFakeEEMc0 = new TH1F("mvaFakeEEMc0", "mvaFakeEEMc0", 54, -3., 15.);  // pT: 1.0-2.0
  TH1F *mvaFakeEEMc1 = new TH1F("mvaFakeEEMc1", "mvaFakeEEMc1", 54, -3., 15.);  // pT: 2.0-5.0
  TH1F *mvaFakeEEMc2 = new TH1F("mvaFakeEEMc2", "mvaFakeEEMc2", 54, -3., 15.);  // >5
  mvaFakeEBMc0->Sumw2();
  mvaFakeEBMc1->Sumw2();
  mvaFakeEBMc2->Sumw2();
  mvaFakeEBMc3->Sumw2();
  mvaFakeEEMc0->Sumw2();
  mvaFakeEEMc1->Sumw2();
  mvaFakeEEMc2->Sumw2();
  
  // -----------------------------------------------------------------------
  // Many pT/eta bins (with weight): ID output
  TH1F *mvaSignalEBMc0WW = new TH1F("mvaSignalEBMc0WW", "mvaSignalEBMc0WW", 54, -3., 15.);  // pT: 1.0-1.5
  TH1F *mvaSignalEBMc1WW = new TH1F("mvaSignalEBMc1WW", "mvaSignalEBMc1WW", 54, -3., 15.);  // 1.5-2.0
  TH1F *mvaSignalEBMc2WW = new TH1F("mvaSignalEBMc2WW", "mvaSignalEBMc2WW", 54, -3., 15.);  // 2.0-5.0
  TH1F *mvaSignalEBMc3WW = new TH1F("mvaSignalEBMc3WW", "mvaSignalEBMc3WW", 54, -3., 15.);  // >5
  TH1F *mvaSignalEEMc0WW = new TH1F("mvaSignalEEMc0WW", "mvaSignalEEMc0WW", 54, -3., 15.);  // pT: 1.0-2.0 
  TH1F *mvaSignalEEMc1WW = new TH1F("mvaSignalEEMc1WW", "mvaSignalEEMc1WW", 54, -3., 15.);  // pT: 2.0-5.0 
  TH1F *mvaSignalEEMc2WW = new TH1F("mvaSignalEEMc2WW", "mvaSignalEEMc2WW", 54, -3., 15.);  // >5
  mvaSignalEBMc0WW->Sumw2();
  mvaSignalEBMc1WW->Sumw2();
  mvaSignalEBMc2WW->Sumw2();
  mvaSignalEBMc3WW->Sumw2();
  mvaSignalEEMc0WW->Sumw2();
  mvaSignalEEMc1WW->Sumw2();
  mvaSignalEEMc2WW->Sumw2();

  // -----------------------------------------------------------------------
  // Many pT/eta bins (no weight): dxy significance
  TH1F *dxysigSignalEBMc0 = new TH1F("dxysigSignalEBMc0", "dxysigSignalEBMc0", 50, -10., 10.);  // pT: 1.0-1.5
  TH1F *dxysigSignalEBMc1 = new TH1F("dxysigSignalEBMc1", "dxysigSignalEBMc1", 50, -10., 10.);  // 1.5-2.0
  TH1F *dxysigSignalEBMc2 = new TH1F("dxysigSignalEBMc2", "dxysigSignalEBMc2", 50, -10., 10.);  // 2.0-5.0
  TH1F *dxysigSignalEBMc3 = new TH1F("dxysigSignalEBMc3", "dxysigSignalEBMc3", 50, -10., 10.);  // >5
  TH1F *dxysigSignalEEMc0 = new TH1F("dxysigSignalEEMc0", "dxysigSignalEEMc0", 50, -10., 10.);  // pT: 1.0-2.0 
  TH1F *dxysigSignalEEMc1 = new TH1F("dxysigSignalEEMc1", "dxysigSignalEEMc1", 50, -10., 10.);  // pT: 2.0-5.0 
  TH1F *dxysigSignalEEMc2 = new TH1F("dxysigSignalEEMc2", "dxysigSignalEEMc2", 50, -10., 10.);  // >5
  dxysigSignalEBMc0->Sumw2();
  dxysigSignalEBMc1->Sumw2();
  dxysigSignalEBMc2->Sumw2();
  dxysigSignalEBMc3->Sumw2();
  dxysigSignalEEMc0->Sumw2();
  dxysigSignalEEMc1->Sumw2();
  dxysigSignalEEMc2->Sumw2();

  // -----------------------------------------------------------------------
  // Many pT/eta bins (no weight): dzTrg
  TH1F *dztrgSignalEBMc0 = new TH1F("dztrgSignalEBMc0", "dztrgSignalEBMc0", 50, -0.5, 0.5);  // pT: 1.0-1.5
  TH1F *dztrgSignalEBMc1 = new TH1F("dztrgSignalEBMc1", "dztrgSignalEBMc1", 50, -0.5, 0.5);  // 1.5-2.0
  TH1F *dztrgSignalEBMc2 = new TH1F("dztrgSignalEBMc2", "dztrgSignalEBMc2", 50, -0.5, 0.5);  // 2.0-5.0
  TH1F *dztrgSignalEBMc3 = new TH1F("dztrgSignalEBMc3", "dztrgSignalEBMc3", 50, -0.5, 0.5);  // >5
  TH1F *dztrgSignalEEMc0 = new TH1F("dztrgSignalEEMc0", "dztrgSignalEEMc0", 50, -0.5, 0.5);  // pT: 1.0-2.0 
  TH1F *dztrgSignalEEMc1 = new TH1F("dztrgSignalEEMc1", "dztrgSignalEEMc1", 50, -0.5, 0.5);  // pT: 2.0-5.0 
  TH1F *dztrgSignalEEMc2 = new TH1F("dztrgSignalEEMc2", "dztrgSignalEEMc2", 50, -0.5, 0.5);  // >5
  dztrgSignalEBMc0->Sumw2();
  dztrgSignalEBMc1->Sumw2();
  dztrgSignalEBMc2->Sumw2();
  dztrgSignalEBMc3->Sumw2();
  dztrgSignalEEMc0->Sumw2();
  dztrgSignalEEMc1->Sumw2();
  dztrgSignalEEMc2->Sumw2();
  //
  TH1F *dztrgFakeEBMc0 = new TH1F("dztrgFakeEBMc0", "dztrgFakeEBMc0", 50, -0.5, 0.5);  // pT: 0.5-1.5
  TH1F *dztrgFakeEBMc1 = new TH1F("dztrgFakeEBMc1", "dztrgFakeEBMc1", 50, -0.5, 0.5);  // 1.5-2.0
  TH1F *dztrgFakeEBMc2 = new TH1F("dztrgFakeEBMc2", "dztrgFakeEBMc2", 50, -0.5, 0.5);  // 2.0-5.0
  TH1F *dztrgFakeEBMc3 = new TH1F("dztrgFakeEBMc3", "dztrgFakeEBMc3", 50, -0.5, 0.5);  // >5
  TH1F *dztrgFakeEEMc0 = new TH1F("dztrgFakeEEMc0", "dztrgFakeEEMc0", 50, -0.5, 0.5);  // pT: 1.0-2.0
  TH1F *dztrgFakeEEMc1 = new TH1F("dztrgFakeEEMc1", "dztrgFakeEEMc1", 50, -0.5, 0.5);  // pT: 2.0-5.0
  TH1F *dztrgFakeEEMc2 = new TH1F("dztrgFakeEEMc2", "dztrgFakeEEMc2", 50, -0.5, 0.5);  // >5
  dztrgFakeEBMc0->Sumw2();
  dztrgFakeEBMc1->Sumw2();
  dztrgFakeEBMc2->Sumw2();
  dztrgFakeEBMc3->Sumw2();
  dztrgFakeEEMc0->Sumw2();
  dztrgFakeEEMc1->Sumw2();
  dztrgFakeEEMc2->Sumw2();

  // -----------------------------------------------------------------------
  // Many pT/eta bins (no weight): iso04rel
  TH1F *iso04relSignalEBMc0 = new TH1F("iso04relSignalEBMc0", "iso04relSignalEBMc0", 50, 0., 10.);  // pT: 1.0-1.5
  TH1F *iso04relSignalEBMc1 = new TH1F("iso04relSignalEBMc1", "iso04relSignalEBMc1", 50, 0., 10.);  // 1.5-2.0
  TH1F *iso04relSignalEBMc2 = new TH1F("iso04relSignalEBMc2", "iso04relSignalEBMc2", 50, 0., 10.);  // 2.0-5.0
  TH1F *iso04relSignalEBMc3 = new TH1F("iso04relSignalEBMc3", "iso04relSignalEBMc3", 50, 0., 10.);  // >5
  TH1F *iso04relSignalEEMc0 = new TH1F("iso04relSignalEEMc0", "iso04relSignalEEMc0", 50, 0., 10.);  // pT: 1.0-2.0 
  TH1F *iso04relSignalEEMc1 = new TH1F("iso04relSignalEEMc1", "iso04relSignalEEMc1", 50, 0., 10.);  // pT: 2.0-5.0 
  TH1F *iso04relSignalEEMc2 = new TH1F("iso04relSignalEEMc2", "iso04relSignalEEMc2", 50, 0., 10.);  // >5
  iso04relSignalEBMc0->Sumw2();
  iso04relSignalEBMc1->Sumw2();
  iso04relSignalEBMc2->Sumw2();
  iso04relSignalEBMc3->Sumw2();
  iso04relSignalEEMc0->Sumw2();
  iso04relSignalEEMc1->Sumw2();
  iso04relSignalEEMc2->Sumw2();
  //
  TH1F *iso04relFakeEBMc0 = new TH1F("iso04relFakeEBMc0", "iso04relFakeEBMc0", 50, 0., 10.);  // pT: 0.5-1.5
  TH1F *iso04relFakeEBMc1 = new TH1F("iso04relFakeEBMc1", "iso04relFakeEBMc1", 50, 0., 10.);  // 1.5-2.0
  TH1F *iso04relFakeEBMc2 = new TH1F("iso04relFakeEBMc2", "iso04relFakeEBMc2", 50, 0., 10.);  // 2.0-5.0
  TH1F *iso04relFakeEBMc3 = new TH1F("iso04relFakeEBMc3", "iso04relFakeEBMc3", 50, 0., 10.);  // >5
  TH1F *iso04relFakeEEMc0 = new TH1F("iso04relFakeEEMc0", "iso04relFakeEEMc0", 50, 0., 10.);  // pT: 1.0-2.0
  TH1F *iso04relFakeEEMc1 = new TH1F("iso04relFakeEEMc1", "iso04relFakeEEMc1", 50, 0., 10.);  // pT: 2.0-5.0
  TH1F *iso04relFakeEEMc2 = new TH1F("iso04relFakeEEMc2", "iso04relFakeEEMc2", 50, 0., 10.);  // >5
  iso04relFakeEBMc0->Sumw2();
  iso04relFakeEBMc1->Sumw2();
  iso04relFakeEBMc2->Sumw2();
  iso04relFakeEBMc3->Sumw2();
  iso04relFakeEEMc0->Sumw2();
  iso04relFakeEEMc1->Sumw2();
  iso04relFakeEEMc2->Sumw2();

  // --------------------------
  // Full pT range: distributions for data/MC checking overlap between LPT and PF
  TH1F *probePtSignalMc_LptPfOverlap  = new TH1F("probePtSignalMc_LptPfOverlap",  "probePtSignalMc_LptPfOverlap",  60,  0.,  15.);
  TH1F *probeEtaSignalMc_LptPfOverlap = new TH1F("probeEtaSignalMc_LptPfOverlap", "probeEtaSignalMc_LptPfOverlap", 40, -2.4, 2.4);
  TH1F *probeMvaSignalMc_LptPfOverlap = new TH1F("probeMvaSignalMc_LptPfOverlap", "probeMvaSignalMc_LptPfOverlap", 54, -3., 15.); 
  probePtSignalMc_LptPfOverlap->Sumw2();
  probeEtaSignalMc_LptPfOverlap->Sumw2();
  probeMvaSignalMc_LptPfOverlap->Sumw2(); 

  TH1F *probePtSignalMc_LptNotPfOverlap  = new TH1F("probePtSignalMc_LptNotPfOverlap",  "probePtSignalMc_LptNotPfOverlap",  60,  0.,  15.);
  TH1F *probeEtaSignalMc_LptNotPfOverlap = new TH1F("probeEtaSignalMc_LptNotPfOverlap", "probeEtaSignalMc_LptNotPfOverlap", 40, -2.4, 2.4);
  TH1F *probeMvaSignalMc_LptNotPfOverlap = new TH1F("probeMvaSignalMc_LptNotPfOverlap", "probeMvaSignalMc_LptNotPfOverlap", 54, -3., 15.); 
  probePtSignalMc_LptNotPfOverlap->Sumw2();
  probeEtaSignalMc_LptNotPfOverlap->Sumw2();
  probeMvaSignalMc_LptNotPfOverlap->Sumw2(); 

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
    TFile fileWeight("files_marchNoReg/probeLowPt/weightFile_tnpVsNani.root"); 
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

    // Sanity check 
    if (studyOverlap==1 && testLPT!=1) { cout << "should not happen!" << endl; continue; }

    // HLT
    if (hlt_9ip6==0) continue;
    
    // Acceptance
    if (fabs(probeEta)>2.4) continue;
    if (probePt<1.0)        continue;

    // LowPt vs PF selection already applied in formatted-tnp ntuples

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
      probePtSignalMc->Fill(probePt, weight);
      probeEtaSignalMc->Fill(probeEta, weight);
      probePtVsEtaSignalMc->Fill(probeEta,probePt, weight);
      probeDxysigSignalMc   ->Fill(probeDxySig, weight);
      probeDztrgSignalMc    ->Fill(probeDzTrg, weight);
      probeIso04relSignalMc ->Fill(probeIso04Rel, weight);

      if (studyOverlap==1) {
	if (probeIsPFOverlap==1) {
	  probePtSignalMc_LptPfOverlap->Fill(probePt, weight);
	  probeEtaSignalMc_LptPfOverlap->Fill(probeEta, weight);
	  probeMvaSignalMc_LptPfOverlap->Fill(probeMvaId, weight);
	} else {
	  probePtSignalMc_LptNotPfOverlap->Fill(probePt, weight);
	  probeEtaSignalMc_LptNotPfOverlap->Fill(probeEta, weight);
	  probeMvaSignalMc_LptNotPfOverlap->Fill(probeMvaId, weight);
	}
      }
    }
    
    // matching mc-truth in MC: full pT range after weight
    if (probeMatchMcFromJPsi==1) { 
      if (testLPT==1) probeMvaSignalMcWW->Fill(probeMvaId, ptEtaSignalWeight*weight);       
      if (testLPT==0) probeMvaSignalMcWW->Fill(probePfmvaId, ptEtaSignalWeight*weight);       
      probePtSignalMcWW->Fill(probePt, ptEtaSignalWeight*weight);
      probeEtaSignalMcWW->Fill(probeEta, ptEtaSignalWeight*weight);
      probePtVsEtaSignalMcWW->Fill(probeEta,probePt,ptEtaSignalWeight*weight);
      probeDxysigSignalMcWW   ->Fill(probeDxySig, ptEtaSignalWeight*weight);
      probeDztrgSignalMcWW    ->Fill(probeDzTrg, ptEtaSignalWeight*weight);
      probeIso04relSignalMcWW ->Fill(probeIso04Rel, ptEtaSignalWeight*weight);
    }

    // not matching mc-truth in MC: full pT range before weight
    if (probeMatchMcFromJPsi==0) { 
      if (testLPT==1) probeMvaFakeMc->Fill(probeMvaId, weight);
      if (testLPT==0) probeMvaFakeMc->Fill(probePfmvaId, weight);
      probePtFakeMc->Fill(probePt, weight);
      probeEtaFakeMc->Fill(probeEta, weight);
      probeDxysigFakeMc   ->Fill(probeDxySig, weight);
      probeDztrgFakeMc    ->Fill(probeDzTrg, weight);
      probeIso04relFakeMc ->Fill(probeIso04Rel, weight);
    }

    // matching mc-truth in MC: eta/pT bins
    if (probeMatchMcFromJPsi==1) {  // signal

      float theId=-1;
      if (testLPT==1) theId=probeMvaId;
      if (testLPT==0) theId=probePfmvaId;

      if (fabs(probeEta)<1.5) {  // barrel
	if (probePt>=1.0 && probePt<1.5) {
	  mvaSignalEBMc0->Fill(theId, weight);
	  mvaSignalEBMc0WW->Fill(theId,ptEtaSignalWeight*weight);
	  dxysigSignalEBMc0->Fill(probeDxySig, weight);
	  dztrgSignalEBMc0->Fill(probeDzTrg, weight);
	  iso04relSignalEBMc0->Fill(probeIso04Rel, weight);
	}
	if (probePt>=1.5 && probePt<2.0) {
	  mvaSignalEBMc1->Fill(theId, weight);
	  mvaSignalEBMc1WW->Fill(theId,ptEtaSignalWeight*weight);
	  dxysigSignalEBMc1->Fill(probeDxySig, weight);
	  dztrgSignalEBMc1->Fill(probeDzTrg, weight);
	  iso04relSignalEBMc1->Fill(probeIso04Rel, weight);
	}
	if (probePt>=2.0 && probePt<5.0) {
	  mvaSignalEBMc2->Fill(theId, weight);
	  mvaSignalEBMc2WW->Fill(theId,ptEtaSignalWeight*weight);
	  dxysigSignalEBMc2->Fill(probeDxySig, weight);
	  dztrgSignalEBMc2->Fill(probeDzTrg, weight);
	  iso04relSignalEBMc2->Fill(probeIso04Rel, weight);
	}
	if (probePt>=5.0) {
	  mvaSignalEBMc3->Fill(theId, weight);
	  mvaSignalEBMc3WW->Fill(theId,ptEtaSignalWeight*weight);
	  dxysigSignalEBMc3->Fill(probeDxySig, weight);
	  dztrgSignalEBMc3->Fill(probeDzTrg, weight);
	  iso04relSignalEBMc3->Fill(probeIso04Rel, weight);
	}

      } else {  // endcap

	if (probePt>=1.0 && probePt<2.0) {
	  mvaSignalEEMc0->Fill(theId, weight);
	  mvaSignalEEMc0WW->Fill(theId,ptEtaSignalWeight*weight);
	  dxysigSignalEEMc0->Fill(probeDxySig, weight);
	  dztrgSignalEEMc0->Fill(probeDzTrg, weight);
	  iso04relSignalEEMc0->Fill(probeIso04Rel, weight);
	}
	if (probePt>=2.0 && probePt<5.0) {
	  mvaSignalEEMc1->Fill(theId, weight);
	  mvaSignalEEMc1WW->Fill(theId,ptEtaSignalWeight*weight);
	  dxysigSignalEEMc1->Fill(probeDxySig, weight);
	  dztrgSignalEEMc1->Fill(probeDzTrg, weight);
	  iso04relSignalEEMc1->Fill(probeIso04Rel, weight);
	}
	if (probePt>=5.0) {
	  mvaSignalEEMc2->Fill(theId, weight);
	  mvaSignalEEMc2WW->Fill(theId,ptEtaSignalWeight*weight);
	  dxysigSignalEEMc2->Fill(probeDxySig, weight);
	  dztrgSignalEEMc2->Fill(probeDzTrg, weight);
	  iso04relSignalEEMc2->Fill(probeIso04Rel, weight);
	}
      }
    }

    // not matching mc-truth in MC: eta/pT bins
    if (probeMatchMcFromJPsi==0) { 

      float theId=-1;
      if (testLPT==1) theId=probeMvaId;
      if (testLPT==0) theId=probePfmvaId;

      if (fabs(probeEta)<1.5) {  // barrel
	if (probePt>=1.0 && probePt<1.5) mvaFakeEBMc0->Fill(theId, weight);
	if (probePt>=1.5 && probePt<2.0) mvaFakeEBMc1->Fill(theId, weight);
	if (probePt>=2.0 && probePt<5.0) mvaFakeEBMc2->Fill(theId, weight);
	if (probePt>=5.0) mvaFakeEBMc3->Fill(theId, weight);

      } else {  // endcap
	if (probePt>=1.0 && probePt<2.0) mvaFakeEEMc0->Fill(theId, weight);
	if (probePt>=2.0 && probePt<5.0) mvaFakeEEMc1->Fill(theId, weight);
	if (probePt>=5.0) mvaFakeEEMc2->Fill(theId, weight);
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
  mvaSignalEBMc2 -> SetLineWidth(2); 
  mvaSignalEBMc3 -> SetLineWidth(2); 
  mvaSignalEBMc0 -> SetLineColor(6); 
  mvaSignalEBMc1 -> SetLineColor(6); 
  mvaSignalEBMc2 -> SetLineColor(6); 
  mvaSignalEBMc3 -> SetLineColor(6); 
  //
  mvaSignalEBMc0WW -> SetLineWidth(2); 
  mvaSignalEBMc1WW -> SetLineWidth(2); 
  mvaSignalEBMc2WW -> SetLineWidth(2); 
  mvaSignalEBMc3WW -> SetLineWidth(2); 
  mvaSignalEBMc0WW -> SetLineColor(6); 
  mvaSignalEBMc1WW -> SetLineColor(6); 
  mvaSignalEBMc2WW -> SetLineColor(6); 
  mvaSignalEBMc3WW -> SetLineColor(6); 
  //
  mvaSignalEEMc0 -> SetLineWidth(2); 
  mvaSignalEEMc1 -> SetLineWidth(2); 
  mvaSignalEEMc2 -> SetLineWidth(2); 
  mvaSignalEEMc0 -> SetLineColor(6); 
  mvaSignalEEMc1 -> SetLineColor(6); 
  mvaSignalEEMc2 -> SetLineColor(6); 
  //
  mvaSignalEEMc0WW -> SetLineWidth(2); 
  mvaSignalEEMc1WW -> SetLineWidth(2); 
  mvaSignalEEMc2WW -> SetLineWidth(2); 
  mvaSignalEEMc0WW -> SetLineColor(6); 
  mvaSignalEEMc1WW -> SetLineColor(6); 
  mvaSignalEEMc2WW -> SetLineColor(6); 
  //
  mvaFakeEBMc0 -> SetLineWidth(2); 
  mvaFakeEBMc1 -> SetLineWidth(2); 
  mvaFakeEBMc2 -> SetLineWidth(2); 
  mvaFakeEBMc3 -> SetLineWidth(2); 
  mvaFakeEBMc0 -> SetLineColor(4); 
  mvaFakeEBMc1 -> SetLineColor(4); 
  mvaFakeEBMc2 -> SetLineColor(4); 
  mvaFakeEBMc3 -> SetLineColor(4); 
  //
  mvaFakeEEMc0 -> SetLineWidth(2); 
  mvaFakeEEMc1 -> SetLineWidth(2); 
  mvaFakeEEMc2 -> SetLineWidth(2); 
  mvaFakeEEMc0 -> SetLineColor(4); 
  mvaFakeEEMc1 -> SetLineColor(4); 
  mvaFakeEEMc2 -> SetLineColor(4); 


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
  if (studyOverlap==1) {
    probePtSignalMc_LptNotPfOverlap->Write();
    probeEtaSignalMc_LptNotPfOverlap->Write();
    probeMvaSignalMc_LptNotPfOverlap->Write();
    probePtSignalMc_LptPfOverlap->Write();
    probeEtaSignalMc_LptPfOverlap->Write();
    probeMvaSignalMc_LptPfOverlap->Write();
  }
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
  mvaSignalEBMc2->Write();
  mvaSignalEBMc3->Write();
  mvaSignalEEMc0->Write();
  mvaSignalEEMc1->Write();
  mvaSignalEEMc2->Write();
  //
  mvaSignalEBMc0WW->Write();
  mvaSignalEBMc1WW->Write();
  mvaSignalEBMc2WW->Write();
  mvaSignalEBMc3WW->Write();
  mvaSignalEEMc0WW->Write();
  mvaSignalEEMc1WW->Write();
  mvaSignalEEMc2WW->Write();
  //
  mvaFakeEBMc0->Write();
  mvaFakeEBMc1->Write();
  mvaFakeEBMc2->Write();
  mvaFakeEBMc3->Write();
  mvaFakeEEMc0->Write();
  mvaFakeEEMc1->Write();
  mvaFakeEEMc2->Write();
  //
  dxysigSignalEBMc0->Write();
  dxysigSignalEBMc1->Write();
  dxysigSignalEBMc2->Write();
  dxysigSignalEBMc3->Write();
  dxysigSignalEEMc0->Write();
  dxysigSignalEEMc1->Write();
  dxysigSignalEEMc2->Write();
  //
  dztrgSignalEBMc0->Write();
  dztrgSignalEBMc1->Write();
  dztrgSignalEBMc2->Write();
  dztrgSignalEBMc3->Write();
  dztrgSignalEEMc0->Write();
  dztrgSignalEEMc1->Write();
  dztrgSignalEEMc2->Write();
  //
  iso04relSignalEBMc0->Write();
  iso04relSignalEBMc1->Write();
  iso04relSignalEBMc2->Write();
  iso04relSignalEBMc3->Write();
  iso04relSignalEEMc0->Write();
  iso04relSignalEEMc1->Write();
  iso04relSignalEEMc2->Write();
  //

  myFile.Close();


  // -----------------------------------------------------------------------
  // Rebin
  mvaSignalEBMc0->Rebin(4);
  mvaSignalEBMc1->Rebin(4);
  mvaSignalEBMc2->Rebin(4);
  mvaSignalEBMc3->Rebin(4);
  mvaSignalEEMc0->Rebin(4);
  mvaSignalEEMc1->Rebin(4);
  mvaSignalEEMc2->Rebin(4);
  mvaFakeEBMc0->Rebin(4);
  mvaFakeEBMc1->Rebin(4);
  mvaFakeEBMc2->Rebin(4);
  mvaFakeEBMc3->Rebin(4);
  mvaFakeEEMc0->Rebin(4);
  mvaFakeEEMc1->Rebin(4);
  mvaFakeEEMc2->Rebin(4);


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
  mvaSignalEBMc0 -> SetTitle("EB, 0.5 < pT < 1.5");
  mvaFakeEBMc0   -> SetTitle("EB, 0.5 < pT < 1.5");
  mvaSignalEBMc0 -> GetXaxis()->SetTitle("Id BDT");
  mvaFakeEBMc0   -> GetXaxis()->SetTitle("Id BDT");
  mvaFakeEBMc0->DrawNormalized("hist");
  mvaSignalEBMc0->DrawNormalized("samehist");
  legB->Draw();
  cmvaeb0.SaveAs("outputBDT_matchVsFake_withTnP_EB0.png");  
  //
  TCanvas cmvaeb1("cmvaeb1","cmvaeb1",1);
  mvaSignalEBMc1 -> SetTitle("EB, 1.5 < pT < 2.0");
  mvaFakeEBMc1   -> SetTitle("EB, 1.5 < pT < 2.0");
  mvaSignalEBMc1 -> GetXaxis()->SetTitle("Id BDT");
  mvaFakeEBMc1   -> GetXaxis()->SetTitle("Id BDT");
  mvaFakeEBMc1->DrawNormalized("hist");
  mvaSignalEBMc1->DrawNormalized("samehist");
  legB->Draw();
  cmvaeb1.SaveAs("outputBDT_matchVsFake_withTnP_EB1.png");  
  //
  TCanvas cmvaeb2("cmvaeb2","cmvaeb2",1);
  mvaSignalEBMc2 -> SetTitle("EB, 2.0 < pT < 5.0");
  mvaFakeEBMc2   -> SetTitle("EB, 2.0 < pT < 5.0");
  mvaSignalEBMc2 -> GetXaxis()->SetTitle("Id BDT");
  mvaFakeEBMc2   -> GetXaxis()->SetTitle("Id BDT");
  mvaFakeEBMc2->DrawNormalized("hist");
  mvaSignalEBMc2->DrawNormalized("samehist");
  legB->Draw();
  cmvaeb2.SaveAs("outputBDT_matchVsFake_withTnP_EB2.png");  
  //
  TCanvas cmvaeb3("cmvaeb3","cmvaeb3",1);
  mvaSignalEBMc3 -> SetTitle("EB, pT >= 5");
  mvaFakeEBMc3   -> SetTitle("EB, pT >= 5");
  mvaSignalEBMc3 -> GetXaxis()->SetTitle("Id BDT");
  mvaFakeEBMc3   -> GetXaxis()->SetTitle("Id BDT");
  mvaSignalEBMc3->DrawNormalized("hist");
  mvaFakeEBMc3->DrawNormalized("samehist");
  legB->Draw();
  cmvaeb3.SaveAs("outputBDT_matchVsFake_withTnP_EB3.png");  

  TCanvas cmvaee0("cmvaee0","cmvaee0",1);
  mvaSignalEEMc0 -> SetTitle("EE, 1.0 < pT < 2.0");
  mvaFakeEEMc0   -> SetTitle("EE, 1.0 < pT < 2.0");
  mvaSignalEEMc0 -> GetXaxis()->SetTitle("Id BDT");
  mvaFakeEEMc0   -> GetXaxis()->SetTitle("Id BDT");
  mvaFakeEEMc0->DrawNormalized("hist");
  mvaSignalEEMc0->DrawNormalized("samehist");
  legB->Draw();
  cmvaee0.SaveAs("outputBDT_matchVsFake_withTnP_EE0.png");  
  //
  TCanvas cmvaee1("cmvaee1","cmvaee1",1);
  mvaSignalEEMc1 -> SetTitle("EE, 2.0 < pT < 5.0");
  mvaFakeEEMc1   -> SetTitle("EE, 2.0 < pT < 5.0");
  mvaSignalEEMc1 -> GetXaxis()->SetTitle("Id BDT");
  mvaFakeEEMc1   -> GetXaxis()->SetTitle("Id BDT");
  mvaFakeEEMc1->DrawNormalized("hist");
  mvaSignalEEMc1->DrawNormalized("samehist");
  legB->Draw();
  cmvaee1.SaveAs("outputBDT_matchVsFake_withTnP_EE1.png");  
  //
  TCanvas cmvaee2("cmvaee2","cmvaee2",1);
  mvaSignalEEMc2 -> SetTitle("EE, pT >= 5.0");
  mvaFakeEEMc2   -> SetTitle("EE, pT >= 5.0");
  mvaSignalEEMc2 -> GetXaxis()->SetTitle("Id BDT");
  mvaFakeEEMc2   -> GetXaxis()->SetTitle("Id BDT");
  mvaSignalEEMc2->DrawNormalized("hist");
  mvaFakeEEMc2->DrawNormalized("samehist");
  legB->Draw();
  cmvaee2.SaveAs("outputBDT_matchVsFake_withTnP_EE2.png");  
}
