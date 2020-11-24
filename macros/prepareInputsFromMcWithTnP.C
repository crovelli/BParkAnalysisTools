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
  // No weight
  TH1F *probePtSignalMc  = new TH1F("probePtSignalMc",  "probePtSignalMc",  60,  0.,  15.);
  TH1F *probePtFakeMc    = new TH1F("probePtFakeMc",    "probePtFakeMc",    60,  0.,  15.);
  TH1F *probeEtaSignalMc = new TH1F("probeEtaSignalMc", "probeEtaSignalMc", 40, -2.4, 2.4);
  TH1F *probeEtaFakeMc   = new TH1F("probeEtaFakeMc",   "probeEtaFakeMc",   40, -2.4, 2.4);
  TH1F *probeMvaSignalMc = new TH1F("probeMvaSignalMc", "probeMvaSignalMc", 60, -10., 10.); 
  TH1F *probeMvaFakeMc   = new TH1F("probeMvaFakeMc",   "probeMvaFakeMc",   60, -10., 10.);
  TH1F *numVtxMc         = new TH1F("numVtxMc",         "numVtxMc",         30,   0., 30.);  
  probePtSignalMc->Sumw2();
  probePtFakeMc->Sumw2();
  probeEtaSignalMc->Sumw2();
  probeEtaFakeMc->Sumw2();
  probeMvaSignalMc->Sumw2(); 
  probeMvaFakeMc->Sumw2();  
  numVtxMc->Sumw2();     
  //
  // To check eta in pT bins, wo weights
  TH1F *probeEtaSignalMc0 = new TH1F("probeEtaSignalMc0", "probeEtaSignalMc0", 40, -2.4, 2.4);
  TH1F *probeEtaFakeMc0   = new TH1F("probeEtaFakeMc0",   "probeEtaFakeMc0",   40, -2.4, 2.4);
  TH1F *probeEtaSignalMc1 = new TH1F("probeEtaSignalMc1", "probeEtaSignalMc1", 40, -2.4, 2.4);
  TH1F *probeEtaFakeMc1   = new TH1F("probeEtaFakeMc1",   "probeEtaFakeMc1",   40, -2.4, 2.4);
  TH1F *probeEtaSignalMc2 = new TH1F("probeEtaSignalMc2", "probeEtaSignalMc2", 40, -2.4, 2.4);
  TH1F *probeEtaFakeMc2   = new TH1F("probeEtaFakeMc2",   "probeEtaFakeMc2",   40, -2.4, 2.4);
  TH1F *probeEtaSignalMc3 = new TH1F("probeEtaSignalMc3", "probeEtaSignalMc3", 40, -2.4, 2.4);
  TH1F *probeEtaFakeMc3   = new TH1F("probeEtaFakeMc3",   "probeEtaFakeMc3",   40, -2.4, 2.4);
  TH1F *probeEtaSignalMc4 = new TH1F("probeEtaSignalMc4", "probeEtaSignalMc4", 40, -2.4, 2.4);
  TH1F *probeEtaFakeMc4   = new TH1F("probeEtaFakeMc4",   "probeEtaFakeMc4",   40, -2.4, 2.4);
  probeEtaSignalMc0->Sumw2();
  probeEtaFakeMc0->Sumw2();
  probeEtaSignalMc1->Sumw2();
  probeEtaFakeMc1->Sumw2();
  probeEtaSignalMc2->Sumw2();
  probeEtaFakeMc2->Sumw2();
  probeEtaSignalMc3->Sumw2();
  probeEtaFakeMc3->Sumw2();
  probeEtaSignalMc4->Sumw2();
  probeEtaFakeMc4->Sumw2();
  //
  // With weight
  TH1F *probePtSignalMcWW  = new TH1F("probePtSignalMcWW",  "probePtSignalMcWW",  60,  0.,  15.);
  TH1F *probeEtaSignalMcWW = new TH1F("probeEtaSignalMcWW", "probeEtaSignalMcWW", 40, -2.4, 2.4);
  TH1F *probeMvaSignalMcWW = new TH1F("probeMvaSignalMcWW", "probeMvaSignalMcWW", 60, -10., 10.); 
  probePtSignalMcWW->Sumw2();
  probeEtaSignalMcWW->Sumw2();
  probeMvaSignalMcWW->Sumw2(); 
  //
  // To check eta in pT bins, with weight
  TH1F *probeEtaSignalMcWW0 = new TH1F("probeEtaSignalMcWW0", "probeEtaSignalMcWW0", 40, -2.4, 2.4);
  TH1F *probeEtaSignalMcWW1 = new TH1F("probeEtaSignalMcWW1", "probeEtaSignalMcWW1", 40, -2.4, 2.4);
  TH1F *probeEtaSignalMcWW2 = new TH1F("probeEtaSignalMcWW2", "probeEtaSignalMcWW2", 40, -2.4, 2.4);
  TH1F *probeEtaSignalMcWW3 = new TH1F("probeEtaSignalMcWW3", "probeEtaSignalMcWW3", 40, -2.4, 2.4);
  TH1F *probeEtaSignalMcWW4 = new TH1F("probeEtaSignalMcWW4", "probeEtaSignalMcWW4", 40, -2.4, 2.4);
  probeEtaSignalMcWW0->Sumw2();
  probeEtaSignalMcWW1->Sumw2();
  probeEtaSignalMcWW2->Sumw2();
  probeEtaSignalMcWW3->Sumw2();
  probeEtaSignalMcWW4->Sumw2();


  // -----------------------------------------------------------------------
  // Full pT range: distributions for data/MC (no weight)
  TH1F *probeFBremSignalMc  = new TH1F("probeFBremSignalMc",  "probeFBremSignalMc",  50,   0.,  1.);
  TH1F *probeFBremFakeMc    = new TH1F("probeFBremFakeMc",    "probeFBremFakeMc",    50,   0.,  1.);
  TH1F *probeDxySigSignalMc = new TH1F("probeDxySigSignalMc", "probeDxySigSignalMc", 50, -50., 50.);
  TH1F *probeDxySigFakeMc   = new TH1F("probeDxySigFakeMc",   "probeDxySigFakeMc",   50, -50., 50.);
  TH1F *probeDzSigSignalMc  = new TH1F("probeDzSigSignalMc",  "probeDzSigSignalMc",  50, -50., 50.);
  TH1F *probeDzSigFakeMc    = new TH1F("probeDzSigFakeMc",    "probeDzSigFakeMc",    50, -50., 50.);
  probeFBremSignalMc->Sumw2();
  probeFBremFakeMc->Sumw2();
  probeDxySigSignalMc->Sumw2();
  probeDxySigFakeMc->Sumw2();
  probeDzSigSignalMc->Sumw2();
  probeDzSigFakeMc->Sumw2();

  // -----------------------------------------------------------------------
  // Many pT/eta bins (no weight): ID output
  TH1F *mvaSignalEBMc0 = new TH1F("mvaSignalEBMc0", "mvaSignalEBMc0", 60, -10., 10.);  // pT: 0.5-1.0
  TH1F *mvaSignalEBMc1 = new TH1F("mvaSignalEBMc1", "mvaSignalEBMc1", 60, -10., 10.);  // 1.0-1.5
  TH1F *mvaSignalEBMc2 = new TH1F("mvaSignalEBMc2", "mvaSignalEBMc2", 60, -10., 10.);  // 1.5-2.0
  TH1F *mvaSignalEBMc3 = new TH1F("mvaSignalEBMc3", "mvaSignalEBMc3", 60, -10., 10.);  // 2.0-5.0
  TH1F *mvaSignalEBMc4 = new TH1F("mvaSignalEBMc4", "mvaSignalEBMc4", 60, -10., 10.);  // >5
  TH1F *mvaSignalEEMc0 = new TH1F("mvaSignalEEMc0", "mvaSignalEEMc0", 60, -10., 10.);  // pT: 0.5-1.0
  TH1F *mvaSignalEEMc1 = new TH1F("mvaSignalEEMc1", "mvaSignalEEMc1", 60, -10., 10.);  // 1.0-1.5
  TH1F *mvaSignalEEMc2 = new TH1F("mvaSignalEEMc2", "mvaSignalEEMc2", 60, -10., 10.);  // 1.5-2.0
  TH1F *mvaSignalEEMc3 = new TH1F("mvaSignalEEMc3", "mvaSignalEEMc3", 60, -10., 10.);  // 2.0-5.0
  TH1F *mvaSignalEEMc4 = new TH1F("mvaSignalEEMc4", "mvaSignalEEMc4", 60, -10., 10.);  // >5
  mvaSignalEBMc0->Sumw2();
  mvaSignalEBMc1->Sumw2();
  mvaSignalEBMc2->Sumw2();
  mvaSignalEBMc3->Sumw2();
  mvaSignalEBMc4->Sumw2();
  mvaSignalEEMc0->Sumw2();
  mvaSignalEEMc1->Sumw2();
  mvaSignalEEMc2->Sumw2();
  mvaSignalEEMc3->Sumw2();
  mvaSignalEEMc4->Sumw2();
  //
  TH1F *mvaFakeEBMc0 = new TH1F("mvaFakeEBMc0", "mvaFakeEBMc0", 60, -10., 10.);  // pT: 0.5-1.0
  TH1F *mvaFakeEBMc1 = new TH1F("mvaFakeEBMc1", "mvaFakeEBMc1", 60, -10., 10.);  // 1.0-1.5
  TH1F *mvaFakeEBMc2 = new TH1F("mvaFakeEBMc2", "mvaFakeEBMc2", 60, -10., 10.);  // 1.5-2.0
  TH1F *mvaFakeEBMc3 = new TH1F("mvaFakeEBMc3", "mvaFakeEBMc3", 60, -10., 10.);  // 2.0-5.0
  TH1F *mvaFakeEBMc4 = new TH1F("mvaFakeEBMc4", "mvaFakeEBMc4", 60, -10., 10.);  // >5
  TH1F *mvaFakeEEMc0 = new TH1F("mvaFakeEEMc0", "mvaFakeEEMc0", 60, -10., 10.);  // pT: 0.5-1.0
  TH1F *mvaFakeEEMc1 = new TH1F("mvaFakeEEMc1", "mvaFakeEEMc1", 60, -10., 10.);  // 1.0-1.5
  TH1F *mvaFakeEEMc2 = new TH1F("mvaFakeEEMc2", "mvaFakeEEMc2", 60, -10., 10.);  // 1.5-2.0
  TH1F *mvaFakeEEMc3 = new TH1F("mvaFakeEEMc3", "mvaFakeEEMc3", 60, -10., 10.);  // 2.0-5.0
  TH1F *mvaFakeEEMc4 = new TH1F("mvaFakeEEMc4", "mvaFakeEEMc4", 60, -10., 10.);  // >5
  mvaFakeEBMc0->Sumw2();
  mvaFakeEBMc1->Sumw2();
  mvaFakeEBMc2->Sumw2();
  mvaFakeEBMc3->Sumw2();
  mvaFakeEBMc4->Sumw2();
  mvaFakeEEMc0->Sumw2();
  mvaFakeEEMc1->Sumw2();
  mvaFakeEEMc2->Sumw2();
  mvaFakeEEMc3->Sumw2();
  mvaFakeEEMc4->Sumw2();

  // -----------------------------------------------------------------------
  // Many pT/eta bins (with weight): ID output
  TH1F *mvaSignalEBMc0WW = new TH1F("mvaSignalEBMc0WW", "mvaSignalEBMc0WW", 60, -10., 10.);  // pT: 0.5-1.0
  TH1F *mvaSignalEBMc1WW = new TH1F("mvaSignalEBMc1WW", "mvaSignalEBMc1WW", 60, -10., 10.);  // 1.0-1.5
  TH1F *mvaSignalEBMc2WW = new TH1F("mvaSignalEBMc2WW", "mvaSignalEBMc2WW", 60, -10., 10.);  // 1.5-2.0
  TH1F *mvaSignalEBMc3WW = new TH1F("mvaSignalEBMc3WW", "mvaSignalEBMc3WW", 60, -10., 10.);  // 2.0-5.0
  TH1F *mvaSignalEBMc4WW = new TH1F("mvaSignalEBMc4WW", "mvaSignalEBMc4WW", 60, -10., 10.);  // >5
  TH1F *mvaSignalEEMc0WW = new TH1F("mvaSignalEEMc0WW", "mvaSignalEEMc0WW", 60, -10., 10.);  // pT: 0.5-1.0
  TH1F *mvaSignalEEMc1WW = new TH1F("mvaSignalEEMc1WW", "mvaSignalEEMc1WW", 60, -10., 10.);  // 1.0-1.5
  TH1F *mvaSignalEEMc2WW = new TH1F("mvaSignalEEMc2WW", "mvaSignalEEMc2WW", 60, -10., 10.);  // 1.5-2.0
  TH1F *mvaSignalEEMc3WW = new TH1F("mvaSignalEEMc3WW", "mvaSignalEEMc3WW", 60, -10., 10.);  // 2.0-5.0
  TH1F *mvaSignalEEMc4WW = new TH1F("mvaSignalEEMc4WW", "mvaSignalEEMc4WW", 60, -10., 10.);  // >5
  mvaSignalEBMc0WW->Sumw2();
  mvaSignalEBMc1WW->Sumw2();
  mvaSignalEBMc2WW->Sumw2();
  mvaSignalEBMc3WW->Sumw2();
  mvaSignalEBMc4WW->Sumw2();
  mvaSignalEEMc0WW->Sumw2();
  mvaSignalEEMc1WW->Sumw2();
  mvaSignalEEMc2WW->Sumw2();
  mvaSignalEEMc3WW->Sumw2();
  mvaSignalEEMc4WW->Sumw2();

  // -----------------------------------------------------------------------
  // Many pT/eta bins (no weight): F-brem
  TH1F *fbremSignalEBMc0 = new TH1F("fbremSignalEBMc0", "fbremSignalEBMc0", 50, 0., 1.);  // pT: 0.5-1.0
  TH1F *fbremSignalEBMc1 = new TH1F("fbremSignalEBMc1", "fbremSignalEBMc1", 50, 0., 1.);  // 1.0-1.5
  TH1F *fbremSignalEBMc2 = new TH1F("fbremSignalEBMc2", "fbremSignalEBMc2", 50, 0., 1.);  // 1.5-2.0
  TH1F *fbremSignalEBMc3 = new TH1F("fbremSignalEBMc3", "fbremSignalEBMc3", 50, 0., 1.);  // 2.0-5.0
  TH1F *fbremSignalEBMc4 = new TH1F("fbremSignalEBMc4", "fbremSignalEBMc4", 50, 0., 1.);  // >5
  TH1F *fbremSignalEEMc0 = new TH1F("fbremSignalEEMc0", "fbremSignalEEMc0", 50, 0., 1.);  // pT: 0.5-1.0
  TH1F *fbremSignalEEMc1 = new TH1F("fbremSignalEEMc1", "fbremSignalEEMc1", 50, 0., 1.);  // 1.0-1.5
  TH1F *fbremSignalEEMc2 = new TH1F("fbremSignalEEMc2", "fbremSignalEEMc2", 50, 0., 1.);  // 1.5-2.0
  TH1F *fbremSignalEEMc3 = new TH1F("fbremSignalEEMc3", "fbremSignalEEMc3", 50, 0., 1.);  // 2.0-5.0
  TH1F *fbremSignalEEMc4 = new TH1F("fbremSignalEEMc4", "fbremSignalEEMc4", 50, 0., 1.);  // >5
  fbremSignalEBMc0->Sumw2();
  fbremSignalEBMc1->Sumw2();
  fbremSignalEBMc2->Sumw2();
  fbremSignalEBMc3->Sumw2();
  fbremSignalEBMc4->Sumw2();
  fbremSignalEEMc0->Sumw2();
  fbremSignalEEMc1->Sumw2();
  fbremSignalEEMc2->Sumw2();
  fbremSignalEEMc3->Sumw2();
  fbremSignalEEMc4->Sumw2();
  //
  TH1F *fbremFakeEBMc0 = new TH1F("fbremFakeEBMc0", "fbremFakeEBMc0", 50, 0., 1.);  // pT: 0.5-1.0
  TH1F *fbremFakeEBMc1 = new TH1F("fbremFakeEBMc1", "fbremFakeEBMc1", 50, 0., 1.);  // 1.0-1.5
  TH1F *fbremFakeEBMc2 = new TH1F("fbremFakeEBMc2", "fbremFakeEBMc2", 50, 0., 1.);  // 1.5-2.0
  TH1F *fbremFakeEBMc3 = new TH1F("fbremFakeEBMc3", "fbremFakeEBMc3", 50, 0., 1.);  // 2.0-5.0
  TH1F *fbremFakeEBMc4 = new TH1F("fbremFakeEBMc4", "fbremFakeEBMc4", 50, 0., 1.);  // >5
  TH1F *fbremFakeEEMc0 = new TH1F("fbremFakeEEMc0", "fbremFakeEEMc0", 50, 0., 1.);  // pT: 0.5-1.0
  TH1F *fbremFakeEEMc1 = new TH1F("fbremFakeEEMc1", "fbremFakeEEMc1", 50, 0., 1.);  // 1.0-1.5
  TH1F *fbremFakeEEMc2 = new TH1F("fbremFakeEEMc2", "fbremFakeEEMc2", 50, 0., 1.);  // 1.5-2.0
  TH1F *fbremFakeEEMc3 = new TH1F("fbremFakeEEMc3", "fbremFakeEEMc3", 50, 0., 1.);  // 2.0-5.0
  TH1F *fbremFakeEEMc4 = new TH1F("fbremFakeEEMc4", "fbremFakeEEMc4", 50, 0., 1.);  // >5
  fbremFakeEBMc0->Sumw2();
  fbremFakeEBMc1->Sumw2();
  fbremFakeEBMc2->Sumw2();
  fbremFakeEBMc3->Sumw2();
  fbremFakeEBMc4->Sumw2();
  fbremFakeEEMc0->Sumw2();
  fbremFakeEEMc1->Sumw2();
  fbremFakeEEMc2->Sumw2();
  fbremFakeEEMc3->Sumw2();
  fbremFakeEEMc4->Sumw2();

  // -----------------------------------------------------------------------
  // Many pT/eta bins (no weight): significance dXY
  TH1F *dxysigSignalEBMc0 = new TH1F("dxysigSignalEBMc0", "dxysigSignalEBMc0", 50, -50., 50.);  // pT: 0.5-1.0
  TH1F *dxysigSignalEBMc1 = new TH1F("dxysigSignalEBMc1", "dxysigSignalEBMc1", 50, -50., 50.);  // 1.0-1.5
  TH1F *dxysigSignalEBMc2 = new TH1F("dxysigSignalEBMc2", "dxysigSignalEBMc2", 50, -50., 50.);  // 1.5-2.0
  TH1F *dxysigSignalEBMc3 = new TH1F("dxysigSignalEBMc3", "dxysigSignalEBMc3", 50, -50., 50.);  // 2.0-5.0
  TH1F *dxysigSignalEBMc4 = new TH1F("dxysigSignalEBMc4", "dxysigSignalEBMc4", 50, -50., 50.);  // >5
  TH1F *dxysigSignalEEMc0 = new TH1F("dxysigSignalEEMc0", "dxysigSignalEEMc0", 50, -50., 50.);  // pT: 0.5-1.0
  TH1F *dxysigSignalEEMc1 = new TH1F("dxysigSignalEEMc1", "dxysigSignalEEMc1", 50, -50., 50.);  // 1.0-1.5
  TH1F *dxysigSignalEEMc2 = new TH1F("dxysigSignalEEMc2", "dxysigSignalEEMc2", 50, -50., 50.);  // 1.5-2.0
  TH1F *dxysigSignalEEMc3 = new TH1F("dxysigSignalEEMc3", "dxysigSignalEEMc3", 50, -50., 50.);  // 2.0-5.0
  TH1F *dxysigSignalEEMc4 = new TH1F("dxysigSignalEEMc4", "dxysigSignalEEMc4", 50, -50., 50.);  // >5
  dxysigSignalEBMc0->Sumw2();
  dxysigSignalEBMc1->Sumw2();
  dxysigSignalEBMc2->Sumw2();
  dxysigSignalEBMc3->Sumw2();
  dxysigSignalEBMc4->Sumw2();
  dxysigSignalEEMc0->Sumw2();
  dxysigSignalEEMc1->Sumw2();
  dxysigSignalEEMc2->Sumw2();
  dxysigSignalEEMc3->Sumw2();
  dxysigSignalEEMc4->Sumw2();
  //
  TH1F *dxysigFakeEBMc0 = new TH1F("dxysigFakeEBMc0", "dxysigFakeEBMc0", 50, -50., 50.);  // pT: 0.5-1.0
  TH1F *dxysigFakeEBMc1 = new TH1F("dxysigFakeEBMc1", "dxysigFakeEBMc1", 50, -50., 50.);  // 1.0-1.5
  TH1F *dxysigFakeEBMc2 = new TH1F("dxysigFakeEBMc2", "dxysigFakeEBMc2", 50, -50., 50.);  // 1.5-2.0
  TH1F *dxysigFakeEBMc3 = new TH1F("dxysigFakeEBMc3", "dxysigFakeEBMc3", 50, -50., 50.);  // 2.0-5.0
  TH1F *dxysigFakeEBMc4 = new TH1F("dxysigFakeEBMc4", "dxysigFakeEBMc4", 50, -50., 50.);  // >5
  TH1F *dxysigFakeEEMc0 = new TH1F("dxysigFakeEEMc0", "dxysigFakeEEMc0", 50, -50., 50.);  // pT: 0.5-1.0
  TH1F *dxysigFakeEEMc1 = new TH1F("dxysigFakeEEMc1", "dxysigFakeEEMc1", 50, -50., 50.);  // 1.0-1.5
  TH1F *dxysigFakeEEMc2 = new TH1F("dxysigFakeEEMc2", "dxysigFakeEEMc2", 50, -50., 50.);  // 1.5-2.0
  TH1F *dxysigFakeEEMc3 = new TH1F("dxysigFakeEEMc3", "dxysigFakeEEMc3", 50, -50., 50.);  // 2.0-5.0
  TH1F *dxysigFakeEEMc4 = new TH1F("dxysigFakeEEMc4", "dxysigFakeEEMc4", 50, -50., 50.);  // >5
  dxysigFakeEBMc0->Sumw2();
  dxysigFakeEBMc1->Sumw2();
  dxysigFakeEBMc2->Sumw2();
  dxysigFakeEBMc3->Sumw2();
  dxysigFakeEBMc4->Sumw2();
  dxysigFakeEEMc0->Sumw2();
  dxysigFakeEEMc1->Sumw2();
  dxysigFakeEEMc2->Sumw2();
  dxysigFakeEEMc3->Sumw2();
  dxysigFakeEEMc4->Sumw2();

  // -----------------------------------------------------------------------
  // Many pT/eta bins (no weight): significance dZ
  TH1F *dzsigSignalEBMc0 = new TH1F("dzsigSignalEBMc0", "dzsigSignalEBMc0", 50, -50., 50.);  // pT: 0.5-1.0
  TH1F *dzsigSignalEBMc1 = new TH1F("dzsigSignalEBMc1", "dzsigSignalEBMc1", 50, -50., 50.);  // 1.0-1.5
  TH1F *dzsigSignalEBMc2 = new TH1F("dzsigSignalEBMc2", "dzsigSignalEBMc2", 50, -50., 50.);  // 1.5-2.0
  TH1F *dzsigSignalEBMc3 = new TH1F("dzsigSignalEBMc3", "dzsigSignalEBMc3", 50, -50., 50.);  // 2.0-5.0
  TH1F *dzsigSignalEBMc4 = new TH1F("dzsigSignalEBMc4", "dzsigSignalEBMc4", 50, -50., 50.);  // >5
  TH1F *dzsigSignalEEMc0 = new TH1F("dzsigSignalEEMc0", "dzsigSignalEEMc0", 50, -50., 50.);  // pT: 0.5-1.0
  TH1F *dzsigSignalEEMc1 = new TH1F("dzsigSignalEEMc1", "dzsigSignalEEMc1", 50, -50., 50.);  // 1.0-1.5
  TH1F *dzsigSignalEEMc2 = new TH1F("dzsigSignalEEMc2", "dzsigSignalEEMc2", 50, -50., 50.);  // 1.5-2.0
  TH1F *dzsigSignalEEMc3 = new TH1F("dzsigSignalEEMc3", "dzsigSignalEEMc3", 50, -50., 50.);  // 2.0-5.0
  TH1F *dzsigSignalEEMc4 = new TH1F("dzsigSignalEEMc4", "dzsigSignalEEMc4", 50, -50., 50.);  // >5
  dzsigSignalEBMc0->Sumw2();
  dzsigSignalEBMc1->Sumw2();
  dzsigSignalEBMc2->Sumw2();
  dzsigSignalEBMc3->Sumw2();
  dzsigSignalEBMc4->Sumw2();
  dzsigSignalEEMc0->Sumw2();
  dzsigSignalEEMc1->Sumw2();
  dzsigSignalEEMc2->Sumw2();
  dzsigSignalEEMc3->Sumw2();
  dzsigSignalEEMc4->Sumw2();
  //
  TH1F *dzsigFakeEBMc0 = new TH1F("dzsigFakeEBMc0", "dzsigFakeEBMc0", 50, -50., 50.);  // pT: 0.5-1.0
  TH1F *dzsigFakeEBMc1 = new TH1F("dzsigFakeEBMc1", "dzsigFakeEBMc1", 50, -50., 50.);  // 1.0-1.5
  TH1F *dzsigFakeEBMc2 = new TH1F("dzsigFakeEBMc2", "dzsigFakeEBMc2", 50, -50., 50.);  // 1.5-2.0
  TH1F *dzsigFakeEBMc3 = new TH1F("dzsigFakeEBMc3", "dzsigFakeEBMc3", 50, -50., 50.);  // 2.0-5.0
  TH1F *dzsigFakeEBMc4 = new TH1F("dzsigFakeEBMc4", "dzsigFakeEBMc4", 50, -50., 50.);  // >5
  TH1F *dzsigFakeEEMc0 = new TH1F("dzsigFakeEEMc0", "dzsigFakeEEMc0", 50, -50., 50.);  // pT: 0.5-1.0
  TH1F *dzsigFakeEEMc1 = new TH1F("dzsigFakeEEMc1", "dzsigFakeEEMc1", 50, -50., 50.);  // 1.0-1.5
  TH1F *dzsigFakeEEMc2 = new TH1F("dzsigFakeEEMc2", "dzsigFakeEEMc2", 50, -50., 50.);  // 1.5-2.0
  TH1F *dzsigFakeEEMc3 = new TH1F("dzsigFakeEEMc3", "dzsigFakeEEMc3", 50, -50., 50.);  // 2.0-5.0
  TH1F *dzsigFakeEEMc4 = new TH1F("dzsigFakeEEMc4", "dzsigFakeEEMc4", 50, -50., 50.);  // >5
  dzsigFakeEBMc0->Sumw2();
  dzsigFakeEBMc1->Sumw2();
  dzsigFakeEBMc2->Sumw2();
  dzsigFakeEBMc3->Sumw2();
  dzsigFakeEBMc4->Sumw2();
  dzsigFakeEEMc0->Sumw2();
  dzsigFakeEEMc1->Sumw2();
  dzsigFakeEEMc2->Sumw2();
  dzsigFakeEEMc3->Sumw2();
  dzsigFakeEEMc4->Sumw2();


  // --------------------------
  // Full pT range: distributions for data/MC checking overlap between LPT and PF
  TH1F *probePtSignalMc_LptPfOverlap  = new TH1F("probePtSignalMc_LptPfOverlap",  "probePtSignalMc_LptPfOverlap",  60,  0.,  15.);
  TH1F *probeEtaSignalMc_LptPfOverlap = new TH1F("probeEtaSignalMc_LptPfOverlap", "probeEtaSignalMc_LptPfOverlap", 40, -2.4, 2.4);
  TH1F *probeMvaSignalMc_LptPfOverlap = new TH1F("probeMvaSignalMc_LptPfOverlap", "probeMvaSignalMc_LptPfOverlap", 60, -10., 10.); 
  probePtSignalMc_LptPfOverlap->Sumw2();
  probeEtaSignalMc_LptPfOverlap->Sumw2();
  probeMvaSignalMc_LptPfOverlap->Sumw2(); 

  TH1F *probePtSignalMc_LptNotPfOverlap  = new TH1F("probePtSignalMc_LptNotPfOverlap",  "probePtSignalMc_LptNotPfOverlap",  60,  0.,  15.);
  TH1F *probeEtaSignalMc_LptNotPfOverlap = new TH1F("probeEtaSignalMc_LptNotPfOverlap", "probeEtaSignalMc_LptNotPfOverlap", 40, -2.4, 2.4);
  TH1F *probeMvaSignalMc_LptNotPfOverlap = new TH1F("probeMvaSignalMc_LptNotPfOverlap", "probeMvaSignalMc_LptNotPfOverlap", 60, -10., 10.); 
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
    TFile fileWeight("weightFile_tnpVsNani.root"); 
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
    if (hlt_9==0) continue;
    
    // Acceptance
    if (fabs(probeEta)>2.4) continue;
    if (probePt<0.5)        continue;

    // LowPt vs PF selection already applied in formatted-tnp ntuples

    // pT vs eta weight
    float ptEtaSignalWeight=1.;

    if (applyWeight) {
      for (int iBinEta=0; iBinEta<=nBinsWEta; iBinEta++) {
	for (int iBinPt=0; iBinPt<=nBinsWPt; iBinPt++) {
	  bool thisBin = false;
	  if (probePt>=minBPt[iBinPt] && probePt<maxBPt[iBinPt] && probeEta>=minBEta[iBinEta] && probeEta<maxBEta[iBinEta]) thisBin = true;
	  if (probeMatchMc==1 && thisBin) ptEtaSignalWeight = weightsSignal[iBinEta][iBinPt];
	}
      }
    }

    // Vertex distribution
    numVtxMc->Fill(numvtx);

    // matching mc-truth in MC: full pT range before weight
    if (probeMatchMc==1) { 
      if (testLPT==1) probeMvaSignalMc->Fill(probeMvaId);
      if (testLPT==0) probeMvaSignalMc->Fill(probePfmvaId);
      probePtSignalMc->Fill(probePt);
      probeEtaSignalMc->Fill(probeEta);
      probePtVsEtaSignalMc->Fill(probeEta,probePt);
      if (probeFBrem>=0 && probeFBrem<=1) probeFBremSignalMc->Fill(probeFBrem);
      probeDxySigSignalMc->Fill(probeDxySig);
      probeDzSigSignalMc->Fill(probeDzSig);
      if (studyOverlap==1) {
	if (probeIsPFOverlap==1) {
	  probePtSignalMc_LptPfOverlap->Fill(probePt);
	  probeEtaSignalMc_LptPfOverlap->Fill(probeEta);
	  probeMvaSignalMc_LptPfOverlap->Fill(probeMvaId);
	} else {
	  probePtSignalMc_LptNotPfOverlap->Fill(probePt);
	  probeEtaSignalMc_LptNotPfOverlap->Fill(probeEta);
	  probeMvaSignalMc_LptNotPfOverlap->Fill(probeMvaId);
	}
      }
    }
    
    // splitting in pT bins to check eta
    if (probeMatchMc==1) {
      if (probePt<1.0) probeEtaSignalMc0->Fill(probeEta);
      if (probePt>=1.0 && probePt<1.5) probeEtaSignalMc1->Fill(probeEta);
      if (probePt>=1.5 && probePt<2.0) probeEtaSignalMc2->Fill(probeEta);
      if (probePt>=2.0 && probePt<5.0) probeEtaSignalMc3->Fill(probeEta);
      if (probePt>=5) probeEtaSignalMc4->Fill(probeEta);
    }

    // matching mc-truth in MC: full pT range after weight
    if (probeMatchMc==1) { 
      if (testLPT==1) probeMvaSignalMcWW->Fill(probeMvaId, ptEtaSignalWeight);       
      if (testLPT==0) probeMvaSignalMcWW->Fill(probePfmvaId, ptEtaSignalWeight);       
      probePtSignalMcWW->Fill(probePt, ptEtaSignalWeight);
      probeEtaSignalMcWW->Fill(probeEta, ptEtaSignalWeight);
      probePtVsEtaSignalMcWW->Fill(probeEta,probePt,ptEtaSignalWeight);
    }

    // splitting in pT bins to check eta
    if (probeMatchMc==1) {
      if (probePt<1.0) probeEtaSignalMcWW0->Fill(probeEta,ptEtaSignalWeight);
      if (probePt>=1.0 && probePt<1.5) probeEtaSignalMcWW1->Fill(probeEta,ptEtaSignalWeight);
      if (probePt>=1.5 && probePt<2.0) probeEtaSignalMcWW2->Fill(probeEta,ptEtaSignalWeight);
      if (probePt>=2.0 && probePt<5.0) probeEtaSignalMcWW3->Fill(probeEta,ptEtaSignalWeight);
      if (probePt>=5.0) probeEtaSignalMcWW4->Fill(probeEta,ptEtaSignalWeight);
    }

    // not matching mc-truth in MC: full pT range before weight
    if (probeMatchMc==0) { 
      if (testLPT==1) probeMvaFakeMc->Fill(probeMvaId);
      if (testLPT==0) probeMvaFakeMc->Fill(probePfmvaId);
      probePtFakeMc->Fill(probePt);
      probeEtaFakeMc->Fill(probeEta);
      if (probeFBrem>=0 && probeFBrem<=1) probeFBremFakeMc->Fill(probeFBrem);
      probeDxySigFakeMc->Fill(probeDxySig);
      probeDzSigFakeMc->Fill(probeDzSig);
    }
    // splitting in pT bins to check eta
    if (probeMatchMc==0) {
      if (probePt<1) probeEtaFakeMc0->Fill(probeEta);
      if (probePt>=1.0 && probePt<1.5) probeEtaFakeMc1->Fill(probeEta);
      if (probePt>=1.5 && probePt<2.0) probeEtaFakeMc2->Fill(probeEta);
      if (probePt>=2.0 && probePt<5.0) probeEtaFakeMc3->Fill(probeEta);
      if (probePt>=5.0) probeEtaFakeMc4->Fill(probeEta);
    }

    // matching mc-truth in MC: eta/pT bins
    if (probeMatchMc==1) {  // signal

      float theId=-1;
      if (testLPT==1) theId=probeMvaId;
      if (testLPT==0) theId=probePfmvaId;

      if (fabs(probeEta)<1.5) {  // barrel
	if (probePt>=0.5 && probePt<1.0) {
	  mvaSignalEBMc0->Fill(theId);
	  mvaSignalEBMc0WW->Fill(theId,ptEtaSignalWeight);
	  if (probeFBrem>=0 && probeFBrem<=1) fbremSignalEBMc0->Fill(probeFBrem);
	  dxysigSignalEBMc0->Fill(probeDxySig); 
	  dzsigSignalEBMc0->Fill(probeDzSig); 
	}
	if (probePt>=1.0 && probePt<1.5) {
	  mvaSignalEBMc1->Fill(theId);
	  mvaSignalEBMc1WW->Fill(theId,ptEtaSignalWeight);
	  if (probeFBrem>=0 && probeFBrem<=1) fbremSignalEBMc1->Fill(probeFBrem);
	  dxysigSignalEBMc1->Fill(probeDxySig); 
	  dzsigSignalEBMc1->Fill(probeDzSig); 
	}
	if (probePt>=1.5 && probePt<2.0) {
	  mvaSignalEBMc2->Fill(theId);
	  mvaSignalEBMc2WW->Fill(theId,ptEtaSignalWeight);
	  if (probeFBrem>=0 && probeFBrem<=1) fbremSignalEBMc2->Fill(probeFBrem);
	  dxysigSignalEBMc2->Fill(probeDxySig); 
	  dzsigSignalEBMc2->Fill(probeDzSig); 
	}
	if (probePt>=2.0 && probePt<5.0) {
	  mvaSignalEBMc3->Fill(theId);
	  mvaSignalEBMc3WW->Fill(theId,ptEtaSignalWeight);
	  if (probeFBrem>=0 && probeFBrem<=1) fbremSignalEBMc3->Fill(probeFBrem);
	  dxysigSignalEBMc3->Fill(probeDxySig); 
	  dzsigSignalEBMc3->Fill(probeDzSig); 
	}
	if (probePt>=5.0) {
	  mvaSignalEBMc4->Fill(theId);
	  mvaSignalEBMc4WW->Fill(theId,ptEtaSignalWeight);
	  if (probeFBrem>=0 && probeFBrem<=1) fbremSignalEBMc4->Fill(probeFBrem);
	  dxysigSignalEBMc4->Fill(probeDxySig); 
	  dzsigSignalEBMc4->Fill(probeDzSig); 
	}

      } else {  // endcap

	if (probePt>=0.5 && probePt<1.0) {
	  mvaSignalEEMc0->Fill(theId);
	  mvaSignalEEMc0WW->Fill(theId,ptEtaSignalWeight);
	  if (probeFBrem>=0 && probeFBrem<=1) fbremSignalEEMc0->Fill(probeFBrem);
	  dxysigSignalEEMc0->Fill(probeDxySig); 
	  dzsigSignalEEMc0->Fill(probeDzSig); 
	}
	if (probePt>=1.0 && probePt<1.5) {
	  mvaSignalEEMc1->Fill(theId);
	  mvaSignalEEMc1WW->Fill(theId,ptEtaSignalWeight);
	  if (probeFBrem>=0 && probeFBrem<=1) fbremSignalEEMc1->Fill(probeFBrem);
	  dxysigSignalEEMc1->Fill(probeDxySig); 
	  dzsigSignalEEMc1->Fill(probeDzSig); 
	}
	if (probePt>=1.5 && probePt<2.0) {
	  mvaSignalEEMc2->Fill(theId);
	  mvaSignalEEMc2WW->Fill(theId,ptEtaSignalWeight);
	  if (probeFBrem>=0 && probeFBrem<=1) fbremSignalEEMc2->Fill(probeFBrem);
	  dxysigSignalEEMc2->Fill(probeDxySig); 
	  dzsigSignalEEMc2->Fill(probeDzSig); 
	}
	if (probePt>=2.0 && probePt<5.0) {
	  mvaSignalEEMc3->Fill(theId);
	  mvaSignalEEMc3WW->Fill(theId,ptEtaSignalWeight);
	  if (probeFBrem>=0 && probeFBrem<=1) fbremSignalEEMc3->Fill(probeFBrem);
	  dxysigSignalEEMc3->Fill(probeDxySig); 
	  dzsigSignalEEMc3->Fill(probeDzSig); 
	}
	if (probePt>=5.0) {
	  mvaSignalEEMc4->Fill(theId);
	  mvaSignalEEMc4WW->Fill(theId,ptEtaSignalWeight);
	  if (probeFBrem>=0 && probeFBrem<=1) fbremSignalEEMc4->Fill(probeFBrem);
	  dxysigSignalEEMc4->Fill(probeDxySig); 
	  dzsigSignalEEMc4->Fill(probeDzSig); 
	}
      }
    }

    // not matching mc-truth in MC: eta/pT bins
    if (probeMatchMc==0) { 

      float theId=-1;
      if (testLPT==1) theId=probeMvaId;
      if (testLPT==0) theId=probePfmvaId;

      if (fabs(probeEta)<1.5) {  // barrel
	if (probePt>=0.5 && probePt<1.0) {
	  mvaFakeEBMc0->Fill(theId);
	  if (probeFBrem>=0 && probeFBrem<=1) fbremFakeEBMc0->Fill(probeFBrem);
	  dxysigFakeEBMc0->Fill(probeDxySig); 
	  dzsigFakeEBMc0->Fill(probeDzSig); 
	}
	if (probePt>=1.0 && probePt<1.5) {
	  mvaFakeEBMc1->Fill(theId);
	  if (probeFBrem>=0 && probeFBrem<=1) fbremFakeEBMc1->Fill(probeFBrem);
	  dxysigFakeEBMc1->Fill(probeDxySig); 
	  dzsigFakeEBMc1->Fill(probeDzSig); 
	}
	if (probePt>=1.5 && probePt<2.0) {
	  mvaFakeEBMc2->Fill(theId);
	  if (probeFBrem>=0 && probeFBrem<=1) fbremFakeEBMc2->Fill(probeFBrem);
	  dxysigFakeEBMc2->Fill(probeDxySig); 
	  dzsigFakeEBMc2->Fill(probeDzSig); 
	}
	if (probePt>=2.0 && probePt<5.0) {
	  mvaFakeEBMc3->Fill(theId);
	  if (probeFBrem>=0 && probeFBrem<=1) fbremFakeEBMc3->Fill(probeFBrem);
	  dxysigFakeEBMc3->Fill(probeDxySig); 
	  dzsigFakeEBMc3->Fill(probeDzSig); 
	}
	if (probePt>=5.0) {
	  mvaFakeEBMc4->Fill(theId);
	  if (probeFBrem>=0 && probeFBrem<=1) fbremFakeEBMc4->Fill(probeFBrem);
	  dxysigFakeEBMc4->Fill(probeDxySig); 
	  dzsigFakeEBMc4->Fill(probeDzSig); 
	}

      } else {  // endcap

	if (probePt>=0.5 && probePt<1.0) {
	  mvaFakeEEMc0->Fill(theId);
	  if (probeFBrem>=0 && probeFBrem<=1) fbremFakeEEMc0->Fill(probeFBrem);
	  dxysigFakeEEMc0->Fill(probeDxySig); 
	  dzsigFakeEEMc0->Fill(probeDzSig); 
	}
	if (probePt>=1.0 && probePt<1.5) {
	  mvaFakeEEMc1->Fill(theId);
	  if (probeFBrem>=0 && probeFBrem<=1) fbremFakeEEMc1->Fill(probeFBrem);
	  dxysigFakeEEMc1->Fill(probeDxySig); 
	  dzsigFakeEEMc1->Fill(probeDzSig); 
	}
	if (probePt>=1.5 && probePt<2.0) {
	  mvaFakeEEMc2->Fill(theId);
	  if (probeFBrem>=0 && probeFBrem<=1) fbremFakeEEMc2->Fill(probeFBrem);
	  dxysigFakeEEMc2->Fill(probeDxySig); 
	  dzsigFakeEEMc2->Fill(probeDzSig); 
	}
	if (probePt>=2.0 && probePt<5.0) {
	  mvaFakeEEMc3->Fill(theId);
	  if (probeFBrem>=0 && probeFBrem<=1) fbremFakeEEMc3->Fill(probeFBrem);
	  dxysigFakeEEMc3->Fill(probeDxySig); 
	  dzsigFakeEEMc3->Fill(probeDzSig); 
	}
	if (probePt>=5.0) {
	  mvaFakeEEMc4->Fill(theId);
	  if (probeFBrem>=0 && probeFBrem<=1) fbremFakeEEMc4->Fill(probeFBrem);
	  dxysigFakeEEMc4->Fill(probeDxySig); 
	  dzsigFakeEEMc4->Fill(probeDzSig); 
	}
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
  probeEtaSignalMc -> SetLineWidth(2);
  probeEtaSignalMc -> SetLineColor(6);
  probeEtaFakeMc   -> SetLineWidth(2);
  probeEtaFakeMc   -> SetLineColor(4);
  //
  mvaSignalEBMc0 -> SetLineWidth(2); 
  mvaSignalEBMc1 -> SetLineWidth(2); 
  mvaSignalEBMc2 -> SetLineWidth(2); 
  mvaSignalEBMc3 -> SetLineWidth(2); 
  mvaSignalEBMc4 -> SetLineWidth(2); 
  mvaSignalEBMc0 -> SetLineColor(6); 
  mvaSignalEBMc1 -> SetLineColor(6); 
  mvaSignalEBMc2 -> SetLineColor(6); 
  mvaSignalEBMc3 -> SetLineColor(6); 
  mvaSignalEBMc4 -> SetLineColor(6); 
  //
  mvaSignalEBMc0WW -> SetLineWidth(2); 
  mvaSignalEBMc1WW -> SetLineWidth(2); 
  mvaSignalEBMc2WW -> SetLineWidth(2); 
  mvaSignalEBMc3WW -> SetLineWidth(2); 
  mvaSignalEBMc4WW -> SetLineWidth(2); 
  mvaSignalEBMc0WW -> SetLineColor(6); 
  mvaSignalEBMc1WW -> SetLineColor(6); 
  mvaSignalEBMc2WW -> SetLineColor(6); 
  mvaSignalEBMc3WW -> SetLineColor(6); 
  mvaSignalEBMc4WW -> SetLineColor(6); 
  //
  mvaSignalEEMc0 -> SetLineWidth(2); 
  mvaSignalEEMc1 -> SetLineWidth(2); 
  mvaSignalEEMc2 -> SetLineWidth(2); 
  mvaSignalEEMc3 -> SetLineWidth(2); 
  mvaSignalEEMc4 -> SetLineWidth(2); 
  mvaSignalEEMc0 -> SetLineColor(6); 
  mvaSignalEEMc1 -> SetLineColor(6); 
  mvaSignalEEMc2 -> SetLineColor(6); 
  mvaSignalEEMc3 -> SetLineColor(6); 
  mvaSignalEEMc4 -> SetLineColor(6); 
  //
  mvaSignalEEMc0WW -> SetLineWidth(2); 
  mvaSignalEEMc1WW -> SetLineWidth(2); 
  mvaSignalEEMc2WW -> SetLineWidth(2); 
  mvaSignalEEMc3WW -> SetLineWidth(2); 
  mvaSignalEEMc4WW -> SetLineWidth(2); 
  mvaSignalEEMc0WW -> SetLineColor(6); 
  mvaSignalEEMc1WW -> SetLineColor(6); 
  mvaSignalEEMc2WW -> SetLineColor(6); 
  mvaSignalEEMc3WW -> SetLineColor(6); 
  mvaSignalEEMc4WW -> SetLineColor(6); 
  //
  mvaFakeEBMc0 -> SetLineWidth(2); 
  mvaFakeEBMc1 -> SetLineWidth(2); 
  mvaFakeEBMc2 -> SetLineWidth(2); 
  mvaFakeEBMc3 -> SetLineWidth(2); 
  mvaFakeEBMc4 -> SetLineWidth(2); 
  mvaFakeEBMc0 -> SetLineColor(4); 
  mvaFakeEBMc1 -> SetLineColor(4); 
  mvaFakeEBMc2 -> SetLineColor(4); 
  mvaFakeEBMc3 -> SetLineColor(4); 
  mvaFakeEBMc4 -> SetLineColor(4); 
  //
  mvaFakeEEMc0 -> SetLineWidth(2); 
  mvaFakeEEMc1 -> SetLineWidth(2); 
  mvaFakeEEMc2 -> SetLineWidth(2); 
  mvaFakeEEMc3 -> SetLineWidth(2); 
  mvaFakeEEMc4 -> SetLineWidth(2); 
  mvaFakeEEMc0 -> SetLineColor(4); 
  mvaFakeEEMc1 -> SetLineColor(4); 
  mvaFakeEEMc2 -> SetLineColor(4); 
  mvaFakeEEMc3 -> SetLineColor(4); 
  mvaFakeEEMc4 -> SetLineColor(4); 


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
  numVtxMc->Write();     
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
  probeEtaSignalMc0->Write();
  probeEtaFakeMc0->Write();
  probeEtaSignalMc1->Write();
  probeEtaFakeMc1->Write();
  probeEtaSignalMc2->Write();
  probeEtaFakeMc2->Write();
  probeEtaSignalMc3->Write();
  probeEtaFakeMc3->Write();
  probeEtaSignalMc4->Write();
  probeEtaFakeMc4->Write();
  //
  probePtSignalMcWW->Write();
  probeEtaSignalMcWW->Write();
  probeMvaSignalMcWW->Write(); 
  //
  probeEtaSignalMcWW0->Write();
  probeEtaSignalMcWW1->Write();
  probeEtaSignalMcWW2->Write();
  probeEtaSignalMcWW3->Write();
  probeEtaSignalMcWW4->Write();
  //
  probePtVsEtaSignalMc->Write();
  probePtVsEtaSignalMcWW->Write();
  //
  probeFBremSignalMc->Write();
  probeFBremFakeMc->Write();
  probeDxySigSignalMc->Write();
  probeDxySigFakeMc->Write();
  probeDzSigSignalMc->Write();
  probeDzSigFakeMc->Write();
  //
  mvaSignalEBMc0->Write();
  mvaSignalEBMc1->Write();
  mvaSignalEBMc2->Write();
  mvaSignalEBMc3->Write();
  mvaSignalEBMc4->Write();
  mvaSignalEEMc0->Write();
  mvaSignalEEMc1->Write();
  mvaSignalEEMc2->Write();
  mvaSignalEEMc3->Write();
  mvaSignalEEMc4->Write();
  //
  mvaSignalEBMc0WW->Write();
  mvaSignalEBMc1WW->Write();
  mvaSignalEBMc2WW->Write();
  mvaSignalEBMc3WW->Write();
  mvaSignalEBMc4WW->Write();
  mvaSignalEEMc0WW->Write();
  mvaSignalEEMc1WW->Write();
  mvaSignalEEMc2WW->Write();
  mvaSignalEEMc3WW->Write();
  mvaSignalEEMc4WW->Write();
  //
  mvaFakeEBMc0->Write();
  mvaFakeEBMc1->Write();
  mvaFakeEBMc2->Write();
  mvaFakeEBMc3->Write();
  mvaFakeEBMc4->Write();
  mvaFakeEEMc0->Write();
  mvaFakeEEMc1->Write();
  mvaFakeEEMc2->Write();
  mvaFakeEEMc3->Write();
  mvaFakeEEMc4->Write();
  //
  fbremSignalEBMc0->Write();
  fbremSignalEBMc1->Write();
  fbremSignalEBMc2->Write();
  fbremSignalEBMc3->Write();
  fbremSignalEBMc4->Write();
  fbremSignalEEMc0->Write();
  fbremSignalEEMc1->Write();
  fbremSignalEEMc2->Write();
  fbremSignalEEMc3->Write();
  fbremSignalEEMc4->Write();
  //
  fbremFakeEBMc0->Write();
  fbremFakeEBMc1->Write();
  fbremFakeEBMc2->Write();
  fbremFakeEBMc3->Write();
  fbremFakeEBMc4->Write();
  fbremFakeEEMc0->Write();
  fbremFakeEEMc1->Write();
  fbremFakeEEMc2->Write();
  fbremFakeEEMc3->Write();
  fbremFakeEEMc4->Write();
  //
  dxysigSignalEBMc0->Write();
  dxysigSignalEBMc1->Write();
  dxysigSignalEBMc2->Write();
  dxysigSignalEBMc3->Write();
  dxysigSignalEBMc4->Write();
  dxysigSignalEEMc0->Write();
  dxysigSignalEEMc1->Write();
  dxysigSignalEEMc2->Write();
  dxysigSignalEEMc3->Write();
  dxysigSignalEEMc4->Write();
  //
  dxysigFakeEBMc0->Write();
  dxysigFakeEBMc1->Write();
  dxysigFakeEBMc2->Write();
  dxysigFakeEBMc3->Write();
  dxysigFakeEBMc4->Write();
  dxysigFakeEEMc0->Write();
  dxysigFakeEEMc1->Write();
  dxysigFakeEEMc2->Write();
  dxysigFakeEEMc3->Write();
  dxysigFakeEEMc4->Write();
  //
  dzsigSignalEBMc0->Write();
  dzsigSignalEBMc1->Write();
  dzsigSignalEBMc2->Write();
  dzsigSignalEBMc3->Write();
  dzsigSignalEBMc4->Write();
  dzsigSignalEEMc0->Write();
  dzsigSignalEEMc1->Write();
  dzsigSignalEEMc2->Write();
  dzsigSignalEEMc3->Write();
  dzsigSignalEEMc4->Write();
  //
  dzsigFakeEBMc0->Write();
  dzsigFakeEBMc1->Write();
  dzsigFakeEBMc2->Write();
  dzsigFakeEBMc3->Write();
  dzsigFakeEBMc4->Write();
  dzsigFakeEEMc0->Write();
  dzsigFakeEEMc1->Write();
  dzsigFakeEEMc2->Write();
  dzsigFakeEEMc3->Write();
  dzsigFakeEEMc4->Write();

  myFile.Close();


  // -----------------------------------------------------------------------
  // Rebin
  mvaSignalEBMc0->Rebin();
  mvaSignalEBMc1->Rebin();
  mvaSignalEBMc2->Rebin();
  mvaSignalEBMc3->Rebin();
  mvaSignalEBMc4->Rebin();
  mvaSignalEEMc0->Rebin();
  mvaSignalEEMc1->Rebin();
  mvaSignalEEMc2->Rebin();
  mvaSignalEEMc3->Rebin();
  mvaSignalEEMc4->Rebin();
  mvaFakeEBMc0->Rebin();
  mvaFakeEBMc1->Rebin();
  mvaFakeEBMc2->Rebin();
  mvaFakeEBMc3->Rebin();
  mvaFakeEBMc4->Rebin();
  mvaFakeEEMc0->Rebin();
  mvaFakeEEMc1->Rebin();
  mvaFakeEEMc2->Rebin();
  mvaFakeEEMc3->Rebin();
  mvaFakeEEMc4->Rebin();


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
  leg->AddEntry(probeMvaSignalMc, "Matching truth", "lp");
  leg->AddEntry(probeMvaFakeMc,   "Not matching truth", "lp");
  //
  TLegend *legB;
  legB = new TLegend(0.40,0.65,0.75,0.90);
  legB->SetFillStyle(0);
  legB->SetBorderSize(0);
  legB->SetTextSize(0.05);
  legB->SetFillColor(0);
  legB->AddEntry(probeMvaSignalMc, "Matching truth", "lp");
  legB->AddEntry(probeMvaFakeMc,   "Not matching truth", "lp");
  //
  TLegend *legC;
  legC = new TLegend(0.25,0.15,0.65,0.30);
  legC->SetFillStyle(0);
  legC->SetBorderSize(0);
  legC->SetTextSize(0.05);
  legC->SetFillColor(0);
  legC->AddEntry(probeMvaSignalMc, "Matching truth", "lp");
  legC->AddEntry(probeMvaFakeMc,   "Not matching truth", "lp");

  TCanvas cmvamc("cmvamc","cmvamc",1);
  probeMvaSignalMc->SetTitle("");
  probeMvaFakeMc->SetTitle("");
  probeMvaSignalMc->GetXaxis()->SetTitle("Id BDT");
  probeMvaFakeMc->GetXaxis()->SetTitle("Id BDT");
  probeMvaSignalMc->DrawNormalized("hist");
  probeMvaFakeMc->DrawNormalized("samehist");
  leg->Draw();
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

  TCanvas cmvaeb0("cmvaeb0","cmvaeb0",1);
  mvaSignalEBMc0 -> SetTitle("EB, 0.5 < pT < 1");
  mvaFakeEBMc0   -> SetTitle("EB, 0.5 < pT < 1");
  mvaSignalEBMc0 -> GetXaxis()->SetTitle("Id BDT");
  mvaFakeEBMc0   -> GetXaxis()->SetTitle("Id BDT");
  mvaSignalEBMc0->DrawNormalized("hist");
  mvaFakeEBMc0->DrawNormalized("samehist");
  leg->Draw();
  cmvaeb0.SaveAs("outputBDT_matchVsFake_withTnP_EB0.png");  
  //
  TCanvas cmvaeb1("cmvaeb1","cmvaeb1",1);
  mvaSignalEBMc1 -> SetTitle("EB, 1.0 < pT < 1.5");
  mvaFakeEBMc1   -> SetTitle("EB, 1.0 < pT < 1.5");
  mvaSignalEBMc1 -> GetXaxis()->SetTitle("Id BDT");
  mvaFakeEBMc1   -> GetXaxis()->SetTitle("Id BDT");
  mvaSignalEBMc1->DrawNormalized("hist");
  mvaFakeEBMc1->DrawNormalized("samehist");
  leg->Draw();
  cmvaeb1.SaveAs("outputBDT_matchVsFake_withTnP_EB1.png");  
  //
  TCanvas cmvaeb2("cmvaeb2","cmvaeb2",1);
  mvaSignalEBMc2 -> SetTitle("EB, 1.5 < pT < 2.0");
  mvaFakeEBMc2   -> SetTitle("EB, 1.5 < pT < 2.0");
  mvaSignalEBMc2 -> GetXaxis()->SetTitle("Id BDT");
  mvaFakeEBMc2   -> GetXaxis()->SetTitle("Id BDT");
  mvaSignalEBMc2->DrawNormalized("hist");
  mvaFakeEBMc2->DrawNormalized("samehist");
  leg->Draw();
  cmvaeb2.SaveAs("outputBDT_matchVsFake_withTnP_EB2.png");  
  //
  TCanvas cmvaeb3("cmvaeb3","cmvaeb3",1);
  mvaSignalEBMc3 -> SetTitle("EB, 2.0 < pT < 5.0");
  mvaFakeEBMc3   -> SetTitle("EB, 2.0 < pT < 5.0");
  mvaSignalEBMc3 -> GetXaxis()->SetTitle("Id BDT");
  mvaFakeEBMc3   -> GetXaxis()->SetTitle("Id BDT");
  mvaSignalEBMc3->DrawNormalized("hist");
  mvaFakeEBMc3->DrawNormalized("samehist");
  leg->Draw();
  cmvaeb3.SaveAs("outputBDT_matchVsFake_withTnP_EB3.png");  
  //
  TCanvas cmvaeb4("cmvaeb4","cmvaeb4",1);
  mvaSignalEBMc4 -> SetTitle("EB, pT > 5");
  mvaFakeEBMc4   -> SetTitle("EB, pT > 5");
  mvaSignalEBMc4 -> GetXaxis()->SetTitle("Id BDT");
  mvaFakeEBMc4   -> GetXaxis()->SetTitle("Id BDT");
  mvaSignalEBMc4->DrawNormalized("hist");
  mvaFakeEBMc4->DrawNormalized("samehist");
  leg->Draw();
  cmvaeb4.SaveAs("outputBDT_matchVsFake_withTnP_EB4.png");  

  TCanvas cmvaee0("cmvaee0","cmvaee0",1);
  mvaSignalEEMc0 -> SetTitle("EE, 0.5 < pT < 1");
  mvaFakeEEMc0   -> SetTitle("EE, 0.5 < pT < 1");
  mvaSignalEEMc0 -> GetXaxis()->SetTitle("Id BDT");
  mvaFakeEEMc0   -> GetXaxis()->SetTitle("Id BDT");
  mvaFakeEEMc0->DrawNormalized("hist");
  mvaSignalEEMc0->DrawNormalized("samehist");
  leg->Draw();
  cmvaee0.SaveAs("outputBDT_matchVsFake_withTnP_EE0.png");  
  //
  TCanvas cmvaee1("cmvaee1","cmvaee1",1);
  mvaSignalEEMc1 -> SetTitle("EE, 1.0 < pT < 1.5");
  mvaFakeEEMc1   -> SetTitle("EE, 1.0 < pT < 1.5");
  mvaSignalEEMc1 -> GetXaxis()->SetTitle("Id BDT");
  mvaFakeEEMc1   -> GetXaxis()->SetTitle("Id BDT");
  mvaSignalEEMc1->DrawNormalized("hist");
  mvaFakeEEMc1->DrawNormalized("samehist");
  leg->Draw();
  cmvaee1.SaveAs("outputBDT_matchVsFake_withTnP_EE1.png");  
  //
  TCanvas cmvaee2("cmvaee2","cmvaee2",1);
  mvaSignalEEMc2 -> SetTitle("EE, 1.5 < pT < 2.0");
  mvaFakeEEMc2   -> SetTitle("EE, 1.5 < pT < 2.0");
  mvaSignalEEMc2 -> GetXaxis()->SetTitle("Id BDT");
  mvaFakeEEMc2   -> GetXaxis()->SetTitle("Id BDT");
  mvaSignalEEMc2->DrawNormalized("hist");
  mvaFakeEEMc2->DrawNormalized("samehist");
  leg->Draw();
  cmvaee2.SaveAs("outputBDT_matchVsFake_withTnP_EE2.png");  
  //
  TCanvas cmvaee3("cmvaee3","cmvaee3",1);
  mvaSignalEEMc3 -> SetTitle("EE, 2.0 < pT < 5.0");
  mvaFakeEEMc3   -> SetTitle("EE, 2.0 < pT < 5.0");
  mvaSignalEEMc3 -> GetXaxis()->SetTitle("Id BDT");
  mvaFakeEEMc3   -> GetXaxis()->SetTitle("Id BDT");
  mvaSignalEEMc3->DrawNormalized("hist");
  mvaFakeEEMc3->DrawNormalized("samehist");
  leg->Draw();
  cmvaee3.SaveAs("outputBDT_matchVsFake_withTnP_EE3.png");  
  //
  TCanvas cmvaee4("cmvaee4","cmvaee4",1);
  mvaSignalEEMc4 -> SetTitle("EE, pT > 5");
  mvaFakeEEMc4   -> SetTitle("EE, pT > 5");
  mvaSignalEEMc4 -> GetXaxis()->SetTitle("Id BDT");
  mvaFakeEEMc4   -> GetXaxis()->SetTitle("Id BDT");
  mvaSignalEEMc4->DrawNormalized("hist");
  mvaFakeEEMc4->DrawNormalized("samehist");
  leg->Draw();
  cmvaee4.SaveAs("outputBDT_matchVsFake_withTnP_EE4.png");  
}
