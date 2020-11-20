#define prepareInputsFromNaniInMc_cxx
#include "prepareInputsFromNaniInMc.h"

// ROOT includes 
#include "TVector3.h"
#include "TLorentzVector.h"
#include <TROOT.h>
#include <TStyle.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <iostream>

using namespace std;

// To be run on nanoAODs MC to make basic distributions (eta, pT, Id-Output)
// mainly to compare the kinematics with and without TnP selection
// Select signal and background based on match with MC-truth

void prepareInputsFromNaniInMc::Loop(int testLowPt, bool studyOverlap)
{
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;

  // -------------------------------------------------------------------------
  // Histos
  TH1F *mvaSignalMc = new TH1F("mvaSignalMc", "mvaSignalMc", 60, -10., 10.);
  TH1F *mvaFakeMc   = new TH1F("mvaFakeMc",   "mvaFakeMc",   60, -10., 10.);
  TH1F *ptSignalMc  = new TH1F("ptSignalMc",  "ptSignalMc",  60,   0., 15.);
  TH1F *ptFakeMc    = new TH1F("ptFakeMc",    "ptFakeMc",    60,   0., 15.);
  TH1F *etaSignalMc = new TH1F("etaSignalMc", "etaSignalMc", 40, -2.4, 2.4);
  TH1F *etaFakeMc   = new TH1F("etaFakeMc",   "etaFakeMc",   40, -2.4, 2.4);
  mvaSignalMc -> Sumw2();
  mvaFakeMc   -> Sumw2();
  ptSignalMc  -> Sumw2();    
  ptFakeMc    -> Sumw2();    
  etaSignalMc -> Sumw2();    
  etaFakeMc   -> Sumw2();    
  // 
  TH1F *mvaSignalMc_LptPfOverlap    = new TH1F("mvaSignalMc_LptPfOverlap",    "mvaSignalMc_LptPfOverlap",    60, -10., 10.); 
  TH1F *ptSignalMc_LptPfOverlap     = new TH1F("ptSignalMc_LptPfOverlap",     "ptSignalMc_LptPfOverlap",     60,   0., 15.);
  TH1F *etaSignalMc_LptPfOverlap    = new TH1F("etaSignalMc_LptPfOverlap",    "etaSignalMc_LptPfOverlap",    40, -2.4, 2.4);
  TH1F *mvaSignalMc_LptNotPfOverlap = new TH1F("mvaSignalMc_LptNotPfOverlap", "mvaSignalMc_LptNotPfOverlap", 60, -10., 10.); 
  TH1F *ptSignalMc_LptNotPfOverlap  = new TH1F("ptSignalMc_LptNotPfOverlap",  "ptSignalMc_LptNotPfOverlap",  60,   0., 15.);
  TH1F *etaSignalMc_LptNotPfOverlap = new TH1F("etaSignalMc_LptNotPfOverlap", "etaSignalMc_LptNotPfOverlap", 40, -2.4, 2.4);
  TH1F *mvaFakeMc_LptPfOverlap      = new TH1F("mvaFakeMc_LptPfOverlap",      "mvaFakeMc_LptPfOverlap",      60, -10., 10.); 
  TH1F *ptFakeMc_LptPfOverlap       = new TH1F("ptFakeMc_LptPfOverlap",       "ptFakeMc_LptPfOverlap",       60,   0., 15.);
  TH1F *etaFakeMc_LptPfOverlap      = new TH1F("etaFakeMc_LptPfOverlap",      "etaFakeMc_LptPfOverlap",      40, -2.4, 2.4);
  TH1F *mvaFakeMc_LptNotPfOverlap   = new TH1F("mvaFakeMc_LptNotPfOverlap",   "mvaFakeMc_LptNotPfOverlap",   60, -10., 10.); 
  TH1F *ptFakeMc_LptNotPfOverlap    = new TH1F("ptFakeMc_LptNotPfOverlap",    "ptFakeMc_LptNotPfOverlap",    60,   0., 15.);
  TH1F *etaFakeMc_LptNotPfOverlap   = new TH1F("etaFakeMc_LptNotPfOverlap",   "etaFakeMc_LptNotPfOverlap",   40, -2.4, 2.4);
  mvaSignalMc_LptPfOverlap->Sumw2();  
  ptSignalMc_LptPfOverlap->Sumw2();
  etaSignalMc_LptPfOverlap->Sumw2();  
  mvaSignalMc_LptNotPfOverlap->Sumw2();  
  ptSignalMc_LptNotPfOverlap->Sumw2();
  etaSignalMc_LptNotPfOverlap->Sumw2();  
  mvaFakeMc_LptPfOverlap->Sumw2();  
  ptFakeMc_LptPfOverlap->Sumw2();
  etaFakeMc_LptPfOverlap->Sumw2();  
  mvaFakeMc_LptNotPfOverlap->Sumw2();  
  ptFakeMc_LptNotPfOverlap->Sumw2();
  etaFakeMc_LptNotPfOverlap->Sumw2();  
  //
  TH1F *etaSignalMc0 = new TH1F("etaSignalMc0", "etaSignalMc0", 40, -2.4, 2.4);
  TH1F *etaFakeMc0   = new TH1F("etaFakeMc0",   "etaFakeMc0",   40, -2.4, 2.4);
  TH1F *etaSignalMc1 = new TH1F("etaSignalMc1", "etaSignalMc1", 40, -2.4, 2.4);
  TH1F *etaFakeMc1   = new TH1F("etaFakeMc1",   "etaFakeMc1",   40, -2.4, 2.4);
  TH1F *etaSignalMc2 = new TH1F("etaSignalMc2", "etaSignalMc2", 40, -2.4, 2.4);
  TH1F *etaFakeMc2   = new TH1F("etaFakeMc2",   "etaFakeMc2",   40, -2.4, 2.4);
  TH1F *etaSignalMc3 = new TH1F("etaSignalMc3", "etaSignalMc3", 40, -2.4, 2.4);
  TH1F *etaFakeMc3   = new TH1F("etaFakeMc3",   "etaFakeMc3",   40, -2.4, 2.4);
  TH1F *etaSignalMc4 = new TH1F("etaSignalMc4", "etaSignalMc4", 40, -2.4, 2.4);
  TH1F *etaFakeMc4   = new TH1F("etaFakeMc4",   "etaFakeMc4",   40, -2.4, 2.4);
  etaSignalMc0 -> Sumw2();    
  etaFakeMc0   -> Sumw2();    
  etaSignalMc1 -> Sumw2();    
  etaFakeMc1   -> Sumw2();    
  etaSignalMc2 -> Sumw2();    
  etaFakeMc2   -> Sumw2();    
  etaSignalMc3 -> Sumw2();    
  etaFakeMc3   -> Sumw2();    
  etaSignalMc4 -> Sumw2();    
  etaFakeMc4   -> Sumw2();    

  TH2F *ptVsEtaSignalMc = new TH2F("ptVsEtaSignalMc","ptVsEtaSignalMc", 40, -2.4, 2.4, 60, 0., 15.);
  TH2F *ptVsEtaFakeMc   = new TH2F("ptVsEtaFakeMc",  "ptVsEtaFakeMc",   40, -2.4, 2.4, 60, 0., 15.);
  ptVsEtaSignalMc->Sumw2();
  ptVsEtaFakeMc->Sumw2();

  // Many pT/eta bins
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


  // -------------------------------------------------------------------------
  // Loop over entries
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    // Look for mc truth e-/e+
    int mcPos = -999;
    int mcEle = -999;
    for (int iGen=0; iGen<nGenPart; iGen++) {
      if ( GenPart_pdgId[iGen]==11 && GenPart_status[iGen]==1 )  mcPos = iGen;
      if ( GenPart_pdgId[iGen]==-11 && GenPart_status[iGen]==1 ) mcEle = iGen;
    } 
    TVector3 mcEleTV3(0,0,0);
    if (mcEle>=0) mcEleTV3.SetPtEtaPhi(GenPart_pt[mcEle], GenPart_eta[mcEle], GenPart_phi[mcEle]);
    TVector3 mcPosTV3(0,0,0);
    if (mcPos>=0) mcPosTV3.SetPtEtaPhi(GenPart_pt[mcPos], GenPart_eta[mcPos], GenPart_phi[mcPos]);
    if (mcEle<0 || mcPos<0) cout << "Error: mc-truth not found" << endl;


    // Loop over electrons
    for (int iEle=0; iEle<nElectron; iEle++) {

      // HLT - not needed in MC
      // int iHLT_Mu9_IP6 = (int)HLT_Mu9_IP6;
      // if (iHLT_Mu9_IP6==0) continue;

      // LowPt or PF only    
      if (testLowPt==1 && Electron_isLowPt[iEle]==0) continue;
      if (testLowPt==0 && Electron_isPF[iEle]==0) continue;

      // Sanity check  
      if (testLowPt==0 && studyOverlap==1) { cout << "should not happen!" << endl; continue; }

      // Acceptance
      if (fabs(Electron_eta[iEle])>2.4) continue;
      if (Electron_pt[iEle]<0.5)        continue;

      // Electron candidate
      TVector3 recoEleTV3(0,0,0); 
      recoEleTV3.SetPtEtaPhi(Electron_pt[iEle], Electron_eta[iEle], Electron_phi[iEle]);
      
      float minDR = 999.;
      float deltaPtOverPt = 999.;
      if (mcEle>=0 && (mcEleTV3.DeltaR(recoEleTV3))<minDR) { 
	minDR = mcEleTV3.DeltaR(recoEleTV3);
	deltaPtOverPt = fabs(Electron_pt[iEle]-GenPart_pt[mcEle])/GenPart_pt[mcEle];
      }
      if (mcPos>=0 && (mcPosTV3.DeltaR(recoEleTV3))<minDR) { 
	minDR = mcPosTV3.DeltaR(recoEleTV3);
	deltaPtOverPt = fabs(Electron_pt[iEle]-GenPart_pt[mcPos])/GenPart_pt[mcPos];
      }
      
      int trueEle = -1;
      if (minDR<0.03 && deltaPtOverPt<0.5) trueEle=1;

      // Full eta/pT range
      if (trueEle==1) { 
	if (testLowPt==1) mvaSignalMc->Fill(Electron_mvaId[iEle]);    
	if (testLowPt==0) mvaSignalMc->Fill(Electron_pfmvaId[iEle]);    
	ptSignalMc->Fill(Electron_pt[iEle]);
	etaSignalMc->Fill(Electron_eta[iEle]);
	ptVsEtaSignalMc->Fill(Electron_eta[iEle],Electron_pt[iEle]);
	// 
	// LPT / PF overlap study
	if (studyOverlap==1) {
	  if (Electron_isLowPt[iEle]==1 && Electron_isPFoverlap[iEle]==1) {
	    mvaSignalMc_LptPfOverlap -> Fill(Electron_mvaId[iEle]);
	    ptSignalMc_LptPfOverlap  -> Fill(Electron_pt[iEle]);
	    etaSignalMc_LptPfOverlap -> Fill(Electron_eta[iEle]);
	  } 
	  if (Electron_isLowPt[iEle]==1 && Electron_isPFoverlap[iEle]==0) {
	    mvaSignalMc_LptNotPfOverlap -> Fill(Electron_mvaId[iEle]);
	    ptSignalMc_LptNotPfOverlap  -> Fill(Electron_pt[iEle]);
	    etaSignalMc_LptNotPfOverlap -> Fill(Electron_eta[iEle]);
	  }
	}
	//
	if (Electron_pt[iEle]<1.0) etaSignalMc0->Fill(Electron_eta[iEle]); 
	if (Electron_pt[iEle]>=1.0 && Electron_pt[iEle]<1.5) etaSignalMc1->Fill(Electron_eta[iEle]); 
	if (Electron_pt[iEle]>=1.5 && Electron_pt[iEle]<2.0) etaSignalMc2->Fill(Electron_eta[iEle]); 
	if (Electron_pt[iEle]>=2.0 && Electron_pt[iEle]<5.0) etaSignalMc3->Fill(Electron_eta[iEle]); 
	if (Electron_pt[iEle]>=5.0) etaSignalMc4->Fill(Electron_eta[iEle]); 
      } else {
	if (testLowPt==1) mvaFakeMc->Fill(Electron_mvaId[iEle]);
	if (testLowPt==0) mvaFakeMc->Fill(Electron_pfmvaId[iEle]);
	ptFakeMc->Fill(Electron_pt[iEle]);
	etaFakeMc->Fill(Electron_eta[iEle]);
	ptVsEtaFakeMc->Fill(Electron_eta[iEle],Electron_pt[iEle]);
	//
	// LPT / PF overlap study
	if (studyOverlap==1) {
	  if (Electron_isLowPt[iEle]==1 && Electron_isPFoverlap[iEle]==1) {
	    mvaFakeMc_LptPfOverlap -> Fill(Electron_mvaId[iEle]);
	    ptFakeMc_LptPfOverlap  -> Fill(Electron_pt[iEle]);
	    etaFakeMc_LptPfOverlap -> Fill(Electron_eta[iEle]);
	  } 
	  if (Electron_isLowPt[iEle]==1 && Electron_isPFoverlap[iEle]==0) {
	    mvaFakeMc_LptNotPfOverlap -> Fill(Electron_mvaId[iEle]);
	    ptFakeMc_LptNotPfOverlap  -> Fill(Electron_pt[iEle]);
	    etaFakeMc_LptNotPfOverlap -> Fill(Electron_eta[iEle]);
	  }
	}
	//
	if (Electron_pt[iEle]<1.0) etaFakeMc0->Fill(Electron_eta[iEle]); 
	if (Electron_pt[iEle]>=1.0 && Electron_pt[iEle]<1.5) etaFakeMc1->Fill(Electron_eta[iEle]); 
	if (Electron_pt[iEle]>=1.5 && Electron_pt[iEle]<2.0) etaFakeMc2->Fill(Electron_eta[iEle]); 
	if (Electron_pt[iEle]>=2.0 && Electron_pt[iEle]<5.0) etaFakeMc3->Fill(Electron_eta[iEle]); 
	if (Electron_pt[iEle]>=5.0) etaFakeMc4->Fill(Electron_eta[iEle]); 
      }

      // Eta/pT bins
      if (trueEle==1) { 
	
	float theId = -1;
	if (testLowPt==1) theId=Electron_mvaId[iEle];
	if (testLowPt==0) theId=Electron_pfmvaId[iEle];

	if (fabs(Electron_eta[iEle])<1.5) {
	  if (Electron_pt[iEle]>=0.5 && Electron_pt[iEle]<1.0) mvaSignalEBMc0->Fill(theId);
	  if (Electron_pt[iEle]>=1.0 && Electron_pt[iEle]<1.5) mvaSignalEBMc1->Fill(theId);
	  if (Electron_pt[iEle]>=1.5 && Electron_pt[iEle]<2.0) mvaSignalEBMc2->Fill(theId);
	  if (Electron_pt[iEle]>=2.0 && Electron_pt[iEle]<5.0) mvaSignalEBMc3->Fill(theId);
	  if (Electron_pt[iEle]>=5.0) mvaSignalEBMc4->Fill(theId);
	} else {
	  if (Electron_pt[iEle]>=0.5 && Electron_pt[iEle]<1.0) mvaSignalEEMc0->Fill(theId);
	  if (Electron_pt[iEle]>=1.0 && Electron_pt[iEle]<1.5) mvaSignalEEMc1->Fill(theId);
	  if (Electron_pt[iEle]>=1.5 && Electron_pt[iEle]<2.0) mvaSignalEEMc2->Fill(theId);
	  if (Electron_pt[iEle]>=2.0 && Electron_pt[iEle]<5.0) mvaSignalEEMc3->Fill(theId);
	  if (Electron_pt[iEle]>=5.0) mvaSignalEEMc4->Fill(theId);
	}
      } else {

	float theId = -1;
	if (testLowPt==1) theId=Electron_mvaId[iEle];
	if (testLowPt==0) theId=Electron_pfmvaId[iEle];

	if (fabs(Electron_eta[iEle])<1.5) {
	  if (Electron_pt[iEle]>=0.5 && Electron_pt[iEle]<1.0) mvaFakeEBMc0->Fill(theId);
	  if (Electron_pt[iEle]>=1.0 && Electron_pt[iEle]<1.5) mvaFakeEBMc1->Fill(theId);
	  if (Electron_pt[iEle]>=1.5 && Electron_pt[iEle]<2.0) mvaFakeEBMc2->Fill(theId);
	  if (Electron_pt[iEle]>=2.0 && Electron_pt[iEle]<5.0) mvaFakeEBMc3->Fill(theId);
	  if (Electron_pt[iEle]>=5.0) mvaFakeEBMc4->Fill(theId);
	} else {
	  if (Electron_pt[iEle]>=0.5 && Electron_pt[iEle]<1.0) mvaFakeEEMc0->Fill(theId);
	  if (Electron_pt[iEle]>=1.0 && Electron_pt[iEle]<1.5) mvaFakeEEMc1->Fill(theId);
	  if (Electron_pt[iEle]>=1.5 && Electron_pt[iEle]<2.0) mvaFakeEEMc2->Fill(theId);
	  if (Electron_pt[iEle]>=2.0 && Electron_pt[iEle]<5.0) mvaFakeEEMc3->Fill(theId);
	  if (Electron_pt[iEle]>=5.0) mvaFakeEEMc4->Fill(theId);
	}
      }

    } // Loop over electrons

  } // Loop over entries


  // -------------------------------------------------------------------------
  // Cosmetics
  mvaSignalMc -> SetLineWidth(2);
  mvaSignalMc -> SetLineColor(6);
  mvaFakeMc   -> SetLineWidth(2);
  mvaFakeMc   -> SetLineColor(4);
  ptSignalMc  -> SetLineWidth(2);
  ptSignalMc  -> SetLineColor(6);
  ptFakeMc    -> SetLineWidth(2);
  ptFakeMc    -> SetLineColor(4);
  etaSignalMc -> SetLineWidth(2);
  etaSignalMc -> SetLineColor(6);
  etaFakeMc   -> SetLineWidth(2);
  etaFakeMc   -> SetLineColor(4);
  //
  mvaSignalMc_LptPfOverlap -> SetLineWidth(2);  
  mvaSignalMc_LptPfOverlap -> SetLineColor(6);  
  mvaFakeMc_LptPfOverlap   -> SetLineWidth(2);  
  mvaFakeMc_LptPfOverlap   -> SetLineColor(6);  
  ptSignalMc_LptPfOverlap  -> SetLineWidth(2);  
  ptSignalMc_LptPfOverlap  -> SetLineColor(6);  
  ptFakeMc_LptPfOverlap    -> SetLineWidth(2);  
  ptFakeMc_LptPfOverlap    -> SetLineColor(6);  
  etaSignalMc_LptPfOverlap -> SetLineWidth(2);  
  etaSignalMc_LptPfOverlap -> SetLineColor(6);  
  etaFakeMc_LptPfOverlap   -> SetLineWidth(2);  
  etaFakeMc_LptPfOverlap   -> SetLineColor(6);  
  mvaSignalMc_LptNotPfOverlap -> SetLineWidth(2);  
  mvaSignalMc_LptNotPfOverlap -> SetLineColor(4);  
  mvaFakeMc_LptNotPfOverlap   -> SetLineWidth(2);  
  mvaFakeMc_LptNotPfOverlap   -> SetLineColor(4);  
  ptSignalMc_LptNotPfOverlap  -> SetLineWidth(2);  
  ptSignalMc_LptNotPfOverlap  -> SetLineColor(4);  
  ptFakeMc_LptNotPfOverlap    -> SetLineWidth(2);  
  ptFakeMc_LptNotPfOverlap    -> SetLineColor(4);  
  etaSignalMc_LptNotPfOverlap -> SetLineWidth(2);  
  etaSignalMc_LptNotPfOverlap -> SetLineColor(4);  
  etaFakeMc_LptNotPfOverlap   -> SetLineWidth(2);  
  etaFakeMc_LptNotPfOverlap   -> SetLineColor(4);  
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


  // -------------------------------------------------------------------------
  // Save Outputs
  TFile myfile("myFileFromNani.root","RECREATE");
  myfile.cd();
  mvaSignalMc    -> Write();
  mvaFakeMc      -> Write();
  ptSignalMc     -> Write();    
  ptFakeMc       -> Write();    
  etaSignalMc    -> Write();    
  etaSignalMc0    -> Write();    
  etaSignalMc1    -> Write();    
  etaSignalMc2    -> Write();    
  etaSignalMc3    -> Write();    
  etaSignalMc4    -> Write();    
  etaFakeMc      -> Write();    
  etaFakeMc0      -> Write();    
  etaFakeMc1      -> Write();    
  etaFakeMc2      -> Write();    
  etaFakeMc3      -> Write();    
  etaFakeMc4      -> Write();    
  ptVsEtaSignalMc -> Write();
  ptVsEtaFakeMc   -> Write();
  mvaSignalEBMc0 -> Write();
  mvaSignalEBMc1 -> Write();
  mvaSignalEBMc2 -> Write();
  mvaSignalEBMc3 -> Write();
  mvaSignalEBMc4 -> Write();
  mvaSignalEEMc0 -> Write();
  mvaSignalEEMc1 -> Write();
  mvaSignalEEMc2 -> Write();
  mvaSignalEEMc3 -> Write();
  mvaSignalEEMc4 -> Write();
  mvaFakeEBMc0   -> Write();
  mvaFakeEBMc1   -> Write();
  mvaFakeEBMc2   -> Write();
  mvaFakeEBMc3   -> Write();
  mvaFakeEBMc4   -> Write();
  mvaFakeEEMc0   -> Write();
  mvaFakeEEMc1   -> Write();
  mvaFakeEEMc2   -> Write();
  mvaFakeEEMc3   -> Write();
  mvaFakeEEMc4   -> Write();
  mvaSignalMc_LptPfOverlap -> Write();
  ptSignalMc_LptPfOverlap  -> Write();
  etaSignalMc_LptPfOverlap -> Write();
  mvaSignalMc_LptNotPfOverlap -> Write();
  ptSignalMc_LptNotPfOverlap  -> Write();
  etaSignalMc_LptNotPfOverlap -> Write();
  mvaFakeMc_LptPfOverlap -> Write();
  ptFakeMc_LptPfOverlap  -> Write();
  etaFakeMc_LptPfOverlap -> Write();
  mvaFakeMc_LptNotPfOverlap -> Write();
  ptFakeMc_LptNotPfOverlap  -> Write();
  etaFakeMc_LptNotPfOverlap -> Write();
  //
  myfile.Close();


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

  // -------------------------------------------------------------------------
  // Plots
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);

  TLegend *leg;
  leg = new TLegend(0.10,0.65,0.45,0.90);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.05);
  leg->SetFillColor(0);
  leg->AddEntry(mvaSignalMc, "Matching truth", "lp");
  leg->AddEntry(mvaFakeMc,   "Not matching truth", "lp");
  //
  TLegend *legB;
  legB = new TLegend(0.40,0.65,0.75,0.90);
  legB->SetFillStyle(0);
  legB->SetBorderSize(0);
  legB->SetTextSize(0.05);
  legB->SetFillColor(0);
  legB->AddEntry(mvaSignalMc, "Matching truth", "lp");
  legB->AddEntry(mvaFakeMc,   "Not matching truth", "lp");
  //
  TLegend *legC;
  legC = new TLegend(0.25,0.15,0.65,0.30);
  legC->SetFillStyle(0);
  legC->SetBorderSize(0);
  legC->SetTextSize(0.05);
  legC->SetFillColor(0);
  legC->AddEntry(mvaSignalMc, "Matching truth", "lp");
  legC->AddEntry(mvaFakeMc,   "Not matching truth", "lp");

  TCanvas cmvamc("cmvamc","cmvamc",1);
  mvaSignalMc->SetTitle("");
  mvaFakeMc->SetTitle("");
  mvaSignalMc->GetXaxis()->SetTitle("Id BDT");
  mvaFakeMc->GetXaxis()->SetTitle("Id BDT");
  mvaSignalMc->DrawNormalized("hist");
  mvaFakeMc->DrawNormalized("samehist");
  leg->Draw();
  cmvamc.SaveAs("outputBDT_matchVsFake_noTnP.png");

  TCanvas cptmc("cptmc","cptmc",1);
  ptSignalMc->SetTitle("");
  ptFakeMc->SetTitle("");
  ptSignalMc->GetXaxis()->SetTitle("pT [GeV]");
  ptFakeMc->GetXaxis()->SetTitle("pT [GeV]");
  ptFakeMc->DrawNormalized("hist");
  ptSignalMc->DrawNormalized("samehist");
  legB->Draw();
  cptmc.SaveAs("pt_matchVsFake_noTnP.png");

  TCanvas cetamc("cetamc","cetamc",1);
  etaSignalMc->SetTitle("");
  etaFakeMc->SetTitle("");
  etaSignalMc->GetXaxis()->SetTitle("#eta");
  etaFakeMc->GetXaxis()->SetTitle("#eta");
  etaSignalMc->DrawNormalized("hist");
  etaFakeMc->DrawNormalized("samehist");
  legC->Draw();
  cetamc.SaveAs("eta_matchVsFake_noTnP.png");

  TCanvas cmvaeb0("cmvaeb0","cmvaeb0",1);
  mvaSignalEBMc0 -> SetTitle("EB, 0.5 < pT < 1");
  mvaFakeEBMc0   -> SetTitle("EB, 0.5 < pT < 1");
  mvaSignalEBMc0 -> GetXaxis()->SetTitle("Id BDT");
  mvaFakeEBMc0   -> GetXaxis()->SetTitle("Id BDT");
  mvaFakeEBMc0->DrawNormalized("hist");
  mvaSignalEBMc0->DrawNormalized("samehist");
  leg->Draw();
  cmvaeb0.SaveAs("outputBDT_matchVsFake_noTnP_EB0.png");  
  //
  TCanvas cmvaeb1("cmvaeb1","cmvaeb1",1);
  mvaSignalEBMc1 -> SetTitle("EB, 1.0 < pT < 1.5");
  mvaFakeEBMc1   -> SetTitle("EB, 1.0 < pT < 1.5");
  mvaSignalEBMc1 -> GetXaxis()->SetTitle("Id BDT");
  mvaFakeEBMc1   -> GetXaxis()->SetTitle("Id BDT");
  mvaFakeEBMc1->DrawNormalized("hist");
  mvaSignalEBMc1->DrawNormalized("samehist");
  leg->Draw();
  cmvaeb1.SaveAs("outputBDT_matchVsFake_noTnP_EB1.png");  
  //
  TCanvas cmvaeb2("cmvaeb2","cmvaeb2",1);
  mvaSignalEBMc2 -> SetTitle("EB, 1.5 < pT < 2.0");
  mvaFakeEBMc2   -> SetTitle("EB, 1.5 < pT < 2.0");
  mvaSignalEBMc2 -> GetXaxis()->SetTitle("Id BDT");
  mvaFakeEBMc2   -> GetXaxis()->SetTitle("Id BDT");
  mvaSignalEBMc2->DrawNormalized("hist");
  mvaFakeEBMc2->DrawNormalized("samehist");
  leg->Draw();
  cmvaeb2.SaveAs("outputBDT_matchVsFake_noTnP_EB2.png");  
  //
  TCanvas cmvaeb3("cmvaeb3","cmvaeb3",1);
  mvaSignalEBMc3 -> SetTitle("EB, 2.0 < pT < 5.0");
  mvaFakeEBMc3   -> SetTitle("EB, 2.0 < pT < 5.0");
  mvaSignalEBMc3 -> GetXaxis()->SetTitle("Id BDT");
  mvaFakeEBMc3   -> GetXaxis()->SetTitle("Id BDT");
  mvaSignalEBMc3->DrawNormalized("hist");
  mvaFakeEBMc3->DrawNormalized("samehist");
  leg->Draw();
  cmvaeb3.SaveAs("outputBDT_matchVsFake_noTnP_EB3.png");  
  //
  TCanvas cmvaeb4("cmvaeb4","cmvaeb4",1);
  mvaSignalEBMc4 -> SetTitle("EB, pT > 5");
  mvaFakeEBMc4   -> SetTitle("EB, pT > 5");
  mvaSignalEBMc4 -> GetXaxis()->SetTitle("Id BDT");
  mvaFakeEBMc4   -> GetXaxis()->SetTitle("Id BDT");
  mvaSignalEBMc4->DrawNormalized("hist");
  mvaFakeEBMc4->DrawNormalized("samehist");
  leg->Draw();
  cmvaeb4.SaveAs("outputBDT_matchVsFake_noTnP_EB4.png");  

  TCanvas cmvaee0("cmvaee0","cmvaee0",1);
  mvaSignalEEMc0 -> SetTitle("EE, 0.5 < pT < 1");
  mvaFakeEEMc0   -> SetTitle("EE, 0.5 < pT < 1");
  mvaSignalEEMc0 -> GetXaxis()->SetTitle("Id BDT");
  mvaFakeEEMc0   -> GetXaxis()->SetTitle("Id BDT");
  mvaFakeEEMc0->DrawNormalized("hist");
  mvaSignalEEMc0->DrawNormalized("samehist");
  leg->Draw();
  cmvaee0.SaveAs("outputBDT_matchVsFake_noTnP_EE0.png");  
  //
  TCanvas cmvaee1("cmvaee1","cmvaee1",1);
  mvaSignalEEMc1 -> SetTitle("EE, 1.0 < pT < 1.5");
  mvaFakeEEMc1   -> SetTitle("EE, 1.0 < pT < 1.5");
  mvaSignalEEMc1 -> GetXaxis()->SetTitle("Id BDT");
  mvaFakeEEMc1   -> GetXaxis()->SetTitle("Id BDT");
  mvaFakeEEMc1->DrawNormalized("hist");
  mvaSignalEEMc1->DrawNormalized("samehist");
  leg->Draw();
  cmvaee1.SaveAs("outputBDT_matchVsFake_noTnP_EE1.png");  
  //
  TCanvas cmvaee2("cmvaee2","cmvaee2",1);
  mvaSignalEEMc2 -> SetTitle("EE, 1.5 < pT < 2.0");
  mvaFakeEEMc2   -> SetTitle("EE, 1.5 < pT < 2.0");
  mvaSignalEEMc2 -> GetXaxis()->SetTitle("Id BDT");
  mvaFakeEEMc2   -> GetXaxis()->SetTitle("Id BDT");
  mvaSignalEEMc2->DrawNormalized("hist");
  mvaFakeEEMc2->DrawNormalized("samehist");
  leg->Draw();
  cmvaee2.SaveAs("outputBDT_matchVsFake_noTnP_EE2.png");  
  //
  TCanvas cmvaee3("cmvaee3","cmvaee3",1);
  mvaSignalEEMc3 -> SetTitle("EE, 2.0 < pT < 5.0");
  mvaFakeEEMc3   -> SetTitle("EE, 2.0 < pT < 5.0");
  mvaSignalEEMc3 -> GetXaxis()->SetTitle("Id BDT");
  mvaFakeEEMc3   -> GetXaxis()->SetTitle("Id BDT");
  mvaSignalEEMc3->DrawNormalized("hist");
  mvaFakeEEMc3->DrawNormalized("samehist");
  leg->Draw();
  cmvaee3.SaveAs("outputBDT_matchVsFake_noTnP_EE3.png");  
  //
  TCanvas cmvaee4("cmvaee4","cmvaee4",1);
  mvaSignalEEMc4 -> SetTitle("EE, pT > 5");
  mvaFakeEEMc4   -> SetTitle("EE, pT > 5");
  mvaSignalEEMc4 -> GetXaxis()->SetTitle("Id BDT");
  mvaFakeEEMc4   -> GetXaxis()->SetTitle("Id BDT");
  mvaSignalEEMc4->DrawNormalized("hist");
  mvaFakeEEMc4->DrawNormalized("samehist");
  leg->Draw();
  cmvaee4.SaveAs("outputBDT_matchVsFake_noTnP_EE4.png");  
  //
  //
  if (studyOverlap==1) {
    
    TLegend *legA1;
    legA1 = new TLegend(0.10,0.65,0.45,0.90);
    legA1->SetFillStyle(0);
    legA1->SetBorderSize(0);
    legA1->SetTextSize(0.05);
    legA1->SetFillColor(0);
    legA1->AddEntry(mvaFakeMc_LptNotPfOverlap, "No PF overlap", "lp");
    legA1->AddEntry(mvaFakeMc_LptPfOverlap,    "PF overlap", "lp");
    //
    TLegend *legB1;
    legB1 = new TLegend(0.40,0.65,0.75,0.90);
    legB1->SetFillStyle(0);
    legB1->SetBorderSize(0);
    legB1->SetTextSize(0.05);
    legB1->SetFillColor(0);
    legB1->AddEntry(mvaFakeMc_LptNotPfOverlap, "No PF overlap", "lp");
    legB1->AddEntry(mvaFakeMc_LptPfOverlap,    "PF overlap", "lp");
    //
    TLegend *legC1;
    legC1 = new TLegend(0.25,0.15,0.65,0.30);
    legC1->SetFillStyle(0);
    legC1->SetBorderSize(0);
    legC1->SetTextSize(0.05);
    legC1->SetFillColor(0);
    legC1->AddEntry(mvaFakeMc_LptNotPfOverlap, "No PF overlap", "lp");
    legC1->AddEntry(mvaFakeMc_LptPfOverlap,    "PF overlap", "lp");

    TCanvas ca1("ca1","ca1",1);
    mvaFakeMc_LptNotPfOverlap -> SetTitle(""); 
    mvaFakeMc_LptPfOverlap    -> SetTitle(""); 
    mvaFakeMc_LptNotPfOverlap -> GetXaxis()->SetTitle("Id BDT");
    mvaFakeMc_LptPfOverlap    -> GetXaxis()->SetTitle("Id BDT");
    mvaFakeMc_LptNotPfOverlap -> DrawNormalized("hist");
    mvaFakeMc_LptPfOverlap    -> DrawNormalized("samehist");
    legA1->Draw();
    ca1.SaveAs("outputBDT_fakeMC_w-wo-overlap.png");
    //
    TCanvas ca1b("ca1b","ca1b",1);
    mvaSignalMc_LptNotPfOverlap -> SetTitle(""); 
    mvaSignalMc_LptPfOverlap    -> SetTitle(""); 
    mvaSignalMc_LptNotPfOverlap -> GetXaxis()->SetTitle("Id BDT");
    mvaSignalMc_LptPfOverlap    -> GetXaxis()->SetTitle("Id BDT");
    mvaSignalMc_LptPfOverlap    -> DrawNormalized("hist");
    mvaSignalMc_LptNotPfOverlap -> DrawNormalized("samehist");
    legA1->Draw();
    ca1b.SaveAs("outputBDT_matchMC_w-wo-overlap.png");

    TCanvas ca2("ca2","ca2",1);
    ptFakeMc_LptNotPfOverlap -> SetTitle(""); 
    ptFakeMc_LptPfOverlap    -> SetTitle(""); 
    ptFakeMc_LptNotPfOverlap -> GetXaxis()->SetTitle("pt [GeV]");
    ptFakeMc_LptPfOverlap    -> GetXaxis()->SetTitle("pt [GeV]");
    ptFakeMc_LptNotPfOverlap -> DrawNormalized("hist");
    ptFakeMc_LptPfOverlap    -> DrawNormalized("samehist");
    legB1->Draw();
    ca2.SaveAs("pt_fakeMC_w-wo-overlap.png");
    //
    TCanvas ca2b("ca2b","ca2b",1);
    ptSignalMc_LptNotPfOverlap -> SetTitle(""); 
    ptSignalMc_LptPfOverlap    -> SetTitle(""); 
    ptSignalMc_LptNotPfOverlap -> GetXaxis()->SetTitle("pt [GeV]");
    ptSignalMc_LptPfOverlap    -> GetXaxis()->SetTitle("pt [GeV]");
    ptSignalMc_LptNotPfOverlap -> DrawNormalized("hist");
    ptSignalMc_LptPfOverlap    -> DrawNormalized("samehist");
    legB1->Draw();
    ca2b.SaveAs("pt_matchMC_w-wo-overlap.png");

    TCanvas ca3("ca3","ca3",1);
    etaFakeMc_LptNotPfOverlap -> Rebin();
    etaFakeMc_LptPfOverlap    -> Rebin();
    etaFakeMc_LptNotPfOverlap -> SetTitle(""); 
    etaFakeMc_LptPfOverlap    -> SetTitle(""); 
    etaFakeMc_LptNotPfOverlap -> GetXaxis()->SetTitle("#eta");
    etaFakeMc_LptPfOverlap    -> GetXaxis()->SetTitle("#eta");
    etaFakeMc_LptNotPfOverlap -> DrawNormalized("hist");
    etaFakeMc_LptPfOverlap    -> DrawNormalized("samehist");
    legC1->Draw();
    ca3.SaveAs("eta_fakeMC_w-wo-overlap.png");
    //
    TCanvas ca3b("ca3b","ca3b",1);
    etaSignalMc_LptNotPfOverlap -> Rebin();
    etaSignalMc_LptPfOverlap    -> Rebin();
    etaSignalMc_LptNotPfOverlap -> SetTitle(""); 
    etaSignalMc_LptPfOverlap    -> SetTitle(""); 
    etaSignalMc_LptNotPfOverlap -> GetXaxis()->SetTitle("#eta");
    etaSignalMc_LptPfOverlap    -> GetXaxis()->SetTitle("#eta");
    etaSignalMc_LptPfOverlap    -> DrawNormalized("hist");
    etaSignalMc_LptNotPfOverlap -> DrawNormalized("samehist");
    legC1->Draw();
    ca3b.SaveAs("eta_matchMC_w-wo-overlap.png");
  }  


}
