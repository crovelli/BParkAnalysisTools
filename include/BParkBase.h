//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Feb  3 14:02:00 2020 by ROOT version 6.18/04
// from TTree Events/Events
// found on file: /eos/cms/store/group/cmst3/group/bpark/BParkingNANO_2020Jan16/BuToKJpsi_Toee_Mufilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/crab_BuToKJpsi_Toee/200116_215618/0000/BParkNANO_mc_2020Jan16_1.root
//////////////////////////////////////////////////////////

#ifndef BParkBase_h
#define BParkBase_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class BParkBase {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   UInt_t          run;
   UInt_t          luminosityBlock;
   ULong64_t       event;
   UInt_t          nBToKEE;
   Float_t         BToKEE_b_iso03[700];   //[nBToKEE]
   Float_t         BToKEE_b_iso04[700];   //[nBToKEE]
   Float_t         BToKEE_cos2D[700];   //[nBToKEE]
   Float_t         BToKEE_eta[700];   //[nBToKEE]
   Float_t         BToKEE_fit_cos2D[700];   //[nBToKEE]
   Float_t         BToKEE_fit_eta[700];   //[nBToKEE]
   Float_t         BToKEE_fit_k_eta[700];   //[nBToKEE]
   Float_t         BToKEE_fit_k_phi[700];   //[nBToKEE]
   Float_t         BToKEE_fit_k_pt[700];   //[nBToKEE]
   Float_t         BToKEE_fit_l1_eta[700];   //[nBToKEE]
   Float_t         BToKEE_fit_l1_phi[700];   //[nBToKEE]
   Float_t         BToKEE_fit_l1_pt[700];   //[nBToKEE]
   Float_t         BToKEE_fit_l2_eta[700];   //[nBToKEE]
   Float_t         BToKEE_fit_l2_phi[700];   //[nBToKEE]
   Float_t         BToKEE_fit_l2_pt[700];   //[nBToKEE]
   Float_t         BToKEE_fit_mass[700];   //[nBToKEE]
   Float_t         BToKEE_fit_massErr[700];   //[nBToKEE]
   Float_t         BToKEE_fit_phi[700];   //[nBToKEE]
   Float_t         BToKEE_fit_pt[700];   //[nBToKEE]
   Float_t         BToKEE_k_iso03[700];   //[nBToKEE]
   Float_t         BToKEE_k_iso04[700];   //[nBToKEE]
   Float_t         BToKEE_l1_iso03[700];   //[nBToKEE]
   Float_t         BToKEE_l1_iso04[700];   //[nBToKEE]
   Float_t         BToKEE_l2_iso03[700];   //[nBToKEE]
   Float_t         BToKEE_l2_iso04[700];   //[nBToKEE]
   Float_t         BToKEE_l_xy[700];   //[nBToKEE]
   Float_t         BToKEE_l_xy_unc[700];   //[nBToKEE]
   Float_t         BToKEE_mass[700];   //[nBToKEE]
   Float_t         BToKEE_maxDR[700];   //[nBToKEE]
   Float_t         BToKEE_minDR[700];   //[nBToKEE]
   Float_t         BToKEE_mllErr_llfit[700];   //[nBToKEE]
   Float_t         BToKEE_mll_fullfit[700];   //[nBToKEE]
   Float_t         BToKEE_mll_llfit[700];   //[nBToKEE]
   Float_t         BToKEE_mll_raw[700];   //[nBToKEE]
   Float_t         BToKEE_phi[700];   //[nBToKEE]
   Float_t         BToKEE_pt[700];   //[nBToKEE]
   Float_t         BToKEE_svprob[700];   //[nBToKEE]
   Float_t         BToKEE_vtx_ex[700];   //[nBToKEE]
   Float_t         BToKEE_vtx_ey[700];   //[nBToKEE]
   Float_t         BToKEE_vtx_ez[700];   //[nBToKEE]
   Float_t         BToKEE_vtx_x[700];   //[nBToKEE]
   Float_t         BToKEE_vtx_y[700];   //[nBToKEE]
   Float_t         BToKEE_vtx_z[700];   //[nBToKEE]
   Int_t           BToKEE_charge[700];   //[nBToKEE]
   Int_t           BToKEE_kIdx[700];   //[nBToKEE]
   Int_t           BToKEE_l1Idx[700];   //[nBToKEE]
   Int_t           BToKEE_l2Idx[700];   //[nBToKEE]
   Int_t           BToKEE_n_k_used[700];   //[nBToKEE]
   Int_t           BToKEE_n_l1_used[700];   //[nBToKEE]
   Int_t           BToKEE_n_l2_used[700];   //[nBToKEE]
   Int_t           BToKEE_pdgId[700];   //[nBToKEE]
   UInt_t          nBToKMuMu;
   Float_t         BToKMuMu_b_iso03[100];   //[nBToKMuMu]
   Float_t         BToKMuMu_b_iso04[100];   //[nBToKMuMu]
   Float_t         BToKMuMu_cos2D[100];   //[nBToKMuMu]
   Float_t         BToKMuMu_eta[100];   //[nBToKMuMu]
   Float_t         BToKMuMu_fit_cos2D[100];   //[nBToKMuMu]
   Float_t         BToKMuMu_fit_eta[100];   //[nBToKMuMu]
   Float_t         BToKMuMu_fit_k_eta[100];   //[nBToKMuMu]
   Float_t         BToKMuMu_fit_k_phi[100];   //[nBToKMuMu]
   Float_t         BToKMuMu_fit_k_pt[100];   //[nBToKMuMu]
   Float_t         BToKMuMu_fit_l1_eta[100];   //[nBToKMuMu]
   Float_t         BToKMuMu_fit_l1_phi[100];   //[nBToKMuMu]
   Float_t         BToKMuMu_fit_l1_pt[100];   //[nBToKMuMu]
   Float_t         BToKMuMu_fit_l2_eta[100];   //[nBToKMuMu]
   Float_t         BToKMuMu_fit_l2_phi[100];   //[nBToKMuMu]
   Float_t         BToKMuMu_fit_l2_pt[100];   //[nBToKMuMu]
   Float_t         BToKMuMu_fit_mass[100];   //[nBToKMuMu]
   Float_t         BToKMuMu_fit_massErr[100];   //[nBToKMuMu]
   Float_t         BToKMuMu_fit_phi[100];   //[nBToKMuMu]
   Float_t         BToKMuMu_fit_pt[100];   //[nBToKMuMu]
   Float_t         BToKMuMu_k_iso03[100];   //[nBToKMuMu]
   Float_t         BToKMuMu_k_iso04[100];   //[nBToKMuMu]
   Float_t         BToKMuMu_l1_iso03[100];   //[nBToKMuMu]
   Float_t         BToKMuMu_l1_iso04[100];   //[nBToKMuMu]
   Float_t         BToKMuMu_l2_iso03[100];   //[nBToKMuMu]
   Float_t         BToKMuMu_l2_iso04[100];   //[nBToKMuMu]
   Float_t         BToKMuMu_l_xy[100];   //[nBToKMuMu]
   Float_t         BToKMuMu_l_xy_unc[100];   //[nBToKMuMu]
   Float_t         BToKMuMu_mass[100];   //[nBToKMuMu]
   Float_t         BToKMuMu_maxDR[100];   //[nBToKMuMu]
   Float_t         BToKMuMu_minDR[100];   //[nBToKMuMu]
   Float_t         BToKMuMu_mllErr_llfit[100];   //[nBToKMuMu]
   Float_t         BToKMuMu_mll_fullfit[100];   //[nBToKMuMu]
   Float_t         BToKMuMu_mll_llfit[100];   //[nBToKMuMu]
   Float_t         BToKMuMu_mll_raw[100];   //[nBToKMuMu]
   Float_t         BToKMuMu_phi[100];   //[nBToKMuMu]
   Float_t         BToKMuMu_pt[100];   //[nBToKMuMu]
   Float_t         BToKMuMu_svprob[100];   //[nBToKMuMu]
   Float_t         BToKMuMu_vtx_ex[100];   //[nBToKMuMu]
   Float_t         BToKMuMu_vtx_ey[100];   //[nBToKMuMu]
   Float_t         BToKMuMu_vtx_ez[100];   //[nBToKMuMu]
   Float_t         BToKMuMu_vtx_x[100];   //[nBToKMuMu]
   Float_t         BToKMuMu_vtx_y[100];   //[nBToKMuMu]
   Float_t         BToKMuMu_vtx_z[100];   //[nBToKMuMu]
   Int_t           BToKMuMu_charge[100];   //[nBToKMuMu]
   Int_t           BToKMuMu_kIdx[100];   //[nBToKMuMu]
   Int_t           BToKMuMu_l1Idx[100];   //[nBToKMuMu]
   Int_t           BToKMuMu_l2Idx[100];   //[nBToKMuMu]
   Int_t           BToKMuMu_n_k_used[100];   //[nBToKMuMu]
   Int_t           BToKMuMu_n_l1_used[100];   //[nBToKMuMu]
   Int_t           BToKMuMu_n_l2_used[100];   //[nBToKMuMu]
   Int_t           BToKMuMu_pdgId[100];   //[nBToKMuMu]
   UInt_t          nElectron;
   Float_t         Electron_dxy[500];   //[nElectron]
   Float_t         Electron_dxyErr[500];   //[nElectron]
   Float_t         Electron_dz[500];   //[nElectron]
   Float_t         Electron_dzErr[500];   //[nElectron]
   Float_t         Electron_eta[500];   //[nElectron]
   Float_t         Electron_fBrem[500];   //[nElectron]
   Float_t         Electron_ip3d[500];   //[nElectron]
   Float_t         Electron_mass[500];   //[nElectron]
   Float_t         Electron_mvaId[500];   //[nElectron]
   Float_t         Electron_pfRelIso[500];   //[nElectron]
   Float_t         Electron_pfmvaId[500];   //[nElectron]
   Float_t         Electron_phi[500];   //[nElectron]
   Float_t         Electron_pt[500];   //[nElectron]
   Float_t         Electron_ptBiased[500];   //[nElectron]
   Float_t         Electron_sip3d[500];   //[nElectron]
   Float_t         Electron_trkRelIso[500];   //[nElectron]
   Float_t         Electron_unBiased[500];   //[nElectron]
   Float_t         Electron_vx[500];   //[nElectron]
   Float_t         Electron_vy[500];   //[nElectron]
   Float_t         Electron_vz[500];   //[nElectron]
   Int_t           Electron_charge[500];   //[nElectron]
   Int_t           Electron_pdgId[500];   //[nElectron]
   Bool_t          Electron_convVeto[500];   //[nElectron]
   Bool_t          Electron_isLowPt[500];   //[nElectron]
   Bool_t          Electron_isPF[500];   //[nElectron]
   Bool_t          Electron_isPFoverlap[500];   //[nElectron]
   UInt_t          nGenPart;
   Float_t         GenPart_eta[600];   //[nGenPart]
   Float_t         GenPart_mass[600];   //[nGenPart]
   Float_t         GenPart_phi[600];   //[nGenPart]
   Float_t         GenPart_pt[600];   //[nGenPart]
   Float_t         GenPart_vx[600];   //[nGenPart]
   Float_t         GenPart_vy[600];   //[nGenPart]
   Float_t         GenPart_vz[600];   //[nGenPart]
   Int_t           GenPart_genPartIdxMother[600];   //[nGenPart]
   Int_t           GenPart_pdgId[600];   //[nGenPart]
   Int_t           GenPart_status[600];   //[nGenPart]
   Int_t           GenPart_statusFlags[600];   //[nGenPart]
   Float_t         Generator_binvar;
   Float_t         Generator_scalePDF;
   Float_t         Generator_weight;
   Float_t         Generator_x1;
   Float_t         Generator_x2;
   Float_t         Generator_xpdf1;
   Float_t         Generator_xpdf2;
   Int_t           Generator_id1;
   Int_t           Generator_id2;
   Float_t         genWeight;
   UInt_t          nPSWeight;
   Float_t         PSWeight[1];   //[nPSWeight]
   UInt_t          nMuon;
   Float_t         Muon_dxy[100];   //[nMuon]
   Float_t         Muon_dxyErr[100];   //[nMuon]
   Float_t         Muon_dz[100];   //[nMuon]
   Float_t         Muon_dzErr[100];   //[nMuon]
   Float_t         Muon_eta[100];   //[nMuon]
   Float_t         Muon_ip3d[100];   //[nMuon]
   Float_t         Muon_mass[100];   //[nMuon]
   Float_t         Muon_pfRelIso03_all[100];   //[nMuon]
   Float_t         Muon_pfRelIso04_all[100];   //[nMuon]
   Float_t         Muon_phi[100];   //[nMuon]
   Float_t         Muon_pt[100];   //[nMuon]
   Float_t         Muon_ptErr[100];   //[nMuon]
   Float_t         Muon_sip3d[100];   //[nMuon]
   Float_t         Muon_vx[100];   //[nMuon]
   Float_t         Muon_vy[100];   //[nMuon]
   Float_t         Muon_vz[100];   //[nMuon]
   Int_t           Muon_charge[100];   //[nMuon]
   Int_t           Muon_isTriggering[100];   //[nMuon]
   Int_t           Muon_pdgId[100];   //[nMuon]
   Bool_t          Muon_isGlobal[100];   //[nMuon]
   Bool_t          Muon_isPFcand[100];   //[nMuon]
   Bool_t          Muon_isTracker[100];   //[nMuon]
   Bool_t          Muon_mediumId[100];   //[nMuon]
   UChar_t         Muon_pfIsoId[100];   //[nMuon]
   Bool_t          Muon_softId[100];   //[nMuon]
   Bool_t          Muon_tightId[100];   //[nMuon]
   UChar_t         Muon_tkIsoId[100];   //[nMuon]
   Bool_t          Muon_triggerIdLoose[100];   //[nMuon]
   UInt_t          nTriggerMuon;
   Float_t         TriggerMuon_eta[20];   //[nTriggerMuon]
   Float_t         TriggerMuon_mass[20];   //[nTriggerMuon]
   Float_t         TriggerMuon_phi[20];   //[nTriggerMuon]
   Float_t         TriggerMuon_pt[20];   //[nTriggerMuon]
   Float_t         TriggerMuon_vx[20];   //[nTriggerMuon]
   Float_t         TriggerMuon_vy[20];   //[nTriggerMuon]
   Float_t         TriggerMuon_vz[20];   //[nTriggerMuon]
   Int_t           TriggerMuon_charge[20];   //[nTriggerMuon]
   Int_t           TriggerMuon_pdgId[20];   //[nTriggerMuon]
   Int_t           TriggerMuon_trgMuonIndex[20];   //[nTriggerMuon]
   Float_t         Pileup_nTrueInt;
   Float_t         Pileup_pudensity;
   Float_t         Pileup_gpudensity;
   Int_t           Pileup_nPU;
   Int_t           Pileup_sumEOOT;
   Int_t           Pileup_sumLOOT;
   Float_t         fixedGridRhoFastjetAll;
   Float_t         fixedGridRhoFastjetCentral;
   Float_t         fixedGridRhoFastjetCentralCalo;
   Float_t         fixedGridRhoFastjetCentralChargedPileUp;
   Float_t         fixedGridRhoFastjetCentralNeutral;
   UInt_t          nProbeTracks;
   Float_t         ProbeTracks_DCASig[1000];   //[nProbeTracks]
   Float_t         ProbeTracks_dxy[1000];   //[nProbeTracks]
   Float_t         ProbeTracks_dxyS[1000];   //[nProbeTracks]
   Float_t         ProbeTracks_dz[1000];   //[nProbeTracks]
   Float_t         ProbeTracks_dzS[1000];   //[nProbeTracks]
   Float_t         ProbeTracks_eta[1000];   //[nProbeTracks]
   Float_t         ProbeTracks_mass[1000];   //[nProbeTracks]
   Float_t         ProbeTracks_phi[1000];   //[nProbeTracks]
   Float_t         ProbeTracks_pt[1000];   //[nProbeTracks]
   Float_t         ProbeTracks_vx[1000];   //[nProbeTracks]
   Float_t         ProbeTracks_vy[1000];   //[nProbeTracks]
   Float_t         ProbeTracks_vz[1000];   //[nProbeTracks]
   Int_t           ProbeTracks_charge[1000];   //[nProbeTracks]
   Int_t           ProbeTracks_isLostTrk[1000];   //[nProbeTracks]
   Int_t           ProbeTracks_isPacked[1000];   //[nProbeTracks]
   Int_t           ProbeTracks_nValidHits[1000];   //[nProbeTracks]
   Int_t           ProbeTracks_pdgId[1000];   //[nProbeTracks]
   Bool_t          ProbeTracks_isMatchedToEle[1000];   //[nProbeTracks]
   Bool_t          ProbeTracks_isMatchedToLooseMuon[1000];   //[nProbeTracks]
   Bool_t          ProbeTracks_isMatchedToMediumMuon[1000];   //[nProbeTracks]
   Bool_t          ProbeTracks_isMatchedToMuon[1000];   //[nProbeTracks]
   Bool_t          ProbeTracks_isMatchedToSoftMuon[1000];   //[nProbeTracks]
   UChar_t         HLT_Mu7_IP4;
   UChar_t         HLT_Mu8_IP6;
   UChar_t         HLT_Mu8_IP5;
   UChar_t         HLT_Mu8_IP3;
   UChar_t         HLT_Mu8p5_IP3p5;
   UChar_t         HLT_Mu9_IP6;
   UChar_t         HLT_Mu9_IP5;
   UChar_t         HLT_Mu9_IP4;
   UChar_t         HLT_Mu10p5_IP3p5;
   UChar_t         HLT_Mu12_IP6;
   UChar_t         L1_SingleMu7er1p5;
   UChar_t         L1_SingleMu8er1p5;
   UChar_t         L1_SingleMu9er1p5;
   UChar_t         L1_SingleMu10er1p5;
   UChar_t         L1_SingleMu12er1p5;
   UChar_t         L1_SingleMu22;
   UInt_t          nTrigObj;
   Float_t         TrigObj_pt[20];   //[nTrigObj]
   Float_t         TrigObj_eta[20];   //[nTrigObj]
   Float_t         TrigObj_phi[20];   //[nTrigObj]
   Float_t         TrigObj_l1pt[20];   //[nTrigObj]
   Float_t         TrigObj_l1pt_2[20];   //[nTrigObj]
   Float_t         TrigObj_l2pt[20];   //[nTrigObj]
   Int_t           TrigObj_id[20];   //[nTrigObj]
   Int_t           TrigObj_l1iso[20];   //[nTrigObj]
   Int_t           TrigObj_l1charge[20];   //[nTrigObj]
   Int_t           TrigObj_filterBits[20];   //[nTrigObj]
   UInt_t          nOtherPV;
   Float_t         OtherPV_z[3];   //[nOtherPV]
   Float_t         PV_ndof;
   Float_t         PV_x;
   Float_t         PV_y;
   Float_t         PV_z;
   Float_t         PV_chi2;
   Float_t         PV_score;
   Int_t           PV_npvs;
   Int_t           PV_npvsGood;
   UInt_t          nSV;
   Float_t         SV_dlen[30];   //[nSV]
   Float_t         SV_dlenSig[30];   //[nSV]
   Float_t         SV_pAngle[30];   //[nSV]
   Int_t           Electron_genPartIdx[500];   //[nElectron]
   Int_t           Electron_genPartFlav[500];   //[nElectron]
   Int_t           Muon_genPartIdx[30];   //[nMuon]
   Int_t           Muon_genPartFlav[30];   //[nMuon]
   Float_t         SV_chi2[30];   //[nSV]
   Float_t         SV_eta[30];   //[nSV]
   Float_t         SV_mass[30];   //[nSV]
   Float_t         SV_ndof[30];   //[nSV]
   Float_t         SV_phi[30];   //[nSV]
   Float_t         SV_pt[30];   //[nSV]
   Float_t         SV_x[30];   //[nSV]
   Float_t         SV_y[30];   //[nSV]
   Float_t         SV_z[30];   //[nSV]
   Int_t           ProbeTracks_genPartIdx[1000];   //[nProbeTracks]
   Int_t           ProbeTracks_genPartFlav[1000];   //[nProbeTracks]

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_luminosityBlock;   //!
   TBranch        *b_event;   //!
   TBranch        *b_nBToKEE;   //!
   TBranch        *b_BToKEE_b_iso03;   //!
   TBranch        *b_BToKEE_b_iso04;   //!
   TBranch        *b_BToKEE_cos2D;   //!
   TBranch        *b_BToKEE_eta;   //!
   TBranch        *b_BToKEE_fit_cos2D;   //!
   TBranch        *b_BToKEE_fit_eta;   //!
   TBranch        *b_BToKEE_fit_k_eta;   //!
   TBranch        *b_BToKEE_fit_k_phi;   //!
   TBranch        *b_BToKEE_fit_k_pt;   //!
   TBranch        *b_BToKEE_fit_l1_eta;   //!
   TBranch        *b_BToKEE_fit_l1_phi;   //!
   TBranch        *b_BToKEE_fit_l1_pt;   //!
   TBranch        *b_BToKEE_fit_l2_eta;   //!
   TBranch        *b_BToKEE_fit_l2_phi;   //!
   TBranch        *b_BToKEE_fit_l2_pt;   //!
   TBranch        *b_BToKEE_fit_mass;   //!
   TBranch        *b_BToKEE_fit_massErr;   //!
   TBranch        *b_BToKEE_fit_phi;   //!
   TBranch        *b_BToKEE_fit_pt;   //!
   TBranch        *b_BToKEE_k_iso03;   //!
   TBranch        *b_BToKEE_k_iso04;   //!
   TBranch        *b_BToKEE_l1_iso03;   //!
   TBranch        *b_BToKEE_l1_iso04;   //!
   TBranch        *b_BToKEE_l2_iso03;   //!
   TBranch        *b_BToKEE_l2_iso04;   //!
   TBranch        *b_BToKEE_l_xy;   //!
   TBranch        *b_BToKEE_l_xy_unc;   //!
   TBranch        *b_BToKEE_mass;   //!
   TBranch        *b_BToKEE_maxDR;   //!
   TBranch        *b_BToKEE_minDR;   //!
   TBranch        *b_BToKEE_mllErr_llfit;   //!
   TBranch        *b_BToKEE_mll_fullfit;   //!
   TBranch        *b_BToKEE_mll_llfit;   //!
   TBranch        *b_BToKEE_mll_raw;   //!
   TBranch        *b_BToKEE_phi;   //!
   TBranch        *b_BToKEE_pt;   //!
   TBranch        *b_BToKEE_svprob;   //!
   TBranch        *b_BToKEE_vtx_ex;   //!
   TBranch        *b_BToKEE_vtx_ey;   //!
   TBranch        *b_BToKEE_vtx_ez;   //!
   TBranch        *b_BToKEE_vtx_x;   //!
   TBranch        *b_BToKEE_vtx_y;   //!
   TBranch        *b_BToKEE_vtx_z;   //!
   TBranch        *b_BToKEE_charge;   //!
   TBranch        *b_BToKEE_kIdx;   //!
   TBranch        *b_BToKEE_l1Idx;   //!
   TBranch        *b_BToKEE_l2Idx;   //!
   TBranch        *b_BToKEE_n_k_used;   //!
   TBranch        *b_BToKEE_n_l1_used;   //!
   TBranch        *b_BToKEE_n_l2_used;   //!
   TBranch        *b_BToKEE_pdgId;   //!
   TBranch        *b_nBToKMuMu;   //!
   TBranch        *b_BToKMuMu_b_iso03;   //!
   TBranch        *b_BToKMuMu_b_iso04;   //!
   TBranch        *b_BToKMuMu_cos2D;   //!
   TBranch        *b_BToKMuMu_eta;   //!
   TBranch        *b_BToKMuMu_fit_cos2D;   //!
   TBranch        *b_BToKMuMu_fit_eta;   //!
   TBranch        *b_BToKMuMu_fit_k_eta;   //!
   TBranch        *b_BToKMuMu_fit_k_phi;   //!
   TBranch        *b_BToKMuMu_fit_k_pt;   //!
   TBranch        *b_BToKMuMu_fit_l1_eta;   //!
   TBranch        *b_BToKMuMu_fit_l1_phi;   //!
   TBranch        *b_BToKMuMu_fit_l1_pt;   //!
   TBranch        *b_BToKMuMu_fit_l2_eta;   //!
   TBranch        *b_BToKMuMu_fit_l2_phi;   //!
   TBranch        *b_BToKMuMu_fit_l2_pt;   //!
   TBranch        *b_BToKMuMu_fit_mass;   //!
   TBranch        *b_BToKMuMu_fit_massErr;   //!
   TBranch        *b_BToKMuMu_fit_phi;   //!
   TBranch        *b_BToKMuMu_fit_pt;   //!
   TBranch        *b_BToKMuMu_k_iso03;   //!
   TBranch        *b_BToKMuMu_k_iso04;   //!
   TBranch        *b_BToKMuMu_l1_iso03;   //!
   TBranch        *b_BToKMuMu_l1_iso04;   //!
   TBranch        *b_BToKMuMu_l2_iso03;   //!
   TBranch        *b_BToKMuMu_l2_iso04;   //!
   TBranch        *b_BToKMuMu_l_xy;   //!
   TBranch        *b_BToKMuMu_l_xy_unc;   //!
   TBranch        *b_BToKMuMu_mass;   //!
   TBranch        *b_BToKMuMu_maxDR;   //!
   TBranch        *b_BToKMuMu_minDR;   //!
   TBranch        *b_BToKMuMu_mllErr_llfit;   //!
   TBranch        *b_BToKMuMu_mll_fullfit;   //!
   TBranch        *b_BToKMuMu_mll_llfit;   //!
   TBranch        *b_BToKMuMu_mll_raw;   //!
   TBranch        *b_BToKMuMu_phi;   //!
   TBranch        *b_BToKMuMu_pt;   //!
   TBranch        *b_BToKMuMu_svprob;   //!
   TBranch        *b_BToKMuMu_vtx_ex;   //!
   TBranch        *b_BToKMuMu_vtx_ey;   //!
   TBranch        *b_BToKMuMu_vtx_ez;   //!
   TBranch        *b_BToKMuMu_vtx_x;   //!
   TBranch        *b_BToKMuMu_vtx_y;   //!
   TBranch        *b_BToKMuMu_vtx_z;   //!
   TBranch        *b_BToKMuMu_charge;   //!
   TBranch        *b_BToKMuMu_kIdx;   //!
   TBranch        *b_BToKMuMu_l1Idx;   //!
   TBranch        *b_BToKMuMu_l2Idx;   //!
   TBranch        *b_BToKMuMu_n_k_used;   //!
   TBranch        *b_BToKMuMu_n_l1_used;   //!
   TBranch        *b_BToKMuMu_n_l2_used;   //!
   TBranch        *b_BToKMuMu_pdgId;   //!
   TBranch        *b_nElectron;   //!
   TBranch        *b_Electron_dxy;   //!
   TBranch        *b_Electron_dxyErr;   //!
   TBranch        *b_Electron_dz;   //!
   TBranch        *b_Electron_dzErr;   //!
   TBranch        *b_Electron_eta;   //!
   TBranch        *b_Electron_fBrem;   //!
   TBranch        *b_Electron_ip3d;   //!
   TBranch        *b_Electron_mass;   //!
   TBranch        *b_Electron_mvaId;   //!
   TBranch        *b_Electron_pfRelIso;   //!
   TBranch        *b_Electron_pfmvaId;   //!
   TBranch        *b_Electron_phi;   //!
   TBranch        *b_Electron_pt;   //!
   TBranch        *b_Electron_ptBiased;   //!
   TBranch        *b_Electron_sip3d;   //!
   TBranch        *b_Electron_trkRelIso;   //!
   TBranch        *b_Electron_unBiased;   //!
   TBranch        *b_Electron_vx;   //!
   TBranch        *b_Electron_vy;   //!
   TBranch        *b_Electron_vz;   //!
   TBranch        *b_Electron_charge;   //!
   TBranch        *b_Electron_pdgId;   //!
   TBranch        *b_Electron_convVeto;   //!
   TBranch        *b_Electron_isLowPt;   //!
   TBranch        *b_Electron_isPF;   //!
   TBranch        *b_Electron_isPFoverlap;   //!
   TBranch        *b_nGenPart;   //!
   TBranch        *b_GenPart_eta;   //!
   TBranch        *b_GenPart_mass;   //!
   TBranch        *b_GenPart_phi;   //!
   TBranch        *b_GenPart_pt;   //!
   TBranch        *b_GenPart_vx;   //!
   TBranch        *b_GenPart_vy;   //!
   TBranch        *b_GenPart_vz;   //!
   TBranch        *b_GenPart_genPartIdxMother;   //!
   TBranch        *b_GenPart_pdgId;   //!
   TBranch        *b_GenPart_status;   //!
   TBranch        *b_GenPart_statusFlags;   //!
   TBranch        *b_Generator_binvar;   //!
   TBranch        *b_Generator_scalePDF;   //!
   TBranch        *b_Generator_weight;   //!
   TBranch        *b_Generator_x1;   //!
   TBranch        *b_Generator_x2;   //!
   TBranch        *b_Generator_xpdf1;   //!
   TBranch        *b_Generator_xpdf2;   //!
   TBranch        *b_Generator_id1;   //!
   TBranch        *b_Generator_id2;   //!
   TBranch        *b_genWeight;   //!
   TBranch        *b_nPSWeight;   //!
   TBranch        *b_PSWeight;   //!
   TBranch        *b_nMuon;   //!
   TBranch        *b_Muon_dxy;   //!
   TBranch        *b_Muon_dxyErr;   //!
   TBranch        *b_Muon_dz;   //!
   TBranch        *b_Muon_dzErr;   //!
   TBranch        *b_Muon_eta;   //!
   TBranch        *b_Muon_ip3d;   //!
   TBranch        *b_Muon_mass;   //!
   TBranch        *b_Muon_pfRelIso03_all;   //!
   TBranch        *b_Muon_pfRelIso04_all;   //!
   TBranch        *b_Muon_phi;   //!
   TBranch        *b_Muon_pt;   //!
   TBranch        *b_Muon_ptErr;   //!
   TBranch        *b_Muon_sip3d;   //!
   TBranch        *b_Muon_vx;   //!
   TBranch        *b_Muon_vy;   //!
   TBranch        *b_Muon_vz;   //!
   TBranch        *b_Muon_charge;   //!
   TBranch        *b_Muon_isTriggering;   //!
   TBranch        *b_Muon_pdgId;   //!
   TBranch        *b_Muon_isGlobal;   //!
   TBranch        *b_Muon_isPFcand;   //!
   TBranch        *b_Muon_isTracker;   //!
   TBranch        *b_Muon_mediumId;   //!
   TBranch        *b_Muon_pfIsoId;   //!
   TBranch        *b_Muon_softId;   //!
   TBranch        *b_Muon_tightId;   //!
   TBranch        *b_Muon_tkIsoId;   //!
   TBranch        *b_Muon_triggerIdLoose;   //!
   TBranch        *b_nTriggerMuon;   //!
   TBranch        *b_TriggerMuon_eta;   //!
   TBranch        *b_TriggerMuon_mass;   //!
   TBranch        *b_TriggerMuon_phi;   //!
   TBranch        *b_TriggerMuon_pt;   //!
   TBranch        *b_TriggerMuon_vx;   //!
   TBranch        *b_TriggerMuon_vy;   //!
   TBranch        *b_TriggerMuon_vz;   //!
   TBranch        *b_TriggerMuon_charge;   //!
   TBranch        *b_TriggerMuon_pdgId;   //!
   TBranch        *b_TriggerMuon_trgMuonIndex;   //!
   TBranch        *b_Pileup_nTrueInt;   //!
   TBranch        *b_Pileup_pudensity;   //!
   TBranch        *b_Pileup_gpudensity;   //!
   TBranch        *b_Pileup_nPU;   //!
   TBranch        *b_Pileup_sumEOOT;   //!
   TBranch        *b_Pileup_sumLOOT;   //!
   TBranch        *b_fixedGridRhoFastjetAll;   //!
   TBranch        *b_fixedGridRhoFastjetCentral;   //!
   TBranch        *b_fixedGridRhoFastjetCentralCalo;   //!
   TBranch        *b_fixedGridRhoFastjetCentralChargedPileUp;   //!
   TBranch        *b_fixedGridRhoFastjetCentralNeutral;   //!
   TBranch        *b_nProbeTracks;   //!
   TBranch        *b_ProbeTracks_DCASig;   //!
   TBranch        *b_ProbeTracks_dxy;   //!
   TBranch        *b_ProbeTracks_dxyS;   //!
   TBranch        *b_ProbeTracks_dz;   //!
   TBranch        *b_ProbeTracks_dzS;   //!
   TBranch        *b_ProbeTracks_eta;   //!
   TBranch        *b_ProbeTracks_mass;   //!
   TBranch        *b_ProbeTracks_phi;   //!
   TBranch        *b_ProbeTracks_pt;   //!
   TBranch        *b_ProbeTracks_vx;   //!
   TBranch        *b_ProbeTracks_vy;   //!
   TBranch        *b_ProbeTracks_vz;   //!
   TBranch        *b_ProbeTracks_charge;   //!
   TBranch        *b_ProbeTracks_isLostTrk;   //!
   TBranch        *b_ProbeTracks_isPacked;   //!
   TBranch        *b_ProbeTracks_nValidHits;   //!
   TBranch        *b_ProbeTracks_pdgId;   //!
   TBranch        *b_ProbeTracks_isMatchedToEle;   //!
   TBranch        *b_ProbeTracks_isMatchedToLooseMuon;   //!
   TBranch        *b_ProbeTracks_isMatchedToMediumMuon;   //!
   TBranch        *b_ProbeTracks_isMatchedToMuon;   //!
   TBranch        *b_ProbeTracks_isMatchedToSoftMuon;   //!
   TBranch        *b_HLT_Mu7_IP4;   //!
   TBranch        *b_HLT_Mu8_IP6;   //!
   TBranch        *b_HLT_Mu8_IP5;   //!
   TBranch        *b_HLT_Mu8_IP3;   //!
   TBranch        *b_HLT_Mu8p5_IP3p5;   //!
   TBranch        *b_HLT_Mu9_IP6;   //!
   TBranch        *b_HLT_Mu9_IP5;   //!
   TBranch        *b_HLT_Mu9_IP4;   //!
   TBranch        *b_HLT_Mu10p5_IP3p5;   //!
   TBranch        *b_HLT_Mu12_IP6;   //!
   TBranch        *b_L1_SingleMu7er1p5;   //!
   TBranch        *b_L1_SingleMu8er1p5;   //!
   TBranch        *b_L1_SingleMu9er1p5;   //!
   TBranch        *b_L1_SingleMu10er1p5;   //!
   TBranch        *b_L1_SingleMu12er1p5;   //!
   TBranch        *b_L1_SingleMu22;   //!
   TBranch        *b_nTrigObj;   //!
   TBranch        *b_TrigObj_pt;   //!
   TBranch        *b_TrigObj_eta;   //!
   TBranch        *b_TrigObj_phi;   //!
   TBranch        *b_TrigObj_l1pt;   //!
   TBranch        *b_TrigObj_l1pt_2;   //!
   TBranch        *b_TrigObj_l2pt;   //!
   TBranch        *b_TrigObj_id;   //!
   TBranch        *b_TrigObj_l1iso;   //!
   TBranch        *b_TrigObj_l1charge;   //!
   TBranch        *b_TrigObj_filterBits;   //!
   TBranch        *b_nOtherPV;   //!
   TBranch        *b_OtherPV_z;   //!
   TBranch        *b_PV_ndof;   //!
   TBranch        *b_PV_x;   //!
   TBranch        *b_PV_y;   //!
   TBranch        *b_PV_z;   //!
   TBranch        *b_PV_chi2;   //!
   TBranch        *b_PV_score;   //!
   TBranch        *b_PV_npvs;   //!
   TBranch        *b_PV_npvsGood;   //!
   TBranch        *b_nSV;   //!
   TBranch        *b_SV_dlen;   //!
   TBranch        *b_SV_dlenSig;   //!
   TBranch        *b_SV_pAngle;   //!
   TBranch        *b_Electron_genPartIdx;   //!
   TBranch        *b_Electron_genPartFlav;   //!
   TBranch        *b_Muon_genPartIdx;   //!
   TBranch        *b_Muon_genPartFlav;   //!
   TBranch        *b_SV_chi2;   //!
   TBranch        *b_SV_eta;   //!
   TBranch        *b_SV_mass;   //!
   TBranch        *b_SV_ndof;   //!
   TBranch        *b_SV_phi;   //!
   TBranch        *b_SV_pt;   //!
   TBranch        *b_SV_x;   //!
   TBranch        *b_SV_y;   //!
   TBranch        *b_SV_z;   //!
   TBranch        *b_ProbeTracks_genPartIdx;   //!
   TBranch        *b_ProbeTracks_genPartFlav;   //!

   BParkBase(TTree *tree=0);
   virtual ~BParkBase();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef BParkBase_cxx
BParkBase::BParkBase(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/eos/cms/store/group/cmst3/group/bpark/BParkingNANO_2020Jan16/BuToKJpsi_Toee_Mufilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/crab_BuToKJpsi_Toee/200116_215618/0000/BParkNANO_mc_2020Jan16_1.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/eos/cms/store/group/cmst3/group/bpark/BParkingNANO_2020Jan16/BuToKJpsi_Toee_Mufilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/crab_BuToKJpsi_Toee/200116_215618/0000/BParkNANO_mc_2020Jan16_1.root");
      }
      f->GetObject("Events",tree);

   }
   Init(tree);
}

BParkBase::~BParkBase()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t BParkBase::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t BParkBase::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void BParkBase::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("luminosityBlock", &luminosityBlock, &b_luminosityBlock);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("nBToKEE", &nBToKEE, &b_nBToKEE);
   fChain->SetBranchAddress("BToKEE_b_iso03", BToKEE_b_iso03, &b_BToKEE_b_iso03);
   fChain->SetBranchAddress("BToKEE_b_iso04", BToKEE_b_iso04, &b_BToKEE_b_iso04);
   fChain->SetBranchAddress("BToKEE_cos2D", BToKEE_cos2D, &b_BToKEE_cos2D);
   fChain->SetBranchAddress("BToKEE_eta", BToKEE_eta, &b_BToKEE_eta);
   fChain->SetBranchAddress("BToKEE_fit_cos2D", BToKEE_fit_cos2D, &b_BToKEE_fit_cos2D);
   fChain->SetBranchAddress("BToKEE_fit_eta", BToKEE_fit_eta, &b_BToKEE_fit_eta);
   fChain->SetBranchAddress("BToKEE_fit_k_eta", BToKEE_fit_k_eta, &b_BToKEE_fit_k_eta);
   fChain->SetBranchAddress("BToKEE_fit_k_phi", BToKEE_fit_k_phi, &b_BToKEE_fit_k_phi);
   fChain->SetBranchAddress("BToKEE_fit_k_pt", BToKEE_fit_k_pt, &b_BToKEE_fit_k_pt);
   fChain->SetBranchAddress("BToKEE_fit_l1_eta", BToKEE_fit_l1_eta, &b_BToKEE_fit_l1_eta);
   fChain->SetBranchAddress("BToKEE_fit_l1_phi", BToKEE_fit_l1_phi, &b_BToKEE_fit_l1_phi);
   fChain->SetBranchAddress("BToKEE_fit_l1_pt", BToKEE_fit_l1_pt, &b_BToKEE_fit_l1_pt);
   fChain->SetBranchAddress("BToKEE_fit_l2_eta", BToKEE_fit_l2_eta, &b_BToKEE_fit_l2_eta);
   fChain->SetBranchAddress("BToKEE_fit_l2_phi", BToKEE_fit_l2_phi, &b_BToKEE_fit_l2_phi);
   fChain->SetBranchAddress("BToKEE_fit_l2_pt", BToKEE_fit_l2_pt, &b_BToKEE_fit_l2_pt);
   fChain->SetBranchAddress("BToKEE_fit_mass", BToKEE_fit_mass, &b_BToKEE_fit_mass);
   fChain->SetBranchAddress("BToKEE_fit_massErr", BToKEE_fit_massErr, &b_BToKEE_fit_massErr);
   fChain->SetBranchAddress("BToKEE_fit_phi", BToKEE_fit_phi, &b_BToKEE_fit_phi);
   fChain->SetBranchAddress("BToKEE_fit_pt", BToKEE_fit_pt, &b_BToKEE_fit_pt);
   fChain->SetBranchAddress("BToKEE_k_iso03", BToKEE_k_iso03, &b_BToKEE_k_iso03);
   fChain->SetBranchAddress("BToKEE_k_iso04", BToKEE_k_iso04, &b_BToKEE_k_iso04);
   fChain->SetBranchAddress("BToKEE_l1_iso03", BToKEE_l1_iso03, &b_BToKEE_l1_iso03);
   fChain->SetBranchAddress("BToKEE_l1_iso04", BToKEE_l1_iso04, &b_BToKEE_l1_iso04);
   fChain->SetBranchAddress("BToKEE_l2_iso03", BToKEE_l2_iso03, &b_BToKEE_l2_iso03);
   fChain->SetBranchAddress("BToKEE_l2_iso04", BToKEE_l2_iso04, &b_BToKEE_l2_iso04);
   fChain->SetBranchAddress("BToKEE_l_xy", BToKEE_l_xy, &b_BToKEE_l_xy);
   fChain->SetBranchAddress("BToKEE_l_xy_unc", BToKEE_l_xy_unc, &b_BToKEE_l_xy_unc);
   fChain->SetBranchAddress("BToKEE_mass", BToKEE_mass, &b_BToKEE_mass);
   fChain->SetBranchAddress("BToKEE_maxDR", BToKEE_maxDR, &b_BToKEE_maxDR);
   fChain->SetBranchAddress("BToKEE_minDR", BToKEE_minDR, &b_BToKEE_minDR);
   fChain->SetBranchAddress("BToKEE_mllErr_llfit", BToKEE_mllErr_llfit, &b_BToKEE_mllErr_llfit);
   fChain->SetBranchAddress("BToKEE_mll_fullfit", BToKEE_mll_fullfit, &b_BToKEE_mll_fullfit);
   fChain->SetBranchAddress("BToKEE_mll_llfit", BToKEE_mll_llfit, &b_BToKEE_mll_llfit);
   fChain->SetBranchAddress("BToKEE_mll_raw", BToKEE_mll_raw, &b_BToKEE_mll_raw);
   fChain->SetBranchAddress("BToKEE_phi", BToKEE_phi, &b_BToKEE_phi);
   fChain->SetBranchAddress("BToKEE_pt", BToKEE_pt, &b_BToKEE_pt);
   fChain->SetBranchAddress("BToKEE_svprob", BToKEE_svprob, &b_BToKEE_svprob);
   fChain->SetBranchAddress("BToKEE_vtx_ex", BToKEE_vtx_ex, &b_BToKEE_vtx_ex);
   fChain->SetBranchAddress("BToKEE_vtx_ey", BToKEE_vtx_ey, &b_BToKEE_vtx_ey);
   fChain->SetBranchAddress("BToKEE_vtx_ez", BToKEE_vtx_ez, &b_BToKEE_vtx_ez);
   fChain->SetBranchAddress("BToKEE_vtx_x", BToKEE_vtx_x, &b_BToKEE_vtx_x);
   fChain->SetBranchAddress("BToKEE_vtx_y", BToKEE_vtx_y, &b_BToKEE_vtx_y);
   fChain->SetBranchAddress("BToKEE_vtx_z", BToKEE_vtx_z, &b_BToKEE_vtx_z);
   fChain->SetBranchAddress("BToKEE_charge", BToKEE_charge, &b_BToKEE_charge);
   fChain->SetBranchAddress("BToKEE_kIdx", BToKEE_kIdx, &b_BToKEE_kIdx);
   fChain->SetBranchAddress("BToKEE_l1Idx", BToKEE_l1Idx, &b_BToKEE_l1Idx);
   fChain->SetBranchAddress("BToKEE_l2Idx", BToKEE_l2Idx, &b_BToKEE_l2Idx);
   fChain->SetBranchAddress("BToKEE_n_k_used", BToKEE_n_k_used, &b_BToKEE_n_k_used);
   fChain->SetBranchAddress("BToKEE_n_l1_used", BToKEE_n_l1_used, &b_BToKEE_n_l1_used);
   fChain->SetBranchAddress("BToKEE_n_l2_used", BToKEE_n_l2_used, &b_BToKEE_n_l2_used);
   fChain->SetBranchAddress("BToKEE_pdgId", BToKEE_pdgId, &b_BToKEE_pdgId);
   fChain->SetBranchAddress("nBToKMuMu", &nBToKMuMu, &b_nBToKMuMu);
   fChain->SetBranchAddress("BToKMuMu_b_iso03", BToKMuMu_b_iso03, &b_BToKMuMu_b_iso03);
   fChain->SetBranchAddress("BToKMuMu_b_iso04", BToKMuMu_b_iso04, &b_BToKMuMu_b_iso04);
   fChain->SetBranchAddress("BToKMuMu_cos2D", BToKMuMu_cos2D, &b_BToKMuMu_cos2D);
   fChain->SetBranchAddress("BToKMuMu_eta", BToKMuMu_eta, &b_BToKMuMu_eta);
   fChain->SetBranchAddress("BToKMuMu_fit_cos2D", BToKMuMu_fit_cos2D, &b_BToKMuMu_fit_cos2D);
   fChain->SetBranchAddress("BToKMuMu_fit_eta", BToKMuMu_fit_eta, &b_BToKMuMu_fit_eta);
   fChain->SetBranchAddress("BToKMuMu_fit_k_eta", BToKMuMu_fit_k_eta, &b_BToKMuMu_fit_k_eta);
   fChain->SetBranchAddress("BToKMuMu_fit_k_phi", BToKMuMu_fit_k_phi, &b_BToKMuMu_fit_k_phi);
   fChain->SetBranchAddress("BToKMuMu_fit_k_pt", BToKMuMu_fit_k_pt, &b_BToKMuMu_fit_k_pt);
   fChain->SetBranchAddress("BToKMuMu_fit_l1_eta", BToKMuMu_fit_l1_eta, &b_BToKMuMu_fit_l1_eta);
   fChain->SetBranchAddress("BToKMuMu_fit_l1_phi", BToKMuMu_fit_l1_phi, &b_BToKMuMu_fit_l1_phi);
   fChain->SetBranchAddress("BToKMuMu_fit_l1_pt", BToKMuMu_fit_l1_pt, &b_BToKMuMu_fit_l1_pt);
   fChain->SetBranchAddress("BToKMuMu_fit_l2_eta", BToKMuMu_fit_l2_eta, &b_BToKMuMu_fit_l2_eta);
   fChain->SetBranchAddress("BToKMuMu_fit_l2_phi", BToKMuMu_fit_l2_phi, &b_BToKMuMu_fit_l2_phi);
   fChain->SetBranchAddress("BToKMuMu_fit_l2_pt", BToKMuMu_fit_l2_pt, &b_BToKMuMu_fit_l2_pt);
   fChain->SetBranchAddress("BToKMuMu_fit_mass", BToKMuMu_fit_mass, &b_BToKMuMu_fit_mass);
   fChain->SetBranchAddress("BToKMuMu_fit_massErr", BToKMuMu_fit_massErr, &b_BToKMuMu_fit_massErr);
   fChain->SetBranchAddress("BToKMuMu_fit_phi", BToKMuMu_fit_phi, &b_BToKMuMu_fit_phi);
   fChain->SetBranchAddress("BToKMuMu_fit_pt", BToKMuMu_fit_pt, &b_BToKMuMu_fit_pt);
   fChain->SetBranchAddress("BToKMuMu_k_iso03", BToKMuMu_k_iso03, &b_BToKMuMu_k_iso03);
   fChain->SetBranchAddress("BToKMuMu_k_iso04", BToKMuMu_k_iso04, &b_BToKMuMu_k_iso04);
   fChain->SetBranchAddress("BToKMuMu_l1_iso03", BToKMuMu_l1_iso03, &b_BToKMuMu_l1_iso03);
   fChain->SetBranchAddress("BToKMuMu_l1_iso04", BToKMuMu_l1_iso04, &b_BToKMuMu_l1_iso04);
   fChain->SetBranchAddress("BToKMuMu_l2_iso03", BToKMuMu_l2_iso03, &b_BToKMuMu_l2_iso03);
   fChain->SetBranchAddress("BToKMuMu_l2_iso04", BToKMuMu_l2_iso04, &b_BToKMuMu_l2_iso04);
   fChain->SetBranchAddress("BToKMuMu_l_xy", BToKMuMu_l_xy, &b_BToKMuMu_l_xy);
   fChain->SetBranchAddress("BToKMuMu_l_xy_unc", BToKMuMu_l_xy_unc, &b_BToKMuMu_l_xy_unc);
   fChain->SetBranchAddress("BToKMuMu_mass", BToKMuMu_mass, &b_BToKMuMu_mass);
   fChain->SetBranchAddress("BToKMuMu_maxDR", BToKMuMu_maxDR, &b_BToKMuMu_maxDR);
   fChain->SetBranchAddress("BToKMuMu_minDR", BToKMuMu_minDR, &b_BToKMuMu_minDR);
   fChain->SetBranchAddress("BToKMuMu_mllErr_llfit", BToKMuMu_mllErr_llfit, &b_BToKMuMu_mllErr_llfit);
   fChain->SetBranchAddress("BToKMuMu_mll_fullfit", BToKMuMu_mll_fullfit, &b_BToKMuMu_mll_fullfit);
   fChain->SetBranchAddress("BToKMuMu_mll_llfit", BToKMuMu_mll_llfit, &b_BToKMuMu_mll_llfit);
   fChain->SetBranchAddress("BToKMuMu_mll_raw", BToKMuMu_mll_raw, &b_BToKMuMu_mll_raw);
   fChain->SetBranchAddress("BToKMuMu_phi", BToKMuMu_phi, &b_BToKMuMu_phi);
   fChain->SetBranchAddress("BToKMuMu_pt", BToKMuMu_pt, &b_BToKMuMu_pt);
   fChain->SetBranchAddress("BToKMuMu_svprob", BToKMuMu_svprob, &b_BToKMuMu_svprob);
   fChain->SetBranchAddress("BToKMuMu_vtx_ex", BToKMuMu_vtx_ex, &b_BToKMuMu_vtx_ex);
   fChain->SetBranchAddress("BToKMuMu_vtx_ey", BToKMuMu_vtx_ey, &b_BToKMuMu_vtx_ey);
   fChain->SetBranchAddress("BToKMuMu_vtx_ez", BToKMuMu_vtx_ez, &b_BToKMuMu_vtx_ez);
   fChain->SetBranchAddress("BToKMuMu_vtx_x", BToKMuMu_vtx_x, &b_BToKMuMu_vtx_x);
   fChain->SetBranchAddress("BToKMuMu_vtx_y", BToKMuMu_vtx_y, &b_BToKMuMu_vtx_y);
   fChain->SetBranchAddress("BToKMuMu_vtx_z", BToKMuMu_vtx_z, &b_BToKMuMu_vtx_z);
   fChain->SetBranchAddress("BToKMuMu_charge", BToKMuMu_charge, &b_BToKMuMu_charge);
   fChain->SetBranchAddress("BToKMuMu_kIdx", BToKMuMu_kIdx, &b_BToKMuMu_kIdx);
   fChain->SetBranchAddress("BToKMuMu_l1Idx", BToKMuMu_l1Idx, &b_BToKMuMu_l1Idx);
   fChain->SetBranchAddress("BToKMuMu_l2Idx", BToKMuMu_l2Idx, &b_BToKMuMu_l2Idx);
   fChain->SetBranchAddress("BToKMuMu_n_k_used", BToKMuMu_n_k_used, &b_BToKMuMu_n_k_used);
   fChain->SetBranchAddress("BToKMuMu_n_l1_used", BToKMuMu_n_l1_used, &b_BToKMuMu_n_l1_used);
   fChain->SetBranchAddress("BToKMuMu_n_l2_used", BToKMuMu_n_l2_used, &b_BToKMuMu_n_l2_used);
   fChain->SetBranchAddress("BToKMuMu_pdgId", BToKMuMu_pdgId, &b_BToKMuMu_pdgId);
   fChain->SetBranchAddress("nElectron", &nElectron, &b_nElectron);
   fChain->SetBranchAddress("Electron_dxy", Electron_dxy, &b_Electron_dxy);
   fChain->SetBranchAddress("Electron_dxyErr", Electron_dxyErr, &b_Electron_dxyErr);
   fChain->SetBranchAddress("Electron_dz", Electron_dz, &b_Electron_dz);
   fChain->SetBranchAddress("Electron_dzErr", Electron_dzErr, &b_Electron_dzErr);
   fChain->SetBranchAddress("Electron_eta", Electron_eta, &b_Electron_eta);
   fChain->SetBranchAddress("Electron_fBrem", Electron_fBrem, &b_Electron_fBrem);
   fChain->SetBranchAddress("Electron_ip3d", Electron_ip3d, &b_Electron_ip3d);
   fChain->SetBranchAddress("Electron_mass", Electron_mass, &b_Electron_mass);
   fChain->SetBranchAddress("Electron_mvaId", Electron_mvaId, &b_Electron_mvaId);
   fChain->SetBranchAddress("Electron_pfRelIso", Electron_pfRelIso, &b_Electron_pfRelIso);
   fChain->SetBranchAddress("Electron_pfmvaId", Electron_pfmvaId, &b_Electron_pfmvaId);
   fChain->SetBranchAddress("Electron_phi", Electron_phi, &b_Electron_phi);
   fChain->SetBranchAddress("Electron_pt", Electron_pt, &b_Electron_pt);
   fChain->SetBranchAddress("Electron_ptBiased", Electron_ptBiased, &b_Electron_ptBiased);
   fChain->SetBranchAddress("Electron_sip3d", Electron_sip3d, &b_Electron_sip3d);
   fChain->SetBranchAddress("Electron_trkRelIso", Electron_trkRelIso, &b_Electron_trkRelIso);
   fChain->SetBranchAddress("Electron_unBiased", Electron_unBiased, &b_Electron_unBiased);
   fChain->SetBranchAddress("Electron_vx", Electron_vx, &b_Electron_vx);
   fChain->SetBranchAddress("Electron_vy", Electron_vy, &b_Electron_vy);
   fChain->SetBranchAddress("Electron_vz", Electron_vz, &b_Electron_vz);
   fChain->SetBranchAddress("Electron_charge", Electron_charge, &b_Electron_charge);
   fChain->SetBranchAddress("Electron_pdgId", Electron_pdgId, &b_Electron_pdgId);
   fChain->SetBranchAddress("Electron_convVeto", Electron_convVeto, &b_Electron_convVeto);
   fChain->SetBranchAddress("Electron_isLowPt", Electron_isLowPt, &b_Electron_isLowPt);
   fChain->SetBranchAddress("Electron_isPF", Electron_isPF, &b_Electron_isPF);
   fChain->SetBranchAddress("Electron_isPFoverlap", Electron_isPFoverlap, &b_Electron_isPFoverlap);
   fChain->SetBranchAddress("nGenPart", &nGenPart, &b_nGenPart);
   fChain->SetBranchAddress("GenPart_eta", GenPart_eta, &b_GenPart_eta);
   fChain->SetBranchAddress("GenPart_mass", GenPart_mass, &b_GenPart_mass);
   fChain->SetBranchAddress("GenPart_phi", GenPart_phi, &b_GenPart_phi);
   fChain->SetBranchAddress("GenPart_pt", GenPart_pt, &b_GenPart_pt);
   fChain->SetBranchAddress("GenPart_vx", GenPart_vx, &b_GenPart_vx);
   fChain->SetBranchAddress("GenPart_vy", GenPart_vy, &b_GenPart_vy);
   fChain->SetBranchAddress("GenPart_vz", GenPart_vz, &b_GenPart_vz);
   fChain->SetBranchAddress("GenPart_genPartIdxMother", GenPart_genPartIdxMother, &b_GenPart_genPartIdxMother);
   fChain->SetBranchAddress("GenPart_pdgId", GenPart_pdgId, &b_GenPart_pdgId);
   fChain->SetBranchAddress("GenPart_status", GenPart_status, &b_GenPart_status);
   fChain->SetBranchAddress("GenPart_statusFlags", GenPart_statusFlags, &b_GenPart_statusFlags);
   fChain->SetBranchAddress("Generator_binvar", &Generator_binvar, &b_Generator_binvar);
   fChain->SetBranchAddress("Generator_scalePDF", &Generator_scalePDF, &b_Generator_scalePDF);
   fChain->SetBranchAddress("Generator_weight", &Generator_weight, &b_Generator_weight);
   fChain->SetBranchAddress("Generator_x1", &Generator_x1, &b_Generator_x1);
   fChain->SetBranchAddress("Generator_x2", &Generator_x2, &b_Generator_x2);
   fChain->SetBranchAddress("Generator_xpdf1", &Generator_xpdf1, &b_Generator_xpdf1);
   fChain->SetBranchAddress("Generator_xpdf2", &Generator_xpdf2, &b_Generator_xpdf2);
   fChain->SetBranchAddress("Generator_id1", &Generator_id1, &b_Generator_id1);
   fChain->SetBranchAddress("Generator_id2", &Generator_id2, &b_Generator_id2);
   fChain->SetBranchAddress("genWeight", &genWeight, &b_genWeight);
   fChain->SetBranchAddress("nPSWeight", &nPSWeight, &b_nPSWeight);
   fChain->SetBranchAddress("PSWeight", PSWeight, &b_PSWeight);
   fChain->SetBranchAddress("nMuon", &nMuon, &b_nMuon);
   fChain->SetBranchAddress("Muon_dxy", Muon_dxy, &b_Muon_dxy);
   fChain->SetBranchAddress("Muon_dxyErr", Muon_dxyErr, &b_Muon_dxyErr);
   fChain->SetBranchAddress("Muon_dz", Muon_dz, &b_Muon_dz);
   fChain->SetBranchAddress("Muon_dzErr", Muon_dzErr, &b_Muon_dzErr);
   fChain->SetBranchAddress("Muon_eta", Muon_eta, &b_Muon_eta);
   fChain->SetBranchAddress("Muon_ip3d", Muon_ip3d, &b_Muon_ip3d);
   fChain->SetBranchAddress("Muon_mass", Muon_mass, &b_Muon_mass);
   fChain->SetBranchAddress("Muon_pfRelIso03_all", Muon_pfRelIso03_all, &b_Muon_pfRelIso03_all);
   fChain->SetBranchAddress("Muon_pfRelIso04_all", Muon_pfRelIso04_all, &b_Muon_pfRelIso04_all);
   fChain->SetBranchAddress("Muon_phi", Muon_phi, &b_Muon_phi);
   fChain->SetBranchAddress("Muon_pt", Muon_pt, &b_Muon_pt);
   fChain->SetBranchAddress("Muon_ptErr", Muon_ptErr, &b_Muon_ptErr);
   fChain->SetBranchAddress("Muon_sip3d", Muon_sip3d, &b_Muon_sip3d);
   fChain->SetBranchAddress("Muon_vx", Muon_vx, &b_Muon_vx);
   fChain->SetBranchAddress("Muon_vy", Muon_vy, &b_Muon_vy);
   fChain->SetBranchAddress("Muon_vz", Muon_vz, &b_Muon_vz);
   fChain->SetBranchAddress("Muon_charge", Muon_charge, &b_Muon_charge);
   fChain->SetBranchAddress("Muon_isTriggering", Muon_isTriggering, &b_Muon_isTriggering);
   fChain->SetBranchAddress("Muon_pdgId", Muon_pdgId, &b_Muon_pdgId);
   fChain->SetBranchAddress("Muon_isGlobal", Muon_isGlobal, &b_Muon_isGlobal);
   fChain->SetBranchAddress("Muon_isPFcand", Muon_isPFcand, &b_Muon_isPFcand);
   fChain->SetBranchAddress("Muon_isTracker", Muon_isTracker, &b_Muon_isTracker);
   fChain->SetBranchAddress("Muon_mediumId", Muon_mediumId, &b_Muon_mediumId);
   fChain->SetBranchAddress("Muon_pfIsoId", Muon_pfIsoId, &b_Muon_pfIsoId);
   fChain->SetBranchAddress("Muon_softId", Muon_softId, &b_Muon_softId);
   fChain->SetBranchAddress("Muon_tightId", Muon_tightId, &b_Muon_tightId);
   fChain->SetBranchAddress("Muon_tkIsoId", Muon_tkIsoId, &b_Muon_tkIsoId);
   fChain->SetBranchAddress("Muon_triggerIdLoose", Muon_triggerIdLoose, &b_Muon_triggerIdLoose);
   fChain->SetBranchAddress("nTriggerMuon", &nTriggerMuon, &b_nTriggerMuon);
   fChain->SetBranchAddress("TriggerMuon_eta", TriggerMuon_eta, &b_TriggerMuon_eta);
   fChain->SetBranchAddress("TriggerMuon_mass", TriggerMuon_mass, &b_TriggerMuon_mass);
   fChain->SetBranchAddress("TriggerMuon_phi", TriggerMuon_phi, &b_TriggerMuon_phi);
   fChain->SetBranchAddress("TriggerMuon_pt", TriggerMuon_pt, &b_TriggerMuon_pt);
   fChain->SetBranchAddress("TriggerMuon_vx", TriggerMuon_vx, &b_TriggerMuon_vx);
   fChain->SetBranchAddress("TriggerMuon_vy", TriggerMuon_vy, &b_TriggerMuon_vy);
   fChain->SetBranchAddress("TriggerMuon_vz", TriggerMuon_vz, &b_TriggerMuon_vz);
   fChain->SetBranchAddress("TriggerMuon_charge", TriggerMuon_charge, &b_TriggerMuon_charge);
   fChain->SetBranchAddress("TriggerMuon_pdgId", TriggerMuon_pdgId, &b_TriggerMuon_pdgId);
   fChain->SetBranchAddress("TriggerMuon_trgMuonIndex", TriggerMuon_trgMuonIndex, &b_TriggerMuon_trgMuonIndex);
   fChain->SetBranchAddress("Pileup_nTrueInt", &Pileup_nTrueInt, &b_Pileup_nTrueInt);
   fChain->SetBranchAddress("Pileup_pudensity", &Pileup_pudensity, &b_Pileup_pudensity);
   fChain->SetBranchAddress("Pileup_gpudensity", &Pileup_gpudensity, &b_Pileup_gpudensity);
   fChain->SetBranchAddress("Pileup_nPU", &Pileup_nPU, &b_Pileup_nPU);
   fChain->SetBranchAddress("Pileup_sumEOOT", &Pileup_sumEOOT, &b_Pileup_sumEOOT);
   fChain->SetBranchAddress("Pileup_sumLOOT", &Pileup_sumLOOT, &b_Pileup_sumLOOT);
   fChain->SetBranchAddress("fixedGridRhoFastjetAll", &fixedGridRhoFastjetAll, &b_fixedGridRhoFastjetAll);
   fChain->SetBranchAddress("fixedGridRhoFastjetCentral", &fixedGridRhoFastjetCentral, &b_fixedGridRhoFastjetCentral);
   fChain->SetBranchAddress("fixedGridRhoFastjetCentralCalo", &fixedGridRhoFastjetCentralCalo, &b_fixedGridRhoFastjetCentralCalo);
   fChain->SetBranchAddress("fixedGridRhoFastjetCentralChargedPileUp", &fixedGridRhoFastjetCentralChargedPileUp, &b_fixedGridRhoFastjetCentralChargedPileUp);
   fChain->SetBranchAddress("fixedGridRhoFastjetCentralNeutral", &fixedGridRhoFastjetCentralNeutral, &b_fixedGridRhoFastjetCentralNeutral);
   fChain->SetBranchAddress("nProbeTracks", &nProbeTracks, &b_nProbeTracks);
   fChain->SetBranchAddress("ProbeTracks_DCASig", ProbeTracks_DCASig, &b_ProbeTracks_DCASig);
   fChain->SetBranchAddress("ProbeTracks_dxy", ProbeTracks_dxy, &b_ProbeTracks_dxy);
   fChain->SetBranchAddress("ProbeTracks_dxyS", ProbeTracks_dxyS, &b_ProbeTracks_dxyS);
   fChain->SetBranchAddress("ProbeTracks_dz", ProbeTracks_dz, &b_ProbeTracks_dz);
   fChain->SetBranchAddress("ProbeTracks_dzS", ProbeTracks_dzS, &b_ProbeTracks_dzS);
   fChain->SetBranchAddress("ProbeTracks_eta", ProbeTracks_eta, &b_ProbeTracks_eta);
   fChain->SetBranchAddress("ProbeTracks_mass", ProbeTracks_mass, &b_ProbeTracks_mass);
   fChain->SetBranchAddress("ProbeTracks_phi", ProbeTracks_phi, &b_ProbeTracks_phi);
   fChain->SetBranchAddress("ProbeTracks_pt", ProbeTracks_pt, &b_ProbeTracks_pt);
   fChain->SetBranchAddress("ProbeTracks_vx", ProbeTracks_vx, &b_ProbeTracks_vx);
   fChain->SetBranchAddress("ProbeTracks_vy", ProbeTracks_vy, &b_ProbeTracks_vy);
   fChain->SetBranchAddress("ProbeTracks_vz", ProbeTracks_vz, &b_ProbeTracks_vz);
   fChain->SetBranchAddress("ProbeTracks_charge", ProbeTracks_charge, &b_ProbeTracks_charge);
   fChain->SetBranchAddress("ProbeTracks_isLostTrk", ProbeTracks_isLostTrk, &b_ProbeTracks_isLostTrk);
   fChain->SetBranchAddress("ProbeTracks_isPacked", ProbeTracks_isPacked, &b_ProbeTracks_isPacked);
   fChain->SetBranchAddress("ProbeTracks_nValidHits", ProbeTracks_nValidHits, &b_ProbeTracks_nValidHits);
   fChain->SetBranchAddress("ProbeTracks_pdgId", ProbeTracks_pdgId, &b_ProbeTracks_pdgId);
   fChain->SetBranchAddress("ProbeTracks_isMatchedToEle", ProbeTracks_isMatchedToEle, &b_ProbeTracks_isMatchedToEle);
   fChain->SetBranchAddress("ProbeTracks_isMatchedToLooseMuon", ProbeTracks_isMatchedToLooseMuon, &b_ProbeTracks_isMatchedToLooseMuon);
   fChain->SetBranchAddress("ProbeTracks_isMatchedToMediumMuon", ProbeTracks_isMatchedToMediumMuon, &b_ProbeTracks_isMatchedToMediumMuon);
   fChain->SetBranchAddress("ProbeTracks_isMatchedToMuon", ProbeTracks_isMatchedToMuon, &b_ProbeTracks_isMatchedToMuon);
   fChain->SetBranchAddress("ProbeTracks_isMatchedToSoftMuon", ProbeTracks_isMatchedToSoftMuon, &b_ProbeTracks_isMatchedToSoftMuon);
   fChain->SetBranchAddress("HLT_Mu7_IP4", &HLT_Mu7_IP4, &b_HLT_Mu7_IP4);
   fChain->SetBranchAddress("HLT_Mu8_IP6", &HLT_Mu8_IP6, &b_HLT_Mu8_IP6);
   fChain->SetBranchAddress("HLT_Mu8_IP5", &HLT_Mu8_IP5, &b_HLT_Mu8_IP5);
   fChain->SetBranchAddress("HLT_Mu8_IP3", &HLT_Mu8_IP3, &b_HLT_Mu8_IP3);
   fChain->SetBranchAddress("HLT_Mu8p5_IP3p5", &HLT_Mu8p5_IP3p5, &b_HLT_Mu8p5_IP3p5);
   fChain->SetBranchAddress("HLT_Mu9_IP6", &HLT_Mu9_IP6, &b_HLT_Mu9_IP6);
   fChain->SetBranchAddress("HLT_Mu9_IP5", &HLT_Mu9_IP5, &b_HLT_Mu9_IP5);
   fChain->SetBranchAddress("HLT_Mu9_IP4", &HLT_Mu9_IP4, &b_HLT_Mu9_IP4);
   fChain->SetBranchAddress("HLT_Mu10p5_IP3p5", &HLT_Mu10p5_IP3p5, &b_HLT_Mu10p5_IP3p5);
   fChain->SetBranchAddress("HLT_Mu12_IP6", &HLT_Mu12_IP6, &b_HLT_Mu12_IP6);
   fChain->SetBranchAddress("L1_SingleMu7er1p5", &L1_SingleMu7er1p5, &b_L1_SingleMu7er1p5);
   fChain->SetBranchAddress("L1_SingleMu8er1p5", &L1_SingleMu8er1p5, &b_L1_SingleMu8er1p5);
   fChain->SetBranchAddress("L1_SingleMu9er1p5", &L1_SingleMu9er1p5, &b_L1_SingleMu9er1p5);
   fChain->SetBranchAddress("L1_SingleMu10er1p5", &L1_SingleMu10er1p5, &b_L1_SingleMu10er1p5);
   fChain->SetBranchAddress("L1_SingleMu12er1p5", &L1_SingleMu12er1p5, &b_L1_SingleMu12er1p5);
   fChain->SetBranchAddress("L1_SingleMu22", &L1_SingleMu22, &b_L1_SingleMu22);
   fChain->SetBranchAddress("nTrigObj", &nTrigObj, &b_nTrigObj);
   fChain->SetBranchAddress("TrigObj_pt", TrigObj_pt, &b_TrigObj_pt);
   fChain->SetBranchAddress("TrigObj_eta", TrigObj_eta, &b_TrigObj_eta);
   fChain->SetBranchAddress("TrigObj_phi", TrigObj_phi, &b_TrigObj_phi);
   fChain->SetBranchAddress("TrigObj_l1pt", TrigObj_l1pt, &b_TrigObj_l1pt);
   fChain->SetBranchAddress("TrigObj_l1pt_2", TrigObj_l1pt_2, &b_TrigObj_l1pt_2);
   fChain->SetBranchAddress("TrigObj_l2pt", TrigObj_l2pt, &b_TrigObj_l2pt);
   fChain->SetBranchAddress("TrigObj_id", TrigObj_id, &b_TrigObj_id);
   fChain->SetBranchAddress("TrigObj_l1iso", TrigObj_l1iso, &b_TrigObj_l1iso);
   fChain->SetBranchAddress("TrigObj_l1charge", TrigObj_l1charge, &b_TrigObj_l1charge);
   fChain->SetBranchAddress("TrigObj_filterBits", TrigObj_filterBits, &b_TrigObj_filterBits);
   fChain->SetBranchAddress("nOtherPV", &nOtherPV, &b_nOtherPV);
   fChain->SetBranchAddress("OtherPV_z", OtherPV_z, &b_OtherPV_z);
   fChain->SetBranchAddress("PV_ndof", &PV_ndof, &b_PV_ndof);
   fChain->SetBranchAddress("PV_x", &PV_x, &b_PV_x);
   fChain->SetBranchAddress("PV_y", &PV_y, &b_PV_y);
   fChain->SetBranchAddress("PV_z", &PV_z, &b_PV_z);
   fChain->SetBranchAddress("PV_chi2", &PV_chi2, &b_PV_chi2);
   fChain->SetBranchAddress("PV_score", &PV_score, &b_PV_score);
   fChain->SetBranchAddress("PV_npvs", &PV_npvs, &b_PV_npvs);
   fChain->SetBranchAddress("PV_npvsGood", &PV_npvsGood, &b_PV_npvsGood);
   fChain->SetBranchAddress("nSV", &nSV, &b_nSV);
   fChain->SetBranchAddress("SV_dlen", SV_dlen, &b_SV_dlen);
   fChain->SetBranchAddress("SV_dlenSig", SV_dlenSig, &b_SV_dlenSig);
   fChain->SetBranchAddress("SV_pAngle", SV_pAngle, &b_SV_pAngle);
   fChain->SetBranchAddress("Electron_genPartIdx", Electron_genPartIdx, &b_Electron_genPartIdx);
   fChain->SetBranchAddress("Electron_genPartFlav", Electron_genPartFlav, &b_Electron_genPartFlav);
   fChain->SetBranchAddress("Muon_genPartIdx", Muon_genPartIdx, &b_Muon_genPartIdx);
   fChain->SetBranchAddress("Muon_genPartFlav", Muon_genPartFlav, &b_Muon_genPartFlav);
   fChain->SetBranchAddress("SV_chi2", SV_chi2, &b_SV_chi2);
   fChain->SetBranchAddress("SV_eta", SV_eta, &b_SV_eta);
   fChain->SetBranchAddress("SV_mass", SV_mass, &b_SV_mass);
   fChain->SetBranchAddress("SV_ndof", SV_ndof, &b_SV_ndof);
   fChain->SetBranchAddress("SV_phi", SV_phi, &b_SV_phi);
   fChain->SetBranchAddress("SV_pt", SV_pt, &b_SV_pt);
   fChain->SetBranchAddress("SV_x", SV_x, &b_SV_x);
   fChain->SetBranchAddress("SV_y", SV_y, &b_SV_y);
   fChain->SetBranchAddress("SV_z", SV_z, &b_SV_z);
   fChain->SetBranchAddress("ProbeTracks_genPartIdx", ProbeTracks_genPartIdx, &b_ProbeTracks_genPartIdx);
   fChain->SetBranchAddress("ProbeTracks_genPartFlav", ProbeTracks_genPartFlav, &b_ProbeTracks_genPartFlav);
   Notify();
}

Bool_t BParkBase::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void BParkBase::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t BParkBase::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef BParkBase_cxx
