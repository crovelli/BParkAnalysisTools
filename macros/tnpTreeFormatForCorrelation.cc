#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TH1.h>
#include <iostream>
#include <vector>
#include <TRandom.h>
#include "TLorentzVector.h"

// to produce formatted TnP ntuples, to check the correlation 
// between variables, separately for signal and background

using namespace std;

void tnpTreeFormatForCorrelation(const char* filename, int isProbeLpt, int isSignal) {

  cout << "Formatting " << filename << endl;  

  // Original ntuple   
  TFile *fileOrig = 0;
  TTree *treeOrig = 0;

  fileOrig = TFile::Open(TString("/tmp/crovelli/")+TString(filename));
  if( fileOrig ) {
    fileOrig->cd();
    treeOrig = (TTree*)fileOrig->Get("TaPtree");
  } else {
    cout << "File " << filename << " not existing !" << endl;
    return;
  }
  
  fileOrig->cd();
  if (!treeOrig) {
    cout << "TaPTree not existing !" << endl; 
    return;    
  }

  treeOrig->SetMakeClass(0);
  cout << "TreeOrig->Size = "<< treeOrig->GetEntries() << endl;

  // number of entries saved in the first tree
  int nentriesOrig = treeOrig->GetEntries();   

  // Tree for the final format
  TFile *fileNew = TFile::Open(TString("/tmp/crovelli/Formatted_")+TString(filename),"recreate");
  fileNew->ls();
  fileNew->cd();
  TDirectory *myDir = (TDirectory*)fileNew->mkdir("tnpAna");
  myDir->cd();
  TTree *treeNew = new TTree("fitter_tree","reduced tree for T&P");
  treeNew->SetAutoSave(-99999999999);
  treeNew->SetAutoFlush(-99999999999);

  std::vector<TTree*> trees; 
  trees.push_back(treeNew);

  // original tree leaves
  Int_t           theRun = 0;
  Int_t           theEvent = 0;
  Int_t           nvtx = 0;
  Int_t           sampleID = 0;
  Float_t         rho = 0;
  Float_t         pu_weight = 0;
  Float_t         pu_n = 0;
  Float_t         perEveW = 0;
  Int_t           hlt12ip6 = 0;
  Int_t           hlt9ip6 = 0;
  Int_t           hlt9ip5 = 0;
  Int_t           hlt9ip4 = 0;
  Int_t           hlt7ip4 = 0;
  Int_t           hlt8ip6 = 0;
  Int_t           hlt8ip5 = 0;
  Int_t           hlt8ip3 = 0;
  Int_t           hlt8d5ip3d5 = 0;
  Int_t           hlt10d5ip3d5 = 0;
  Int_t           selectedBSize = 0;
  vector<float>   *tag_pt = nullptr;
  vector<float>   *tag_eta = nullptr;
  vector<float>   *tag_phi = nullptr;
  vector<bool>    *tag_isPF = nullptr;
  vector<bool>    *tag_isPFOverlap = nullptr;
  vector<bool>    *tag_isLowPt = nullptr;
  vector<float>   *tag_mvaId = nullptr;
  vector<float>   *tag_pfmvaId = nullptr;
  vector<int>     *tag_convveto = nullptr;
  vector<bool>    *tag_matchMc = nullptr;
  vector<bool>    *tag_matchMcFromJPsi = nullptr;
  vector<float>   *probe_Bmass = nullptr;
  vector<float>   *probe_Bpt = nullptr;
  vector<float>   *probe_Bcos2D = nullptr;
  vector<float>   *probe_Bsvprob = nullptr;
  vector<float>   *probe_Bxysig = nullptr;
  vector<bool>    *probe_BmatchMC = nullptr;
  vector<float>   *probe_Kpt = nullptr;
  vector<float>   *probe_pt = nullptr;
  vector<float>   *probe_eta = nullptr;
  vector<float>   *probe_phi = nullptr;
  vector<bool>    *probe_isPF = nullptr;
  vector<bool>    *probe_isPFOverlap = nullptr;
  vector<bool>    *probe_isLowPt = nullptr;
  vector<float>   *probe_mvaId = nullptr;
  vector<float>   *probe_pfmvaId = nullptr;
  vector<float>   *probe_dxySig = nullptr;
  vector<float>   *probe_dzSig = nullptr;
  vector<float>   *probe_fBrem = nullptr;
  vector<float>   *probe_unBiased = nullptr;
  vector<float>   *probe_ptBiased = nullptr;
  vector<float>   *probe_normpt = nullptr;  
  vector<float>   *probe_dzTrg = nullptr;  
  vector<float>   *probe_iso04_rel = nullptr;
  vector<float>   *probe_invMass = nullptr;
  vector<int>     *probe_convveto = nullptr;
  vector<bool>    *probe_matchMc = nullptr;
  vector<bool>    *probe_matchMcFromJPsi = nullptr;
  Int_t           selectedPairsSize = 0;
  
  // List of branches - original tree
  TBranch        *b_theRun;   //!
  TBranch        *b_theEvent;   //!
  TBranch        *b_nvtx;   //!
  TBranch        *b_sampleID;   //!
  TBranch        *b_rho;   //!
  TBranch        *b_pu_weight;   //!
  TBranch        *b_pu_n;   //!
  TBranch        *b_perEveW;   //!
  TBranch        *b_hlt12ip6;
  TBranch        *b_hlt9ip6;
  TBranch        *b_hlt9ip5;
  TBranch        *b_hlt9ip4;
  TBranch        *b_hlt7ip4;
  TBranch        *b_hlt8ip6;
  TBranch        *b_hlt8ip5;
  TBranch        *b_hlt8ip3;
  TBranch        *b_hlt8d5ip3d5;
  TBranch        *b_hlt10d5ip3d5;
  TBranch        *b_selectedBSize;   //!
  TBranch        *b_tag_pt;   //!
  TBranch        *b_tag_eta;   //!
  TBranch        *b_tag_phi;   //!
  TBranch        *b_tag_isPF;   //!
  TBranch        *b_tag_isPFOverlap;   //!
  TBranch        *b_tag_isLowPt;   //!
  TBranch        *b_tag_mvaId;   //!
  TBranch        *b_tag_pfmvaId;   //!
  TBranch        *b_tag_convveto;   //!
  TBranch        *b_tag_matchMc;   //!
  TBranch        *b_tag_matchMcFromJPsi;   //!
  TBranch        *b_probe_Bmass;   //!
  TBranch        *b_probe_Bpt;   //!
  TBranch        *b_probe_Bcos2D;   //!
  TBranch        *b_probe_Bsvprob;   //!
  TBranch        *b_probe_Bxysig;   //!
  TBranch        *b_probe_BmatchMC;   //!
  TBranch        *b_probe_Kpt;   //!
  TBranch        *b_probe_pt;   //!
  TBranch        *b_probe_eta;   //!
  TBranch        *b_probe_phi;   //!
  TBranch        *b_probe_isPF;   //!
  TBranch        *b_probe_isPFOverlap;   //!
  TBranch        *b_probe_isLowPt;   //!
  TBranch        *b_probe_mvaId;   //!
  TBranch        *b_probe_pfmvaId;   //!
  TBranch        *b_probe_dxySig;   //!
  TBranch        *b_probe_dzSig;   //!
  TBranch        *b_probe_fBrem;   //!
  TBranch        *b_probe_unBiased;   //!
  TBranch        *b_probe_ptBiased;   //!
  TBranch        *b_probe_normpt; //!
  TBranch        *b_probe_dzTrg;  //!
  TBranch        *b_probe_iso04_rel; //!
  TBranch        *b_probe_invMass;   //!
  TBranch        *b_probe_convveto;   //!
  TBranch        *b_probe_matchMc;   //!
  TBranch        *b_probe_matchMcFromJPsi;   //!
  TBranch        *b_selectedPairsSize;   //!
  
  // Set branch addresses and branch pointers 
  treeOrig->SetBranchAddress("theRun", &theRun, &b_theRun);
  treeOrig->SetBranchAddress("theEvent", &theEvent, &b_theEvent);
  treeOrig->SetBranchAddress("nvtx", &nvtx, &b_nvtx);
  treeOrig->SetBranchAddress("sampleID", &sampleID, &b_sampleID);
  treeOrig->SetBranchAddress("rho", &rho, &b_rho);
  treeOrig->SetBranchAddress("pu_weight", &pu_weight, &b_pu_weight);
  treeOrig->SetBranchAddress("pu_n", &pu_n, &b_pu_n);
  treeOrig->SetBranchAddress("perEveW", &perEveW, &b_perEveW);
  treeOrig->SetBranchAddress("hlt12ip6", &hlt12ip6, &b_hlt12ip6);
  treeOrig->SetBranchAddress("hlt9ip6", &hlt9ip6, &b_hlt9ip6);
  treeOrig->SetBranchAddress("hlt9ip5", &hlt9ip5, &b_hlt9ip5);
  treeOrig->SetBranchAddress("hlt9ip4", &hlt9ip4, &b_hlt9ip4);
  treeOrig->SetBranchAddress("hlt7ip4", &hlt7ip4, &b_hlt7ip4);
  treeOrig->SetBranchAddress("hlt8ip6", &hlt8ip6, &b_hlt8ip6);
  treeOrig->SetBranchAddress("hlt8ip5", &hlt8ip5, &b_hlt8ip5);
  treeOrig->SetBranchAddress("hlt8ip3", &hlt8ip3, &b_hlt8ip3);
  treeOrig->SetBranchAddress("hlt8d5ip3d5", &hlt8d5ip3d5, &b_hlt8d5ip3d5);
  treeOrig->SetBranchAddress("hlt10d5ip3d5", &hlt10d5ip3d5, &b_hlt10d5ip3d5);
  treeOrig->SetBranchAddress("selectedBSize", &selectedBSize, &b_selectedBSize);
  treeOrig->SetBranchAddress("tag_pt", &tag_pt, &b_tag_pt);
  treeOrig->SetBranchAddress("tag_eta", &tag_eta, &b_tag_eta);
  treeOrig->SetBranchAddress("tag_phi", &tag_phi, &b_tag_phi);
  treeOrig->SetBranchAddress("tag_isPF", &tag_isPF, &b_tag_isPF);
  treeOrig->SetBranchAddress("tag_isPFOverlap", &tag_isPFOverlap, &b_tag_isPFOverlap);
  treeOrig->SetBranchAddress("tag_isLowPt", &tag_isLowPt, &b_tag_isLowPt);
  treeOrig->SetBranchAddress("tag_mvaId", &tag_mvaId, &b_tag_mvaId);
  treeOrig->SetBranchAddress("tag_pfmvaId", &tag_pfmvaId, &b_tag_pfmvaId);
  treeOrig->SetBranchAddress("tag_convveto", &tag_convveto, &b_tag_convveto);
  treeOrig->SetBranchAddress("tag_matchMc", &tag_matchMc, &b_tag_matchMc);
  treeOrig->SetBranchAddress("tag_matchMcFromJPsi", &tag_matchMcFromJPsi, &b_tag_matchMcFromJPsi);
  treeOrig->SetBranchAddress("probe_Bmass", &probe_Bmass, &b_probe_Bmass);
  treeOrig->SetBranchAddress("probe_Bpt", &probe_Bpt, &b_probe_Bpt);
  treeOrig->SetBranchAddress("probe_Bcos2D", &probe_Bcos2D, &b_probe_Bcos2D);
  treeOrig->SetBranchAddress("probe_Bsvprob", &probe_Bsvprob, &b_probe_Bsvprob);
  treeOrig->SetBranchAddress("probe_Bxysig", &probe_Bxysig, &b_probe_Bxysig);
  treeOrig->SetBranchAddress("probe_BmatchMC", &probe_BmatchMC, &b_probe_BmatchMC);
  treeOrig->SetBranchAddress("probe_Kpt", &probe_Kpt, &b_probe_Kpt);
  treeOrig->SetBranchAddress("probe_pt", &probe_pt, &b_probe_pt);
  treeOrig->SetBranchAddress("probe_eta", &probe_eta, &b_probe_eta);
  treeOrig->SetBranchAddress("probe_phi", &probe_phi, &b_probe_phi);
  treeOrig->SetBranchAddress("probe_isPF", &probe_isPF, &b_probe_isPF);
  treeOrig->SetBranchAddress("probe_isPFOverlap", &probe_isPFOverlap, &b_probe_isPFOverlap);
  treeOrig->SetBranchAddress("probe_isLowPt", &probe_isLowPt, &b_probe_isLowPt);
  treeOrig->SetBranchAddress("probe_mvaId", &probe_mvaId, &b_probe_mvaId);
  treeOrig->SetBranchAddress("probe_pfmvaId", &probe_pfmvaId, &b_probe_pfmvaId);
  treeOrig->SetBranchAddress("probe_dxySig", &probe_dxySig, &b_probe_dxySig);
  treeOrig->SetBranchAddress("probe_dzSig", &probe_dzSig, &b_probe_dzSig);
  treeOrig->SetBranchAddress("probe_fBrem", &probe_fBrem, &b_probe_fBrem);
  treeOrig->SetBranchAddress("probe_unBiased", &probe_unBiased, &b_probe_unBiased);
  treeOrig->SetBranchAddress("probe_ptBiased", &probe_ptBiased, &b_probe_ptBiased);
  treeOrig->SetBranchAddress("probe_normpt", &probe_normpt, &b_probe_normpt);
  treeOrig->SetBranchAddress("probe_dzTrg", &probe_dzTrg, &b_probe_dzTrg);
  treeOrig->SetBranchAddress("probe_iso04_rel", &probe_iso04_rel, &b_probe_iso04_rel);
  treeOrig->SetBranchAddress("probe_invMass", &probe_invMass, &b_probe_invMass);
  treeOrig->SetBranchAddress("probe_convveto", &probe_convveto, &b_probe_convveto);
  treeOrig->SetBranchAddress("probe_matchMc", &probe_matchMc, &b_probe_matchMc);
  treeOrig->SetBranchAddress("probe_matchMcFromJPsi", &probe_matchMcFromJPsi, &b_probe_matchMcFromJPsi);
  treeOrig->SetBranchAddress("selectedPairsSize", &selectedPairsSize, &b_selectedPairsSize);

  // New variables
  Float_t   probePt;
  Float_t   probeMvaId;
  Float_t   probePfmvaId;
  Float_t   probeDzTrg;
  Float_t   probeIso04Rel;
  Float_t   probeDxySig;
  Float_t   probeDzSig;

  // New branches
  for(int i=0; i<(int)trees.size();i++) {
    TTree *theTreeNew = trees[i];
    theTreeNew->Branch("probePt",&probePt,"probePt/F");
    theTreeNew->Branch("probeMvaId",&probeMvaId,"probeMvaId/F");
    theTreeNew->Branch("probePfmvaId",&probePfmvaId,"probePfmvaId/F");
    theTreeNew->Branch("probeDzTrg",&probeDzTrg,"probeDzTrg/F");
    theTreeNew->Branch("probeIso04Rel",&probeIso04Rel,"probeIso04Rel/F");
    theTreeNew->Branch("probeDxySig",&probeDxySig,"probeDxySig/F");
    theTreeNew->Branch("probeDzSig",&probeDzSig,"probeDzSig/F");
  }

  cout << "Now preparing the new tree" << endl;
  for(int i=0; i<nentriesOrig; i++) {
    
    if (i%10000 == 0) std::cout << ">>> Event # " << i << " / " << nentriesOrig << " entries" << std::endl; 
    treeOrig->GetEntry(i);

    // Trigger
    if (hlt9ip6==0) continue;

    // Loop over electrons
    for (unsigned int ii=0; ii<probe_Bmass->size(); ii++) {

      // conversion veto
      if ( tag_convveto->at(ii)==0)   continue;      
      if ( probe_convveto->at(ii)==0) continue;      

      // further selection on tag: low pt probes (loose selection)
      if (isProbeLpt==1) { 
	if ( tag_isPF->at(ii)==0 ) continue;
	if ( probe_isLowPt->at(ii)==0) continue;                   
	if ( tag_pfmvaId->at(ii)<-1. ) continue;        
      }

      // further selection on tag: low pt probes
      if (isProbeLpt==0) { 
	if ( probe_isPF->at(ii)==0) continue;  
	if ( tag_isLowPt->at(ii)==1 && tag_mvaId->at(ii)<-1) continue;
      }

      // Selection on probe (usually moved to analysis)
      if (fabs(probe_eta->at(ii))>2.4) continue;
      if (probe_pt->at(ii)<1.0)        continue;
      if (probe_pt->at(ii)>1000) continue;
      if (fabs(probe_dxySig->at(ii))>200) continue;   
      if (fabs(probe_dzTrg->at(ii))>10) continue;   
      if (probe_iso04_rel->at(ii)<-1 || probe_iso04_rel->at(ii)>5000) continue; 

      // e+e- invariant mass selection
      if (probe_invMass->at(ii)<2 || probe_invMass->at(ii)>4) continue;  

      // signal or background
      if (isSignal==1 && probe_matchMcFromJPsi->at(ii)==0) continue;
      if (isSignal==0 && probe_matchMcFromJPsi->at(ii)==1) continue;
      
      // save new variables, making flat tree
      probePt = probe_pt->at(ii);
      probeMvaId = probe_mvaId->at(ii);
      probePfmvaId = probe_pfmvaId->at(ii);
      probeDzTrg = probe_dzTrg->at(ii);
      probeIso04Rel = probe_iso04_rel->at(ii);
      probeDxySig = probe_dxySig->at(ii);
      probeDzSig = probe_dzSig->at(ii);

      treeNew->Fill();
    }
  }

  // new format
  treeNew->Write();
  fileNew->Close();
  fileNew->ls();
  
  fileOrig->cd();
  fileOrig->Close();  
}

