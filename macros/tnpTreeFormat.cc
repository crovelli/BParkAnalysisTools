#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TH1.h>
#include <iostream>
#include <vector>
#include <TRandom.h>

using namespace std;

void tnpTreeFormat(const char* filename, float lumiForW) {

  cout << "Formatting " << filename << endl;  

  // Original ntuple   
  TFile *fileOrig = 0;
  TTree *treeOrig = 0;

  fileOrig = TFile::Open(filename);
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
  Int_t           theLumi = 0;
  Int_t           nvtx = 0;
  Int_t           sampleID = 0;
  Float_t         rho = 0;
  Float_t         pu_weight = 0;
  Float_t         pu_n = 0;
  Float_t         perEveW = 0;
  Int_t           hlt9 = 0;
  Int_t           hlt12 = 0;
  Int_t           selectedBSize = 0;
  vector<float>   *tag_pt = nullptr;
  vector<float>   *tag_eta = nullptr;
  vector<float>   *tag_phi = nullptr;
  vector<bool>    *tag_isPF = nullptr;
  vector<bool>    *tag_isPFOverlap = nullptr;
  vector<bool>    *tag_isLowPt = nullptr;
  vector<float>   *tag_mvaId = nullptr;
  vector<float>   *tag_pfmvaId = nullptr;
  vector<float>   *tag_unBiased = nullptr;
  vector<float>   *tag_pfRelIso = nullptr;
  vector<bool>    *tag_matchMC = nullptr;
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
  vector<float>   *probe_mvaIdExtra = nullptr;
  vector<float>   *probe_pfmvaId = nullptr;
  vector<float>   *probe_dxySig = nullptr;
  vector<float>   *probe_dzSig = nullptr;
  vector<float>   *probe_pfRelIso = nullptr;
  vector<float>   *probe_trkRelIso = nullptr;
  vector<float>   *probe_fBrem = nullptr;
  vector<float>   *probe_unBiased = nullptr;
  vector<float>   *probe_ptBiased = nullptr;
  vector<bool>    *probe_isTag = nullptr;
  vector<float>   *probe_invMass = nullptr;
  vector<bool>    *probe_matchMC = nullptr;
  Int_t           selectedPairsSize = 0;
  
  // List of branches - original tree
  TBranch        *b_theRun;   //!
  TBranch        *b_theEvent;   //!
  TBranch        *b_theLumi;   //!
  TBranch        *b_nvtx;   //!
  TBranch        *b_sampleID;   //!
  TBranch        *b_rho;   //!
  TBranch        *b_pu_weight;   //!
  TBranch        *b_pu_n;   //!
  TBranch        *b_perEveW;   //!
  TBranch        *b_hlt9;   //!
  TBranch        *b_hlt12;   //!
  TBranch        *b_selectedBSize;   //!
  TBranch        *b_tag_pt;   //!
  TBranch        *b_tag_eta;   //!
  TBranch        *b_tag_phi;   //!
  TBranch        *b_tag_isPF;   //!
  TBranch        *b_tag_isPFOverlap;   //!
  TBranch        *b_tag_isLowPt;   //!
  TBranch        *b_tag_mvaId;   //!
  TBranch        *b_tag_pfmvaId;   //!
  TBranch        *b_tag_unBiased;   //!
  TBranch        *b_tag_pfRelIso;   //!
  TBranch        *b_tag_matchMC;   //!
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
  TBranch        *b_probe_mvaIdExtra;   //!
  TBranch        *b_probe_pfmvaId;   //!
  TBranch        *b_probe_dxySig;   //!
  TBranch        *b_probe_dzSig;   //!
  TBranch        *b_probe_pfRelIso;   //!
  TBranch        *b_probe_trkRelIso;   //!
  TBranch        *b_probe_fBrem;   //!
  TBranch        *b_probe_unBiased;   //!
  TBranch        *b_probe_ptBiased;   //!
  TBranch        *b_probe_isTag;   //!
  TBranch        *b_probe_invMass;   //!
  TBranch        *b_probe_matchMC;   //!
  TBranch        *b_selectedPairsSize;   //!
  
  // Set branch addresses and branch pointers 
  treeOrig->SetBranchAddress("theRun", &theRun, &b_theRun);
  treeOrig->SetBranchAddress("theEvent", &theEvent, &b_theEvent);
  treeOrig->SetBranchAddress("theLumi", &theLumi, &b_theLumi);
  treeOrig->SetBranchAddress("nvtx", &nvtx, &b_nvtx);
  treeOrig->SetBranchAddress("sampleID", &sampleID, &b_sampleID);
  treeOrig->SetBranchAddress("rho", &rho, &b_rho);
  treeOrig->SetBranchAddress("pu_weight", &pu_weight, &b_pu_weight);
  treeOrig->SetBranchAddress("pu_n", &pu_n, &b_pu_n);
  treeOrig->SetBranchAddress("perEveW", &perEveW, &b_perEveW);
  treeOrig->SetBranchAddress("hlt9", &hlt9, &b_hlt9);
  treeOrig->SetBranchAddress("hlt12", &hlt12, &b_hlt12);
  treeOrig->SetBranchAddress("selectedBSize", &selectedBSize, &b_selectedBSize);
  treeOrig->SetBranchAddress("tag_pt", &tag_pt, &b_tag_pt);
  treeOrig->SetBranchAddress("tag_eta", &tag_eta, &b_tag_eta);
  treeOrig->SetBranchAddress("tag_phi", &tag_phi, &b_tag_phi);
  treeOrig->SetBranchAddress("tag_isPF", &tag_isPF, &b_tag_isPF);
  treeOrig->SetBranchAddress("tag_isPFOverlap", &tag_isPFOverlap, &b_tag_isPFOverlap);
  treeOrig->SetBranchAddress("tag_isLowPt", &tag_isLowPt, &b_tag_isLowPt);
  treeOrig->SetBranchAddress("tag_mvaId", &tag_mvaId, &b_tag_mvaId);
  treeOrig->SetBranchAddress("tag_pfmvaId", &tag_pfmvaId, &b_tag_pfmvaId);
  treeOrig->SetBranchAddress("tag_unBiased", &tag_unBiased, &b_tag_unBiased);
  treeOrig->SetBranchAddress("tag_pfRelIso", &tag_pfRelIso, &b_tag_pfRelIso);
  treeOrig->SetBranchAddress("tag_matchMC", &tag_matchMC, &b_tag_matchMC);
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
  treeOrig->SetBranchAddress("probe_mvaIdExtra", &probe_mvaIdExtra, &b_probe_mvaIdExtra);
  treeOrig->SetBranchAddress("probe_pfmvaId", &probe_pfmvaId, &b_probe_pfmvaId);
  treeOrig->SetBranchAddress("probe_dxySig", &probe_dxySig, &b_probe_dxySig);
  treeOrig->SetBranchAddress("probe_dzSig", &probe_dzSig, &b_probe_dzSig);
  treeOrig->SetBranchAddress("probe_pfRelIso", &probe_pfRelIso, &b_probe_pfRelIso);
  treeOrig->SetBranchAddress("probe_trkRelIso", &probe_trkRelIso, &b_probe_trkRelIso);
  treeOrig->SetBranchAddress("probe_fBrem", &probe_fBrem, &b_probe_fBrem);
  treeOrig->SetBranchAddress("probe_unBiased", &probe_unBiased, &b_probe_unBiased);
  treeOrig->SetBranchAddress("probe_ptBiased", &probe_ptBiased, &b_probe_ptBiased);
  treeOrig->SetBranchAddress("probe_isTag", &probe_isTag, &b_probe_isTag);
  treeOrig->SetBranchAddress("probe_invMass", &probe_invMass, &b_probe_invMass);
  treeOrig->SetBranchAddress("probe_matchMC", &probe_matchMC, &b_probe_matchMC);
  treeOrig->SetBranchAddress("selectedPairsSize", &selectedPairsSize, &b_selectedPairsSize);

  // New variables
  Int_t     hlt_9;
  Int_t     hlt_12;
  Int_t     tagMatchMC;
  Float_t   probePt;
  Float_t   probeAbsEta;
  Float_t   probeEta;
  Int_t     probeIsPF;
  Int_t     probeIsLowPt;
  Float_t   probeMvaId;
  Float_t   probeMvaIdExtra;
  Float_t   probePfmvaId;
  Float_t   probeUnBiased;
  Int_t     probeMatchMC;
  Float_t   B_mass;
  Float_t   pair_mass;
  Float_t   weight;

  // New branches
  for(int i=0; i<(int)trees.size();i++) {
    TTree *theTreeNew = trees[i];
    theTreeNew->Branch("hlt_9", &hlt_9, "hlt_9/I");
    theTreeNew->Branch("hlt_12", &hlt_12, "hlt_12/I");
    theTreeNew->Branch("tagMatchMC",&tagMatchMC,"tagMatchMC/I");
    theTreeNew->Branch("probePt",&probePt,"probePt/F");
    theTreeNew->Branch("probeAbsEta",&probeAbsEta,"probeAbsEta/F");
    theTreeNew->Branch("probeEta",&probeEta,"probeEta/F");
    theTreeNew->Branch("probeIsPF",&probeIsPF,"probeIsPF/I");
    theTreeNew->Branch("probeIsLowPt",&probeIsLowPt,"probeIsLowPt/I");
    theTreeNew->Branch("probeMvaId",&probeMvaId,"probeMvaId/F");
    theTreeNew->Branch("probePfmvaId",&probePfmvaId,"probePfmvaId/F");
    theTreeNew->Branch("probeUnBiased",&probeUnBiased,"probeUnBiased/F");
    theTreeNew->Branch("probeMatchMC",&probeMatchMC,"probeMatchMC/I");
    theTreeNew->Branch("B_mass", &B_mass, "B_mass/F");
    theTreeNew->Branch("pair_mass", &pair_mass, "pair_mass/F");
    theTreeNew->Branch("weight", &weight, "weight/F");
  }

  cout << "Now preparing the new tree" << endl;
  for(int i=0; i<nentriesOrig; i++) {
    
    if (i%10000 == 0) std::cout << ">>> Event # " << i << " / " << nentriesOrig << " entries" << std::endl; 
    treeOrig->GetEntry(i);

    for (unsigned int ii=0; ii<probe_Bmass->size(); ii++) {

      // further selection on tag 
      if ( tag_isPF->at(ii)==0 ) continue;
      if ( tag_pt->at(ii)<5 && tag_pfmvaId->at(ii)<1. )  continue;
      if ( tag_pt->at(ii)>=5 && tag_pfmvaId->at(ii)<2. ) continue;
      // chiara: add conv. veto

      // further selection on probe 
      // chiara: add conv. veto

      // e+e- invariant mass selection
      if (probe_invMass->at(ii)<2 || probe_invMass->at(ii)>4) continue;

      // save new variables, making flat tree
      hlt_9  = hlt9;
      hlt_12 = hlt12;
      
      B_mass = probe_Bmass->at(ii);
      pair_mass = probe_invMass->at(ii);
      
      tagMatchMC = tag_matchMC->at(ii);
      
      probePt = probe_pt->at(ii);
      probeAbsEta = fabs(probe_eta->at(ii));
      probeEta = probe_eta->at(ii);
      probeIsPF = probe_isPF->at(ii);
      probeIsLowPt = probe_isLowPt->at(ii);
      probeMvaId = probe_mvaId->at(ii);
      probeMvaIdExtra = probe_mvaIdExtra->at(ii);
      probePfmvaId = probe_pfmvaId->at(ii);
      probeUnBiased = probe_unBiased->at(ii);
      probeMatchMC = probe_matchMC->at(ii);

      // weights
      if (theRun==1) {   // MC                                                                                                                   
	//weight = perEveW * lumiForW * lumiWeight;     // chiara
	weight = perEveW;
      } else {
	weight     = 1.;
      }

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

