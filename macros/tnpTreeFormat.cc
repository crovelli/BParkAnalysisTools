#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TH1.h>
#include <iostream>
#include <vector>
#include <TRandom.h>

using namespace std;

void tnpTreeFormat(const char* filename, float lumiForW, int isProbeLpt) {

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
  vector<float>   *tag_pfRelIso = nullptr;
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
  vector<float>   *probe_pfRelIso = nullptr;
  vector<float>   *probe_trkRelIso = nullptr;
  vector<float>   *probe_fBrem = nullptr;
  vector<float>   *probe_unBiased = nullptr;
  vector<float>   *probe_ptBiased = nullptr;
  vector<bool>    *probe_isTag = nullptr;
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
  TBranch        *b_tag_pfRelIso;   //!
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
  TBranch        *b_probe_pfRelIso;   //!
  TBranch        *b_probe_trkRelIso;   //!
  TBranch        *b_probe_fBrem;   //!
  TBranch        *b_probe_unBiased;   //!
  TBranch        *b_probe_ptBiased;   //!
  TBranch        *b_probe_isTag;   //!
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
  treeOrig->SetBranchAddress("tag_pfRelIso", &tag_pfRelIso, &b_tag_pfRelIso);
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
  treeOrig->SetBranchAddress("probe_pfRelIso", &probe_pfRelIso, &b_probe_pfRelIso);
  treeOrig->SetBranchAddress("probe_trkRelIso", &probe_trkRelIso, &b_probe_trkRelIso);
  treeOrig->SetBranchAddress("probe_fBrem", &probe_fBrem, &b_probe_fBrem);
  treeOrig->SetBranchAddress("probe_unBiased", &probe_unBiased, &b_probe_unBiased);
  treeOrig->SetBranchAddress("probe_ptBiased", &probe_ptBiased, &b_probe_ptBiased);
  treeOrig->SetBranchAddress("probe_isTag", &probe_isTag, &b_probe_isTag);
  treeOrig->SetBranchAddress("probe_invMass", &probe_invMass, &b_probe_invMass);
  treeOrig->SetBranchAddress("probe_convveto", &probe_convveto, &b_probe_convveto);
  treeOrig->SetBranchAddress("probe_matchMc", &probe_matchMc, &b_probe_matchMc);
  treeOrig->SetBranchAddress("probe_matchMcFromJPsi", &probe_matchMcFromJPsi, &b_probe_matchMcFromJPsi);
  treeOrig->SetBranchAddress("selectedPairsSize", &selectedPairsSize, &b_selectedPairsSize);

  // New variables
  Int_t     hlt_9;
  Int_t     hlt_12;
  Int_t     numvtx;
  Int_t     tagMatchMc;
  Int_t     tagMatchMcFromJPsi;
  Float_t   tagPt;
  Float_t   tagEta;
  Float_t   tagRelIso;
  Float_t   tagPfmvaId;
  Float_t   probePt;
  Float_t   probeAbsEta;
  Float_t   probeEta;
  Int_t     probeIsPF;
  Int_t     probeIsPFOverlap;
  Int_t     probeIsLowPt;
  Float_t   probeMvaId;
  Float_t   probePfmvaId;
  Float_t   probeUnBiased;
  Float_t   probePtBiased;
  Int_t     probeConvVeto;
  Int_t     probeMatchMc;
  Int_t     probeMatchMcFromJPsi;
  Float_t   probeTrkRelIso;
  Float_t   probePFRelIso;
  Float_t   probeFBrem;
  Float_t   probeDxySig;
  Float_t   probeDzSig;
  Float_t   K_pt;
  Float_t   B_mass;
  Float_t   pair_mass;
  Float_t   weight;

  // New branches
  for(int i=0; i<(int)trees.size();i++) {
    TTree *theTreeNew = trees[i];
    theTreeNew->Branch("hlt_9", &hlt_9, "hlt_9/I");
    theTreeNew->Branch("hlt_12", &hlt_12, "hlt_12/I");
    theTreeNew->Branch("numvtx", &numvtx, "numvtx/I");
    theTreeNew->Branch("tagMatchMc",&tagMatchMc,"tagMatchMc/I");
    theTreeNew->Branch("tagMatchMcFromJPsi",&tagMatchMcFromJPsi,"tagMatchMcFromJPsi/I");
    theTreeNew->Branch("tagPt",&tagPt,"tagPt/F");
    theTreeNew->Branch("tagEta",&tagEta,"tagEta/F");
    theTreeNew->Branch("tagRelIso",&tagRelIso,"tagRelIso/F");
    theTreeNew->Branch("tagPfmvaId",&tagPfmvaId,"tagPfmvaId/F");
    theTreeNew->Branch("probePt",&probePt,"probePt/F");
    theTreeNew->Branch("probeAbsEta",&probeAbsEta,"probeAbsEta/F");
    theTreeNew->Branch("probeEta",&probeEta,"probeEta/F");
    theTreeNew->Branch("probeIsPF",&probeIsPF,"probeIsPF/I");
    theTreeNew->Branch("probeIsPFOverlap",&probeIsPFOverlap,"probeIsPFOverlap/I");
    theTreeNew->Branch("probeIsLowPt",&probeIsLowPt,"probeIsLowPt/I");
    theTreeNew->Branch("probeMvaId",&probeMvaId,"probeMvaId/F");
    theTreeNew->Branch("probePfmvaId",&probePfmvaId,"probePfmvaId/F");
    theTreeNew->Branch("probeUnBiased",&probeUnBiased,"probeUnBiased/F");
    theTreeNew->Branch("probePtBiased",&probePtBiased,"probePtBiased/F");
    theTreeNew->Branch("probeConvVeto",&probeConvVeto,"probeConvVeto/I");
    theTreeNew->Branch("probeMatchMc",&probeMatchMc,"probeMatchMc/I");
    theTreeNew->Branch("probeMatchMcFromJPsi",&probeMatchMcFromJPsi,"probeMatchMcFromJPsi/I");
    theTreeNew->Branch("probeTrkRelIso",&probeTrkRelIso,"probeTrkRelIso/F");
    theTreeNew->Branch("probePFRelIso",&probePFRelIso,"probePFRelIso/F");
    theTreeNew->Branch("probeFBrem",&probeFBrem,"probeFBrem/F");
    theTreeNew->Branch("probeDxySig",&probeDxySig,"probeDxySig/F");
    theTreeNew->Branch("probeDzSig",&probeDzSig,"probeDzSig/F");
    theTreeNew->Branch("K_pt",&K_pt,"K_pt/F");
    theTreeNew->Branch("B_mass", &B_mass, "B_mass/F");
    theTreeNew->Branch("pair_mass", &pair_mass, "pair_mass/F");
    theTreeNew->Branch("weight", &weight, "weight/F");
  }

  cout << "Now preparing the new tree" << endl;
  for(int i=0; i<nentriesOrig; i++) {
    
    if (i%10000 == 0) std::cout << ">>> Event # " << i << " / " << nentriesOrig << " entries" << std::endl; 
    treeOrig->GetEntry(i);

    for (unsigned int ii=0; ii<probe_Bmass->size(); ii++) {

      if ( tag_convveto->at(ii)==0)   continue;      
      if ( probe_convveto->at(ii)==0) continue;      

      // further selection on tag and probe
      if (isProbeLpt==1) { // probe=lowPT; tag=PF
	if ( tag_isPF->at(ii)==0 ) continue;
	if ( probe_isLowPt->at(ii)==0) continue;                   
	///////////if ( tag_pt->at(ii)<5 && tag_pfmvaId->at(ii)<1. )  continue;
	///////////if ( tag_pt->at(ii)>=5 && tag_pfmvaId->at(ii)<2. ) continue;
	if ( tag_pfmvaId->at(ii)<0. )  continue;
      }

      if (isProbeLpt==0) { // probe=PF; tag=PF or good lowPT
	if ( probe_isPF->at(ii)==0) continue;  
	if ( tag_isPF->at(ii)==0 && tag_mvaId->at(ii)<2) continue;
      }

      // e+e- invariant mass selection
      ////////////// if (probe_invMass->at(ii)<2 || probe_invMass->at(ii)>4) continue;  // comment for ROCs, uncomment for TnP

      // save new variables, making flat tree
      hlt_9  = hlt9;
      hlt_12 = hlt12;

      numvtx = nvtx;     

      B_mass = probe_Bmass->at(ii);
      pair_mass = probe_invMass->at(ii);
      
      tagMatchMc = tag_matchMc->at(ii);
      tagMatchMcFromJPsi = tag_matchMcFromJPsi->at(ii);
      tagPt  = tag_pt->at(ii);
      tagEta = tag_eta->at(ii);
      tagRelIso = tag_pfRelIso->at(ii);
      tagPfmvaId = tag_pfmvaId->at(ii);
      probePt = probe_pt->at(ii);
      probeAbsEta = fabs(probe_eta->at(ii));
      probeEta = probe_eta->at(ii);
      probeIsPF = probe_isPF->at(ii);
      probeIsPFOverlap = probe_isPFOverlap->at(ii);
      probeIsLowPt = probe_isLowPt->at(ii);
      probeMvaId = probe_mvaId->at(ii);
      probePfmvaId = probe_pfmvaId->at(ii);
      probeUnBiased = probe_unBiased->at(ii);
      probePtBiased = probe_ptBiased->at(ii);
      probeConvVeto = (int)probe_convveto->at(ii);
      probeMatchMc = probe_matchMc->at(ii);
      probeMatchMcFromJPsi = probe_matchMcFromJPsi->at(ii);
      probeTrkRelIso = probe_trkRelIso->at(ii);
      probePFRelIso = probe_pfRelIso->at(ii);
      probeFBrem = probe_fBrem->at(ii);
      probeDxySig = probe_dxySig->at(ii);
      probeDzSig = probe_dzSig->at(ii);
      K_pt = probe_Kpt->at(ii);

      // weights
      if (theRun==1) {   // MC                                                                                                                   
	//weight = perEveW * lumiForW * lumiWeight;     // chiara
	//weight = perEveW;
	weight = perEveW*pu_weight;
      } else {
	weight = 1.;
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

