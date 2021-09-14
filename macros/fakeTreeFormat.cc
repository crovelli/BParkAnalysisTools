#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include "TLorentzVector.h"
#include <iostream>
#include <vector>

using namespace std;

void fakeTreeFormat(const char* filename, float lumiForW, int isProbeLpt=-1) {

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
  Float_t         pu_weight;  
  Float_t         perEveW = 0;
  Int_t           hlt9 = 0;
  Int_t           selectedElesSize = 0; 
  vector<float>   *Bmass = nullptr;  
  vector<float>   *Bpt = nullptr;  
  vector<float>   *Bcos2D = nullptr;  
  vector<float>   *Bsvprob = nullptr;  
  vector<float>   *Bxysig = nullptr;  
  vector<int>     *BmatchMC = nullptr;  
  vector<float>   *deltaR_ele_mu1 = nullptr;
  vector<float>   *deltaR_ele_mu2 = nullptr;
  vector<float>   *deltaR_ele_k = nullptr;
  vector<float>   *k_pt = nullptr;
  vector<float>   *k_eta = nullptr;
  vector<float>   *k_phi = nullptr;
  vector<float>   *mu1_pt = nullptr;
  vector<float>   *mu1_eta = nullptr;
  vector<float>   *mu1_phi = nullptr;
  vector<int>     *mu1_isTriggering = nullptr;
  vector<int>     *mu1_matchMcFromJPsi = nullptr;
  vector<float>   *mu2_pt = nullptr;
  vector<float>   *mu2_eta = nullptr;
  vector<float>   *mu2_phi = nullptr;
  vector<int>     *mu2_isTriggering = nullptr;
  vector<int>     *mu2_matchMcFromJPsi = nullptr;
  vector<float>   *ele_pt = nullptr;
  vector<float>   *ele_eta = nullptr;
  vector<float>   *ele_phi = nullptr;
  vector<int>     *ele_isPF = nullptr;
  vector<int>     *ele_isPFOverlap = nullptr;
  vector<int>     *ele_isLowPt = nullptr;
  vector<float>   *ele_mvaId = nullptr;
  vector<float>   *ele_pfmvaId = nullptr;
  vector<float>   *ele_unBiased = nullptr;
  vector<float>   *ele_ptBiased = nullptr;
  vector<int>     *ele_convveto = nullptr;
  vector<float>   *probe_closeToEle_pt  = nullptr;
  vector<float>   *probe_closeToEle_eta = nullptr;
  vector<float>   *probe_closeToEle_phi = nullptr;

  // List of branches - original tree
  TBranch        *b_theRun;   //!
  TBranch        *b_theEvent;   //!
  TBranch        *b_nvtx;   //!
  TBranch        *b_sampleID;   //!
  TBranch        *b_rho;   //!
  TBranch        *b_pu_weight;   //!    
  TBranch        *b_perEveW;   //!
  TBranch        *b_hlt9;   //!
  TBranch        *b_selectedElesSize;   //!  
  TBranch        *b_Bmass; //!
  TBranch        *b_Bpt; //!
  TBranch        *b_Bcos2D; //!
  TBranch        *b_Bsvprob; //!
  TBranch        *b_Bxysig; //!
  TBranch        *b_BmatchMC; //!
  TBranch        *b_deltaR_ele_mu1; //!
  TBranch        *b_deltaR_ele_mu2; //!
  TBranch        *b_deltaR_ele_k; //!
  TBranch        *b_k_pt;   //!
  TBranch        *b_k_eta;   //!
  TBranch        *b_k_phi;   //!
  TBranch        *b_mu1_pt;   //!
  TBranch        *b_mu1_eta;   //!
  TBranch        *b_mu1_phi;   //!
  TBranch        *b_mu1_isTriggering;   //!
  TBranch        *b_mu1_matchMcFromJPsi;   //!
  TBranch        *b_mu2_pt;   //!
  TBranch        *b_mu2_eta;   //!
  TBranch        *b_mu2_phi;   //!
  TBranch        *b_mu2_isTriggering;   //!
  TBranch        *b_mu2_matchMcFromJPsi;   //!
  TBranch        *b_ele_pt;   //!
  TBranch        *b_ele_eta;   //!
  TBranch        *b_ele_phi;   //!
  TBranch        *b_ele_isPF;   //!
  TBranch        *b_ele_isPFOverlap;   //!
  TBranch        *b_ele_isLowPt;   //!
  TBranch        *b_ele_mvaId;   //!
  TBranch        *b_ele_pfmvaId;   //!
  TBranch        *b_ele_unBiased;   //!
  TBranch        *b_ele_ptBiased;   //!
  TBranch        *b_ele_convveto;   //!
  TBranch        *b_probe_closeToEle_pt;   //!
  TBranch        *b_probe_closeToEle_eta;   //!
  TBranch        *b_probe_closeToEle_phi;   //!

  // Set branch addresses and branch pointers 
  treeOrig->SetBranchAddress("theRun", &theRun, &b_theRun);
  treeOrig->SetBranchAddress("theEvent", &theEvent, &b_theEvent);
  treeOrig->SetBranchAddress("nvtx", &nvtx, &b_nvtx);
  treeOrig->SetBranchAddress("sampleID", &sampleID, &b_sampleID);
  treeOrig->SetBranchAddress("rho", &rho, &b_rho);
  treeOrig->SetBranchAddress("pu_weight", &pu_weight, &b_pu_weight);  
  treeOrig->SetBranchAddress("perEveW", &perEveW, &b_perEveW);
  treeOrig->SetBranchAddress("hlt9", &hlt9, &b_hlt9);
  treeOrig->SetBranchAddress("selectedElesSize", &selectedElesSize, &b_selectedElesSize);
  treeOrig->SetBranchAddress("Bmass", &Bmass, &b_Bmass);
  treeOrig->SetBranchAddress("Bpt", &Bpt, &b_Bpt);
  treeOrig->SetBranchAddress("Bcos2D", &Bcos2D, &b_Bcos2D);
  treeOrig->SetBranchAddress("Bsvprob", &Bsvprob, &b_Bsvprob);
  treeOrig->SetBranchAddress("Bxysig", &Bxysig, &b_Bxysig);
  treeOrig->SetBranchAddress("BmatchMC", &BmatchMC, &b_BmatchMC);
  treeOrig->SetBranchAddress("deltaR_ele_mu1", &deltaR_ele_mu1, &b_deltaR_ele_mu1);
  treeOrig->SetBranchAddress("deltaR_ele_mu2", &deltaR_ele_mu2, &b_deltaR_ele_mu2);
  treeOrig->SetBranchAddress("deltaR_ele_k", &deltaR_ele_k, &b_deltaR_ele_k);
  treeOrig->SetBranchAddress("k_pt", &k_pt, &b_k_pt);
  treeOrig->SetBranchAddress("k_eta", &k_eta, &b_k_eta);
  treeOrig->SetBranchAddress("k_phi", &k_phi, &b_k_phi);
  treeOrig->SetBranchAddress("mu1_pt", &mu1_pt, &b_mu1_pt);
  treeOrig->SetBranchAddress("mu1_eta", &mu1_eta, &b_mu1_eta);
  treeOrig->SetBranchAddress("mu1_phi", &mu1_phi, &b_mu1_phi);
  treeOrig->SetBranchAddress("mu1_isTriggering", &mu1_isTriggering, &b_mu1_isTriggering);
  treeOrig->SetBranchAddress("mu1_matchMcFromJPsi", &mu1_matchMcFromJPsi, &b_mu1_matchMcFromJPsi);
  treeOrig->SetBranchAddress("mu2_pt", &mu2_pt, &b_mu2_pt);
  treeOrig->SetBranchAddress("mu2_eta", &mu2_eta, &b_mu2_eta);
  treeOrig->SetBranchAddress("mu2_phi", &mu2_phi, &b_mu2_phi);
  treeOrig->SetBranchAddress("mu2_isTriggering", &mu2_isTriggering, &b_mu2_isTriggering);
  treeOrig->SetBranchAddress("mu2_matchMcFromJPsi", &mu2_matchMcFromJPsi, &b_mu2_matchMcFromJPsi);
  treeOrig->SetBranchAddress("ele_pt", &ele_pt, &b_ele_pt);
  treeOrig->SetBranchAddress("ele_eta", &ele_eta, &b_ele_eta);
  treeOrig->SetBranchAddress("ele_phi", &ele_phi, &b_ele_phi);
  treeOrig->SetBranchAddress("ele_isPF", &ele_isPF, &b_ele_isPF);
  treeOrig->SetBranchAddress("ele_isPFOverlap", &ele_isPFOverlap, &b_ele_isPFOverlap);
  treeOrig->SetBranchAddress("ele_isLowPt", &ele_isLowPt, &b_ele_isLowPt);
  treeOrig->SetBranchAddress("ele_mvaId", &ele_mvaId, &b_ele_mvaId);
  treeOrig->SetBranchAddress("ele_pfmvaId", &ele_pfmvaId, &b_ele_pfmvaId);
  treeOrig->SetBranchAddress("ele_unBiased", &ele_unBiased, &b_ele_unBiased);
  treeOrig->SetBranchAddress("ele_ptBiased", &ele_ptBiased, &b_ele_ptBiased);
  treeOrig->SetBranchAddress("ele_convveto", &ele_convveto, &b_ele_convveto);
  treeOrig->SetBranchAddress("probe_closeToEle_pt",  &probe_closeToEle_pt,  &b_probe_closeToEle_pt);
  treeOrig->SetBranchAddress("probe_closeToEle_eta", &probe_closeToEle_eta, &b_probe_closeToEle_eta);
  treeOrig->SetBranchAddress("probe_closeToEle_phi", &probe_closeToEle_phi, &b_probe_closeToEle_phi);

  // New variables
  Int_t     hlt_9;
  Int_t     numvtx;
  Float_t   bMass;
  Float_t   bPt;
  Float_t   bCos2D;
  Float_t   bSvprob;
  Float_t   bXysig;
  Int_t     bMatchMc;
  Float_t   mumuMass;
  Float_t   dR_ele_mu1;
  Float_t   dR_ele_mu2;
  Float_t   dR_ele_k;
  Float_t   dR_ele_mu1mu2;
  Float_t   kPt;
  Float_t   kEta;
  Float_t   kPhi;
  Float_t   mu1Pt;
  Float_t   mu1Eta;
  Float_t   mu1Phi;
  Int_t     mu1IsTriggering;
  Int_t     mu1MatchMcFromJPsi;
  Float_t   mu2Pt;
  Float_t   mu2Eta;
  Float_t   mu2Phi;
  Int_t     mu2IsTriggering;
  Int_t     mu2MatchMcFromJPsi;
  Float_t   elePt;
  Float_t   eleEta;
  Float_t   elePhi;
  Int_t     eleIsPF;
  Int_t     eleIsPFOverlap;
  Int_t     eleIsLowPt;
  Float_t   eleMvaId;
  Float_t   elePfmvaId;
  Int_t     eleConvVeto;
  Float_t   eleTrkRelIso;
  Float_t   elePFRelIso;
  Float_t   eleFBrem;
  Float_t   eleDxySig;
  Float_t   eleDzSig;
  Float_t   weight;
  Float_t   probeCloseToElePt;
  Float_t   probeCloseToEleEta;
  Float_t   probeCloseToElePhi;

  // New branches
  for(int i=0; i<(int)trees.size();i++) {
    TTree *theTreeNew = trees[i];
    theTreeNew->Branch("hlt_9",  &hlt_9,  "hlt_9/I");
    theTreeNew->Branch("numvtx", &numvtx, "numvtx/I");
    theTreeNew->Branch("bMass", &bMass, "bMass/F");
    theTreeNew->Branch("bPt", &bPt, "bPt/F");
    theTreeNew->Branch("bCos2D", &bCos2D, "bCos2D/F");
    theTreeNew->Branch("bSvprob", &bSvprob, "bSvprob/F");
    theTreeNew->Branch("bXysig", &bXysig, "bXysig/F");
    theTreeNew->Branch("bMatchMc", &bMatchMc, "bMatchMc/I");
    theTreeNew->Branch("mumuMass", &mumuMass, "mumuMass/F");
    theTreeNew->Branch("dR_ele_mu1", &dR_ele_mu1, "dR_ele_mu1/F");
    theTreeNew->Branch("dR_ele_mu2", &dR_ele_mu2, "dR_ele_mu2/F");
    theTreeNew->Branch("dR_ele_k",   &dR_ele_k,   "dR_ele_k/F");
    theTreeNew->Branch("dR_ele_mu1mu2", &dR_ele_mu1mu2, "dR_ele_mu1mu2/F");  
    theTreeNew->Branch("kPt",&kPt,"kPt/F");
    theTreeNew->Branch("kEta",&kEta,"kEta/F");
    theTreeNew->Branch("kPhi",&kPhi,"kPhi/F");
    theTreeNew->Branch("mu1Pt",&mu1Pt,"mu1Pt/F");
    theTreeNew->Branch("mu1Eta",&mu1Eta,"mu1Eta/F");
    theTreeNew->Branch("mu1Phi",&mu1Phi,"mu1Phi/F");
    theTreeNew->Branch("mu1IsTriggering",&mu1IsTriggering,"mu1IsTriggering/I");
    theTreeNew->Branch("mu1MatchMcFromJPsi",&mu1MatchMcFromJPsi,"mu1MatchMcFromJPsi/I");
    theTreeNew->Branch("mu2Pt",&mu2Pt,"mu2Pt/F");
    theTreeNew->Branch("mu2Eta",&mu2Eta,"mu2Eta/F");
    theTreeNew->Branch("mu2Phi",&mu2Phi,"mu2Phi/F");
    theTreeNew->Branch("mu2IsTriggering",&mu2IsTriggering,"mu2IsTriggering/I");
    theTreeNew->Branch("mu2MatchMcFromJPsi",&mu2MatchMcFromJPsi,"mu2MatchMcFromJPsi/I");
    theTreeNew->Branch("elePt",&elePt,"elePt/F");
    theTreeNew->Branch("eleEta",&eleEta,"eleEta/F");
    theTreeNew->Branch("elePhi",&elePhi,"elePhi/F");
    theTreeNew->Branch("eleIsPF",&eleIsPF,"eleIsPF/I");
    theTreeNew->Branch("eleIsPFOverlap",&eleIsPFOverlap,"eleIsPFOverlap/I");
    theTreeNew->Branch("eleIsLowPt",&eleIsLowPt,"eleIsLowPt/I");
    theTreeNew->Branch("eleMvaId",&eleMvaId,"eleMvaId/F");
    theTreeNew->Branch("elePfmvaId",&elePfmvaId,"elePfmvaId/F");
    theTreeNew->Branch("eleConvVeto",&eleConvVeto,"eleConvVeto/I");
    theTreeNew->Branch("eleTrkRelIso",&eleTrkRelIso,"eleTrkRelIso/F");
    theTreeNew->Branch("elePFRelIso",&elePFRelIso,"elePFRelIso/F");
    theTreeNew->Branch("eleFBrem",&eleFBrem,"eleFBrem/F");
    theTreeNew->Branch("eleDxySig",&eleDxySig,"eleDxySig/F");
    theTreeNew->Branch("eleDzSig",&eleDzSig,"eleDzSig/F");
    theTreeNew->Branch("weight", &weight, "weight/F");
    theTreeNew->Branch("probeCloseToElePt",  &probeCloseToElePt,  "probeCloseToElePt/F");
    theTreeNew->Branch("probeCloseToEleEta", &probeCloseToEleEta, "probeCloseToEleEta/F");
    theTreeNew->Branch("probeCloseToElePhi", &probeCloseToElePhi, "probeCloseToElePhi/F");
  }

  cout << "Now preparing the new tree" << endl;
  for(int i=0; i<nentriesOrig; i++) {
    
    if (i%10000 == 0) std::cout << ">>> Event # " << i << " / " << nentriesOrig << " entries" << std::endl; 
    treeOrig->GetEntry(i);

    // further selection: trigger 
    if (hlt9==0) continue; 

    // Loop on electrons
    for (unsigned int ii=0; ii<ele_pt->size(); ii++) {
      
      // further selection: remove tag side
      if (mu1_isTriggering->at(ii) || mu2_isTriggering->at(ii)) continue;

      // Build cone around mu-mu
      TLorentzVector mu1TLV(0,0,0,0);     
      TLorentzVector mu2TLV(0,0,0,0);     
      mu1TLV.SetPtEtaPhiM(mu1_pt->at(ii), mu1_eta->at(ii), mu1_phi->at(ii), 0);
      mu2TLV.SetPtEtaPhiM(mu2_pt->at(ii), mu2_eta->at(ii), mu2_phi->at(ii), 0);
      TLorentzVector mumuTLV = mu1TLV + mu2TLV;
      float theMumuMass = mumuTLV.M(); 
      mumuMass = theMumuMass;

      // Ele: kinematics
      if (ele_convveto->at(ii)==0)   continue;
      if (ele_pt->at(ii)<0.5)        continue;
      if (fabs(ele_eta->at(ii))>2.4) continue;

      // LowPt or PF
      if (isProbeLpt==1 && ele_isLowPt->at(ii)==0) continue;
      if (isProbeLpt==0 && ele_isPF->at(ii)==0)    continue;

      // Tight selection
      // further selection on B: smaller mass range 
      if (Bmass->at(ii)<4.7 || Bmass->at(ii)>5.6) continue; 

      // J/Psi region
      if (mumuMass>3.4 || mumuMass<2.8) continue;

      // Check dR ele-mumu
      TLorentzVector eleTLV(0,0,0,0);
      eleTLV.SetPtEtaPhiM(ele_pt->at(ii), ele_eta->at(ii), ele_phi->at(ii), 0);
      float deltaR_ele_mu1mu2 = eleTLV.DeltaR(mumuTLV);

      // save new variables, making flat tree
      hlt_9  = hlt9;
      numvtx = nvtx;     

      bMass = Bmass->at(ii); 
      bPt = Bpt->at(ii); 
      bCos2D = Bcos2D->at(ii); 
      bSvprob = Bsvprob->at(ii);  
      bXysig = Bxysig->at(ii);   
      bMatchMc = BmatchMC->at(ii); 

      dR_ele_mu1    = deltaR_ele_mu1->at(ii); 
      dR_ele_mu2    = deltaR_ele_mu2->at(ii); 
      dR_ele_k      = deltaR_ele_k->at(ii); 
      dR_ele_mu1mu2 = deltaR_ele_mu1mu2; 

      kPt  = k_pt->at(ii);
      kEta = k_eta->at(ii);
      kPhi = k_phi->at(ii);

      mu1Pt  = mu1_pt->at(ii);
      mu1Eta = mu1_eta->at(ii);
      mu1Phi = mu1_phi->at(ii);
      mu1IsTriggering = mu1_isTriggering->at(ii);
      mu1MatchMcFromJPsi = mu1_matchMcFromJPsi->at(ii);
      mu2Pt  = mu2_pt->at(ii);
      mu2Eta = mu2_eta->at(ii);
      mu2Phi = mu2_phi->at(ii);
      mu2IsTriggering = mu2_isTriggering->at(ii);
      mu2MatchMcFromJPsi = mu2_matchMcFromJPsi->at(ii);

      elePt = ele_pt->at(ii);
      eleEta = ele_eta->at(ii);
      elePhi = ele_phi->at(ii);
      eleIsPF = ele_isPF->at(ii);
      eleIsLowPt = ele_isLowPt->at(ii);
      eleIsPFOverlap = ele_isPFOverlap->at(ii);
      eleMvaId = ele_mvaId->at(ii);
      elePfmvaId = ele_pfmvaId->at(ii);
      eleConvVeto = (int)ele_convveto->at(ii);

      probeCloseToElePt  = probe_closeToEle_pt->at(ii);
      probeCloseToEleEta = probe_closeToEle_eta->at(ii);
      probeCloseToElePhi = probe_closeToEle_phi->at(ii);

      // weights
      if (theRun==1) {   // MC                                                                           
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

