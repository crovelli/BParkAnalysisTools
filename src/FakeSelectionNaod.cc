#include "TMath.h"
#include "TTree.h"
#include "TLorentzVector.h"

#include <iostream> 

#include "../include/FakeSelectionNaod.hh"    

#define MAX_PU_REWEIGHT 60

using namespace std;

FakeSelectionNaod::FakeSelectionNaod(TTree *tree)     
  : BParkBase(tree) {        

  // Chiara: to be set by hand   
  sampleID = 1;           // 0 = data, >=1 MC
  donvtxreweight_ = 0;    // 
  nvtxWFileName_  = "/afs/cern.ch/user/c/crovelli/public/bphys/nvtxWeights__bin1.root"; 
}

FakeSelectionNaod::~FakeSelectionNaod() {

  // output
  outFile_ -> cd();
  h_entries -> Write();
  h_selection -> Write();
  outTree_ -> Write();
  outFile_ -> Close();
}     

void FakeSelectionNaod::Loop() {

  if (fChain == 0) return;

  // Loop over events
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  cout << "entries : " <<  nentries << endl;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {

    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (jentry%10000==0) cout << jentry << endl;

    // To keep track of the total number of events 
    h_entries->Fill(5);
    
    // Event info     
    theRun = run;
    theEvent = event;
    theSampleID = sampleID;

    // # Vertices
    nvtx = PV_npvs;

    // Energy density
    rho = fixedGridRhoFastjetAll;
    
    // PU weight (for MC only and if requested)
    pu_weight = 1.;
    if (sampleID>0 && donvtxreweight_==1) {     // MC
      pu_weight = GetNvtxWeight(nvtx);         
    }

    // other weights for the dataset
    perEveW = 1.;
    if (sampleID>0) perEveW = Generator_weight;

    // Events breakdown  
    h_selection->Fill(0.,perEveW);
    

    // ----------------------------------------------------
    // PV must be there
    bool goodPV = false;
    if (PV_x>-999) goodPV = true;
    if (!goodPV) continue;
    h_selection->Fill(1.,perEveW);

    // Trigger 
    int iHLT_Mu12_IP6 = (int)HLT_Mu12_IP6;
    int iHLT_Mu9_IP6  = (int)HLT_Mu9_IP6;
    hlt9  = iHLT_Mu9_IP6;
    hlt12 = iHLT_Mu12_IP6;
    
    // B candidates
    if (nBToKMuMu<=0) continue;
    h_selection->Fill(2.,perEveW);

    // Minimal Bcandidate requirements
    vector<int> goodBs;
    for (u_int iB=0; iB<nBToKMuMu; iB++) {

      // preparing variables
      float mu1_eta = BToKMuMu_fit_l1_eta[iB];     
      float mu2_eta = BToKMuMu_fit_l2_eta[iB];    
      float k_pt    = BToKMuMu_fit_k_pt[iB];    
      float k_eta   = BToKMuMu_fit_k_eta[iB];    

      float b_xySig = BToKMuMu_l_xy[iB]/BToKMuMu_l_xy_unc[iB];

      // B selection (standard cut)
      // bool vtxFitSel = BToKMuMu_fit_pt[iB]>10.0 && b_xySig>6.0 && BToKMuMu_svprob[iB]>0.1 && BToKMuMu_fit_cos2D[iB]>0.999;
      // bool mu1Sel = fabs(mu1_eta)<2.4;  
      // bool mu2Sel = fabs(mu2_eta)<2.4;  
      // bool kSel   = k_pt>0.7 && fabs(k_eta)<2.4; 
      // bool additionalSel = BToKMuMu_fit_mass[iB]>4.5 && BToKMuMu_fit_mass[iB]<6.0;
      // bool isBsel = vtxFitSel && mu1Sel && mu2Sel && kSel && additionalSel;

      // B selection (relaxed)
      bool vtxFitSel = BToKMuMu_fit_pt[iB]>5.0 && BToKMuMu_svprob[iB]>0.1 && BToKMuMu_fit_cos2D[iB]>0.99;
      bool mu1Sel = fabs(mu1_eta)<2.4;  
      bool mu2Sel = fabs(mu2_eta)<2.4;  
      bool kSel   = k_pt>0.7 && fabs(k_eta)<2.4; 
      bool additionalSel = BToKMuMu_fit_mass[iB]>4.5 && BToKMuMu_fit_mass[iB]<6.0;
      bool isBsel = vtxFitSel && mu1Sel && mu2Sel && kSel && additionalSel;

      if (!isBsel) continue;
      goodBs.push_back(iB);
    }
    if (goodBs.size()>0) h_selection->Fill(3.,perEveW);

    
    // ------------------------------------------------------------------
    // Looking for electrons in B->Kmumu events
    for (u_int iB=0; iB<goodBs.size(); iB++) {

      int thisB = goodBs[iB];

      float my_mu1_pt  = BToKMuMu_fit_l1_pt[thisB]; 
      float my_mu1_eta = BToKMuMu_fit_l1_eta[thisB]; 
      float my_mu1_phi = BToKMuMu_fit_l1_phi[thisB]; 
      TLorentzVector mu1TLV(0,0,0,0);
      mu1TLV.SetPtEtaPhiM(my_mu1_pt,my_mu1_eta,my_mu1_phi,0);

      float my_mu2_pt  = BToKMuMu_fit_l2_pt[thisB]; 
      float my_mu2_eta = BToKMuMu_fit_l2_eta[thisB]; 
      float my_mu2_phi = BToKMuMu_fit_l2_phi[thisB]; 
      TLorentzVector mu2TLV(0,0,0,0);
      mu2TLV.SetPtEtaPhiM(my_mu2_pt,my_mu2_eta,my_mu2_phi,0);
  
      float my_k_pt  = BToKMuMu_fit_k_pt[thisB]; 
      float my_k_eta = BToKMuMu_fit_k_eta[thisB]; 
      float my_k_phi = BToKMuMu_fit_k_phi[thisB]; 
      TLorentzVector kTLV(0,0,0,0);
      kTLV.SetPtEtaPhiM(my_k_pt,my_k_eta,my_k_phi,0);

      int mu1_idx = BToKMuMu_l1Idx[thisB];
      int mu2_idx = BToKMuMu_l2Idx[thisB];

      // Loop over fake electrons
      for (u_int iEle=0; iEle<nElectron; iEle++) {
	
	float my_ele_pt  = Electron_pt[iEle];      
	float my_ele_eta = Electron_eta[iEle];   
	float my_ele_phi = Electron_phi[iEle];          
	TLorentzVector eleTLV(0,0,0,0);
	eleTLV.SetPtEtaPhiM(my_ele_pt,my_ele_eta,my_ele_phi,0);

	float my_deltaR_ele_mu1 = eleTLV.DeltaR(mu1TLV);
	float my_deltaR_ele_mu2 = eleTLV.DeltaR(mu2TLV);
	float my_deltaR_ele_k   = eleTLV.DeltaR(kTLV);

	// Filling tree

	// B-infos for further cuts offline
	Bmass.push_back(BToKMuMu_fit_mass[thisB]);
	Bpt.push_back(BToKMuMu_fit_pt[thisB]);
	Bcos2D.push_back(BToKMuMu_fit_cos2D[thisB]);
	Bsvprob.push_back(BToKMuMu_svprob[thisB]);
	Bxysig.push_back(BToKMuMu_l_xy[thisB]/BToKMuMu_l_xy_unc[thisB]);

	// Match with MC-truth
	bool isThisAMcB = -1;
	if (sampleID>0) isThisAMcB = isMcB(thisB);
	BmatchMC.push_back(isThisAMcB);

	bool isMu1McFromJPsi = -1;
	bool isMu2McFromJPsi = -1;
	if (sampleID>0) {
	  isMu1McFromJPsi = isMcMuFromJPsi(mu1_idx);
	  isMu2McFromJPsi = isMcMuFromJPsi(mu2_idx);
	}

	deltaR_ele_mu1.push_back(my_deltaR_ele_mu1);
	deltaR_ele_mu2.push_back(my_deltaR_ele_mu2);
	deltaR_ele_k.push_back(my_deltaR_ele_k);

	// K-infos for further cuts offline
	k_pt.push_back(my_k_pt);
	k_eta.push_back(my_k_eta);
	k_phi.push_back(my_k_phi);

	// Muons
	mu1_pt.push_back(my_mu1_pt);
	mu1_eta.push_back(my_mu1_eta);
	mu1_phi.push_back(my_mu1_phi);
	mu1_isTriggering.push_back(Muon_isTriggering[mu1_idx]);
	mu1_matchMcFromJPsi.push_back(isMu1McFromJPsi);
	//
	mu2_pt.push_back(my_mu2_pt);
	mu2_eta.push_back(my_mu2_eta);
	mu2_phi.push_back(my_mu2_phi);
	mu2_isTriggering.push_back(Muon_isTriggering[mu2_idx]);
	mu2_matchMcFromJPsi.push_back(isMu2McFromJPsi);	

	// Electrons
	ele_pt.push_back(Electron_pt[iEle]);
	ele_eta.push_back(Electron_eta[iEle]);
	ele_phi.push_back(Electron_phi[iEle]);
	ele_isPF.push_back(Electron_isPF[iEle]);
	ele_isPFOverlap.push_back(Electron_isPFoverlap[iEle]); 
	ele_isLowPt.push_back(Electron_isLowPt[iEle]); 
	ele_mvaId.push_back(Electron_mvaId[iEle]);
	ele_pfmvaId.push_back(Electron_pfmvaId[iEle]);
	ele_dxySig.push_back(Electron_dxy[iEle]/Electron_dxyErr[iEle]);  
	ele_dzSig.push_back(Electron_dz[iEle]/Electron_dzErr[iEle]);  
	ele_pfRelIso.push_back(Electron_pfRelIso[iEle]);  
	ele_trkRelIso.push_back(Electron_trkRelIso[iEle]);  
	ele_fBrem.push_back(Electron_fBrem[iEle]);  
	ele_unBiased.push_back(Electron_unBiased[iEle]);  
	ele_ptBiased.push_back(Electron_ptBiased[iEle]);  
	ele_convveto.push_back(Electron_convVeto[iEle]);

	// Probe track closest to electron 
	float minProbeDrEle = 999.;
	float probeCloseToEle_pt  = -999.;
	float probeCloseToEle_eta = -999.;
	float probeCloseToEle_phi = -999.;

	for (int ipt=0; ipt<nProbeTracks; ipt++) {
	  float pt_pt  = ProbeTracks_pt[ipt];
	  float pt_eta = ProbeTracks_eta[ipt];   
	  float pt_phi = ProbeTracks_phi[ipt];   
	  TLorentzVector ptTLV(0,0,0,0);
	  ptTLV.SetPtEtaPhiM(pt_pt,pt_eta,pt_phi,0);
	  float dRele = ptTLV.DeltaR(eleTLV);
	  if ( dRele<minProbeDrEle ) {
	    minProbeDrEle=dRele;
	    probeCloseToEle_pt=pt_pt;
	    probeCloseToEle_eta=pt_eta;
	    probeCloseToEle_phi=pt_phi;
	  } 
	}

	probe_closeToEle_pt.push_back(probeCloseToEle_pt);
	probe_closeToEle_eta.push_back(probeCloseToEle_eta);
	probe_closeToEle_phi.push_back(probeCloseToEle_phi);

      } // loop over electrons

    } // Loop over good Bs

    // At least one tag and one probe
    selectedElesSize = deltaR_ele_mu1.size();
    if (selectedElesSize<=0) continue;
    h_selection->Fill(4.,perEveW);

    // Filling the output tree
    outTree_->Fill();


    // Cleaning all vectors used for the output tree, ready for a new entry
    Bmass.clear();      
    Bpt.clear();
    Bcos2D.clear();
    Bsvprob.clear();
    Bxysig.clear();
    BmatchMC.clear();      
    deltaR_ele_mu1.clear();
    deltaR_ele_mu2.clear();
    deltaR_ele_k.clear();
    k_pt.clear();   
    k_eta.clear();   
    k_phi.clear();   
    mu1_pt.clear();   
    mu1_eta.clear();   
    mu1_phi.clear();   
    mu1_isTriggering.clear();   
    mu1_matchMcFromJPsi.clear();      
    mu2_pt.clear();   
    mu2_eta.clear();   
    mu2_phi.clear();   
    mu2_isTriggering.clear();   
    mu2_matchMcFromJPsi.clear();      
    ele_pt.clear();
    ele_eta.clear();
    ele_phi.clear();
    ele_isPF.clear();
    ele_isPFOverlap.clear();
    ele_isLowPt.clear();
    ele_mvaId.clear();
    ele_pfmvaId.clear();
    ele_dxySig.clear();
    ele_dzSig.clear();
    ele_pfRelIso.clear();
    ele_trkRelIso.clear();
    ele_fBrem.clear();
    ele_unBiased.clear();
    ele_ptBiased.clear();
    ele_convveto.clear();
    probe_closeToEle_pt.clear();
    probe_closeToEle_eta.clear();
    probe_closeToEle_phi.clear();

  } // loop over entries

}

void FakeSelectionNaod::SetNvtxWeights(std::string nvtxWeightFile) {

  if (nvtxWeightFile == "") {
    std::cout << "you need a weights file to use this function" << std::endl;
    return;
  }
  std::cout << "PU REWEIGHTING Based on #vertices:: Using file " << nvtxWeightFile << std::endl;
  
  TFile *f_nvtx = new TFile(nvtxWeightFile.c_str(),"READ");
  f_nvtx->cd();
  
  TH1F *nvtxweights = 0;
  TH1F *mc_nvtx = 0;
  mc_nvtx     = (TH1F*) f_nvtx->Get("mcNvtx");
  nvtxweights = (TH1F*) f_nvtx->Get("weights");
  
  if (!nvtxweights || !mc_nvtx) {
    std::cout << "weights histograms not found in file " << nvtxWeightFile << std::endl;
    return;
  }
  TH1F* weightedNvtx= (TH1F*)mc_nvtx->Clone("weightedNvtx");
  weightedNvtx->Multiply(nvtxweights);
  
  // Rescaling weights in order to preserve same integral of events     
  TH1F* weights = (TH1F*)nvtxweights->Clone("rescaledWeights");
  weights->Scale( mc_nvtx->Integral() / weightedNvtx->Integral() );
  
  float sumNvtxweights=0.;
  for (int i = 0; i<nvtxweights->GetNbinsX(); i++) {
    float weight=1.;
    weight=weights->GetBinContent(i+1);
    sumNvtxweights+=weight;
    nvtxweights_.push_back(weight);
    float lowedge=weights->GetBinLowEdge(i+1);
    nvtxlowedge_.push_back(lowedge);
  }
}

float FakeSelectionNaod::GetNvtxWeight(float nvtx) {

  int thesize   = nvtxlowedge_.size();
  int thesizem1 = nvtxlowedge_.size()-1;
  float weight=1;

  if (sampleID!=0 && thesize>0 && donvtxreweight_) {
    for (int i = 0; i<thesizem1; i++) {   
      if (nvtxlowedge_[i]<=nvtx && nvtxlowedge_[i+1]>nvtx) { 
	weight = nvtxweights_[i];
      }
    }
    if (nvtxlowedge_[thesizem1]<=nvtx) { 
      weight = nvtxweights_[thesizem1];
    }
  }

  return weight;
}

void FakeSelectionNaod::PrepareOutputs(std::string filename) 
{
  _datasetName=filename;

  std::string outname = _datasetName+".root";    
  cout << "output: " << outname << endl;
  outFile_ = new TFile(outname.c_str(),"RECREATE");

  bookOutputTree();
  bookOutputHistos();

  // loading weights for pileup if needed  
  if (donvtxreweight_) SetNvtxWeights(nvtxWFileName_);
};


void FakeSelectionNaod::bookOutputTree() 
{
  outTree_ = new TTree("TaPtree", "TaPtree");
  
  cout << "Booking tree" << endl;
  
  outTree_->Branch("theRun",   &theRun,   "theRun/I");    
  outTree_->Branch("theEvent", &theEvent, "theEvent/I");    
  outTree_->Branch("nvtx",     &nvtx,     "nvtx/I");    
  outTree_->Branch("sampleID", &sampleID, "sampleID/I");    
  outTree_->Branch("rho",      &rho,      "rho/F");    
  outTree_->Branch("pu_weight", &pu_weight, "pu_weight/F");
  outTree_->Branch("perEveW",  &perEveW,  "perEveW/F");    

  outTree_->Branch("hlt9",     &hlt9,     "hlt9/I");
  outTree_->Branch("hlt12",    &hlt12,    "hlt12/I");

  outTree_->Branch("selectedElesSize",  &selectedElesSize,  "selectedElesSize/I");  

  outTree_->Branch("Bmass",   "std::vector<float>", &Bmass);   
  outTree_->Branch("Bpt",     "std::vector<float>", &Bpt);   
  outTree_->Branch("Bcos2D",  "std::vector<float>", &Bcos2D);   
  outTree_->Branch("Bsvprob", "std::vector<float>", &Bsvprob);   
  outTree_->Branch("Bxysig",  "std::vector<float>", &Bxysig);   

  outTree_->Branch("BmatchMC", "std::vector<int>",   &BmatchMC);   

  outTree_->Branch("deltaR_ele_mu1",   "std::vector<float>",  &deltaR_ele_mu1);   
  outTree_->Branch("deltaR_ele_mu2",   "std::vector<float>",  &deltaR_ele_mu2);   
  outTree_->Branch("deltaR_ele_k",     "std::vector<float>",  &deltaR_ele_k);   

  outTree_->Branch("k_pt",             "std::vector<float>",  &k_pt);     
  outTree_->Branch("k_eta",            "std::vector<float>",  &k_eta);     
  outTree_->Branch("k_phi",            "std::vector<float>",  &k_phi);       

  outTree_->Branch("mu1_pt",              "std::vector<float>",  &mu1_pt);     
  outTree_->Branch("mu1_eta",             "std::vector<float>",  &mu1_eta);     
  outTree_->Branch("mu1_phi",             "std::vector<float>",  &mu1_phi);       
  outTree_->Branch("mu1_isTriggering",    "std::vector<int>",    &mu1_isTriggering);  
  outTree_->Branch("mu1_matchMcFromJPsi", "std::vector<int>",    &mu1_matchMcFromJPsi);    

  outTree_->Branch("mu2_pt",             "std::vector<float>",   &mu2_pt);     
  outTree_->Branch("mu2_eta",             "std::vector<float>",  &mu2_eta);     
  outTree_->Branch("mu2_phi",             "std::vector<float>",  &mu2_phi);       
  outTree_->Branch("mu2_isTriggering",    "std::vector<int>",    &mu2_isTriggering);
  outTree_->Branch("mu2_matchMcFromJPsi", "std::vector<int>",    &mu2_matchMcFromJPsi);     

  outTree_->Branch("ele_pt",            "std::vector<float>",  &ele_pt);     
  outTree_->Branch("ele_eta",           "std::vector<float>",  &ele_eta);     
  outTree_->Branch("ele_phi",           "std::vector<float>",  &ele_phi);     
  outTree_->Branch("ele_isPF",          "std::vector<bool>",   &ele_isPF);    
  outTree_->Branch("ele_isPFOverlap",   "std::vector<bool>",   &ele_isPFOverlap); 
  outTree_->Branch("ele_isLowPt",       "std::vector<bool>",   &ele_isLowPt);    
  outTree_->Branch("ele_mvaId",         "std::vector<float>",  &ele_mvaId);     
  outTree_->Branch("ele_pfmvaId",       "std::vector<float>",  &ele_pfmvaId);     
  outTree_->Branch("ele_dxySig",        "std::vector<float>",  &ele_dxySig);     
  outTree_->Branch("ele_dzSig",         "std::vector<float>",  &ele_dzSig);    
  outTree_->Branch("ele_pfRelIso",      "std::vector<float>",  &ele_pfRelIso);     
  outTree_->Branch("ele_trkRelIso",     "std::vector<float>",  &ele_trkRelIso);      
  outTree_->Branch("ele_fBrem",         "std::vector<float>",  &ele_fBrem);     
  outTree_->Branch("ele_unBiased",      "std::vector<float>",  &ele_unBiased);     
  outTree_->Branch("ele_ptBiased",      "std::vector<float>",  &ele_ptBiased);     
  outTree_->Branch("ele_convveto",      "std::vector<bool>",   &ele_convveto);   

  outTree_->Branch("probe_closeToEle_pt",  "std::vector<float>", &probe_closeToEle_pt);
  outTree_->Branch("probe_closeToEle_eta", "std::vector<float>", &probe_closeToEle_eta);
  outTree_->Branch("probe_closeToEle_phi", "std::vector<float>", &probe_closeToEle_phi);
}

void FakeSelectionNaod::bookOutputHistos() 
{
  cout << "Booking histos" << endl;
  //
  h_entries   = new TH1F("h_entries",  "Number of entries",   3,  3.5, 6.5);
  h_selection = new TH1F("h_selection","Selection breakdown", 8, -0.5, 7.5);
}

bool FakeSelectionNaod::isMcB( int theB ) {
  
  // taking index
  int mu1_idx = BToKMuMu_l1Idx[theB];
  int mu2_idx = BToKMuMu_l2Idx[theB];
  int k_idx   = BToKMuMu_kIdx[theB];

  // Gen tree
  int k_genPartIdx      = ProbeTracks_genPartIdx[k_idx];  
  int k_genMotherIdx    = GenPart_genPartIdxMother[k_genPartIdx];
  int k_genGMotherIdx   = GenPart_genPartIdxMother[k_genMotherIdx];
  int k_genPdgId        = GenPart_pdgId[k_genPartIdx];
  int k_genMotherPdgId  = GenPart_pdgId[k_genMotherIdx];
  int k_genGMotherPdgId = GenPart_pdgId[k_genGMotherIdx];

  int mu1_genPartIdx      = Muon_genPartIdx[mu1_idx];  
  int mu1_genMotherIdx    = GenPart_genPartIdxMother[mu1_genPartIdx];
  int mu1_genGMotherIdx   = GenPart_genPartIdxMother[mu1_genMotherIdx];
  int mu1_genPdgId        = GenPart_pdgId[mu1_genPartIdx];
  int mu1_genMotherPdgId  = GenPart_pdgId[mu1_genMotherIdx];
  int mu1_genGMotherPdgId = GenPart_pdgId[mu1_genGMotherIdx];

  int mu2_genPartIdx      = Muon_genPartIdx[mu2_idx];  
  int mu2_genMotherIdx    = GenPart_genPartIdxMother[mu2_genPartIdx];
  int mu2_genGMotherIdx   = GenPart_genPartIdxMother[mu2_genMotherIdx];
  int mu2_genPdgId        = GenPart_pdgId[mu2_genPartIdx];
  int mu2_genMotherPdgId  = GenPart_pdgId[mu2_genMotherIdx];
  int mu2_genGMotherPdgId = GenPart_pdgId[mu2_genGMotherIdx];

  // B -> K J/psi(ll) at gen level
  bool okMatch = (mu1_genPartIdx>-0.5 && mu2_genPartIdx>-0.5 && k_genPartIdx>-0.5);
  bool RK_res1 = abs(mu1_genMotherPdgId)==443 && abs(k_genMotherPdgId)==521;
  bool RK_res2 = (mu1_genMotherPdgId==mu2_genMotherPdgId) && (k_genMotherPdgId==mu1_genGMotherPdgId) && (k_genMotherPdgId==mu2_genGMotherPdgId);
  bool RK_res = okMatch && RK_res1 && RK_res2;

  return RK_res;
}

bool FakeSelectionNaod::isMcMuFromJPsi( int mu_idx ) {

  // Gen tree
  int mu_genPartIdx      = Muon_genPartIdx[mu_idx];  
  int mu_genMotherIdx    = GenPart_genPartIdxMother[mu_genPartIdx];
  int mu_genGMotherIdx   = GenPart_genPartIdxMother[mu_genMotherIdx];
  int mu_genPdgId        = GenPart_pdgId[mu_genPartIdx];
  int mu_genMotherPdgId  = GenPart_pdgId[mu_genMotherIdx];
  int mu_genGMotherPdgId = GenPart_pdgId[mu_genGMotherIdx];

  // B -> K J/psi(ll) at gen level
  bool okMatch = (mu_genPartIdx>-0.5) && (abs(mu_genMotherPdgId)==443);

  return okMatch;
}
