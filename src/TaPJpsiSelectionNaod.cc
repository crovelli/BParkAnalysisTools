#include "TMath.h"
#include "TTree.h"
#include "TLorentzVector.h"

#include <iostream> 

#include "../include/TaPJpsiSelectionNaod.hh"    

using namespace std;

TaPJpsiSelectionNaod::TaPJpsiSelectionNaod(TTree *tree)     
  : BParkBase(tree) {        

  // Chiara: to be set by hand   
  sampleID = 1;           // 0 = data, >=1 MC
  donvtxreweight_ = 0;    // 
  nvtxWFileName_ = "/afs/cern.ch/user/c/crovelli/public/bphys/march/nvtxWeights_run2018D.root"; 
}

TaPJpsiSelectionNaod::~TaPJpsiSelectionNaod() {

  // output
  outFile_->cd();
  h_entries -> Write();
  h_selection -> Write();
  outTree_->Write();
  outFile_->Close();
}     

void TaPJpsiSelectionNaod::Loop() {

  if (fChain == 0) return;

  // Loop over events
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  cout << "entries : " <<  nentries << endl;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {

    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (jentry%50000==0) cout << jentry << endl;

    // To keep track of the total number of events 
    h_entries->Fill(5);
    
    // Event info     
    theRun   = run;
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
    int iHLT_Mu9_IP5  = (int)HLT_Mu9_IP5;
    int iHLT_Mu9_IP4  = (int)HLT_Mu9_IP4;
    int iHLT_Mu7_IP4  = (int)HLT_Mu7_IP4;
    int iHLT_Mu8_IP6  = (int)HLT_Mu8_IP6;
    int iHLT_Mu8_IP5  = (int)HLT_Mu8_IP5;
    int iHLT_Mu8_IP3  = (int)HLT_Mu8_IP3;
    int iHLT_Mu8p5_IP3p5  = (int)HLT_Mu8p5_IP3p5;
    int iHLT_Mu10p5_IP3p5 = (int)HLT_Mu10p5_IP3p5;
    hlt12ip6 = iHLT_Mu12_IP6;
    hlt9ip6  = iHLT_Mu9_IP6;
    hlt9ip5  = iHLT_Mu9_IP5;
    hlt9ip4  = iHLT_Mu9_IP4;
    hlt7ip4  = iHLT_Mu7_IP4;
    hlt8ip6  = iHLT_Mu8_IP6;
    hlt8ip5  = iHLT_Mu8_IP5;
    hlt8ip3  = iHLT_Mu8_IP3;
    hlt8d5ip3d5  = iHLT_Mu8p5_IP3p5;
    hlt10d5ip3d5 = iHLT_Mu10p5_IP3p5;
    
    // B candidates
    if (nBToKEE<=0) continue;
    h_selection->Fill(2.,perEveW);


    // Minimal Bcandidate requirements
    vector<int> goodBs;
    for (u_int iB=0; iB<nBToKEE; iB++) {

      // preparing variables
      int ele1_idx = BToKEE_l1Idx[iB];
      int ele2_idx = BToKEE_l2Idx[iB];

      float ele1_pt = BToKEE_fit_l1_pt[iB];
      float ele2_pt = BToKEE_fit_l2_pt[iB];
      float k_pt    = BToKEE_fit_k_pt[iB];

      float ele1_eta = BToKEE_fit_l1_eta[iB];     
      float ele2_eta = BToKEE_fit_l2_eta[iB];    
      float k_eta    = BToKEE_fit_k_eta[iB];    
      float ele1_phi = BToKEE_fit_l1_phi[iB];    
      float ele2_phi = BToKEE_fit_l2_phi[iB];    

      bool ele1_convveto = Electron_convVeto[ele1_idx];     
      bool ele2_convveto = Electron_convVeto[ele2_idx];     

      float b_xySig = BToKEE_l_xy[iB]/BToKEE_l_xy_unc[iB];

      // B selection (standard cut)
      // bool vtxFitSel = BToKEE_fit_pt[iB]>10.0 && b_xySig>6.0 && BToKEE_svprob[iB]>0.1 && BToKEE_fit_cos2D[iB]>0.999;
      // bool ele1Sel = ele1_pt>1.5 && fabs(ele1_eta)<2.4;  
      // bool ele2Sel = ele2_pt>0.5 && fabs(ele2_eta)<2.4;  
      // bool kSel = k_pt>0.7 && fabs(k_eta)<2.4; 
      
      // B selection (relaxed) 
      bool vtxFitSel = BToKEE_fit_pt[iB]>5.0 && BToKEE_svprob[iB]>0.1 && BToKEE_fit_cos2D[iB]>0.99;
      bool ele1Sel = fabs(ele1_eta)<2.4;  
      bool ele2Sel = fabs(ele2_eta)<2.4;  
      bool kSel = k_pt>0.7 && fabs(k_eta)<2.4; 

      bool additionalSel = BToKEE_fit_mass[iB]>4.5 && BToKEE_fit_mass[iB]<6.0;
      bool isBsel = vtxFitSel && ele1Sel && ele2Sel && kSel && additionalSel;

      if (!isBsel) continue;

      goodBs.push_back(iB);
    }
    if (goodBs.size()>0) h_selection->Fill(3.,perEveW);

    
    // ------------------------------------------------------------------
    // Selecting just one B among those with ee in common
    vector<int> cleanGoodBs;
    if (goodBs.size()>0) cleanGoodBs.push_back(goodBs[0]);

    for (u_int iB=1; iB<goodBs.size(); iB++) {
      int thisB    = goodBs[iB];
      int ele1_idx = BToKEE_l1Idx[thisB];
      int ele2_idx = BToKEE_l2Idx[thisB];
      float theXySig = BToKEE_l_xy[thisB]/BToKEE_l_xy_unc[thisB];

      int toadd = -1;
      vector<int> toBeErased;
      for (u_int iBC=0; iBC<cleanGoodBs.size(); iBC++) {
	int thisBC    = cleanGoodBs[iBC];
	int ele1_idxC = BToKEE_l1Idx[thisBC];
	int ele2_idxC = BToKEE_l2Idx[thisBC];
	float theXySigC = BToKEE_l_xy[thisBC]/BToKEE_l_xy_unc[thisBC];
	
	// already there - eventually to be replaced
	if ( ele1_idxC==ele1_idx && ele2_idxC==ele2_idx && theXySig>theXySigC ){
	  toBeErased.push_back(iBC);
	  cleanGoodBs.push_back(thisB);
	  // this B is already in, not to be added
	  toadd = -100;      

	} else if ( ele1_idxC==ele1_idx && ele2_idxC==ele2_idx && theXySig<=theXySigC ) {  // there's already a better B candidate with same leptons -> not to be added
	  toadd = -100;

	} else if ( (ele1_idxC!=ele1_idx || ele2_idxC!=ele2_idx) && toadd>-50 ) {   // new leptons and B not yet added -> to be added
	  toadd = thisB; 
	}
      }
      if (toadd>-1) cleanGoodBs.push_back(thisB);
    
      for (u_int iBTBE=0; iBTBE<toBeErased.size(); iBTBE++) {
	int thisTBE = toBeErased[iBTBE];
	if (thisTBE==0) cleanGoodBs.erase(cleanGoodBs.begin());  
	if (thisTBE>0)  cleanGoodBs.erase(cleanGoodBs.begin()+thisTBE);  
      }	
      toBeErased.clear();
    }
    goodBs.clear();

    selectedBSize = cleanGoodBs.size();


    // ------------------------------------------------------------------
    // At least one good B candidate
    if (cleanGoodBs.size()<=0) continue;
    h_selection->Fill(4.,perEveW);


    // ------------------------------------------------------------------
    // Tag and probe using all the J/Psi candidates and the 2 possible combinations
    for (u_int iB=0; iB<cleanGoodBs.size(); iB++) {

      // B-infos for further cuts offline
      int thisB = cleanGoodBs[iB];
      float thisBmass   = BToKEE_fit_mass[thisB];
      float thisBpt     = BToKEE_fit_pt[thisB];
      float thisBcos    = BToKEE_fit_cos2D[thisB];
      float thisBsvprob = BToKEE_svprob[thisB];
      float thisBxysig  = BToKEE_l_xy[thisB]/BToKEE_l_xy_unc[thisB];
      bool isThisAMcB = -1;
      if (sampleID>0) isThisAMcB = isMcB(thisB);      

      // K-infos for further cuts offline
      int thisK = BToKEE_kIdx[thisB];
      float thisKpt = BToKEE_fit_k_pt[thisB];   

      // Electrons
      int ele1_idx = BToKEE_l1Idx[thisB];
      int ele2_idx = BToKEE_l2Idx[thisB];

      float ele1_pt  = BToKEE_fit_l1_pt[thisB];    
      float ele2_pt  = BToKEE_fit_l2_pt[thisB];      
      float ele1_eta = BToKEE_fit_l1_eta[thisB]; 
      float ele2_eta = BToKEE_fit_l2_eta[thisB];     
      float ele1_phi = BToKEE_fit_l1_phi[thisB];     
      float ele2_phi = BToKEE_fit_l2_phi[thisB]; 

      TLorentzVector ele1TLV(0,0,0,0);
      ele1TLV.SetPtEtaPhiM(ele1_pt,ele1_eta,ele1_phi,0);
      TLorentzVector ele2TLV(0,0,0,0);
      ele2TLV.SetPtEtaPhiM(ele2_pt,ele2_eta,ele2_phi,0);
      float mee = (ele1TLV+ele2TLV).M();

      if ( 1 ) { 
	//
	tag_pt.push_back(ele1_pt);
	tag_eta.push_back(ele1_eta);
	tag_phi.push_back(ele1_phi);
	tag_isPF.push_back(Electron_isPF[ele1_idx]);           
	tag_isPFOverlap.push_back(Electron_isPFoverlap[ele1_idx]); 
	tag_isLowPt.push_back(Electron_isLowPt[ele1_idx]);
	tag_mvaId.push_back(Electron_mvaId[ele1_idx]);
	tag_pfmvaId.push_back(Electron_pfmvaId[ele1_idx]);
	tag_convveto.push_back(Electron_convVeto[ele1_idx]);
	tag_pfRelIso.push_back(Electron_pfRelIso[ele1_idx]);  
	//
	probe_pt.push_back(ele2_pt);
	probe_eta.push_back(ele2_eta);
	probe_phi.push_back(ele2_phi);
	probe_isPF.push_back(Electron_isPF[ele2_idx]);
	probe_isPFOverlap.push_back(Electron_isPFoverlap[ele2_idx]); 
	probe_isLowPt.push_back(Electron_isLowPt[ele2_idx]); 
	probe_mvaId.push_back(Electron_mvaId[ele2_idx]);
	probe_pfmvaId.push_back(Electron_pfmvaId[ele2_idx]);
	probe_dxySig.push_back(Electron_dxy[ele2_idx]/Electron_dxyErr[ele2_idx]);  
	probe_dzSig.push_back(Electron_dz[ele2_idx]/Electron_dzErr[ele2_idx]);  
	probe_pfRelIso.push_back(Electron_pfRelIso[ele2_idx]);  
	probe_trkRelIso.push_back(Electron_trkRelIso[ele2_idx]);  
	probe_fBrem.push_back(Electron_fBrem[ele2_idx]);  
	probe_unBiased.push_back(Electron_unBiased[ele2_idx]);  
	probe_ptBiased.push_back(Electron_ptBiased[ele2_idx]);  
	probe_convveto.push_back(Electron_convVeto[ele2_idx]);
	// 
	probe_Bmass.push_back(thisBmass);
	probe_Bpt.push_back(thisBpt);
	probe_Bcos2D.push_back(thisBcos);
	probe_Bsvprob.push_back(thisBsvprob);
	probe_Bxysig.push_back(thisBxysig);
	probe_BmatchMC.push_back(isThisAMcB);
	probe_Kpt.push_back(thisKpt);
	//
	probe_invMass.push_back(mee);  
	//
	if (sampleID>0) {     // MC      
	  bool isTagAMcEleFromJPsi   = isMcEleFromJPsi(ele1_idx);
	  bool isProbeAMcEleFromJPsi = isMcEleFromJPsi(ele2_idx);
	  tag_matchMcFromJPsi.push_back(isTagAMcEleFromJPsi);
	  probe_matchMcFromJPsi.push_back(isProbeAMcEleFromJPsi);
	  //
	  bool isTagAMcEle   = (Electron_genPartIdx[ele1_idx]>-0.5);
	  bool isProbeAMcEle = (Electron_genPartIdx[ele2_idx]>-0.5);
	  tag_matchMc.push_back(isTagAMcEle);
	  probe_matchMc.push_back(isProbeAMcEle);

	} else {
	  tag_matchMcFromJPsi.push_back(0);  
	  probe_matchMcFromJPsi.push_back(0);  
	  tag_matchMc.push_back(0);  
	  probe_matchMc.push_back(0);  
	}
      }
      //
      //
      if ( 1 ) { 
	//
	tag_pt.push_back(ele2_pt);
	tag_eta.push_back(ele2_eta);
	tag_phi.push_back(ele2_phi);
	tag_isPF.push_back(Electron_isPF[ele2_idx]);
	tag_isPFOverlap.push_back(Electron_isPFoverlap[ele2_idx]); 
	tag_isLowPt.push_back(Electron_isLowPt[ele2_idx]);
	tag_mvaId.push_back(Electron_mvaId[ele2_idx]);
	tag_pfmvaId.push_back(Electron_pfmvaId[ele2_idx]);
	tag_convveto.push_back(Electron_convVeto[ele2_idx]);
	tag_pfRelIso.push_back(Electron_pfRelIso[ele2_idx]);  
	//
	probe_pt.push_back(ele1_pt);    
	probe_eta.push_back(ele1_eta); 
	probe_phi.push_back(ele1_phi);    
	probe_isPF.push_back(Electron_isPF[ele1_idx]);
	probe_isPFOverlap.push_back(Electron_isPFoverlap[ele1_idx]); 
	probe_isLowPt.push_back(Electron_isLowPt[ele1_idx]);
	probe_mvaId.push_back(Electron_mvaId[ele1_idx]);
	probe_pfmvaId.push_back(Electron_pfmvaId[ele1_idx]);
	probe_dxySig.push_back(Electron_dxy[ele1_idx]/Electron_dxyErr[ele1_idx]);  
	probe_dzSig.push_back(Electron_dz[ele1_idx]/Electron_dzErr[ele1_idx]);  
	probe_pfRelIso.push_back(Electron_pfRelIso[ele1_idx]);  
	probe_trkRelIso.push_back(Electron_trkRelIso[ele1_idx]);  
	probe_fBrem.push_back(Electron_fBrem[ele1_idx]);  
	probe_unBiased.push_back(Electron_unBiased[ele1_idx]);  
	probe_ptBiased.push_back(Electron_ptBiased[ele1_idx]);  
	probe_convveto.push_back(Electron_convVeto[ele1_idx]);
	//
	probe_Bmass.push_back(thisBmass);
	probe_Bpt.push_back(thisBpt);
	probe_Bcos2D.push_back(thisBcos);
	probe_Bsvprob.push_back(thisBsvprob);
	probe_Bxysig.push_back(thisBxysig);
	probe_BmatchMC.push_back(isThisAMcB);
	probe_Kpt.push_back(thisKpt);
	//
	probe_invMass.push_back(mee);  
	//
	if (sampleID>0) {     // MC      
	  bool isTagAMcEleFromJPsi   = isMcEleFromJPsi(ele2_idx);
	  bool isProbeAMcEleFromJPsi = isMcEleFromJPsi(ele1_idx);
	  tag_matchMcFromJPsi.push_back(isTagAMcEleFromJPsi);
	  probe_matchMcFromJPsi.push_back(isProbeAMcEleFromJPsi);
	  //
	  bool isTagAMcEle   = (Electron_genPartIdx[ele2_idx]>-0.5);
	  bool isProbeAMcEle = (Electron_genPartIdx[ele1_idx]>-0.5);
	  tag_matchMc.push_back(isTagAMcEle);
	  probe_matchMc.push_back(isProbeAMcEle);
	} else {
	  tag_matchMcFromJPsi.push_back(0);  
	  probe_matchMcFromJPsi.push_back(0);  
	  tag_matchMc.push_back(0);  
	  probe_matchMc.push_back(0);  
	}
      }

    } // Loop over good Bs
      
    // At least one tag and one probe
    selectedPairsSize = tag_pt.size();
    if (selectedPairsSize<=0) continue;
    h_selection->Fill(5.,perEveW);

    // Filling the output tree
    outTree_->Fill();


    // Cleaning all vectors used for the selection
    cleanGoodBs.clear();
    
    // Cleaning all vectors used for the output tree, ready for a new entry
    tag_pt.clear();  
    tag_eta.clear();  
    tag_phi.clear();  
    tag_isPF.clear();  
    tag_isPFOverlap.clear();
    tag_isLowPt.clear();  
    tag_mvaId.clear();  
    tag_pfmvaId.clear();  
    tag_convveto.clear();
    tag_pfRelIso.clear();
    tag_matchMcFromJPsi.clear();  
    tag_matchMc.clear();  
    //
    probe_pt.clear();  
    probe_eta.clear();  
    probe_phi.clear();  
    probe_isPF.clear();  
    probe_isPFOverlap.clear();
    probe_isLowPt.clear();  
    probe_mvaId.clear();  
    probe_pfmvaId.clear();  
    probe_dxySig.clear();  
    probe_dzSig.clear();  
    probe_pfRelIso.clear();  
    probe_trkRelIso.clear();  
    probe_fBrem.clear();  
    probe_unBiased.clear();  
    probe_ptBiased.clear();  
    probe_convveto.clear();
    probe_matchMcFromJPsi.clear();  
    probe_matchMc.clear();  
    //
    probe_Bmass.clear();
    probe_Bpt.clear();
    probe_Bcos2D.clear();
    probe_Bsvprob.clear();    
    probe_Bxysig.clear();    
    probe_BmatchMC.clear();
    probe_Kpt.clear();
    //
    probe_invMass.clear();  
    //
  }
}

void TaPJpsiSelectionNaod::SetNvtxWeights(std::string nvtxWeightFile) {

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

float TaPJpsiSelectionNaod::GetNvtxWeight(float nvtx) {

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

bool TaPJpsiSelectionNaod::isMcB( int theB ) {
  
  // taking index
  int ele1_idx = BToKEE_l1Idx[theB];
  int ele2_idx = BToKEE_l2Idx[theB];
  int k_idx    = BToKEE_kIdx[theB];

  // Gen tree
  int k_genPartIdx      = ProbeTracks_genPartIdx[k_idx];  
  int k_genMotherIdx    = GenPart_genPartIdxMother[k_genPartIdx];
  int k_genGMotherIdx   = GenPart_genPartIdxMother[k_genMotherIdx];
  int k_genPdgId        = GenPart_pdgId[k_genPartIdx];
  int k_genMotherPdgId  = GenPart_pdgId[k_genMotherIdx];
  int k_genGMotherPdgId = GenPart_pdgId[k_genGMotherIdx];

  int ele1_genPartIdx      = Electron_genPartIdx[ele1_idx];  
  int ele1_genMotherIdx    = GenPart_genPartIdxMother[ele1_genPartIdx];
  int ele1_genGMotherIdx   = GenPart_genPartIdxMother[ele1_genMotherIdx];
  int ele1_genPdgId        = GenPart_pdgId[ele1_genPartIdx];
  int ele1_genMotherPdgId  = GenPart_pdgId[ele1_genMotherIdx];
  int ele1_genGMotherPdgId = GenPart_pdgId[ele1_genGMotherIdx];

  int ele2_genPartIdx      = Electron_genPartIdx[ele2_idx];  
  int ele2_genMotherIdx    = GenPart_genPartIdxMother[ele2_genPartIdx];
  int ele2_genGMotherIdx   = GenPart_genPartIdxMother[ele2_genMotherIdx];
  int ele2_genPdgId        = GenPart_pdgId[ele2_genPartIdx];
  int ele2_genMotherPdgId  = GenPart_pdgId[ele2_genMotherIdx];
  int ele2_genGMotherPdgId = GenPart_pdgId[ele2_genGMotherIdx];

  // B -> K J/psi(ll) at gen level
  bool okMatch = (ele1_genPartIdx>-0.5 && ele2_genPartIdx>-0.5 && k_genPartIdx>-0.5);
  bool RK_res1 = abs(ele1_genMotherPdgId)==443 && abs(k_genMotherPdgId)==521;
  bool RK_res2 = (ele1_genMotherPdgId==ele2_genMotherPdgId) && (k_genMotherPdgId==ele1_genGMotherPdgId) && (k_genMotherPdgId==ele2_genGMotherPdgId);
  bool RK_res = okMatch && RK_res1 && RK_res2;

  return RK_res;
}

bool TaPJpsiSelectionNaod::isMcEleFromJPsi( int ele_idx ) {

  // Gen tree
  int ele_genPartIdx      = Electron_genPartIdx[ele_idx];  
  int ele_genMotherIdx    = GenPart_genPartIdxMother[ele_genPartIdx];
  int ele_genGMotherIdx   = GenPart_genPartIdxMother[ele_genMotherIdx];
  int ele_genPdgId        = GenPart_pdgId[ele_genPartIdx];
  int ele_genMotherPdgId  = GenPart_pdgId[ele_genMotherIdx];
  int ele_genGMotherPdgId = GenPart_pdgId[ele_genGMotherIdx];

  // B -> K J/psi(ll) at gen level
  bool okMatch = (ele_genPartIdx>-0.5) && (abs(ele_genMotherPdgId)==443);

  return okMatch;
}

void TaPJpsiSelectionNaod::PrepareOutputs(std::string filename) 
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


void TaPJpsiSelectionNaod::bookOutputTree() 
{
  outTree_ = new TTree("TaPtree", "TaPtree");
  
  cout << "Booking tree" << endl;
  
  outTree_->Branch("theRun", &theRun, "theRun/I");    
  outTree_->Branch("theEvent", &theEvent, "theEvent/I");    
  outTree_->Branch("nvtx", &nvtx, "nvtx/I");    
  outTree_->Branch("sampleID", &sampleID, "sampleID/I");    
  outTree_->Branch("rho", &rho, "rho/F");    
  outTree_->Branch("pu_weight", &pu_weight, "pu_weight/F");    
  outTree_->Branch("pu_n", &pu_n, "pu_n/F");    
  outTree_->Branch("perEveW", &perEveW, "perEveW/F");    

  outTree_->Branch("hlt12ip6", &hlt12ip6, "hlt12ip6/I");
  outTree_->Branch("hlt9ip6", &hlt9ip6, "hlt9ip6/I");
  outTree_->Branch("hlt9ip5", &hlt9ip5, "hlt9ip5/I");
  outTree_->Branch("hlt9ip4", &hlt9ip4, "hlt9ip4/I");
  outTree_->Branch("hlt7ip4", &hlt7ip4, "hlt7ip4/I");
  outTree_->Branch("hlt8ip6", &hlt8ip6, "hlt8ip6/I");
  outTree_->Branch("hlt8ip5", &hlt8ip5, "hlt8ip5/I");
  outTree_->Branch("hlt8ip3", &hlt8ip3, "hlt8ip3/I");
  outTree_->Branch("hlt8d5ip3d5", &hlt8d5ip3d5, "hlt8d5ip3d5/I");
  outTree_->Branch("hlt10d5ip3d5", &hlt10d5ip3d5, "hlt10d5ip3d5/I");

  outTree_->Branch("selectedBSize",  &selectedBSize,  "selectedBSize/I");   
  
  outTree_->Branch("tag_pt", "std::vector<float>", &tag_pt);  
  outTree_->Branch("tag_eta", "std::vector<float>", &tag_eta);  
  outTree_->Branch("tag_phi", "std::vector<float>", &tag_phi);  
  outTree_->Branch("tag_isPF", "std::vector<bool>", &tag_isPF);  
  outTree_->Branch("tag_isPFOverlap", "std::vector<bool>", &tag_isPFOverlap);  
  outTree_->Branch("tag_isLowPt", "std::vector<bool>", &tag_isLowPt);  
  outTree_->Branch("tag_mvaId", "std::vector<float>", &tag_mvaId);  
  outTree_->Branch("tag_pfmvaId", "std::vector<float>", &tag_pfmvaId);  
  outTree_->Branch("tag_convveto", "std::vector<bool>", &tag_convveto);  
  outTree_->Branch("tag_pfRelIso", "std::vector<float>", &tag_pfRelIso);  
  outTree_->Branch("tag_matchMcFromJPsi", "std::vector<bool>", &tag_matchMcFromJPsi);  
  outTree_->Branch("tag_matchMc", "std::vector<bool>", &tag_matchMc);  

  outTree_->Branch("probe_pt", "std::vector<float>", &probe_pt);  
  outTree_->Branch("probe_eta", "std::vector<float>", &probe_eta);  
  outTree_->Branch("probe_phi", "std::vector<float>", &probe_phi);  
  outTree_->Branch("probe_isPF", "std::vector<bool>", &probe_isPF);  
  outTree_->Branch("probe_isPFOverlap", "std::vector<bool>", &probe_isPFOverlap);  
  outTree_->Branch("probe_isLowPt", "std::vector<bool>", &probe_isLowPt);  
  outTree_->Branch("probe_mvaId", "std::vector<float>", &probe_mvaId);  
  outTree_->Branch("probe_pfmvaId", "std::vector<float>", &probe_pfmvaId);  
  outTree_->Branch("probe_dxySig", "std::vector<float>", &probe_dxySig);  
  outTree_->Branch("probe_dzSig", "std::vector<float>", &probe_dzSig);  
  outTree_->Branch("probe_pfRelIso", "std::vector<float>", &probe_pfRelIso);  
  outTree_->Branch("probe_trkRelIso", "std::vector<float>", &probe_trkRelIso);  
  outTree_->Branch("probe_fBrem", "std::vector<float>", &probe_fBrem);  
  outTree_->Branch("probe_unBiased", "std::vector<float>", &probe_unBiased);  
  outTree_->Branch("probe_ptBiased", "std::vector<float>", &probe_ptBiased);  
  outTree_->Branch("probe_convveto", "std::vector<bool>", &probe_convveto);  
  outTree_->Branch("probe_invMass", "std::vector<float>", &probe_invMass);  
  outTree_->Branch("probe_matchMcFromJPsi", "std::vector<bool>", &probe_matchMcFromJPsi);  
  outTree_->Branch("probe_matchMc", "std::vector<bool>", &probe_matchMc);  

  outTree_->Branch("probe_Bmass",   "std::vector<float>", &probe_Bmass);  
  outTree_->Branch("probe_Bpt",     "std::vector<float>", &probe_Bpt);  
  outTree_->Branch("probe_Bcos2D",  "std::vector<float>", &probe_Bcos2D);  
  outTree_->Branch("probe_Bsvprob", "std::vector<float>", &probe_Bsvprob);  
  outTree_->Branch("probe_Bxysig",  "std::vector<float>", &probe_Bxysig);  
  outTree_->Branch("probe_BmatchMC", "std::vector<bool>", &probe_BmatchMC);  
  outTree_->Branch("probe_Kpt", "std::vector<float>", &probe_Kpt);  

  outTree_->Branch("selectedPairsSize",  &selectedPairsSize,  "selectedPairsSize/I");   
}

void TaPJpsiSelectionNaod::bookOutputHistos() 
{
  cout << "Booking histos" << endl;
  //
  h_entries   = new TH1F("h_entries",  "Number of entries",   3,  3.5, 6.5);
  h_selection = new TH1F("h_selection","Selection breakdown", 8, -0.5, 7.5);
}
