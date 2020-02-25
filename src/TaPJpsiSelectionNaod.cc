#include "TMath.h"
#include "TTree.h"
#include "TLorentzVector.h"

#include <iostream> 

#include "../include/TaPJpsiSelectionNaod.hh"    

#define MAX_PU_REWEIGHT 60

using namespace std;

TaPJpsiSelectionNaod::TaPJpsiSelectionNaod(TTree *tree)     
  : BParkBase(tree) {        

  // Chiara: To be set by hand   
  dopureweight_ = false;
  sampleID = 0;    // 0 = data
  puWFileName_ = "puWFileName";   
  lumiWeight_  = 1.;
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
    if (jentry%10000==0) cout << jentry << endl;

    // To keep track of the total number of events 
    h_entries->Fill(5);
    
    // Event info     
    theRun   = run;
    theLumi  = luminosityBlock;
    theEvent = event;
    theSampleID = sampleID;
    theLumiWeight = lumiWeight_;

    // # Vertices
    nvtx = nOtherPV+1;   // chiara, controlla che sia ok

    // Energy density
    rho = fixedGridRhoFastjetAll;
    
    // PU weight (for MC only and if requested)
    pu_weight = 1.;
    pu_n      = -1.;
    /*
    if (sampleID>0) {     // MC
      pu_n = Pileup_nTrueInt;           // chiara: verifica che sia questo
      if (dopureweight_) pu_weight = GetPUWeight(pu_n);         
    }
    */

    // other weights for the dataset
    perEveW = 1.;
    //if (sampleID>0) { 
    //const auto & eveWeights = genInfo->weights();
    //if(!eveWeights.empty()) perEveW = eveWeights[0];
    //}


    // Events breakdown  
    h_selection->Fill(0.,perEveW);
    

    // ----------------------------------------------------
    // save events only if:
    // 1) good vertex
    // 2) triggering muon
    // 3) reconstructed B + JPsi
    // 2) at least one tag and one probe
    // ----------------------------------------------------
    

    // PV must be there - chiara: should always be the case
    bool goodPV = false;
    if (PV_x>-999) goodPV = true;
    if (!goodPV) continue;
    h_selection->Fill(1.,perEveW);

    // Trigger - chiara: capire se basta cosi' o se bisogna controllare altro (L1?)
    bool okTrigger = false;
    if (HLT_Mu7_IP4 || HLT_Mu8_IP6 || HLT_Mu8_IP5 || HLT_Mu8_IP3 || HLT_Mu8p5_IP3p5 || HLT_Mu9_IP6 || HLT_Mu9_IP5 || HLT_Mu9_IP4 || HLT_Mu10p5_IP3p5 || HLT_Mu12_IP6) okTrigger = true;
    if (!okTrigger) continue;

    // Triggering muons - chiara
    if (nTriggerMuon<=0) continue;
    h_selection->Fill(2.,perEveW);

    // B candidates
    if (nBToKEE<=0) continue;
    h_selection->Fill(3.,perEveW);

    // Minimal Bcandidate requirements
    vector<int> goodBs;
    for (u_int iB=0; iB<nBToKEE; iB++) {

      // preparing variables
      int ele1_idx = BToKEE_l1Idx[iB];
      int ele2_idx = BToKEE_l2Idx[iB];
      int k_idx    = BToKEE_kIdx[iB];

      float ele1_pt = Electron_pt[ele1_idx];
      float ele2_pt = Electron_pt[ele2_idx];
      float k_pt    = ProbeTracks_pt[k_idx];

      float ele1_eta = Electron_eta[ele1_idx];
      float ele2_eta = Electron_eta[ele2_idx];
      float k_eta    = ProbeTracks_eta[k_idx];
      float ele1_phi = Electron_phi[ele1_idx];
      float ele2_phi = Electron_phi[ele2_idx];

      bool ele1_convveto = Electron_convVeto[ele1_idx];     
      bool ele2_convveto = Electron_convVeto[ele2_idx];     

      float b_xySig = BToKEE_l_xy[iB]/BToKEE_l_xy_unc[iB];

      // B selection
      // standard cut based selection:  BToKEE_fit_pt[iB]>10.0, cos>0.999
      bool vtxFitSel = BToKEE_fit_pt[iB]>3.0 && b_xySig>6.0 && BToKEE_svprob[iB]>0.1 && BToKEE_fit_cos2D[iB]>0.99;
      bool ele1Sel = ele1_convveto && BToKEE_fit_l1_pt[iB]>0.5 && abs(ele1_eta)<2.4;  
      bool ele2Sel = ele2_convveto && BToKEE_fit_l2_pt[iB]>0.5 && abs(ele2_eta)<2.4;  
      // standard cut based selection: pt>1.5
      bool kSel = k_pt>0.8 && fabs(k_eta)<2.4; 
      bool additionalSel = BToKEE_fit_mass[iB]>4.5 && BToKEE_fit_mass[iB]<6.0;
      bool isBsel = vtxFitSel && ele1Sel && ele2Sel && kSel && additionalSel;

      // JPsi selection - capire se quantita pre-post fit, chiara
      TLorentzVector ele1TLV(0,0,0,0);
      ele1TLV.SetPtEtaPhiM(ele1_pt,ele1_eta,ele1_phi,0);
      TLorentzVector ele2TLV(0,0,0,0);
      ele2TLV.SetPtEtaPhiM(ele2_pt,ele2_eta,ele2_phi,0);
      float mee = (ele1TLV+ele2TLV).M();
      if (mee>3.4 || mee<2.8) continue;  // chiara, mio

      if (!isBsel) continue;

      goodBs.push_back(iB);
    }
    if (goodBs.size()>0) h_selection->Fill(4.,perEveW);

    
    // ------------------------------------------------------------------
    // Selecting just one B among those with ee in common
    vector<int> cleanGoodBs;
    if (goodBs.size()>0) cleanGoodBs.push_back(goodBs[0]);

    for (u_int iB=1; iB<goodBs.size(); iB++) {
      int thisB    = goodBs[iB];
      int ele1_idx = BToKEE_l1Idx[thisB];
      int ele2_idx = BToKEE_l2Idx[thisB];
      float BFitprob = BToKEE_svprob[thisB];

      int toadd = -1;
      vector<int> toBeErased;
      for (u_int iBC=0; iBC<cleanGoodBs.size(); iBC++) {
	int thisBC    = cleanGoodBs[iBC];
	int ele1_idxC = BToKEE_l1Idx[thisBC];
	int ele2_idxC = BToKEE_l2Idx[thisBC];
	float BCFitprob = BToKEE_svprob[thisBC];
	
	// already there - eventually to be replaced
	if ( ele1_idxC==ele1_idx && ele2_idxC==ele2_idx && BFitprob>BCFitprob ){
	  toBeErased.push_back(iBC);
	  cleanGoodBs.push_back(thisB);
	  // this B is already in, not to be added
	  toadd = -100;      

	} else if ( ele1_idxC==ele1_idx && ele2_idxC==ele2_idx && BFitprob<=BCFitprob ) {  // there's already a better B candidate with same leptons -> not to be added
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
    h_selection->Fill(5.,perEveW);


    // ------------------------------------------------------------------
    // Tag and probe using all the J/Psi candidates and the 2 possible combinations
    for (u_int iB=0; iB<cleanGoodBs.size(); iB++) {

      // B-infos for further cuts offline
      int thisB = cleanGoodBs[iB];
      float thisBmass = BToKEE_fit_mass[thisB];
      float thisBpt   = BToKEE_fit_pt[thisB];
      float thisBcos  = BToKEE_fit_cos2D[thisB];
      bool isThisAMcB = -1;
      if (sampleID>0) isThisAMcB = isMcB(thisB);

      // K-infos for further cuts offline
      int thisK = BToKEE_kIdx[thisB];
      float thisKpt = ProbeTracks_pt[thisK];

      // Electrons
      int ele1_idx = BToKEE_l1Idx[thisB];
      int ele2_idx = BToKEE_l2Idx[thisB];

      float ele1_pt  = Electron_pt[ele1_idx];
      float ele2_pt  = Electron_pt[ele2_idx];
      float ele1_eta = Electron_eta[ele1_idx];
      float ele2_eta = Electron_eta[ele2_idx];
      float ele1_phi = Electron_phi[ele1_idx];
      float ele2_phi = Electron_phi[ele2_idx];

      TLorentzVector ele1TLV(0,0,0,0);
      ele1TLV.SetPtEtaPhiM(ele1_pt,ele1_eta,ele1_phi,0);
      TLorentzVector ele2TLV(0,0,0,0);
      ele2TLV.SetPtEtaPhiM(ele2_pt,ele2_eta,ele2_phi,0);
      float mee = (ele1TLV+ele2TLV).M();

      bool is1tag   = isTag(ele1_idx);
      bool is2tag   = isTag(ele2_idx);
      bool is1probe = isProbe(ele1_idx);
      bool is2probe = isProbe(ele2_idx);

      if ( is1tag && is2probe ) { 
	//
	tag_pt.push_back(Electron_pt[ele1_idx]);
	tag_eta.push_back(Electron_eta[ele1_idx]);
	tag_phi.push_back(Electron_phi[ele1_idx]);
	tag_isPF.push_back(Electron_isPF[ele1_idx]);
	tag_isLowPt.push_back(Electron_isLowPt[ele1_idx]);
	tag_mvaId.push_back(Electron_mvaId[ele1_idx]);
	tag_unBiased.push_back(Electron_unBiased[ele1_idx]);  
	tag_pfmvaId.push_back(Electron_pfmvaId[ele1_idx]);
	//
	probe_Bmass.push_back(thisBmass);
	probe_Bpt.push_back(thisBpt);
	probe_Bcos2D.push_back(thisBcos);
	probe_BmatchMC.push_back(isThisAMcB);
	probe_pt.push_back(Electron_pt[ele2_idx]);
	probe_eta.push_back(Electron_eta[ele2_idx]);
	probe_phi.push_back(Electron_phi[ele2_idx]);
	probe_isPF.push_back(Electron_isPF[ele2_idx]);
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
	probe_isTag.push_back( is2tag );
	probe_invMass.push_back(mee);  
	//
	if (sampleID>0) {     // MC      
	  bool isTagAMcEle   = isMcEleFromJPsi(ele1_idx);
	  bool isProbeAMcEle = isMcEleFromJPsi(ele2_idx);
	  tag_matchMC.push_back(isTagAMcEle);
	  probe_matchMC.push_back(isProbeAMcEle);
	} else {
	  tag_matchMC.push_back(0);  
	  probe_matchMC.push_back(0);  
	}
      }
      //
      //
      if ( is2tag && is1probe ) { 
	//
	tag_pt.push_back(Electron_pt[ele2_idx]);
	tag_eta.push_back(Electron_eta[ele2_idx]);
	tag_phi.push_back(Electron_phi[ele2_idx]);
	tag_isPF.push_back(Electron_isPF[ele2_idx]);
	tag_isLowPt.push_back(Electron_isLowPt[ele2_idx]);
	tag_mvaId.push_back(Electron_mvaId[ele2_idx]);
	tag_unBiased.push_back(Electron_unBiased[ele2_idx]);  
	tag_pfmvaId.push_back(Electron_pfmvaId[ele2_idx]);
	//
	probe_Bmass.push_back(thisBmass);
	probe_Bpt.push_back(thisBpt);
	probe_Bcos2D.push_back(thisBcos);
	probe_BmatchMC.push_back(isThisAMcB);
	probe_pt.push_back(Electron_pt[ele1_idx]);
	probe_eta.push_back(Electron_eta[ele1_idx]);
	probe_phi.push_back(Electron_phi[ele1_idx]);
	probe_isPF.push_back(Electron_isPF[ele1_idx]);
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
	probe_isTag.push_back( is1tag );
	probe_invMass.push_back(mee);  
	//
	if (sampleID>0) {     // MC      
	  bool isTagAMcEle   = isMcEleFromJPsi(ele2_idx);
	  bool isProbeAMcEle = isMcEleFromJPsi(ele1_idx);
	  tag_matchMC.push_back(isTagAMcEle);
	  probe_matchMC.push_back(isProbeAMcEle);
	} else {
	  tag_matchMC.push_back(0);  
	  probe_matchMC.push_back(0);  
	}
      }

    } // Loop over good Bs
      
    // At least one tag and one probe
    selectedPairsSize = tag_pt.size();
    if (selectedPairsSize<=0) continue;
    h_selection->Fill(6.,perEveW);

    // Filling the output tree
    outTree_->Fill();


    // Cleaning all vectors used for the selection
    cleanGoodBs.clear();
    
    // Cleaning all vectors used for the output tree, ready for a new entry
    tag_pt.clear();  
    tag_eta.clear();  
    tag_phi.clear();  
    tag_isPF.clear();  
    tag_isLowPt.clear();  
    tag_mvaId.clear();  
    tag_pfmvaId.clear();  
    tag_unBiased.clear();  
    tag_matchMC.clear();  
    //
    probe_Bmass.clear();
    probe_Bpt.clear();
    probe_Bcos2D.clear();
    probe_BmatchMC.clear();
    probe_pt.clear();  
    probe_eta.clear();  
    probe_phi.clear();  
    probe_isPF.clear();  
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
    probe_isTag.clear();  
    probe_invMass.clear();  
    probe_matchMC.clear();  
  }
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

bool TaPJpsiSelectionNaod::isTag( int theEle ) {

  bool isTag = true;
  
  // chiara

  return isTag;  
}

bool TaPJpsiSelectionNaod::isProbe( int theEle ) {

  bool isProbe = true;
  
  // chiara

  return isProbe;  
}

void TaPJpsiSelectionNaod::PrepareOutputs(std::string filename) 
{
  _datasetName=filename;

  std::string outname = _datasetName+".root";    
  cout << "output: " << outname << endl;
  outFile_ = new TFile(outname.c_str(),"RECREATE");

  bookOutputTree();
  bookOutputHistos();
};


void TaPJpsiSelectionNaod::bookOutputTree() 
{
  outTree_ = new TTree("TaPtree", "TaPtree");
  
  cout << "Booking tree" << endl;
  
  outTree_->Branch("theRun", &theRun, "theRun/I");    
  outTree_->Branch("theEvent", &theEvent, "theEvent/I");    
  outTree_->Branch("theLumi", &theLumi, "theLumi/I");    
  outTree_->Branch("nvtx", &nvtx, "nvtx/I");    
  outTree_->Branch("sampleID", &sampleID, "sampleID/I");    
  outTree_->Branch("rho", &rho, "rho/F");    
  outTree_->Branch("theLumiWeight", &theLumiWeight, "theLumiWeight/F");    
  outTree_->Branch("pu_weight", &pu_weight, "pu_weight/F");    
  outTree_->Branch("pu_n", &pu_n, "pu_n/F");    
  outTree_->Branch("perEveW", &perEveW, "perEveW/F");    

  outTree_->Branch("selectedBSize",  &selectedBSize,  "selectedBSize/I");   
  
  outTree_->Branch("tag_pt", "std::vector<float>", &tag_pt);  
  outTree_->Branch("tag_eta", "std::vector<float>", &tag_eta);  
  outTree_->Branch("tag_phi", "std::vector<float>", &tag_phi);  
  outTree_->Branch("tag_isPF", "std::vector<bool>", &tag_isPF);  
  outTree_->Branch("tag_isLowPt", "std::vector<bool>", &tag_isLowPt);  
  outTree_->Branch("tag_mvaId", "std::vector<float>", &tag_mvaId);  
  outTree_->Branch("tag_pfmvaId", "std::vector<float>", &tag_pfmvaId);  
  outTree_->Branch("tag_unBiased", "std::vector<float>", &tag_unBiased);  
  outTree_->Branch("tag_matchMC", "std::vector<bool>", &tag_matchMC);  

  outTree_->Branch("probe_Bmass", "std::vector<float>", &probe_Bmass);  
  outTree_->Branch("probe_Bpt", "std::vector<float>", &probe_Bpt);  
  outTree_->Branch("probe_Bcos2D", "std::vector<float>", &probe_Bcos2D);  
  outTree_->Branch("probe_BmatchMC", "std::vector<bool>", &probe_BmatchMC);  
  outTree_->Branch("probe_pt", "std::vector<float>", &probe_pt);  
  outTree_->Branch("probe_eta", "std::vector<float>", &probe_eta);  
  outTree_->Branch("probe_phi", "std::vector<float>", &probe_phi);  
  outTree_->Branch("probe_isPF", "std::vector<bool>", &probe_isPF);  
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
  outTree_->Branch("probe_isTag", "std::vector<bool>", &probe_isTag);  
  outTree_->Branch("probe_invMass", "std::vector<float>", &probe_invMass);  
  outTree_->Branch("probe_matchMC", "std::vector<bool>", &probe_matchMC);  

  outTree_->Branch("selectedPairsSize",  &selectedPairsSize,  "selectedPairsSize/I");   
}

void TaPJpsiSelectionNaod::bookOutputHistos() 
{
  cout << "Booking histos" << endl;
  //
  h_entries   = new TH1F("h_entries",  "Number of entries",   3,  3.5, 6.5);
  h_selection = new TH1F("h_selection","Selection breakdown", 8, -0.5, 7.5);
}
