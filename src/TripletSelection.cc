#include "TMath.h"
#include "TTree.h"
#include "TLorentzVector.h"

#include <iostream> 

#include "../include/TripletSelection.hh"    

#define MAX_PU_REWEIGHT 60

using namespace std;

TripletSelection::TripletSelection(TTree *tree)     
  : BParkBase(tree) { }

TripletSelection::~TripletSelection() {

  // output
  outFile_->cd();
  h_entries -> Write();
  h_selection -> Write();
  outTree_->Write();
  outFile_->Close();
}     

void TripletSelection::Loop() {

  if (fChain == 0) return;

  // Loop over events
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  int theE = fChain->GetEntriesFast();
  cout << "entries : " <<  theE << endl;

  // Extra checks -----------
  int myAllGoodBs=0;
  int myAllGoodPreselBs=0;
  int allTriplets=0;
  int nTriplets[20];
  for (int nn=0; nn<20; nn++) nTriplets[nn]=0;
  // Extra checks -----------


  for (Long64_t jentry=0; jentry<nentries;jentry++) {

    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (jentry%10000==0) cout << jentry << endl;

    float perEveW = 1.;

    // To keep track of the total number of events 
    h_entries->Fill(5);

    // Event info 
    theEvent = event;
    
    // Energy density
    rho = fixedGridRhoFastjetAll;
    
    // Events breakdown  
    h_selection->Fill(0.,perEveW);
    
    // PV must be there
    bool goodPV = false;
    if (PV_x>-999) goodPV = true;
    if (!goodPV) continue;
    h_selection->Fill(1.,perEveW);

    // Trigger
    iHLT_Mu12_IP6 = (int)HLT_Mu12_IP6;
    iHLT_Mu9_IP6  = (int)HLT_Mu9_IP6;

    bool okTrigger = false;
    if (HLT_Mu7_IP4 || HLT_Mu8_IP6 || HLT_Mu8_IP5 || HLT_Mu8_IP3 || HLT_Mu8p5_IP3p5 || HLT_Mu9_IP6 || HLT_Mu9_IP5 || HLT_Mu9_IP4 || HLT_Mu10p5_IP3p5 || HLT_Mu12_IP6) okTrigger = true;
    if (!okTrigger) continue;

    // Triggering muons
    if (nTriggerMuon<=0) continue;
    h_selection->Fill(2.,perEveW);

    // B candidates
    if (nBToKEE<=0) continue;
    h_selection->Fill(3.,perEveW);


    // ------------------ chiara -----------------------------------------
    // Check: efficiency - per true B - of my further preselection
    for (u_int iB=0; iB<nBToKEE; iB++) {

      // Only good Bs
      if (isMcB(iB)==0) continue;
      
      // preparing variables
      int ele1_idx = BToKEE_l1Idx[iB];
      int ele2_idx = BToKEE_l2Idx[iB];
      int k_idx    = BToKEE_kIdx[iB];

      // Removing B from low-pT electrons with PF overlap
      if (Electron_isPFoverlap[ele1_idx]==1) continue;
      if (Electron_isPFoverlap[ele2_idx]==1) continue;

      myAllGoodBs++;

      // Further preselection
      float k_pt     = ProbeTracks_pt[k_idx];
      float ele1_eta = Electron_eta[ele1_idx];
      float ele2_eta = Electron_eta[ele2_idx];
      float k_eta    = ProbeTracks_eta[k_idx];
      bool ele1_convveto = Electron_convVeto[ele1_idx];     
      bool ele2_convveto = Electron_convVeto[ele2_idx];     
      float b_xySig = BToKEE_l_xy[iB]/BToKEE_l_xy_unc[iB];

      bool vtxFitSel = BToKEE_fit_pt[iB]>3.0 && b_xySig>6.0 && BToKEE_svprob[iB]>0.1 && BToKEE_fit_cos2D[iB]>0.99;
      bool ele1Sel = ele1_convveto && BToKEE_fit_l1_pt[iB]>0.5 && abs(ele1_eta)<2.4;  
      bool ele2Sel = ele2_convveto && BToKEE_fit_l2_pt[iB]>0.5 && abs(ele2_eta)<2.4;  
      bool kSel = k_pt>0.7 && fabs(k_eta)<2.4; 
      bool additionalSel = BToKEE_fit_mass[iB]>4.5 && BToKEE_fit_mass[iB]<6.0;
      bool isBsel = vtxFitSel && ele1Sel && ele2Sel && kSel && additionalSel;
      if (!isBsel) continue;

      //bool islplp=true;
      //if (Electron_pfmvaId[ele1_idx]<20) islplp=false;
      //if (Electron_pfmvaId[ele2_idx]<20) islplp=false;
      //if (!islplp) continue;

      myAllGoodPreselBs++;
    }
    
    // -----------------------------------------------
    // Minimal Bcandidate requirements
    vector<int> goodBs;
    vector<int> goodTrueBs;
    vector<int> goodCombBs;
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
      float k_phi    = ProbeTracks_phi[k_idx];  

      bool ele1_convveto = Electron_convVeto[ele1_idx];     
      bool ele2_convveto = Electron_convVeto[ele2_idx];     

      float b_xySig = BToKEE_l_xy[iB]/BToKEE_l_xy_unc[iB];

      // B selection - usata presentazione 1 aprile
      bool vtxFitSel = BToKEE_fit_pt[iB]>3.0 && b_xySig>6.0 && BToKEE_svprob[iB]>0.1 && BToKEE_fit_cos2D[iB]>0.99;
      bool ele1Sel = ele1_convveto && BToKEE_fit_l1_pt[iB]>0.5 && abs(ele1_eta)<2.4;  
      bool ele2Sel = ele2_convveto && BToKEE_fit_l2_pt[iB]>0.5 && abs(ele2_eta)<2.4;  
      bool kSel = k_pt>0.7 && fabs(k_eta)<2.4; 
      bool additionalSel = BToKEE_fit_mass[iB]>4.5 && BToKEE_fit_mass[iB]<6.0;
      bool isBsel = vtxFitSel && ele1Sel && ele2Sel && kSel && additionalSel;

      // Extra to be tested
      //if (ele1_pt<=5 && Electron_pfmvaId[ele1_idx]<-2.4) continue;   
      //if (ele2_pt<=5 && Electron_pfmvaId[ele2_idx]<-2.4) continue;   
      //if (ele1_pt>5 && Electron_pfmvaId[ele1_idx]<-3.2)  continue;   
      //if (ele2_pt>5 && Electron_pfmvaId[ele2_idx]<-3.2)  continue;   
      //if (Electron_mvaId[ele1_idx]<0)  continue;
      //if (Electron_mvaId[ele2_idx]<-3) continue;
      //if (k_pt<=BToKEE_fit_l2_pt[iB] && k_pt<0.8) continue;
      //if (k_pt>BToKEE_fit_l2_pt[iB] && BToKEE_fit_l2_pt[iB]<0.8) continue;
      // Extra to be tested

      //if (!isBsel) continue;

      // Da tenere sempre, a meno di restringermi solo a lowPt-lowPt
      // Removing B from low-pT electrons with PF overlap
      if (Electron_isPFoverlap[ele1_idx]==1) continue;
      if (Electron_isPFoverlap[ele2_idx]==1) continue;

      // Only lowPt-lowPt, without PF overlap removal
      // bool islplp=true;
      // if (Electron_pfmvaId[ele1_idx]<20) islplp=false;
      // if (Electron_pfmvaId[ele2_idx]<20) islplp=false;
      // if (!islplp) continue;

      // how many "good" Bs, real or from combinatorics
      goodBs.push_back(iB);
      if (isMcB(iB)==1) goodTrueBs.push_back(iB);
      if (isMcB(iB)==0) goodCombBs.push_back(iB);
    }


    // Counting, per event
    goodBSize       = goodBs.size();
    goodTrueBSize   = goodTrueBs.size();
    goodCombBSize   = goodCombBs.size();

    // At least one good B candidate
    if (goodBs.size()<=0) continue;
    h_selection->Fill(4.,perEveW);


    // -------------------------------------------------
    // Mc truth
    int myGenEle = -1;
    int myGenPos = -1;
    int myGenK   = -1;
    for (int iGen=0; iGen<nGenPart; iGen++) {
      int theId    = GenPart_pdgId[iGen];
      int theMum   = GenPart_genPartIdxMother[iGen];
      int theMumId = GenPart_pdgId[theMum];
      if (abs(theMumId)==521 && theId==11)       myGenPos = iGen;
      if (abs(theMumId)==521 && theId==-11)      myGenEle = iGen;
      if (abs(theMumId)==521 && abs(theId)==321) myGenK   = iGen;
    }
    TVector3 genEle(0.,0.,0.);
    TVector3 genPos(0.,0.,0.);
    TVector3 genK(0.,0.,0.);
    genEle.SetPtEtaPhi(GenPart_pt[myGenEle], GenPart_eta[myGenEle], GenPart_phi[myGenEle]);
    genPos.SetPtEtaPhi(GenPart_pt[myGenPos], GenPart_eta[myGenPos], GenPart_phi[myGenPos]);
    genK.SetPtEtaPhi(GenPart_pt[myGenK], GenPart_eta[myGenK], GenPart_phi[myGenK]);
    // to verify that all events have the correct MC chain
    rightMcTruth = 1;
    if (myGenEle<0 || myGenPos<0 || myGenK<0) rightMcTruth = 0;


    // -------------------------------------------------
    // Extra: Hand-made triplets
    handMadeB = -1;
    handMadeBmass = -999.;
    int myBestRecoEle = -1;
    int myBestRecoPos = -1;
    int myBestRecoK   = -1;
    float mindREle = 999.;
    float mindRPos = 999.;
    float mindRK   = 999.;
    for (int iEle=0; iEle<nElectron; iEle++) {
      TVector3 thisEle(0.,0.,0.);
      thisEle.SetPtEtaPhi(Electron_pt[iEle],Electron_eta[iEle],Electron_phi[iEle]);
      float dREle = thisEle.DeltaR(genEle);
      float dRPos = thisEle.DeltaR(genPos);
      float dPtOverPtEle = fabs(Electron_pt[iEle]-GenPart_pt[myGenEle])/Electron_pt[iEle];
      float dPtOverPtPos = fabs(Electron_pt[iEle]-GenPart_pt[myGenPos])/Electron_pt[iEle];
      if (dREle<mindREle && dPtOverPtEle<0.5) {
	mindREle = dREle;
	myBestRecoEle = iEle;
      }
      if (dRPos<mindRPos && dPtOverPtPos<0.5) { 
	mindRPos = dRPos;
	myBestRecoPos = iEle;
      }
    }
    for (int iK=0; iK<nProbeTracks; iK++) {
      TVector3 thisK(0.,0.,0.); 
      thisK.SetPtEtaPhi(ProbeTracks_pt[iK],ProbeTracks_eta[iK],ProbeTracks_phi[iK]);
      float dR = thisK.DeltaR(genK);
      float dPtOverPt = fabs(ProbeTracks_pt[iK]-GenPart_pt[myGenK])/ProbeTracks_pt[iK];
      if (dR<mindRK && dPtOverPt<0.5) {
	mindRK = dR;
	myBestRecoK = iK;
      }
    }
    if (mindREle<0.03 && mindRPos<0.03 && mindRK<0.03) {
      handMadeB = 1; 
      TLorentzVector bestEle1(0.,0.,0.,0.);
      TLorentzVector bestEle2(0.,0.,0.,0.);
      TLorentzVector bestK(0.,0.,0.,0.);
      bestEle1.SetPtEtaPhiM(Electron_pt[myBestRecoEle], Electron_eta[myBestRecoEle], Electron_phi[myBestRecoEle], 0.000511);
      bestEle2.SetPtEtaPhiM(Electron_pt[myBestRecoPos], Electron_eta[myBestRecoPos], Electron_phi[myBestRecoPos], 0.000511);
      bestK.SetPtEtaPhiM(ProbeTracks_pt[myBestRecoK],ProbeTracks_eta[myBestRecoK],ProbeTracks_phi[myBestRecoK], 0.000494);    
      TLorentzVector bestHandMadeB = bestEle1+bestEle2+bestK;
      handMadeBmass = bestHandMadeB.M();
    }


    // -----------------------------------------------------------
    // Extra: reduction when keeping just 1,2,3... N triplets / event
    for (u_int iB=0; iB<goodBs.size(); iB++) {        
      allTriplets++;
      for (int nn=0; nn<20; nn++) {
	int nnp1=nn+1;
	if (iB<nnp1) nTriplets[nn]++;
      }}


    // ------------------------------------
    // Debug SV prob: flat or not?
    if (goodBs.size()>=1) {
      int myBestSvprob = -1;
      float mySvprobMax  = -999.;
      for (u_int iB=0; iB<goodBs.size(); iB++) {        
	int thisB = goodBs[iB];
	if (BToKEE_svprob[thisB]>mySvprobMax) {
	  mySvprobMax  = BToKEE_svprob[thisB];
	  myBestSvprob = thisB; 
	}
      }
      int debug_ele1_idx = BToKEE_l1Idx[myBestSvprob];
      int debug_ele2_idx = BToKEE_l2Idx[myBestSvprob];
      debug_svprob.push_back(mySvprobMax);      
      if ( (Electron_isPF[debug_ele1_idx])==1 && (Electron_isPF[debug_ele2_idx])==1 ) debug_pf_svprob.push_back(mySvprobMax);
      if (isMcB(myBestSvprob)==1) {
	debug_svprob_match.push_back(1);  
	if ( (Electron_isPF[debug_ele1_idx])==1 && (Electron_isPF[debug_ele2_idx])==1 ) debug_pf_svprob_match.push_back(1);
      } else {
	debug_svprob_match.push_back(0);  
	if ( (Electron_isPF[debug_ele1_idx])==1 && (Electron_isPF[debug_ele2_idx])==1 ) debug_pf_svprob_match.push_back(0);
      }
    }
    // Debug SV prob: flat or not?
    // ------------------------------------


    // ------------------------------------------------------------
    // "True" Bs: duplication of electrons and tracks for cases with at least 2 "true" Bs
    if (goodTrueBs.size()>1) { 
      for (u_int iB0=0; iB0<(goodTrueBs.size()-1); iB0++) {        
	
	float maxMinDREle = -999.;
	float maxMinDREle_dEta = -999.;
	float maxMinDREle_dPhi = -999.;
	float maxMinDREle_dPtOverPt = -999.;
	float maxDRTrack = -999.;
	
	int thisB0      = goodTrueBs[iB0];
	int ele01_idx   = BToKEE_l1Idx[thisB0];
	int ele02_idx   = BToKEE_l2Idx[thisB0];
	int k0_idx      = BToKEE_kIdx[thisB0];
	float ele01_pt  = Electron_pt[ele01_idx];
	float ele02_pt  = Electron_pt[ele02_idx];
	float ele01_eta = Electron_eta[ele01_idx];
	float ele02_eta = Electron_eta[ele02_idx];
	float ele01_phi = Electron_phi[ele01_idx];
	float ele02_phi = Electron_phi[ele02_idx];
	float k0_pt     = ProbeTracks_pt[k0_idx];
	float k0_eta    = ProbeTracks_eta[k0_idx];
	float k0_phi    = ProbeTracks_phi[k0_idx];
	float BFitprob0 = BToKEE_svprob[thisB0];
	
	for (u_int iB1=iB0+1; iB1<goodTrueBs.size(); iB1++) {      
	  
	  int thisB1      = goodTrueBs[iB1];
	  int ele11_idx   = BToKEE_l1Idx[thisB1];
	  int ele12_idx   = BToKEE_l2Idx[thisB1];
	  int k1_idx      = BToKEE_kIdx[thisB1];
	  float ele11_pt  = Electron_pt[ele11_idx];
	  float ele12_pt  = Electron_pt[ele12_idx];
	  float ele11_eta = Electron_eta[ele11_idx];
	  float ele12_eta = Electron_eta[ele12_idx];
	  float ele11_phi = Electron_phi[ele11_idx];
	  float ele12_phi = Electron_phi[ele12_idx];
	  float k1_pt     = ProbeTracks_pt[k1_idx];
	  float k1_eta    = ProbeTracks_eta[k1_idx];
	  float k1_phi    = ProbeTracks_phi[k1_idx];
	  float BFitprob1 = BToKEE_svprob[thisB1];
	  
	  TVector3 ele01(0.,0.,0.);
	  TVector3 ele02(0.,0.,0.);
	  TVector3 ele11(0.,0.,0.);
	  TVector3 ele12(0.,0.,0.);
	  ele01.SetPtEtaPhi(ele01_pt, ele01_eta, ele01_phi);
	  ele02.SetPtEtaPhi(ele02_pt, ele02_eta, ele02_phi);
	  ele11.SetPtEtaPhi(ele11_pt, ele11_eta, ele11_phi);
	  ele12.SetPtEtaPhi(ele12_pt, ele12_eta, ele12_phi);
	  float deltaR_01_11 = ele01.DeltaR(ele11); 
	  float deltaR_01_12 = ele01.DeltaR(ele12); 
	  float deltaR_02_11 = ele02.DeltaR(ele11); 
	  float deltaR_02_12 = ele02.DeltaR(ele12); 
	  
	  TVector3 track0(0.,0.,0.);
	  TVector3 track1(0.,0.,0.);
	  track0.SetPtEtaPhi(k0_pt, k0_eta, k0_phi);
	  track1.SetPtEtaPhi(k1_pt, k1_eta, k1_phi);
	  float deltaR_tracks = track0.DeltaR(track1);
	  
	  // electrons
	  float min1 = -999;
	  float min2 = -999;
	  float theDEta1 = -999.;
	  float theDEta2 = -999.;
	  float theDPhi1 = -999.;
	  float theDPhi2 = -999.;
	  float theDPtOverPt1 = -999.;
	  float theDPtOverPt2 = -999.;
	  if ( deltaR_01_11<deltaR_01_12 ) { 
	    min1 = deltaR_01_11;
	    theDEta1 = fabs(ele01_eta-ele11_eta);
	    theDPhi1 = ele01.DeltaPhi(ele11);
	    theDPtOverPt1 = (ele01_pt-ele11_pt)/(ele01_pt+ele11_pt);
	  } else {
	    min1 = deltaR_01_12;
	    theDEta1 = fabs(ele01_eta-ele12_eta);
	    theDPhi1 = ele01.DeltaPhi(ele12);	  
	    theDPtOverPt1 = (ele01_pt-ele12_pt)/(ele01_pt+ele12_pt);
	  }
	  if ( deltaR_02_11<deltaR_02_12 ) {
	    min2 = deltaR_02_11;
	    theDEta2 = fabs(ele02_eta-ele11_eta);
	    theDPhi2 = ele02.DeltaPhi(ele11);
	    theDPtOverPt2 = (ele02_pt-ele11_pt)/(ele02_pt+ele11_pt);
	  } else { 
	    min2 = deltaR_02_12;
	    theDEta2 = fabs(ele02_eta-ele12_eta);
	    theDPhi2 = ele02.DeltaPhi(ele12);
	    theDPtOverPt2 = (ele02_pt-ele12_pt)/(ele02_pt+ele12_pt);
	  }
	  
	  if (min1>=min2 && min1>maxMinDREle) { 
	    maxMinDREle = min1;
	    maxMinDREle_dEta = theDEta1;
	    maxMinDREle_dPhi = theDPhi1;
	    maxMinDREle_dPtOverPt = theDPtOverPt1;
	  } else if (min2>min1 && min2>maxMinDREle) {   
	    maxMinDREle = min2;
	    maxMinDREle_dEta = theDEta2;
	    maxMinDREle_dPhi = theDPhi2;
	    maxMinDREle_dPtOverPt = theDPtOverPt2;
	  }
	  
	  // tracks
	  if (deltaR_tracks>maxDRTrack) maxDRTrack = deltaR_tracks;
	}

	float drgen = dRgen(thisB0); 
	goodTrueB_dRgen.push_back(drgen);      
	goodTrueB_maxMinDREle.push_back(maxMinDREle);
	goodTrueB_maxMinDREle_dEta.push_back(maxMinDREle_dEta);
	goodTrueB_maxMinDREle_dPhi.push_back(maxMinDREle_dPhi);
	goodTrueB_maxMinDREle_dPtOverPt.push_back(maxMinDREle_dPtOverPt);
	goodTrueB_maxDRTrack.push_back(maxDRTrack);
      }
    }
    

    // ----------------------------------------------------------------
    // Bs from combinatorics: duplication of electrons and tracks for events with at least 2 comb Bs
    if (goodCombBs.size()>1) { 
      for (u_int iB0=0; iB0<(goodCombBs.size()-1); iB0++) {        

	float maxMinDREle = -999.;
	float maxMinDREle_dEta = -999.;
	float maxMinDREle_dPhi = -999.;
	float maxMinDREle_dPtOverPt = -999.;
	float maxDRTrack = -999.;
	
	int thisB0      = goodCombBs[iB0];
	int ele01_idx   = BToKEE_l1Idx[thisB0];
	int ele02_idx   = BToKEE_l2Idx[thisB0];
	int k0_idx      = BToKEE_kIdx[thisB0];
	float ele01_pt  = Electron_pt[ele01_idx];
	float ele02_pt  = Electron_pt[ele02_idx];
	float ele01_eta = Electron_eta[ele01_idx];
	float ele02_eta = Electron_eta[ele02_idx];
	float ele01_phi = Electron_phi[ele01_idx];
	float ele02_phi = Electron_phi[ele02_idx];
	float k0_pt     = ProbeTracks_pt[k0_idx];
	float k0_eta    = ProbeTracks_eta[k0_idx];
	float k0_phi    = ProbeTracks_phi[k0_idx];
	float BFitprob0 = BToKEE_svprob[thisB0];
	
	for (u_int iB1=iB0+1; iB1<goodCombBs.size(); iB1++) {      
	  
	  int thisB1      = goodCombBs[iB1];
	  int ele11_idx   = BToKEE_l1Idx[thisB1];
	  int ele12_idx   = BToKEE_l2Idx[thisB1];
	  int k1_idx      = BToKEE_kIdx[thisB1];
	  float ele11_pt  = Electron_pt[ele11_idx];
	  float ele12_pt  = Electron_pt[ele12_idx];
	  float ele11_eta = Electron_eta[ele11_idx];
	  float ele12_eta = Electron_eta[ele12_idx];
	  float ele11_phi = Electron_phi[ele11_idx];
	  float ele12_phi = Electron_phi[ele12_idx];
	  float k1_pt     = ProbeTracks_pt[k1_idx];
	  float k1_eta    = ProbeTracks_eta[k1_idx];
	  float k1_phi    = ProbeTracks_phi[k1_idx];
	  float BFitprob1 = BToKEE_svprob[thisB1];
	  
	  TVector3 ele01(0.,0.,0.);
	  TVector3 ele02(0.,0.,0.);
	  TVector3 ele11(0.,0.,0.);
	  TVector3 ele12(0.,0.,0.);
	  ele01.SetPtEtaPhi(ele01_pt, ele01_eta, ele01_phi);
	  ele02.SetPtEtaPhi(ele02_pt, ele02_eta, ele02_phi);
	  ele11.SetPtEtaPhi(ele11_pt, ele11_eta, ele11_phi);
	  ele12.SetPtEtaPhi(ele12_pt, ele12_eta, ele12_phi);
	  float deltaR_01_11 = ele01.DeltaR(ele11); 
	  float deltaR_01_12 = ele01.DeltaR(ele12); 
	  float deltaR_02_11 = ele02.DeltaR(ele11); 
	  float deltaR_02_12 = ele02.DeltaR(ele12); 
	  
	  TVector3 track0(0.,0.,0.);
	  TVector3 track1(0.,0.,0.);
	  track0.SetPtEtaPhi(k0_pt, k0_eta, k0_phi);
	  track1.SetPtEtaPhi(k1_pt, k1_eta, k1_phi);
	  float deltaR_tracks = track0.DeltaR(track1);
	  
	  // electrons
	  float min1 = -999;
	  float min2 = -999;
	  float theDEta1 = -999.;
	  float theDEta2 = -999.;
	  float theDPhi1 = -999.;
	  float theDPhi2 = -999.;
	  float theDPtOverPt1 = -999.;
	  float theDPtOverPt2 = -999.;
	  if ( deltaR_01_11<deltaR_01_12 ) { 
	    min1 = deltaR_01_11;
	    theDEta1 = fabs(ele01_eta-ele11_eta);
	    theDPhi1 = ele01.DeltaPhi(ele11);
	    theDPtOverPt1 = (ele01_pt-ele11_pt)/(ele01_pt+ele11_pt);
	  } else {
	    min1 = deltaR_01_12;
	    theDEta1 = fabs(ele01_eta-ele12_eta);
	    theDPhi1 = ele01.DeltaPhi(ele12);	  
	    theDPtOverPt1 = (ele01_pt-ele12_pt)/(ele01_pt+ele12_pt);
	  }
	  if ( deltaR_02_11<deltaR_02_12 ) {
	    min2 = deltaR_02_11;
	    theDEta2 = fabs(ele02_eta-ele11_eta);
	    theDPhi2 = ele02.DeltaPhi(ele11);
	    theDPtOverPt2 = (ele02_pt-ele11_pt)/(ele02_pt+ele11_pt);
	  } else { 
	    min2 = deltaR_02_12;
	    theDEta2 = fabs(ele02_eta-ele12_eta);
	    theDPhi2 = ele02.DeltaPhi(ele12);
	    theDPtOverPt2 = (ele02_pt-ele12_pt)/(ele02_pt+ele12_pt);
	  }
	  
	  if (min1>=min2 && min1>maxMinDREle) { 
	    maxMinDREle = min1;
	    maxMinDREle_dEta = theDEta1;
	    maxMinDREle_dPhi = theDPhi1;
	    maxMinDREle_dPtOverPt = theDPtOverPt1;
	  } else if (min2>min1 && min2>maxMinDREle) {   
	    maxMinDREle = min2;
	    maxMinDREle_dEta = theDEta2;
	    maxMinDREle_dPhi = theDPhi2;
	    maxMinDREle_dPtOverPt = theDPtOverPt2;
	  }
	  
	  // tracks
	  if (deltaR_tracks>maxDRTrack) maxDRTrack = deltaR_tracks;
	}
	
	goodCombB_maxMinDREle.push_back(maxMinDREle);
	goodCombB_maxMinDREle_dEta.push_back(maxMinDREle_dEta);
	goodCombB_maxMinDREle_dPhi.push_back(maxMinDREle_dPhi);
	goodCombB_maxMinDREle_dPtOverPt.push_back(maxMinDREle_dPtOverPt);
	goodCombB_maxDRTrack.push_back(maxDRTrack);
      }
    }


    // -----------------------------------------------------------------------------
    // 
    // Choice of the best electrons for the triplet -> the electrons in good Bs are ~always the same. 
    // Q: Which are the better reconstructed ones? 
    // 
    // a) Best reco-gen DR, as a reference.   bestMatchedB = -1 ==> no "true" B 
    int bestMatchedB = -1;
    float dRtotMin=999.;
    for (u_int iB=0; iB<goodTrueBs.size(); iB++) {        
      int thisB    = goodTrueBs[iB];
      int ele1_idx = BToKEE_l1Idx[thisB];
      int ele2_idx = BToKEE_l2Idx[thisB];
      float dR1 = dRRecoGenEle(ele1_idx);
      float dR2 = dRRecoGenEle(ele2_idx);
      float dRtot = dR1+dR2;
      if (dRtot<dRtotMin) {          // if more than 1 found with same dRtot keep the first 
	dRtotMin     = dRtot;
	bestMatchedB = thisB;
      }
    }
    if (bestMatchedB==-1 && goodTrueBs.size()>0) cout << "Error! bestMatchedB not found" << endl;

    // b) Distributions related to the best matched B
    bestMatch_Bmass  = -999.;
    bestMatch_SvProb = -999.;
    bestMatch_XYSig  = -999.;
    bestMatch_Cos2D  = -999.;
    bestMatch_PtSum  = -999.;
    bestMatch_KPt    = -999.;
    bestMatch_KEta    = -999.;
    bestMatch_Ele1Pt  = -999.;
    bestMatch_Ele2Pt  = -999.;
    bestMatch_MinPt   = -999.;
    bestMatch_Ele1Eta = -999.;
    bestMatch_Ele2Eta = -999.;
    bestMatch_Ele1pfmva  = -999.;
    bestMatch_Ele2pfmva  = -999.;
    bestMatch_Ele1lptmva = -999.;
    bestMatch_Ele2lptmva = -999.;
    bestMatch_maxDrRecoGen = -999.; 
    bestMatch_minDrRecoGen = -999.;
    bestMatch_drRecoGenK   = -999.;
    bestMatch_maxDrRecoGenFromB = -999.; 
    bestMatch_minDrRecoGenFromB = -999.;
    bestMatch_drRecoGenFromBK   = -999.;
    //
    if (goodTrueBs.size()>0) {
      int ele1BM_idx   = BToKEE_l1Idx[bestMatchedB];
      int ele2BM_idx   = BToKEE_l2Idx[bestMatchedB];
      int kBM_idx      = BToKEE_kIdx[bestMatchedB];
      bestMatch_Bmass  = BToKEE_fit_mass[bestMatchedB];
      bestMatch_SvProb = BToKEE_svprob[bestMatchedB];
      bestMatch_XYSig  = BToKEE_l_xy[bestMatchedB]/BToKEE_l_xy_unc[bestMatchedB];
      bestMatch_Cos2D  = BToKEE_fit_cos2D[bestMatchedB]; 
      bestMatch_PtSum  = Electron_pt[ele1BM_idx] + Electron_pt[ele2BM_idx];
      bestMatch_KPt    = ProbeTracks_pt[kBM_idx];
      bestMatch_KEta   = ProbeTracks_eta[kBM_idx];
      float bestMatch_KPhi   = ProbeTracks_phi[kBM_idx];
      bestMatch_Ele1Pt = Electron_pt[ele1BM_idx];
      bestMatch_Ele2Pt = Electron_pt[ele2BM_idx];
      if (bestMatch_Ele2Pt>bestMatch_KPt)
	bestMatch_MinPt = bestMatch_KPt;
      else
	bestMatch_MinPt = bestMatch_Ele2Pt;
      bestMatch_Ele1Eta = Electron_eta[ele1BM_idx];
      bestMatch_Ele2Eta = Electron_eta[ele2BM_idx];
      float bestMatch_Ele1Phi = Electron_phi[ele1BM_idx];
      float bestMatch_Ele2Phi = Electron_phi[ele2BM_idx];
      bestMatch_Ele1pfmva  = Electron_pfmvaId[ele1BM_idx];
      bestMatch_Ele2pfmva  = Electron_pfmvaId[ele2BM_idx];
      bestMatch_Ele1lptmva = Electron_mvaId[ele1BM_idx];
      bestMatch_Ele2lptmva = Electron_mvaId[ele2BM_idx];

      // To verify match with MC-truth using the best matched gen particles 
      float dRRecoGenEle1BM = dRRecoGenEle(ele1BM_idx);
      float dRRecoGenEle2BM = dRRecoGenEle(ele2BM_idx);
      if (dRRecoGenEle1BM>dRRecoGenEle2BM) {
	bestMatch_maxDrRecoGen = dRRecoGenEle1BM;
	bestMatch_minDrRecoGen = dRRecoGenEle2BM;
      } else {
	bestMatch_maxDrRecoGen = dRRecoGenEle2BM;
	bestMatch_minDrRecoGen = dRRecoGenEle1BM;
      }
      bestMatch_drRecoGenK = dRRecoGenK(kBM_idx);

      // To verify match with MC-truth using the right particles from B
      TVector3 recoEle1(0.,0.,0.);
      TVector3 recoEle2(0.,0.,0.);
      TVector3 recoK(0.,0.,0.);
      recoEle1.SetPtEtaPhi(bestMatch_Ele1Pt,bestMatch_Ele1Eta,bestMatch_Ele1Phi);
      recoEle2.SetPtEtaPhi(bestMatch_Ele2Pt,bestMatch_Ele2Eta,bestMatch_Ele2Phi);
      recoK.SetPtEtaPhi(bestMatch_KPt,bestMatch_KEta,bestMatch_KPhi);
      float dR_GenEleFromB_RecoEle1 = recoEle1.DeltaR(genEle);
      float dR_GenPosFromB_RecoEle1 = recoEle1.DeltaR(genPos);
      float dR_GenEleFromB_RecoEle2 = recoEle2.DeltaR(genEle);
      float dR_GenPosFromB_RecoEle2 = recoEle2.DeltaR(genPos);
      float dR_GenFromB_RecoEle1    = dR_GenEleFromB_RecoEle1;
      if (dR_GenPosFromB_RecoEle1<dR_GenEleFromB_RecoEle1) dR_GenFromB_RecoEle1 = dR_GenPosFromB_RecoEle1;
      float dR_GenFromB_RecoEle2    = dR_GenEleFromB_RecoEle2;
      if (dR_GenPosFromB_RecoEle2<dR_GenEleFromB_RecoEle2) dR_GenFromB_RecoEle2 = dR_GenPosFromB_RecoEle2;
      if (dR_GenFromB_RecoEle1>dR_GenFromB_RecoEle2) {
	bestMatch_maxDrRecoGenFromB = dR_GenFromB_RecoEle1;
	bestMatch_minDrRecoGenFromB = dR_GenFromB_RecoEle2;
      } else {
	bestMatch_maxDrRecoGenFromB = dR_GenFromB_RecoEle2;
	bestMatch_minDrRecoGenFromB = dR_GenFromB_RecoEle1;
      }
      bestMatch_drRecoGenFromBK = recoK.DeltaR(genK);
    }

    // c) Try different reco criteria on matched ("true") Bs and compare with best dR with MC truth
    //    Counting per event, only in events with >=2 true Bs.
    goodTrueBs_SvProbMatch=-1;
    goodTrueBs_xySigMatch=-1;
    goodTrueBs_cos2DMatch=-1;
    goodTrueBs_ptSumMatch=-1;
    goodTrueBs_kptMatch=-1;
    //
    goodTrueBs_bestSvProb_dRmax=-1.;
    goodTrueBs_bestSvProb_dRmin=-1.;
    goodTrueBs_bestXYSig_dRmax=-1.;
    goodTrueBs_bestXYSig_dRmin=-1.;
    goodTrueBs_bestCos2d_dRmax=-1.;
    goodTrueBs_bestCos2d_dRmin=-1.;
    goodTrueBs_bestPtsum_dRmax=-1.;
    goodTrueBs_bestPtsum_dRmin=-1.;
    goodTrueBs_bestKpt_dRmax=-1.;
    goodTrueBs_bestKpt_dRmin=-1.;

    if (goodTrueBs.size()>=2) {

      int bestSvprob_onlyTrue = -1;
      int bestXYsig_onlyTrue  = -1;
      int bestCos2D_onlyTrue  = -1;
      int bestPtSum_onlyTrue  = -1;
      int bestKPt_onlyTrue    = -1;
      float svprobMax_onlyTrue  = -999.;
      float xysigMax_onlyTrue   = -999.;
      float cos2DMax_onlyTrue   = -999.;
      float ptsumMax_onlyTrue   = -999.;
      float kptMax_onlyTrue     = -999.;
    
      for (u_int iB=0; iB<goodTrueBs.size(); iB++) {        

	int thisB    = goodTrueBs[iB];
	int ele1_idx = BToKEE_l1Idx[thisB];
	int ele2_idx = BToKEE_l2Idx[thisB];
	int k_idx    = BToKEE_kIdx[thisB];
	
	float ele1_pt = Electron_pt[ele1_idx];
	float ele2_pt = Electron_pt[ele2_idx];
	float k_pt    = ProbeTracks_pt[k_idx];

	// best fit
	if (BToKEE_svprob[thisB]>svprobMax_onlyTrue) {
	  svprobMax_onlyTrue  = BToKEE_svprob[thisB];
	  bestSvprob_onlyTrue = thisB; 
	}
	
	// best XYsign
	float theXySig = BToKEE_l_xy[thisB]/BToKEE_l_xy_unc[thisB];
	if (theXySig>xysigMax_onlyTrue) {
	  xysigMax_onlyTrue  = theXySig;
	  bestXYsig_onlyTrue = thisB;
	} 
	
	// highest cos2D
	if (BToKEE_fit_cos2D[thisB]>cos2DMax_onlyTrue) {
	  cos2DMax_onlyTrue  = BToKEE_fit_cos2D[thisB];
	  bestCos2D_onlyTrue = thisB;
	} 
	
	// highest sum_pT for electrons
	float ptsum = Electron_pt[ele1_idx] + Electron_pt[ele2_idx];
	if (ptsum>ptsumMax_onlyTrue) {
	  ptsumMax_onlyTrue  = ptsum;
	  bestPtSum_onlyTrue = thisB;
	}
	
	// highest K pT
	if (k_pt>kptMax_onlyTrue) {
	  kptMax_onlyTrue  = k_pt;
	  bestKPt_onlyTrue = thisB;
	}
      } // Loop over true Bs

      // 
      goodTrueBs_SvProbMatch=0;
      goodTrueBs_xySigMatch=0;
      goodTrueBs_cos2DMatch=0;
      goodTrueBs_ptSumMatch=0;
      goodTrueBs_kptMatch=0;
      //
      // index can be the same
      if (bestSvprob_onlyTrue==bestMatchedB) goodTrueBs_SvProbMatch=1;
      if (bestXYsig_onlyTrue==bestMatchedB)  goodTrueBs_xySigMatch=1;
      if (bestCos2D_onlyTrue==bestMatchedB)  goodTrueBs_cos2DMatch=1;
      if (bestPtSum_onlyTrue==bestMatchedB)  goodTrueBs_ptSumMatch=1;
      if (bestKPt_onlyTrue==bestMatchedB)    goodTrueBs_kptMatch=1;
      //
      // index can be different, but the reco variable be the same
      if (fabs(svprobMax_onlyTrue-bestMatch_SvProb)<0.000001) { goodTrueBs_SvProbMatch=1; }
      if (fabs(xysigMax_onlyTrue-bestMatch_XYSig)<0.000001)   { goodTrueBs_xySigMatch=1;  }
      if (fabs(cos2DMax_onlyTrue-bestMatch_Cos2D)<0.000001)   { goodTrueBs_cos2DMatch=1;  }
      if (fabs(ptsumMax_onlyTrue-bestMatch_PtSum)<0.000001)   { goodTrueBs_ptSumMatch=1;  }
      if (fabs(kptMax_onlyTrue-bestMatch_KPt)<0.000001)       { goodTrueBs_kptMatch=1;    }
      //
      // index can be different, but dR be the same (i.e. when an electrons is low-pt and PF: eta/phi the same, differences in pt -> in vtx variables)
      // correcting for this, and saving dR to check the cases which still do not match
      float drtotbm = bestMatch_maxDrRecoGen+bestMatch_minDrRecoGen;
      if (goodTrueBs_SvProbMatch==0) { 
	int idx1  = BToKEE_l1Idx[bestSvprob_onlyTrue];
	int idx2  = BToKEE_l2Idx[bestSvprob_onlyTrue];
	float dr1 = dRRecoGenEle(idx1);
	float dr2 = dRRecoGenEle(idx2);
	if ( (dr1+dr2)==drtotbm ) goodTrueBs_SvProbMatch=1;
	else { 
	  if (dr1>dr2) { goodTrueBs_bestSvProb_dRmax = dr1; goodTrueBs_bestSvProb_dRmin = dr2; }
	  else { goodTrueBs_bestSvProb_dRmax = dr2; goodTrueBs_bestSvProb_dRmin =dr1; }
	}
      }
      if (goodTrueBs_xySigMatch==0) {
	int idx1  = BToKEE_l1Idx[bestXYsig_onlyTrue];
	int idx2  = BToKEE_l2Idx[bestXYsig_onlyTrue];
	float dr1 = dRRecoGenEle(idx1);
	float dr2 = dRRecoGenEle(idx2);
	if ( (dr1+dr2)==drtotbm ) goodTrueBs_xySigMatch=1;
	else { 
	  if (dr1>dr2) { goodTrueBs_bestXYSig_dRmax = dr1; goodTrueBs_bestXYSig_dRmin = dr2; }
	  else { goodTrueBs_bestXYSig_dRmax = dr2; goodTrueBs_bestXYSig_dRmin =dr1; }
	}
      }
      if (goodTrueBs_cos2DMatch==0) { 
	int idx1  = BToKEE_l1Idx[bestCos2D_onlyTrue];
	int idx2  = BToKEE_l2Idx[bestCos2D_onlyTrue];
	float dr1 = dRRecoGenEle(idx1);
	float dr2 = dRRecoGenEle(idx2);
	if ( (dr1+dr2)==drtotbm ) goodTrueBs_cos2DMatch=1;
	else {
	  if (dr1>dr2) { goodTrueBs_bestCos2d_dRmax = dr1; goodTrueBs_bestCos2d_dRmin = dr2; }
	  else { goodTrueBs_bestCos2d_dRmax = dr2; goodTrueBs_bestCos2d_dRmin =dr1; }
	}
      }
      if (goodTrueBs_ptSumMatch==0) {
	int idx1  = BToKEE_l1Idx[bestPtSum_onlyTrue];
	int idx2  = BToKEE_l2Idx[bestPtSum_onlyTrue];
	float dr1 = dRRecoGenEle(idx1);
	float dr2 = dRRecoGenEle(idx2);
	if ( (dr1+dr2)==drtotbm ) goodTrueBs_ptSumMatch=1;
	else {
	  if (dr1>dr2) { goodTrueBs_bestPtsum_dRmax = dr1; goodTrueBs_bestPtsum_dRmin = dr2; }
	  else { goodTrueBs_bestPtsum_dRmax = dr2; goodTrueBs_bestPtsum_dRmin = dr1; }
	}
      }
      if (goodTrueBs_kptMatch==0) {
	int idx1  = BToKEE_l1Idx[bestKPt_onlyTrue];
	int idx2  = BToKEE_l2Idx[bestKPt_onlyTrue];
	float dr1 = dRRecoGenEle(idx1);
	float dr2 = dRRecoGenEle(idx2);
	if ( (dr1+dr2)==drtotbm ) goodTrueBs_kptMatch=1;
	else {
	  if (dr1<dr2) { goodTrueBs_bestKpt_dRmax = dr1; goodTrueBs_bestKpt_dRmin = dr2; }
	  else { goodTrueBs_bestKpt_dRmax = dr2; goodTrueBs_bestKpt_dRmin =dr1; }
	}
      }

    }  // at least 2 true Bs

    // d) fill histos for distributions related to matched Bs: only not best matched
    if (goodTrueBs.size()>0) {

      for (u_int iB=0; iB<goodTrueBs.size(); iB++) {        
	int thisB = goodTrueBs[iB];
	
	// Remove best match
	if (thisB==bestMatchedB) continue;
	
	// Remove cases where dR is the same as for best match
	int ele1_idx = BToKEE_l1Idx[thisB];
	int ele2_idx = BToKEE_l2Idx[thisB];
	float drGenRecoEle1 = dRRecoGenEle(ele1_idx);
	float drGenRecoEle2 = dRRecoGenEle(ele2_idx);
	if ( (drGenRecoEle1+drGenRecoEle2)==(bestMatch_maxDrRecoGen+bestMatch_minDrRecoGen) ) continue;
	
	// Compute variables
	int k_idx  = BToKEE_kIdx[thisB];
	float k_pt = ProbeTracks_pt[k_idx];
	float ele_ptsum = Electron_pt[ele1_idx] + Electron_pt[ele2_idx];
	float theXySig  = BToKEE_l_xy[thisB]/BToKEE_l_xy_unc[thisB];
	
	// Check if the variable is the same as for best match. If not fill histos
	if (fabs(BToKEE_svprob[thisB]-bestMatch_SvProb)>0.000001)   goodTrueB_svProb_notBestMatch.push_back(BToKEE_svprob[thisB]);
	if (fabs(theXySig-bestMatch_XYSig)>0.000001)                goodTrueB_xySig_notBestMatch.push_back(theXySig);
	if (fabs(BToKEE_fit_cos2D[thisB]-bestMatch_Cos2D)>0.000001) goodTrueB_cos2D_notBestMatch.push_back(BToKEE_fit_cos2D[thisB]);
	if (fabs(ele_ptsum-bestMatch_PtSum)>0.000001)               goodTrueB_ptsum_notBestMatch.push_back(ele_ptsum);  
	if (fabs(k_pt-bestMatch_KPt)>0.000001)                      goodTrueB_kpt_notBestMatch.push_back(ProbeTracks_pt[k_idx]);
      }
    }

    // e) fill histos for distributions related to NOT matched Bs
    for (u_int iB=0; iB<goodCombBs.size(); iB++) {        
      int thisB    = goodCombBs[iB];
      int ele1_idx = BToKEE_l1Idx[thisB];
      int ele2_idx = BToKEE_l2Idx[thisB];
      int k_idx    = BToKEE_kIdx[thisB];

      int k_genPartIdx    = ProbeTracks_genPartIdx[k_idx];  
      int ele1_genPartIdx = Electron_genPartIdx[ele1_idx];  
      int ele2_genPartIdx = Electron_genPartIdx[ele2_idx];  
      float ptsum    = Electron_pt[ele1_idx] + Electron_pt[ele2_idx];
      float theXySig = BToKEE_l_xy[thisB]/BToKEE_l_xy_unc[thisB];
      float dRRecoGenEle1 = dRRecoGenEle(ele1_idx);
      float dRRecoGenEle2 = dRRecoGenEle(ele2_idx);
      float maxDrRecoGenEle, minDrRecoGenEle;
      float maxDrRecoGenFromBEle, minDrRecoGenFromBEle;
      float drRecoGenK = dRRecoGenK(k_idx);
      if (dRRecoGenEle1>=dRRecoGenEle2) {
	maxDrRecoGenEle = dRRecoGenEle1;
	minDrRecoGenEle = dRRecoGenEle2;
      } else {
	maxDrRecoGenEle = dRRecoGenEle2;
	minDrRecoGenEle = dRRecoGenEle1;
      }
      goodCombB_svProb.push_back(BToKEE_svprob[thisB]);	
      goodCombB_xySig.push_back(theXySig);
      goodCombB_cos2D.push_back(BToKEE_fit_cos2D[thisB]);
      goodCombB_ptsum.push_back(ptsum);
      goodCombB_kpt.push_back(ProbeTracks_pt[k_idx]);
      goodCombB_keta.push_back(ProbeTracks_eta[k_idx]);
      goodCombB_ele1pt.push_back(Electron_pt[ele1_idx]);
      goodCombB_ele2pt.push_back(Electron_pt[ele2_idx]);
      if (Electron_pt[ele2_idx]>ProbeTracks_pt[k_idx])
	goodCombB_minpt.push_back(ProbeTracks_pt[k_idx]);
      else 
	goodCombB_minpt.push_back(Electron_pt[ele2_idx]);
      goodCombB_ele1eta.push_back(Electron_eta[ele1_idx]);
      goodCombB_ele2eta.push_back(Electron_eta[ele2_idx]);
      goodCombB_ele1pfmva.push_back(Electron_pfmvaId[ele1_idx]);
      goodCombB_ele2pfmva.push_back(Electron_pfmvaId[ele2_idx]);
      goodCombB_ele1lptmva.push_back(Electron_mvaId[ele1_idx]);
      goodCombB_ele2lptmva.push_back(Electron_mvaId[ele2_idx]);
      if (ele1_genPartIdx>-0.5) goodCombB_causeEle1.push_back(0);
      else goodCombB_causeEle1.push_back(1);
      if (ele2_genPartIdx>-0.5) goodCombB_causeEle2.push_back(0);
      else goodCombB_causeEle2.push_back(1);
      if (k_genPartIdx>-0.5) goodCombB_causeK.push_back(0);
      else goodCombB_causeK.push_back(1);
      goodCombB_notmatching.push_back(whichMcB(thisB));



      // To verify match with MC-truth using the right particles from B
      TVector3 recoEle1(0.,0.,0.);
      TVector3 recoEle2(0.,0.,0.);
      TVector3 recoK(0.,0.,0.);
      recoEle1.SetPtEtaPhi(Electron_pt[ele1_idx], Electron_eta[ele1_idx], Electron_phi[ele1_idx]);
      recoEle2.SetPtEtaPhi(Electron_pt[ele2_idx], Electron_eta[ele2_idx], Electron_phi[ele2_idx]);
      recoK.SetPtEtaPhi(ProbeTracks_pt[k_idx],ProbeTracks_eta[k_idx],ProbeTracks_phi[k_idx]);
      float dR_GenFromBEle_RecoEle1 = recoEle1.DeltaR(genEle);
      float dR_GenFromBPos_RecoEle1 = recoEle1.DeltaR(genPos);
      float dR_GenFromBEle_RecoEle2 = recoEle2.DeltaR(genEle);
      float dR_GenFromBPos_RecoEle2 = recoEle2.DeltaR(genPos);
      float dR_GenFromB_RecoEle1    = dR_GenFromBEle_RecoEle1;
      if (dR_GenFromBPos_RecoEle1<dR_GenFromBEle_RecoEle1) dR_GenFromB_RecoEle1 = dR_GenFromBPos_RecoEle1;
      float dR_GenFromB_RecoEle2    = dR_GenFromBEle_RecoEle2;
      if (dR_GenFromBPos_RecoEle2<dR_GenFromBEle_RecoEle2) dR_GenFromB_RecoEle2 = dR_GenFromBPos_RecoEle2;
      if (dR_GenFromB_RecoEle1>dR_GenFromB_RecoEle2) {
        maxDrRecoGenFromBEle = dR_GenFromB_RecoEle1;
        minDrRecoGenFromBEle = dR_GenFromB_RecoEle2;
      } else {
	maxDrRecoGenFromBEle = dR_GenFromB_RecoEle2;
	minDrRecoGenFromBEle = dR_GenFromB_RecoEle1;
      }
      float drRecoGenFromBK = recoK.DeltaR(genK);

      // with the associated gen-particle
      goodCombB_maxDrRecoGen.push_back(maxDrRecoGenEle);
      goodCombB_minDrRecoGen.push_back(minDrRecoGenEle);
      goodCombB_drRecoGenK.push_back(drRecoGenK);
      // with the gen-particle from B
      goodCombB_maxDrRecoGenFromB.push_back(maxDrRecoGenFromBEle);
      goodCombB_minDrRecoGenFromB.push_back(minDrRecoGenFromBEle);
      goodCombB_drRecoGenFromBK.push_back(drRecoGenFromBK);
    }

    // f) how often the best B according to reco criteria is NOT matched (not only not the best, but not matched at all)
    //    only cases with at least 1 true and 1 combinatorics selected 
    bestSvProbMatch=-999;
    bestSvProbMatch_notmatching=-999;
    bestSvProbMatchCat0=-999;
    bestSvProbMatchCat1=-999;
    bestSvProbMatchCat2=-999;
    bestSvProbMatchCatNew0=-999;
    bestSvProbMatchCatNew1=-999;
    bestSvProbMatchCatNew2=-999;
    bestSvProbMatch_causeEle1=-999;
    bestSvProbMatch_causeEle2=-999;
    bestSvProbMatch_causeK=-999;
    bestSvProbMatch_notok_ele1pt=-999;
    bestSvProbMatch_notok_ele2pt=-999;
    bestSvProbMatch_notok_kpt=-999;
    bestSvProbMatch_notok_ele1eta=-999;
    bestSvProbMatch_notok_ele2eta=-999;
    bestSvProbMatch_notok_keta=-999;
    bestSvProbMatch_ok_ele1pt=-999;
    bestSvProbMatch_ok_ele2pt=-999;
    bestSvProbMatch_ok_kpt=-999;
    bestSvProbMatch_ok_ele1eta=-999;
    bestSvProbMatch_ok_ele2eta=-999;
    bestSvProbMatch_ok_keta=-999;
    numberBetterSvProbTriplets=-999;
    bestSvProbMatch_second=-999;
    bestSvProbMatchCat0_second=-999;
    bestSvProbMatchCat1_second=-999;
    bestSvProbMatchCat2_second=-999;
    bestSvProbMatchCatNew0_second=-999;
    bestSvProbMatchCatNew1_second=-999;
    bestSvProbMatchCatNew2_second=-999;

    bestXYsigMatch=-999;
    bestXYsigMatch_notmatching=-999;
    bestXYsigMatchCat0=-999;
    bestXYsigMatchCat1=-999;
    bestXYsigMatchCat2=-999;
    bestXYsigMatchCatNew0=-999;
    bestXYsigMatchCatNew1=-999;
    bestXYsigMatchCatNew2=-999;
    bestXYsigMatch_causeEle1=-999;
    bestXYsigMatch_causeEle2=-999;
    bestXYsigMatch_causeK=-999;
    bestXYsigMatch_notok_ele1pt=-999;
    bestXYsigMatch_notok_ele2pt=-999;
    bestXYsigMatch_notok_kpt=-999;
    bestXYsigMatch_notok_minpt=-999;
    bestXYsigMatch_notok_ele1eta=-999;
    bestXYsigMatch_notok_ele2eta=-999;
    bestXYsigMatch_notok_keta=-999;
    bestXYsigMatch_notok_pfmva1=-999;
    bestXYsigMatch_notok_pfmva2=-999;
    bestXYsigMatch_notok_lptmva1=-999;
    bestXYsigMatch_notok_lptmva2=-999;
    bestXYsigMatch_notok_costhetaSK=-999;
    bestXYsigMatch_notok_costhetaSKCS=-999;
    bestXYsigMatch_notok_costhetaL=-999;
    bestXYsigMatch_ok_ele1pt=-999;
    bestXYsigMatch_ok_ele2pt=-999;
    bestXYsigMatch_ok_kpt=-999;
    bestXYsigMatch_ok_minpt=-999;
    bestXYsigMatch_ok_ele1eta=-999;
    bestXYsigMatch_ok_ele2eta=-999;
    bestXYsigMatch_ok_keta=-999;
    bestXYsigMatch_ok_pfmva1=-999;
    bestXYsigMatch_ok_pfmva2=-999;
    bestXYsigMatch_ok_lptmva1=-999;
    bestXYsigMatch_ok_lptmva2=-999;
    bestXYsigMatch_ok_costhetaSK=-999;
    bestXYsigMatch_ok_costhetaSKCS=-999;
    bestXYsigMatch_ok_costhetaL=-999;
    bestXYsigMatch_ok_costhetaSK_gen=-999;
    numberBetterXYsigTriplets=-999;
    bestXYsigMatch_second=-999;
    bestXYsigMatchCat0_second=-999;
    bestXYsigMatchCat1_second=-999;
    bestXYsigMatchCat2_second=-999;
    bestXYsigMatchCatNew0_second=-999;
    bestXYsigMatchCatNew1_second=-999;
    bestXYsigMatchCatNew2_second=-999;

    bestCos2DMatch=-999;
    bestCos2DMatchCat0=-999;
    bestCos2DMatchCat1=-999;
    bestCos2DMatchCat2=-999;
    bestCos2DMatchCatNew0=-999;
    bestCos2DMatchCatNew1=-999;
    bestCos2DMatchCatNew2=-999;
    bestCos2DMatch_causeEle1=-999;
    bestCos2DMatch_causeEle2=-999;
    bestCos2DMatch_causeK=-999;
    bestCos2DMatch_notok_ele1pt=-999;
    bestCos2DMatch_notok_ele2pt=-999;
    bestCos2DMatch_notok_kpt=-999;
    bestCos2DMatch_notok_ele1eta=-999;
    bestCos2DMatch_notok_ele2eta=-999;
    bestCos2DMatch_notok_keta=-999;
    bestCos2DMatch_ok_ele1pt=-999;
    bestCos2DMatch_ok_ele2pt=-999;
    bestCos2DMatch_ok_kpt=-999;
    bestCos2DMatch_ok_ele1eta=-999;
    bestCos2DMatch_ok_ele2eta=-999;
    bestCos2DMatch_ok_keta=-999;
    numberBetterCos2DTriplets=-999;
    bestCos2DMatch_second=-999;
    bestCos2DMatchCat0_second=-999;
    bestCos2DMatchCat1_second=-999;
    bestCos2DMatchCat2_second=-999;
    bestCos2DMatchCatNew0_second=-999;
    bestCos2DMatchCatNew1_second=-999;
    bestCos2DMatchCatNew2_second=-999;

    bestAllPtSumMatch=-999;
    bestAllPtSumMatchCat0=-999;
    bestAllPtSumMatchCat1=-999;
    bestAllPtSumMatchCat2=-999;
    bestAllPtSumMatchCatNew0=-999;
    bestAllPtSumMatchCatNew1=-999;
    bestAllPtSumMatchCatNew2=-999;
    bestAllPtSumMatch_causeEle1=-999;
    bestAllPtSumMatch_causeEle2=-999;
    bestAllPtSumMatch_causeK=-999;
    numberBetterAllPtSumTriplets=-999;  
    bestAllPtSumMatch_second=-999;
    bestAllPtSumMatchCat0_second=-999;
    bestAllPtSumMatchCat1_second=-999;
    bestAllPtSumMatchCat2_second=-999;
    bestAllPtSumMatchCatNew0_second=-999;
    bestAllPtSumMatchCatNew1_second=-999;
    bestAllPtSumMatchCatNew2_second=-999;

    //if (goodTrueBs.size()>0 && goodCombBs.size()>0) {
    if (1) {
      
      int bestSvprob_all  = -1;
      int bestSvprob_cat0 = -1;
      int bestSvprob_cat1 = -1;
      int bestSvprob_cat2 = -1;
      int bestSvprob_catNew0 = -1;
      int bestSvprob_catNew1 = -1;
      int bestSvprob_catNew2 = -1;
      int bestXYsig_all   = -1;
      int bestXYsig_cat0  = -1;
      int bestXYsig_cat1  = -1;
      int bestXYsig_cat2  = -1;
      int bestXYsig_catNew0  = -1;
      int bestXYsig_catNew1  = -1;
      int bestXYsig_catNew2  = -1;
      int bestCos2D_all   = -1;
      int bestCos2D_cat0 = -1;
      int bestCos2D_cat1 = -1;
      int bestCos2D_cat2 = -1;
      int bestCos2D_catNew0 = -1;
      int bestCos2D_catNew1 = -1;
      int bestCos2D_catNew2 = -1;
      int bestAllPtSum_all   = -1;
      int bestAllPtSum_cat0 = -1;
      int bestAllPtSum_cat1 = -1;
      int bestAllPtSum_cat2 = -1;
      int bestAllPtSum_catNew0 = -1;
      int bestAllPtSum_catNew1 = -1;
      int bestAllPtSum_catNew2 = -1;
      //
      int secondSvprob_all  = -1;
      int secondSvprob_cat0 = -1;
      int secondSvprob_cat1 = -1;
      int secondSvprob_cat2 = -1;
      int secondSvprob_catNew0 = -1;
      int secondSvprob_catNew1 = -1;
      int secondSvprob_catNew2 = -1;
      int secondXYsig_all  = -1;
      int secondXYsig_cat0 = -1;
      int secondXYsig_cat1 = -1;
      int secondXYsig_cat2 = -1;
      int secondXYsig_catNew0 = -1;
      int secondXYsig_catNew1 = -1;
      int secondXYsig_catNew2 = -1;
      int secondCos2D_all  = -1;
      int secondCos2D_cat0 = -1;
      int secondCos2D_cat1 = -1;
      int secondCos2D_cat2 = -1;
      int secondCos2D_catNew0 = -1;
      int secondCos2D_catNew1 = -1;
      int secondCos2D_catNew2 = -1;
      int secondAllPtSum_all  = -1;
      int secondAllPtSum_cat0 = -1;
      int secondAllPtSum_cat1 = -1;
      int secondAllPtSum_cat2 = -1;
      int secondAllPtSum_catNew0 = -1;
      int secondAllPtSum_catNew1 = -1;
      int secondAllPtSum_catNew2 = -1;
      //
      float svprobBest_all  = -999.;
      float svprobBest_cat0 = -999.;
      float svprobBest_cat1 = -999.;
      float svprobBest_cat2 = -999.;
      float svprobBest_catNew0 = -999.;
      float svprobBest_catNew1 = -999.;
      float svprobBest_catNew2 = -999.;
      float xysigBest_all   = -999.;
      float xysigBest_cat0  = -999.;
      float xysigBest_cat1  = -999.;
      float xysigBest_cat2  = -999.;
      float xysigBest_catNew0 = -999.;
      float xysigBest_catNew1 = -999.;
      float xysigBest_catNew2 = -999.;
      float cos2DBest_all   = -999.;
      float cos2DBest_cat0  = -999.;
      float cos2DBest_cat1  = -999.;
      float cos2DBest_cat2  = -999.;
      float cos2DBest_catNew0 = -999.;
      float cos2DBest_catNew1 = -999.;
      float cos2DBest_catNew2 = -999.;
      float allptsumBest_all   = -999.;
      float allptsumBest_cat0  = -999.;
      float allptsumBest_cat1  = -999.;
      float allptsumBest_cat2  = -999.;
      float allptsumBest_catNew0 = -999.;
      float allptsumBest_catNew1 = -999.;
      float allptsumBest_catNew2 = -999.;
      //
      float svprobSecond_all  = -999.;
      float svprobSecond_cat0 = -999.;
      float svprobSecond_cat1 = -999.;
      float svprobSecond_cat2 = -999.;
      float svprobSecond_catNew0 = -999.;
      float svprobSecond_catNew1 = -999.;
      float svprobSecond_catNew2 = -999.;
      float xysigSecond_all     = -999.;
      float xysigSecond_cat0    = -999.;
      float xysigSecond_cat1    = -999.;
      float xysigSecond_cat2    = -999.;
      float xysigSecond_catNew0 = -999.;
      float xysigSecond_catNew1 = -999.;
      float xysigSecond_catNew2 = -999.;
      float cos2DSecond_all     = -999.;
      float cos2DSecond_cat0    = -999.;
      float cos2DSecond_cat1    = -999.;
      float cos2DSecond_cat2    = -999.;
      float cos2DSecond_catNew0 = -999.;
      float cos2DSecond_catNew1 = -999.;
      float cos2DSecond_catNew2 = -999.;
      float allptsumSecond_all   = -999.;
      float allptsumSecond_cat0  = -999.;
      float allptsumSecond_cat1  = -999.;
      float allptsumSecond_cat2  = -999.;
      float allptsumSecond_catNew0 = -999.;
      float allptsumSecond_catNew1 = -999.;
      float allptsumSecond_catNew2 = -999.;
      //
      //
      float svprobBest_all_matched   = -999.;
      float xysigBest_all_matched    = -999.;
      float cos2DBest_all_matched    = -999.;
      float allPtSumBest_all_matched = -999.;

      for (u_int iB=0; iB<goodBs.size(); iB++) {        
	int thisB    = goodBs[iB];
	int ele1_idx = BToKEE_l1Idx[thisB];
	int ele2_idx = BToKEE_l2Idx[thisB];
	int k_idx    = BToKEE_kIdx[thisB];
	
	int whichCat = -1;
	if (Electron_pfmvaId[ele1_idx]<20 && Electron_pfmvaId[ele2_idx]<20)  whichCat = 0;
	else if (Electron_mvaId[ele1_idx]<20 && Electron_mvaId[ele2_idx]<20) whichCat = 2;
	else if ( (Electron_pfmvaId[ele1_idx]<20 && Electron_mvaId[ele2_idx]<20) || (Electron_pfmvaId[ele2_idx]<20 && Electron_mvaId[ele1_idx]<20) ) whichCat = 1;

	int whichCatNew = -1;
	if (Electron_pfmvaId[ele2_idx]<20)  whichCatNew = 0;
	else if (Electron_mvaId[ele2_idx]<20 && Electron_mvaId[ele2_idx]>2) whichCatNew = 1;
	else if (Electron_mvaId[ele2_idx]<20 && Electron_mvaId[ele2_idx]<2) whichCatNew = 2;
	if (whichCatNew<0) cout << "problem with new cat" << endl;

	// best fit
	if (BToKEE_svprob[thisB]>=svprobBest_all) {
	  svprobSecond_all = svprobBest_all;
	  secondSvprob_all = bestSvprob_all;
	  svprobBest_all   = BToKEE_svprob[thisB];
	  bestSvprob_all   = thisB; 
	} else if (BToKEE_svprob[thisB]>svprobSecond_all) {
	  svprobSecond_all = BToKEE_svprob[thisB];
	  secondSvprob_all = thisB; 
	}

	if (whichCat==0 && BToKEE_svprob[thisB]>svprobBest_cat0) {
	  svprobSecond_cat0 = svprobBest_cat0;
	  secondSvprob_cat0 = bestSvprob_cat0;
	  svprobBest_cat0   = BToKEE_svprob[thisB];
	  bestSvprob_cat0   = thisB; 
	} else if (whichCat==0 && BToKEE_svprob[thisB]>svprobSecond_cat0) {
	  svprobSecond_cat0 = BToKEE_svprob[thisB];
	  secondSvprob_cat0 = thisB; 
	}

	if (whichCat==1 && BToKEE_svprob[thisB]>svprobBest_cat1) {
	  svprobSecond_cat1 = svprobBest_cat1;
	  secondSvprob_cat1 = bestSvprob_cat1;
	  svprobBest_cat1   = BToKEE_svprob[thisB];
	  bestSvprob_cat1   = thisB; 
	} else if (whichCat==1 && BToKEE_svprob[thisB]>svprobSecond_cat1) {
	  svprobSecond_cat1 = BToKEE_svprob[thisB];
	  secondSvprob_cat1 = thisB; 
	}

	if (whichCat==2 && BToKEE_svprob[thisB]>svprobBest_cat2) {
	  svprobSecond_cat2 = svprobBest_cat2;
	  secondSvprob_cat2 = bestSvprob_cat2;
	  svprobBest_cat2   = BToKEE_svprob[thisB];
	  bestSvprob_cat2   = thisB; 
	} else if (whichCat==2 && BToKEE_svprob[thisB]>svprobSecond_cat2) {
	  svprobSecond_cat2 = BToKEE_svprob[thisB];
	  secondSvprob_cat2 = thisB; 
	}

	if (whichCatNew==0 && BToKEE_svprob[thisB]>svprobBest_catNew0) {
	  svprobSecond_catNew0 = svprobBest_catNew0;
	  secondSvprob_catNew0 = bestSvprob_catNew0;
	  svprobBest_catNew0   = BToKEE_svprob[thisB];
	  bestSvprob_catNew0   = thisB; 
	} else if (whichCatNew==0 && BToKEE_svprob[thisB]>svprobSecond_catNew0) {
	  svprobSecond_catNew0 = BToKEE_svprob[thisB];
	  secondSvprob_catNew0 = thisB; 
	}

	if (whichCatNew==1 && BToKEE_svprob[thisB]>svprobBest_catNew1) {
	  svprobSecond_catNew1 = svprobBest_catNew1;
	  secondSvprob_catNew1 = bestSvprob_catNew1;
	  svprobBest_catNew1   = BToKEE_svprob[thisB];
	  bestSvprob_catNew1   = thisB; 
	} else if (whichCatNew==1 && BToKEE_svprob[thisB]>svprobSecond_catNew1) {
	  svprobSecond_catNew1 = BToKEE_svprob[thisB];
	  secondSvprob_catNew1 = thisB; 
	}

	if (whichCatNew==2 && BToKEE_svprob[thisB]>svprobBest_catNew2) {
	  svprobSecond_catNew2 = svprobBest_catNew2;
	  secondSvprob_catNew2 = bestSvprob_catNew2;
	  svprobBest_catNew2   = BToKEE_svprob[thisB];
	  bestSvprob_catNew2   = thisB; 
	} else if (whichCatNew==2 && BToKEE_svprob[thisB]>svprobSecond_catNew2) {
	  svprobSecond_catNew2 = BToKEE_svprob[thisB];
	  secondSvprob_catNew2 = thisB; 
	}

	// best XYsign
	float theXySig = BToKEE_l_xy[thisB]/BToKEE_l_xy_unc[thisB];
	if (theXySig>=xysigBest_all) {
	  xysigSecond_all = xysigBest_all;
	  secondXYsig_all = bestXYsig_all;
	  xysigBest_all   = theXySig;
	  bestXYsig_all   = thisB; 
	} else if (theXySig>xysigSecond_all) {
	  xysigSecond_all = theXySig;
	  secondXYsig_all = thisB; 
	}

	if (whichCat==0 && theXySig>=xysigBest_cat0) {
	  xysigSecond_cat0 = xysigBest_cat0;
	  secondXYsig_cat0 = bestXYsig_cat0;
	  xysigBest_cat0   = theXySig;
	  bestXYsig_cat0   = thisB; 
	} else if (whichCat==0 && theXySig>xysigSecond_cat0) {
	  xysigSecond_cat0 = theXySig;
	  secondXYsig_cat0 = thisB; 
	}

	if (whichCat==1 && theXySig>=xysigBest_cat1) {
	  xysigSecond_cat1 = xysigBest_cat1;
	  secondXYsig_cat1 = bestXYsig_cat1;
	  xysigBest_cat1   = theXySig;
	  bestXYsig_cat1   = thisB; 
	} else if (whichCat==1 && theXySig>xysigSecond_cat1) {
	  xysigSecond_cat1 = theXySig;
	  secondXYsig_cat1 = thisB; 
	}

	if (whichCat==2 && theXySig>=xysigBest_cat2) {
	  xysigSecond_cat2 = xysigBest_cat2;
	  secondXYsig_cat2 = bestXYsig_cat2;
	  xysigBest_cat2   = theXySig;
	  bestXYsig_cat2   = thisB; 
	} else if (whichCat==2 && theXySig>xysigSecond_cat2) {
	  xysigSecond_cat2 = theXySig;
	  secondXYsig_cat2 = thisB; 
	}

	if (whichCatNew==0 && theXySig>=xysigBest_catNew0) {
	  xysigSecond_catNew0 = xysigBest_catNew0;
	  secondXYsig_catNew0 = bestXYsig_catNew0;
	  xysigBest_catNew0   = theXySig;
	  bestXYsig_catNew0   = thisB; 
	} else if (whichCatNew==0 && theXySig>xysigSecond_catNew0) {
	  xysigSecond_catNew0 = theXySig;
	  secondXYsig_catNew0 = thisB; 
	}

	if (whichCatNew==1 && theXySig>=xysigBest_catNew1) {
	  xysigSecond_catNew1 = xysigBest_catNew1;
	  secondXYsig_catNew1 = bestXYsig_catNew1;
	  xysigBest_catNew1   = theXySig;
	  bestXYsig_catNew1   = thisB; 
	} else if (whichCatNew==1 && theXySig>xysigSecond_catNew1) {
	  xysigSecond_catNew1 = theXySig;
	  secondXYsig_catNew1 = thisB; 
	}

	if (whichCatNew==2 && theXySig>=xysigBest_catNew2) {
	  xysigSecond_catNew2 = xysigBest_catNew2;
	  secondXYsig_catNew2 = bestXYsig_catNew2;
	  xysigBest_catNew2   = theXySig;
	  bestXYsig_catNew2   = thisB; 
	} else if (whichCatNew==2 && theXySig>xysigSecond_catNew2) {
	  xysigSecond_catNew2 = theXySig;
	  secondXYsig_catNew2 = thisB; 
	}

	// highest cos2D
	if (BToKEE_fit_cos2D[thisB]>=cos2DBest_all) {
	  cos2DSecond_all = cos2DBest_all;
	  secondCos2D_all = bestCos2D_all;
	  cos2DBest_all   = BToKEE_fit_cos2D[thisB];
	  bestCos2D_all   = thisB; 
	} else if (BToKEE_fit_cos2D[thisB]>cos2DSecond_all) {
	  cos2DSecond_all = BToKEE_fit_cos2D[thisB];
	  secondCos2D_all = thisB; 
	}

	if (whichCat==0 && BToKEE_fit_cos2D[thisB]>=cos2DBest_cat0) {
	  cos2DSecond_cat0 = cos2DBest_cat0;
	  secondCos2D_cat0 = bestCos2D_cat0;
	  cos2DBest_cat0   = BToKEE_fit_cos2D[thisB];
	  bestCos2D_cat0   = thisB; 
	} else if (whichCat==0 && BToKEE_fit_cos2D[thisB]>cos2DSecond_cat0) {
	  cos2DSecond_cat0 = BToKEE_fit_cos2D[thisB];
	  secondCos2D_cat0 = thisB; 
	}

	if (whichCat==1 && BToKEE_fit_cos2D[thisB]>=cos2DBest_cat1) {
	  cos2DSecond_cat1 = cos2DBest_cat1;
	  secondCos2D_cat1 = bestCos2D_cat1;
	  cos2DBest_cat1   = BToKEE_fit_cos2D[thisB];
	  bestCos2D_cat1   = thisB; 
	} else if (whichCat==1 && BToKEE_fit_cos2D[thisB]>cos2DSecond_cat1) {
	  cos2DSecond_cat1 = BToKEE_fit_cos2D[thisB];
	  secondCos2D_cat1 = thisB; 
	}
	
	if (whichCat==2 && BToKEE_fit_cos2D[thisB]>=cos2DBest_cat2) {
	  cos2DSecond_cat2 = cos2DBest_cat2;
	  secondCos2D_cat2 = bestCos2D_cat2;
	  cos2DBest_cat2   = BToKEE_fit_cos2D[thisB];
	  bestCos2D_cat2   = thisB; 
	} else if (whichCat==2 && BToKEE_fit_cos2D[thisB]>cos2DSecond_cat2) {
	  cos2DSecond_cat2 = BToKEE_fit_cos2D[thisB];
	  secondCos2D_cat2 = thisB; 
	}

	if (whichCatNew==0 && BToKEE_fit_cos2D[thisB]>=cos2DBest_catNew0) {
	  cos2DSecond_catNew0 = cos2DBest_catNew0;
	  secondCos2D_catNew0 = bestCos2D_catNew0;
	  cos2DBest_catNew0   = BToKEE_fit_cos2D[thisB];
	  bestCos2D_catNew0   = thisB; 
	} else if (whichCatNew==0 && BToKEE_fit_cos2D[thisB]>cos2DSecond_catNew0) {
	  cos2DSecond_catNew0 = BToKEE_fit_cos2D[thisB];
	  secondCos2D_catNew0 = thisB; 
	}

	if (whichCatNew==1 && BToKEE_fit_cos2D[thisB]>=cos2DBest_catNew1) {
	  cos2DSecond_catNew1 = cos2DBest_catNew1;
	  secondCos2D_catNew1 = bestCos2D_catNew1;
	  cos2DBest_catNew1   = BToKEE_fit_cos2D[thisB];
	  bestCos2D_catNew1   = thisB; 
	} else if (whichCatNew==1 && BToKEE_fit_cos2D[thisB]>cos2DSecond_catNew1) {
	  cos2DSecond_catNew1 = BToKEE_fit_cos2D[thisB];
	  secondCos2D_catNew1 = thisB; 
	}

	if (whichCatNew==2 && BToKEE_fit_cos2D[thisB]>=cos2DBest_catNew2) {
	  cos2DSecond_catNew2 = cos2DBest_catNew2;
	  secondCos2D_catNew2 = bestCos2D_catNew2;
	  cos2DBest_catNew2   = BToKEE_fit_cos2D[thisB];
	  bestCos2D_catNew2   = thisB; 
	} else if (whichCatNew==2 && BToKEE_fit_cos2D[thisB]>cos2DSecond_catNew2) {
	  cos2DSecond_catNew2 = BToKEE_fit_cos2D[thisB];
	  secondCos2D_catNew2 = thisB; 
	}

	// ele1+ele2+ele3 sumpt
	float allptsum = Electron_pt[ele1_idx]+Electron_pt[ele2_idx]+ProbeTracks_pt[k_idx];
	if (allptsum>allptsumBest_all) {
	  allptsumSecond_all = allptsumBest_all;
	  secondAllPtSum_all = bestAllPtSum_all;
	  allptsumBest_all   = allptsum;
	  bestAllPtSum_all   = thisB; 
	} else if (allptsum>allptsumSecond_all) {
	  allptsumSecond_all = allptsum;
	  secondAllPtSum_all = thisB; 
	}

	if (whichCat==0 && allptsum>allptsumBest_cat0) {
	  allptsumSecond_cat0 = allptsumBest_cat0;
	  secondAllPtSum_cat0 = bestAllPtSum_cat0;
	  allptsumBest_cat0   = allptsum;
	  bestAllPtSum_cat0   = thisB; 
	} else if (whichCat==0 && allptsum>allptsumSecond_cat0) {
	  allptsumSecond_cat0 = allptsum;
	  secondAllPtSum_cat0 = thisB; 
	}

	if (whichCat==1 && allptsum>allptsumBest_cat1) {
	  allptsumSecond_cat1 = allptsumBest_cat1;
	  secondAllPtSum_cat1 = bestAllPtSum_cat1;
	  allptsumBest_cat1   = allptsum;
	  bestAllPtSum_cat1   = thisB; 
	} else if (whichCat==1 && allptsum>allptsumSecond_cat1) {
	  allptsumSecond_cat1 = allptsum;
	  secondAllPtSum_cat1 = thisB; 
	}

	if (whichCat==2 && allptsum>allptsumBest_cat2) {
	  allptsumSecond_cat2 = allptsumBest_cat2;
	  secondAllPtSum_cat2 = bestAllPtSum_cat2;
	  allptsumBest_cat2   = allptsum;
	  bestAllPtSum_cat2   = thisB; 
	} else if (whichCat==2 && allptsum>allptsumSecond_cat2) {
	  allptsumSecond_cat2 = allptsum;
	  secondAllPtSum_cat2 = thisB; 
	}

	if (whichCatNew==0 && allptsum>allptsumBest_catNew0) {
	  allptsumSecond_catNew0 = allptsumBest_catNew0;
	  secondAllPtSum_catNew0 = bestAllPtSum_catNew0;
	  allptsumBest_catNew0   = allptsum;
	  bestAllPtSum_catNew0   = thisB; 
	} else if (whichCatNew==0 && allptsum>allptsumSecond_catNew0) {
	  allptsumSecond_catNew0 = allptsum;
	  secondAllPtSum_catNew0 = thisB; 
	}

	if (whichCatNew==1 && allptsum>allptsumBest_catNew1) {
	  allptsumSecond_catNew1 = allptsumBest_catNew1;
	  secondAllPtSum_catNew1 = bestAllPtSum_catNew1;
	  allptsumBest_catNew1   = allptsum;
	  bestAllPtSum_catNew1   = thisB; 
	} else if (whichCatNew==1 && allptsum>allptsumSecond_catNew1) {
	  allptsumSecond_catNew1 = allptsum;
	  secondAllPtSum_catNew1 = thisB; 
	}

	if (whichCatNew==2 && allptsum>allptsumBest_catNew2) {
	  allptsumSecond_catNew2 = allptsumBest_catNew2;
	  secondAllPtSum_catNew2 = bestAllPtSum_catNew2;
	  allptsumBest_catNew2   = allptsum;
	  bestAllPtSum_catNew2   = thisB; 
	} else if (whichCatNew==2 && allptsum>allptsumSecond_catNew2) {
	  allptsumSecond_catNew2 = allptsum;
	  secondAllPtSum_catNew2 = thisB; 
	}

	// matched only 
	if (isMcB(thisB)) {
	  if (BToKEE_svprob[thisB]>svprobBest_all_matched) svprobBest_all_matched = BToKEE_svprob[thisB];
	  if (theXySig>xysigBest_all_matched) xysigBest_all_matched = theXySig;
	  if (BToKEE_fit_cos2D[thisB]>cos2DBest_all_matched) cos2DBest_all_matched = BToKEE_fit_cos2D[thisB];
	  if (allptsum>allPtSumBest_all_matched) allPtSumBest_all_matched = allptsum;
	}

      } // Loop over good Bs

      // Position of first matching B: how many combinatorics B there are before the good one
      numberBetterSvProbTriplets=0;
      numberBetterXYsigTriplets=0;
      numberBetterCos2DTriplets=0;
      numberBetterAllPtSumTriplets=0;

      for (u_int iB=0; iB<goodBs.size(); iB++) {        
	int thisB = goodBs[iB];
	int ele1_idx = BToKEE_l1Idx[thisB];
	int ele2_idx = BToKEE_l2Idx[thisB];
	int k_idx    = BToKEE_kIdx[thisB];
	float theXySig = BToKEE_l_xy[thisB]/BToKEE_l_xy_unc[thisB];
	float eleptsum = Electron_pt[ele1_idx] + Electron_pt[ele2_idx];
	float allptsum = eleptsum+ProbeTracks_pt[k_idx];

	if (svprobBest_all_matched>-900) {
	  if (BToKEE_svprob[thisB]>svprobBest_all_matched) numberBetterSvProbTriplets++;
	} else {
	  numberBetterSvProbTriplets=-999;
	}
	if (xysigBest_all_matched>-900) {
	  if (theXySig>xysigBest_all_matched) numberBetterXYsigTriplets++;
	} else {
	  numberBetterXYsigTriplets=-999;
	}
	if (cos2DBest_all_matched>-900) {
	  if (BToKEE_fit_cos2D[thisB]>cos2DBest_all_matched) numberBetterCos2DTriplets++;
	} else {
	  numberBetterCos2DTriplets=-999;
	}
	if (allPtSumBest_all_matched>-900) {
	  if (allptsum>allPtSumBest_all_matched) numberBetterAllPtSumTriplets++;
	} else {
	  numberBetterAllPtSumTriplets=-999;
	}

      } // Loop over good Bs

    
      // Counting, per event
      bestSvProbMatch=0;
      bestSvProbMatch_second=0;
      bestSvProbMatch_notmatching=-999;
      bestSvProbMatchCat0=-999;
      bestSvProbMatchCat1=-999;
      bestSvProbMatchCat2=-999;
      bestSvProbMatchCat0_second=-999;
      bestSvProbMatchCat1_second=-999;
      bestSvProbMatchCat2_second=-999;
      bestSvProbMatchCatNew0=-999;
      bestSvProbMatchCatNew1=-999;
      bestSvProbMatchCatNew2=-999;
      bestSvProbMatchCatNew0_second=-999;
      bestSvProbMatchCatNew1_second=-999;
      bestSvProbMatchCatNew2_second=-999;
      //
      bestXYsigMatch=0;
      bestXYsigMatch_second=0;
      bestXYsigMatch_notmatching=-999;
      bestXYsigMatchCat0=-999;     // PFPF
      bestXYsigMatchCat1=-999;     // PFLP
      bestXYsigMatchCat2=-999;     // LPLP
      bestXYsigMatchCat0_second=-999;  
      bestXYsigMatchCat1_second=-999;  
      bestXYsigMatchCat2_second=-999;  
      bestXYsigMatchCatNew0=-999;     // ele2 PF
      bestXYsigMatchCatNew1=-999;     // ele2 LP, ID>2
      bestXYsigMatchCatNew2=-999;     // ele2 LP, ID<2
      bestXYsigMatchCatNew0_second=-999;  
      bestXYsigMatchCatNew1_second=-999;  
      bestXYsigMatchCatNew2_second=-999;  
      //
      bestCos2DMatch=0;
      bestCos2DMatch_second=0;
      bestCos2DMatchCat0=-999;
      bestCos2DMatchCat1=-999;
      bestCos2DMatchCat2=-999;
      bestCos2DMatchCat0_second=-999;
      bestCos2DMatchCat1_second=-999;
      bestCos2DMatchCat2_second=-999;
      bestCos2DMatchCatNew0=-999;
      bestCos2DMatchCatNew1=-999;
      bestCos2DMatchCatNew2=-999;
      bestCos2DMatchCatNew0_second=-999;
      bestCos2DMatchCatNew1_second=-999;
      bestCos2DMatchCatNew2_second=-999;
      //
      bestAllPtSumMatch=0;
      bestAllPtSumMatch_second=0;
      bestAllPtSumMatchCat0=-999;
      bestAllPtSumMatchCat1=-999;
      bestAllPtSumMatchCat2=-999;
      bestAllPtSumMatchCat0_second=-999;
      bestAllPtSumMatchCat1_second=-999;
      bestAllPtSumMatchCat2_second=-999;
      bestAllPtSumMatchCatNew0=-999;
      bestAllPtSumMatchCatNew1=-999;
      bestAllPtSumMatchCatNew2=-999;
      bestAllPtSumMatchCatNew0_second=-999;
      bestAllPtSumMatchCatNew1_second=-999;
      bestAllPtSumMatchCatNew2_second=-999;
      //
      //
      if (isMcB(bestSvprob_all)) { 
	bestSvProbMatch=1;
	int ele1_idx = BToKEE_l1Idx[bestSvprob_all];
	int ele2_idx = BToKEE_l2Idx[bestSvprob_all];
	int k_idx    = BToKEE_kIdx[bestSvprob_all];
	bestSvProbMatch_ok_ele1pt  = Electron_pt[ele1_idx];
	bestSvProbMatch_ok_ele2pt  = Electron_pt[ele2_idx];
	bestSvProbMatch_ok_kpt     = ProbeTracks_pt[k_idx];
	bestSvProbMatch_ok_ele1eta = Electron_eta[ele1_idx];
	bestSvProbMatch_ok_ele2eta = Electron_eta[ele2_idx];
	bestSvProbMatch_ok_keta    = ProbeTracks_eta[k_idx];
      } else {
	bestSvProbMatch=0;
	int ele1_idx = BToKEE_l1Idx[bestSvprob_all];
	int ele2_idx = BToKEE_l2Idx[bestSvprob_all];
	int k_idx    = BToKEE_kIdx[bestSvprob_all];
	int k_genPartIdx    = ProbeTracks_genPartIdx[k_idx];  
	int ele1_genPartIdx = Electron_genPartIdx[ele1_idx];  
	int ele2_genPartIdx = Electron_genPartIdx[ele2_idx];  
	if (ele1_genPartIdx>-0.5) bestSvProbMatch_causeEle1=0;
	else bestSvProbMatch_causeEle1=1;
	if (ele2_genPartIdx>-0.5) bestSvProbMatch_causeEle2=0;
	else bestSvProbMatch_causeEle2=1;
	if (k_genPartIdx>-0.5) bestSvProbMatch_causeK=0;
	else bestSvProbMatch_causeK=1;
	bestSvProbMatch_notmatching = whichMcB(bestSvprob_all);
	//
	bestSvProbMatch_notok_ele1pt  = Electron_pt[ele1_idx];
	bestSvProbMatch_notok_ele2pt  = Electron_pt[ele2_idx];
	bestSvProbMatch_notok_kpt     = ProbeTracks_pt[k_idx];
	bestSvProbMatch_notok_ele1eta = Electron_eta[ele1_idx];
	bestSvProbMatch_notok_ele2eta = Electron_eta[ele2_idx];
	bestSvProbMatch_notok_keta    = ProbeTracks_eta[k_idx];
      }

      if (isMcB(secondSvprob_all)) bestSvProbMatch_second=1;
      else bestSvProbMatch_second=0;

      if (bestSvprob_cat0>=0) {
	if (isMcB(bestSvprob_cat0)) bestSvProbMatchCat0=1;
	else bestSvProbMatchCat0=0;
	if (isMcB(secondSvprob_cat0)) bestSvProbMatchCat0_second=1;
	else bestSvProbMatchCat0_second=0;
      }

      else if (bestSvprob_cat1>=0) {                         // with else this becomes exclusive: order events in 0-1-2
	if (isMcB(bestSvprob_cat1)) bestSvProbMatchCat1=1;
	else bestSvProbMatchCat1=0;
	if (isMcB(secondSvprob_cat1)) bestSvProbMatchCat1_second=1;
	else bestSvProbMatchCat1_second=0;
      }

      else if (bestSvprob_cat2>=0) {
	if (isMcB(bestSvprob_cat2)) bestSvProbMatchCat2=1;
	else bestSvProbMatchCat2=0;
	if (isMcB(secondSvprob_cat2)) bestSvProbMatchCat2_second=1;
	else bestSvProbMatchCat2_second=0;
      }

      if (bestSvprob_catNew0>=0) {
	if (isMcB(bestSvprob_catNew0)) bestSvProbMatchCatNew0=1;
	else bestSvProbMatchCatNew0=0;
	if (isMcB(secondSvprob_catNew0)) bestSvProbMatchCatNew0_second=1;
	else bestSvProbMatchCatNew0_second=0;
      }

      else if (bestSvprob_catNew1>=0) {                         // with else this becomes exclusive: order events in 0-1-2
	if (isMcB(bestSvprob_catNew1)) bestSvProbMatchCatNew1=1;
	else bestSvProbMatchCatNew1=0;
	if (isMcB(secondSvprob_catNew1)) bestSvProbMatchCatNew1_second=1;
	else bestSvProbMatchCatNew1_second=0;
      }

      else if (bestSvprob_catNew2>=0) {
	if (isMcB(bestSvprob_catNew2)) bestSvProbMatchCatNew2=1;
	else bestSvProbMatchCatNew2=0;
	if (isMcB(secondSvprob_catNew2)) bestSvProbMatchCatNew2_second=1;
	else bestSvProbMatchCatNew2_second=0;
      }

      if (isMcB(bestXYsig_all)) { 
	bestXYsigMatch=1;
	int ele1_idx = BToKEE_l1Idx[bestXYsig_all];
	int ele2_idx = BToKEE_l2Idx[bestXYsig_all];
	int k_idx    = BToKEE_kIdx[bestXYsig_all];
	bestXYsigMatch_ok_ele1pt  = Electron_pt[ele1_idx];
	bestXYsigMatch_ok_ele2pt  = Electron_pt[ele2_idx];
	bestXYsigMatch_ok_kpt     = ProbeTracks_pt[k_idx];
	if (Electron_pt[ele2_idx]<=ProbeTracks_pt[k_idx])
	  bestXYsigMatch_ok_minpt = Electron_pt[ele2_idx];
	else 
	  bestXYsigMatch_ok_minpt = ProbeTracks_pt[k_idx];
	bestXYsigMatch_ok_ele1eta = Electron_eta[ele1_idx];
	bestXYsigMatch_ok_ele2eta = Electron_eta[ele2_idx];
	bestXYsigMatch_ok_keta    = ProbeTracks_eta[k_idx];
	bestXYsigMatch_ok_pfmva1  = Electron_pfmvaId[ele1_idx]; 
	bestXYsigMatch_ok_pfmva2  = Electron_pfmvaId[ele2_idx]; 
	bestXYsigMatch_ok_lptmva1 = Electron_mvaId[ele1_idx]; 
	bestXYsigMatch_ok_lptmva2 = Electron_mvaId[ele2_idx]; 
	bestXYsigMatch_ok_costhetaSK     = cosThetaStarK(bestXYsig_all);
	bestXYsigMatch_ok_costhetaSKCS   = cosThetaStarKCS(bestXYsig_all);
	bestXYsigMatch_ok_costhetaL      = cosThetaL(bestXYsig_all);
	bestXYsigMatch_ok_costhetaSK_gen = cosThetaStarKGen(bestXYsig_all);

      } else {
	bestXYsigMatch=0;
	int ele1_idx = BToKEE_l1Idx[bestXYsig_all];
	int ele2_idx = BToKEE_l2Idx[bestXYsig_all];
	int k_idx    = BToKEE_kIdx[bestXYsig_all];
	int k_genPartIdx    = ProbeTracks_genPartIdx[k_idx];  
	int ele1_genPartIdx = Electron_genPartIdx[ele1_idx];  
	int ele2_genPartIdx = Electron_genPartIdx[ele2_idx];  
	if (ele1_genPartIdx>-0.5) bestXYsigMatch_causeEle1=0;
	else bestXYsigMatch_causeEle1=1;
	if (ele2_genPartIdx>-0.5) bestXYsigMatch_causeEle2=0;
	else bestXYsigMatch_causeEle2=1;
	if (k_genPartIdx>-0.5) bestXYsigMatch_causeK=0;
	else bestXYsigMatch_causeK=1;
	//
	bestXYsigMatch_notok_ele1pt  = Electron_pt[ele1_idx];
	bestXYsigMatch_notok_ele2pt  = Electron_pt[ele2_idx];
	bestXYsigMatch_notok_kpt     = ProbeTracks_pt[k_idx];
	if (Electron_pt[ele2_idx]<=ProbeTracks_pt[k_idx])
	  bestXYsigMatch_notok_minpt = Electron_pt[ele2_idx];
	else 
	  bestXYsigMatch_notok_minpt = ProbeTracks_pt[k_idx];
	bestXYsigMatch_notok_ele1eta = Electron_eta[ele1_idx];
	bestXYsigMatch_notok_ele2eta = Electron_eta[ele2_idx];
	bestXYsigMatch_notok_keta    = ProbeTracks_eta[k_idx];
	bestXYsigMatch_notok_pfmva1  = Electron_pfmvaId[ele1_idx]; 
	bestXYsigMatch_notok_pfmva2  = Electron_pfmvaId[ele2_idx]; 
	bestXYsigMatch_notok_lptmva1 = Electron_mvaId[ele1_idx]; 
	bestXYsigMatch_notok_lptmva2 = Electron_mvaId[ele2_idx]; 	
	bestXYsigMatch_notok_costhetaSK   = cosThetaStarK(bestXYsig_all);
	bestXYsigMatch_notok_costhetaSKCS = cosThetaStarKCS(bestXYsig_all);
	bestXYsigMatch_notok_costhetaL    = cosThetaL(bestXYsig_all);
	bestXYsigMatch_notmatching = whichMcB(bestXYsig_all);
      }

      if (isMcB(secondXYsig_all)) bestXYsigMatch_second=1;
      else bestXYsigMatch_second=0;

      if (bestXYsig_cat0>=0) {
	if (isMcB(bestXYsig_cat0)) bestXYsigMatchCat0=1;
	else bestXYsigMatchCat0=0;
	if (isMcB(secondXYsig_cat0)) bestXYsigMatchCat0_second=1;
	else bestXYsigMatchCat0_second=0;
      }

      else if (bestXYsig_cat1>=0) {                         // with else this becomes exclusive: order events in 0-1-2
	if (isMcB(bestXYsig_cat1)) bestXYsigMatchCat1=1;
	else bestXYsigMatchCat1=0;
	if (isMcB(secondXYsig_cat1)) bestXYsigMatchCat1_second=1;
	else bestXYsigMatchCat1_second=0;
      }

      else if (bestXYsig_cat2>=0) {
	if (isMcB(bestXYsig_cat2)) bestXYsigMatchCat2=1;
	else bestXYsigMatchCat2=0;
	if (isMcB(secondXYsig_cat2)) bestXYsigMatchCat2_second=1;
	else bestXYsigMatchCat2_second=0;
      }

      if (bestXYsig_catNew0>=0) {
	if (isMcB(bestXYsig_catNew0)) bestXYsigMatchCatNew0=1;
	else bestXYsigMatchCatNew0=0;
	if (isMcB(secondXYsig_catNew0)) bestXYsigMatchCatNew0_second=1;
	else bestXYsigMatchCatNew0_second=0;
      }

      else if (bestXYsig_catNew1>=0) {                         // with else this becomes exclusive: order events in 0-1-2
	if (isMcB(bestXYsig_catNew1)) bestXYsigMatchCatNew1=1;
	else bestXYsigMatchCatNew1=0;
	if (isMcB(secondXYsig_catNew1)) bestXYsigMatchCatNew1_second=1;
	else bestXYsigMatchCatNew1_second=0;
      }

      else if (bestXYsig_catNew2>=0) {
	if (isMcB(bestXYsig_catNew2)) bestXYsigMatchCatNew2=1;
	else bestXYsigMatchCatNew2=0;
	if (isMcB(secondXYsig_catNew2)) bestXYsigMatchCatNew2_second=1;
	else bestXYsigMatchCatNew2_second=0;
      }

      if (isMcB(bestCos2D_all)) {
	bestCos2DMatch=1;
	int ele1_idx = BToKEE_l1Idx[bestCos2D_all];
	int ele2_idx = BToKEE_l2Idx[bestCos2D_all];
	int k_idx    = BToKEE_kIdx[bestCos2D_all];
	bestCos2DMatch_ok_ele1pt  = Electron_pt[ele1_idx];
	bestCos2DMatch_ok_ele2pt  = Electron_pt[ele2_idx];
	bestCos2DMatch_ok_kpt     = ProbeTracks_pt[k_idx];
	bestCos2DMatch_ok_ele1eta = Electron_eta[ele1_idx];
	bestCos2DMatch_ok_ele2eta = Electron_eta[ele2_idx];
	bestCos2DMatch_ok_keta    = ProbeTracks_eta[k_idx];
      } else {
	bestCos2DMatch=0;
	int ele1_idx = BToKEE_l1Idx[bestCos2D_all];
	int ele2_idx = BToKEE_l2Idx[bestCos2D_all];
	int k_idx    = BToKEE_kIdx[bestCos2D_all];
	int k_genPartIdx    = ProbeTracks_genPartIdx[k_idx];  
	int ele1_genPartIdx = Electron_genPartIdx[ele1_idx];  
	int ele2_genPartIdx = Electron_genPartIdx[ele2_idx];  
	if (ele1_genPartIdx>-0.5) bestCos2DMatch_causeEle1=0;
	else bestCos2DMatch_causeEle1=1;
	if (ele2_genPartIdx>-0.5) bestCos2DMatch_causeEle2=0;
	else bestCos2DMatch_causeEle2=1;
	if (k_genPartIdx>-0.5) bestCos2DMatch_causeK=0;
	else bestCos2DMatch_causeK=1;
	//
	bestCos2DMatch_notok_ele1pt  = Electron_pt[ele1_idx];
	bestCos2DMatch_notok_ele2pt  = Electron_pt[ele2_idx];
	bestCos2DMatch_notok_kpt     = ProbeTracks_pt[k_idx];
	bestCos2DMatch_notok_ele1eta = Electron_eta[ele1_idx];
	bestCos2DMatch_notok_ele2eta = Electron_eta[ele2_idx];
	bestCos2DMatch_notok_keta    = ProbeTracks_eta[k_idx];	
      }

      if (isMcB(secondCos2D_all)) bestCos2DMatch_second=1;
      else bestCos2DMatch_second=0;

      if (bestCos2D_cat0>=0) {
	if (isMcB(bestCos2D_cat0)) bestCos2DMatchCat0=1;
	else bestCos2DMatchCat0=0;
	if (isMcB(secondCos2D_cat0)) bestCos2DMatchCat0_second=1;
	else bestCos2DMatchCat0_second=0;
      }

      else if (bestCos2D_cat1>=0) {
	if (isMcB(bestCos2D_cat1)) bestCos2DMatchCat1=1;
	else bestCos2DMatchCat1=0;
	if (isMcB(secondCos2D_cat1)) bestCos2DMatchCat1_second=1;
	else bestCos2DMatchCat1_second=0;
      }

      else if (bestCos2D_cat2>=0) {
	if (isMcB(bestCos2D_cat2)) bestCos2DMatchCat2=1;
	else bestCos2DMatchCat2=0;
	if (isMcB(secondCos2D_cat2)) bestCos2DMatchCat2_second=1;
	else bestCos2DMatchCat2_second=0;
      }

      if (bestCos2D_catNew0>=0) {
	if (isMcB(bestCos2D_catNew0)) bestCos2DMatchCatNew0=1;
	else bestCos2DMatchCatNew0=0;
	if (isMcB(secondCos2D_catNew0)) bestCos2DMatchCatNew0_second=1;
	else bestCos2DMatchCatNew0_second=0;
      }

      else if (bestCos2D_catNew1>=0) {
	if (isMcB(bestCos2D_catNew1)) bestCos2DMatchCatNew1=1;
	else bestCos2DMatchCatNew1=0;
	if (isMcB(secondCos2D_catNew1)) bestCos2DMatchCatNew1_second=1;
	else bestCos2DMatchCatNew1_second=0;
      }

      else if (bestCos2D_catNew2>=0) {
	if (isMcB(bestCos2D_catNew2)) bestCos2DMatchCatNew2=1;
	else bestCos2DMatchCatNew2=0;
	if (isMcB(secondCos2D_catNew2)) bestCos2DMatchCatNew2_second=1;
	else bestCos2DMatchCatNew2_second=0;
      }


      if (isMcB(bestAllPtSum_all)) bestAllPtSumMatch=1;
      else {
	bestAllPtSumMatch=0;
	int ele1_idx = BToKEE_l1Idx[bestAllPtSum_all];
	int ele2_idx = BToKEE_l2Idx[bestAllPtSum_all];
	int k_idx    = BToKEE_kIdx[bestAllPtSum_all];
	int k_genPartIdx    = ProbeTracks_genPartIdx[k_idx];  
	int ele1_genPartIdx = Electron_genPartIdx[ele1_idx];  
	int ele2_genPartIdx = Electron_genPartIdx[ele2_idx];  
	if (ele1_genPartIdx>-0.5) bestAllPtSumMatch_causeEle1=0;
	else bestAllPtSumMatch_causeEle1=1;
	if (ele2_genPartIdx>-0.5) bestAllPtSumMatch_causeEle2=0;
	else bestAllPtSumMatch_causeEle2=1;
	if (k_genPartIdx>-0.5) bestAllPtSumMatch_causeK=0;
	else bestAllPtSumMatch_causeK=1;
      }

      if (isMcB(secondAllPtSum_all)) bestAllPtSumMatch_second=1;
      else bestAllPtSumMatch_second=0;

      if (bestAllPtSum_cat0>=0) {
	if (isMcB(bestAllPtSum_cat0)) bestAllPtSumMatchCat0=1;
	else bestAllPtSumMatchCat0=0;
	if (isMcB(secondAllPtSum_cat0)) bestAllPtSumMatchCat0_second=1;
	else bestAllPtSumMatchCat0_second=0;
      }

      else if (bestAllPtSum_cat1>=0) {
	if (isMcB(bestAllPtSum_cat1)) bestAllPtSumMatchCat1=1;
	else bestAllPtSumMatchCat1=0;
	if (isMcB(secondAllPtSum_cat1)) bestAllPtSumMatchCat1_second=1;
	else bestAllPtSumMatchCat1_second=0;
      }

      else if (bestAllPtSum_cat2>=0) {
	if (isMcB(bestAllPtSum_cat2)) bestAllPtSumMatchCat2=1;
	else bestAllPtSumMatchCat2=0;
	if (isMcB(secondAllPtSum_cat2)) bestAllPtSumMatchCat2_second=1;
	else bestAllPtSumMatchCat2_second=0;
      }

      if (bestAllPtSum_catNew0>=0) {
	if (isMcB(bestAllPtSum_catNew0)) bestAllPtSumMatchCatNew0=1;
	else bestAllPtSumMatchCatNew0=0;
	if (isMcB(secondAllPtSum_catNew0)) bestAllPtSumMatchCatNew0_second=1;
	else bestAllPtSumMatchCatNew0_second=0;
      }

      else if (bestAllPtSum_catNew1>=0) {
	if (isMcB(bestAllPtSum_catNew1)) bestAllPtSumMatchCatNew1=1;
	else bestAllPtSumMatchCatNew1=0;
	if (isMcB(secondAllPtSum_catNew1)) bestAllPtSumMatchCatNew1_second=1;
	else bestAllPtSumMatchCatNew1_second=0;
      }

      else if (bestAllPtSum_catNew2>=0) {
	if (isMcB(bestAllPtSum_catNew2)) bestAllPtSumMatchCatNew2=1;
	else bestAllPtSumMatchCatNew2=0;
	if (isMcB(secondAllPtSum_catNew2)) bestAllPtSumMatchCatNew2_second=1;
	else bestAllPtSumMatchCatNew2_second=0;
      }

    } // >=1 1 true good B, >= 1 comb good B

    // Filling the output tree
    outTree_->Fill();


    // --------------------------------------------------------------------------

    // Cleaning all vectors used for the selection
    goodBs.clear();
    goodTrueBs.clear();
    goodCombBs.clear();

    // Debug
    debug_svprob.clear();
    debug_svprob_match.clear();
    debug_pf_svprob.clear();
    debug_pf_svprob_match.clear();

    // Cleaning all vectors used for the output tree, ready for a new entry
    goodTrueB_maxMinDREle.clear();  
    goodTrueB_maxMinDREle_dEta.clear();  
    goodTrueB_maxMinDREle_dPhi.clear();  
    goodTrueB_maxMinDREle_dPtOverPt.clear();
    goodTrueB_dRgen.clear();
    goodTrueB_maxDRTrack.clear();
    //
    goodCombB_maxMinDREle.clear();  
    goodCombB_maxMinDREle_dEta.clear();  
    goodCombB_maxMinDREle_dPhi.clear();  
    goodCombB_maxMinDREle_dPtOverPt.clear();
    goodCombB_maxDRTrack.clear();
    //
    goodTrueB_svProb_notBestMatch.clear();
    goodTrueB_xySig_notBestMatch.clear();
    goodTrueB_cos2D_notBestMatch.clear();
    goodTrueB_ptsum_notBestMatch.clear();
    goodTrueB_kpt_notBestMatch.clear();
    //
    goodCombB_svProb.clear();
    goodCombB_xySig.clear();
    goodCombB_cos2D.clear();
    goodCombB_ptsum.clear();
    goodCombB_kpt.clear();

    goodCombB_keta.clear();
    goodCombB_ele1pt.clear();
    goodCombB_ele2pt.clear();
    goodCombB_minpt.clear();
    goodCombB_ele1eta.clear();
    goodCombB_ele2eta.clear();
    goodCombB_ele1pfmva.clear();
    goodCombB_ele2pfmva.clear();
    goodCombB_ele1lptmva.clear();
    goodCombB_ele2lptmva.clear();
    goodCombB_causeEle1.clear();
    goodCombB_causeEle2.clear();
    goodCombB_causeK.clear();
    goodCombB_notmatching.clear();
    goodCombB_maxDrRecoGen.clear();
    goodCombB_minDrRecoGen.clear();
    goodCombB_drRecoGenK.clear();
    goodCombB_maxDrRecoGenFromB.clear();
    goodCombB_minDrRecoGenFromB.clear();
    goodCombB_drRecoGenFromBK.clear();
  }

  // Summary
  cout << "myAllGoodPreselBs = " << myAllGoodPreselBs << " over myAllGoodBs = " << myAllGoodBs << endl;
  cout << endl;
  cout << "Keeping all triplets, I've " << allTriplets << " triplets to analyze" << endl;
  for (int nn=0; nn<20; nn++) cout << "Keeping " << nn+1 << " ==> " << nTriplets[nn] << endl;
}



// for B -> K J/Psi -> Kee (resonant)
/*
bool TripletSelection::isMcB( int theB ) {
  
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
  // 443 = J/Psi; 521 = B+
  bool okMatch = (ele1_genPartIdx>-0.5 && ele2_genPartIdx>-0.5 && k_genPartIdx>-0.5);
  bool RK_res1 = abs(ele1_genMotherPdgId)==443 && abs(k_genMotherPdgId)==521;
  bool RK_res2 = (ele1_genMotherPdgId==ele2_genMotherPdgId) && (k_genMotherPdgId==ele1_genGMotherPdgId) && (k_genMotherPdgId==ele2_genGMotherPdgId);
  bool RK_res = okMatch && RK_res1 && RK_res2;

  return RK_res;
}
*/

// for B -> Kee (non resonant)
bool TripletSelection::isMcB( int theB ) {
  
  // taking index
  int ele1_idx = BToKEE_l1Idx[theB];
  int ele2_idx = BToKEE_l2Idx[theB];
  int k_idx    = BToKEE_kIdx[theB];

  // Gen tree
  int k_genPartIdx     = ProbeTracks_genPartIdx[k_idx];  
  int k_genMotherIdx   = GenPart_genPartIdxMother[k_genPartIdx];
  int k_genMotherPdgId = GenPart_pdgId[k_genMotherIdx];

  int ele1_genPartIdx     = Electron_genPartIdx[ele1_idx];  
  int ele1_genMotherIdx   = GenPart_genPartIdxMother[ele1_genPartIdx];
  int ele1_genMotherPdgId = GenPart_pdgId[ele1_genMotherIdx];

  int ele2_genPartIdx     = Electron_genPartIdx[ele2_idx];  
  int ele2_genMotherIdx   = GenPart_genPartIdxMother[ele2_genPartIdx];
  int ele2_genMotherPdgId = GenPart_pdgId[ele2_genMotherIdx];

  // B -> K J/psi(ll) at gen level
  // 521 = B+
  bool okMatch  = (ele1_genPartIdx>-0.5 && ele2_genPartIdx>-0.5 && k_genPartIdx>-0.5);
  bool RK_nres1 = abs(ele1_genMotherPdgId)==521 && abs(ele2_genMotherPdgId)==521 && abs(k_genMotherPdgId)==521;
  bool RK_nres2 = (ele1_genMotherPdgId==ele2_genMotherPdgId) && (k_genMotherPdgId==ele1_genMotherPdgId) && (k_genMotherPdgId==ele2_genMotherPdgId);
  bool RK_nres  = okMatch && RK_nres1 && RK_nres2;

  return RK_nres;
}

// for B -> Kee (non resonant)
int TripletSelection::whichMcB( int theB ) {
  
  // taking index
  int ele1_idx = BToKEE_l1Idx[theB];
  int ele2_idx = BToKEE_l2Idx[theB];
  int k_idx    = BToKEE_kIdx[theB];

  // Gen tree
  int k_genPartIdx     = ProbeTracks_genPartIdx[k_idx];  
  int k_genMotherIdx   = GenPart_genPartIdxMother[k_genPartIdx];
  int k_genMotherPdgId = GenPart_pdgId[k_genMotherIdx];

  int ele1_genPartIdx     = Electron_genPartIdx[ele1_idx];  
  int ele1_genMotherIdx   = GenPart_genPartIdxMother[ele1_genPartIdx];
  int ele1_genMotherPdgId = GenPart_pdgId[ele1_genMotherIdx];

  int ele2_genPartIdx     = Electron_genPartIdx[ele2_idx];  
  int ele2_genMotherIdx   = GenPart_genPartIdxMother[ele2_genPartIdx];
  int ele2_genMotherPdgId = GenPart_pdgId[ele2_genMotherIdx];

  // B -> K ee at gen level
  // 521 = B+
  bool bad1 = false;
  bool bad2 = false;
  bool badK = false;
  if (ele1_genPartIdx<=-0.5) bad1 = true;
  if (ele2_genPartIdx<=-0.5) bad2 = true;
  if (k_genPartIdx<=-0.5)    badK = true;
  //
  if (abs(ele1_genMotherPdgId)!=521) bad1 = true;
  if (abs(ele2_genMotherPdgId)!=521) bad2 = true;
  if (abs(k_genMotherPdgId)!=521)    badK = true;
  //
  //if (ele1_genMotherPdgId!=ele2_genMotherPdgId && ele1_genMotherPdgId!=k_genMotherPdgId) bad1 = true;
  //if (ele2_genMotherPdgId!=ele1_genMotherPdgId && ele2_genMotherPdgId!=k_genMotherPdgId) bad2 = true;
  //if (k_genMotherPdgId!=ele1_genMotherPdgId && k_genMotherPdgId!=ele2_genMotherPdgId)    badK = true;

  if (bad1 && !bad2 && !badK) return 1;   // only ele1 bad
  if (bad2 && !bad1 && !badK) return 2;   // only ele2 bad
  if (badK && !bad1 && !bad2) return 3;   // only K bad
  if (bad1 && bad2 && !badK)  return 4;   // ele1 and ele2 bad; K good
  if (bad1 && badK && !bad2)  return 5;   // ele1 and K bad; ele2 good
  if (bad2 && badK && !bad1)  return 6;   // ele2 and K bad; ele1 good
  if (bad1 && bad2 && badK)   return 7;   // all bad
  else return -1;
}

float TripletSelection::dRgen( int theB ) {
  
  // taking index
  int ele1_idx = BToKEE_l1Idx[theB];
  int ele2_idx = BToKEE_l2Idx[theB];

  // Gen tree
  int ele1_genPartIdx = Electron_genPartIdx[ele1_idx];  
  int ele2_genPartIdx = Electron_genPartIdx[ele2_idx];  

  float ele1_genEta = GenPart_eta[ele1_genPartIdx];
  float ele2_genEta = GenPart_eta[ele2_genPartIdx];
  float ele1_genPhi = GenPart_phi[ele1_genPartIdx];
  float ele2_genPhi = GenPart_phi[ele2_genPartIdx];
  float ele1_genPt  = GenPart_pt[ele1_genPartIdx];
  float ele2_genPt  = GenPart_pt[ele2_genPartIdx];

  TVector3 ele1V(0.,0.,0.);
  TVector3 ele2V(0.,0.,0.);
  ele1V.SetPtEtaPhi(ele1_genPt, ele1_genEta, ele1_genPhi);
  ele2V.SetPtEtaPhi(ele2_genPt, ele2_genEta, ele2_genPhi);

  float deltaR = ele1V.DeltaR(ele2V);
  
  return deltaR;
}

float TripletSelection::dRRecoGenEle( int theRecoEle ) {
  
  // Reco ele
  float reco_eta = Electron_eta[theRecoEle];
  float reco_phi = Electron_phi[theRecoEle];
  float reco_pt  = Electron_pt[theRecoEle];
  TVector3 reco_ele(0.,0.,0.);
  reco_ele.SetPtEtaPhi(reco_pt, reco_eta, reco_phi);

  // Gen ele
  int theGenEle = Electron_genPartIdx[theRecoEle];  
  float gen_eta = GenPart_eta[theGenEle];
  float gen_phi = GenPart_phi[theGenEle];
  float gen_pt  = GenPart_pt[theGenEle];
  TVector3 gen_ele(0.,0.,0.);
  gen_ele.SetPtEtaPhi(gen_pt, gen_eta, gen_phi);

  float deltaR = reco_ele.DeltaR(gen_ele);
  return deltaR;
}

float TripletSelection::dRRecoGenK( int theRecoK ) {
  
  // Reco K
  float reco_eta = ProbeTracks_eta[theRecoK];
  float reco_phi = ProbeTracks_phi[theRecoK];
  float reco_pt  = ProbeTracks_pt[theRecoK];
  TVector3 reco_k(0.,0.,0.);
  reco_k.SetPtEtaPhi(reco_pt, reco_eta, reco_phi);

  // Gen K
  int theGenK = ProbeTracks_genPartIdx[theRecoK];  
  float gen_eta = GenPart_eta[theGenK];
  float gen_phi = GenPart_phi[theGenK];
  float gen_pt  = GenPart_pt[theGenK];
  TVector3 gen_k(0.,0.,0.);
  gen_k.SetPtEtaPhi(gen_pt, gen_eta, gen_phi);

  float deltaR = reco_k.DeltaR(gen_k);
  return deltaR;
}

// Reco level: cos(theta*) = cos angle between B flight direction (defined from B momentum) and K direction in B rest frame
float TripletSelection::cosThetaStarK( int theB ) {

  int ele1_idx = BToKEE_l1Idx[theB];
  int ele2_idx = BToKEE_l2Idx[theB];
  int k_idx    = BToKEE_kIdx[theB];

  TLorentzVector ele1V3(0.,0.,0.,0.);
  TLorentzVector ele2V3(0.,0.,0.,0.);
  TLorentzVector kV3(0.,0.,0.,0.);
  ele1V3.SetPtEtaPhiM(Electron_pt[ele1_idx], Electron_eta[ele1_idx], Electron_phi[ele1_idx], 0.000511);
  ele2V3.SetPtEtaPhiM(Electron_pt[ele2_idx], Electron_eta[ele2_idx], Electron_phi[ele2_idx], 0.000511);
  kV3.SetPtEtaPhiM(ProbeTracks_pt[k_idx],ProbeTracks_eta[k_idx],ProbeTracks_phi[k_idx], 0.000494);
  
  TLorentzVector B = ele1V3 + ele2V3 + kV3;
  TLorentzVector K_Bstar(kV3);
  K_Bstar.Boost(-B.BoostVector());

  TVector3 K_Bstar_perp(0.,0.,0.); 
  K_Bstar_perp.SetPtEtaPhi(K_Bstar.Pt(),K_Bstar.Eta(),K_Bstar.Phi());
  K_Bstar_perp.SetZ(0);

  TVector3 B_perp(0.,0.,0.);
  B_perp.SetPtEtaPhi(B.Pt(),B.Eta(),B.Phi());
  B_perp.SetZ(0);

  double den = (B_perp.Mag() * K_Bstar_perp.Mag());
  if (den!= 0.) return B_perp.Dot(K_Bstar_perp)/den;
  else return -2.;
}

// Gen level: cos(theta*) = cos angle between B flight direction (defined from B momentum) and K direction in B rest frame
float TripletSelection::cosThetaStarKGen( int theB ) {

  int ele1_idx = BToKEE_l1Idx[theB];
  int ele2_idx = BToKEE_l2Idx[theB];
  int k_idx    = BToKEE_kIdx[theB];

  int ele1_genPartIdx = Electron_genPartIdx[ele1_idx];  
  int ele2_genPartIdx = Electron_genPartIdx[ele2_idx];  
  int k_genPartIdx    = ProbeTracks_genPartIdx[k_idx];  
  
  TLorentzVector ele1V3(0.,0.,0.,0.);
  TLorentzVector ele2V3(0.,0.,0.,0.);
  TLorentzVector kV3(0.,0.,0.,0.);
  ele1V3.SetPtEtaPhiM(GenPart_pt[ele1_genPartIdx], GenPart_eta[ele1_genPartIdx], GenPart_phi[ele1_genPartIdx], 0.000511);
  ele2V3.SetPtEtaPhiM(GenPart_pt[ele2_genPartIdx], GenPart_eta[ele2_genPartIdx], GenPart_phi[ele2_genPartIdx], 0.000511);
  kV3.SetPtEtaPhiM(GenPart_pt[k_genPartIdx], GenPart_eta[k_genPartIdx], GenPart_phi[k_genPartIdx], 0.000494);
  
  TLorentzVector B = ele1V3 + ele2V3 + kV3;
  TLorentzVector K_Bstar(kV3);
  K_Bstar.Boost(-B.BoostVector());

  TVector3 K_Bstar_perp(0.,0.,0.); 
  K_Bstar_perp.SetPtEtaPhi(K_Bstar.Pt(),K_Bstar.Eta(),K_Bstar.Phi());
  K_Bstar_perp.SetZ(0);

  TVector3 B_perp(0.,0.,0.);
  B_perp.SetPtEtaPhi(B.Pt(),B.Eta(),B.Phi());
  B_perp.SetZ(0);

  double den = (B_perp.Mag() * K_Bstar_perp.Mag());
  if (den!= 0.) return B_perp.Dot(K_Bstar_perp)/den;
  else return -2.;
}

// Reco level: cos(thetaL) = cos angle between ele1 momentum and the direction opposite to the B momentum, in e+e- rest frame
float TripletSelection::cosThetaL( int theB ) {

  int ele1_idx = BToKEE_l1Idx[theB];
  int ele2_idx = BToKEE_l2Idx[theB];
  int k_idx    = BToKEE_kIdx[theB];

  TLorentzVector ele1V3(0.,0.,0.,0.);
  TLorentzVector ele2V3(0.,0.,0.,0.);
  TLorentzVector kV3(0.,0.,0.,0.);
  ele1V3.SetPtEtaPhiM(Electron_pt[ele1_idx], Electron_eta[ele1_idx], Electron_phi[ele1_idx], 0.000511);
  ele2V3.SetPtEtaPhiM(Electron_pt[ele2_idx], Electron_eta[ele2_idx], Electron_phi[ele2_idx], 0.000511);
  kV3.SetPtEtaPhiM(ProbeTracks_pt[k_idx],ProbeTracks_eta[k_idx],ProbeTracks_phi[k_idx], 0.000494);
  
  TLorentzVector ll = ele1V3 + ele2V3;

  TLorentzVector ele1_LLstar(ele1V3);
  ele1_LLstar.Boost(-ll.BoostVector());

  TLorentzVector minusB = -(ele1V3 + ele2V3 + kV3);
  TLorentzVector minusB_LLstar(minusB);
  minusB_LLstar.Boost(-ll.BoostVector());

  TVector3 ele1_LLstar_perp(0.,0.,0.); 
  ele1_LLstar_perp.SetPtEtaPhi(ele1_LLstar.Pt(),ele1_LLstar.Eta(),ele1_LLstar.Phi());
  ele1_LLstar_perp.SetZ(0);

  TVector3 minusB_LLstar_perp(0.,0.,0.); 
  minusB_LLstar_perp.SetPtEtaPhi(minusB_LLstar.Pt(),minusB_LLstar.Eta(),minusB_LLstar.Phi());
  minusB_LLstar_perp.SetZ(0);

  double den = (ele1_LLstar_perp.Mag() * minusB_LLstar_perp.Mag());
  if (den!= 0.) return ele1_LLstar_perp.Dot(minusB_LLstar_perp)/den;
  else return -2.;
}

// Reco level: cos(theta*) in Collins Sopper frame = cos angle between K and the line that bisects the acute angle between the two
// colliding protons in B rest frame
float TripletSelection::cosThetaStarKCS( int theB ) {

  int ele1_idx = BToKEE_l1Idx[theB];
  int ele2_idx = BToKEE_l2Idx[theB];
  int k_idx    = BToKEE_kIdx[theB];

  TLorentzVector ele1V3(0.,0.,0.,0.);
  TLorentzVector ele2V3(0.,0.,0.,0.);
  TLorentzVector kV3(0.,0.,0.,0.);
  ele1V3.SetPtEtaPhiM(Electron_pt[ele1_idx], Electron_eta[ele1_idx], Electron_phi[ele1_idx], 0.000511);
  ele2V3.SetPtEtaPhiM(Electron_pt[ele2_idx], Electron_eta[ele2_idx], Electron_phi[ele2_idx], 0.000511);
  kV3.SetPtEtaPhiM(ProbeTracks_pt[k_idx],ProbeTracks_eta[k_idx],ProbeTracks_phi[k_idx], 0.000494);
  
  TLorentzVector B = ele1V3 + ele2V3 + kV3;
  TLorentzVector K_Bstar(kV3);
  K_Bstar.Boost(-B.BoostVector());

  return K_Bstar.CosTheta();
}

void TripletSelection::PrepareOutputs(std::string filename) 
{
  _datasetName=filename;

  std::string outname = _datasetName+".root";    
  cout << "output: " << outname << endl;
  outFile_ = new TFile(outname.c_str(),"RECREATE");

  bookOutputTree();
  bookOutputHistos();
};


void TripletSelection::bookOutputTree() 
{
  outTree_ = new TTree("TaPtree", "TaPtree");
  
  cout << "Booking tree" << endl;

  outTree_->Branch("debug_svprob", "std::vector<float>", &debug_svprob);
  outTree_->Branch("debug_svprob_match", "std::vector<int>", &debug_svprob_match);
  outTree_->Branch("debug_pf_svprob", "std::vector<float>", &debug_pf_svprob);
  outTree_->Branch("debug_pf_svprob_match", "std::vector<int>", &debug_pf_svprob_match);

  outTree_->Branch("iHLT_Mu12_IP6", &iHLT_Mu12_IP6, "iHLT_Mu12_IP6/I");
  outTree_->Branch("iHLT_Mu9_IP6",  &iHLT_Mu9_IP6,  "iHLT_Mu9_IP6/I");

  outTree_->Branch("theEvent", &theEvent, "theEvent/I");
  
  outTree_->Branch("rho", &rho, "rho/F");    

  outTree_->Branch("rightMcTruth",     &rightMcTruth,     "rightMcTruth/I");   
  outTree_->Branch("goodBSize",        &goodBSize,        "goodBSize/I");   
  outTree_->Branch("goodTrueBSize",    &goodTrueBSize,    "goodTrueBSize/I");   
  outTree_->Branch("goodCombBSize",    &goodCombBSize,    "goodCombBSize/I");   

  outTree_->Branch("handMadeB",        &handMadeB,        "handMadeB/I");   
  outTree_->Branch("handMadeBmass",    &handMadeBmass,    "handMadeBmass/F");   

  outTree_->Branch("goodTrueB_maxMinDREle",           "std::vector<float>",   &goodTrueB_maxMinDREle); 
  outTree_->Branch("goodTrueB_maxMinDREle_dEta",      "std::vector<float>",   &goodTrueB_maxMinDREle_dEta); 
  outTree_->Branch("goodTrueB_maxMinDREle_dPhi",      "std::vector<float>",   &goodTrueB_maxMinDREle_dPhi); 
  outTree_->Branch("goodTrueB_maxMinDREle_dPtOverPt", "std::vector<float>",   &goodTrueB_maxMinDREle_dPtOverPt); 
  outTree_->Branch("goodTrueB_dRgen",                 "std::vector<float>",   &goodTrueB_dRgen); 
  outTree_->Branch("goodTrueB_maxDRTrack",            "std::vector<float>",   &goodTrueB_maxDRTrack);

  outTree_->Branch("goodCombB_maxMinDREle",           "std::vector<float>",   &goodCombB_maxMinDREle); 
  outTree_->Branch("goodCombB_maxMinDREle_dEta",      "std::vector<float>",   &goodCombB_maxMinDREle_dEta); 
  outTree_->Branch("goodCombB_maxMinDREle_dPhi",      "std::vector<float>",   &goodCombB_maxMinDREle_dPhi); 
  outTree_->Branch("goodCombB_maxMinDREle_dPtOverPt", "std::vector<float>",   &goodCombB_maxMinDREle_dPtOverPt); 
  outTree_->Branch("goodCombB_maxDRTrack",            "std::vector<float>",   &goodCombB_maxDRTrack);

  outTree_->Branch("goodTrueBs_SvProbMatch", &goodTrueBs_SvProbMatch, "goodTrueBs_SvProbMatch/I");
  outTree_->Branch("goodTrueBs_xySigMatch",  &goodTrueBs_xySigMatch,  "goodTrueBs_xySigMatch/I");
  outTree_->Branch("goodTrueBs_cos2DMatch",  &goodTrueBs_cos2DMatch,  "goodTrueBs_cos2DMatch/I");
  outTree_->Branch("goodTrueBs_ptSumMatch",  &goodTrueBs_ptSumMatch,  "goodTrueBs_ptSumMatch/I");
  outTree_->Branch("goodTrueBs_kptMatch",    &goodTrueBs_kptMatch,    "goodTrueBs_kptMatch/I");

  outTree_->Branch("goodTrueBs_bestSvProb_dRmax", &goodTrueBs_bestSvProb_dRmax, "goodTrueBs_bestSvProb_dRmax/F");
  outTree_->Branch("goodTrueBs_bestXYSig_dRmax",  &goodTrueBs_bestXYSig_dRmax,  "goodTrueBs_bestXYSig_dRmax/F");
  outTree_->Branch("goodTrueBs_bestCos2d_dRmax",  &goodTrueBs_bestCos2d_dRmax,  "goodTrueBs_bestCos2d_dRmax/F");
  outTree_->Branch("goodTrueBs_bestPtsum_dRmax",  &goodTrueBs_bestPtsum_dRmax,  "goodTrueBs_bestPtsum_dRmax/F");
  outTree_->Branch("goodTrueBs_bestKpt_dRmax",    &goodTrueBs_bestKpt_dRmax,    "goodTrueBs_bestKpt_dRmax/F");
  outTree_->Branch("goodTrueBs_bestSvProb_dRmin", &goodTrueBs_bestSvProb_dRmin, "goodTrueBs_bestSvProb_dRmin/F");
  outTree_->Branch("goodTrueBs_bestXYSig_dRmin",  &goodTrueBs_bestXYSig_dRmin,  "goodTrueBs_bestXYSig_dRmin/F");
  outTree_->Branch("goodTrueBs_bestCos2d_dRmin",  &goodTrueBs_bestCos2d_dRmin,  "goodTrueBs_bestCos2d_dRmin/F");
  outTree_->Branch("goodTrueBs_bestPtsum_dRmin",  &goodTrueBs_bestPtsum_dRmin,  "goodTrueBs_bestPtsum_dRmin/F");
  outTree_->Branch("goodTrueBs_bestKpt_dRmin",    &goodTrueBs_bestKpt_dRmin,    "goodTrueBs_bestKpt_dRmin/F");

  outTree_->Branch("bestMatch_Bmass",        &bestMatch_Bmass,        "bestMatch_Bmass/F");
  outTree_->Branch("bestMatch_SvProb",       &bestMatch_SvProb,       "bestMatch_SvProb/F");
  outTree_->Branch("bestMatch_XYSig",        &bestMatch_XYSig,        "bestMatch_XYSig/F");
  outTree_->Branch("bestMatch_Cos2D",        &bestMatch_Cos2D,        "bestMatch_Cos2D/F");
  outTree_->Branch("bestMatch_PtSum",        &bestMatch_PtSum,        "bestMatch_PtSum/F");
  outTree_->Branch("bestMatch_KPt",          &bestMatch_KPt,          "bestMatch_KPt/F");
  outTree_->Branch("bestMatch_KEta",         &bestMatch_KEta,         "bestMatch_KEta/F");
  outTree_->Branch("bestMatch_Ele1Pt",       &bestMatch_Ele1Pt,       "bestMatch_Ele1Pt/F");
  outTree_->Branch("bestMatch_Ele2Pt",       &bestMatch_Ele2Pt,       "bestMatch_Ele2Pt/F");
  outTree_->Branch("bestMatch_MinPt",        &bestMatch_MinPt,        "bestMatch_MinPt/F");
  outTree_->Branch("bestMatch_Ele1Eta",      &bestMatch_Ele1Eta,      "bestMatch_Ele1Eta/F");
  outTree_->Branch("bestMatch_Ele2Eta",      &bestMatch_Ele2Eta,      "bestMatch_Ele2Eta/F");
  outTree_->Branch("bestMatch_Ele1pfmva",    &bestMatch_Ele1pfmva,    "bestMatch_Ele1pfmva/F");
  outTree_->Branch("bestMatch_Ele2pfmva",    &bestMatch_Ele2pfmva,    "bestMatch_Ele2pfmva/F");
  outTree_->Branch("bestMatch_Ele1lptmva",   &bestMatch_Ele1lptmva,   "bestMatch_Ele1lptmva/F");
  outTree_->Branch("bestMatch_Ele2lptmva",   &bestMatch_Ele2lptmva,   "bestMatch_Ele2lptmva/F");
  outTree_->Branch("bestMatch_maxDrRecoGen", &bestMatch_maxDrRecoGen, "bestMatch_maxDrRecoGen/F");
  outTree_->Branch("bestMatch_minDrRecoGen", &bestMatch_minDrRecoGen, "bestMatch_minDrRecoGen/F");
  outTree_->Branch("bestMatch_drRecoGenK",   &bestMatch_drRecoGenK,   "bestMatch_drRecoGenK/F");
  outTree_->Branch("bestMatch_maxDrRecoGenFromB", &bestMatch_maxDrRecoGenFromB, "bestMatch_maxDrRecoGenFromB/F");
  outTree_->Branch("bestMatch_minDrRecoGenFromB", &bestMatch_minDrRecoGenFromB, "bestMatch_minDrRecoGenFromB/F");
  outTree_->Branch("bestMatch_drRecoGenFromBK",   &bestMatch_drRecoGenFromBK,   "bestMatch_drRecoGenFromBK/F");

  outTree_->Branch("goodTrueB_svProb_notBestMatch",       "std::vector<float>", &goodTrueB_svProb_notBestMatch);
  outTree_->Branch("goodTrueB_xySig_notBestMatch",        "std::vector<float>", &goodTrueB_xySig_notBestMatch);
  outTree_->Branch("goodTrueB_cos2D_notBestMatch",        "std::vector<float>", &goodTrueB_cos2D_notBestMatch);
  outTree_->Branch("goodTrueB_ptsum_notBestMatch",        "std::vector<float>", &goodTrueB_ptsum_notBestMatch);
  outTree_->Branch("goodTrueB_kpt_notBestMatch",          "std::vector<float>", &goodTrueB_kpt_notBestMatch);

  outTree_->Branch("goodCombB_svProb",       "std::vector<float>", &goodCombB_svProb);
  outTree_->Branch("goodCombB_xySig"  ,      "std::vector<float>", &goodCombB_xySig);
  outTree_->Branch("goodCombB_cos2D",        "std::vector<float>", &goodCombB_cos2D);
  outTree_->Branch("goodCombB_ptsum",        "std::vector<float>", &goodCombB_ptsum);
  outTree_->Branch("goodCombB_kpt",          "std::vector<float>", &goodCombB_kpt);
  outTree_->Branch("goodCombB_keta",         "std::vector<float>", &goodCombB_keta);
  outTree_->Branch("goodCombB_ele1pt",       "std::vector<float>", &goodCombB_ele1pt);
  outTree_->Branch("goodCombB_ele2pt",       "std::vector<float>", &goodCombB_ele2pt);
  outTree_->Branch("goodCombB_minpt",        "std::vector<float>", &goodCombB_minpt);
  outTree_->Branch("goodCombB_ele1eta",      "std::vector<float>", &goodCombB_ele1eta);
  outTree_->Branch("goodCombB_ele2eta",      "std::vector<float>", &goodCombB_ele2eta);
  outTree_->Branch("goodCombB_ele1pfmva",    "std::vector<float>", &goodCombB_ele1pfmva);
  outTree_->Branch("goodCombB_ele2pfmva",    "std::vector<float>", &goodCombB_ele2pfmva);
  outTree_->Branch("goodCombB_ele1lptmva",   "std::vector<float>", &goodCombB_ele1lptmva);
  outTree_->Branch("goodCombB_ele2lptmva",   "std::vector<float>", &goodCombB_ele2lptmva);
  outTree_->Branch("goodCombB_causeEle1",    "std::vector<float>", &goodCombB_causeEle1);
  outTree_->Branch("goodCombB_causeEle2",    "std::vector<float>", &goodCombB_causeEle2);
  outTree_->Branch("goodCombB_causeK",       "std::vector<float>", &goodCombB_causeK);
  outTree_->Branch("goodCombB_notmatching",  "std::vector<int>",   &goodCombB_notmatching);
  outTree_->Branch("goodCombB_maxDrRecoGen", "std::vector<float>", &goodCombB_maxDrRecoGen);
  outTree_->Branch("goodCombB_minDrRecoGen", "std::vector<float>", &goodCombB_minDrRecoGen);
  outTree_->Branch("goodCombB_drRecoGenK",   "std::vector<float>", &goodCombB_drRecoGenK);
  outTree_->Branch("goodCombB_maxDrRecoGenFromB", "std::vector<float>", &goodCombB_maxDrRecoGenFromB);
  outTree_->Branch("goodCombB_minDrRecoGenFromB", "std::vector<float>", &goodCombB_minDrRecoGenFromB);
  outTree_->Branch("goodCombB_drRecoGenFromBK",   "std::vector<float>", &goodCombB_drRecoGenFromBK);

  outTree_->Branch("bestSvProbMatch",               &bestSvProbMatch,               "bestSvProbMatch/I");
  outTree_->Branch("bestSvProbMatch_second",        &bestSvProbMatch_second,        "bestSvProbMatch_second/I");
  outTree_->Branch("bestSvProbMatchCat0",           &bestSvProbMatchCat0,           "bestSvProbMatchCat0/I");
  outTree_->Branch("bestSvProbMatchCat1",           &bestSvProbMatchCat1,           "bestSvProbMatchCat1/I");
  outTree_->Branch("bestSvProbMatchCat2",           &bestSvProbMatchCat2,           "bestSvProbMatchCat2/I");
  outTree_->Branch("bestSvProbMatchCat0_second",    &bestSvProbMatchCat0_second,    "bestSvProbMatchCat0_second/I");
  outTree_->Branch("bestSvProbMatchCat1_second",    &bestSvProbMatchCat1_second,    "bestSvProbMatchCat1_second/I");
  outTree_->Branch("bestSvProbMatchCat2_second",    &bestSvProbMatchCat2_second,    "bestSvProbMatchCat2_second/I");
  outTree_->Branch("bestSvProbMatchCatNew0",        &bestSvProbMatchCatNew0,        "bestSvProbMatchCatNew0/I");
  outTree_->Branch("bestSvProbMatchCatNew1",        &bestSvProbMatchCatNew1,        "bestSvProbMatchCatNew1/I");
  outTree_->Branch("bestSvProbMatchCatNew2",        &bestSvProbMatchCatNew2,        "bestSvProbMatchCatNew2/I");
  outTree_->Branch("bestSvProbMatchCatNew0_second", &bestSvProbMatchCatNew0_second, "bestSvProbMatchCatNew0_second/I");
  outTree_->Branch("bestSvProbMatchCatNew1_second", &bestSvProbMatchCatNew1_second, "bestSvProbMatchCatNew1_second/I");
  outTree_->Branch("bestSvProbMatchCatNew2_second", &bestSvProbMatchCatNew2_second, "bestSvProbMatchCatNew2_second/I");
  outTree_->Branch("bestSvProbMatch_causeEle1",     &bestSvProbMatch_causeEle1,     "bestSvProbMatch_causeEle1/I");
  outTree_->Branch("bestSvProbMatch_causeEle2",     &bestSvProbMatch_causeEle2,     "bestSvProbMatch_causeEle2/I");
  outTree_->Branch("bestSvProbMatch_notmatching",   &bestSvProbMatch_notmatching,   "bestSvProbMatch_notmatching/I");
  outTree_->Branch("bestSvProbMatch_causeK",        &bestSvProbMatch_causeK,        "bestSvProbMatch_causeK/I");
  outTree_->Branch("bestSvProbMatch_notok_ele1pt",  &bestSvProbMatch_notok_ele1pt,  "bestSvProbMatch_notok_ele1pt/F");
  outTree_->Branch("bestSvProbMatch_notok_ele2pt",  &bestSvProbMatch_notok_ele2pt,  "bestSvProbMatch_notok_ele2pt/F");
  outTree_->Branch("bestSvProbMatch_notok_kpt",     &bestSvProbMatch_notok_kpt,     "bestSvProbMatch_notok_kpt/F");
  outTree_->Branch("bestSvProbMatch_notok_ele1eta", &bestSvProbMatch_notok_ele1eta, "bestSvProbMatch_notok_ele1eta/F");
  outTree_->Branch("bestSvProbMatch_notok_ele2eta", &bestSvProbMatch_notok_ele2eta, "bestSvProbMatch_notok_ele2eta/F");
  outTree_->Branch("bestSvProbMatch_notok_keta",    &bestSvProbMatch_notok_keta,    "bestSvProbMatch_notok_keta/F");
  outTree_->Branch("bestSvProbMatch_ok_ele1pt",     &bestSvProbMatch_ok_ele1pt,     "bestSvProbMatch_ok_ele1pt/F");
  outTree_->Branch("bestSvProbMatch_ok_ele2pt",     &bestSvProbMatch_ok_ele2pt,     "bestSvProbMatch_ok_ele2pt/F");
  outTree_->Branch("bestSvProbMatch_ok_kpt",        &bestSvProbMatch_ok_kpt,        "bestSvProbMatch_ok_kpt/F");
  outTree_->Branch("bestSvProbMatch_ok_ele1eta",    &bestSvProbMatch_ok_ele1eta,    "bestSvProbMatch_ok_ele1eta/F");
  outTree_->Branch("bestSvProbMatch_ok_ele2eta",    &bestSvProbMatch_ok_ele2eta,    "bestSvProbMatch_ok_ele2eta/F");
  outTree_->Branch("bestSvProbMatch_ok_keta",       &bestSvProbMatch_ok_keta,       "bestSvProbMatch_ok_keta/F");

  outTree_->Branch("bestXYsigMatch",                &bestXYsigMatch,                "bestXYsigMatch/I");
  outTree_->Branch("bestXYsigMatch_second",         &bestXYsigMatch_second,         "bestXYsigMatch_second/I");
  outTree_->Branch("bestXYsigMatchCat0",            &bestXYsigMatchCat0,            "bestXYsigMatchCat0/I");
  outTree_->Branch("bestXYsigMatchCat1",            &bestXYsigMatchCat1,            "bestXYsigMatchCat1/I");
  outTree_->Branch("bestXYsigMatchCat2",            &bestXYsigMatchCat2,            "bestXYsigMatchCat2/I");
  outTree_->Branch("bestXYsigMatchCat0_second",     &bestXYsigMatchCat0_second,     "bestXYsigMatchCat0_second/I");
  outTree_->Branch("bestXYsigMatchCat1_second",     &bestXYsigMatchCat1_second,     "bestXYsigMatchCat1_second/I");
  outTree_->Branch("bestXYsigMatchCat2_second",     &bestXYsigMatchCat2_second,     "bestXYsigMatchCat2_second/I");
  outTree_->Branch("bestXYsigMatchCatNew0",         &bestXYsigMatchCatNew0,         "bestXYsigMatchCatNew0/I");
  outTree_->Branch("bestXYsigMatchCatNew1",         &bestXYsigMatchCatNew1,         "bestXYsigMatchCatNew1/I");
  outTree_->Branch("bestXYsigMatchCatNew2",         &bestXYsigMatchCatNew2,         "bestXYsigMatchCatNew2/I");
  outTree_->Branch("bestXYsigMatchCatNew0_second",  &bestXYsigMatchCatNew0_second,  "bestXYsigMatchCatNew0_second/I");
  outTree_->Branch("bestXYsigMatchCatNew1_second",  &bestXYsigMatchCatNew1_second,  "bestXYsigMatchCatNew1_second/I");
  outTree_->Branch("bestXYsigMatchCatNew2_second",  &bestXYsigMatchCatNew2_second,  "bestXYsigMatchCatNew2_second/I");
  outTree_->Branch("bestXYsigMatch_causeEle1",      &bestXYsigMatch_causeEle1,      "bestXYsigMatch_causeEle1/I");
  outTree_->Branch("bestXYsigMatch_causeEle2",      &bestXYsigMatch_causeEle2,      "bestXYsigMatch_causeEle2/I");
  outTree_->Branch("bestXYsigMatch_notmatching",    &bestXYsigMatch_notmatching,    "bestXYsigMatch_notmatching/I");
  outTree_->Branch("bestXYsigMatch_causeK",         &bestXYsigMatch_causeK,         "bestXYsigMatch_causeK/I");
  outTree_->Branch("bestXYsigMatch_notok_ele1pt",   &bestXYsigMatch_notok_ele1pt,   "bestXYsigMatch_notok_ele1pt/F");
  outTree_->Branch("bestXYsigMatch_notok_ele2pt",   &bestXYsigMatch_notok_ele2pt,   "bestXYsigMatch_notok_ele2pt/F");
  outTree_->Branch("bestXYsigMatch_notok_kpt",      &bestXYsigMatch_notok_kpt,      "bestXYsigMatch_notok_kpt/F");
  outTree_->Branch("bestXYsigMatch_notok_minpt",    &bestXYsigMatch_notok_minpt,    "bestXYsigMatch_notok_minpt/F");
  outTree_->Branch("bestXYsigMatch_notok_ele1eta",  &bestXYsigMatch_notok_ele1eta,  "bestXYsigMatch_notok_ele1eta/F");
  outTree_->Branch("bestXYsigMatch_notok_ele2eta",  &bestXYsigMatch_notok_ele2eta,  "bestXYsigMatch_notok_ele2eta/F");
  outTree_->Branch("bestXYsigMatch_notok_keta",     &bestXYsigMatch_notok_keta,     "bestXYsigMatch_notok_keta/F");
  outTree_->Branch("bestXYsigMatch_notok_pfmva1",   &bestXYsigMatch_notok_pfmva1,   "bestXYsigMatch_notok_pfmva1/F");
  outTree_->Branch("bestXYsigMatch_notok_pfmva2",   &bestXYsigMatch_notok_pfmva2,   "bestXYsigMatch_notok_pfmva2/F");
  outTree_->Branch("bestXYsigMatch_notok_lptmva1",  &bestXYsigMatch_notok_lptmva1,  "bestXYsigMatch_notok_lptmva1/F");
  outTree_->Branch("bestXYsigMatch_notok_lptmva2",  &bestXYsigMatch_notok_lptmva2,  "bestXYsigMatch_notok_lptmva2/F");
  outTree_->Branch("bestXYsigMatch_notok_costhetaSK",   &bestXYsigMatch_notok_costhetaSK,   "bestXYsigMatch_notok_costhetaSK/F");
  outTree_->Branch("bestXYsigMatch_notok_costhetaSKCS", &bestXYsigMatch_notok_costhetaSKCS, "bestXYsigMatch_notok_costhetaSKCS/F");
  outTree_->Branch("bestXYsigMatch_notok_costhetaL",    &bestXYsigMatch_notok_costhetaL,    "bestXYsigMatch_notok_costhetaL/F");
  outTree_->Branch("bestXYsigMatch_ok_ele1pt",      &bestXYsigMatch_ok_ele1pt,      "bestXYsigMatch_ok_ele1pt/F");
  outTree_->Branch("bestXYsigMatch_ok_ele2pt",      &bestXYsigMatch_ok_ele2pt,      "bestXYsigMatch_ok_ele2pt/F");
  outTree_->Branch("bestXYsigMatch_ok_kpt",         &bestXYsigMatch_ok_kpt,         "bestXYsigMatch_ok_kpt/F");
  outTree_->Branch("bestXYsigMatch_ok_minpt",       &bestXYsigMatch_ok_minpt,       "bestXYsigMatch_ok_minpt/F");
  outTree_->Branch("bestXYsigMatch_ok_ele1eta",     &bestXYsigMatch_ok_ele1eta,     "bestXYsigMatch_ok_ele1eta/F");
  outTree_->Branch("bestXYsigMatch_ok_ele2eta",     &bestXYsigMatch_ok_ele2eta,     "bestXYsigMatch_ok_ele2eta/F");
  outTree_->Branch("bestXYsigMatch_ok_keta",        &bestXYsigMatch_ok_keta,        "bestXYsigMatch_ok_keta/F");
  outTree_->Branch("bestXYsigMatch_ok_pfmva1",      &bestXYsigMatch_ok_pfmva1,      "bestXYsigMatch_ok_pfmva1/F");
  outTree_->Branch("bestXYsigMatch_ok_pfmva2",      &bestXYsigMatch_ok_pfmva2,      "bestXYsigMatch_ok_pfmva2/F");
  outTree_->Branch("bestXYsigMatch_ok_lptmva1",     &bestXYsigMatch_ok_lptmva1,     "bestXYsigMatch_ok_lptmva1/F");
  outTree_->Branch("bestXYsigMatch_ok_lptmva2",     &bestXYsigMatch_ok_lptmva2,     "bestXYsigMatch_ok_lptmva2/F");
  outTree_->Branch("bestXYsigMatch_ok_costhetaSK",     &bestXYsigMatch_ok_costhetaSK,     "bestXYsigMatch_ok_costhetaSK/F");
  outTree_->Branch("bestXYsigMatch_ok_costhetaSKCS",   &bestXYsigMatch_ok_costhetaSKCS,   "bestXYsigMatch_ok_costhetaSKCS/F");
  outTree_->Branch("bestXYsigMatch_ok_costhetaL",      &bestXYsigMatch_ok_costhetaL,      "bestXYsigMatch_ok_costhetaL/F");
  outTree_->Branch("bestXYsigMatch_ok_costhetaSK_gen", &bestXYsigMatch_ok_costhetaSK_gen, "bestXYsigMatch_ok_costhetaSK_gen/F");

  outTree_->Branch("numberBetterSvProbTriplets",   &numberBetterSvProbTriplets,   "numberBetterSvProbTriplets/I");
  outTree_->Branch("numberBetterXYsigTriplets",    &numberBetterXYsigTriplets,    "numberBetterXYsigTriplets/I");
  outTree_->Branch("numberBetterCos2DTriplets",    &numberBetterCos2DTriplets,    "numberBetterCos2DTriplets/I");
  outTree_->Branch("numberBetterAllPtSumTriplets", &numberBetterAllPtSumTriplets, "numberBetterAllPtSumTriplets/I");

  outTree_->Branch("bestCos2DMatch",                &bestCos2DMatch,                "bestCos2DMatch/I");
  outTree_->Branch("bestCos2DMatch_second",         &bestCos2DMatch_second,         "bestCos2DMatch_second/I");
  outTree_->Branch("bestCos2DMatchCat0",            &bestCos2DMatchCat0,            "bestCos2DMatchCat0/I");
  outTree_->Branch("bestCos2DMatchCat1",            &bestCos2DMatchCat1,            "bestCos2DMatchCat1/I");
  outTree_->Branch("bestCos2DMatchCat2",            &bestCos2DMatchCat2,            "bestCos2DMatchCat2/I");
  outTree_->Branch("bestCos2DMatchCat0_second",     &bestCos2DMatchCat0_second,     "bestCos2DMatchCat0_second/I");
  outTree_->Branch("bestCos2DMatchCat1_second",     &bestCos2DMatchCat1_second,     "bestCos2DMatchCat1_second/I");
  outTree_->Branch("bestCos2DMatchCat2_second",     &bestCos2DMatchCat2_second,     "bestCos2DMatchCat2_second/I");
  outTree_->Branch("bestCos2DMatchCatNew0",         &bestCos2DMatchCatNew0,         "bestCos2DMatchCatNew0/I");
  outTree_->Branch("bestCos2DMatchCatNew1",         &bestCos2DMatchCatNew1,         "bestCos2DMatchCatNew1/I");
  outTree_->Branch("bestCos2DMatchCatNew2",         &bestCos2DMatchCatNew2,         "bestCos2DMatchCatNew2/I");
  outTree_->Branch("bestCos2DMatchCatNew0_second",  &bestCos2DMatchCatNew0_second,  "bestCos2DMatchCatNew0_second/I");
  outTree_->Branch("bestCos2DMatchCatNew1_second",  &bestCos2DMatchCatNew1_second,  "bestCos2DMatchCatNew1_second/I");
  outTree_->Branch("bestCos2DMatchCatNew2_second",  &bestCos2DMatchCatNew2_second,  "bestCos2DMatchCatNew2_second/I");
  outTree_->Branch("bestCos2DMatch_causeEle1",      &bestCos2DMatch_causeEle1,      "bestCos2DMatch_causeEle1/I");
  outTree_->Branch("bestCos2DMatch_causeEle2",      &bestCos2DMatch_causeEle2,      "bestCos2DMatch_causeEle2/I");
  outTree_->Branch("bestCos2DMatch_causeK",         &bestCos2DMatch_causeK,         "bestCos2DMatch_causeK/I");
  outTree_->Branch("bestCos2DMatch_notok_ele1pt",   &bestCos2DMatch_notok_ele1pt,   "bestCos2DMatch_notok_ele1pt/F");
  outTree_->Branch("bestCos2DMatch_notok_ele2pt",   &bestCos2DMatch_notok_ele2pt,   "bestCos2DMatch_notok_ele2pt/F");
  outTree_->Branch("bestCos2DMatch_notok_kpt",      &bestCos2DMatch_notok_kpt,      "bestCos2DMatch_notok_kpt/F");
  outTree_->Branch("bestCos2DMatch_notok_ele1eta",  &bestCos2DMatch_notok_ele1eta,  "bestCos2DMatch_notok_ele1eta/F");
  outTree_->Branch("bestCos2DMatch_notok_ele2eta",  &bestCos2DMatch_notok_ele2eta,  "bestCos2DMatch_notok_ele2eta/F");
  outTree_->Branch("bestCos2DMatch_notok_keta",     &bestCos2DMatch_notok_keta,     "bestCos2DMatch_notok_keta/F");
  outTree_->Branch("bestCos2DMatch_ok_ele1pt",      &bestCos2DMatch_ok_ele1pt,      "bestCos2DMatch_ok_ele1pt/F");
  outTree_->Branch("bestCos2DMatch_ok_ele2pt",      &bestCos2DMatch_ok_ele2pt,      "bestCos2DMatch_ok_ele2pt/F");
  outTree_->Branch("bestCos2DMatch_ok_kpt",         &bestCos2DMatch_ok_kpt,         "bestCos2DMatch_ok_kpt/F");
  outTree_->Branch("bestCos2DMatch_ok_ele1eta",     &bestCos2DMatch_ok_ele1eta,     "bestCos2DMatch_ok_ele1eta/F");
  outTree_->Branch("bestCos2DMatch_ok_ele2eta",     &bestCos2DMatch_ok_ele2eta,     "bestCos2DMatch_ok_ele2eta/F");
  outTree_->Branch("bestCos2DMatch_ok_keta",        &bestCos2DMatch_ok_keta,        "bestCos2DMatch_ok_keta/F");

  outTree_->Branch("bestAllPtSumMatch",            &bestAllPtSumMatch,            "bestAllPtSumMatch/I");
  outTree_->Branch("bestAllPtSumMatch_second",     &bestAllPtSumMatch_second,     "bestAllPtSumMatch_second/I");
  outTree_->Branch("bestAllPtSumMatchCat0",        &bestAllPtSumMatchCat0,        "bestAllPtSumMatchCat0/I");
  outTree_->Branch("bestAllPtSumMatchCat1",        &bestAllPtSumMatchCat1,        "bestAllPtSumMatchCat1/I");
  outTree_->Branch("bestAllPtSumMatchCat2",        &bestAllPtSumMatchCat2,        "bestAllPtSumMatchCat2/I");
  outTree_->Branch("bestAllPtSumMatchCatNew0",     &bestAllPtSumMatchCatNew0,     "bestAllPtSumMatchCatNew0/I");
  outTree_->Branch("bestAllPtSumMatchCatNew1",     &bestAllPtSumMatchCatNew1,     "bestAllPtSumMatchCatNew1/I");
  outTree_->Branch("bestAllPtSumMatchCatNew2",     &bestAllPtSumMatchCatNew2,     "bestAllPtSumMatchCatNew2/I");
  outTree_->Branch("bestAllPtSumMatchCat0_second",     &bestAllPtSumMatchCat0_second,     "bestAllPtSumMatchCat0_second/I");
  outTree_->Branch("bestAllPtSumMatchCat1_second",     &bestAllPtSumMatchCat1_second,     "bestAllPtSumMatchCat1_second/I");
  outTree_->Branch("bestAllPtSumMatchCat2_second",     &bestAllPtSumMatchCat2_second,     "bestAllPtSumMatchCat2_second/I");
  outTree_->Branch("bestAllPtSumMatchCatNew0_second",  &bestAllPtSumMatchCatNew0_second,  "bestAllPtSumMatchCatNew0_second/I");
  outTree_->Branch("bestAllPtSumMatchCatNew1_second",  &bestAllPtSumMatchCatNew1_second,  "bestAllPtSumMatchCatNew1_second/I");
  outTree_->Branch("bestAllPtSumMatchCatNew2_second",  &bestAllPtSumMatchCatNew2_second,  "bestAllPtSumMatchCatNew2_second/I");
  outTree_->Branch("bestAllPtSumMatch_causeEle1",  &bestAllPtSumMatch_causeEle1,  "bestAllPtSumMatch_causeEle1/I");
  outTree_->Branch("bestAllPtSumMatch_causeEle2",  &bestAllPtSumMatch_causeEle2,  "bestAllPtSumMatch_causeEle2/I");
  outTree_->Branch("bestAllPtSumMatch_causeK",     &bestAllPtSumMatch_causeK,     "bestAllPtSumMatch_causeK/I");
}

void TripletSelection::bookOutputHistos() 
{
  cout << "Booking histos" << endl;
  //
  h_entries   = new TH1F("h_entries",  "Number of entries",   3,  3.5, 6.5);
  h_selection = new TH1F("h_selection","Selection breakdown", 8, -0.5, 7.5);
}
