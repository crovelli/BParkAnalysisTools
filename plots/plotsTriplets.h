//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Mar 16 18:26:08 2020 by ROOT version 6.14/00
// from TTree TaPtree/TaPtree
// found on file: BuToKJpsi_Toee_BParkNANO_mc_2020Jan16.txt.root
//////////////////////////////////////////////////////////

#ifndef plotsTriplets_h
#define plotsTriplets_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"

class plotsTriplets {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   vector<float>   *debug_svprob;
   vector<int>     *debug_svprob_match;
   vector<float>   *debug_pf_svprob;
   vector<int>     *debug_pf_svprob_match;
   Int_t           theEvent;
   Float_t         rho;
   Int_t           goodBSize;
   Int_t           goodTrueBSize;
   Int_t           goodCombBSize;
   vector<float>   *goodTrueB_maxMinDREle;
   vector<float>   *goodTrueB_maxMinDREle_dEta;
   vector<float>   *goodTrueB_maxMinDREle_dPhi;
   vector<float>   *goodTrueB_maxMinDREle_dPtOverPt;
   vector<float>   *goodTrueB_dRgen;
   vector<float>   *goodTrueB_maxDRTrack;
   vector<float>   *goodCombB_maxMinDREle;
   vector<float>   *goodCombB_maxMinDREle_dEta;
   vector<float>   *goodCombB_maxMinDREle_dPhi;
   vector<float>   *goodCombB_maxMinDREle_dPtOverPt;
   vector<float>   *goodCombB_maxDRTrack;
   Int_t           goodTrueBs_SvProbMatch;
   Int_t           goodTrueBs_xySigMatch;
   Int_t           goodTrueBs_cos2DMatch;
   Int_t           goodTrueBs_ptSumMatch;
   Int_t           goodTrueBs_kptMatch;
   Float_t         goodTrueBs_bestSvProb_dRmax;
   Float_t         goodTrueBs_bestXYSig_dRmax;
   Float_t         goodTrueBs_bestCos2d_dRmax;
   Float_t         goodTrueBs_bestPtsum_dRmax;
   Float_t         goodTrueBs_bestKpt_dRmax;
   Float_t         goodTrueBs_bestSvProb_dRmin;
   Float_t         goodTrueBs_bestXYSig_dRmin;
   Float_t         goodTrueBs_bestCos2d_dRmin;
   Float_t         goodTrueBs_bestPtsum_dRmin;
   Float_t         goodTrueBs_bestKpt_dRmin;
   Float_t         bestMatch_SvProb;
   Float_t         bestMatch_XYSig;
   Float_t         bestMatch_Cos2D;
   Float_t         bestMatch_PtSum;
   Float_t         bestMatch_KPt;
   Float_t         bestMatch_maxDrRecoGen;
   Float_t         bestMatch_minDrRecoGen;
   Float_t         bestMatch_drRecoGenK;
   vector<float>   *goodTrueB_svProb_notBestMatch;
   vector<float>   *goodTrueB_xySig_notBestMatch;
   vector<float>   *goodTrueB_cos2D_notBestMatch;
   vector<float>   *goodTrueB_ptsum_notBestMatch;
   vector<float>   *goodTrueB_kpt_notBestMatch;
   vector<float>   *goodCombB_svProb;
   vector<float>   *goodCombB_xySig;
   vector<float>   *goodCombB_cos2D;
   vector<float>   *goodCombB_ptsum;
   vector<float>   *goodCombB_kpt;
   vector<float>   *goodCombB_maxDrRecoGen;
   vector<float>   *goodCombB_minDrRecoGen;
   vector<float>   *goodCombB_drRecoGenK;
   vector<float>   *combToBestB_svProb;
   vector<float>   *combToBestB_xySig;
   vector<float>   *combToBestB_cos2D;
   vector<float>   *combToBestB_ptsum;
   vector<float>   *combToBestB_kp;
   vector<float>   *combToBestB_maxDrRecoGen;
   vector<float>   *combToBestB_minDrRecoGen;
   vector<float>   *combToBestB_drRecoGenK;
   Int_t           bestSvProbMatch;
   Int_t           bestSvProbMatch_causeEle1;
   Int_t           bestSvProbMatch_causeEle2;
   Int_t           bestSvProbMatch_causeK;
   Float_t         bestSvProbMatch_notok_ele1pt;
   Float_t         bestSvProbMatch_notok_ele2pt;
   Float_t         bestSvProbMatch_notok_kpt;
   Float_t         bestSvProbMatch_notok_ele1eta;
   Float_t         bestSvProbMatch_notok_ele2eta;
   Float_t         bestSvProbMatch_notok_keta;
   Float_t         bestSvProbMatch_ok_ele1pt;
   Float_t         bestSvProbMatch_ok_ele2pt;
   Float_t         bestSvProbMatch_ok_kpt;
   Float_t         bestSvProbMatch_ok_ele1eta;
   Float_t         bestSvProbMatch_ok_ele2eta;
   Float_t         bestSvProbMatch_ok_keta;
   Int_t           bestXYsigMatch;
   Int_t           bestXYsigMatch_causeEle1;
   Int_t           bestXYsigMatch_causeEle2;
   Int_t           bestXYsigMatch_causeK;
   Float_t         bestXYsigMatch_notok_ele1pt;
   Float_t         bestXYsigMatch_notok_ele2pt;
   Float_t         bestXYsigMatch_notok_kpt;
   Float_t         bestXYsigMatch_notok_ele1eta;
   Float_t         bestXYsigMatch_notok_ele2eta;
   Float_t         bestXYsigMatch_notok_keta;
   Float_t         bestXYsigMatch_notok_pfmva1;
   Float_t         bestXYsigMatch_notok_pfmva2;
   Float_t         bestXYsigMatch_notok_lptmva1;
   Float_t         bestXYsigMatch_notok_lptmva2;
   Float_t         bestXYsigMatch_ok_ele1pt;
   Float_t         bestXYsigMatch_ok_ele2pt;
   Float_t         bestXYsigMatch_ok_kpt;
   Float_t         bestXYsigMatch_ok_ele1eta;
   Float_t         bestXYsigMatch_ok_ele2eta;
   Float_t         bestXYsigMatch_ok_keta;
   Float_t         bestXYsigMatch_ok_pfmva1;
   Float_t         bestXYsigMatch_ok_pfmva2;
   Float_t         bestXYsigMatch_ok_lptmva1;
   Float_t         bestXYsigMatch_ok_lptmva2;
   Int_t           bestCos2DMatch;
   Int_t           bestCos2DMatch_causeEle1;
   Int_t           bestCos2DMatch_causeEle2;
   Int_t           bestCos2DMatch_causeK;
   Float_t         bestCos2DMatch_notok_ele1pt;
   Float_t         bestCos2DMatch_notok_ele2pt;
   Float_t         bestCos2DMatch_notok_kpt;
   Float_t         bestCos2DMatch_notok_ele1eta;
   Float_t         bestCos2DMatch_notok_ele2eta;
   Float_t         bestCos2DMatch_notok_keta;
   Float_t         bestCos2DMatch_ok_ele1pt;
   Float_t         bestCos2DMatch_ok_ele2pt;
   Float_t         bestCos2DMatch_ok_kpt;
   Float_t         bestCos2DMatch_ok_ele1eta;
   Float_t         bestCos2DMatch_ok_ele2eta;
   Float_t         bestCos2DMatch_ok_keta;
   Int_t           bestPtSumMatch;
   Int_t           bestPtSumMatch_causeEle1;
   Int_t           bestPtSumMatch_causeEle2;
   Int_t           bestPtSumMatch_causeK;
   Int_t           bestKPtMatch;
   Int_t           bestKPtMatch_causeEle1;
   Int_t           bestKPtMatch_causeEle2;
   Int_t           bestKPtMatch_causeK;

   // List of branches
   TBranch        *b_debug_svprob;   //!
   TBranch        *b_debug_svprob_match;   //!
   TBranch        *b_debug_pf_svprob;   //!
   TBranch        *b_debug_pf_svprob_match;   //!
   TBranch        *b_theEvent;   //!
   TBranch        *b_rho;   //!
   TBranch        *b_goodBSize;   //!
   TBranch        *b_goodTrueBSize;   //!
   TBranch        *b_goodCombBSize;   //!
   TBranch        *b_goodTrueB_maxMinDREle;   //!
   TBranch        *b_goodTrueB_maxMinDREle_dEta;   //!
   TBranch        *b_goodTrueB_maxMinDREle_dPhi;   //!
   TBranch        *b_goodTrueB_maxMinDREle_dPtOverPt;   //!
   TBranch        *b_goodTrueB_dRgen;   //!
   TBranch        *b_goodTrueB_maxDRTrack;   //!
   TBranch        *b_goodCombB_maxMinDREle;   //!
   TBranch        *b_goodCombB_maxMinDREle_dEta;   //!
   TBranch        *b_goodCombB_maxMinDREle_dPhi;   //!
   TBranch        *b_goodCombB_maxMinDREle_dPtOverPt;   //!
   TBranch        *b_goodCombB_maxDRTrack;   //!
   TBranch        *b_goodTrueBs_SvProbMatch;   //!
   TBranch        *b_goodTrueBs_xySigMatch;   //!
   TBranch        *b_goodTrueBs_cos2DMatch;   //!
   TBranch        *b_goodTrueBs_ptSumMatch;   //!
   TBranch        *b_goodTrueBs_kptMatch;   //!
   TBranch        *b_goodTrueBs_bestSvProb_dRmax;   //!
   TBranch        *b_goodTrueBs_bestXYSig_dRmax;   //!
   TBranch        *b_goodTrueBs_bestCos2d_dRmax;   //!
   TBranch        *b_goodTrueBs_bestPtsum_dRmax;   //!
   TBranch        *b_goodTrueBs_bestKpt_dRmax;   //!
   TBranch        *b_goodTrueBs_bestSvProb_dRmin;   //!
   TBranch        *b_goodTrueBs_bestXYSig_dRmin;   //!
   TBranch        *b_goodTrueBs_bestCos2d_dRmin;   //!
   TBranch        *b_goodTrueBs_bestPtsum_dRmin;   //!
   TBranch        *b_goodTrueBs_bestKpt_dRmin;   //!
   TBranch        *b_bestMatch_SvProb;   //!
   TBranch        *b_bestMatch_XYSig;   //!
   TBranch        *b_bestMatch_Cos2D;   //!
   TBranch        *b_bestMatch_PtSum;   //!
   TBranch        *b_bestMatch_KPt;   //!
   TBranch        *b_bestMatch_maxDrRecoGen;   //!
   TBranch        *b_bestMatch_minDrRecoGen;   //!
   TBranch        *b_bestMatch_drRecoGenK;   //!
   TBranch        *b_goodTrueB_svProb_notBestMatch;   //!
   TBranch        *b_goodTrueB_xySig_notBestMatch;   //!
   TBranch        *b_goodTrueB_cos2D_notBestMatch;   //!
   TBranch        *b_goodTrueB_ptsum_notBestMatch;   //!
   TBranch        *b_goodTrueB_kpt_notBestMatch;   //!
   TBranch        *b_goodCombB_svProb;   //!
   TBranch        *b_goodCombB_xySig;   //!
   TBranch        *b_goodCombB_cos2D;   //!
   TBranch        *b_goodCombB_ptsum;   //!
   TBranch        *b_goodCombB_kpt;   //!
   TBranch        *b_goodCombB_maxDrRecoGen;   //!
   TBranch        *b_goodCombB_minDrRecoGen;   //!
   TBranch        *b_goodCombB_drRecoGenK;   //!
   TBranch        *b_combToBestB_svProb;   //!
   TBranch        *b_combToBestB_xySig;   //!
   TBranch        *b_combToBestB_cos2D;   //!
   TBranch        *b_combToBestB_ptsum;   //!
   TBranch        *b_combToBestB_kp;   //!
   TBranch        *b_combToBestB_maxDrRecoGen;   //!
   TBranch        *b_combToBestB_minDrRecoGen;   //!
   TBranch        *b_combToBestB_drRecoGenK;   //!
   TBranch        *b_bestSvProbMatch;   //!
   TBranch        *b_bestSvProbMatch_causeEle1;   //!
   TBranch        *b_bestSvProbMatch_causeEle2;   //!
   TBranch        *b_bestSvProbMatch_causeK;   //!
   TBranch        *b_bestSvProbMatch_notok_ele1pt;   //!
   TBranch        *b_bestSvProbMatch_notok_ele2pt;   //!
   TBranch        *b_bestSvProbMatch_notok_kpt;   //!
   TBranch        *b_bestSvProbMatch_notok_ele1eta;   //!
   TBranch        *b_bestSvProbMatch_notok_ele2eta;   //!
   TBranch        *b_bestSvProbMatch_notok_keta;   //!
   TBranch        *b_bestSvProbMatch_ok_ele1pt;   //!
   TBranch        *b_bestSvProbMatch_ok_ele2pt;   //!
   TBranch        *b_bestSvProbMatch_ok_kpt;   //!
   TBranch        *b_bestSvProbMatch_ok_ele1eta;   //!
   TBranch        *b_bestSvProbMatch_ok_ele2eta;   //!
   TBranch        *b_bestSvProbMatch_ok_keta;   //!
   TBranch        *b_bestXYsigMatch;   //!
   TBranch        *b_bestXYsigMatch_causeEle1;   //!
   TBranch        *b_bestXYsigMatch_causeEle2;   //!
   TBranch        *b_bestXYsigMatch_causeK;   //!
   TBranch        *b_bestXYsigMatch_notok_ele1pt;   //!
   TBranch        *b_bestXYsigMatch_notok_ele2pt;   //!
   TBranch        *b_bestXYsigMatch_notok_kpt;   //!
   TBranch        *b_bestXYsigMatch_notok_ele1eta;   //!
   TBranch        *b_bestXYsigMatch_notok_ele2eta;   //!
   TBranch        *b_bestXYsigMatch_notok_keta;   //!
   TBranch        *b_bestXYsigMatch_notok_pfmva1;   //!
   TBranch        *b_bestXYsigMatch_notok_pfmva2;   //!
   TBranch        *b_bestXYsigMatch_notok_lptmva1;   //!
   TBranch        *b_bestXYsigMatch_notok_lptmva2;   //!
   TBranch        *b_bestXYsigMatch_ok_ele1pt;   //!
   TBranch        *b_bestXYsigMatch_ok_ele2pt;   //!
   TBranch        *b_bestXYsigMatch_ok_kpt;   //!
   TBranch        *b_bestXYsigMatch_ok_ele1eta;   //!
   TBranch        *b_bestXYsigMatch_ok_ele2eta;   //!
   TBranch        *b_bestXYsigMatch_ok_keta;   //!
   TBranch        *b_bestXYsigMatch_ok_pfmva1;   //!
   TBranch        *b_bestXYsigMatch_ok_pfmva2;   //!
   TBranch        *b_bestXYsigMatch_ok_lptmva1;   //!
   TBranch        *b_bestXYsigMatch_ok_lptmva2;   //!
   TBranch        *b_bestCos2DMatch;   //!
   TBranch        *b_bestCos2DMatch_causeEle1;   //!
   TBranch        *b_bestCos2DMatch_causeEle2;   //!
   TBranch        *b_bestCos2DMatch_causeK;   //!
   TBranch        *b_bestCos2DMatch_notok_ele1pt;   //!
   TBranch        *b_bestCos2DMatch_notok_ele2pt;   //!
   TBranch        *b_bestCos2DMatch_notok_kpt;   //!
   TBranch        *b_bestCos2DMatch_notok_ele1eta;   //!
   TBranch        *b_bestCos2DMatch_notok_ele2eta;   //!
   TBranch        *b_bestCos2DMatch_notok_keta;   //!
   TBranch        *b_bestCos2DMatch_ok_ele1pt;   //!
   TBranch        *b_bestCos2DMatch_ok_ele2pt;   //!
   TBranch        *b_bestCos2DMatch_ok_kpt;   //!
   TBranch        *b_bestCos2DMatch_ok_ele1eta;   //!
   TBranch        *b_bestCos2DMatch_ok_ele2eta;   //!
   TBranch        *b_bestCos2DMatch_ok_keta;   //!
   TBranch        *b_bestPtSumMatch;   //!
   TBranch        *b_bestPtSumMatch_causeEle1;   //!
   TBranch        *b_bestPtSumMatch_causeEle2;   //!
   TBranch        *b_bestPtSumMatch_causeK;   //!
   TBranch        *b_bestKPtMatch;   //!
   TBranch        *b_bestKPtMatch_causeEle1;   //!
   TBranch        *b_bestKPtMatch_causeEle2;   //!
   TBranch        *b_bestKPtMatch_causeK;   //!

   plotsTriplets(TTree *tree=0);
   virtual ~plotsTriplets();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef plotsTriplets_cxx
plotsTriplets::plotsTriplets(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("BuToKJpsi_Toee_BParkNANO_mc_2020Jan16.txt.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("BuToKJpsi_Toee_BParkNANO_mc_2020Jan16.txt.root");
      }
      f->GetObject("TaPtree",tree);

   }
   Init(tree);
}

plotsTriplets::~plotsTriplets()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t plotsTriplets::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t plotsTriplets::LoadTree(Long64_t entry)
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

void plotsTriplets::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   debug_svprob = 0;
   debug_svprob_match = 0;
   debug_pf_svprob = 0;
   debug_pf_svprob_match = 0;
   goodTrueB_maxMinDREle = 0;
   goodTrueB_maxMinDREle_dEta = 0;
   goodTrueB_maxMinDREle_dPhi = 0;
   goodTrueB_maxMinDREle_dPtOverPt = 0;
   goodTrueB_dRgen = 0;
   goodTrueB_maxDRTrack = 0;
   goodCombB_maxMinDREle = 0;
   goodCombB_maxMinDREle_dEta = 0;
   goodCombB_maxMinDREle_dPhi = 0;
   goodCombB_maxMinDREle_dPtOverPt = 0;
   goodCombB_maxDRTrack = 0;
   goodTrueB_svProb_notBestMatch = 0;
   goodTrueB_xySig_notBestMatch = 0;
   goodTrueB_cos2D_notBestMatch = 0;
   goodTrueB_ptsum_notBestMatch = 0;
   goodTrueB_kpt_notBestMatch = 0;
   goodCombB_svProb = 0;
   goodCombB_xySig = 0;
   goodCombB_cos2D = 0;
   goodCombB_ptsum = 0;
   goodCombB_kpt = 0;
   goodCombB_maxDrRecoGen = 0;
   goodCombB_minDrRecoGen = 0;
   goodCombB_drRecoGenK = 0;
   combToBestB_svProb = 0;
   combToBestB_xySig = 0;
   combToBestB_cos2D = 0;
   combToBestB_ptsum = 0;
   combToBestB_kp = 0;
   combToBestB_maxDrRecoGen = 0;
   combToBestB_minDrRecoGen = 0;
   combToBestB_drRecoGenK = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("debug_svprob", &debug_svprob, &b_debug_svprob);
   fChain->SetBranchAddress("debug_svprob_match", &debug_svprob_match, &b_debug_svprob_match);
   fChain->SetBranchAddress("debug_pf_svprob", &debug_pf_svprob, &b_debug_pf_svprob);
   fChain->SetBranchAddress("debug_pf_svprob_match", &debug_pf_svprob_match, &b_debug_pf_svprob_match);
   fChain->SetBranchAddress("theEvent", &theEvent, &b_theEvent);
   fChain->SetBranchAddress("rho", &rho, &b_rho);
   fChain->SetBranchAddress("goodBSize", &goodBSize, &b_goodBSize);
   fChain->SetBranchAddress("goodTrueBSize", &goodTrueBSize, &b_goodTrueBSize);
   fChain->SetBranchAddress("goodCombBSize", &goodCombBSize, &b_goodCombBSize);
   fChain->SetBranchAddress("goodTrueB_maxMinDREle", &goodTrueB_maxMinDREle, &b_goodTrueB_maxMinDREle);
   fChain->SetBranchAddress("goodTrueB_maxMinDREle_dEta", &goodTrueB_maxMinDREle_dEta, &b_goodTrueB_maxMinDREle_dEta);
   fChain->SetBranchAddress("goodTrueB_maxMinDREle_dPhi", &goodTrueB_maxMinDREle_dPhi, &b_goodTrueB_maxMinDREle_dPhi);
   fChain->SetBranchAddress("goodTrueB_maxMinDREle_dPtOverPt", &goodTrueB_maxMinDREle_dPtOverPt, &b_goodTrueB_maxMinDREle_dPtOverPt);
   fChain->SetBranchAddress("goodTrueB_dRgen", &goodTrueB_dRgen, &b_goodTrueB_dRgen);
   fChain->SetBranchAddress("goodTrueB_maxDRTrack", &goodTrueB_maxDRTrack, &b_goodTrueB_maxDRTrack);
   fChain->SetBranchAddress("goodCombB_maxMinDREle", &goodCombB_maxMinDREle, &b_goodCombB_maxMinDREle);
   fChain->SetBranchAddress("goodCombB_maxMinDREle_dEta", &goodCombB_maxMinDREle_dEta, &b_goodCombB_maxMinDREle_dEta);
   fChain->SetBranchAddress("goodCombB_maxMinDREle_dPhi", &goodCombB_maxMinDREle_dPhi, &b_goodCombB_maxMinDREle_dPhi);
   fChain->SetBranchAddress("goodCombB_maxMinDREle_dPtOverPt", &goodCombB_maxMinDREle_dPtOverPt, &b_goodCombB_maxMinDREle_dPtOverPt);
   fChain->SetBranchAddress("goodCombB_maxDRTrack", &goodCombB_maxDRTrack, &b_goodCombB_maxDRTrack);
   fChain->SetBranchAddress("goodTrueBs_SvProbMatch", &goodTrueBs_SvProbMatch, &b_goodTrueBs_SvProbMatch);
   fChain->SetBranchAddress("goodTrueBs_xySigMatch", &goodTrueBs_xySigMatch, &b_goodTrueBs_xySigMatch);
   fChain->SetBranchAddress("goodTrueBs_cos2DMatch", &goodTrueBs_cos2DMatch, &b_goodTrueBs_cos2DMatch);
   fChain->SetBranchAddress("goodTrueBs_ptSumMatch", &goodTrueBs_ptSumMatch, &b_goodTrueBs_ptSumMatch);
   fChain->SetBranchAddress("goodTrueBs_kptMatch", &goodTrueBs_kptMatch, &b_goodTrueBs_kptMatch);
   fChain->SetBranchAddress("goodTrueBs_bestSvProb_dRmax", &goodTrueBs_bestSvProb_dRmax, &b_goodTrueBs_bestSvProb_dRmax);
   fChain->SetBranchAddress("goodTrueBs_bestXYSig_dRmax", &goodTrueBs_bestXYSig_dRmax, &b_goodTrueBs_bestXYSig_dRmax);
   fChain->SetBranchAddress("goodTrueBs_bestCos2d_dRmax", &goodTrueBs_bestCos2d_dRmax, &b_goodTrueBs_bestCos2d_dRmax);
   fChain->SetBranchAddress("goodTrueBs_bestPtsum_dRmax", &goodTrueBs_bestPtsum_dRmax, &b_goodTrueBs_bestPtsum_dRmax);
   fChain->SetBranchAddress("goodTrueBs_bestKpt_dRmax", &goodTrueBs_bestKpt_dRmax, &b_goodTrueBs_bestKpt_dRmax);
   fChain->SetBranchAddress("goodTrueBs_bestSvProb_dRmin", &goodTrueBs_bestSvProb_dRmin, &b_goodTrueBs_bestSvProb_dRmin);
   fChain->SetBranchAddress("goodTrueBs_bestXYSig_dRmin", &goodTrueBs_bestXYSig_dRmin, &b_goodTrueBs_bestXYSig_dRmin);
   fChain->SetBranchAddress("goodTrueBs_bestCos2d_dRmin", &goodTrueBs_bestCos2d_dRmin, &b_goodTrueBs_bestCos2d_dRmin);
   fChain->SetBranchAddress("goodTrueBs_bestPtsum_dRmin", &goodTrueBs_bestPtsum_dRmin, &b_goodTrueBs_bestPtsum_dRmin);
   fChain->SetBranchAddress("goodTrueBs_bestKpt_dRmin", &goodTrueBs_bestKpt_dRmin, &b_goodTrueBs_bestKpt_dRmin);
   fChain->SetBranchAddress("bestMatch_SvProb", &bestMatch_SvProb, &b_bestMatch_SvProb);
   fChain->SetBranchAddress("bestMatch_XYSig", &bestMatch_XYSig, &b_bestMatch_XYSig);
   fChain->SetBranchAddress("bestMatch_Cos2D", &bestMatch_Cos2D, &b_bestMatch_Cos2D);
   fChain->SetBranchAddress("bestMatch_PtSum", &bestMatch_PtSum, &b_bestMatch_PtSum);
   fChain->SetBranchAddress("bestMatch_KPt", &bestMatch_KPt, &b_bestMatch_KPt);
   fChain->SetBranchAddress("bestMatch_maxDrRecoGen", &bestMatch_maxDrRecoGen, &b_bestMatch_maxDrRecoGen);
   fChain->SetBranchAddress("bestMatch_minDrRecoGen", &bestMatch_minDrRecoGen, &b_bestMatch_minDrRecoGen);
   fChain->SetBranchAddress("bestMatch_drRecoGenK", &bestMatch_drRecoGenK, &b_bestMatch_drRecoGenK);
   fChain->SetBranchAddress("goodTrueB_svProb_notBestMatch", &goodTrueB_svProb_notBestMatch, &b_goodTrueB_svProb_notBestMatch);
   fChain->SetBranchAddress("goodTrueB_xySig_notBestMatch", &goodTrueB_xySig_notBestMatch, &b_goodTrueB_xySig_notBestMatch);
   fChain->SetBranchAddress("goodTrueB_cos2D_notBestMatch", &goodTrueB_cos2D_notBestMatch, &b_goodTrueB_cos2D_notBestMatch);
   fChain->SetBranchAddress("goodTrueB_ptsum_notBestMatch", &goodTrueB_ptsum_notBestMatch, &b_goodTrueB_ptsum_notBestMatch);
   fChain->SetBranchAddress("goodTrueB_kpt_notBestMatch", &goodTrueB_kpt_notBestMatch, &b_goodTrueB_kpt_notBestMatch);
   fChain->SetBranchAddress("goodCombB_svProb", &goodCombB_svProb, &b_goodCombB_svProb);
   fChain->SetBranchAddress("goodCombB_xySig", &goodCombB_xySig, &b_goodCombB_xySig);
   fChain->SetBranchAddress("goodCombB_cos2D", &goodCombB_cos2D, &b_goodCombB_cos2D);
   fChain->SetBranchAddress("goodCombB_ptsum", &goodCombB_ptsum, &b_goodCombB_ptsum);
   fChain->SetBranchAddress("goodCombB_kpt", &goodCombB_kpt, &b_goodCombB_kpt);
   fChain->SetBranchAddress("goodCombB_maxDrRecoGen", &goodCombB_maxDrRecoGen, &b_goodCombB_maxDrRecoGen);
   fChain->SetBranchAddress("goodCombB_minDrRecoGen", &goodCombB_minDrRecoGen, &b_goodCombB_minDrRecoGen);
   fChain->SetBranchAddress("goodCombB_drRecoGenK", &goodCombB_drRecoGenK, &b_goodCombB_drRecoGenK);
   fChain->SetBranchAddress("combToBestB_svProb", &combToBestB_svProb, &b_combToBestB_svProb);
   fChain->SetBranchAddress("combToBestB_xySig", &combToBestB_xySig, &b_combToBestB_xySig);
   fChain->SetBranchAddress("combToBestB_cos2D", &combToBestB_cos2D, &b_combToBestB_cos2D);
   fChain->SetBranchAddress("combToBestB_ptsum", &combToBestB_ptsum, &b_combToBestB_ptsum);
   fChain->SetBranchAddress("combToBestB_kp", &combToBestB_kp, &b_combToBestB_kp);
   fChain->SetBranchAddress("combToBestB_maxDrRecoGen", &combToBestB_maxDrRecoGen, &b_combToBestB_maxDrRecoGen);
   fChain->SetBranchAddress("combToBestB_minDrRecoGen", &combToBestB_minDrRecoGen, &b_combToBestB_minDrRecoGen);
   fChain->SetBranchAddress("combToBestB_drRecoGenK", &combToBestB_drRecoGenK, &b_combToBestB_drRecoGenK);
   fChain->SetBranchAddress("bestSvProbMatch", &bestSvProbMatch, &b_bestSvProbMatch);
   fChain->SetBranchAddress("bestSvProbMatch_causeEle1", &bestSvProbMatch_causeEle1, &b_bestSvProbMatch_causeEle1);
   fChain->SetBranchAddress("bestSvProbMatch_causeEle2", &bestSvProbMatch_causeEle2, &b_bestSvProbMatch_causeEle2);
   fChain->SetBranchAddress("bestSvProbMatch_causeK", &bestSvProbMatch_causeK, &b_bestSvProbMatch_causeK);
   fChain->SetBranchAddress("bestSvProbMatch_notok_ele1pt", &bestSvProbMatch_notok_ele1pt, &b_bestSvProbMatch_notok_ele1pt);
   fChain->SetBranchAddress("bestSvProbMatch_notok_ele2pt", &bestSvProbMatch_notok_ele2pt, &b_bestSvProbMatch_notok_ele2pt);
   fChain->SetBranchAddress("bestSvProbMatch_notok_kpt", &bestSvProbMatch_notok_kpt, &b_bestSvProbMatch_notok_kpt);
   fChain->SetBranchAddress("bestSvProbMatch_notok_ele1eta", &bestSvProbMatch_notok_ele1eta, &b_bestSvProbMatch_notok_ele1eta);
   fChain->SetBranchAddress("bestSvProbMatch_notok_ele2eta", &bestSvProbMatch_notok_ele2eta, &b_bestSvProbMatch_notok_ele2eta);
   fChain->SetBranchAddress("bestSvProbMatch_notok_keta", &bestSvProbMatch_notok_keta, &b_bestSvProbMatch_notok_keta);
   fChain->SetBranchAddress("bestSvProbMatch_ok_ele1pt", &bestSvProbMatch_ok_ele1pt, &b_bestSvProbMatch_ok_ele1pt);
   fChain->SetBranchAddress("bestSvProbMatch_ok_ele2pt", &bestSvProbMatch_ok_ele2pt, &b_bestSvProbMatch_ok_ele2pt);
   fChain->SetBranchAddress("bestSvProbMatch_ok_kpt", &bestSvProbMatch_ok_kpt, &b_bestSvProbMatch_ok_kpt);
   fChain->SetBranchAddress("bestSvProbMatch_ok_ele1eta", &bestSvProbMatch_ok_ele1eta, &b_bestSvProbMatch_ok_ele1eta);
   fChain->SetBranchAddress("bestSvProbMatch_ok_ele2eta", &bestSvProbMatch_ok_ele2eta, &b_bestSvProbMatch_ok_ele2eta);
   fChain->SetBranchAddress("bestSvProbMatch_ok_keta", &bestSvProbMatch_ok_keta, &b_bestSvProbMatch_ok_keta);
   fChain->SetBranchAddress("bestXYsigMatch", &bestXYsigMatch, &b_bestXYsigMatch);
   fChain->SetBranchAddress("bestXYsigMatch_causeEle1", &bestXYsigMatch_causeEle1, &b_bestXYsigMatch_causeEle1);
   fChain->SetBranchAddress("bestXYsigMatch_causeEle2", &bestXYsigMatch_causeEle2, &b_bestXYsigMatch_causeEle2);
   fChain->SetBranchAddress("bestXYsigMatch_causeK", &bestXYsigMatch_causeK, &b_bestXYsigMatch_causeK);
   fChain->SetBranchAddress("bestXYsigMatch_notok_ele1pt", &bestXYsigMatch_notok_ele1pt, &b_bestXYsigMatch_notok_ele1pt);
   fChain->SetBranchAddress("bestXYsigMatch_notok_ele2pt", &bestXYsigMatch_notok_ele2pt, &b_bestXYsigMatch_notok_ele2pt);
   fChain->SetBranchAddress("bestXYsigMatch_notok_kpt", &bestXYsigMatch_notok_kpt, &b_bestXYsigMatch_notok_kpt);
   fChain->SetBranchAddress("bestXYsigMatch_notok_ele1eta", &bestXYsigMatch_notok_ele1eta, &b_bestXYsigMatch_notok_ele1eta);
   fChain->SetBranchAddress("bestXYsigMatch_notok_ele2eta", &bestXYsigMatch_notok_ele2eta, &b_bestXYsigMatch_notok_ele2eta);
   fChain->SetBranchAddress("bestXYsigMatch_notok_keta", &bestXYsigMatch_notok_keta, &b_bestXYsigMatch_notok_keta);
   fChain->SetBranchAddress("bestXYsigMatch_notok_pfmva1", &bestXYsigMatch_notok_pfmva1, &b_bestXYsigMatch_notok_pfmva1);
   fChain->SetBranchAddress("bestXYsigMatch_notok_pfmva2", &bestXYsigMatch_notok_pfmva2, &b_bestXYsigMatch_notok_pfmva2);
   fChain->SetBranchAddress("bestXYsigMatch_notok_lptmva1", &bestXYsigMatch_notok_lptmva1, &b_bestXYsigMatch_notok_lptmva1);
   fChain->SetBranchAddress("bestXYsigMatch_notok_lptmva2", &bestXYsigMatch_notok_lptmva2, &b_bestXYsigMatch_notok_lptmva2);
   fChain->SetBranchAddress("bestXYsigMatch_ok_ele1pt", &bestXYsigMatch_ok_ele1pt, &b_bestXYsigMatch_ok_ele1pt);
   fChain->SetBranchAddress("bestXYsigMatch_ok_ele2pt", &bestXYsigMatch_ok_ele2pt, &b_bestXYsigMatch_ok_ele2pt);
   fChain->SetBranchAddress("bestXYsigMatch_ok_kpt", &bestXYsigMatch_ok_kpt, &b_bestXYsigMatch_ok_kpt);
   fChain->SetBranchAddress("bestXYsigMatch_ok_ele1eta", &bestXYsigMatch_ok_ele1eta, &b_bestXYsigMatch_ok_ele1eta);
   fChain->SetBranchAddress("bestXYsigMatch_ok_ele2eta", &bestXYsigMatch_ok_ele2eta, &b_bestXYsigMatch_ok_ele2eta);
   fChain->SetBranchAddress("bestXYsigMatch_ok_keta", &bestXYsigMatch_ok_keta, &b_bestXYsigMatch_ok_keta);
   fChain->SetBranchAddress("bestXYsigMatch_ok_pfmva1", &bestXYsigMatch_ok_pfmva1, &b_bestXYsigMatch_ok_pfmva1);
   fChain->SetBranchAddress("bestXYsigMatch_ok_pfmva2", &bestXYsigMatch_ok_pfmva2, &b_bestXYsigMatch_ok_pfmva2);
   fChain->SetBranchAddress("bestXYsigMatch_ok_lptmva1", &bestXYsigMatch_ok_lptmva1, &b_bestXYsigMatch_ok_lptmva1);
   fChain->SetBranchAddress("bestXYsigMatch_ok_lptmva2", &bestXYsigMatch_ok_lptmva2, &b_bestXYsigMatch_ok_lptmva2);
   fChain->SetBranchAddress("bestCos2DMatch", &bestCos2DMatch, &b_bestCos2DMatch);
   fChain->SetBranchAddress("bestCos2DMatch_causeEle1", &bestCos2DMatch_causeEle1, &b_bestCos2DMatch_causeEle1);
   fChain->SetBranchAddress("bestCos2DMatch_causeEle2", &bestCos2DMatch_causeEle2, &b_bestCos2DMatch_causeEle2);
   fChain->SetBranchAddress("bestCos2DMatch_causeK", &bestCos2DMatch_causeK, &b_bestCos2DMatch_causeK);
   fChain->SetBranchAddress("bestCos2DMatch_notok_ele1pt", &bestCos2DMatch_notok_ele1pt, &b_bestCos2DMatch_notok_ele1pt);
   fChain->SetBranchAddress("bestCos2DMatch_notok_ele2pt", &bestCos2DMatch_notok_ele2pt, &b_bestCos2DMatch_notok_ele2pt);
   fChain->SetBranchAddress("bestCos2DMatch_notok_kpt", &bestCos2DMatch_notok_kpt, &b_bestCos2DMatch_notok_kpt);
   fChain->SetBranchAddress("bestCos2DMatch_notok_ele1eta", &bestCos2DMatch_notok_ele1eta, &b_bestCos2DMatch_notok_ele1eta);
   fChain->SetBranchAddress("bestCos2DMatch_notok_ele2eta", &bestCos2DMatch_notok_ele2eta, &b_bestCos2DMatch_notok_ele2eta);
   fChain->SetBranchAddress("bestCos2DMatch_notok_keta", &bestCos2DMatch_notok_keta, &b_bestCos2DMatch_notok_keta);
   fChain->SetBranchAddress("bestCos2DMatch_ok_ele1pt", &bestCos2DMatch_ok_ele1pt, &b_bestCos2DMatch_ok_ele1pt);
   fChain->SetBranchAddress("bestCos2DMatch_ok_ele2pt", &bestCos2DMatch_ok_ele2pt, &b_bestCos2DMatch_ok_ele2pt);
   fChain->SetBranchAddress("bestCos2DMatch_ok_kpt", &bestCos2DMatch_ok_kpt, &b_bestCos2DMatch_ok_kpt);
   fChain->SetBranchAddress("bestCos2DMatch_ok_ele1eta", &bestCos2DMatch_ok_ele1eta, &b_bestCos2DMatch_ok_ele1eta);
   fChain->SetBranchAddress("bestCos2DMatch_ok_ele2eta", &bestCos2DMatch_ok_ele2eta, &b_bestCos2DMatch_ok_ele2eta);
   fChain->SetBranchAddress("bestCos2DMatch_ok_keta", &bestCos2DMatch_ok_keta, &b_bestCos2DMatch_ok_keta);
   fChain->SetBranchAddress("bestPtSumMatch", &bestPtSumMatch, &b_bestPtSumMatch);
   fChain->SetBranchAddress("bestPtSumMatch_causeEle1", &bestPtSumMatch_causeEle1, &b_bestPtSumMatch_causeEle1);
   fChain->SetBranchAddress("bestPtSumMatch_causeEle2", &bestPtSumMatch_causeEle2, &b_bestPtSumMatch_causeEle2);
   fChain->SetBranchAddress("bestPtSumMatch_causeK", &bestPtSumMatch_causeK, &b_bestPtSumMatch_causeK);
   fChain->SetBranchAddress("bestKPtMatch", &bestKPtMatch, &b_bestKPtMatch);
   fChain->SetBranchAddress("bestKPtMatch_causeEle1", &bestKPtMatch_causeEle1, &b_bestKPtMatch_causeEle1);
   fChain->SetBranchAddress("bestKPtMatch_causeEle2", &bestKPtMatch_causeEle2, &b_bestKPtMatch_causeEle2);
   fChain->SetBranchAddress("bestKPtMatch_causeK", &bestKPtMatch_causeK, &b_bestKPtMatch_causeK);
   Notify();
}

Bool_t plotsTriplets::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void plotsTriplets::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t plotsTriplets::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef plotsTriplets_cxx
