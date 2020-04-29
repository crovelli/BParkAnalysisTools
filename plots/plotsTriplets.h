//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat Apr 18 12:17:47 2020 by ROOT version 6.20/02
// from TTree TaPtree/TaPtree
// found on file: ../data/Jan16prod/BuToKee_Toee_BParkNANO_mc_2020Jan16.txt.root
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
   Int_t           iHLT_Mu12_IP6;
   Int_t           iHLT_Mu9_IP6;
   Int_t           theEvent;
   Float_t         rho;
   Int_t           rightMcTruth;
   Float_t         gen_costhetaSKRad;
   Float_t         gen_costhetaSKCS;
   Float_t         gen_costhetaSKHel;
   Int_t           goodBSize;
   Int_t           goodTrueBSize;
   Int_t           goodCombBSize;
   Int_t           handMadeB;
   Float_t         handMadeBmass;
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
   Float_t         bestMatch_Bmass;
   Float_t         bestMatch_SvProb;
   Float_t         bestMatch_XYSig;
   Float_t         bestMatch_Cos2D;
   Float_t         bestMatch_PtSum;
   Float_t         bestMatch_KPt;
   Float_t         bestMatch_KEta;
   Float_t         bestMatch_Ele1Pt;
   Float_t         bestMatch_Ele2Pt;
   Float_t         bestMatch_costhetaSKRad;
   Float_t         bestMatch_costhetaSKCS;
   Float_t         bestMatch_costhetaSKHel;
   Float_t         bestMatch_MinPt;
   Float_t         bestMatch_Ele1Eta;
   Float_t         bestMatch_Ele2Eta;
   Float_t         bestMatch_Ele1pfmva;
   Float_t         bestMatch_Ele2pfmva;
   Float_t         bestMatch_Ele1lptmva;
   Float_t         bestMatch_Ele2lptmva;
   Float_t         bestMatch_maxDrRecoGen;
   Float_t         bestMatch_minDrRecoGen;
   Float_t         bestMatch_drRecoGenK;
   Float_t         bestMatch_maxDrRecoGenFromB;
   Float_t         bestMatch_minDrRecoGenFromB;
   Float_t         bestMatch_drRecoGenFromBK;
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
   vector<float>   *goodCombB_keta;
   vector<float>   *goodCombB_ele1pt;
   vector<float>   *goodCombB_ele2pt;
   vector<float>   *goodCombB_costhetaSKRad;
   vector<float>   *goodCombB_costhetaSKCS;
   vector<float>   *goodCombB_costhetaSKHel;
   vector<float>   *goodCombB_minpt;
   vector<float>   *goodCombB_ele1eta;
   vector<float>   *goodCombB_ele2eta;
   vector<float>   *goodCombB_ele1pfmva;
   vector<float>   *goodCombB_ele2pfmva;
   vector<float>   *goodCombB_ele1lptmva;
   vector<float>   *goodCombB_ele2lptmva;
   vector<float>   *goodCombB_causeEle1;
   vector<float>   *goodCombB_causeEle2;
   vector<float>   *goodCombB_causeK;
   vector<int>     *goodCombB_notmatching;
   vector<float>   *goodCombB_maxDrRecoGen;
   vector<float>   *goodCombB_minDrRecoGen;
   vector<float>   *goodCombB_drRecoGenK;
   vector<float>   *goodCombB_maxDrRecoGenFromB;
   vector<float>   *goodCombB_minDrRecoGenFromB;
   vector<float>   *goodCombB_drRecoGenFromBK;
   Int_t           bestSvProbMatch;
   Int_t           bestSvProbMatch_second;
   Int_t           bestSvProbMatchCat0;
   Int_t           bestSvProbMatchCat1;
   Int_t           bestSvProbMatchCat2;
   Int_t           bestSvProbMatchCat0_second;
   Int_t           bestSvProbMatchCat1_second;
   Int_t           bestSvProbMatchCat2_second;
   Int_t           bestSvProbMatchCatNew0;
   Int_t           bestSvProbMatchCatNew1;
   Int_t           bestSvProbMatchCatNew2;
   Int_t           bestSvProbMatchCatNew0_second;
   Int_t           bestSvProbMatchCatNew1_second;
   Int_t           bestSvProbMatchCatNew2_second;
   Int_t           bestSvProbMatch_causeEle1;
   Int_t           bestSvProbMatch_causeEle2;
   Int_t           bestSvProbMatch_notmatching;
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
   Int_t           bestXYsigMatch_second;
   Int_t           bestXYsigMatchCat0;
   Int_t           bestXYsigMatchCat1;
   Int_t           bestXYsigMatchCat2;
   Int_t           bestXYsigMatchCat0_second;
   Int_t           bestXYsigMatchCat1_second;
   Int_t           bestXYsigMatchCat2_second;
   Int_t           bestXYsigMatchCatNew0;
   Int_t           bestXYsigMatchCatNew1;
   Int_t           bestXYsigMatchCatNew2;
   Int_t           bestXYsigMatchCatNew0_second;
   Int_t           bestXYsigMatchCatNew1_second;
   Int_t           bestXYsigMatchCatNew2_second;
   Int_t           bestXYsigMatch_causeEle1;
   Int_t           bestXYsigMatch_causeEle2;
   Int_t           bestXYsigMatch_notmatching;
   Int_t           bestXYsigMatch_causeK;
   Float_t         bestXYsigMatch_notok_ele1pt;
   Float_t         bestXYsigMatch_notok_ele2pt;
   Float_t         bestXYsigMatch_notok_kpt;
   Float_t         bestXYsigMatch_notok_minpt;
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
   Float_t         bestXYsigMatch_ok_minpt;
   Float_t         bestXYsigMatch_ok_ele1eta;
   Float_t         bestXYsigMatch_ok_ele2eta;
   Float_t         bestXYsigMatch_ok_keta;
   Float_t         bestXYsigMatch_ok_pfmva1;
   Float_t         bestXYsigMatch_ok_pfmva2;
   Float_t         bestXYsigMatch_ok_lptmva1;
   Float_t         bestXYsigMatch_ok_lptmva2;
   Int_t           numberBetterSvProbTriplets;
   Int_t           numberBetterXYsigTriplets;
   Int_t           numberBetterCos2DTriplets;
   Int_t           numberBetterAllPtSumTriplets;
   Int_t           numberBetterBPtTriplets;
   Int_t           bestCos2DMatch;
   Int_t           bestCos2DMatch_second;
   Int_t           bestCos2DMatchCat0;
   Int_t           bestCos2DMatchCat1;
   Int_t           bestCos2DMatchCat2;
   Int_t           bestCos2DMatchCat0_second;
   Int_t           bestCos2DMatchCat1_second;
   Int_t           bestCos2DMatchCat2_second;
   Int_t           bestCos2DMatchCatNew0;
   Int_t           bestCos2DMatchCatNew1;
   Int_t           bestCos2DMatchCatNew2;
   Int_t           bestCos2DMatchCatNew0_second;
   Int_t           bestCos2DMatchCatNew1_second;
   Int_t           bestCos2DMatchCatNew2_second;
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
   Int_t           bestAllPtSumMatch;
   Int_t           bestAllPtSumMatch_second;
   Int_t           bestAllPtSumMatchCat0;
   Int_t           bestAllPtSumMatchCat1;
   Int_t           bestAllPtSumMatchCat2;
   Int_t           bestAllPtSumMatchCatNew0;
   Int_t           bestAllPtSumMatchCatNew1;
   Int_t           bestAllPtSumMatchCatNew2;
   Int_t           bestAllPtSumMatchCat0_second;
   Int_t           bestAllPtSumMatchCat1_second;
   Int_t           bestAllPtSumMatchCat2_second;
   Int_t           bestAllPtSumMatchCatNew0_second;
   Int_t           bestAllPtSumMatchCatNew1_second;
   Int_t           bestAllPtSumMatchCatNew2_second;
   Int_t           bestAllPtSumMatch_causeEle1;
   Int_t           bestAllPtSumMatch_causeEle2;
   Int_t           bestAllPtSumMatch_causeK;
   Int_t           bestBPtMatch;
   Int_t           bestBPtMatch_second;
   Int_t           bestBPtMatchCat0;
   Int_t           bestBPtMatchCat1;
   Int_t           bestBPtMatchCat2;
   Int_t           bestBPtMatchCatNew0;
   Int_t           bestBPtMatchCatNew1;
   Int_t           bestBPtMatchCatNew2;
   Int_t           bestBPtMatchCat0_second;
   Int_t           bestBPtMatchCat1_second;
   Int_t           bestBPtMatchCat2_second;
   Int_t           bestBPtMatchCatNew0_second;
   Int_t           bestBPtMatchCatNew1_second;
   Int_t           bestBPtMatchCatNew2_second;
   Int_t           bestBPtMatch_causeEle1;
   Int_t           bestBPtMatch_causeEle2;
   Int_t           bestBPtMatch_causeK;

   // List of branches
   TBranch        *b_debug_svprob;   //!
   TBranch        *b_debug_svprob_match;   //!
   TBranch        *b_debug_pf_svprob;   //!
   TBranch        *b_debug_pf_svprob_match;   //!
   TBranch        *b_iHLT_Mu12_IP6;   //!
   TBranch        *b_iHLT_Mu9_IP6;   //!
   TBranch        *b_theEvent;   //!
   TBranch        *b_rho;   //!
   TBranch        *b_rightMcTruth;   //!
   TBranch        *b_gen_costhetaSKRad; 
   TBranch        *b_gen_costhetaSKCS; 
   TBranch        *b_gen_costhetaSKHel; 
   TBranch        *b_goodBSize;   //!
   TBranch        *b_goodTrueBSize;   //!
   TBranch        *b_goodCombBSize;   //!
   TBranch        *b_handMadeB;   //!
   TBranch        *b_handMadeBmass;   //!
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
   TBranch        *b_bestMatch_Bmass;   //!
   TBranch        *b_bestMatch_SvProb;   //!
   TBranch        *b_bestMatch_XYSig;   //!
   TBranch        *b_bestMatch_Cos2D;   //!
   TBranch        *b_bestMatch_PtSum;   //!
   TBranch        *b_bestMatch_KPt;   //!
   TBranch        *b_bestMatch_KEta;   //!
   TBranch        *b_bestMatch_Ele1Pt;   //!
   TBranch        *b_bestMatch_Ele2Pt;   //!
   TBranch        *b_bestMatch_costhetaSKRad;   //!
   TBranch        *b_bestMatch_costhetaSKCS;   //!
   TBranch        *b_bestMatch_costhetaSKHel;   //!
   TBranch        *b_bestMatch_MinPt;   //!
   TBranch        *b_bestMatch_Ele1Eta;   //!
   TBranch        *b_bestMatch_Ele2Eta;   //!
   TBranch        *b_bestMatch_Ele1pfmva;   //!
   TBranch        *b_bestMatch_Ele2pfmva;   //!
   TBranch        *b_bestMatch_Ele1lptmva;   //!
   TBranch        *b_bestMatch_Ele2lptmva;   //!
   TBranch        *b_bestMatch_maxDrRecoGen;   //!
   TBranch        *b_bestMatch_minDrRecoGen;   //!
   TBranch        *b_bestMatch_drRecoGenK;   //!
   TBranch        *b_bestMatch_maxDrRecoGenFromB;   //!
   TBranch        *b_bestMatch_minDrRecoGenFromB;   //!
   TBranch        *b_bestMatch_drRecoGenFromBK;   //!
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
   TBranch        *b_goodCombB_keta;   //!
   TBranch        *b_goodCombB_ele1pt;   //!
   TBranch        *b_goodCombB_ele2pt;   //!
   TBranch        *b_goodCombB_costhetaSKRad;   //!
   TBranch        *b_goodCombB_costhetaSKCS;   //!
   TBranch        *b_goodCombB_costhetaSKHel;   //!
   TBranch        *b_goodCombB_minpt;   //!
   TBranch        *b_goodCombB_ele1eta;   //!
   TBranch        *b_goodCombB_ele2eta;   //!
   TBranch        *b_goodCombB_ele1pfmva;   //!
   TBranch        *b_goodCombB_ele2pfmva;   //!
   TBranch        *b_goodCombB_ele1lptmva;   //!
   TBranch        *b_goodCombB_ele2lptmva;   //!
   TBranch        *b_goodCombB_causeEle1;   //!
   TBranch        *b_goodCombB_causeEle2;   //!
   TBranch        *b_goodCombB_causeK;   //!
   TBranch        *b_goodCombB_notmatching;   //!
   TBranch        *b_goodCombB_maxDrRecoGen;   //!
   TBranch        *b_goodCombB_minDrRecoGen;   //!
   TBranch        *b_goodCombB_drRecoGenK;   //!
   TBranch        *b_goodCombB_maxDrRecoGenFromB;   //!
   TBranch        *b_goodCombB_minDrRecoGenFromB;   //!
   TBranch        *b_goodCombB_drRecoGenFromBK;   //!
   TBranch        *b_bestSvProbMatch;   //!
   TBranch        *b_bestSvProbMatch_second;   //!
   TBranch        *b_bestSvProbMatchCat0;   //!
   TBranch        *b_bestSvProbMatchCat1;   //!
   TBranch        *b_bestSvProbMatchCat2;   //!
   TBranch        *b_bestSvProbMatchCat0_second;   //!
   TBranch        *b_bestSvProbMatchCat1_second;   //!
   TBranch        *b_bestSvProbMatchCat2_second;   //!
   TBranch        *b_bestSvProbMatchCatNew0;   //!
   TBranch        *b_bestSvProbMatchCatNew1;   //!
   TBranch        *b_bestSvProbMatchCatNew2;   //!
   TBranch        *b_bestSvProbMatchCatNew0_second;   //!
   TBranch        *b_bestSvProbMatchCatNew1_second;   //!
   TBranch        *b_bestSvProbMatchCatNew2_second;   //!
   TBranch        *b_bestSvProbMatch_causeEle1;   //!
   TBranch        *b_bestSvProbMatch_causeEle2;   //!
   TBranch        *b_bestSvProbMatch_notmatching;   //!
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
   TBranch        *b_bestXYsigMatch_second;   //!
   TBranch        *b_bestXYsigMatchCat0;   //!
   TBranch        *b_bestXYsigMatchCat1;   //!
   TBranch        *b_bestXYsigMatchCat2;   //!
   TBranch        *b_bestXYsigMatchCat0_second;   //!
   TBranch        *b_bestXYsigMatchCat1_second;   //!
   TBranch        *b_bestXYsigMatchCat2_second;   //!
   TBranch        *b_bestXYsigMatchCatNew0;   //!
   TBranch        *b_bestXYsigMatchCatNew1;   //!
   TBranch        *b_bestXYsigMatchCatNew2;   //!
   TBranch        *b_bestXYsigMatchCatNew0_second;   //!
   TBranch        *b_bestXYsigMatchCatNew1_second;   //!
   TBranch        *b_bestXYsigMatchCatNew2_second;   //!
   TBranch        *b_bestXYsigMatch_causeEle1;   //!
   TBranch        *b_bestXYsigMatch_causeEle2;   //!
   TBranch        *b_bestXYsigMatch_notmatching;   //!
   TBranch        *b_bestXYsigMatch_causeK;   //!
   TBranch        *b_bestXYsigMatch_notok_ele1pt;   //!
   TBranch        *b_bestXYsigMatch_notok_ele2pt;   //!
   TBranch        *b_bestXYsigMatch_notok_kpt;   //!
   TBranch        *b_bestXYsigMatch_notok_minpt;   //!
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
   TBranch        *b_bestXYsigMatch_ok_minpt;   //!
   TBranch        *b_bestXYsigMatch_ok_ele1eta;   //!
   TBranch        *b_bestXYsigMatch_ok_ele2eta;   //!
   TBranch        *b_bestXYsigMatch_ok_keta;   //!
   TBranch        *b_bestXYsigMatch_ok_pfmva1;   //!
   TBranch        *b_bestXYsigMatch_ok_pfmva2;   //!
   TBranch        *b_bestXYsigMatch_ok_lptmva1;   //!
   TBranch        *b_bestXYsigMatch_ok_lptmva2;   //!
   TBranch        *b_numberBetterSvProbTriplets;   //!
   TBranch        *b_numberBetterXYsigTriplets;   //!
   TBranch        *b_numberBetterCos2DTriplets;   //!
   TBranch        *b_numberBetterAllPtSumTriplets;   //!
   TBranch        *b_numberBetterBPtTriplets;   //!
   TBranch        *b_bestCos2DMatch;   //!
   TBranch        *b_bestCos2DMatch_second;   //!
   TBranch        *b_bestCos2DMatchCat0;   //!
   TBranch        *b_bestCos2DMatchCat1;   //!
   TBranch        *b_bestCos2DMatchCat2;   //!
   TBranch        *b_bestCos2DMatchCat0_second;   //!
   TBranch        *b_bestCos2DMatchCat1_second;   //!
   TBranch        *b_bestCos2DMatchCat2_second;   //!
   TBranch        *b_bestCos2DMatchCatNew0;   //!
   TBranch        *b_bestCos2DMatchCatNew1;   //!
   TBranch        *b_bestCos2DMatchCatNew2;   //!
   TBranch        *b_bestCos2DMatchCatNew0_second;   //!
   TBranch        *b_bestCos2DMatchCatNew1_second;   //!
   TBranch        *b_bestCos2DMatchCatNew2_second;   //!
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
   TBranch        *b_bestAllPtSumMatch;   //!
   TBranch        *b_bestAllPtSumMatch_second;   //!
   TBranch        *b_bestAllPtSumMatchCat0;   //!
   TBranch        *b_bestAllPtSumMatchCat1;   //!
   TBranch        *b_bestAllPtSumMatchCat2;   //!
   TBranch        *b_bestAllPtSumMatchCatNew0;   //!
   TBranch        *b_bestAllPtSumMatchCatNew1;   //!
   TBranch        *b_bestAllPtSumMatchCatNew2;   //!
   TBranch        *b_bestAllPtSumMatchCat0_second;   //!
   TBranch        *b_bestAllPtSumMatchCat1_second;   //!
   TBranch        *b_bestAllPtSumMatchCat2_second;   //!
   TBranch        *b_bestAllPtSumMatchCatNew0_second;   //!
   TBranch        *b_bestAllPtSumMatchCatNew1_second;   //!
   TBranch        *b_bestAllPtSumMatchCatNew2_second;   //!
   TBranch        *b_bestAllPtSumMatch_causeEle1;   //!
   TBranch        *b_bestAllPtSumMatch_causeEle2;   //!
   TBranch        *b_bestAllPtSumMatch_causeK;   //!
   TBranch        *b_bestBPtMatch;   //!
   TBranch        *b_bestBPtMatch_second;   //!
   TBranch        *b_bestBPtMatchCat0;   //!
   TBranch        *b_bestBPtMatchCat1;   //!
   TBranch        *b_bestBPtMatchCat2;   //!
   TBranch        *b_bestBPtMatchCatNew0;   //!
   TBranch        *b_bestBPtMatchCatNew1;   //!
   TBranch        *b_bestBPtMatchCatNew2;   //!
   TBranch        *b_bestBPtMatchCat0_second;   //!
   TBranch        *b_bestBPtMatchCat1_second;   //!
   TBranch        *b_bestBPtMatchCat2_second;   //!
   TBranch        *b_bestBPtMatchCatNew0_second;   //!
   TBranch        *b_bestBPtMatchCatNew1_second;   //!
   TBranch        *b_bestBPtMatchCatNew2_second;   //!
   TBranch        *b_bestBPtMatch_causeEle1;   //!
   TBranch        *b_bestBPtMatch_causeEle2;   //!
   TBranch        *b_bestBPtMatch_causeK;   //!

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
    //TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../data/Jan16prod/BuToKee_Toee_BParkNANO_mc_2020Jan16.txt.root");
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../data/Jan16prod/BuToKJpsi_Toee_BParkNANO_mc_2020Jan16.txt.root");
    if (!f || !f->IsOpen()) {
      //f = new TFile("../data/Jan16prod/BuToKee_Toee_BParkNANO_mc_2020Jan16.txt.root");
      f = new TFile("../data/Jan16prod/BuToKJpsi_Toee_BParkNANO_mc_2020Jan16.txt.root");
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
   goodCombB_keta = 0;
   goodCombB_ele1pt = 0;
   goodCombB_ele2pt = 0;
   goodCombB_costhetaSKRad = 0;
   goodCombB_costhetaSKCS = 0;
   goodCombB_costhetaSKHel = 0;
   goodCombB_minpt = 0;
   goodCombB_ele1eta = 0;
   goodCombB_ele2eta = 0;
   goodCombB_ele1pfmva = 0;
   goodCombB_ele2pfmva = 0;
   goodCombB_ele1lptmva = 0;
   goodCombB_ele2lptmva = 0;
   goodCombB_causeEle1 = 0;
   goodCombB_causeEle2 = 0;
   goodCombB_causeK = 0;
   goodCombB_notmatching = 0;
   goodCombB_maxDrRecoGen = 0;
   goodCombB_minDrRecoGen = 0;
   goodCombB_drRecoGenK = 0;
   goodCombB_maxDrRecoGenFromB = 0;
   goodCombB_minDrRecoGenFromB = 0;
   goodCombB_drRecoGenFromBK = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("debug_svprob", &debug_svprob, &b_debug_svprob);
   fChain->SetBranchAddress("debug_svprob_match", &debug_svprob_match, &b_debug_svprob_match);
   fChain->SetBranchAddress("debug_pf_svprob", &debug_pf_svprob, &b_debug_pf_svprob);
   fChain->SetBranchAddress("debug_pf_svprob_match", &debug_pf_svprob_match, &b_debug_pf_svprob_match);
   fChain->SetBranchAddress("iHLT_Mu12_IP6", &iHLT_Mu12_IP6, &b_iHLT_Mu12_IP6);
   fChain->SetBranchAddress("iHLT_Mu9_IP6", &iHLT_Mu9_IP6, &b_iHLT_Mu9_IP6);
   fChain->SetBranchAddress("theEvent", &theEvent, &b_theEvent);
   fChain->SetBranchAddress("rho", &rho, &b_rho);
   fChain->SetBranchAddress("rightMcTruth", &rightMcTruth, &b_rightMcTruth);
   fChain->SetBranchAddress("gen_costhetaSKRad", &gen_costhetaSKRad, &b_gen_costhetaSKRad);
   fChain->SetBranchAddress("gen_costhetaSKCS",  &gen_costhetaSKCS,  &b_gen_costhetaSKCS);
   fChain->SetBranchAddress("gen_costhetaSKHel", &gen_costhetaSKHel, &b_gen_costhetaSKHel);
   fChain->SetBranchAddress("goodBSize", &goodBSize, &b_goodBSize);
   fChain->SetBranchAddress("goodTrueBSize", &goodTrueBSize, &b_goodTrueBSize);
   fChain->SetBranchAddress("goodCombBSize", &goodCombBSize, &b_goodCombBSize);
   fChain->SetBranchAddress("handMadeB", &handMadeB, &b_handMadeB);
   fChain->SetBranchAddress("handMadeBmass", &handMadeBmass, &b_handMadeBmass);
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
   fChain->SetBranchAddress("bestMatch_Bmass", &bestMatch_Bmass, &b_bestMatch_Bmass);
   fChain->SetBranchAddress("bestMatch_SvProb", &bestMatch_SvProb, &b_bestMatch_SvProb);
   fChain->SetBranchAddress("bestMatch_XYSig", &bestMatch_XYSig, &b_bestMatch_XYSig);
   fChain->SetBranchAddress("bestMatch_Cos2D", &bestMatch_Cos2D, &b_bestMatch_Cos2D);
   fChain->SetBranchAddress("bestMatch_PtSum", &bestMatch_PtSum, &b_bestMatch_PtSum);
   fChain->SetBranchAddress("bestMatch_KPt", &bestMatch_KPt, &b_bestMatch_KPt);
   fChain->SetBranchAddress("bestMatch_KEta", &bestMatch_KEta, &b_bestMatch_KEta);
   fChain->SetBranchAddress("bestMatch_Ele1Pt", &bestMatch_Ele1Pt, &b_bestMatch_Ele1Pt);
   fChain->SetBranchAddress("bestMatch_Ele2Pt", &bestMatch_Ele2Pt, &b_bestMatch_Ele2Pt);
   fChain->SetBranchAddress("bestMatch_costhetaSKRad", &bestMatch_costhetaSKRad, &b_bestMatch_costhetaSKRad);
   fChain->SetBranchAddress("bestMatch_costhetaSKCS", &bestMatch_costhetaSKCS, &b_bestMatch_costhetaSKCS);
   fChain->SetBranchAddress("bestMatch_costhetaSKHel", &bestMatch_costhetaSKHel, &b_bestMatch_costhetaSKHel);
   fChain->SetBranchAddress("bestMatch_MinPt", &bestMatch_MinPt, &b_bestMatch_MinPt);
   fChain->SetBranchAddress("bestMatch_Ele1Eta", &bestMatch_Ele1Eta, &b_bestMatch_Ele1Eta);
   fChain->SetBranchAddress("bestMatch_Ele2Eta", &bestMatch_Ele2Eta, &b_bestMatch_Ele2Eta);
   fChain->SetBranchAddress("bestMatch_Ele1pfmva", &bestMatch_Ele1pfmva, &b_bestMatch_Ele1pfmva);
   fChain->SetBranchAddress("bestMatch_Ele2pfmva", &bestMatch_Ele2pfmva, &b_bestMatch_Ele2pfmva);
   fChain->SetBranchAddress("bestMatch_Ele1lptmva", &bestMatch_Ele1lptmva, &b_bestMatch_Ele1lptmva);
   fChain->SetBranchAddress("bestMatch_Ele2lptmva", &bestMatch_Ele2lptmva, &b_bestMatch_Ele2lptmva);
   fChain->SetBranchAddress("bestMatch_maxDrRecoGen", &bestMatch_maxDrRecoGen, &b_bestMatch_maxDrRecoGen);
   fChain->SetBranchAddress("bestMatch_minDrRecoGen", &bestMatch_minDrRecoGen, &b_bestMatch_minDrRecoGen);
   fChain->SetBranchAddress("bestMatch_drRecoGenK", &bestMatch_drRecoGenK, &b_bestMatch_drRecoGenK);
   fChain->SetBranchAddress("bestMatch_maxDrRecoGenFromB", &bestMatch_maxDrRecoGenFromB, &b_bestMatch_maxDrRecoGenFromB);
   fChain->SetBranchAddress("bestMatch_minDrRecoGenFromB", &bestMatch_minDrRecoGenFromB, &b_bestMatch_minDrRecoGenFromB);
   fChain->SetBranchAddress("bestMatch_drRecoGenFromBK", &bestMatch_drRecoGenFromBK, &b_bestMatch_drRecoGenFromBK);
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
   fChain->SetBranchAddress("goodCombB_keta", &goodCombB_keta, &b_goodCombB_keta);
   fChain->SetBranchAddress("goodCombB_ele1pt", &goodCombB_ele1pt, &b_goodCombB_ele1pt);
   fChain->SetBranchAddress("goodCombB_ele2pt", &goodCombB_ele2pt, &b_goodCombB_ele2pt);
   fChain->SetBranchAddress("goodCombB_costhetaSKRad", &goodCombB_costhetaSKRad, &b_goodCombB_costhetaSKRad);
   fChain->SetBranchAddress("goodCombB_costhetaSKCS", &goodCombB_costhetaSKCS, &b_goodCombB_costhetaSKCS);
   fChain->SetBranchAddress("goodCombB_costhetaSKHel", &goodCombB_costhetaSKHel, &b_goodCombB_costhetaSKHel);
   fChain->SetBranchAddress("goodCombB_minpt", &goodCombB_minpt, &b_goodCombB_minpt);
   fChain->SetBranchAddress("goodCombB_ele1eta", &goodCombB_ele1eta, &b_goodCombB_ele1eta);
   fChain->SetBranchAddress("goodCombB_ele2eta", &goodCombB_ele2eta, &b_goodCombB_ele2eta);
   fChain->SetBranchAddress("goodCombB_ele1pfmva", &goodCombB_ele1pfmva, &b_goodCombB_ele1pfmva);
   fChain->SetBranchAddress("goodCombB_ele2pfmva", &goodCombB_ele2pfmva, &b_goodCombB_ele2pfmva);
   fChain->SetBranchAddress("goodCombB_ele1lptmva", &goodCombB_ele1lptmva, &b_goodCombB_ele1lptmva);
   fChain->SetBranchAddress("goodCombB_ele2lptmva", &goodCombB_ele2lptmva, &b_goodCombB_ele2lptmva);
   fChain->SetBranchAddress("goodCombB_causeEle1", &goodCombB_causeEle1, &b_goodCombB_causeEle1);
   fChain->SetBranchAddress("goodCombB_causeEle2", &goodCombB_causeEle2, &b_goodCombB_causeEle2);
   fChain->SetBranchAddress("goodCombB_causeK", &goodCombB_causeK, &b_goodCombB_causeK);
   fChain->SetBranchAddress("goodCombB_notmatching", &goodCombB_notmatching, &b_goodCombB_notmatching);
   fChain->SetBranchAddress("goodCombB_maxDrRecoGen", &goodCombB_maxDrRecoGen, &b_goodCombB_maxDrRecoGen);
   fChain->SetBranchAddress("goodCombB_minDrRecoGen", &goodCombB_minDrRecoGen, &b_goodCombB_minDrRecoGen);
   fChain->SetBranchAddress("goodCombB_drRecoGenK", &goodCombB_drRecoGenK, &b_goodCombB_drRecoGenK);
   fChain->SetBranchAddress("goodCombB_maxDrRecoGenFromB", &goodCombB_maxDrRecoGenFromB, &b_goodCombB_maxDrRecoGenFromB);
   fChain->SetBranchAddress("goodCombB_minDrRecoGenFromB", &goodCombB_minDrRecoGenFromB, &b_goodCombB_minDrRecoGenFromB);
   fChain->SetBranchAddress("goodCombB_drRecoGenFromBK", &goodCombB_drRecoGenFromBK, &b_goodCombB_drRecoGenFromBK);
   fChain->SetBranchAddress("bestSvProbMatch", &bestSvProbMatch, &b_bestSvProbMatch);
   fChain->SetBranchAddress("bestSvProbMatch_second", &bestSvProbMatch_second, &b_bestSvProbMatch_second);
   fChain->SetBranchAddress("bestSvProbMatchCat0", &bestSvProbMatchCat0, &b_bestSvProbMatchCat0);
   fChain->SetBranchAddress("bestSvProbMatchCat1", &bestSvProbMatchCat1, &b_bestSvProbMatchCat1);
   fChain->SetBranchAddress("bestSvProbMatchCat2", &bestSvProbMatchCat2, &b_bestSvProbMatchCat2);
   fChain->SetBranchAddress("bestSvProbMatchCat0_second", &bestSvProbMatchCat0_second, &b_bestSvProbMatchCat0_second);
   fChain->SetBranchAddress("bestSvProbMatchCat1_second", &bestSvProbMatchCat1_second, &b_bestSvProbMatchCat1_second);
   fChain->SetBranchAddress("bestSvProbMatchCat2_second", &bestSvProbMatchCat2_second, &b_bestSvProbMatchCat2_second);
   fChain->SetBranchAddress("bestSvProbMatchCatNew0", &bestSvProbMatchCatNew0, &b_bestSvProbMatchCatNew0);
   fChain->SetBranchAddress("bestSvProbMatchCatNew1", &bestSvProbMatchCatNew1, &b_bestSvProbMatchCatNew1);
   fChain->SetBranchAddress("bestSvProbMatchCatNew2", &bestSvProbMatchCatNew2, &b_bestSvProbMatchCatNew2);
   fChain->SetBranchAddress("bestSvProbMatchCatNew0_second", &bestSvProbMatchCatNew0_second, &b_bestSvProbMatchCatNew0_second);
   fChain->SetBranchAddress("bestSvProbMatchCatNew1_second", &bestSvProbMatchCatNew1_second, &b_bestSvProbMatchCatNew1_second);
   fChain->SetBranchAddress("bestSvProbMatchCatNew2_second", &bestSvProbMatchCatNew2_second, &b_bestSvProbMatchCatNew2_second);
   fChain->SetBranchAddress("bestSvProbMatch_causeEle1", &bestSvProbMatch_causeEle1, &b_bestSvProbMatch_causeEle1);
   fChain->SetBranchAddress("bestSvProbMatch_causeEle2", &bestSvProbMatch_causeEle2, &b_bestSvProbMatch_causeEle2);
   fChain->SetBranchAddress("bestSvProbMatch_notmatching", &bestSvProbMatch_notmatching, &b_bestSvProbMatch_notmatching);
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
   fChain->SetBranchAddress("bestXYsigMatch_second", &bestXYsigMatch_second, &b_bestXYsigMatch_second);
   fChain->SetBranchAddress("bestXYsigMatchCat0", &bestXYsigMatchCat0, &b_bestXYsigMatchCat0);
   fChain->SetBranchAddress("bestXYsigMatchCat1", &bestXYsigMatchCat1, &b_bestXYsigMatchCat1);
   fChain->SetBranchAddress("bestXYsigMatchCat2", &bestXYsigMatchCat2, &b_bestXYsigMatchCat2);
   fChain->SetBranchAddress("bestXYsigMatchCat0_second", &bestXYsigMatchCat0_second, &b_bestXYsigMatchCat0_second);
   fChain->SetBranchAddress("bestXYsigMatchCat1_second", &bestXYsigMatchCat1_second, &b_bestXYsigMatchCat1_second);
   fChain->SetBranchAddress("bestXYsigMatchCat2_second", &bestXYsigMatchCat2_second, &b_bestXYsigMatchCat2_second);
   fChain->SetBranchAddress("bestXYsigMatchCatNew0", &bestXYsigMatchCatNew0, &b_bestXYsigMatchCatNew0);
   fChain->SetBranchAddress("bestXYsigMatchCatNew1", &bestXYsigMatchCatNew1, &b_bestXYsigMatchCatNew1);
   fChain->SetBranchAddress("bestXYsigMatchCatNew2", &bestXYsigMatchCatNew2, &b_bestXYsigMatchCatNew2);
   fChain->SetBranchAddress("bestXYsigMatchCatNew0_second", &bestXYsigMatchCatNew0_second, &b_bestXYsigMatchCatNew0_second);
   fChain->SetBranchAddress("bestXYsigMatchCatNew1_second", &bestXYsigMatchCatNew1_second, &b_bestXYsigMatchCatNew1_second);
   fChain->SetBranchAddress("bestXYsigMatchCatNew2_second", &bestXYsigMatchCatNew2_second, &b_bestXYsigMatchCatNew2_second);
   fChain->SetBranchAddress("bestXYsigMatch_causeEle1", &bestXYsigMatch_causeEle1, &b_bestXYsigMatch_causeEle1);
   fChain->SetBranchAddress("bestXYsigMatch_causeEle2", &bestXYsigMatch_causeEle2, &b_bestXYsigMatch_causeEle2);
   fChain->SetBranchAddress("bestXYsigMatch_notmatching", &bestXYsigMatch_notmatching, &b_bestXYsigMatch_notmatching);
   fChain->SetBranchAddress("bestXYsigMatch_causeK", &bestXYsigMatch_causeK, &b_bestXYsigMatch_causeK);
   fChain->SetBranchAddress("bestXYsigMatch_notok_ele1pt", &bestXYsigMatch_notok_ele1pt, &b_bestXYsigMatch_notok_ele1pt);
   fChain->SetBranchAddress("bestXYsigMatch_notok_ele2pt", &bestXYsigMatch_notok_ele2pt, &b_bestXYsigMatch_notok_ele2pt);
   fChain->SetBranchAddress("bestXYsigMatch_notok_kpt", &bestXYsigMatch_notok_kpt, &b_bestXYsigMatch_notok_kpt);
   fChain->SetBranchAddress("bestXYsigMatch_notok_minpt", &bestXYsigMatch_notok_minpt, &b_bestXYsigMatch_notok_minpt);
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
   fChain->SetBranchAddress("bestXYsigMatch_ok_minpt", &bestXYsigMatch_ok_minpt, &b_bestXYsigMatch_ok_minpt);
   fChain->SetBranchAddress("bestXYsigMatch_ok_ele1eta", &bestXYsigMatch_ok_ele1eta, &b_bestXYsigMatch_ok_ele1eta);
   fChain->SetBranchAddress("bestXYsigMatch_ok_ele2eta", &bestXYsigMatch_ok_ele2eta, &b_bestXYsigMatch_ok_ele2eta);
   fChain->SetBranchAddress("bestXYsigMatch_ok_keta", &bestXYsigMatch_ok_keta, &b_bestXYsigMatch_ok_keta);
   fChain->SetBranchAddress("bestXYsigMatch_ok_pfmva1", &bestXYsigMatch_ok_pfmva1, &b_bestXYsigMatch_ok_pfmva1);
   fChain->SetBranchAddress("bestXYsigMatch_ok_pfmva2", &bestXYsigMatch_ok_pfmva2, &b_bestXYsigMatch_ok_pfmva2);
   fChain->SetBranchAddress("bestXYsigMatch_ok_lptmva1", &bestXYsigMatch_ok_lptmva1, &b_bestXYsigMatch_ok_lptmva1);
   fChain->SetBranchAddress("bestXYsigMatch_ok_lptmva2", &bestXYsigMatch_ok_lptmva2, &b_bestXYsigMatch_ok_lptmva2);
   fChain->SetBranchAddress("numberBetterSvProbTriplets", &numberBetterSvProbTriplets, &b_numberBetterSvProbTriplets);
   fChain->SetBranchAddress("numberBetterXYsigTriplets", &numberBetterXYsigTriplets, &b_numberBetterXYsigTriplets);
   fChain->SetBranchAddress("numberBetterCos2DTriplets", &numberBetterCos2DTriplets, &b_numberBetterCos2DTriplets);
   fChain->SetBranchAddress("numberBetterAllPtSumTriplets", &numberBetterAllPtSumTriplets, &b_numberBetterAllPtSumTriplets);
   fChain->SetBranchAddress("numberBetterBPtTriplets", &numberBetterBPtTriplets, &b_numberBetterBPtTriplets);
   fChain->SetBranchAddress("bestCos2DMatch", &bestCos2DMatch, &b_bestCos2DMatch);
   fChain->SetBranchAddress("bestCos2DMatch_second", &bestCos2DMatch_second, &b_bestCos2DMatch_second);
   fChain->SetBranchAddress("bestCos2DMatchCat0", &bestCos2DMatchCat0, &b_bestCos2DMatchCat0);
   fChain->SetBranchAddress("bestCos2DMatchCat1", &bestCos2DMatchCat1, &b_bestCos2DMatchCat1);
   fChain->SetBranchAddress("bestCos2DMatchCat2", &bestCos2DMatchCat2, &b_bestCos2DMatchCat2);
   fChain->SetBranchAddress("bestCos2DMatchCat0_second", &bestCos2DMatchCat0_second, &b_bestCos2DMatchCat0_second);
   fChain->SetBranchAddress("bestCos2DMatchCat1_second", &bestCos2DMatchCat1_second, &b_bestCos2DMatchCat1_second);
   fChain->SetBranchAddress("bestCos2DMatchCat2_second", &bestCos2DMatchCat2_second, &b_bestCos2DMatchCat2_second);
   fChain->SetBranchAddress("bestCos2DMatchCatNew0", &bestCos2DMatchCatNew0, &b_bestCos2DMatchCatNew0);
   fChain->SetBranchAddress("bestCos2DMatchCatNew1", &bestCos2DMatchCatNew1, &b_bestCos2DMatchCatNew1);
   fChain->SetBranchAddress("bestCos2DMatchCatNew2", &bestCos2DMatchCatNew2, &b_bestCos2DMatchCatNew2);
   fChain->SetBranchAddress("bestCos2DMatchCatNew0_second", &bestCos2DMatchCatNew0_second, &b_bestCos2DMatchCatNew0_second);
   fChain->SetBranchAddress("bestCos2DMatchCatNew1_second", &bestCos2DMatchCatNew1_second, &b_bestCos2DMatchCatNew1_second);
   fChain->SetBranchAddress("bestCos2DMatchCatNew2_second", &bestCos2DMatchCatNew2_second, &b_bestCos2DMatchCatNew2_second);
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
   fChain->SetBranchAddress("bestAllPtSumMatch", &bestAllPtSumMatch, &b_bestAllPtSumMatch);
   fChain->SetBranchAddress("bestAllPtSumMatch_second", &bestAllPtSumMatch_second, &b_bestAllPtSumMatch_second);
   fChain->SetBranchAddress("bestAllPtSumMatchCat0", &bestAllPtSumMatchCat0, &b_bestAllPtSumMatchCat0);
   fChain->SetBranchAddress("bestAllPtSumMatchCat1", &bestAllPtSumMatchCat1, &b_bestAllPtSumMatchCat1);
   fChain->SetBranchAddress("bestAllPtSumMatchCat2", &bestAllPtSumMatchCat2, &b_bestAllPtSumMatchCat2);
   fChain->SetBranchAddress("bestAllPtSumMatchCatNew0", &bestAllPtSumMatchCatNew0, &b_bestAllPtSumMatchCatNew0);
   fChain->SetBranchAddress("bestAllPtSumMatchCatNew1", &bestAllPtSumMatchCatNew1, &b_bestAllPtSumMatchCatNew1);
   fChain->SetBranchAddress("bestAllPtSumMatchCatNew2", &bestAllPtSumMatchCatNew2, &b_bestAllPtSumMatchCatNew2);
   fChain->SetBranchAddress("bestAllPtSumMatchCat0_second", &bestAllPtSumMatchCat0_second, &b_bestAllPtSumMatchCat0_second);
   fChain->SetBranchAddress("bestAllPtSumMatchCat1_second", &bestAllPtSumMatchCat1_second, &b_bestAllPtSumMatchCat1_second);
   fChain->SetBranchAddress("bestAllPtSumMatchCat2_second", &bestAllPtSumMatchCat2_second, &b_bestAllPtSumMatchCat2_second);
   fChain->SetBranchAddress("bestAllPtSumMatchCatNew0_second", &bestAllPtSumMatchCatNew0_second, &b_bestAllPtSumMatchCatNew0_second);
   fChain->SetBranchAddress("bestAllPtSumMatchCatNew1_second", &bestAllPtSumMatchCatNew1_second, &b_bestAllPtSumMatchCatNew1_second);
   fChain->SetBranchAddress("bestAllPtSumMatchCatNew2_second", &bestAllPtSumMatchCatNew2_second, &b_bestAllPtSumMatchCatNew2_second);
   fChain->SetBranchAddress("bestAllPtSumMatch_causeEle1", &bestAllPtSumMatch_causeEle1, &b_bestAllPtSumMatch_causeEle1);
   fChain->SetBranchAddress("bestAllPtSumMatch_causeEle2", &bestAllPtSumMatch_causeEle2, &b_bestAllPtSumMatch_causeEle2);
   fChain->SetBranchAddress("bestAllPtSumMatch_causeK", &bestAllPtSumMatch_causeK, &b_bestAllPtSumMatch_causeK);
   fChain->SetBranchAddress("bestBPtMatch", &bestBPtMatch, &b_bestBPtMatch);
   fChain->SetBranchAddress("bestBPtMatch_second", &bestBPtMatch_second, &b_bestBPtMatch_second);
   fChain->SetBranchAddress("bestBPtMatchCat0", &bestBPtMatchCat0, &b_bestBPtMatchCat0);
   fChain->SetBranchAddress("bestBPtMatchCat1", &bestBPtMatchCat1, &b_bestBPtMatchCat1);
   fChain->SetBranchAddress("bestBPtMatchCat2", &bestBPtMatchCat2, &b_bestBPtMatchCat2);
   fChain->SetBranchAddress("bestBPtMatchCatNew0", &bestBPtMatchCatNew0, &b_bestBPtMatchCatNew0);
   fChain->SetBranchAddress("bestBPtMatchCatNew1", &bestBPtMatchCatNew1, &b_bestBPtMatchCatNew1);
   fChain->SetBranchAddress("bestBPtMatchCatNew2", &bestBPtMatchCatNew2, &b_bestBPtMatchCatNew2);
   fChain->SetBranchAddress("bestBPtMatchCat0_second", &bestBPtMatchCat0_second, &b_bestBPtMatchCat0_second);
   fChain->SetBranchAddress("bestBPtMatchCat1_second", &bestBPtMatchCat1_second, &b_bestBPtMatchCat1_second);
   fChain->SetBranchAddress("bestBPtMatchCat2_second", &bestBPtMatchCat2_second, &b_bestBPtMatchCat2_second);
   fChain->SetBranchAddress("bestBPtMatchCatNew0_second", &bestBPtMatchCatNew0_second, &b_bestBPtMatchCatNew0_second);
   fChain->SetBranchAddress("bestBPtMatchCatNew1_second", &bestBPtMatchCatNew1_second, &b_bestBPtMatchCatNew1_second);
   fChain->SetBranchAddress("bestBPtMatchCatNew2_second", &bestBPtMatchCatNew2_second, &b_bestBPtMatchCatNew2_second);
   fChain->SetBranchAddress("bestBPtMatch_causeEle1", &bestBPtMatch_causeEle1, &b_bestBPtMatch_causeEle1);
   fChain->SetBranchAddress("bestBPtMatch_causeEle2", &bestBPtMatch_causeEle2, &b_bestBPtMatch_causeEle2);
   fChain->SetBranchAddress("bestBPtMatch_causeK", &bestBPtMatch_causeK, &b_bestBPtMatch_causeK);
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
