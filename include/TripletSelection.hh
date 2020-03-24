#ifndef TripletSelection_h
#define TripletSelection_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include <vector>
#include <string>
#include <TLorentzVector.h>

#include "./BParkBase.h"

using namespace std;

class TripletSelection : public BParkBase{
public:

  //! constructor
  TripletSelection(TTree *tree=0);
  //! destructor
  virtual ~TripletSelection();
  //! loop over events
  void Loop();
  void PrepareOutputs(std::string filename);             
  
private:
  
  // Analysis methods
  bool isMcB( int myB );
  float dRgen( int myB ); 
  float dRRecoGenEle( int theRecoEle );
  float dRRecoGenK( int theRecoK );
  float cosThetaL( int myB );
  float cosThetaStarK( int myB );
  float cosThetaStarKCS( int myB );
  float cosThetaStarKGen( int myB );
  void bookOutputTree();
  void bookOutputHistos();

  // ---- outputs
  TFile* outFile_;
  TTree* outTree_;
  TH1F* h_entries;
  TH1F* h_selection;

  // dataset name
  std::string _datasetName;      

  //---output tree branches variables
  vector <float> debug_svprob={};
  vector <float> debug_pf_svprob={};
  vector <int> debug_svprob_match={};
  vector <int> debug_pf_svprob_match={};

  int theEvent;
  float  rho;
  //
  int iHLT_Mu12_IP6;
  int iHLT_Mu9_IP6;
  // 
  int goodBSize;
  int goodTrueBSize;
  int goodCombBSize;
  //
  vector <float> goodTrueB_maxMinDREle={};   
  vector <float> goodTrueB_maxMinDREle_dEta={};   
  vector <float> goodTrueB_maxMinDREle_dPhi={};   
  vector <float> goodTrueB_maxMinDREle_dPtOverPt={};   
  vector <float> goodTrueB_dRgen={};   
  vector <float> goodTrueB_maxDRTrack={};   
  //
  vector <float> goodCombB_maxMinDREle={};   
  vector <float> goodCombB_maxMinDREle_dEta={};   
  vector <float> goodCombB_maxMinDREle_dPhi={};   
  vector <float> goodCombB_maxMinDREle_dPtOverPt={};   
  vector <float> goodCombB_maxDRTrack={};   
  //
  int goodTrueBs_SvProbMatch;
  int goodTrueBs_xySigMatch;
  int goodTrueBs_cos2DMatch;
  int goodTrueBs_ptSumMatch;
  int goodTrueBs_kptMatch;
  //
  float goodTrueBs_bestSvProb_dRmax;
  float goodTrueBs_bestSvProb_dRmin;
  float goodTrueBs_bestXYSig_dRmax;
  float goodTrueBs_bestXYSig_dRmin;
  float goodTrueBs_bestCos2d_dRmax;
  float goodTrueBs_bestCos2d_dRmin;
  float goodTrueBs_bestPtsum_dRmax;
  float goodTrueBs_bestPtsum_dRmin;
  float goodTrueBs_bestKpt_dRmax;
  float goodTrueBs_bestKpt_dRmin;
  //
  float bestMatch_SvProb;
  float bestMatch_XYSig;
  float bestMatch_Cos2D;
  float bestMatch_PtSum;
  float bestMatch_KPt;
  float bestMatch_KEta;
  float bestMatch_Ele1Pt;
  float bestMatch_Ele2Pt;
  float bestMatch_MinPt;
  float bestMatch_Ele1Eta;
  float bestMatch_Ele2Eta;
  float bestMatch_Ele1pfmva;
  float bestMatch_Ele2pfmva;
  float bestMatch_Ele1lptmva;
  float bestMatch_Ele2lptmva;
  float bestMatch_maxDrRecoGen;
  float bestMatch_minDrRecoGen;
  float bestMatch_drRecoGenK;
  //            
  vector<float> goodTrueB_svProb_notBestMatch={};
  vector<float> goodTrueB_xySig_notBestMatch={};
  vector<float> goodTrueB_cos2D_notBestMatch={};
  vector<float> goodTrueB_ptsum_notBestMatch={};
  vector<float> goodTrueB_kpt_notBestMatch={};
  //
  vector<float> goodCombB_svProb={};
  vector<float> goodCombB_xySig={};
  vector<float> goodCombB_cos2D={};
  vector<float> goodCombB_ptsum={};
  vector<float> goodCombB_kpt={};
  vector<float> goodCombB_keta={};
  vector<float> goodCombB_ele1pt={};
  vector<float> goodCombB_ele2pt={};
  vector<float> goodCombB_minpt={};
  vector<float> goodCombB_ele1eta={};
  vector<float> goodCombB_ele2eta={};
  vector<float> goodCombB_ele1pfmva={};
  vector<float> goodCombB_ele2pfmva={};
  vector<float> goodCombB_ele1lptmva={};
  vector<float> goodCombB_ele2lptmva={};
  vector<float> goodCombB_maxDrRecoGen={};
  vector<float> goodCombB_minDrRecoGen={};
  vector<float> goodCombB_drRecoGenK={};
  //
  vector<float> combToBestB_svProb={};
  vector<float> combToBestB_xySig={};
  vector<float> combToBestB_cos2D={};
  vector<float> combToBestB_ptsum={};
  vector<float> combToBestB_kpt={};
  vector<float> combToBestB_maxDrRecoGen={};
  vector<float> combToBestB_minDrRecoGen={};
  vector<float> combToBestB_drRecoGenK={};
  //
  int bestSvProbMatch;
  int bestSvProbMatchCat0;
  int bestSvProbMatchCat1;
  int bestSvProbMatchCat2;
  int bestSvProbMatch_causeEle1;
  int bestSvProbMatch_causeEle2;
  int bestSvProbMatch_causeK;
  float bestSvProbMatch_notok_ele1pt;
  float bestSvProbMatch_notok_ele2pt;
  float bestSvProbMatch_notok_kpt;
  float bestSvProbMatch_notok_ele1eta;
  float bestSvProbMatch_notok_ele2eta;
  float bestSvProbMatch_notok_keta;
  float bestSvProbMatch_ok_ele1pt;
  float bestSvProbMatch_ok_ele2pt;
  float bestSvProbMatch_ok_kpt;
  float bestSvProbMatch_ok_ele1eta;
  float bestSvProbMatch_ok_ele2eta;
  float bestSvProbMatch_ok_keta;

  int bestXYsigMatch;
  int bestXYsigMatchCat0;
  int bestXYsigMatchCat1;
  int bestXYsigMatchCat2;
  int bestXYsigMatch_causeEle1;
  int bestXYsigMatch_causeEle2;
  int bestXYsigMatch_causeK;
  float bestXYsigMatch_notok_ele1pt;
  float bestXYsigMatch_notok_ele2pt;
  float bestXYsigMatch_notok_kpt;
  float bestXYsigMatch_notok_minpt;
  float bestXYsigMatch_notok_ele1eta;
  float bestXYsigMatch_notok_ele2eta;
  float bestXYsigMatch_notok_keta;
  float bestXYsigMatch_notok_pfmva1;
  float bestXYsigMatch_notok_pfmva2;
  float bestXYsigMatch_notok_lptmva1;
  float bestXYsigMatch_notok_lptmva2;
  float bestXYsigMatch_notok_costhetaSK;
  float bestXYsigMatch_notok_costhetaSKCS;
  float bestXYsigMatch_notok_costhetaL;
  float bestXYsigMatch_ok_ele1pt;
  float bestXYsigMatch_ok_ele2pt;
  float bestXYsigMatch_ok_kpt;
  float bestXYsigMatch_ok_minpt;
  float bestXYsigMatch_ok_ele1eta;
  float bestXYsigMatch_ok_ele2eta;
  float bestXYsigMatch_ok_keta;
  float bestXYsigMatch_ok_pfmva1;
  float bestXYsigMatch_ok_pfmva2;
  float bestXYsigMatch_ok_lptmva1;
  float bestXYsigMatch_ok_lptmva2;
  float bestXYsigMatch_ok_costhetaSK;
  float bestXYsigMatch_ok_costhetaSKCS;
  float bestXYsigMatch_ok_costhetaL;
  float bestXYsigMatch_ok_costhetaSK_gen;

  int bestCos2DMatch;
  int bestCos2DMatch_causeEle1;
  int bestCos2DMatch_causeEle2;
  int bestCos2DMatch_causeK;
  float bestCos2DMatch_notok_ele1pt;
  float bestCos2DMatch_notok_ele2pt;
  float bestCos2DMatch_notok_kpt;
  float bestCos2DMatch_notok_ele1eta;
  float bestCos2DMatch_notok_ele2eta;
  float bestCos2DMatch_notok_keta;
  float bestCos2DMatch_ok_ele1pt;
  float bestCos2DMatch_ok_ele2pt;
  float bestCos2DMatch_ok_kpt;
  float bestCos2DMatch_ok_ele1eta;
  float bestCos2DMatch_ok_ele2eta;
  float bestCos2DMatch_ok_keta;

  int bestPtSumMatch;
  int bestPtSumMatch_causeEle1;
  int bestPtSumMatch_causeEle2;
  int bestPtSumMatch_causeK;
  int bestKPtMatch;
  int bestKPtMatch_causeEle1;
  int bestKPtMatch_causeEle2;
  int bestKPtMatch_causeK;
};

#endif
