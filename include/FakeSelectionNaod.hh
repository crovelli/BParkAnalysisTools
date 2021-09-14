#ifndef FakeSelectionNaod_h
#define FakeSelectionNaod_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include <vector>
#include <string>
#include <TLorentzVector.h>

#include "./BParkBase.h"

using namespace std;

class FakeSelectionNaod : public BParkBase{
public:

  //! constructor
  FakeSelectionNaod(TTree *tree=0);
  //! destructor
  virtual ~FakeSelectionNaod();
  //! loop over events
  void Loop();
  void PrepareOutputs(std::string filename);             
  
private:
  
  // Analysis methods
  void bookOutputTree();
  void bookOutputHistos();
  bool isMcB(int theB);
  bool isMcMuFromJPsi(int mu_idx);
  void SetNvtxWeights(std::string nvtxWeightFile);
  float GetNvtxWeight(float nvtx);
  
  // settings
  bool donvtxreweight_;
  int sampleID;
  string nvtxWFileName_;

  // to compute weights for pileup   
  std::vector<Double_t> nvtxweights_;
  std::vector<Double_t> nvtxlowedge_;

  // ---- outputs
  TFile* outFile_;
  TTree* outTree_;
  TH1F* h_entries;
  TH1F* h_selection;

  // dataset name
  std::string _datasetName;      

  //---output tree branches variables
  int    theRun;
  int    theEvent;
  int    nvtx;
  int    theSampleID;
  float  rho;
  float  pu_weight;
  float perEveW;
  //
  int hlt9;
  // 
  int selectedElesSize;
  //
  vector <float> Bmass={};
  vector <float> Bpt={};
  vector <float> Bcos2D={};
  vector <float> Bsvprob={};
  vector <float> Bxysig={};
  vector <int> BmatchMC={};
  vector <float> deltaR_ele_mu1={};   
  vector <float> deltaR_ele_mu2={};   
  vector <float> deltaR_ele_k={};   
  vector <float> k_pt={};
  vector <float> k_eta={};
  vector <float> k_phi={};
  vector <float> mu1_pt={};
  vector <float> mu1_eta={};
  vector <float> mu1_phi={};
  vector <int> mu1_isTriggering={};
  vector <int> mu1_matchMcFromJPsi={};
  vector <float> mu2_pt={};
  vector <float> mu2_eta={};
  vector <float> mu2_phi={};
  vector <int> mu2_isTriggering={};
  vector <int> mu2_matchMcFromJPsi={};
  vector <float> ele_pt={};
  vector <float> ele_eta={};
  vector <float> ele_phi={};
  vector <bool> ele_isPF={};
  vector <bool> ele_isPFOverlap={};
  vector <bool> ele_isLowPt={};
  vector <float> ele_mvaId={};
  vector <float> ele_pfmvaId={};
  vector <float> ele_dxySig={};
  vector <float> ele_dzSig={};
  vector <float> ele_pfRelIso={};
  vector <float> ele_trkRelIso={};
  vector <float> ele_fBrem={};
  vector <float> ele_unBiased={};
  vector <float> ele_ptBiased={};
  vector <bool> ele_convveto={};
  vector <float> probe_closeToEle_pt={};
  vector <float> probe_closeToEle_eta={};
  vector <float> probe_closeToEle_phi={};
};

#endif
