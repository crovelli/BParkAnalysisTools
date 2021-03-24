#ifndef TaPJpsiSelectionNaod_h
#define TaPJpsiSelectionNaod_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include <vector>
#include <string>
#include <TLorentzVector.h>

#include "./BParkBase.h"

using namespace std;

class TaPJpsiSelectionNaod : public BParkBase{
public:

  //! constructor
  TaPJpsiSelectionNaod(TTree *tree=0);
  //! destructor
  virtual ~TaPJpsiSelectionNaod();
  //! loop over events
  void Loop();
  void PrepareOutputs(std::string filename);             
  
private:
  
  // Analysis methods
  bool isMcB( int myB );
  bool isMcEleFromJPsi (int myEle );
  void bookOutputTree();
  void bookOutputHistos();
  void SetNvtxWeights(std::string nvtxWeightFile);
  float GetNvtxWeight(float nvtx);

  // to compute weights for pileup
  std::vector<Double_t> nvtxweights_;
  std::vector<Double_t> nvtxlowedge_;

  // settings
  bool donvtxreweight_;
  int sampleID;
  string nvtxWFileName_;
  float lumiWeight_;

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
  float  pu_n;
  float perEveW;
  //
  int hlt9;
  int hlt12;
  // 
  int selectedBSize;
  //
  vector <float> tag_pt={};   
  vector <float> tag_eta={};
  vector <float> tag_phi={};
  vector <bool> tag_isPF={};
  vector <bool> tag_isPFOverlap={};
  vector <bool> tag_isLowPt={};
  vector <float> tag_mvaId={};
  vector <float> tag_pfmvaId={};
  vector <bool>  tag_matchMcFromJPsi={};
  vector <bool>  tag_matchMc={};
  vector <bool> tag_convveto={};
  vector <float> tag_pfRelIso={};
  //
  vector <float> probe_Bmass={};   
  vector <float> probe_Bpt={};   
  vector <float> probe_Bcos2D={};   
  vector <float> probe_Bsvprob={};
  vector <float> probe_Bxysig={};
  vector <bool> probe_BmatchMC={};   
  vector <float> probe_Kpt={}; 
  vector <float> probe_pt={}; 
  vector <float> probe_eta ={};
  vector <float> probe_phi={};
  vector <bool> probe_isPF ={};
  vector <bool> probe_isPFOverlap ={};
  vector <bool> probe_isLowPt ={};
  vector <float> probe_mvaId={};        
  vector <float> probe_pfmvaId={};        
  vector <float> probe_dxySig={};        
  vector <float> probe_dzSig={};        
  vector <float> probe_pfRelIso ={};
  vector <float> probe_trkRelIso ={};
  vector <float> probe_fBrem ={};
  vector <float> probe_unBiased ={};
  vector <float> probe_ptBiased ={};
  vector <bool> probe_convveto={};
  vector <float> probe_invMass={};
  vector <bool>  probe_matchMcFromJPsi={};
  vector <bool>  probe_matchMc={};
  //
  int selectedPairsSize;
};

#endif
