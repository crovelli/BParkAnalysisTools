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
  bool isTag( int myEle );     
  bool isProbe( int myEle );     
  void bookOutputTree();
  void bookOutputHistos();

  // to compute weights for pileup
  std::vector<Double_t> puweights_;

  // settings
  bool dopureweight_;
  int sampleID;
  string puWFileName_;
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
  int    theLumi;
  int    nvtx;
  int    theSampleID;
  float  rho;
  float  theLumiWeight;
  float  pu_weight;
  float  pu_n;
  float perEveW;
  // 
  int selectedBSize;
  //
  vector <float> tag_pt={};   
  vector <float> tag_eta={};
  vector <float> tag_phi={};
  vector <bool> tag_isPF={};
  vector <bool> tag_isLowPt={};
  vector <float> tag_mvaId={};
  vector <float> tag_pfmvaId={};
  vector <float> tag_unBiased ={};
  vector <bool>  tag_matchMC={};
  //
  vector <float> probe_Bmass={};   
  vector <float> probe_Bpt={};   
  vector <float> probe_Bcos2D={};   
  vector <bool> probe_BmatchMC={};   
  vector <float> probe_pt={}; 
  vector <float> probe_eta ={};
  vector <float> probe_phi={};
  vector <bool> probe_isPF ={};
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
  vector <bool> probe_isTag ={};
  vector <float> probe_invMass={};
  vector <bool>  probe_matchMC={};
  //
  int selectedPairsSize;
};

#endif
