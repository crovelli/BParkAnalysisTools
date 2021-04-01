#ifndef Vertices_h
#define Vertices_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include <vector>
#include <string>

#include "./BParkBase.h"

using namespace std;

class Vertices : public BParkBase{
public:

  //! constructor
  Vertices(TTree *tree=0);
  //! destructor
  virtual ~Vertices();
  //! loop over events
  void Loop();
  void PrepareOutputs(std::string filename);             
  
private:
  
  // Analysis methods
  void bookOutputTree();
  void SetNvtxWeights(std::string nvtxWeightFile);
  float GetNvtxWeight(float nvtx);

  // to compute weights for pileup
  std::vector<Double_t> nvtxweights_;
  std::vector<Double_t> nvtxlowedge_;

  // settings
  bool donvtxreweight_;
  string nvtxWFileName_;

  // ---- outputs
  TFile* outFile_;
  TTree* outTree_;

  // dataset name
  std::string _datasetName;      

  //---output tree branches variables
  int nvtx;
  float pu_weight;
  int hlt9;
  int hlt12;
};

#endif
