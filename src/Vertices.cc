#include "TMath.h"
#include "TTree.h"

#include <iostream> 

#include "../include/Vertices.hh"    

using namespace std;

Vertices::Vertices(TTree *tree)     
  : BParkBase(tree) {        

  // Chiara: to be set by hand   
  donvtxreweight_ = 0;    // 
  nvtxWFileName_ = "/afs/cern.ch/user/c/crovelli/public/bphys/nvtxWeights__bin1.root";
}

Vertices::~Vertices() {

  // output
  outFile_->cd();
  outTree_->Write();
  outFile_->Close();
}     

void Vertices::Loop() {

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

    // Trigger 
    int iHLT_Mu12_IP6 = (int)HLT_Mu12_IP6;
    int iHLT_Mu9_IP6  = (int)HLT_Mu9_IP6;
    hlt9  = iHLT_Mu9_IP6;
    hlt12 = iHLT_Mu12_IP6;
    
    // at least HLT_9 required
    if (hlt9==0) continue;

    // # Vertices
    nvtx = PV_npvs;

    // PU weight (only if requested)
    pu_weight = 1.;
    if (donvtxreweight_==1) {     
      pu_weight = GetNvtxWeight(nvtx);         
    }

    // Filling the output tree
    outTree_->Fill();
  }
}

void Vertices::SetNvtxWeights(std::string nvtxWeightFile) {

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
  
  for (int i = 0; i<nvtxweights->GetNbinsX(); i++) {
    float weight=weights->GetBinContent(i+1);
    float lowedge=weights->GetBinLowEdge(i+1);
    nvtxweights_.push_back(weight);        // weight for bins from 1 to N at position 0 to N-1 => assume no over/under-flow
    nvtxlowedge_.push_back(lowedge);       // low-edge of bins from from 1 to N at position 0 to N-1
  }
}

float Vertices::GetNvtxWeight(float nvtx) {

  int thesize   = nvtxlowedge_.size();    
  int thesizem1 = nvtxlowedge_.size()-1;  
  float weight=1;

  if (thesize>0 && donvtxreweight_) {
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

void Vertices::PrepareOutputs(std::string filename) 
{
  _datasetName=filename;

  std::string outname = _datasetName+".root";    
  cout << "output: " << outname << endl;
  outFile_ = new TFile(outname.c_str(),"RECREATE");

  bookOutputTree();

  // loading weights for pileup if needed
  if (donvtxreweight_) SetNvtxWeights(nvtxWFileName_);
};


void Vertices::bookOutputTree() 
{
  outTree_ = new TTree("TaPtree", "TaPtree");
  
  cout << "Booking tree" << endl;
  
  outTree_->Branch("nvtx", &nvtx, "nvtx/I");    
  outTree_->Branch("pu_weight", &pu_weight, "pu_weight/F");    

  outTree_->Branch("hlt9", &hlt9, "hlt9/I");
  outTree_->Branch("hlt12", &hlt12, "hlt12/I");
}

