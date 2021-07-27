#define nvtxInput_cxx
#include "nvtxInput.h"
#include <TH1.h>
#include <iostream>

// For pileup reweighting using number of reconstructed vertices
// To be run on dedicated trees produced with Vertices.cc
// Create distribution with number of vertices, input to nvtxWeights

using namespace std; 

void nvtxInput::Loop()
{
  if (fChain == 0) return;
  
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  
  //int nbins = 36;
  //float edges[37] = { 0, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 50 };
  int nbins = 27;
  float edges[28] = { 0, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 50 };

  TH1F *H_nvtx = new TH1F("H_nvtx", "H_nvtx", nbins, edges);
  
  cout << nentries << endl;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    if (jentry%500000==0) cout << jentry << endl;

    if (hlt9==0) continue;

    /////if (nvtx>37) continue; // never in MC

    for (int i = 0; i <nbins; i++) {   
      if (edges[i]<=nvtx && edges[i+1]>nvtx) {  
	int ip1 = i+1;
	float thecontent = (edges[i]+edges[i+1])/2.;
	H_nvtx->Fill(thecontent, pu_weight); 
      }
    }
    if (nvtx>=edges[nbins]) { 
      float thecontent = edges[nbins-1]+0.01;
      H_nvtx->Fill(thecontent, pu_weight); 
    }
  }

  TFile fileNvtx("fileNumberVtx.root","RECREATE");
  H_nvtx->Write();
}
