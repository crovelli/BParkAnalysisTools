#define correlationFromMcWithTnP_cxx
#include <TROOT.h>
#include <TH2.h>
#include <TFile.h>
#include <TTree.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TColor.h>
#include <iostream>

// To be run on MC dedicated formatted ntuples to get correlations between variables
// Separately for signal and background

using namespace std;

/** Get a value from a TTree by event id and name */
float get_ntuple_entry(TTree* nt, int i, std::string field) {
  float v;
  nt->SetBranchAddress(field.c_str(), &v);
  nt->GetEvent(i);
  nt->ResetBranchAddresses();
  return v;
}

// Prepare and plot the correlation matrix
void correlationFromMcWithTnP(TFile myfile, bool testLPT) {

  TTree *fChain = (TTree*)myfile.Get("tnpAna/fitter_tree");
  if (fChain == 0) cout << "problem with tree" << endl;
  Long64_t nentries = fChain->GetEntriesFast();   
  
  // get list of variable names
  std::vector<std::string> names;
  if (testLPT==1) names.push_back("probeMvaId");
  if (testLPT==0) names.push_back("probePfmvaId");
  names.push_back("probePt");
  names.push_back("probeDxySig");
  names.push_back("probeDzTrg");
  names.push_back("probeIso04Rel");

  // Matrix
  std::vector<float> matrix(names.size() * names.size());
  TH2F *matrixH = new TH2F("matrixH","correlation matrix",names.size(),1,names.size(),names.size(),1,names.size());
  matrixH->Sumw2();
  matrixH->SetTitle("");
  
  for (size_t i=1; i<=names.size(); i++) matrixH->GetXaxis()->SetBinLabel (i, (names[i-1]).c_str());
  for (size_t j=1; j<=names.size(); j++) matrixH->GetYaxis()->SetBinLabel (j, (names[j-1]).c_str());
  
  // Convert the ntuple to a vector, calculating means as we go
  std::vector<float> table(names.size() * nentries);
  std::vector<float> means(names.size(), 0);
  for (size_t i=0; i<nentries; i++) {
    for (size_t j=0; j<names.size(); j++) {
      if (i==0) cout << "conversion for variable = " <<  names[j] << endl;
      float v = get_ntuple_entry(fChain, i, names[j]);
      table[j + i * names.size()] = v;
      means[j] += v;
    }
  }

  // sums to means
  for (size_t i=0; i<names.size(); i++) means[i] /= nentries;

  // compute correlations
  cout << "correlation matrix computation" << endl;
  for (size_t i=0; i<names.size(); i++) {
    for (size_t j=i; j<names.size(); j++) {
      float t = 0;
      float dx2 = 0;
      float dy2 = 0;
      for (int k=0; k<nentries; k++) {
	float x1 = table[k * names.size() + i] - means[i];
	float x2 = table[k * names.size() + j] - means[j];
	t += x1 * x2;
	dx2 += x1 * x1;
	dy2 += x2 * x2;
      }
      matrix[i * names.size() + j] = t / sqrt(dx2 * dy2);
      matrixH->SetBinContent(i+1, j+1, matrix[i * names.size() + j] );
      matrixH->SetBinContent(j+1, i+1, matrix[i * names.size() + j] );
    }
  }

  
  // Save output
  TFile myOutFile("myOutFile.root","RECREATE");
  myOutFile.cd();
  matrixH->Write();
  myOutFile.Close();
  
  // Plot
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);

  TCanvas c("c","c",1);
  matrixH->Draw("colz");
  c.SaveAs("matrix.png");
 
}
