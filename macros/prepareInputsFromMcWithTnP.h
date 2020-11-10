//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Oct 23 09:27:22 2020 by ROOT version 6.22/02
// from TTree fitter_tree/reduced tree for T&P
// found on file: Formatted_ParkingBPH1_Run2018D_ALL_probeLowPt.root
//////////////////////////////////////////////////////////

#ifndef prepareInputsFromMcWithTnP_h
#define prepareInputsFromMcWithTnP_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class prepareInputsFromMcWithTnP {
 public :
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain
  
  // Fixed size dimensions of array or collections stored in the TTree if any.
  
  // Declaration of leaf types
  Int_t           hlt_9;
  Int_t           hlt_12;
  Int_t           numvtx;
  Int_t           tagMatchMc;
  Int_t           tagMatchMcFromJPsi;
  Float_t         tagPt;
  Float_t         tagEta;
  Float_t         tagRelIso;
  Float_t         probePt;
  Float_t         probeAbsEta;
  Float_t         probeEta;
  Int_t           probeIsPF;
  Int_t           probeIsLowPt;
  Float_t         probeMvaId;
  Float_t         probePfmvaId;
  Float_t         probeUnBiased;
  Float_t         probePtBiased;
  Int_t           probeConvVeto;
  Int_t           probeMatchMc;
  Int_t           probeMatchMcFromJPsi;
  Float_t         probeTrkRelIso;
  Float_t         probePFRelIso;
  Float_t         probeFBrem;
  Float_t         probeDxySig;
  Float_t         probeDzSig;
  Float_t         K_pt;
  Float_t         B_mass;
  Float_t         pair_mass;
  Float_t         weight;
  
  // List of branches
  TBranch        *b_hlt_9;   //!
  TBranch        *b_hlt_12;   //!
  TBranch        *b_numvtx;   //!
  TBranch        *b_tagMatchMc;   //!
  TBranch        *b_tagMatchMcFromJPsi;   //!
  TBranch        *b_tagPt;   //!
  TBranch        *b_tagEta;   //!
  TBranch        *b_tagRelIso;   //!
  TBranch        *b_probePt;   //!
  TBranch        *b_probeAbsEta;   //!
  TBranch        *b_probeEta;   //!
  TBranch        *b_probeIsPF;   //!
  TBranch        *b_probeIsLowPt;   //!
  TBranch        *b_probeMvaId;   //!
  TBranch        *b_probePfmvaId;   //!
  TBranch        *b_probeUnBiased;   //!
  TBranch        *b_probePtBiased;   //!
  TBranch        *b_probeConvVeto;   //!
  TBranch        *b_probeMatchMc;   //!
  TBranch        *b_probeMatchMcFromJPsi;   //!
  TBranch        *b_probeTrkRelIso;   //!
  TBranch        *b_probePFRelIso;   //!
  TBranch        *b_probeFBrem;   //!
  TBranch        *b_probeDxySig;   //!
  TBranch        *b_probeDzSig;   //!
  TBranch        *b_K_pt;   //!
  TBranch        *b_B_mass;   //!
  TBranch        *b_pair_mass;   //!
  TBranch        *b_weight;   //!
  
  prepareInputsFromMcWithTnP(TTree *tree=0);
  virtual ~prepareInputsFromMcWithTnP();
  virtual Int_t    Cut(Long64_t entry);
  virtual Int_t    GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void     Init(TTree *tree);
  virtual void     Loop(bool applyWeight, bool testLPT);
  virtual Bool_t   Notify();
  virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef prepareInputsFromMcWithTnP_cxx
prepareInputsFromMcWithTnP::prepareInputsFromMcWithTnP(TTree *tree) : fChain(0) 
{
  // if parameter tree is not specified (or zero), connect the file
  // used to generate this class and read the Tree.
  if (tree == 0) {
    // TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Formatted_BuToKJpsi_Toee_BParkNANO_mc_2020May16_ext_probeLowPt.root");
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Formatted_BuToKJpsi_Toee_BParkNANO_mc_2020May16_ext_probePF.root");
    if (!f || !f->IsOpen()) {
      // f = new TFile("Formatted_BuToKJpsi_Toee_BParkNANO_mc_2020May16_ext_probeLowPt.root");
      f = new TFile("Formatted_BuToKJpsi_Toee_BParkNANO_mc_2020May16_ext_probePF.root");
    }
    // TDirectory * dir = (TDirectory*)f->Get("Formatted_BuToKJpsi_Toee_BParkNANO_mc_2020May16_ext_probeLowPt.root:/tnpAna");
    TDirectory * dir = (TDirectory*)f->Get("Formatted_BuToKJpsi_Toee_BParkNANO_mc_2020May16_ext_probePF.root:/tnpAna");
    dir->GetObject("fitter_tree",tree);
    
  }
  Init(tree);
}

prepareInputsFromMcWithTnP::~prepareInputsFromMcWithTnP()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

Int_t prepareInputsFromMcWithTnP::GetEntry(Long64_t entry)
{
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}

Long64_t prepareInputsFromMcWithTnP::LoadTree(Long64_t entry)
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

void prepareInputsFromMcWithTnP::Init(TTree *tree)
{
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses and branch
  // pointers of the tree will be set.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).
  
  // Set branch addresses and branch pointers
  if (!tree) return;
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);
  
  fChain->SetBranchAddress("hlt_9", &hlt_9, &b_hlt_9);
  fChain->SetBranchAddress("hlt_12", &hlt_12, &b_hlt_12);
  fChain->SetBranchAddress("numvtx", &numvtx, &b_numvtx);
  fChain->SetBranchAddress("tagMatchMc", &tagMatchMc, &b_tagMatchMc);
  fChain->SetBranchAddress("tagMatchMcFromJPsi", &tagMatchMcFromJPsi, &b_tagMatchMcFromJPsi);
  fChain->SetBranchAddress("tagPt", &tagPt, &b_tagPt);
  fChain->SetBranchAddress("tagEta", &tagEta, &b_tagEta);
  fChain->SetBranchAddress("tagRelIso", &tagRelIso, &b_tagRelIso);
  fChain->SetBranchAddress("probePt", &probePt, &b_probePt);
  fChain->SetBranchAddress("probeAbsEta", &probeAbsEta, &b_probeAbsEta);
  fChain->SetBranchAddress("probeEta", &probeEta, &b_probeEta);
  fChain->SetBranchAddress("probeIsPF", &probeIsPF, &b_probeIsPF);
  fChain->SetBranchAddress("probeIsLowPt", &probeIsLowPt, &b_probeIsLowPt);
  fChain->SetBranchAddress("probeMvaId", &probeMvaId, &b_probeMvaId);
  fChain->SetBranchAddress("probePfmvaId", &probePfmvaId, &b_probePfmvaId);
  fChain->SetBranchAddress("probeUnBiased", &probeUnBiased, &b_probeUnBiased);
  fChain->SetBranchAddress("probePtBiased", &probePtBiased, &b_probePtBiased);
  fChain->SetBranchAddress("probeConvVeto", &probeConvVeto, &b_probeConvVeto);
  fChain->SetBranchAddress("probeMatchMc", &probeMatchMc, &b_probeMatchMc);
  fChain->SetBranchAddress("probeMatchMcFromJPsi", &probeMatchMcFromJPsi, &b_probeMatchMcFromJPsi);
  fChain->SetBranchAddress("probeTrkRelIso", &probeTrkRelIso, &b_probeTrkRelIso);
  fChain->SetBranchAddress("probePFRelIso", &probePFRelIso, &b_probePFRelIso);
  fChain->SetBranchAddress("probeFBrem", &probeFBrem, &b_probeFBrem);
  fChain->SetBranchAddress("probeDxySig", &probeDxySig, &b_probeDxySig);
  fChain->SetBranchAddress("probeDzSig", &probeDzSig, &b_probeDzSig);
  fChain->SetBranchAddress("K_pt", &K_pt, &b_K_pt);
  fChain->SetBranchAddress("B_mass", &B_mass, &b_B_mass);
  fChain->SetBranchAddress("pair_mass", &pair_mass, &b_pair_mass);
  fChain->SetBranchAddress("weight", &weight, &b_weight);
  Notify();
}

Bool_t prepareInputsFromMcWithTnP::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.
  
  return kTRUE;
}

void prepareInputsFromMcWithTnP::Show(Long64_t entry)
{
  // Print contents of entry.
  // If entry is not specified, print current entry
  if (!fChain) return;
  fChain->Show(entry);
}

Int_t prepareInputsFromMcWithTnP::Cut(Long64_t entry)
{
  // This function may be called from Loop.
  // returns  1 if entry is accepted.
  // returns -1 otherwise.
  return 1;
}
#endif // #ifdef prepareInputsFromMcWithTnP_cxx
