#ifndef nvtxInput_h
#define nvtxInput_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class nvtxInput {
 public :
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain
  
  // Fixed size dimensions of array or collections stored in the TTree if any.
  
  // Declaration of leaf types
  Int_t           nvtx;
  Float_t         pu_weight;
  Int_t           hlt9;
  Int_t           hlt12;
  
  // List of branches
  TBranch        *b_nvtx;   //!
  TBranch        *b_pu_weight;   //!
  TBranch        *b_hlt9;   //!
  TBranch        *b_hlt12;   //!
  
  nvtxInput(TTree *tree=0);
  virtual ~nvtxInput();
  virtual Int_t    Cut(Long64_t entry);
  virtual Int_t    GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void     Init(TTree *tree);
  virtual void     Loop();
  virtual Bool_t   Notify();
  virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef nvtxInput_cxx
nvtxInput::nvtxInput(TTree *tree) : fChain(0) 
{
  if (tree == 0) {
    //TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/eos/cms/store/user/crovelli/LowPtEle/Vertices/Vertices__BuToKJpsi_Toee_BParkNANO_mc_2020May16_ext__withVtxWeight_bin1.root");
    //TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/eos/cms/store/user/crovelli/LowPtEle/Vertices/Vertices__BuToKJpsi_Toee_BParkNANO_mc_2020May16_ext.root");
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/eos/cms/store/user/crovelli/LowPtEle/Vertices/Vertices__ParkingBPH1_Run2018D_part1_0000.root");
    if (!f || !f->IsOpen()) {
      //f = new TFile("/eos/cms/store/user/crovelli/LowPtEle/Vertices/Vertices__BuToKJpsi_Toee_BParkNANO_mc_2020May16_ext__withVtxWeight_bin1.root");
      //f = new TFile("/eos/cms/store/user/crovelli/LowPtEle/Vertices/Vertices__BuToKJpsi_Toee_BParkNANO_mc_2020May16_ext.root");
      f = new TFile("/eos/cms/store/user/crovelli/LowPtEle/Vertices/Vertices__ParkingBPH1_Run2018D_part1_0000.root");
    }
    f->GetObject("TaPtree",tree);
    
  }
  Init(tree);
}

nvtxInput::~nvtxInput()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

Int_t nvtxInput::GetEntry(Long64_t entry)
{
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}

Long64_t nvtxInput::LoadTree(Long64_t entry)
{
  if (!fChain) return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0) return centry;
  if (fChain->GetTreeNumber() != fCurrent) {
    fCurrent = fChain->GetTreeNumber();
    Notify();
  }
  return centry;
}

void nvtxInput::Init(TTree *tree)
{
  if (!tree) return;
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);
  
  fChain->SetBranchAddress("nvtx", &nvtx, &b_nvtx);
  fChain->SetBranchAddress("pu_weight", &pu_weight, &b_pu_weight);
  fChain->SetBranchAddress("hlt9", &hlt9, &b_hlt9);
  fChain->SetBranchAddress("hlt12", &hlt12, &b_hlt12);
  Notify();
}

Bool_t nvtxInput::Notify()
{
  return kTRUE;
}

void nvtxInput::Show(Long64_t entry)
{
  if (!fChain) return;
  fChain->Show(entry);
}

Int_t nvtxInput::Cut(Long64_t entry)
{
  return 1;
}

#endif // #ifdef nvtxInput_cxx
