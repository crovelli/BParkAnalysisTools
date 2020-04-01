//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Mar 30 18:28:53 2020 by ROOT version 6.18/04
// from TTree TaPtree/TaPtree
// found on file: ../data/Jan16prod/BuToKee_Toee_BParkNANO_mc_2020Jan16__okForAna.root
//////////////////////////////////////////////////////////

#ifndef plotsTripletsExtra_h
#define plotsTripletsExtra_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"

class plotsTripletsExtra {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           iHLT_Mu12_IP6;
   Int_t           iHLT_Mu9_IP6;
   Int_t           theEvent;
   Float_t         rho;
   Int_t           goodBSize;
   Int_t           goodTrueBSize;
   Int_t           goodCombBSize;
   Float_t         bestMatch_Bmass;
   Float_t         bestMatch_SvProb;
   Float_t         bestMatch_XYSig;
   Float_t         bestMatch_Cos2D;
   Float_t         bestMatch_KPt;
   Float_t         bestMatch_KEta;
   Float_t         bestMatch_Ele1Pt;
   Float_t         bestMatch_Ele2Pt;
   Float_t         bestMatch_MinPt;
   Float_t         bestMatch_Ele1Eta;
   Float_t         bestMatch_Ele2Eta;
   Float_t         bestMatch_Ele1pfmva;
   Float_t         bestMatch_Ele2pfmva;
   Float_t         bestMatch_Ele1lptmva;
   Float_t         bestMatch_Ele2lptmva;
   vector<float>   *goodCombB_svProb;
   vector<float>   *goodCombB_xySig;
   vector<float>   *goodCombB_cos2D;
   vector<float>   *goodCombB_kpt;
   vector<float>   *goodCombB_keta;
   vector<float>   *goodCombB_ele1pt;
   vector<float>   *goodCombB_ele2pt;
   vector<float>   *goodCombB_minpt;
   vector<float>   *goodCombB_ele1eta;
   vector<float>   *goodCombB_ele2eta;
   vector<float>   *goodCombB_ele1pfmva;
   vector<float>   *goodCombB_ele2pfmva;
   vector<float>   *goodCombB_ele1lptmva;
   vector<float>   *goodCombB_ele2lptmva;
   vector<int>     *goodCombB_notmatching;
   Int_t           bestSvProbMatch;
   Int_t           bestSvProbMatch_notmatching;
   Int_t           bestXYsigMatch;
   Int_t           bestXYsigMatch_notmatching;
   Int_t           numberBetterSvProbTriplets;
   Int_t           numberBetterXYsigTriplets;

   // List of branches
   TBranch        *b_iHLT_Mu12_IP6;   //!
   TBranch        *b_iHLT_Mu9_IP6;   //!
   TBranch        *b_theEvent;   //!
   TBranch        *b_rho;   //!
   TBranch        *b_goodBSize;   //!
   TBranch        *b_goodTrueBSize;   //!
   TBranch        *b_goodCombBSize;   //!
   TBranch        *b_bestMatch_Bmass;   //!
   TBranch        *b_bestMatch_SvProb;   //!
   TBranch        *b_bestMatch_XYSig;   //!
   TBranch        *b_bestMatch_Cos2D;   //!
   TBranch        *b_bestMatch_KPt;   //!
   TBranch        *b_bestMatch_KEta;   //!
   TBranch        *b_bestMatch_Ele1Pt;   //!
   TBranch        *b_bestMatch_Ele2Pt;   //!
   TBranch        *b_bestMatch_MinPt;   //!
   TBranch        *b_bestMatch_Ele1Eta;   //!
   TBranch        *b_bestMatch_Ele2Eta;   //!
   TBranch        *b_bestMatch_Ele1pfmva;   //!
   TBranch        *b_bestMatch_Ele2pfmva;   //!
   TBranch        *b_bestMatch_Ele1lptmva;   //!
   TBranch        *b_bestMatch_Ele2lptmva;   //!
   TBranch        *b_goodCombB_svProb;   //!
   TBranch        *b_goodCombB_xySig;   //!
   TBranch        *b_goodCombB_cos2D;   //!
   TBranch        *b_goodCombB_kpt;   //!
   TBranch        *b_goodCombB_keta;   //!
   TBranch        *b_goodCombB_ele1pt;   //!
   TBranch        *b_goodCombB_ele2pt;   //!
   TBranch        *b_goodCombB_minpt;   //!
   TBranch        *b_goodCombB_ele1eta;   //!
   TBranch        *b_goodCombB_ele2eta;   //!
   TBranch        *b_goodCombB_ele1pfmva;   //!
   TBranch        *b_goodCombB_ele2pfmva;   //!
   TBranch        *b_goodCombB_ele1lptmva;   //!
   TBranch        *b_goodCombB_ele2lptmva;   //!
   TBranch        *b_goodCombB_notmatching;   //!
   TBranch        *b_bestSvProbMatch;   //!
   TBranch        *b_bestSvProbMatch_notmatching;   //!
   TBranch        *b_bestXYsigMatch;   //!
   TBranch        *b_bestXYsigMatch_notmatching;   //!
   TBranch        *b_numberBetterSvProbTriplets;   //!
   TBranch        *b_numberBetterXYsigTriplets;   //!

   plotsTripletsExtra(TTree *tree=0);
   virtual ~plotsTripletsExtra();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef plotsTripletsExtra_cxx
plotsTripletsExtra::plotsTripletsExtra(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../data/Jan16prod/BuToKee_Toee_BParkNANO_mc_2020Jan16.txt.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("../data/Jan16prod/BuToKee_Toee_BParkNANO_mc_2020Jan16.txt.root");
      }
      f->GetObject("TaPtree",tree);

   }
   Init(tree);
}

plotsTripletsExtra::~plotsTripletsExtra()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t plotsTripletsExtra::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t plotsTripletsExtra::LoadTree(Long64_t entry)
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

void plotsTripletsExtra::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   goodCombB_svProb = 0;
   goodCombB_xySig = 0;
   goodCombB_cos2D = 0;
   goodCombB_kpt = 0;
   goodCombB_keta = 0;
   goodCombB_ele1pt = 0;
   goodCombB_ele2pt = 0;
   goodCombB_minpt = 0;
   goodCombB_ele1eta = 0;
   goodCombB_ele2eta = 0;
   goodCombB_ele1pfmva = 0;
   goodCombB_ele2pfmva = 0;
   goodCombB_ele1lptmva = 0;
   goodCombB_ele2lptmva = 0;
   goodCombB_notmatching = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("iHLT_Mu12_IP6", &iHLT_Mu12_IP6, &b_iHLT_Mu12_IP6);
   fChain->SetBranchAddress("iHLT_Mu9_IP6", &iHLT_Mu9_IP6, &b_iHLT_Mu9_IP6);
   fChain->SetBranchAddress("theEvent", &theEvent, &b_theEvent);
   fChain->SetBranchAddress("rho", &rho, &b_rho);
   fChain->SetBranchAddress("goodBSize", &goodBSize, &b_goodBSize);
   fChain->SetBranchAddress("goodTrueBSize", &goodTrueBSize, &b_goodTrueBSize);
   fChain->SetBranchAddress("goodCombBSize", &goodCombBSize, &b_goodCombBSize);
   fChain->SetBranchAddress("bestMatch_Bmass", &bestMatch_Bmass, &b_bestMatch_Bmass);
   fChain->SetBranchAddress("bestMatch_SvProb", &bestMatch_SvProb, &b_bestMatch_SvProb);
   fChain->SetBranchAddress("bestMatch_XYSig", &bestMatch_XYSig, &b_bestMatch_XYSig);
   fChain->SetBranchAddress("bestMatch_Cos2D", &bestMatch_Cos2D, &b_bestMatch_Cos2D);
   fChain->SetBranchAddress("bestMatch_KPt", &bestMatch_KPt, &b_bestMatch_KPt);
   fChain->SetBranchAddress("bestMatch_KEta", &bestMatch_KEta, &b_bestMatch_KEta);
   fChain->SetBranchAddress("bestMatch_Ele1Pt", &bestMatch_Ele1Pt, &b_bestMatch_Ele1Pt);
   fChain->SetBranchAddress("bestMatch_Ele2Pt", &bestMatch_Ele2Pt, &b_bestMatch_Ele2Pt);
   fChain->SetBranchAddress("bestMatch_MinPt", &bestMatch_MinPt, &b_bestMatch_MinPt);
   fChain->SetBranchAddress("bestMatch_Ele1Eta", &bestMatch_Ele1Eta, &b_bestMatch_Ele1Eta);
   fChain->SetBranchAddress("bestMatch_Ele2Eta", &bestMatch_Ele2Eta, &b_bestMatch_Ele2Eta);
   fChain->SetBranchAddress("bestMatch_Ele1pfmva", &bestMatch_Ele1pfmva, &b_bestMatch_Ele1pfmva);
   fChain->SetBranchAddress("bestMatch_Ele2pfmva", &bestMatch_Ele2pfmva, &b_bestMatch_Ele2pfmva);
   fChain->SetBranchAddress("bestMatch_Ele1lptmva", &bestMatch_Ele1lptmva, &b_bestMatch_Ele1lptmva);
   fChain->SetBranchAddress("bestMatch_Ele2lptmva", &bestMatch_Ele2lptmva, &b_bestMatch_Ele2lptmva);
   fChain->SetBranchAddress("goodCombB_svProb", &goodCombB_svProb, &b_goodCombB_svProb);
   fChain->SetBranchAddress("goodCombB_xySig", &goodCombB_xySig, &b_goodCombB_xySig);
   fChain->SetBranchAddress("goodCombB_cos2D", &goodCombB_cos2D, &b_goodCombB_cos2D);
   fChain->SetBranchAddress("goodCombB_kpt", &goodCombB_kpt, &b_goodCombB_kpt);
   fChain->SetBranchAddress("goodCombB_keta", &goodCombB_keta, &b_goodCombB_keta);
   fChain->SetBranchAddress("goodCombB_ele1pt", &goodCombB_ele1pt, &b_goodCombB_ele1pt);
   fChain->SetBranchAddress("goodCombB_ele2pt", &goodCombB_ele2pt, &b_goodCombB_ele2pt);
   fChain->SetBranchAddress("goodCombB_minpt", &goodCombB_minpt, &b_goodCombB_minpt);
   fChain->SetBranchAddress("goodCombB_ele1eta", &goodCombB_ele1eta, &b_goodCombB_ele1eta);
   fChain->SetBranchAddress("goodCombB_ele2eta", &goodCombB_ele2eta, &b_goodCombB_ele2eta);
   fChain->SetBranchAddress("goodCombB_ele1pfmva", &goodCombB_ele1pfmva, &b_goodCombB_ele1pfmva);
   fChain->SetBranchAddress("goodCombB_ele2pfmva", &goodCombB_ele2pfmva, &b_goodCombB_ele2pfmva);
   fChain->SetBranchAddress("goodCombB_ele1lptmva", &goodCombB_ele1lptmva, &b_goodCombB_ele1lptmva);
   fChain->SetBranchAddress("goodCombB_ele2lptmva", &goodCombB_ele2lptmva, &b_goodCombB_ele2lptmva);
   fChain->SetBranchAddress("goodCombB_notmatching", &goodCombB_notmatching, &b_goodCombB_notmatching);
   fChain->SetBranchAddress("bestSvProbMatch", &bestSvProbMatch, &b_bestSvProbMatch);
   fChain->SetBranchAddress("bestSvProbMatch_notmatching", &bestSvProbMatch_notmatching, &b_bestSvProbMatch_notmatching);
   fChain->SetBranchAddress("bestXYsigMatch", &bestXYsigMatch, &b_bestXYsigMatch);
   fChain->SetBranchAddress("bestXYsigMatch_notmatching", &bestXYsigMatch_notmatching, &b_bestXYsigMatch_notmatching);
   fChain->SetBranchAddress("numberBetterSvProbTriplets", &numberBetterSvProbTriplets, &b_numberBetterSvProbTriplets);
   fChain->SetBranchAddress("numberBetterXYsigTriplets", &numberBetterXYsigTriplets, &b_numberBetterXYsigTriplets);
   Notify();
}

Bool_t plotsTripletsExtra::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void plotsTripletsExtra::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t plotsTripletsExtra::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef plotsTripletsExtra_cxx
