#define plotsTripletsExtra_cxx
#include "plotsTripletsExtra.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <iostream> 

using namespace std;

void plotsTripletsExtra::Loop()
{
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  
  // Histos
  TH1F *h_numB      = new TH1F("h_numB", "     h_numB",      25, -0.5, 24.5);
  TH1F *h_numCombB  = new TH1F("h_numCombB",  "h_numCombB",  25, -0.5, 24.5);
  TH1F *h_numTrueB  = new TH1F("h_numTrueB",  "h_numTrueB",  10, -0.5, 9.5);

  // True Bs
  TH1F *h_trueBestMatch_svProb   = new TH1F("h_trueBestMatch_svProb", "h_trueBestMatch_svProb", 100, -0.01, 1.01);
  TH1F *h_trueBestMatch_xysig    = new TH1F("h_trueBestMatch_xysig",  "h_trueBestMatch_xysig",  50,  0., 100.);
  TH1F *h_trueBestMatch_cos2d    = new TH1F("h_trueBestMatch_cos2d",  "h_trueBestMatch_cos2d",  100,  0.989, 1.001);
  TH1F *h_trueBestMatch_kpt      = new TH1F("h_trueBestMatch_kpt",      "h_trueBestMatch_kpt",      60, 0., 30.);
  TH1F *h_trueBestMatch_keta     = new TH1F("h_trueBestMatch_keta",    "h_trueBestMatch_keta",    50,-2.5,2.5);
  TH1F *h_trueBestMatch_ele1eta  = new TH1F("h_trueBestMatch_ele1eta", "h_trueBestMatch_ele1eta", 50,-2.5,2.5);
  TH1F *h_trueBestMatch_ele2eta  = new TH1F("h_trueBestMatch_ele2eta", "h_trueBestMatch_ele2eta", 50,-2.5,2.5);
  TH1F *h_trueBestMatch_ele1pt   = new TH1F("h_trueBestMatch_ele1pt",  "h_trueBestMatch_ele1pt",  60, 0., 30.);
  TH1F *h_trueBestMatch_ele2pt   = new TH1F("h_trueBestMatch_ele2pt",  "h_trueBestMatch_ele2pt",  60, 0., 30.);
  TH1F *h_trueBestMatch_pfmva1   = new TH1F("h_trueBestMatch_pfmva1",  "h_trueBestMatch_pfmva1",  50,-10.,10.);
  TH1F *h_trueBestMatch_pfmva2   = new TH1F("h_trueBestMatch_pfmva2",  "h_trueBestMatch_pfmva2",  50,-10.,10.);
  TH1F *h_trueBestMatch_lptmva1  = new TH1F("h_trueBestMatch_lptmva1", "h_trueBestMatch_lptmva1", 50,-10.,10.);
  TH1F *h_trueBestMatch_lptmva2  = new TH1F("h_trueBestMatch_lptmva2", "h_trueBestMatch_lptmva2", 50,-10.,10.);
  TH1F *hz_trueBestMatch_kpt     = new TH1F("hz_trueBestMatch_kpt",    "hz_trueBestMatch_kpt",    50, 0., 5.);
  TH1F *hz_trueBestMatch_ele1pt  = new TH1F("hz_trueBestMatch_ele1pt", "hz_trueBestMatch_ele1pt", 50, 0., 10.);
  TH1F *hz_trueBestMatch_ele2pt  = new TH1F("hz_trueBestMatch_ele2pt", "hz_trueBestMatch_ele2pt", 50, 0., 5.);
  TH1F *hz_trueBestMatch_minpt   = new TH1F("hz_trueBestMatch_minpt",  "hz_trueBestMatch_minpt",  50, 0., 5.);
 
  // Combinatorics - all
  TH1F *h_comb_svProb   = new TH1F("h_comb_svProb",  "h_comb_svProb", 100, -0.01, 1.01);
  TH1F *h_comb_xysig    = new TH1F("h_comb_xysig",   "h_comb_xysig",   50, 0., 100.);
  TH1F *h_comb_cos2d    = new TH1F("h_comb_cos2d",   "h_comb_cos2d",  100, 0.989, 1.001);
  TH1F *h_comb_kpt      = new TH1F("h_comb_kpt",     "h_comb_kpt",     60, 0., 30.);
  TH1F *h_comb_keta     = new TH1F("h_comb_keta",    "h_comb_keta",    50,-2.5,2.5);
  TH1F *h_comb_ele1eta  = new TH1F("h_comb_ele1eta", "h_comb_ele1eta", 50,-2.5,2.5);
  TH1F *h_comb_ele2eta  = new TH1F("h_comb_ele2eta", "h_comb_ele2eta", 50,-2.5,2.5);
  TH1F *h_comb_ele1pt   = new TH1F("h_comb_ele1pt",  "h_comb_ele1pt",  60, 0., 30.);
  TH1F *h_comb_ele2pt   = new TH1F("h_comb_ele2pt",  "h_comb_ele2pt",  60, 0., 30.);
  TH1F *h_comb_pfmva1   = new TH1F("h_comb_pfmva1",  "h_comb_pfmva1",  50,-10.,10.);
  TH1F *h_comb_pfmva2   = new TH1F("h_comb_pfmva2",  "h_comb_pfmva2",  50,-10.,10.);
  TH1F *h_comb_lptmva1  = new TH1F("h_comb_lptmva1", "h_comb_lptmva1", 50,-10.,10.);
  TH1F *h_comb_lptmva2  = new TH1F("h_comb_lptmva2", "h_comb_lptmva2", 50,-10.,10.);
  TH1F *h_comb_pfmva1_badele1  = new TH1F("h_comb_pfmva1_badele1",  "h_comb_pfmva1_badele1",  50,-10.,10.);
  TH1F *h_comb_pfmva2_badele2  = new TH1F("h_comb_pfmva2_badele2",  "h_comb_pfmva2_badele2",  50,-10.,10.);
  TH1F *h_comb_lptmva1_badele1 = new TH1F("h_comb_lptmva1_badele1", "h_comb_lptmva1_badele1", 50,-10.,10.);
  TH1F *h_comb_lptmva2_badele2 = new TH1F("h_comb_lptmva2_badele2", "h_comb_lptmva2_badele2", 50,-10.,10.);
  TH1F *hz_comb_kpt    = new TH1F("hz_comb_kpt",    "hz_comb_kpt",    50, 0., 5.);
  TH1F *hz_comb_ele1pt = new TH1F("hz_comb_ele1pt", "hz_comb_ele1pt", 50, 0., 10.);
  TH1F *hz_comb_ele2pt = new TH1F("hz_comb_ele2pt", "hz_comb_ele2pt", 50, 0., 5.);
  TH1F *hz_comb_kpt_badk       = new TH1F("hz_comb_kpt_badk",       "hz_comb_kpt_badk",        50, 0., 5.);
  TH1F *hz_comb_ele1pt_badele1 = new TH1F("hz_comb_ele1pt_badele1", "hz_comb_ele1pt_badele1",  50, 0., 10.);
  TH1F *hz_comb_ele2pt_badele2 = new TH1F("hz_comb_ele2pt_badele2", "hz_comb_ele2pt_badele2",  50, 0., 5.);
  TH1F *hz_comb_minpt  = new TH1F("hz_comb_minpt",  "hz_comb_minpt",  50, 0., 5.);
  TH1F *h_comb_notmatching = new TH1F("h_comb_notmatching", "h_comb_notmatching", 9, -1.5, 7.5);

  // Combinatorics - 1 not matched object
  TH1F *h_comb1_svProb   = new TH1F("h_comb1_svProb",  "h_comb1_svProb", 100, -0.01, 1.01);
  TH1F *h_comb1_xysig    = new TH1F("h_comb1_xysig",   "h_comb1_xysig",   50,  0., 100.);
  TH1F *h_comb1_cos2d    = new TH1F("h_comb1_cos2d",   "h_comb1_cos2d",  100,  0.989, 1.001);
  TH1F *h_comb1_kpt      = new TH1F("h_comb1_kpt",     "h_comb1_kpt",     60, 0., 30.);
  TH1F *h_comb1_keta     = new TH1F("h_comb1_keta",    "h_comb1_keta",    50,-2.5,2.5);
  TH1F *h_comb1_ele1eta  = new TH1F("h_comb1_ele1eta", "h_comb1_ele1eta", 50,-2.5,2.5);
  TH1F *h_comb1_ele2eta  = new TH1F("h_comb1_ele2eta", "h_comb1_ele2eta", 50,-2.5,2.5);
  TH1F *h_comb1_ele1pt   = new TH1F("h_comb1_ele1pt",  "h_comb1_ele1pt",  60, 0., 30.);
  TH1F *h_comb1_ele2pt   = new TH1F("h_comb1_ele2pt",  "h_comb1_ele2pt",  60, 0., 30.);
  TH1F *h_comb1_pfmva1   = new TH1F("h_comb1_pfmva1",  "h_comb1_pfmva1",  50,-10.,10.);
  TH1F *h_comb1_pfmva2   = new TH1F("h_comb1_pfmva2",  "h_comb1_pfmva2",  50,-10.,10.);
  TH1F *h_comb1_lptmva1  = new TH1F("h_comb1_lptmva1", "h_comb1_lptmva1", 50,-10.,10.);
  TH1F *h_comb1_lptmva2  = new TH1F("h_comb1_lptmva2", "h_comb1_lptmva2", 50,-10.,10.);
  TH1F *h_comb1_pfmva1_badele1  = new TH1F("h_comb1_pfmva1_badele1",  "h_comb1_pfmva1_badele1",  50,-10.,10.);
  TH1F *h_comb1_pfmva2_badele2  = new TH1F("h_comb1_pfmva2_badele2",  "h_comb1_pfmva2_badele2",  50,-10.,10.);
  TH1F *h_comb1_lptmva1_badele1 = new TH1F("h_comb1_lptmva1_badele1", "h_comb1_lptmva1_badele1", 50,-10.,10.);
  TH1F *h_comb1_lptmva2_badele2 = new TH1F("h_comb1_lptmva2_badele2", "h_comb1_lptmva2_badele2", 50,-10.,10.);
  TH1F *hz_comb1_kpt    = new TH1F("hz_comb1_kpt",    "hz_comb1_kpt",    50, 0., 5.);
  TH1F *hz_comb1_ele1pt = new TH1F("hz_comb1_ele1pt", "hz_comb1_ele1pt", 50, 0., 10.);
  TH1F *hz_comb1_ele2pt = new TH1F("hz_comb1_ele2pt", "hz_comb1_ele2pt", 50, 0., 5.);
  TH1F *hz_comb1_kpt_badk       = new TH1F("hz_comb1_kpt_badk",       "hz_comb1_kpt_badk",        50, 0., 5.);
  TH1F *hz_comb1_ele1pt_badele1 = new TH1F("hz_comb1_ele1pt_badele1", "hz_comb1_ele1pt_badele1",  50, 0., 10.);
  TH1F *hz_comb1_ele2pt_badele2 = new TH1F("hz_comb1_ele2pt_badele2", "hz_comb1_ele2pt_badele2",  50, 0., 5.);
  TH1F *hz_comb1_minpt  = new TH1F("hz_comb1_minpt",  "hz_comb1_minpt",  50, 0., 5.);
  TH1F *h_comb1_notmatching = new TH1F("h_comb1_notmatching", "h_comb1_notmatching", 9, -1.5, 7.5);

  // Combinatorics - 2 not matched objects
  TH1F *h_comb2_svProb   = new TH1F("h_comb2_svProb",  "h_comb2_svProb", 100, -0.01, 1.01);
  TH1F *h_comb2_xysig    = new TH1F("h_comb2_xysig",   "h_comb2_xysig",   50, 0., 100.);
  TH1F *h_comb2_cos2d    = new TH1F("h_comb2_cos2d",   "h_comb2_cos2d",  100, 0.989, 1.001);
  TH1F *h_comb2_kpt      = new TH1F("h_comb2_kpt",     "h_comb2_kpt",     60, 0., 30.);
  TH1F *h_comb2_keta     = new TH1F("h_comb2_keta",    "h_comb2_keta",    50,-2.5,2.5);
  TH1F *h_comb2_ele1eta  = new TH1F("h_comb2_ele1eta", "h_comb2_ele1eta", 50,-2.5,2.5);
  TH1F *h_comb2_ele2eta  = new TH1F("h_comb2_ele2eta", "h_comb2_ele2eta", 50,-2.5,2.5);
  TH1F *h_comb2_ele1pt   = new TH1F("h_comb2_ele1pt",  "h_comb2_ele1pt",  60, 0., 30.);
  TH1F *h_comb2_ele2pt   = new TH1F("h_comb2_ele2pt",  "h_comb2_ele2pt",  60, 0., 30.);
  TH1F *h_comb2_pfmva1   = new TH1F("h_comb2_pfmva1",  "h_comb2_pfmva1",  50,-10.,10.);
  TH1F *h_comb2_pfmva2   = new TH1F("h_comb2_pfmva2",  "h_comb2_pfmva2",  50,-10.,10.);
  TH1F *h_comb2_lptmva1  = new TH1F("h_comb2_lptmva1", "h_comb2_lptmva1", 50,-10.,10.);
  TH1F *h_comb2_lptmva2  = new TH1F("h_comb2_lptmva2", "h_comb2_lptmva2", 50,-10.,10.);
  TH1F *h_comb2_pfmva1_badele1  = new TH1F("h_comb2_pfmva1_badele1",  "h_comb2_pfmva1_badele1",  50,-10.,10.);
  TH1F *h_comb2_pfmva2_badele2  = new TH1F("h_comb2_pfmva2_badele2",  "h_comb2_pfmva2_badele2",  50,-10.,10.);
  TH1F *h_comb2_lptmva1_badele1 = new TH1F("h_comb2_lptmva1_badele1", "h_comb2_lptmva1_badele1", 50,-10.,10.);
  TH1F *h_comb2_lptmva2_badele2 = new TH1F("h_comb2_lptmva2_badele2", "h_comb2_lptmva2_badele2", 50,-10.,10.);
  TH1F *hz_comb2_kpt    = new TH1F("hz_comb2_kpt",    "hz_comb2_kpt",    50, 0., 5.);
  TH1F *hz_comb2_ele1pt = new TH1F("hz_comb2_ele1pt", "hz_comb2_ele1pt", 50, 0., 10.);
  TH1F *hz_comb2_ele2pt = new TH1F("hz_comb2_ele2pt", "hz_comb2_ele2pt", 50, 0., 5.);
  TH1F *hz_comb2_kpt_badk       = new TH1F("hz_comb2_kpt_badk",       "hz_comb2_kpt_badk",        50, 0., 5.);
  TH1F *hz_comb2_ele1pt_badele1 = new TH1F("hz_comb2_ele1pt_badele1", "hz_comb2_ele1pt_badele1",  50, 0., 10.);
  TH1F *hz_comb2_ele2pt_badele2 = new TH1F("hz_comb2_ele2pt_badele2", "hz_comb2_ele2pt_badele2",  50, 0., 5.);
  TH1F *hz_comb2_minpt  = new TH1F("hz_comb2_minpt",  "hz_comb2_minpt",  50, 0., 5.);
  TH1F *h_comb2_notmatching = new TH1F("h_comb2_notmatching", "h_comb2_notmatching", 9, -1.5, 7.5);

  // Combinatorics - 3 not matched objects
  TH1F *h_comb3_svProb   = new TH1F("h_comb3_svProb",  "h_comb3_svProb", 100, -0.01, 1.01);
  TH1F *h_comb3_xysig    = new TH1F("h_comb3_xysig",   "h_comb3_xysig",   50, 0., 100.);
  TH1F *h_comb3_cos2d    = new TH1F("h_comb3_cos2d",   "h_comb3_cos2d",  100, 0.989, 1.001);
  TH1F *h_comb3_kpt      = new TH1F("h_comb3_kpt",     "h_comb3_kpt",     60, 0., 30.);
  TH1F *h_comb3_keta     = new TH1F("h_comb3_keta",    "h_comb3_keta",    50,-2.5,2.5);
  TH1F *h_comb3_ele1eta  = new TH1F("h_comb3_ele1eta", "h_comb3_ele1eta", 50,-2.5,2.5);
  TH1F *h_comb3_ele2eta  = new TH1F("h_comb3_ele2eta", "h_comb3_ele2eta", 50,-2.5,2.5);
  TH1F *h_comb3_ele1pt   = new TH1F("h_comb3_ele1pt",  "h_comb3_ele1pt",  60, 0., 30.);
  TH1F *h_comb3_ele2pt   = new TH1F("h_comb3_ele2pt",  "h_comb3_ele2pt",  60, 0., 30.);
  TH1F *h_comb3_pfmva1   = new TH1F("h_comb3_pfmva1",  "h_comb3_pfmva1",  50,-10.,10.);
  TH1F *h_comb3_pfmva2   = new TH1F("h_comb3_pfmva2",  "h_comb3_pfmva2",  50,-10.,10.);
  TH1F *h_comb3_lptmva1  = new TH1F("h_comb3_lptmva1", "h_comb3_lptmva1", 50,-10.,10.);
  TH1F *h_comb3_lptmva2  = new TH1F("h_comb3_lptmva2", "h_comb3_lptmva2", 50,-10.,10.);
  TH1F *h_comb3_pfmva1_badele1  = new TH1F("h_comb3_pfmva1_badele1",  "h_comb3_pfmva1_badele1",  50,-10.,10.);
  TH1F *h_comb3_pfmva2_badele2  = new TH1F("h_comb3_pfmva2_badele2",  "h_comb3_pfmva2_badele2",  50,-10.,10.);
  TH1F *h_comb3_lptmva1_badele1 = new TH1F("h_comb3_lptmva1_badele1", "h_comb3_lptmva1_badele1", 50,-10.,10.);
  TH1F *h_comb3_lptmva2_badele2 = new TH1F("h_comb3_lptmva2_badele2", "h_comb3_lptmva2_badele2", 50,-10.,10.);
  TH1F *hz_comb3_kpt    = new TH1F("hz_comb3_kpt",    "hz_comb3_kpt",    50, 0., 5.);
  TH1F *hz_comb3_ele1pt = new TH1F("hz_comb3_ele1pt", "hz_comb3_ele1pt", 50, 0., 10.);
  TH1F *hz_comb3_ele2pt = new TH1F("hz_comb3_ele2pt", "hz_comb3_ele2pt", 50, 0., 5.);
  TH1F *hz_comb3_kpt_badk       = new TH1F("hz_comb3_kpt_badk",       "hz_comb3_kpt_badk",        50, 0., 5.);
  TH1F *hz_comb3_ele1pt_badele1 = new TH1F("hz_comb3_ele1pt_badele1", "hz_comb3_ele1pt_badele1",  50, 0., 10.);
  TH1F *hz_comb3_ele2pt_badele2 = new TH1F("hz_comb3_ele2pt_badele2", "hz_comb3_ele2pt_badele2",  50, 0., 5.);
  TH1F *hz_comb3_minpt  = new TH1F("hz_comb3_minpt",  "hz_comb3_minpt",  50, 0., 5.);
  TH1F *h_comb3_notmatching = new TH1F("h_comb3_notmatching", "h_comb3_notmatching", 9, -1.5, 7.5);

  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    // multiplicity
    h_numB      -> Fill(goodBSize);
    h_numCombB  -> Fill(goodCombBSize);
    h_numTrueB  -> Fill(goodTrueBSize);

    // Best matched B related quantities
    if (bestMatch_SvProb>-999) {

      h_trueBestMatch_svProb   -> Fill(bestMatch_SvProb);
      h_trueBestMatch_xysig    -> Fill(bestMatch_XYSig);
      h_trueBestMatch_cos2d    -> Fill(bestMatch_Cos2D);
      h_trueBestMatch_kpt      -> Fill(bestMatch_KPt);
      h_trueBestMatch_keta     -> Fill(bestMatch_KEta);
      h_trueBestMatch_ele1eta  -> Fill(bestMatch_Ele1Eta);
      h_trueBestMatch_ele2eta  -> Fill(bestMatch_Ele2Eta);
      h_trueBestMatch_ele1pt   -> Fill(bestMatch_Ele1Pt);
      h_trueBestMatch_ele2pt   -> Fill(bestMatch_Ele2Pt);
      h_trueBestMatch_pfmva1   -> Fill(bestMatch_Ele1pfmva);
      h_trueBestMatch_pfmva2   -> Fill(bestMatch_Ele2pfmva);
      h_trueBestMatch_lptmva1  -> Fill(bestMatch_Ele1lptmva);
      h_trueBestMatch_lptmva2  -> Fill(bestMatch_Ele2lptmva);
      hz_trueBestMatch_kpt     -> Fill(bestMatch_KPt);
      hz_trueBestMatch_ele1pt  -> Fill(bestMatch_Ele1Pt);
      hz_trueBestMatch_ele2pt  -> Fill(bestMatch_Ele2Pt);
      hz_trueBestMatch_minpt   -> Fill(bestMatch_MinPt);
    }
  
    // Combinatorics
    int goodCombB_size = goodCombB_svProb->size();
    for (int ii=0; ii<goodCombB_size; ii++) {
      
      // all combinatorics
      h_comb_svProb   -> Fill(goodCombB_svProb->at(ii));
      h_comb_xysig    -> Fill(goodCombB_xySig->at(ii));
      h_comb_cos2d    -> Fill(goodCombB_cos2D->at(ii));
      h_comb_kpt      -> Fill(goodCombB_kpt->at(ii));
      h_comb_keta     -> Fill(goodCombB_keta->at(ii));
      h_comb_ele1eta  -> Fill(goodCombB_ele1eta->at(ii));
      h_comb_ele2eta  -> Fill(goodCombB_ele2eta->at(ii));
      h_comb_ele1pt   -> Fill(goodCombB_ele1pt->at(ii));
      h_comb_ele2pt   -> Fill(goodCombB_ele2pt->at(ii));
      h_comb_pfmva1   -> Fill(goodCombB_ele1pfmva->at(ii));
      h_comb_pfmva2   -> Fill(goodCombB_ele2pfmva->at(ii));
      h_comb_lptmva1  -> Fill(goodCombB_ele1lptmva->at(ii));
      h_comb_lptmva2  -> Fill(goodCombB_ele2lptmva->at(ii));
      hz_comb_ele1pt  -> Fill(goodCombB_ele1pt->at(ii));
      hz_comb_ele2pt  -> Fill(goodCombB_ele2pt->at(ii));
      hz_comb_kpt     -> Fill(goodCombB_kpt->at(ii));
      hz_comb_minpt   -> Fill(goodCombB_minpt->at(ii));
      h_comb_notmatching -> Fill(goodCombB_notmatching->at(ii));

      if (goodCombB_notmatching->at(ii)<=3) { // just one not-matched object
	h_comb1_svProb   -> Fill(goodCombB_svProb->at(ii));
	h_comb1_xysig    -> Fill(goodCombB_xySig->at(ii));
	h_comb1_cos2d    -> Fill(goodCombB_cos2D->at(ii));
	h_comb1_kpt      -> Fill(goodCombB_kpt->at(ii));
	h_comb1_keta     -> Fill(goodCombB_keta->at(ii));
	h_comb1_ele1eta  -> Fill(goodCombB_ele1eta->at(ii));
	h_comb1_ele2eta  -> Fill(goodCombB_ele2eta->at(ii));
	h_comb1_ele1pt   -> Fill(goodCombB_ele1pt->at(ii));
	h_comb1_ele2pt   -> Fill(goodCombB_ele2pt->at(ii));
	h_comb1_pfmva1   -> Fill(goodCombB_ele1pfmva->at(ii));
	h_comb1_pfmva2   -> Fill(goodCombB_ele2pfmva->at(ii));
	h_comb1_lptmva1  -> Fill(goodCombB_ele1lptmva->at(ii));
	h_comb1_lptmva2  -> Fill(goodCombB_ele2lptmva->at(ii));
	hz_comb1_ele1pt  -> Fill(goodCombB_ele1pt->at(ii));
	hz_comb1_ele2pt  -> Fill(goodCombB_ele2pt->at(ii));
	hz_comb1_kpt     -> Fill(goodCombB_kpt->at(ii));
	hz_comb1_minpt   -> Fill(goodCombB_minpt->at(ii));
	h_comb1_notmatching -> Fill(goodCombB_notmatching->at(ii));

      } else if (goodCombB_notmatching->at(ii)<=6) { // 2 not-matched objects   
	h_comb2_svProb   -> Fill(goodCombB_svProb->at(ii));
	h_comb2_xysig    -> Fill(goodCombB_xySig->at(ii));
	h_comb2_cos2d    -> Fill(goodCombB_cos2D->at(ii));
	h_comb2_kpt      -> Fill(goodCombB_kpt->at(ii));
	h_comb2_keta     -> Fill(goodCombB_keta->at(ii));
	h_comb2_ele1eta  -> Fill(goodCombB_ele1eta->at(ii));
	h_comb2_ele2eta  -> Fill(goodCombB_ele2eta->at(ii));
	h_comb2_ele1pt   -> Fill(goodCombB_ele1pt->at(ii));
	h_comb2_ele2pt   -> Fill(goodCombB_ele2pt->at(ii));
	h_comb2_pfmva1   -> Fill(goodCombB_ele1pfmva->at(ii));
	h_comb2_pfmva2   -> Fill(goodCombB_ele2pfmva->at(ii));
	h_comb2_lptmva1  -> Fill(goodCombB_ele1lptmva->at(ii));
	h_comb2_lptmva2  -> Fill(goodCombB_ele2lptmva->at(ii));
	hz_comb2_ele1pt  -> Fill(goodCombB_ele1pt->at(ii));
	hz_comb2_ele2pt  -> Fill(goodCombB_ele2pt->at(ii));
	hz_comb2_kpt     -> Fill(goodCombB_kpt->at(ii));
	hz_comb2_minpt   -> Fill(goodCombB_minpt->at(ii));
	h_comb2_notmatching -> Fill(goodCombB_notmatching->at(ii));

      } else if (goodCombB_notmatching->at(ii)==7) { // 3 not-matched objects   
	h_comb3_svProb   -> Fill(goodCombB_svProb->at(ii));
	h_comb3_xysig    -> Fill(goodCombB_xySig->at(ii));
	h_comb3_cos2d    -> Fill(goodCombB_cos2D->at(ii));
	h_comb3_kpt      -> Fill(goodCombB_kpt->at(ii));
	h_comb3_keta     -> Fill(goodCombB_keta->at(ii));
	h_comb3_ele1eta  -> Fill(goodCombB_ele1eta->at(ii));
	h_comb3_ele2eta  -> Fill(goodCombB_ele2eta->at(ii));
	h_comb3_ele1pt   -> Fill(goodCombB_ele1pt->at(ii));
	h_comb3_ele2pt   -> Fill(goodCombB_ele2pt->at(ii));
	h_comb3_pfmva1   -> Fill(goodCombB_ele1pfmva->at(ii));
	h_comb3_pfmva2   -> Fill(goodCombB_ele2pfmva->at(ii));
	h_comb3_lptmva1  -> Fill(goodCombB_ele1lptmva->at(ii));
	h_comb3_lptmva2  -> Fill(goodCombB_ele2lptmva->at(ii));
	hz_comb3_ele1pt  -> Fill(goodCombB_ele1pt->at(ii));
	hz_comb3_ele2pt  -> Fill(goodCombB_ele2pt->at(ii));
	hz_comb3_kpt     -> Fill(goodCombB_kpt->at(ii));
	hz_comb3_minpt   -> Fill(goodCombB_minpt->at(ii));
	h_comb3_notmatching -> Fill(goodCombB_notmatching->at(ii));
      }

      // all combinatorics
      if (goodCombB_notmatching->at(ii)==1 || goodCombB_notmatching->at(ii)==4 || goodCombB_notmatching->at(ii)==5 || goodCombB_notmatching->at(ii)==7) { // bad ele1
	h_comb_pfmva1_badele1  -> Fill(goodCombB_ele1pfmva->at(ii));
	h_comb_lptmva1_badele1 -> Fill(goodCombB_ele1lptmva->at(ii));
	hz_comb_ele1pt_badele1 -> Fill(goodCombB_ele1pt->at(ii));
      }
      if (goodCombB_notmatching->at(ii)==2 || goodCombB_notmatching->at(ii)==4 || goodCombB_notmatching->at(ii)==6 || goodCombB_notmatching->at(ii)==7) { // bad ele2
	h_comb_pfmva2_badele2  -> Fill(goodCombB_ele2pfmva->at(ii));
	h_comb_lptmva2_badele2 -> Fill(goodCombB_ele2lptmva->at(ii));
	hz_comb_ele2pt_badele2 -> Fill(goodCombB_ele2pt->at(ii));
      }
      if (goodCombB_notmatching->at(ii)==3 || goodCombB_notmatching->at(ii)==5 || goodCombB_notmatching->at(ii)==6 || goodCombB_notmatching->at(ii)==7) { // bad K
	hz_comb_kpt_badk       -> Fill(goodCombB_kpt->at(ii));
      }

      // just one not-matched object
      if (goodCombB_notmatching->at(ii)==1) {  // bad ele1
	h_comb1_pfmva1_badele1  -> Fill(goodCombB_ele1pfmva->at(ii));
	h_comb1_lptmva1_badele1 -> Fill(goodCombB_ele1lptmva->at(ii));
	hz_comb1_ele1pt_badele1 -> Fill(goodCombB_ele1pt->at(ii));
      }
      if (goodCombB_notmatching->at(ii)==2) {  // bad ele2
	h_comb1_pfmva2_badele2  -> Fill(goodCombB_ele2pfmva->at(ii));
	h_comb1_lptmva2_badele2 -> Fill(goodCombB_ele2lptmva->at(ii));
	hz_comb1_ele2pt_badele2 -> Fill(goodCombB_ele2pt->at(ii));
      }
      if (goodCombB_notmatching->at(ii)==3) {  // bad K
	hz_comb1_kpt_badk       -> Fill(goodCombB_kpt->at(ii));
      }

      // two not-matched objects
      if (goodCombB_notmatching->at(ii)==4 || goodCombB_notmatching->at(ii)==5) { // bad ele1
	h_comb2_pfmva1_badele1  -> Fill(goodCombB_ele1pfmva->at(ii));
	h_comb2_lptmva1_badele1 -> Fill(goodCombB_ele1lptmva->at(ii));
	hz_comb2_ele1pt_badele1 -> Fill(goodCombB_ele1pt->at(ii));
      }
      if (goodCombB_notmatching->at(ii)==4 || goodCombB_notmatching->at(ii)==6) { // bad ele2
	h_comb2_pfmva2_badele2  -> Fill(goodCombB_ele2pfmva->at(ii));
	h_comb2_lptmva2_badele2 -> Fill(goodCombB_ele2lptmva->at(ii));
	hz_comb2_ele2pt_badele2 -> Fill(goodCombB_ele2pt->at(ii));
      }
      if (goodCombB_notmatching->at(ii)==5 || goodCombB_notmatching->at(ii)==6) {
	hz_comb2_kpt_badk       -> Fill(goodCombB_kpt->at(ii));
      }

      // three not-matched objects
      if (goodCombB_notmatching->at(ii)==7) { 
	h_comb3_pfmva1_badele1  -> Fill(goodCombB_ele1pfmva->at(ii));
	h_comb3_lptmva1_badele1 -> Fill(goodCombB_ele1lptmva->at(ii));
	hz_comb3_ele1pt_badele1 -> Fill(goodCombB_ele1pt->at(ii));
	h_comb3_pfmva2_badele2  -> Fill(goodCombB_ele2pfmva->at(ii));
	h_comb3_lptmva2_badele2 -> Fill(goodCombB_ele2lptmva->at(ii));
	hz_comb3_ele2pt_badele2 -> Fill(goodCombB_ele2pt->at(ii));
	hz_comb3_kpt_badk       -> Fill(goodCombB_kpt->at(ii));
      }
    }



  } // Loop over entries


  // Plots
  gStyle->SetOptStat(0);

  TCanvas c1("c1","",1);
  h_numB -> SetLineColor(2);
  h_numB -> SetLineWidth(2);
  h_numB -> Draw();
  h_numB -> GetXaxis()->SetTitle("Number of selected Bs");
  h_numB -> SetTitle("");
  c1.SetLogy();
  c1.SaveAs("NumberOfBs.png");

  TCanvas c2("c2","",1);
  h_numTrueB -> SetLineColor(2);
  h_numTrueB -> SetLineWidth(2);
  h_numTrueB -> Draw();
  h_numTrueB -> GetXaxis()->SetTitle("Number of selected Bs matching MC-truth");
  h_numTrueB -> SetTitle("");
  c2.SetLogy();
  c2.SaveAs("NumberOfTrueBs.png");

  TCanvas c3("c3","",1);
  h_numCombB -> SetLineColor(2);
  h_numCombB -> SetLineWidth(2);
  h_numCombB -> Draw();
  h_numCombB -> GetXaxis()->SetTitle("Number of selected Bs not matching MC-truth");
  h_numCombB -> SetTitle("");
  c3.SetLogy();
  c3.SaveAs("NumberOfCombBs.png");


  // True, distributions
  h_trueBestMatch_svProb   -> SetLineColor(2);
  h_trueBestMatch_svProb   -> SetLineWidth(2);
  h_trueBestMatch_svProb   -> GetXaxis()->SetTitle("SV Fit probability");
  h_trueBestMatch_svProb   -> SetTitle("");
  h_trueBestMatch_xysig    -> SetLineColor(2);
  h_trueBestMatch_xysig    -> SetLineWidth(2);
  h_trueBestMatch_xysig    -> GetXaxis()->SetTitle("dxy significance");
  h_trueBestMatch_xysig    -> SetTitle("");
  h_trueBestMatch_cos2d    -> SetLineColor(2);
  h_trueBestMatch_cos2d    -> SetLineWidth(2);
  h_trueBestMatch_cos2d    -> GetXaxis()->SetTitle("cos2D");
  h_trueBestMatch_cos2d    -> SetTitle("");
  h_trueBestMatch_kpt      -> SetLineColor(2);
  h_trueBestMatch_kpt      -> SetLineWidth(2);
  h_trueBestMatch_kpt      -> GetXaxis()->SetTitle("K pT");
  h_trueBestMatch_kpt      -> SetTitle("");
  h_trueBestMatch_keta     -> SetLineColor(2);
  h_trueBestMatch_keta     -> SetLineWidth(2);
  h_trueBestMatch_keta     -> GetXaxis()->SetTitle("K #eta");
  h_trueBestMatch_keta     -> SetTitle("");
  h_trueBestMatch_ele1eta  -> SetLineColor(2);
  h_trueBestMatch_ele1eta  -> SetLineWidth(2);
  h_trueBestMatch_ele1eta  -> GetXaxis()->SetTitle("ele1 #eta");
  h_trueBestMatch_ele1eta  -> SetTitle("");
  h_trueBestMatch_ele2eta  -> SetLineColor(2);
  h_trueBestMatch_ele2eta  -> SetLineWidth(2);
  h_trueBestMatch_ele2eta  -> GetXaxis()->SetTitle("ele2 #eta");
  h_trueBestMatch_ele2eta  -> SetTitle("");
  h_trueBestMatch_ele1pt   -> SetLineColor(2);
  h_trueBestMatch_ele1pt   -> SetLineWidth(2);
  h_trueBestMatch_ele1pt   -> GetXaxis()->SetTitle("ele1 pt");
  h_trueBestMatch_ele1pt   -> SetTitle("");
  h_trueBestMatch_ele2pt   -> SetLineColor(2);
  h_trueBestMatch_ele2pt   -> SetLineWidth(2);
  h_trueBestMatch_ele2pt   -> GetXaxis()->SetTitle("ele2 pt");
  h_trueBestMatch_ele2pt   -> SetTitle("");
  h_trueBestMatch_pfmva1   -> SetLineColor(2);
  h_trueBestMatch_pfmva1   -> SetLineWidth(2);
  h_trueBestMatch_pfmva1   -> GetXaxis()->SetTitle("ele1 PF mva");
  h_trueBestMatch_pfmva1   -> SetTitle("");
  h_trueBestMatch_pfmva2   -> SetLineColor(2);
  h_trueBestMatch_pfmva2   -> SetLineWidth(2);
  h_trueBestMatch_pfmva2   -> GetXaxis()->SetTitle("ele1 PF mva");
  h_trueBestMatch_pfmva2   -> SetTitle("");
  h_trueBestMatch_lptmva1  -> SetLineColor(2);
  h_trueBestMatch_lptmva1  -> SetLineWidth(2);
  h_trueBestMatch_lptmva1  -> GetXaxis()->SetTitle("ele1 LPT mva");
  h_trueBestMatch_lptmva1  -> SetTitle("");
  h_trueBestMatch_lptmva2  -> SetLineColor(2);
  h_trueBestMatch_lptmva2  -> SetLineWidth(2);
  h_trueBestMatch_lptmva2  -> GetXaxis()->SetTitle("ele2 LPT mva");
  h_trueBestMatch_lptmva2  -> SetTitle("");
  hz_trueBestMatch_kpt      -> SetLineColor(2);
  hz_trueBestMatch_kpt      -> SetLineWidth(2);
  hz_trueBestMatch_kpt      -> GetXaxis()->SetTitle("K pT");
  hz_trueBestMatch_kpt      -> SetTitle("");
  hz_trueBestMatch_ele1pt   -> SetLineColor(2);
  hz_trueBestMatch_ele1pt   -> SetLineWidth(2);
  hz_trueBestMatch_ele1pt   -> GetXaxis()->SetTitle("ele1 pt");
  hz_trueBestMatch_ele1pt   -> SetTitle("");
  hz_trueBestMatch_ele2pt   -> SetLineColor(2);
  hz_trueBestMatch_ele2pt   -> SetLineWidth(2);
  hz_trueBestMatch_ele2pt   -> GetXaxis()->SetTitle("ele2 pt");
  hz_trueBestMatch_ele2pt   -> SetTitle("");
  hz_trueBestMatch_minpt   -> SetLineColor(2);
  hz_trueBestMatch_minpt   -> SetLineWidth(2);
  hz_trueBestMatch_minpt   -> GetXaxis()->SetTitle("min(K, ele2) pt");
  hz_trueBestMatch_minpt   -> SetTitle("");
  //
  // Comb, distributions
  h_comb_svProb   -> SetLineColor(4);
  h_comb_svProb   -> SetLineWidth(2);
  h_comb_svProb   -> GetXaxis()->SetTitle("SV Fit probability");
  h_comb_svProb   -> SetTitle("");
  h_comb_xysig    -> SetLineColor(4);
  h_comb_xysig    -> SetLineWidth(2);
  h_comb_xysig    -> GetXaxis()->SetTitle("dxy significance");
  h_comb_xysig    -> SetTitle("");
  h_comb_cos2d    -> SetLineColor(4);
  h_comb_cos2d    -> SetLineWidth(2);
  h_comb_cos2d    -> GetXaxis()->SetTitle("cos2D");
  h_comb_cos2d    -> SetTitle("");
  h_comb_kpt      -> SetLineColor(4);
  h_comb_kpt      -> SetLineWidth(2);
  h_comb_kpt      -> GetXaxis()->SetTitle("K pT");
  h_comb_kpt      -> SetTitle("");
  h_comb_keta      -> SetLineColor(4);
  h_comb_keta      -> SetLineWidth(2);
  h_comb_keta      -> GetXaxis()->SetTitle("K #eta");
  h_comb_keta      -> SetTitle("");
  h_comb_ele1eta      -> SetLineColor(4);
  h_comb_ele1eta      -> SetLineWidth(2);
  h_comb_ele1eta      -> GetXaxis()->SetTitle("ele1 #eta");
  h_comb_ele1eta      -> SetTitle("");
  h_comb_ele2eta      -> SetLineColor(4);
  h_comb_ele2eta      -> SetLineWidth(2);
  h_comb_ele2eta      -> GetXaxis()->SetTitle("ele2 #eta");
  h_comb_ele2eta      -> SetTitle("");
  h_comb_ele1pt      -> SetLineColor(4);
  h_comb_ele1pt      -> SetLineWidth(2);
  h_comb_ele1pt      -> GetXaxis()->SetTitle("ele1 pT");
  h_comb_ele1pt      -> SetTitle("");
  h_comb_ele2pt      -> SetLineColor(4);
  h_comb_ele2pt      -> SetLineWidth(2);
  h_comb_ele2pt      -> GetXaxis()->SetTitle("ele2 pT");
  h_comb_ele2pt      -> SetTitle("");
  h_comb_pfmva1 -> SetLineColor(4);
  h_comb_pfmva1 -> SetLineWidth(2);
  h_comb_pfmva1 -> GetXaxis()->SetTitle("ele1 PF mva");
  h_comb_pfmva1 -> SetTitle("");
  h_comb_pfmva1_badele1 -> SetLineColor(4);
  h_comb_pfmva1_badele1 -> SetLineWidth(2);
  h_comb_pfmva1_badele1 -> GetXaxis()->SetTitle("ele1 PF mva");
  h_comb_pfmva1_badele1 -> SetTitle("");
  h_comb_pfmva2 -> SetLineColor(4);
  h_comb_pfmva2 -> SetLineWidth(2);
  h_comb_pfmva2 -> GetXaxis()->SetTitle("ele2 PF mva");
  h_comb_pfmva2 -> SetTitle("");
  h_comb_pfmva2_badele2 -> SetLineColor(4);
  h_comb_pfmva2_badele2 -> SetLineWidth(2);
  h_comb_pfmva2_badele2 -> GetXaxis()->SetTitle("ele1 PF mva");
  h_comb_pfmva2_badele2 -> SetTitle("");
  h_comb_lptmva1 -> SetLineColor(4);
  h_comb_lptmva1 -> SetLineWidth(2);
  h_comb_lptmva1 -> GetXaxis()->SetTitle("ele1 LPT mva");
  h_comb_lptmva1 -> SetTitle("");
  h_comb_lptmva1_badele1 -> SetLineColor(4);
  h_comb_lptmva1_badele1 -> SetLineWidth(2);
  h_comb_lptmva1_badele1 -> GetXaxis()->SetTitle("ele1 LPT mva");
  h_comb_lptmva1_badele1 -> SetTitle("");
  h_comb_lptmva2 -> SetLineColor(4);
  h_comb_lptmva2 -> SetLineWidth(2);
  h_comb_lptmva2 -> GetXaxis()->SetTitle("ele2 LPT mva");
  h_comb_lptmva2 -> SetTitle("");
  h_comb_lptmva2_badele2 -> SetLineColor(4);
  h_comb_lptmva2_badele2 -> SetLineWidth(2);
  h_comb_lptmva2_badele2 -> GetXaxis()->SetTitle("ele2 LPT mva");
  h_comb_lptmva2_badele2 -> SetTitle("");
  hz_comb_kpt      -> SetLineColor(4);
  hz_comb_kpt      -> SetLineWidth(2);
  hz_comb_kpt      -> GetXaxis()->SetTitle("K pT");
  hz_comb_kpt      -> SetTitle("");
  hz_comb_ele1pt      -> SetLineColor(4);
  hz_comb_ele1pt      -> SetLineWidth(2);
  hz_comb_ele1pt      -> GetXaxis()->SetTitle("ele1 pT");
  hz_comb_ele1pt      -> SetTitle("");
  hz_comb_ele2pt      -> SetLineColor(4);
  hz_comb_ele2pt      -> SetLineWidth(2);
  hz_comb_ele2pt      -> GetXaxis()->SetTitle("ele2 pT");
  hz_comb_ele2pt      -> SetTitle("");
  hz_comb_kpt_badk       -> SetLineColor(4);
  hz_comb_kpt_badk       -> SetLineWidth(2);
  hz_comb_kpt_badk       -> GetXaxis()->SetTitle("K pT");
  hz_comb_kpt_badk       -> SetTitle("");
  hz_comb_ele1pt_badele1 -> SetLineColor(4);
  hz_comb_ele1pt_badele1 -> SetLineWidth(2);
  hz_comb_ele1pt_badele1 -> GetXaxis()->SetTitle("ele1 pT");
  hz_comb_ele1pt_badele1 -> SetTitle("");
  hz_comb_ele2pt_badele2 -> SetLineColor(4);
  hz_comb_ele2pt_badele2 -> SetLineWidth(2);
  hz_comb_ele2pt_badele2 -> GetXaxis()->SetTitle("ele2 pT");
  hz_comb_ele2pt_badele2 -> SetTitle("");
  hz_comb_minpt       -> SetLineColor(4);
  hz_comb_minpt       -> SetLineWidth(2);
  hz_comb_minpt       -> GetXaxis()->SetTitle("min (ele2, K) pT");
  hz_comb_minpt       -> SetTitle("");
  h_comb_notmatching  -> SetLineColor(4);
  h_comb_notmatching  -> SetLineWidth(2);
  h_comb_notmatching  -> GetXaxis()->SetTitle("not matching");
  h_comb_notmatching  -> SetTitle("");;
  //
  h_comb1_svProb   -> SetLineColor(6);
  h_comb1_svProb   -> SetLineWidth(2);
  h_comb1_svProb   -> GetXaxis()->SetTitle("SV Fit probability");
  h_comb1_svProb   -> SetTitle("");
  h_comb1_xysig    -> SetLineColor(6);
  h_comb1_xysig    -> SetLineWidth(2);
  h_comb1_xysig    -> GetXaxis()->SetTitle("dxy significance");
  h_comb1_xysig    -> SetTitle("");
  h_comb1_cos2d    -> SetLineColor(6);
  h_comb1_cos2d    -> SetLineWidth(2);
  h_comb1_cos2d    -> GetXaxis()->SetTitle("cos2D");
  h_comb1_cos2d    -> SetTitle("");
  h_comb1_kpt      -> SetLineColor(6);
  h_comb1_kpt      -> SetLineWidth(2);
  h_comb1_kpt      -> GetXaxis()->SetTitle("K pT");
  h_comb1_kpt      -> SetTitle("");
  h_comb1_keta      -> SetLineColor(6);
  h_comb1_keta      -> SetLineWidth(2);
  h_comb1_keta      -> GetXaxis()->SetTitle("K #eta");
  h_comb1_keta      -> SetTitle("");
  h_comb1_ele1eta      -> SetLineColor(6);
  h_comb1_ele1eta      -> SetLineWidth(2);
  h_comb1_ele1eta      -> GetXaxis()->SetTitle("ele1 #eta");
  h_comb1_ele1eta      -> SetTitle("");
  h_comb1_ele2eta      -> SetLineColor(6);
  h_comb1_ele2eta      -> SetLineWidth(2);
  h_comb1_ele2eta      -> GetXaxis()->SetTitle("ele2 #eta");
  h_comb1_ele2eta      -> SetTitle("");
  h_comb1_ele1pt      -> SetLineColor(6);
  h_comb1_ele1pt      -> SetLineWidth(2);
  h_comb1_ele1pt      -> GetXaxis()->SetTitle("ele1 pT");
  h_comb1_ele1pt      -> SetTitle("");
  h_comb1_ele2pt      -> SetLineColor(6);
  h_comb1_ele2pt      -> SetLineWidth(2);
  h_comb1_ele2pt      -> GetXaxis()->SetTitle("ele2 pT");
  h_comb1_ele2pt      -> SetTitle("");
  h_comb1_pfmva1 -> SetLineColor(6);
  h_comb1_pfmva1 -> SetLineWidth(2);
  h_comb1_pfmva1 -> GetXaxis()->SetTitle("ele1 PF mva");
  h_comb1_pfmva1 -> SetTitle("");
  h_comb1_pfmva1_badele1 -> SetLineColor(6);
  h_comb1_pfmva1_badele1 -> SetLineWidth(2);
  h_comb1_pfmva1_badele1 -> GetXaxis()->SetTitle("ele1 PF mva");
  h_comb1_pfmva1_badele1 -> SetTitle("");
  h_comb1_pfmva2 -> SetLineColor(6);
  h_comb1_pfmva2 -> SetLineWidth(2);
  h_comb1_pfmva2 -> GetXaxis()->SetTitle("ele2 PF mva");
  h_comb1_pfmva2 -> SetTitle("");
  h_comb1_pfmva2_badele2 -> SetLineColor(6);
  h_comb1_pfmva2_badele2 -> SetLineWidth(2);
  h_comb1_pfmva2_badele2 -> GetXaxis()->SetTitle("ele1 PF mva");
  h_comb1_pfmva2_badele2 -> SetTitle("");
  h_comb1_lptmva1 -> SetLineColor(6);
  h_comb1_lptmva1 -> SetLineWidth(2);
  h_comb1_lptmva1 -> GetXaxis()->SetTitle("ele1 LPT mva");
  h_comb1_lptmva1 -> SetTitle("");
  h_comb1_lptmva1_badele1 -> SetLineColor(6);
  h_comb1_lptmva1_badele1 -> SetLineWidth(2);
  h_comb1_lptmva1_badele1 -> GetXaxis()->SetTitle("ele1 LPT mva");
  h_comb1_lptmva1_badele1 -> SetTitle("");
  h_comb1_lptmva2 -> SetLineColor(6);
  h_comb1_lptmva2 -> SetLineWidth(2);
  h_comb1_lptmva2 -> GetXaxis()->SetTitle("ele2 LPT mva");
  h_comb1_lptmva2 -> SetTitle("");
  h_comb1_lptmva2_badele2 -> SetLineColor(6);
  h_comb1_lptmva2_badele2 -> SetLineWidth(2);
  h_comb1_lptmva2_badele2 -> GetXaxis()->SetTitle("ele2 LPT mva");
  h_comb1_lptmva2_badele2 -> SetTitle("");
  hz_comb1_kpt      -> SetLineColor(6);
  hz_comb1_kpt      -> SetLineWidth(2);
  hz_comb1_kpt      -> GetXaxis()->SetTitle("K pT");
  hz_comb1_kpt      -> SetTitle("");
  hz_comb1_ele1pt      -> SetLineColor(6);
  hz_comb1_ele1pt      -> SetLineWidth(2);
  hz_comb1_ele1pt      -> GetXaxis()->SetTitle("ele1 pT");
  hz_comb1_ele1pt      -> SetTitle("");
  hz_comb1_ele2pt      -> SetLineColor(6);
  hz_comb1_ele2pt      -> SetLineWidth(2);
  hz_comb1_ele2pt      -> GetXaxis()->SetTitle("ele2 pT");
  hz_comb1_ele2pt      -> SetTitle("");
  hz_comb1_kpt_badk       -> SetLineColor(6);
  hz_comb1_kpt_badk       -> SetLineWidth(2);
  hz_comb1_kpt_badk       -> GetXaxis()->SetTitle("K pT");
  hz_comb1_kpt_badk       -> SetTitle("");
  hz_comb1_ele1pt_badele1 -> SetLineColor(6);
  hz_comb1_ele1pt_badele1 -> SetLineWidth(2);
  hz_comb1_ele1pt_badele1 -> GetXaxis()->SetTitle("ele1 pT");
  hz_comb1_ele1pt_badele1 -> SetTitle("");
  hz_comb1_ele2pt_badele2 -> SetLineColor(6);
  hz_comb1_ele2pt_badele2 -> SetLineWidth(2);
  hz_comb1_ele2pt_badele2 -> GetXaxis()->SetTitle("ele2 pT");
  hz_comb1_ele2pt_badele2 -> SetTitle("");
  hz_comb1_minpt       -> SetLineColor(6);
  hz_comb1_minpt       -> SetLineWidth(2);
  hz_comb1_minpt       -> GetXaxis()->SetTitle("min (ele2, K) pT");
  hz_comb1_minpt       -> SetTitle("");
  h_comb1_notmatching  -> SetLineColor(6);
  h_comb1_notmatching  -> SetLineWidth(2);
  h_comb1_notmatching  -> GetXaxis()->SetTitle("not matching");
  h_comb1_notmatching  -> SetTitle("");;
  //
  h_comb2_svProb   -> SetLineColor(3);
  h_comb2_svProb   -> SetLineWidth(2);
  h_comb2_svProb   -> GetXaxis()->SetTitle("SV Fit probability");
  h_comb2_svProb   -> SetTitle("");
  h_comb2_xysig    -> SetLineColor(3);
  h_comb2_xysig    -> SetLineWidth(2);
  h_comb2_xysig    -> GetXaxis()->SetTitle("dxy significance");
  h_comb2_xysig    -> SetTitle("");
  h_comb2_cos2d    -> SetLineColor(3);
  h_comb2_cos2d    -> SetLineWidth(2);
  h_comb2_cos2d    -> GetXaxis()->SetTitle("cos2D");
  h_comb2_cos2d    -> SetTitle("");
  h_comb2_kpt      -> SetLineColor(3);
  h_comb2_kpt      -> SetLineWidth(2);
  h_comb2_kpt      -> GetXaxis()->SetTitle("K pT");
  h_comb2_kpt      -> SetTitle("");
  h_comb2_keta      -> SetLineColor(3);
  h_comb2_keta      -> SetLineWidth(2);
  h_comb2_keta      -> GetXaxis()->SetTitle("K #eta");
  h_comb2_keta      -> SetTitle("");
  h_comb2_ele1eta      -> SetLineColor(3);
  h_comb2_ele1eta      -> SetLineWidth(2);
  h_comb2_ele1eta      -> GetXaxis()->SetTitle("ele1 #eta");
  h_comb2_ele1eta      -> SetTitle("");
  h_comb2_ele2eta      -> SetLineColor(3);
  h_comb2_ele2eta      -> SetLineWidth(2);
  h_comb2_ele2eta      -> GetXaxis()->SetTitle("ele2 #eta");
  h_comb2_ele2eta      -> SetTitle("");
  h_comb2_ele1pt      -> SetLineColor(3);
  h_comb2_ele1pt      -> SetLineWidth(2);
  h_comb2_ele1pt      -> GetXaxis()->SetTitle("ele1 pT");
  h_comb2_ele1pt      -> SetTitle("");
  h_comb2_ele2pt      -> SetLineColor(3);
  h_comb2_ele2pt      -> SetLineWidth(2);
  h_comb2_ele2pt      -> GetXaxis()->SetTitle("ele2 pT");
  h_comb2_ele2pt      -> SetTitle("");
  h_comb2_pfmva1 -> SetLineColor(3);
  h_comb2_pfmva1 -> SetLineWidth(2);
  h_comb2_pfmva1 -> GetXaxis()->SetTitle("ele1 PF mva");
  h_comb2_pfmva1 -> SetTitle("");
  h_comb2_pfmva1_badele1 -> SetLineColor(3);
  h_comb2_pfmva1_badele1 -> SetLineWidth(2);
  h_comb2_pfmva1_badele1 -> GetXaxis()->SetTitle("ele1 PF mva");
  h_comb2_pfmva1_badele1 -> SetTitle("");
  h_comb2_pfmva2 -> SetLineColor(3);
  h_comb2_pfmva2 -> SetLineWidth(2);
  h_comb2_pfmva2 -> GetXaxis()->SetTitle("ele2 PF mva");
  h_comb2_pfmva2 -> SetTitle("");
  h_comb2_pfmva2_badele2 -> SetLineColor(3);
  h_comb2_pfmva2_badele2 -> SetLineWidth(2);
  h_comb2_pfmva2_badele2 -> GetXaxis()->SetTitle("ele1 PF mva");
  h_comb2_pfmva2_badele2 -> SetTitle("");
  h_comb2_lptmva1 -> SetLineColor(3);
  h_comb2_lptmva1 -> SetLineWidth(2);
  h_comb2_lptmva1 -> GetXaxis()->SetTitle("ele1 LPT mva");
  h_comb2_lptmva1 -> SetTitle("");
  h_comb2_lptmva1_badele1 -> SetLineColor(3);
  h_comb2_lptmva1_badele1 -> SetLineWidth(2);
  h_comb2_lptmva1_badele1 -> GetXaxis()->SetTitle("ele1 LPT mva");
  h_comb2_lptmva1_badele1 -> SetTitle("");
  h_comb2_lptmva2 -> SetLineColor(3);
  h_comb2_lptmva2 -> SetLineWidth(2);
  h_comb2_lptmva2 -> GetXaxis()->SetTitle("ele2 LPT mva");
  h_comb2_lptmva2 -> SetTitle("");
  h_comb2_lptmva2_badele2 -> SetLineColor(3);
  h_comb2_lptmva2_badele2 -> SetLineWidth(2);
  h_comb2_lptmva2_badele2 -> GetXaxis()->SetTitle("ele2 LPT mva");
  h_comb2_lptmva2_badele2 -> SetTitle("");
  hz_comb2_kpt      -> SetLineColor(3);
  hz_comb2_kpt      -> SetLineWidth(2);
  hz_comb2_kpt      -> GetXaxis()->SetTitle("K pT");
  hz_comb2_kpt      -> SetTitle("");
  hz_comb2_ele1pt      -> SetLineColor(3);
  hz_comb2_ele1pt      -> SetLineWidth(2);
  hz_comb2_ele1pt      -> GetXaxis()->SetTitle("ele1 pT");
  hz_comb2_ele1pt      -> SetTitle("");
  hz_comb2_ele2pt      -> SetLineColor(3);
  hz_comb2_ele2pt      -> SetLineWidth(2);
  hz_comb2_ele2pt      -> GetXaxis()->SetTitle("ele2 pT");
  hz_comb2_ele2pt      -> SetTitle("");
  hz_comb2_kpt_badk       -> SetLineColor(3);
  hz_comb2_kpt_badk       -> SetLineWidth(2);
  hz_comb2_kpt_badk       -> GetXaxis()->SetTitle("K pT");
  hz_comb2_kpt_badk       -> SetTitle("");
  hz_comb2_ele1pt_badele1 -> SetLineColor(3);
  hz_comb2_ele1pt_badele1 -> SetLineWidth(2);
  hz_comb2_ele1pt_badele1 -> GetXaxis()->SetTitle("ele1 pT");
  hz_comb2_ele1pt_badele1 -> SetTitle("");
  hz_comb2_ele2pt_badele2 -> SetLineColor(3);
  hz_comb2_ele2pt_badele2 -> SetLineWidth(2);
  hz_comb2_ele2pt_badele2 -> GetXaxis()->SetTitle("ele2 pT");
  hz_comb2_ele2pt_badele2 -> SetTitle("");
  hz_comb2_minpt       -> SetLineColor(3);
  hz_comb2_minpt       -> SetLineWidth(2);
  hz_comb2_minpt       -> GetXaxis()->SetTitle("min (ele2, K) pT");
  hz_comb2_minpt       -> SetTitle("");
  h_comb2_notmatching  -> SetLineColor(3);
  h_comb2_notmatching  -> SetLineWidth(2);
  h_comb2_notmatching  -> GetXaxis()->SetTitle("not matching");
  h_comb2_notmatching  -> SetTitle("");;
  //
  h_comb3_svProb   -> SetLineColor(1);
  h_comb3_svProb   -> SetLineWidth(2);
  h_comb3_svProb   -> GetXaxis()->SetTitle("SV Fit probability");
  h_comb3_svProb   -> SetTitle("");
  h_comb3_xysig    -> SetLineColor(1);
  h_comb3_xysig    -> SetLineWidth(2);
  h_comb3_xysig    -> GetXaxis()->SetTitle("dxy significance");
  h_comb3_xysig    -> SetTitle("");
  h_comb3_cos2d    -> SetLineColor(1);
  h_comb3_cos2d    -> SetLineWidth(2);
  h_comb3_cos2d    -> GetXaxis()->SetTitle("cos2D");
  h_comb3_cos2d    -> SetTitle("");
  h_comb3_kpt      -> SetLineColor(1);
  h_comb3_kpt      -> SetLineWidth(2);
  h_comb3_kpt      -> GetXaxis()->SetTitle("K pT");
  h_comb3_kpt      -> SetTitle("");
  h_comb3_keta      -> SetLineColor(1);
  h_comb3_keta      -> SetLineWidth(2);
  h_comb3_keta      -> GetXaxis()->SetTitle("K #eta");
  h_comb3_keta      -> SetTitle("");
  h_comb3_ele1eta      -> SetLineColor(1);
  h_comb3_ele1eta      -> SetLineWidth(2);
  h_comb3_ele1eta      -> GetXaxis()->SetTitle("ele1 #eta");
  h_comb3_ele1eta      -> SetTitle("");
  h_comb3_ele2eta      -> SetLineColor(1);
  h_comb3_ele2eta      -> SetLineWidth(2);
  h_comb3_ele2eta      -> GetXaxis()->SetTitle("ele2 #eta");
  h_comb3_ele2eta      -> SetTitle("");
  h_comb3_ele1pt      -> SetLineColor(1);
  h_comb3_ele1pt      -> SetLineWidth(2);
  h_comb3_ele1pt      -> GetXaxis()->SetTitle("ele1 pT");
  h_comb3_ele1pt      -> SetTitle("");
  h_comb3_ele2pt      -> SetLineColor(1);
  h_comb3_ele2pt      -> SetLineWidth(2);
  h_comb3_ele2pt      -> GetXaxis()->SetTitle("ele2 pT");
  h_comb3_ele2pt      -> SetTitle("");
  h_comb3_pfmva1 -> SetLineColor(1);
  h_comb3_pfmva1 -> SetLineWidth(2);
  h_comb3_pfmva1 -> GetXaxis()->SetTitle("ele1 PF mva");
  h_comb3_pfmva1 -> SetTitle("");
  h_comb3_pfmva1_badele1 -> SetLineColor(1);
  h_comb3_pfmva1_badele1 -> SetLineWidth(2);
  h_comb3_pfmva1_badele1 -> GetXaxis()->SetTitle("ele1 PF mva");
  h_comb3_pfmva1_badele1 -> SetTitle("");
  h_comb3_pfmva2 -> SetLineColor(1);
  h_comb3_pfmva2 -> SetLineWidth(2);
  h_comb3_pfmva2 -> GetXaxis()->SetTitle("ele2 PF mva");
  h_comb3_pfmva2 -> SetTitle("");
  h_comb3_pfmva2_badele2 -> SetLineColor(1);
  h_comb3_pfmva2_badele2 -> SetLineWidth(2);
  h_comb3_pfmva2_badele2 -> GetXaxis()->SetTitle("ele1 PF mva");
  h_comb3_pfmva2_badele2 -> SetTitle("");
  h_comb3_lptmva1 -> SetLineColor(1);
  h_comb3_lptmva1 -> SetLineWidth(2);
  h_comb3_lptmva1 -> GetXaxis()->SetTitle("ele1 LPT mva");
  h_comb3_lptmva1 -> SetTitle("");
  h_comb3_lptmva1_badele1 -> SetLineColor(1);
  h_comb3_lptmva1_badele1 -> SetLineWidth(2);
  h_comb3_lptmva1_badele1 -> GetXaxis()->SetTitle("ele1 LPT mva");
  h_comb3_lptmva1_badele1 -> SetTitle("");
  h_comb3_lptmva2 -> SetLineColor(1);
  h_comb3_lptmva2 -> SetLineWidth(2);
  h_comb3_lptmva2 -> GetXaxis()->SetTitle("ele2 LPT mva");
  h_comb3_lptmva2 -> SetTitle("");
  h_comb3_lptmva2_badele2 -> SetLineColor(1);
  h_comb3_lptmva2_badele2 -> SetLineWidth(2);
  h_comb3_lptmva2_badele2 -> GetXaxis()->SetTitle("ele2 LPT mva");
  h_comb3_lptmva2_badele2 -> SetTitle("");
  hz_comb3_kpt      -> SetLineColor(1);
  hz_comb3_kpt      -> SetLineWidth(2);
  hz_comb3_kpt      -> GetXaxis()->SetTitle("K pT");
  hz_comb3_kpt      -> SetTitle("");
  hz_comb3_ele1pt      -> SetLineColor(1);
  hz_comb3_ele1pt      -> SetLineWidth(2);
  hz_comb3_ele1pt      -> GetXaxis()->SetTitle("ele1 pT");
  hz_comb3_ele1pt      -> SetTitle("");
  hz_comb3_ele2pt      -> SetLineColor(1);
  hz_comb3_ele2pt      -> SetLineWidth(2);
  hz_comb3_ele2pt      -> GetXaxis()->SetTitle("ele2 pT");
  hz_comb3_ele2pt      -> SetTitle("");
  hz_comb3_kpt_badk       -> SetLineColor(1);
  hz_comb3_kpt_badk       -> SetLineWidth(2);
  hz_comb3_kpt_badk       -> GetXaxis()->SetTitle("K pT");
  hz_comb3_kpt_badk       -> SetTitle("");
  hz_comb3_ele1pt_badele1 -> SetLineColor(1);
  hz_comb3_ele1pt_badele1 -> SetLineWidth(2);
  hz_comb3_ele1pt_badele1 -> GetXaxis()->SetTitle("ele1 pT");
  hz_comb3_ele1pt_badele1 -> SetTitle("");
  hz_comb3_ele2pt_badele2 -> SetLineColor(1);
  hz_comb3_ele2pt_badele2 -> SetLineWidth(2);
  hz_comb3_ele2pt_badele2 -> GetXaxis()->SetTitle("ele2 pT");
  hz_comb3_ele2pt_badele2 -> SetTitle("");
  hz_comb3_minpt       -> SetLineColor(1); 
  hz_comb3_minpt       -> SetLineWidth(2);
  hz_comb3_minpt       -> GetXaxis()->SetTitle("min (ele2, K) pT");
  hz_comb3_minpt       -> SetTitle("");
  h_comb3_notmatching  -> SetLineColor(1);
  h_comb3_notmatching  -> SetLineWidth(2);
  h_comb3_notmatching  -> GetXaxis()->SetTitle("not matching");
  h_comb3_notmatching  -> SetTitle("");;
  //
  TLegend leg (0.55,0.7,0.95,0.9);
  leg.SetFillColor(0);
  leg.SetFillStyle(0);
  leg.SetBorderSize(0);
  leg.AddEntry(h_trueBestMatch_svProb,"Best match to MC-truth");
  leg.AddEntry(h_comb_svProb,"Combinatorics, all");
  leg.AddEntry(h_comb1_svProb,"Combinatorics, 1 not matched");
  leg.AddEntry(h_comb2_svProb,"Combinatorics, 2 not matched");
  leg.AddEntry(h_comb3_svProb,"Combinatorics, 3 not matched");
  //
  TLegend leg2 (0.15,0.7,0.55,0.9);
  leg2.SetFillColor(0);
  leg2.SetFillStyle(0);
  leg2.SetBorderSize(0);
  leg2.AddEntry(h_trueBestMatch_svProb,"Best match to MC-truth");
  leg2.AddEntry(h_comb_svProb,"Combinatorics, all");
  leg2.AddEntry(h_comb1_svProb,"Combinatorics, 1 not matched");
  leg2.AddEntry(h_comb2_svProb,"Combinatorics, 2 not matched");
  leg2.AddEntry(h_comb3_svProb,"Combinatorics, 3 not matched");
  //
  //
  TCanvas c0a("c0a","",1);
  h_comb_notmatching->DrawNormalized();
  c0a.SaveAs("Combinatorics_all.png");
  //
  TCanvas c0b("c0b","",1);
  h_comb1_notmatching->DrawNormalized();
  c0b.SaveAs("Combinatorics_1.png");
  //
  TCanvas c0c("c0c","",1);
  h_comb2_notmatching->DrawNormalized();
  c0c.SaveAs("Combinatorics_2.png");
  //
  TCanvas c0d("c0d","",1);
  h_comb3_notmatching->DrawNormalized();
  c0d.SaveAs("Combinatorics_3.png");
  //
  //
  TCanvas c10("c10","",1);
  h_comb_svProb->DrawNormalized();
  h_comb1_svProb->DrawNormalized("same");
  h_comb2_svProb->DrawNormalized("same");
  h_comb3_svProb->DrawNormalized("same");
  h_trueBestMatch_svProb->DrawNormalized("same");
  leg.Draw("same");  
  c10.SaveAs("TrueVsComb_svprob_split.png");
  //
  TCanvas c11("c11","",1);
  h_comb_xysig->DrawNormalized();
  h_comb1_xysig->DrawNormalized("same");
  h_comb2_xysig->DrawNormalized("same");
  h_comb3_xysig->DrawNormalized("same");
  h_trueBestMatch_xysig->DrawNormalized("same");
  leg.Draw("same");  
  c11.SetLogy();
  c11.SaveAs("TrueVsComb_dxy_split.png");
  //
  TCanvas c12("c12","",1);
  h_trueBestMatch_cos2d->DrawNormalized();
  h_comb_cos2d->DrawNormalized("same");
  h_comb1_cos2d->DrawNormalized("same");
  h_comb2_cos2d->DrawNormalized("same");
  h_comb3_cos2d->DrawNormalized("same");
  leg.Draw("same");  
  c12.SetLogy();
  c12.SaveAs("TrueVsComb_cos2d_split.png");
  //
  TCanvas c14("c14","",1);
  h_comb_kpt->DrawNormalized();
  h_comb1_kpt->DrawNormalized("same");
  h_comb2_kpt->DrawNormalized("same");
  h_comb3_kpt->DrawNormalized("same");
  h_trueBestMatch_kpt->DrawNormalized("same");
  leg.Draw("same");  
  c14.SetLogy();
  c14.SaveAs("TrueVsComb_kpt_split.png");
  //
  TCanvas c14a("c14a","",1);
  h_trueBestMatch_keta->DrawNormalized();
  h_comb_keta->DrawNormalized("same");
  h_comb1_keta->DrawNormalized("same");
  h_comb2_keta->DrawNormalized("same");
  h_comb3_keta->DrawNormalized("same");
  leg.Draw("same");  
  c14a.SaveAs("TrueVsComb_keta_split.png");
  //
  TCanvas c14b("c14b","",1);
  h_trueBestMatch_ele1eta->DrawNormalized();
  h_comb_ele1eta->DrawNormalized("same");
  h_comb1_ele1eta->DrawNormalized("same");
  h_comb2_ele1eta->DrawNormalized("same");
  h_comb3_ele1eta->DrawNormalized("same");
  leg.Draw("same");  
  c14b.SaveAs("TrueVsComb_ele1eta_split.png");
  //
  TCanvas c14c("c14c","",1);
  h_trueBestMatch_ele2eta->DrawNormalized();
  h_comb_ele2eta->DrawNormalized("same");
  h_comb1_ele2eta->DrawNormalized("same");
  h_comb2_ele2eta->DrawNormalized("same");
  h_comb3_ele2eta->DrawNormalized("same");
  leg.Draw("same");  
  c14c.SaveAs("TrueVsComb_ele2eta_split.png");
  //
  TCanvas c14d("c14d","",1);
  h_comb_ele1pt->DrawNormalized();
  h_comb1_ele1pt->DrawNormalized("same");
  h_comb2_ele1pt->DrawNormalized("same");
  h_comb3_ele1pt->DrawNormalized("same");
  h_trueBestMatch_ele1pt->DrawNormalized("same");
  leg.Draw("same");  
  c14d.SetLogy();
  c14d.SaveAs("TrueVsComb_ele1pt_split.png");
  //
  TCanvas c14e("c14e","",1);
  h_comb_ele2pt->DrawNormalized();
  h_comb1_ele2pt->DrawNormalized("same");
  h_comb2_ele2pt->DrawNormalized("same");
  h_comb3_ele2pt->DrawNormalized("same");
  h_trueBestMatch_ele2pt->DrawNormalized("same");
  leg.Draw("same");  
  c14e.SetLogy();
  c14e.SaveAs("TrueVsComb_ele2pt_split.png");
  //
  TCanvas c14f("c14f","",1);
  h_trueBestMatch_pfmva1->DrawNormalized();
  h_comb_pfmva1->DrawNormalized("same");
  h_comb1_pfmva1->DrawNormalized("same");
  h_comb2_pfmva1->DrawNormalized("same");
  h_comb3_pfmva1->DrawNormalized("same");
  leg2.Draw("same");  
  c14f.SaveAs("TrueVsComb_ele1pfmvaId_split.png");
  //
  TCanvas c14fbis("c14fbis","",1);
  h_trueBestMatch_pfmva1->DrawNormalized();
  h_comb_pfmva1_badele1->DrawNormalized("same");
  h_comb1_pfmva1_badele1->DrawNormalized("same");
  h_comb2_pfmva1_badele1->DrawNormalized("same");
  h_comb3_pfmva1_badele1->DrawNormalized("same");
  leg2.Draw("same");  
  c14fbis.SaveAs("TrueVsComb_ele1pfmvaId_ele1bad_split.png");
  //
  TCanvas c14g("c14g","",1);
  h_trueBestMatch_pfmva2->DrawNormalized();
  h_comb_pfmva2->DrawNormalized("same");
  h_comb1_pfmva2->DrawNormalized("same");
  h_comb2_pfmva2->DrawNormalized("same");
  h_comb3_pfmva2->DrawNormalized("same");
  leg2.Draw("same");  
  c14g.SaveAs("TrueVsComb_ele2pfmvaId_split.png");
  //
  TCanvas c14gbis("c14gbis","",1);
  h_trueBestMatch_pfmva2->DrawNormalized();
  h_comb_pfmva2_badele2->DrawNormalized("same");
  h_comb1_pfmva2_badele2->DrawNormalized("same");
  h_comb2_pfmva2_badele2->DrawNormalized("same");
  h_comb3_pfmva2_badele2->DrawNormalized("same");
  leg2.Draw("same");  
  c14gbis.SaveAs("TrueVsComb_ele2pfmvaId_ele2bad_split.png");
  //
  TCanvas c14h("c14h","",1);
  h_trueBestMatch_lptmva1->DrawNormalized();
  h_comb_lptmva1->DrawNormalized("same");
  h_comb1_lptmva1->DrawNormalized("same");
  h_comb2_lptmva1->DrawNormalized("same");
  h_comb3_lptmva1->DrawNormalized("same");
  leg2.Draw("same");  
  c14h.SaveAs("TrueVsComb_ele1lptmvaId_split.png");
  //
  TCanvas c14hbis("c14hbis","",1);
  h_trueBestMatch_lptmva1->DrawNormalized();
  h_comb_lptmva1_badele1->DrawNormalized("same");
  h_comb1_lptmva1_badele1->DrawNormalized("same");
  h_comb2_lptmva1_badele1->DrawNormalized("same");
  h_comb3_lptmva1_badele1->DrawNormalized("same");
  leg2.Draw("same");  
  c14hbis.SaveAs("TrueVsComb_ele1lptmvaId_ele1bad_split.png");
  //
  TCanvas c14i("c14i","",1);
  h_trueBestMatch_lptmva2->DrawNormalized();
  h_comb_lptmva2->DrawNormalized("same");
  h_comb1_lptmva2->DrawNormalized("same");
  h_comb2_lptmva2->DrawNormalized("same");
  h_comb3_lptmva2->DrawNormalized("same");
  leg2.Draw("same");  
  c14i.SaveAs("TrueVsComb_ele2lptmvaId_split.png");
  //
  TCanvas c14ibis("c14ibis","",1);
  h_trueBestMatch_lptmva2->DrawNormalized();
  h_comb_lptmva2_badele2->DrawNormalized("same");
  h_comb1_lptmva2_badele2->DrawNormalized("same");
  h_comb2_lptmva2_badele2->DrawNormalized("same");
  h_comb3_lptmva2_badele2->DrawNormalized("same");
  leg2.Draw("same");  
  c14ibis.SaveAs("TrueVsComb_ele2lptmvaId_ele2bad_split.png");
  //
  TCanvas c141("c141","",1);
  hz_comb_kpt->DrawNormalized();
  hz_comb1_kpt->DrawNormalized("same");
  hz_comb2_kpt->DrawNormalized("same");
  hz_comb3_kpt->DrawNormalized("same");
  hz_trueBestMatch_kpt->DrawNormalized("same");
  leg.Draw("same");  
  c141.SaveAs("TrueVsComb_kpt_zoom_split.png");
  //
  TCanvas c141d("c141d","",1);
  hz_comb_ele1pt->DrawNormalized();
  hz_comb1_ele1pt->DrawNormalized("same");
  hz_comb2_ele1pt->DrawNormalized("same");
  hz_comb3_ele1pt->DrawNormalized("same");
  hz_trueBestMatch_ele1pt->DrawNormalized("same");
  leg.Draw("same");  
  c141d.SaveAs("TrueVsComb_ele1pt_zoom_split.png");
  //
  TCanvas c141e("c141e","",1);
  hz_comb_ele2pt->DrawNormalized();
  hz_comb1_ele2pt->DrawNormalized("same");
  hz_comb2_ele2pt->DrawNormalized("same");
  hz_comb3_ele2pt->DrawNormalized("same");
  hz_trueBestMatch_ele2pt->DrawNormalized("same");
  leg.Draw("same");  
  c141e.SaveAs("TrueVsComb_ele2pt_zoom_split.png");
  //
  TCanvas c141bis("c141bis","",1);
  hz_comb_kpt_badk->DrawNormalized();
  hz_comb1_kpt_badk->DrawNormalized("same");
  hz_comb2_kpt_badk->DrawNormalized("same");
  hz_comb3_kpt_badk->DrawNormalized("same");
  hz_trueBestMatch_kpt->DrawNormalized("same");
  leg.Draw("same");  
  c141bis.SaveAs("TrueVsComb_kpt_badKonly_zoom_split.png");
  //
  TCanvas c141dbis("c141dbis","",1);
  hz_comb_ele1pt_badele1->DrawNormalized();
  hz_comb1_ele1pt_badele1->DrawNormalized("same");
  hz_comb2_ele1pt_badele1->DrawNormalized("same");
  hz_comb3_ele1pt_badele1->DrawNormalized("same");
  hz_trueBestMatch_ele1pt->DrawNormalized("same");
  leg.Draw("same");  
  c141dbis.SaveAs("TrueVsComb_ele1pt_badele1only_zoom_split.png");
  //
  TCanvas c141ebis("c141ebis","",1);
  hz_comb1_ele2pt_badele2->DrawNormalized();
  hz_comb_ele2pt_badele2->DrawNormalized("same");
  hz_comb2_ele2pt_badele2->DrawNormalized("same");
  hz_comb3_ele2pt_badele2->DrawNormalized("same");
  hz_trueBestMatch_ele2pt->DrawNormalized("same");
  leg.Draw("same");  
  c141ebis.SaveAs("TrueVsComb_ele2pt_badele2only_zoom_split.png");
  //
  TCanvas c141f("c141f","",1);
  hz_comb1_minpt->DrawNormalized("same");
  hz_comb1_minpt->DrawNormalized("same");
  hz_comb2_minpt->DrawNormalized("same");
  hz_comb3_minpt->DrawNormalized("same");
  hz_trueBestMatch_minpt->DrawNormalized("same");
  leg.Draw("same");  
  c141f.SaveAs("TrueVsComb_minpt_zoom_split.png");
  //
}
