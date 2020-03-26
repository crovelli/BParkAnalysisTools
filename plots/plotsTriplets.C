#define plotsTriplets_cxx
#include "plotsTriplets.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <iostream> 

using namespace std;

void plotsTriplets::Loop()
{
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;

  // Histos
  TH1F *h_numB     = new TH1F("h_numB", "h_numB", 25, -0.5, 24.5);
  TH1F *h_numCombB = new TH1F("h_numCombB", "h_numCombB", 25, -0.5, 24.5);
  TH1F *h_numTrueB = new TH1F("h_numTrueB", "h_numTrueB", 10, -0.5, 9.5);

  TH1F *h_trueB_dRminmaxEle = new TH1F("h_trueB_dRminmaxEle", "h_trueB_dRminmaxEle", 100, 0., 0.05);  
  TH1F *h_combB_dRminmaxEle = new TH1F("h_combB_dRminmaxEle", "h_combB_dRminmaxEle", 100, 0., 6.);  
  TH1F *h_trueB_dRmaxTrack  = new TH1F("h_trueB_dRmaxTrack",  "h_trueB_dRmaxTrack",  100, 0., 6.);  
  TH1F *h_combB_dRmaxTrack  = new TH1F("h_combB_dRmaxTrack",  "h_combB_dRmaxTrack",  100, 0., 6.);  
  TH1F *h_trueB_dRGenEle    = new TH1F("h_trueB_dRGenEle",    "h_trueB_dRGenEle",    100, 0., 6.);  

  TH1F *h_debug_svProbMatch      = new TH1F("h_debug_svProbMatch",      "h_debug_svProbMatch",      50, -0.01, 1.01);
  TH1F *h_debug_svProbUnMatch    = new TH1F("h_debug_svProbUnMatch",    "h_debug_svProbUnMatch",    50, -0.01, 1.01);
  TH1F *h_debug_pf_svProbMatch   = new TH1F("h_debug_pf_svProbMatch",   "h_debug_pf_svProbMatch",   50, -0.01, 1.01);
  TH1F *h_debug_pf_svProbUnMatch = new TH1F("h_debug_pf_svProbUnMatch", "h_debug_pf_svProbUnMatch", 50, -0.01, 1.01);

  TH1F *h_trueBestMatch_svProb   = new TH1F("h_trueBestMatch_svProb", "h_trueBestMatch_svProb", 100, -0.01, 1.01);
  TH1F *h_trueBestMatch_xysig    = new TH1F("h_trueBestMatch_xysig",  "h_trueBestMatch_xysig",  50,  0., 100.);
  TH1F *h_trueBestMatch_cos2d    = new TH1F("h_trueBestMatch_cos2d",  "h_trueBestMatch_cos2d",  100,  0.989, 1.001);
  TH1F *h_trueBestMatch_eleptsum = new TH1F("h_trueBestMatch_eleptsum", "h_trueBestMatch_eleptsum", 100, 0., 100.);
  TH1F *h_trueBestMatch_kpt      = new TH1F("h_trueBestMatch_kpt",      "h_trueBestMatch_kpt",      60, 0., 30.);
  TH1F *h_trueBestMatch_keta     = new TH1F("h_trueBestMatch_keta",    "h_trueBestMatch_keta",    50,-2.5,2.5);
  TH1F *h_trueBestMatch_ele1eta  = new TH1F("h_trueBestMatch_ele1eta", "h_trueBestMatch_ele1eta", 50,-2.5,2.5);
  TH1F *h_trueBestMatch_ele2eta  = new TH1F("h_trueBestMatch_ele2eta", "h_trueBestMatch_ele2eta", 50,-2.5,2.5);
  TH1F *h_trueBestMatch_ele2eta_lowpt_ptlt1   = new TH1F("h_trueBestMatch_ele2eta_lowpt_ptlt1",   "h_trueBestMatch_ele2eta_lowpt_ptlt1",   50,-2.5,2.5);
  TH1F *h_trueBestMatch_ele2eta_lowpt_ptlt2   = new TH1F("h_trueBestMatch_ele2eta_lowpt_ptlt2",   "h_trueBestMatch_ele2eta_lowpt_ptlt2",   50,-2.5,2.5);
  TH1F *h_trueBestMatch_ele2eta_lowpt_ptlt1d5 = new TH1F("h_trueBestMatch_ele2eta_lowpt_ptlt1d5", "h_trueBestMatch_ele2eta_lowpt_ptlt1d5", 50,-2.5,2.5);
  TH1F *h_trueBestMatch_ele2eta_lowpt_idlt0   = new TH1F("h_trueBestMatch_ele2eta_lowpt_idlt0",   "h_trueBestMatch_ele2eta_lowpt_idlt0",   50,-2.5,2.5);
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
  TH1F *h_trueBestMatch_maxDrRecoGen    = new TH1F("h_trueBestMatch_maxDrRecoGen",    "h_trueBestMatch_maxDrRecoGen",    100, 0., 6.);
  TH1F *h_trueBestMatch_minDrRecoGen    = new TH1F("h_trueBestMatch_minDrRecoGen",    "h_trueBestMatch_minDrRecoGen",    100, 0., 6.);
  TH1F *h_trueBestMatch_drRecoGenTrack  = new TH1F("h_trueBestMatch_drRecoGenTrack",  "h_trueBestMatch_drRecoGenTrack",  100, 0., 6.);
  TH2F *hz_trueBestMatch_ele2Vskpt = new TH2F("hz_trueBestMatch_ele2Vskpt", "hz_trueBestMatch_ele2Vskpt", 23, 0.7, 3., 25, 0.5, 3.);
 
  TH1F *h_trueBestMatch_Bmass     = new TH1F("h_trueBestMatch_Bmass",      "h_trueBestMatch_Bmass",      50, 4.5, 5.7);
  TH1F *h_trueBestMatch_BmassPFPF = new TH1F("h_trueBestMatch_BmassPFPF",  "h_trueBestMatch_BmassPFPF",  50, 4.5, 5.7);
  TH1F *h_trueBestMatch_BmassPFLP = new TH1F("h_trueBestMatch_BmassPFLP",  "h_trueBestMatch_BmassPFLP",  50, 4.5, 5.7);
  TH1F *h_trueBestMatch_BmassLPLP = new TH1F("h_trueBestMatch_BmassLPLP",  "h_trueBestMatch_BmassLPLP",  50, 4.5, 5.7);
  TH1F *h_trueBestMatch_BmassPFLP_ptgt2 = new TH1F("h_trueBestMatch_BmassPFLP_ptgt2",  "h_trueBestMatch_BmassPFLP_ptgt2",  50, 4.5, 5.7);
  TH1F *h_trueBestMatch_BmassPFLP_ptlt2 = new TH1F("h_trueBestMatch_BmassPFLP_ptlt2",  "h_trueBestMatch_BmassPFLP_ptlt2",  50, 4.5, 5.7);  
  TH1F *h_trueBestMatch_BmassLPLP_ptgt2 = new TH1F("h_trueBestMatch_BmassLPLP_ptgt2",  "h_trueBestMatch_BmassLPLP_ptgt2",  50, 4.5, 5.7);
  TH1F *h_trueBestMatch_BmassLPLP_ptlt2 = new TH1F("h_trueBestMatch_BmassLPLP_ptlt2",  "h_trueBestMatch_BmassLPLP_ptlt2",  50, 4.5, 5.7);
  TH1F *h_trueBestMatch_BmassPFLP_idgt2 = new TH1F("h_trueBestMatch_BmassPFLP_idgt2",  "h_trueBestMatch_BmassPFLP_idgt2",  50, 4.5, 5.7);
  TH1F *h_trueBestMatch_BmassPFLP_idlt2 = new TH1F("h_trueBestMatch_BmassPFLP_idlt2",  "h_trueBestMatch_BmassPFLP_idlt2",  50, 4.5, 5.7);  
  TH1F *h_trueBestMatch_BmassLPLP_idgt2 = new TH1F("h_trueBestMatch_BmassLPLP_idgt2",  "h_trueBestMatch_BmassLPLP_idgt2",  50, 4.5, 5.7);
  TH1F *h_trueBestMatch_BmassLPLP_idlt2 = new TH1F("h_trueBestMatch_BmassLPLP_idlt2",  "h_trueBestMatch_BmassLPLP_idlt2",  50, 4.5, 5.7);
  TH1F *h_trueBestMatch_Bmass_PFPForPFLP_ele2PF    = new TH1F("h_trueBestMatch_Bmass_PFPForPFLP_ele2PF",    "h_trueBestMatch_Bmass_PFPForPFLP_ele2PF",    50, 4.5, 5.7);
  TH1F *h_trueBestMatch_Bmass_PFLPorLPLP_ele2PtGt2 = new TH1F("h_trueBestMatch_Bmass_PFLPorLPLP_ele2PtGt2", "h_trueBestMatch_Bmass_PFLPorLPLP_ele2PtGt2", 50, 4.5, 5.7);
  TH1F *h_trueBestMatch_Bmass_PFLPorLPLP_ele2PtLt2 = new TH1F("h_trueBestMatch_Bmass_PFLPorLPLP_ele2PtLt2", "h_trueBestMatch_Bmass_PFLPorLPLP_ele2PtLt2", 50, 4.5, 5.7);

  TH1F *h_comb_svProb   = new TH1F("h_comb_svProb", "h_comb_svProb", 100, -0.01, 1.01);
  TH1F *h_comb_xysig    = new TH1F("h_comb_xysig",  "h_comb_xysig",  50,  0., 100.);
  TH1F *h_comb_cos2d    = new TH1F("h_comb_cos2d",  "h_comb_cos2d",  100,  0.989, 1.001);
  TH1F *h_comb_eleptsum = new TH1F("h_comb_eleptsum", "h_comb_eleptsum", 100, 0., 100.);
  TH1F *h_comb_kpt      = new TH1F("h_comb_kpt",      "h_comb_kpt",      60, 0., 30.);
  TH1F *h_comb_keta     = new TH1F("h_comb_keta",    "h_comb_keta",    50,-2.5,2.5);
  TH1F *h_comb_ele1eta  = new TH1F("h_comb_ele1eta", "h_comb_ele1eta", 50,-2.5,2.5);
  TH1F *h_comb_ele2eta  = new TH1F("h_comb_ele2eta", "h_comb_ele2eta", 50,-2.5,2.5);
  TH1F *h_comb_ele2eta_lowpt_ptlt1   = new TH1F("h_comb_ele2eta_lowpt_ptlt1",   "h_comb_ele2eta_lowpt_ptlt1",   50,-2.5,2.5);
  TH1F *h_comb_ele2eta_lowpt_ptlt1d5 = new TH1F("h_comb_ele2eta_lowpt_ptlt1d5", "h_comb_ele2eta_lowpt_ptlt1d5", 50,-2.5,2.5);
  TH1F *h_comb_ele2eta_lowpt_ptlt2   = new TH1F("h_comb_ele2eta_lowpt_ptlt2",   "h_comb_ele2eta_lowpt_ptlt2",   50,-2.5,2.5);
  TH1F *h_comb_ele2eta_lowpt_idlt0   = new TH1F("h_comb_ele2eta_lowpt_idlt0",   "h_comb_ele2eta_lowpt_idlt0",   50,-2.5,2.5);
  TH1F *h_comb_ele1pt   = new TH1F("h_comb_ele1pt",  "h_comb_ele1pt",  60, 0., 30.);
  TH1F *h_comb_ele2pt   = new TH1F("h_comb_ele2pt",  "h_comb_ele2pt",  60, 0., 30.);
  TH1F *h_comb_pfmva1   = new TH1F("h_comb_pfmva1",  "h_comb_pfmva1",  50,-10.,10.);
  TH1F *h_comb_pfmva2   = new TH1F("h_comb_pfmva2",  "h_comb_pfmva2",  50,-10.,10.);
  TH1F *h_comb_lptmva1  = new TH1F("h_comb_lptmva1", "h_comb_lptmva1", 50,-10.,10.);
  TH1F *h_comb_lptmva2  = new TH1F("h_comb_lptmva2", "h_comb_lptmva2", 50,-10.,10.);
  TH1F *h_comb_pfmva1_badele1   = new TH1F("h_comb_pfmva1_badele1",  "h_comb_pfmva1_badele1",  50,-10.,10.);
  TH1F *h_comb_pfmva2_badele2   = new TH1F("h_comb_pfmva2_badele2",  "h_comb_pfmva2_badele2",  50,-10.,10.);
  TH1F *h_comb_lptmva1_badele1  = new TH1F("h_comb_lptmva1_badele1", "h_comb_lptmva1_badele1", 50,-10.,10.);
  TH1F *h_comb_lptmva2_badele2  = new TH1F("h_comb_lptmva2_badele2", "h_comb_lptmva2_badele2", 50,-10.,10.);
  TH1F *hz_comb_kpt     = new TH1F("hz_comb_kpt",    "hz_comb_kpt",    50, 0., 5.);
  TH1F *hz_comb_ele1pt  = new TH1F("hz_comb_ele1pt", "hz_comb_ele1pt", 50, 0., 10.);
  TH1F *hz_comb_ele2pt  = new TH1F("hz_comb_ele2pt", "hz_comb_ele2pt", 50, 0., 5.);
  TH1F *hz_comb_kpt_badk       = new TH1F("hz_comb_kpt_badk",       "hz_comb_kpt_badk",        50, 0., 5.);
  TH1F *hz_comb_ele1pt_badele1 = new TH1F("hz_comb_ele1pt_badele1", "hz_comb_ele1pt_badele1",  50, 0., 10.);
  TH1F *hz_comb_ele2pt_badele2 = new TH1F("hz_comb_ele2pt_badele2", "hz_comb_ele2pt_badele2",  50, 0., 5.);
  TH1F *hz_comb_minpt   = new TH1F("hz_comb_minpt",  "hz_comb_minpt",  50, 0., 5.);
  TH1F *h_comb_maxDrRecoGen    = new TH1F("h_comb_maxDrRecoGen",    "h_comb_maxDrRecoGen",    100, 0., 6.);
  TH1F *h_comb_minDrRecoGen    = new TH1F("h_comb_minDrRecoGen",    "h_comb_minDrRecoGen",    100, 0., 6.);
  TH1F *h_comb_drRecoGenTrack  = new TH1F("h_comb_drRecoGenTrack",  "h_comb_drRecoGenTrack",  100, 0., 6.);
  TH2F *hz_comb_ele2Vskpt = new TH2F("hz_comb_ele2Vskpt", "hz_comb_ele2Vskpt", 23, 0.7,3., 25, 0.5,3.);

  TH1F *h_svProbM_notok_ele1pt  = new TH1F("h_svProbM_notok_ele1pt", "h_svProbM_notok_ele1pt",  60, 0.,30.);
  TH1F *h_svProbM_notok_ele2pt  = new TH1F("h_svProbM_notok_ele2pt", "h_svProbM_notok_ele2pt",  60, 0.,30.);
  TH1F *h_svProbM_notok_kpt     = new TH1F("h_svProbM_notok_kpt",    "h_svProbM_notok_kpt",     60, 0.,30.);
  TH1F *h_svProbM_notok_ele1eta = new TH1F("h_svProbM_notok_ele1eta","h_svProbM_notok_ele1eta", 50,-2.5,2.5);
  TH1F *h_svProbM_notok_ele2eta = new TH1F("h_svProbM_notok_ele2eta","h_svProbM_notok_ele2eta", 50,-2.5,2.5);
  TH1F *h_svProbM_notok_keta    = new TH1F("h_svProbM_notok_keta",   "h_svProbM_notok_keta",    50,-2.5,2.5);
  TH1F *h_svProbM_ok_ele1pt     = new TH1F("h_svProbM_ok_ele1pt",    "h_svProbM_ok_ele1pt",     60, 0.,30.);
  TH1F *h_svProbM_ok_ele2pt     = new TH1F("h_svProbM_ok_ele2pt",    "h_svProbM_ok_ele2pt",     60, 0.,30.);
  TH1F *h_svProbM_ok_kpt        = new TH1F("h_svProbM_ok_kpt",       "h_svProbM_ok_kpt",        60, 0.,30.);
  TH1F *h_svProbM_ok_ele1eta    = new TH1F("h_svProbM_ok_ele1eta",   "h_svProbM_ok_ele1eta",    50,-2.5,2.5);
  TH1F *h_svProbM_ok_ele2eta    = new TH1F("h_svProbM_ok_ele2eta",   "h_svProbM_ok_ele2eta",    50,-2.5,2.5);
  TH1F *h_svProbM_ok_keta       = new TH1F("h_svProbM_ok_keta",      "h_svProbM_ok_keta",       50,-2.5,2.5);
  TH1F *hz_svProbM_notok_ele1pt = new TH1F("hz_svProbM_notok_ele1pt","hz_svProbM_notok_ele1pt", 50, 0.,10.);
  TH1F *hz_svProbM_notok_ele2pt = new TH1F("hz_svProbM_notok_ele2pt","hz_svProbM_notok_ele2pt", 50, 0.,5.);
  TH1F *hz_svProbM_notok_kpt    = new TH1F("hz_svProbM_notok_kpt",   "hz_svProbM_notok_kpt",    50, 0.,5.);
  TH1F *hz_svProbM_ok_ele1pt    = new TH1F("hz_svProbM_ok_ele1pt",   "hz_svProbM_ok_ele1pt",    50, 0.,10.);
  TH1F *hz_svProbM_ok_ele2pt    = new TH1F("hz_svProbM_ok_ele2pt",   "hz_svProbM_ok_ele2pt",    50, 0.,5.);
  TH1F *hz_svProbM_ok_kpt       = new TH1F("hz_svProbM_ok_kpt",      "hz_svProbM_ok_kpt",       50, 0.,5.);

  TH1F *h_xySigM_notok_ele1pt   = new TH1F("h_xySigM_notok_ele1pt",  "h_xySigM_notok_ele1pt",  60, 0.,30.);
  TH1F *h_xySigM_notok_ele2pt   = new TH1F("h_xySigM_notok_ele2pt",  "h_xySigM_notok_ele2pt",  60, 0.,30.);
  TH1F *h_xySigM_notok_kpt      = new TH1F("h_xySigM_notok_kpt",     "h_xySigM_notok_kpt",     60, 0.,30.);
  TH1F *h_xySigM_notok_ele1eta  = new TH1F("h_xySigM_notok_ele1eta", "h_xySigM_notok_ele1eta", 50,-2.5,2.5);
  TH1F *h_xySigM_notok_ele2eta  = new TH1F("h_xySigM_notok_ele2eta", "h_xySigM_notok_ele2eta", 50,-2.5,2.5);
  TH1F *h_xySigM_notok_keta     = new TH1F("h_xySigM_notok_keta",    "h_xySigM_notok_keta",    50,-2.5,2.5);
  TH1F *h_xySigM_notok_pfmva1   = new TH1F("h_xySigM_notok_pfmva1",  "h_xySigM_notok_pfmva1",  50,-10.,10.);
  TH1F *h_xySigM_notok_pfmva2   = new TH1F("h_xySigM_notok_pfmva2",  "h_xySigM_notok_pfmva2",  50,-10.,10.);
  TH1F *h_xySigM_notok_lptmva1  = new TH1F("h_xySigM_notok_lptmva1", "h_xySigM_notok_lptmva1", 50,-10.,10.);
  TH1F *h_xySigM_notok_lptmva2  = new TH1F("h_xySigM_notok_lptmva2", "h_xySigM_notok_lptmva2", 50,-10.,10.);
  TH1F *h_xySigM_notok_costhetaSK   = new TH1F("h_xySigM_notok_costhetaSK",   "h_xySigM_notok_costhetaSK",   11,0.,1.1);
  TH1F *h_xySigM_notok_costhetaSKCS = new TH1F("h_xySigM_notok_costhetaSKCS", "h_xySigM_notok_costhetaSKCS", 11,0.,1.1);
  TH1F *h_xySigM_notok_costhetaL    = new TH1F("h_xySigM_notok_costhetaL",    "h_xySigM_notok_costhetaL",    11,0.,1.1);
  TH1F *h_xySigM_ok_ele1pt      = new TH1F("h_xySigM_ok_ele1pt",     "h_xySigM_ok_ele1pt",     60, 0.,30.);
  TH1F *h_xySigM_ok_ele2pt      = new TH1F("h_xySigM_ok_ele2pt",     "h_xySigM_ok_ele2pt",     60, 0.,30.);
  TH1F *h_xySigM_ok_kpt         = new TH1F("h_xySigM_ok_kpt",        "h_xySigM_ok_kpt",        60, 0.,30.);
  TH1F *h_xySigM_ok_ele1eta     = new TH1F("h_xySigM_ok_ele1eta",    "h_xySigM_ok_ele1eta",    50,-2.5,2.5);
  TH1F *h_xySigM_ok_ele2eta     = new TH1F("h_xySigM_ok_ele2eta",    "h_xySigM_ok_ele2eta",    50,-2.5,2.5);
  TH1F *h_xySigM_ok_keta        = new TH1F("h_xySigM_ok_keta",       "h_xySigM_ok_keta",       50,-2.5,2.5);
  TH1F *h_xySigM_ok_pfmva1      = new TH1F("h_xySigM_ok_pfmva1",     "h_xySigM_ok_pfmva1",     50,-10.,10.);
  TH1F *h_xySigM_ok_pfmva2      = new TH1F("h_xySigM_ok_pfmva2",     "h_xySigM_ok_pfmva2",     50,-10.,10.);
  TH1F *h_xySigM_ok_lptmva1     = new TH1F("h_xySigM_ok_lptmva1",    "h_xySigM_ok_lptmva1",    50,-10.,10.);
  TH1F *h_xySigM_ok_lptmva2     = new TH1F("h_xySigM_ok_lptmva2",    "h_xySigM_ok_lptmva2",    50,-10.,10.);
  TH1F *h_xySigM_ok_costhetaSK_gen = new TH1F("h_xySigM_ok_costhetaSK_gen", "h_xySigM_ok_costhetaSK_gen", 11,0.,1.1);
  TH1F *h_xySigM_ok_costhetaSK     = new TH1F("h_xySigM_ok_costhetaSK",     "h_xySigM_ok_costhetaSK",     11,0.,1.1);
  TH1F *h_xySigM_ok_costhetaSKCS   = new TH1F("h_xySigM_ok_costhetaSKCS",   "h_xySigM_ok_costhetaSKCS",   11,0.,1.1);
  TH1F *h_xySigM_ok_costhetaL      = new TH1F("h_xySigM_ok_costhetaL",      "h_xySigM_ok_costhetaL",      11,0.,1.1);
  TH1F *hz_xySigM_notok_ele1pt  = new TH1F("hz_xySigM_notok_ele1pt", "hz_xySigM_notok_ele1pt", 50, 0.,10.);
  TH1F *hz_xySigM_notok_ele2pt  = new TH1F("hz_xySigM_notok_ele2pt", "hz_xySigM_notok_ele2pt", 50, 0.,5.);
  TH1F *hz_xySigM_notok_kpt     = new TH1F("hz_xySigM_notok_kpt",    "hz_xySigM_notok_kpt",    50, 0.,5.);
  TH1F *hz_xySigM_notok_minpt   = new TH1F("hz_xySigM_notok_minpt",  "hz_xySigM_notok_minpt",  50, 0.,5.);
  TH1F *hz_xySigM_ok_ele1pt     = new TH1F("hz_xySigM_ok_ele1pt",    "hz_xySigM_ok_ele1pt",    50, 0.,10.);
  TH1F *hz_xySigM_ok_ele2pt     = new TH1F("hz_xySigM_ok_ele2pt",    "hz_xySigM_ok_ele2pt",    50, 0.,5.);
  TH1F *hz_xySigM_ok_kpt        = new TH1F("hz_xySigM_ok_kpt",       "hz_xySigM_ok_kpt",       50, 0.,5.);
  TH1F *hz_xySigM_ok_minpt      = new TH1F("hz_xySigM_ok_minpt",     "hz_xySigM_ok_minpt",     50, 0.,5.);
  TH2F *hz_xySigM_notok_ele2Vskpt = new TH2F("hz_xySigM_notok_ele2Vskpt", "hz_xySigM_noyok_ele2Vskpt", 50, 0.,5., 50, 0.,5.);
  TH2F *hz_xySigM_ok_ele2Vskpt    = new TH2F("hz_xySigM_ok_ele2Vskpt",    "hz_xySigM_ok_ele2Vskpt",    50, 0.,5., 50, 0.,5.);

  TH1F *h_cos2DM_notok_ele1pt   = new TH1F("h_cos2DM_notok_ele1pt",  "h_cos2DM_notok_ele1pt",  60, 0.,30.);
  TH1F *h_cos2DM_notok_ele2pt   = new TH1F("h_cos2DM_notok_ele2pt",  "h_cos2DM_notok_ele2pt",  60, 0.,30.);
  TH1F *h_cos2DM_notok_kpt      = new TH1F("h_cos2DM_notok_kpt",     "h_cos2DM_notok_kpt",     60, 0.,30.);
  TH1F *h_cos2DM_notok_ele1eta  = new TH1F("h_cos2DM_notok_ele1eta", "h_cos2DM_notok_ele1eta", 50,-2.5,2.5);
  TH1F *h_cos2DM_notok_ele2eta  = new TH1F("h_cos2DM_notok_ele2eta", "h_cos2DM_notok_ele2eta", 50,-2.5,2.5);
  TH1F *h_cos2DM_notok_keta     = new TH1F("h_cos2DM_notok_keta",    "h_cos2DM_notok_keta",    50,-2.5,2.5);
  TH1F *h_cos2DM_ok_ele1pt      = new TH1F("h_cos2DM_ok_ele1pt",     "h_cos2DM_ok_ele1pt",     60, 0.,30.);
  TH1F *h_cos2DM_ok_ele2pt      = new TH1F("h_cos2DM_ok_ele2pt",     "h_cos2DM_ok_ele2pt",     60, 0.,30.);
  TH1F *h_cos2DM_ok_kpt         = new TH1F("h_cos2DM_ok_kpt",        "h_cos2DM_ok_kpt",        60, 0.,30.);
  TH1F *h_cos2DM_ok_ele1eta     = new TH1F("h_cos2DM_ok_ele1eta",    "h_cos2DM_ok_ele1eta",    50,-2.5,2.5);
  TH1F *h_cos2DM_ok_ele2eta     = new TH1F("h_cos2DM_ok_ele2eta",    "h_cos2DM_ok_ele2eta",    50,-2.5,2.5);
  TH1F *h_cos2DM_ok_keta        = new TH1F("h_cos2DM_ok_keta",       "h_cos2DM_ok_keta",       50,-2.5,2.5);
  TH1F *hz_cos2DM_notok_ele1pt  = new TH1F("hz_cos2DM_notok_ele1pt", "hz_cos2DM_notok_ele1pt", 50, 0.,10.);
  TH1F *hz_cos2DM_notok_ele2pt  = new TH1F("hz_cos2DM_notok_ele2pt", "hz_cos2DM_notok_ele2pt", 50, 0.,5.);
  TH1F *hz_cos2DM_notok_kpt     = new TH1F("hz_cos2DM_notok_kpt",    "hz_cos2DM_notok_kpt",    50, 0.,5.);
  TH1F *hz_cos2DM_ok_ele1pt     = new TH1F("hz_cos2DM_ok_ele1pt",    "hz_cos2DM_ok_ele1pt",    50, 0.,10.);
  TH1F *hz_cos2DM_ok_ele2pt     = new TH1F("hz_cos2DM_ok_ele2pt",    "hz_cos2DM_ok_ele2pt",    50, 0.,5.);
  TH1F *hz_cos2DM_ok_kpt        = new TH1F("hz_cos2DM_ok_kpt",       "hz_cos2DM_ok_kpt",       50, 0.,5.);

  // Counters
  float tot_ok_all = 0.;
  float tot_ok_presel = 0.;
  float tot_notok_all = 0.;
  float tot_notok_presel = 0.;

  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    // debug
    int debugsize   = debug_svprob_match->size();
    int debugpfsize = debug_pf_svprob_match->size();
    for (int ii=0; ii<debugsize; ii++) {
      if (debug_svprob_match->at(ii)==1)    h_debug_svProbMatch      -> Fill(debug_svprob->at(ii));
      if (debug_svprob_match->at(ii)==0)    h_debug_svProbUnMatch    -> Fill(debug_svprob->at(ii));
    }
    for (int ii=0; ii<debugpfsize; ii++) {
      if (debug_pf_svprob_match->at(ii)==1) h_debug_pf_svProbMatch   -> Fill(debug_pf_svprob->at(ii));
      if (debug_pf_svprob_match->at(ii)==0) h_debug_pf_svProbUnMatch -> Fill(debug_pf_svprob->at(ii));
    }

    // multiplicity
    h_numB     -> Fill(goodBSize);
    h_numCombB -> Fill(goodCombBSize);
    h_numTrueB -> Fill(goodTrueBSize);

    // which true
    int trueBminmaxsize = goodTrueB_maxMinDREle->size();
    for (int ii=0; ii<trueBminmaxsize; ii++) { 
      h_trueB_dRminmaxEle->Fill(goodTrueB_maxMinDREle->at(ii));
      h_trueB_dRmaxTrack->Fill(goodTrueB_maxDRTrack->at(ii));      
    }
    int combBminmaxsize = goodCombB_maxMinDREle->size();
    for (int ii=0; ii<combBminmaxsize; ii++) { 
      h_combB_dRminmaxEle->Fill(goodCombB_maxMinDREle->at(ii));
      h_combB_dRmaxTrack->Fill(goodCombB_maxDRTrack->at(ii));      
    }
    int trueBdrgensize = goodTrueB_dRgen->size();
    for (int ii=0; ii<trueBdrgensize; ii++) h_trueB_dRGenEle->Fill(goodTrueB_dRgen->at(ii));

    // true vs comb: best match - only when there is at least one good
    if (bestMatch_SvProb>-999) {

      // to study possible categories
      bool is1pf = true;
      bool is2pf = true;
      if (bestMatch_Ele1pfmva<20) is1pf = true;
      else is1pf = false;
      if (bestMatch_Ele2pfmva<20) is2pf = true;
      else is2pf = false;
      // all
      h_trueBestMatch_Bmass -> Fill(bestMatch_Bmass);
      // standard cat
      if (is1pf && is2pf)                           h_trueBestMatch_BmassPFPF -> Fill(bestMatch_Bmass);
      if ( (is1pf && !is2pf) || (is2pf && !is1pf) ) h_trueBestMatch_BmassPFLP -> Fill(bestMatch_Bmass);
      if (!is1pf && !is2pf)                         h_trueBestMatch_BmassLPLP -> Fill(bestMatch_Bmass);
      // various checks
      if ( (is1pf && !is2pf && bestMatch_Ele2Pt<2) || (is2pf && !is1pf && bestMatch_Ele1Pt<2) ) h_trueBestMatch_BmassPFLP_ptlt2 -> Fill(bestMatch_Bmass);
      if ( (is1pf && !is2pf && bestMatch_Ele2Pt>2) || (is2pf && !is1pf && bestMatch_Ele1Pt>2) ) h_trueBestMatch_BmassPFLP_ptgt2 -> Fill(bestMatch_Bmass);
      if ( !is1pf && !is2pf && bestMatch_Ele1Pt<2 && bestMatch_Ele2Pt<2) h_trueBestMatch_BmassLPLP_ptlt2 -> Fill(bestMatch_Bmass);
      if ( !is1pf && !is2pf && bestMatch_Ele1Pt>2 && bestMatch_Ele2Pt>2) h_trueBestMatch_BmassLPLP_ptgt2 -> Fill(bestMatch_Bmass);
      //
      if ( (is1pf && bestMatch_Ele2lptmva<2) || (is2pf && bestMatch_Ele1lptmva<2) ) h_trueBestMatch_BmassPFLP_idlt2 -> Fill(bestMatch_Bmass);
      if ( (is1pf && !is2pf && bestMatch_Ele2lptmva>2) || (is2pf && !is1pf && bestMatch_Ele1lptmva>2) ) h_trueBestMatch_BmassPFLP_idgt2 -> Fill(bestMatch_Bmass);
      if (!is1pf && !is2pf && bestMatch_Ele2lptmva<2) h_trueBestMatch_BmassLPLP_idlt2 -> Fill(bestMatch_Bmass);
      if (!is1pf && !is2pf && bestMatch_Ele2lptmva>2) h_trueBestMatch_BmassLPLP_idgt2 -> Fill(bestMatch_Bmass);
      // PFPF o PFLP con leading lowPT; PFLP/LPLP con subleading LP e pT>2; PFLP/LPLP con subleading LP e pT<2  
      if (is2pf) h_trueBestMatch_Bmass_PFPForPFLP_ele2PF -> Fill(bestMatch_Bmass);                            // PFPF o PFLP con leading lowPT
      if (!is2pf && bestMatch_Ele2Pt>2) h_trueBestMatch_Bmass_PFLPorLPLP_ele2PtGt2 -> Fill(bestMatch_Bmass);  // PFLP/LPLP con subleading LP e pT>2     
      if (!is2pf && bestMatch_Ele2Pt<2) h_trueBestMatch_Bmass_PFLPorLPLP_ele2PtLt2 -> Fill(bestMatch_Bmass);  // PFLP/LPLP con subleading LP e pT<2 
      //
      //
      h_trueBestMatch_svProb   -> Fill(bestMatch_SvProb);
      h_trueBestMatch_xysig    -> Fill(bestMatch_XYSig);
      h_trueBestMatch_cos2d    -> Fill(bestMatch_Cos2D);
      h_trueBestMatch_eleptsum -> Fill(bestMatch_PtSum);
      h_trueBestMatch_kpt      -> Fill(bestMatch_KPt);
      h_trueBestMatch_keta     -> Fill(bestMatch_KEta);
      h_trueBestMatch_ele1eta  -> Fill(bestMatch_Ele1Eta);
      h_trueBestMatch_ele2eta  -> Fill(bestMatch_Ele2Eta);
      if (bestMatch_Ele2lptmva<20 && bestMatch_Ele2Pt<1 )   h_trueBestMatch_ele2eta_lowpt_ptlt1   -> Fill(bestMatch_Ele2Eta);
      if (bestMatch_Ele2lptmva<20 && bestMatch_Ele2Pt<1.5 ) h_trueBestMatch_ele2eta_lowpt_ptlt1d5 -> Fill(bestMatch_Ele2Eta);
      if (bestMatch_Ele2lptmva<20 && bestMatch_Ele2Pt<2 )   h_trueBestMatch_ele2eta_lowpt_ptlt2   -> Fill(bestMatch_Ele2Eta);
      if (bestMatch_Ele2lptmva<0)                           h_trueBestMatch_ele2eta_lowpt_idlt0   -> Fill(bestMatch_Ele2Eta); 
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
      hz_trueBestMatch_ele2Vskpt     -> Fill(bestMatch_KPt, bestMatch_Ele2Pt); 
      h_trueBestMatch_maxDrRecoGen   -> Fill(bestMatch_maxDrRecoGen);
      h_trueBestMatch_minDrRecoGen   -> Fill(bestMatch_minDrRecoGen);
      h_trueBestMatch_drRecoGenTrack -> Fill(bestMatch_drRecoGenK);
    }

    // true vs comb: comb
    int goodCombB_size = goodCombB_svProb->size();
    for (int ii=0; ii<goodCombB_size; ii++) { 
      h_comb_svProb   -> Fill(goodCombB_svProb->at(ii));
      h_comb_xysig    -> Fill(goodCombB_xySig->at(ii));
      h_comb_cos2d    -> Fill(goodCombB_cos2D->at(ii));
      h_comb_eleptsum -> Fill(goodCombB_ptsum->at(ii));
      h_comb_kpt      -> Fill(goodCombB_kpt->at(ii));
      h_comb_keta     -> Fill(goodCombB_keta->at(ii));
      h_comb_ele1eta  -> Fill(goodCombB_ele1eta->at(ii));
      h_comb_ele2eta  -> Fill(goodCombB_ele2eta->at(ii));
      if (goodCombB_ele2lptmva->at(ii)<20 && goodCombB_ele2pt->at(ii)<1 )   h_comb_ele2eta_lowpt_ptlt1   -> Fill(goodCombB_ele2eta->at(ii));
      if (goodCombB_ele2lptmva->at(ii)<20 && goodCombB_ele2pt->at(ii)<1.5 ) h_comb_ele2eta_lowpt_ptlt1d5 -> Fill(goodCombB_ele2eta->at(ii));
      if (goodCombB_ele2lptmva->at(ii)<20 && goodCombB_ele2pt->at(ii)<2 )   h_comb_ele2eta_lowpt_ptlt2   -> Fill(goodCombB_ele2eta->at(ii));
      if (goodCombB_ele2lptmva->at(ii)<0)                                   h_comb_ele2eta_lowpt_idlt0   -> Fill(goodCombB_ele2eta->at(ii)); 
      h_comb_ele1pt   -> Fill(goodCombB_ele1pt->at(ii));
      h_comb_ele2pt   -> Fill(goodCombB_ele2pt->at(ii));
      h_comb_pfmva1   -> Fill(goodCombB_ele1pfmva->at(ii));
      h_comb_pfmva2   -> Fill(goodCombB_ele2pfmva->at(ii));
      h_comb_lptmva1  -> Fill(goodCombB_ele1lptmva->at(ii));
      h_comb_lptmva2  -> Fill(goodCombB_ele2lptmva->at(ii));
      if (goodCombB_causeEle1->at(ii)==1) h_comb_pfmva1_badele1  -> Fill(goodCombB_ele1pfmva->at(ii));
      if (goodCombB_causeEle2->at(ii)==1) h_comb_pfmva2_badele2  -> Fill(goodCombB_ele2pfmva->at(ii));
      if (goodCombB_causeEle1->at(ii)==1) h_comb_lptmva1_badele1 -> Fill(goodCombB_ele1lptmva->at(ii));
      if (goodCombB_causeEle2->at(ii)==1) h_comb_lptmva2_badele2 -> Fill(goodCombB_ele2lptmva->at(ii));
      hz_comb_kpt     -> Fill(goodCombB_kpt->at(ii));
      hz_comb_ele1pt  -> Fill(goodCombB_ele1pt->at(ii));
      hz_comb_ele2pt  -> Fill(goodCombB_ele2pt->at(ii));
      if (goodCombB_causeK->at(ii)==1)    hz_comb_kpt_badk       -> Fill(goodCombB_kpt->at(ii));
      if (goodCombB_causeEle1->at(ii)==1) hz_comb_ele1pt_badele1 -> Fill(goodCombB_ele1pt->at(ii));
      if (goodCombB_causeEle2->at(ii)==1) hz_comb_ele2pt_badele2 -> Fill(goodCombB_ele2pt->at(ii));
      hz_comb_minpt   -> Fill(goodCombB_minpt->at(ii));
      h_comb_maxDrRecoGen   -> Fill(goodCombB_maxDrRecoGen->at(ii));
      h_comb_minDrRecoGen   -> Fill(goodCombB_minDrRecoGen->at(ii));
      h_comb_drRecoGenTrack -> Fill(goodCombB_drRecoGenK->at(ii));
      hz_comb_ele2Vskpt     -> Fill(goodCombB_kpt->at(ii),goodCombB_ele2pt->at(ii));
    }

    // offline criteria: bad vs good
    if (bestSvProbMatch_notok_ele1pt>-990)  h_svProbM_notok_ele1pt  -> Fill(bestSvProbMatch_notok_ele1pt);
    if (bestSvProbMatch_notok_ele2pt>-990)  h_svProbM_notok_ele2pt  -> Fill(bestSvProbMatch_notok_ele2pt);
    if (bestSvProbMatch_notok_kpt>-990)     h_svProbM_notok_kpt     -> Fill(bestSvProbMatch_notok_kpt);
    if (bestSvProbMatch_notok_ele1eta>-990) h_svProbM_notok_ele1eta -> Fill(bestSvProbMatch_notok_ele1eta);
    if (bestSvProbMatch_notok_ele2eta>-990) h_svProbM_notok_ele2eta -> Fill(bestSvProbMatch_notok_ele2eta);
    if (bestSvProbMatch_notok_keta>-990)    h_svProbM_notok_keta    -> Fill(bestSvProbMatch_notok_keta);
    if (bestSvProbMatch_ok_ele1pt>-990)     h_svProbM_ok_ele1pt     -> Fill(bestSvProbMatch_ok_ele1pt);
    if (bestSvProbMatch_ok_ele2pt>-990)     h_svProbM_ok_ele2pt     -> Fill(bestSvProbMatch_ok_ele2pt);
    if (bestSvProbMatch_ok_kpt>-990)        h_svProbM_ok_kpt        -> Fill(bestSvProbMatch_ok_kpt);
    if (bestSvProbMatch_ok_ele1eta>-990)    h_svProbM_ok_ele1eta    -> Fill(bestSvProbMatch_ok_ele1eta);
    if (bestSvProbMatch_ok_ele2eta>-990)    h_svProbM_ok_ele2eta    -> Fill(bestSvProbMatch_ok_ele2eta);
    if (bestSvProbMatch_ok_keta>-990)       h_svProbM_ok_keta       -> Fill(bestSvProbMatch_ok_keta);
    if (bestSvProbMatch_notok_ele1pt>-990)  hz_svProbM_notok_ele1pt -> Fill(bestSvProbMatch_notok_ele1pt);
    if (bestSvProbMatch_notok_ele2pt>-990)  hz_svProbM_notok_ele2pt -> Fill(bestSvProbMatch_notok_ele2pt);
    if (bestSvProbMatch_notok_kpt>-990)     hz_svProbM_notok_kpt    -> Fill(bestSvProbMatch_notok_kpt);
    if (bestSvProbMatch_ok_ele1pt>-990)     hz_svProbM_ok_ele1pt    -> Fill(bestSvProbMatch_ok_ele1pt);
    if (bestSvProbMatch_ok_ele2pt>-990)     hz_svProbM_ok_ele2pt    -> Fill(bestSvProbMatch_ok_ele2pt);
    if (bestSvProbMatch_ok_kpt>-990)        hz_svProbM_ok_kpt       -> Fill(bestSvProbMatch_ok_kpt);

    // offline criteria: bad vs good
    if (bestXYsigMatch_notok_ele1pt>-990)  h_xySigM_notok_ele1pt  -> Fill(bestXYsigMatch_notok_ele1pt);
    if (bestXYsigMatch_notok_ele2pt>-990)  h_xySigM_notok_ele2pt  -> Fill(bestXYsigMatch_notok_ele2pt);
    if (bestXYsigMatch_notok_kpt>-990)     h_xySigM_notok_kpt     -> Fill(bestXYsigMatch_notok_kpt);
    if (bestXYsigMatch_notok_ele1eta>-990) h_xySigM_notok_ele1eta -> Fill(bestXYsigMatch_notok_ele1eta);
    if (bestXYsigMatch_notok_ele2eta>-990) h_xySigM_notok_ele2eta -> Fill(bestXYsigMatch_notok_ele2eta);
    if (bestXYsigMatch_notok_keta>-990)    h_xySigM_notok_keta    -> Fill(bestXYsigMatch_notok_keta);
    if (bestXYsigMatch_notok_pfmva1>-990  && bestXYsigMatch_notok_pfmva1<20)  h_xySigM_notok_pfmva1  -> Fill(bestXYsigMatch_notok_pfmva1);
    if (bestXYsigMatch_notok_pfmva2>-990  && bestXYsigMatch_notok_pfmva2<20)  h_xySigM_notok_pfmva2  -> Fill(bestXYsigMatch_notok_pfmva2);
    if (bestXYsigMatch_notok_lptmva1>-990 && bestXYsigMatch_notok_lptmva1<20) h_xySigM_notok_lptmva1 -> Fill(bestXYsigMatch_notok_lptmva1);
    if (bestXYsigMatch_notok_lptmva2>-990 && bestXYsigMatch_notok_lptmva2<20) h_xySigM_notok_lptmva2 -> Fill(bestXYsigMatch_notok_lptmva2);
    if (bestXYsigMatch_notok_costhetaSK>-990)   h_xySigM_notok_costhetaSK   -> Fill(fabs(bestXYsigMatch_notok_costhetaSK));
    if (bestXYsigMatch_notok_costhetaSKCS>-990) h_xySigM_notok_costhetaSKCS -> Fill(fabs(bestXYsigMatch_notok_costhetaSKCS));
    if (bestXYsigMatch_notok_costhetaL>-990)    h_xySigM_notok_costhetaL    -> Fill(fabs(bestXYsigMatch_notok_costhetaL));

    if (bestXYsigMatch_ok_ele1pt>-990)     h_xySigM_ok_ele1pt     -> Fill(bestXYsigMatch_ok_ele1pt);
    if (bestXYsigMatch_ok_ele2pt>-990)     h_xySigM_ok_ele2pt     -> Fill(bestXYsigMatch_ok_ele2pt);
    if (bestXYsigMatch_ok_kpt>-990)        h_xySigM_ok_kpt        -> Fill(bestXYsigMatch_ok_kpt);
    if (bestXYsigMatch_ok_ele1eta>-990)    h_xySigM_ok_ele1eta    -> Fill(bestXYsigMatch_ok_ele1eta);
    if (bestXYsigMatch_ok_ele2eta>-990)    h_xySigM_ok_ele2eta    -> Fill(bestXYsigMatch_ok_ele2eta);
    if (bestXYsigMatch_ok_keta>-990)       h_xySigM_ok_keta       -> Fill(bestXYsigMatch_ok_keta);
    if (bestXYsigMatch_ok_pfmva1>-990  && bestXYsigMatch_ok_pfmva1<20)  h_xySigM_ok_pfmva1  -> Fill(bestXYsigMatch_ok_pfmva1);
    if (bestXYsigMatch_ok_pfmva2>-990  && bestXYsigMatch_ok_pfmva2<20)  h_xySigM_ok_pfmva2  -> Fill(bestXYsigMatch_ok_pfmva2);
    if (bestXYsigMatch_ok_lptmva1>-990 && bestXYsigMatch_ok_lptmva1<20) h_xySigM_ok_lptmva1 -> Fill(bestXYsigMatch_ok_lptmva1);
    if (bestXYsigMatch_ok_lptmva2>-990 && bestXYsigMatch_ok_lptmva2<20) h_xySigM_ok_lptmva2 -> Fill(bestXYsigMatch_ok_lptmva2);
    if (bestXYsigMatch_ok_costhetaSK_gen>-990) h_xySigM_ok_costhetaSK_gen -> Fill(fabs(bestXYsigMatch_ok_costhetaSK_gen));
    if (bestXYsigMatch_ok_costhetaSK>-990)     h_xySigM_ok_costhetaSK     -> Fill(fabs(bestXYsigMatch_ok_costhetaSK));
    if (bestXYsigMatch_ok_costhetaSKCS>-990)   h_xySigM_ok_costhetaSKCS   -> Fill(fabs(bestXYsigMatch_ok_costhetaSKCS));
    if (bestXYsigMatch_ok_costhetaL>-990)      h_xySigM_ok_costhetaL      -> Fill(fabs(bestXYsigMatch_ok_costhetaL));

    if (bestXYsigMatch_notok_ele1pt>-990)  hz_xySigM_notok_ele1pt -> Fill(bestXYsigMatch_notok_ele1pt);
    if (bestXYsigMatch_notok_ele2pt>-990)  hz_xySigM_notok_ele2pt -> Fill(bestXYsigMatch_notok_ele2pt);
    if (bestXYsigMatch_notok_kpt>-990)     hz_xySigM_notok_kpt    -> Fill(bestXYsigMatch_notok_kpt);
    if (bestXYsigMatch_notok_ele2pt>-990 && bestXYsigMatch_notok_kpt>-990) { 
      hz_xySigM_notok_ele2Vskpt -> Fill(bestXYsigMatch_notok_kpt, bestXYsigMatch_notok_ele2pt);
      if (bestXYsigMatch_notok_ele2pt<bestXYsigMatch_notok_kpt)  hz_xySigM_notok_minpt -> Fill(bestXYsigMatch_notok_ele2pt);
      if (bestXYsigMatch_notok_ele2pt>=bestXYsigMatch_notok_kpt) hz_xySigM_notok_minpt -> Fill(bestXYsigMatch_notok_kpt);
    }
    if (bestXYsigMatch_ok_ele1pt>-990)     hz_xySigM_ok_ele1pt    -> Fill(bestXYsigMatch_ok_ele1pt);
    if (bestXYsigMatch_ok_ele2pt>-990)     hz_xySigM_ok_ele2pt    -> Fill(bestXYsigMatch_ok_ele2pt);
    if (bestXYsigMatch_ok_kpt>-990)        hz_xySigM_ok_kpt       -> Fill(bestXYsigMatch_ok_kpt);
    if (bestXYsigMatch_ok_ele2pt>-990 && bestXYsigMatch_ok_kpt>-990) { 
      hz_xySigM_ok_ele2Vskpt -> Fill(bestXYsigMatch_ok_kpt, bestXYsigMatch_ok_ele2pt);
      if (bestXYsigMatch_ok_ele2pt<bestXYsigMatch_ok_kpt)  hz_xySigM_ok_minpt -> Fill(bestXYsigMatch_ok_ele2pt);
      if (bestXYsigMatch_ok_ele2pt>=bestXYsigMatch_ok_kpt) hz_xySigM_ok_minpt -> Fill(bestXYsigMatch_ok_kpt);
    }

    // offline criteria: bad vs good
    if (bestCos2DMatch_notok_ele1pt>-990)  h_cos2DM_notok_ele1pt  -> Fill(bestCos2DMatch_notok_ele1pt);
    if (bestCos2DMatch_notok_ele2pt>-990)  h_cos2DM_notok_ele2pt  -> Fill(bestCos2DMatch_notok_ele2pt);
    if (bestCos2DMatch_notok_kpt>-990)     h_cos2DM_notok_kpt     -> Fill(bestCos2DMatch_notok_kpt);
    if (bestCos2DMatch_notok_ele1eta>-990) h_cos2DM_notok_ele1eta -> Fill(bestCos2DMatch_notok_ele1eta);
    if (bestCos2DMatch_notok_ele2eta>-990) h_cos2DM_notok_ele2eta -> Fill(bestCos2DMatch_notok_ele2eta);
    if (bestCos2DMatch_notok_keta>-990)    h_cos2DM_notok_keta    -> Fill(bestCos2DMatch_notok_keta);
    if (bestCos2DMatch_ok_ele1pt>-990)     h_cos2DM_ok_ele1pt     -> Fill(bestCos2DMatch_ok_ele1pt);
    if (bestCos2DMatch_ok_ele2pt>-990)     h_cos2DM_ok_ele2pt     -> Fill(bestCos2DMatch_ok_ele2pt);
    if (bestCos2DMatch_ok_kpt>-990)        h_cos2DM_ok_kpt        -> Fill(bestCos2DMatch_ok_kpt);
    if (bestCos2DMatch_ok_ele1eta>-990)    h_cos2DM_ok_ele1eta    -> Fill(bestCos2DMatch_ok_ele1eta);
    if (bestCos2DMatch_ok_ele2eta>-990)    h_cos2DM_ok_ele2eta    -> Fill(bestCos2DMatch_ok_ele2eta);
    if (bestCos2DMatch_ok_keta>-990)       h_cos2DM_ok_keta       -> Fill(bestCos2DMatch_ok_keta);
    if (bestCos2DMatch_notok_ele1pt>-990)  hz_cos2DM_notok_ele1pt -> Fill(bestCos2DMatch_notok_ele1pt);
    if (bestCos2DMatch_notok_ele2pt>-990)  hz_cos2DM_notok_ele2pt -> Fill(bestCos2DMatch_notok_ele2pt);
    if (bestCos2DMatch_notok_kpt>-990)     hz_cos2DM_notok_kpt    -> Fill(bestCos2DMatch_notok_kpt);
    if (bestCos2DMatch_ok_ele1pt>-990)     hz_cos2DM_ok_ele1pt    -> Fill(bestCos2DMatch_ok_ele1pt);
    if (bestCos2DMatch_ok_ele2pt>-990)     hz_cos2DM_ok_ele2pt    -> Fill(bestCos2DMatch_ok_ele2pt);
    if (bestCos2DMatch_ok_kpt>-990)        hz_cos2DM_ok_kpt       -> Fill(bestCos2DMatch_ok_kpt);


    // Minimal preselection  
    if (bestXYsigMatch_notok_lptmva2>-990) {
      tot_notok_all++; 
      if ( bestXYsigMatch_notok_lptmva1<20 && bestXYsigMatch_notok_lptmva1<=-2 ) continue;
      if ( bestXYsigMatch_notok_lptmva2<20 && bestXYsigMatch_notok_lptmva2<=-4 ) continue;
      if ( bestXYsigMatch_notok_pfmva1<20 && bestXYsigMatch_notok_pfmva1<=-4 )   continue;
      if ( bestXYsigMatch_notok_pfmva2<20 && bestXYsigMatch_notok_pfmva2<=-4 )   continue;
      tot_notok_presel++;	  
    }    
    if (bestXYsigMatch_ok_lptmva2>-990) {
      tot_ok_all++; 
      if ( bestXYsigMatch_ok_lptmva1<20 && bestXYsigMatch_ok_lptmva1<=-2 ) continue;
      if ( bestXYsigMatch_ok_lptmva2<20 && bestXYsigMatch_ok_lptmva2<=-4 ) continue;
      if ( bestXYsigMatch_ok_pfmva1<20 && bestXYsigMatch_ok_pfmva1<=-4 )   continue;
      if ( bestXYsigMatch_ok_pfmva2<20 && bestXYsigMatch_ok_pfmva2<=-4 )   continue;
      tot_ok_presel++;	  
    }    

  } // Loop over entries


  // Efficiency studies
  cout << endl;
  cout << "Efficiency" << endl;
  cout << "Matched:   tot_ok_presel/tot_notok_all = " << tot_ok_presel/tot_ok_all << endl; 
  cout << "UNMatched: tot_notok_presel/tot_notok_all = " << tot_notok_presel/tot_notok_all << endl; 
  cout << endl;
  cout << endl;

  // Plots
  gStyle->SetOptStat(0);

  TCanvas c0("c0","",1);
  h_debug_svProbMatch   -> SetLineColor(2);
  h_debug_svProbMatch   -> SetLineWidth(2);
  h_debug_svProbMatch   -> GetXaxis()-> SetTitle("SV Fit probability");
  h_debug_svProbMatch   -> SetTitle("");
  h_debug_svProbUnMatch -> SetLineColor(4);
  h_debug_svProbUnMatch -> SetLineWidth(2);
  h_debug_svProbUnMatch -> GetXaxis()-> SetTitle("SV Fit probability");
  h_debug_svProbUnMatch -> SetTitle("");
  h_debug_svProbUnMatch   -> DrawNormalized();
  h_debug_svProbMatch -> DrawNormalized("same");
  c0.SaveAs("DebugSVProb.png");

  TCanvas c0b("c0b","",1);
  h_debug_pf_svProbMatch   -> SetLineColor(2);
  h_debug_pf_svProbMatch   -> SetLineWidth(2);
  h_debug_pf_svProbMatch   -> GetXaxis()-> SetTitle("SV Fit probability");
  h_debug_pf_svProbMatch   -> SetTitle("");
  h_debug_pf_svProbUnMatch -> SetLineColor(4);
  h_debug_pf_svProbUnMatch -> SetLineWidth(2);
  h_debug_pf_svProbUnMatch -> GetXaxis()-> SetTitle("SV Fit probability");
  h_debug_pf_svProbUnMatch -> SetTitle("");
  h_debug_pf_svProbUnMatch   -> DrawNormalized();
  h_debug_pf_svProbMatch -> DrawNormalized("same");
  c0b.SaveAs("DebugPFSVProb.png");
  
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

  TCanvas c4("c4","",1);
  h_trueB_dRminmaxEle -> SetLineColor(2);
  h_trueB_dRminmaxEle -> SetLineWidth(2);
  h_trueB_dRminmaxEle -> Draw();
  h_trueB_dRminmaxEle -> GetXaxis()->SetTitle("max (min #DeltaR, ele)");
  h_trueB_dRminmaxEle -> SetTitle("");
  c4.SetLogy();
  c4.SaveAs("MaxMinDrEle_trueBs.png");

  TCanvas c5("c5","",1);
  h_trueB_dRmaxTrack -> SetLineColor(2);
  h_trueB_dRmaxTrack -> SetLineWidth(2);
  h_trueB_dRmaxTrack -> Draw();
  h_trueB_dRmaxTrack -> GetXaxis()->SetTitle("max #DeltaR, K");
  h_trueB_dRmaxTrack -> SetTitle("");
  c5.SetLogy();
  c5.SaveAs("MaxDrK_trueBs.png");

  TCanvas c4a("c4a","",1);
  h_combB_dRminmaxEle -> SetLineColor(2);
  h_combB_dRminmaxEle -> SetLineWidth(2);
  h_combB_dRminmaxEle -> Draw();
  h_combB_dRminmaxEle -> GetXaxis()->SetTitle("max (min #DeltaR, ele)");
  h_combB_dRminmaxEle -> SetTitle("");
  c4a.SetLogy();
  c4a.SaveAs("MaxMinDrEle_combBs.png");

  TCanvas c5a("c5a","",1);
  h_combB_dRmaxTrack -> SetLineColor(2);
  h_combB_dRmaxTrack -> SetLineWidth(2);
  h_combB_dRmaxTrack -> Draw();
  h_combB_dRmaxTrack -> GetXaxis()->SetTitle("max #DeltaR, K");
  h_combB_dRmaxTrack -> SetTitle("");
  c5a.SetLogy();
  c5a.SaveAs("MaxDrK_combBs.png");

  TCanvas c6("c6","",1);
  h_trueB_dRGenEle -> SetLineColor(2);
  h_trueB_dRGenEle -> SetLineWidth(2);
  h_trueB_dRGenEle -> Draw();
  h_trueB_dRGenEle -> GetXaxis()->SetTitle("#DeltaR (gen ele1, gen ele2)");
  h_trueB_dRGenEle -> SetTitle("");
  c6.SetLogy();
  c6.SaveAs("DrGenEle_TrueBs.png");


  h_trueBestMatch_Bmass           -> SetLineColor(1);
  h_trueBestMatch_BmassPFPF       -> SetLineColor(2);
  h_trueBestMatch_BmassPFLP       -> SetLineColor(3);
  h_trueBestMatch_BmassLPLP       -> SetLineColor(4);
  h_trueBestMatch_Bmass_PFPForPFLP_ele2PF    -> SetLineColor(2);
  h_trueBestMatch_Bmass_PFLPorLPLP_ele2PtGt2 -> SetLineColor(3);
  h_trueBestMatch_Bmass_PFLPorLPLP_ele2PtLt2 -> SetLineColor(4);
  h_trueBestMatch_BmassPFLP_ptgt2 -> SetLineColor(5);
  h_trueBestMatch_BmassPFLP_ptlt2 -> SetLineColor(6);
  h_trueBestMatch_BmassLPLP_ptgt2 -> SetLineColor(7);
  h_trueBestMatch_BmassLPLP_ptlt2 -> SetLineColor(8);
  h_trueBestMatch_BmassPFLP_idgt2 -> SetLineColor(5);
  h_trueBestMatch_BmassPFLP_idlt2 -> SetLineColor(6);
  h_trueBestMatch_BmassLPLP_idgt2 -> SetLineColor(7);
  h_trueBestMatch_BmassLPLP_idlt2 -> SetLineColor(8);
  h_trueBestMatch_Bmass           -> SetLineWidth(2);
  h_trueBestMatch_BmassPFPF       -> SetLineWidth(2);
  h_trueBestMatch_BmassPFLP       -> SetLineWidth(2);
  h_trueBestMatch_BmassLPLP       -> SetLineWidth(2);
  h_trueBestMatch_Bmass_PFPForPFLP_ele2PF    -> SetLineWidth(2);
  h_trueBestMatch_Bmass_PFLPorLPLP_ele2PtGt2 -> SetLineWidth(2);
  h_trueBestMatch_Bmass_PFLPorLPLP_ele2PtLt2 -> SetLineWidth(2);
  h_trueBestMatch_BmassPFLP_ptgt2 -> SetLineWidth(2);
  h_trueBestMatch_BmassPFLP_ptlt2 -> SetLineWidth(2);
  h_trueBestMatch_BmassLPLP_ptgt2 -> SetLineWidth(2);
  h_trueBestMatch_BmassLPLP_ptlt2 -> SetLineWidth(2);
  h_trueBestMatch_BmassPFLP_idgt2 -> SetLineWidth(2);
  h_trueBestMatch_BmassPFLP_idlt2 -> SetLineWidth(2);
  h_trueBestMatch_BmassLPLP_idgt2 -> SetLineWidth(2);
  h_trueBestMatch_BmassLPLP_idlt2 -> SetLineWidth(2);
  h_trueBestMatch_Bmass           -> GetXaxis()->SetTitle("B mass");
  h_trueBestMatch_BmassPFPF       -> GetXaxis()->SetTitle("B mass");
  h_trueBestMatch_BmassPFLP       -> GetXaxis()->SetTitle("B mass");
  h_trueBestMatch_BmassLPLP       -> GetXaxis()->SetTitle("B mass");
  h_trueBestMatch_Bmass_PFPForPFLP_ele2PF    -> GetXaxis()->SetTitle("B mass");
  h_trueBestMatch_Bmass_PFLPorLPLP_ele2PtGt2 -> GetXaxis()->SetTitle("B mass");
  h_trueBestMatch_Bmass_PFLPorLPLP_ele2PtLt2 -> GetXaxis()->SetTitle("B mass");
  h_trueBestMatch_BmassPFLP_ptgt2 -> GetXaxis()->SetTitle("B mass");
  h_trueBestMatch_BmassPFLP_ptlt2 -> GetXaxis()->SetTitle("B mass");
  h_trueBestMatch_BmassLPLP_ptgt2 -> GetXaxis()->SetTitle("B mass");
  h_trueBestMatch_BmassLPLP_ptlt2 -> GetXaxis()->SetTitle("B mass");
  h_trueBestMatch_BmassPFLP_idgt2 -> GetXaxis()->SetTitle("B mass");
  h_trueBestMatch_BmassPFLP_idlt2 -> GetXaxis()->SetTitle("B mass");
  h_trueBestMatch_BmassLPLP_idgt2 -> GetXaxis()->SetTitle("B mass");
  h_trueBestMatch_BmassLPLP_idlt2 -> GetXaxis()->SetTitle("B mass");
  h_trueBestMatch_Bmass           -> SetTitle("");
  h_trueBestMatch_BmassPFPF       -> SetTitle("");
  h_trueBestMatch_BmassPFLP       -> SetTitle("");
  h_trueBestMatch_BmassLPLP       -> SetTitle("");
  h_trueBestMatch_Bmass_PFPForPFLP_ele2PF    -> SetTitle("");
  h_trueBestMatch_Bmass_PFLPorLPLP_ele2PtGt2 -> SetTitle("");
  h_trueBestMatch_Bmass_PFLPorLPLP_ele2PtLt2 -> SetTitle("");
  h_trueBestMatch_BmassPFLP_ptgt2 -> SetTitle("");
  h_trueBestMatch_BmassPFLP_ptlt2 -> SetTitle("");
  h_trueBestMatch_BmassLPLP_ptgt2 -> SetTitle("");
  h_trueBestMatch_BmassLPLP_ptlt2 -> SetTitle("");
  h_trueBestMatch_BmassPFLP_idgt2 -> SetTitle("");
  h_trueBestMatch_BmassPFLP_idlt2 -> SetTitle("");
  h_trueBestMatch_BmassLPLP_idgt2 -> SetTitle("");
  h_trueBestMatch_BmassLPLP_idlt2 -> SetTitle("");

  TCanvas c00a("c00a","",1);
  h_trueBestMatch_BmassPFPF       -> DrawNormalized(); 
  h_trueBestMatch_BmassPFLP       -> DrawNormalized("same");
  h_trueBestMatch_BmassLPLP       -> DrawNormalized("same");
  TLegend leg00a (0.15,0.7,0.45,0.9);
  leg00a.SetFillColor(0);
  leg00a.SetFillStyle(0);
  leg00a.SetBorderSize(0);
  leg00a.AddEntry(h_trueBestMatch_BmassPFPF,"PF-PF");
  leg00a.AddEntry(h_trueBestMatch_BmassPFLP,"PF-LP");
  leg00a.AddEntry(h_trueBestMatch_BmassLPLP,"LP-LP");
  leg00a.Draw();
  c00a.SaveAs("Bmass_cat.png");

  TCanvas c00aa("c00aa","",1);
  h_trueBestMatch_Bmass_PFPForPFLP_ele2PF    -> DrawNormalized(); 
  h_trueBestMatch_Bmass_PFLPorLPLP_ele2PtGt2 -> DrawNormalized("same");
  h_trueBestMatch_Bmass_PFLPorLPLP_ele2PtLt2 -> DrawNormalized("same");
  TLegend leg00aa (0.15,0.7,0.45,0.9);
  leg00aa.SetFillColor(0);
  leg00aa.SetFillStyle(0);
  leg00aa.SetBorderSize(0);
  leg00aa.AddEntry(h_trueBestMatch_Bmass_PFPForPFLP_ele2PF,   "ele2 PF");
  leg00aa.AddEntry(h_trueBestMatch_Bmass_PFLPorLPLP_ele2PtGt2,"ele2 LP, pT2>2");
  leg00aa.AddEntry(h_trueBestMatch_Bmass_PFLPorLPLP_ele2PtLt2,"ele2 LP, pT2<2");
  leg00aa.Draw();
  c00aa.SaveAs("Bmass_catNew.png");

  TCanvas c00b("c00b","",1);
  h_trueBestMatch_BmassPFPF       -> DrawNormalized(); 
  h_trueBestMatch_BmassPFLP_ptgt2 -> DrawNormalized("same"); 
  h_trueBestMatch_BmassLPLP_ptgt2 -> DrawNormalized("same"); 
  TLegend leg00b (0.15,0.7,0.45,0.9);
  leg00b.SetFillColor(0);
  leg00b.SetFillStyle(0);
  leg00b.SetBorderSize(0);
  leg00b.AddEntry(h_trueBestMatch_BmassPFPF,"PF-PF");
  leg00b.AddEntry(h_trueBestMatch_BmassPFLP_ptgt2,"PF - LPwithPt>2");
  leg00b.AddEntry(h_trueBestMatch_BmassLPLP_ptgt2,"LPwithPt>2 - LPwithPt>2");
  leg00b.Draw();
  c00b.SaveAs("Bmass_ptgt2.png");

  TCanvas c00b2("c00b2","",1);
  h_trueBestMatch_BmassPFPF       -> DrawNormalized(); 
  h_trueBestMatch_BmassPFLP_idgt2 -> DrawNormalized("same"); 
  h_trueBestMatch_BmassLPLP_idgt2 -> DrawNormalized("same"); 
  TLegend leg00b2 (0.15,0.7,0.45,0.9);
  leg00b2.SetFillColor(0);
  leg00b2.SetFillStyle(0);
  leg00b2.SetBorderSize(0);
  leg00b2.AddEntry(h_trueBestMatch_BmassPFPF,"PF-PF");
  leg00b2.AddEntry(h_trueBestMatch_BmassPFLP_idgt2,"PF - LPwithMva>2");
  leg00b2.AddEntry(h_trueBestMatch_BmassLPLP_idgt2,"LPwithMva>2 - LPwithMva>2");
  leg00b2.Draw();
  c00b2.SaveAs("Bmass_idgt2.png");

  TCanvas c00c("c00c","",1);
  h_trueBestMatch_BmassLPLP       -> DrawNormalized(); 
  h_trueBestMatch_BmassLPLP_ptlt2 -> DrawNormalized("same"); 
  TLegend leg00c (0.15,0.7,0.45,0.9);
  leg00c.SetFillColor(0);
  leg00c.SetFillStyle(0);
  leg00c.SetBorderSize(0);
  leg00c.AddEntry(h_trueBestMatch_BmassLPLP,"LP-LP");
  leg00c.AddEntry(h_trueBestMatch_BmassLPLP_ptlt2,"LPwithPt>2 - LPwithPt>2");
  leg00c.Draw();
  c00c.SaveAs("Bmass_ptlt2.png");

  TCanvas c00c2("c00c2","",1);
  h_trueBestMatch_BmassLPLP       -> DrawNormalized(); 
  h_trueBestMatch_BmassLPLP_idlt2 -> DrawNormalized("same"); 
  TLegend leg00c2(0.15,0.7,0.45,0.9);
  leg00c2.SetFillColor(0);
  leg00c2.SetFillStyle(0);
  leg00c2.SetBorderSize(0);
  leg00c2.AddEntry(h_trueBestMatch_BmassLPLP,"LP-LP");
  leg00c2.AddEntry(h_trueBestMatch_BmassLPLP_idlt2,"LPwithMva>2 - LPwithMva>2");
  leg00c2.Draw();
  c00c2.SaveAs("Bmass_idlt2.png");

  TCanvas c00d("c00d","",1);
  h_trueBestMatch_BmassPFLP_ptlt2 -> DrawNormalized(); 
  h_trueBestMatch_BmassPFLP_ptgt2 -> DrawNormalized("same"); 
  TLegend leg00d (0.15,0.7,0.45,0.9);
  leg00d.SetFillColor(0);
  leg00d.SetFillStyle(0);
  leg00d.SetBorderSize(0);
  leg00d.AddEntry(h_trueBestMatch_BmassPFLP_ptgt2,"PF-LP, LPwithPt>2");
  leg00d.AddEntry(h_trueBestMatch_BmassPFLP_ptlt2,"PF-LP, LPwithPt<2");
  leg00d.Draw();
  c00d.SaveAs("Bmass_PFLP_ptsplit.png");

  TCanvas c00d2("c00d2","",1);
  h_trueBestMatch_BmassPFLP_idgt2 -> DrawNormalized(); 
  h_trueBestMatch_BmassPFLP_idlt2 -> DrawNormalized("same"); 
  TLegend leg00d2(0.15,0.7,0.45,0.9);
  leg00d2.SetFillColor(0);
  leg00d2.SetFillStyle(0);
  leg00d2.SetBorderSize(0);
  leg00d2.AddEntry(h_trueBestMatch_BmassPFLP_idgt2,"PF-LP, LPwithId>2");
  leg00d2.AddEntry(h_trueBestMatch_BmassPFLP_idlt2,"PF-LP, LPwithId<2");
  leg00d2.Draw();
  c00d2.SaveAs("Bmass_PFLP_idsplit.png");

  TCanvas c00e("c00e","",1);
  h_trueBestMatch_BmassLPLP_ptlt2 -> DrawNormalized(); 
  h_trueBestMatch_BmassLPLP_ptgt2 -> DrawNormalized("same"); 
  TLegend leg00e (0.15,0.7,0.45,0.9);
  leg00e.SetFillColor(0);
  leg00e.SetFillStyle(0);
  leg00e.SetBorderSize(0);
  leg00e.AddEntry(h_trueBestMatch_BmassLPLP_ptgt2,"LP-LP, LPwithPt>2");
  leg00e.AddEntry(h_trueBestMatch_BmassLPLP_ptlt2,"LP-LP, LPwithPt<2");
  leg00e.Draw();
  c00e.SaveAs("Bmass_LPLP_ptsplit.png");

  TCanvas c00e2("c00e2","",1);
  h_trueBestMatch_BmassLPLP_idgt2 -> DrawNormalized(); 
  h_trueBestMatch_BmassLPLP_idlt2 -> DrawNormalized("same"); 
  TLegend leg00e2(0.15,0.7,0.45,0.9);
  leg00e2.SetFillColor(0);
  leg00e2.SetFillStyle(0);
  leg00e2.SetBorderSize(0);
  leg00e2.AddEntry(h_trueBestMatch_BmassLPLP_idgt2,"LP-LP, LPwithId>2");
  leg00e2.AddEntry(h_trueBestMatch_BmassLPLP_idlt2,"LP-LP, LPwithId<2");
  leg00e2.Draw();
  c00e2.SaveAs("Bmass_LPLP_idsplit.png");

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
  h_trueBestMatch_eleptsum -> SetLineColor(2);
  h_trueBestMatch_eleptsum -> SetLineWidth(2);
  h_trueBestMatch_eleptsum -> GetXaxis()->SetTitle("ele1 pT + ele2 pT");
  h_trueBestMatch_eleptsum -> SetTitle("");
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
  h_trueBestMatch_ele2eta_lowpt_ptlt1 -> SetLineColor(2);
  h_trueBestMatch_ele2eta_lowpt_ptlt1 -> SetLineWidth(2);
  h_trueBestMatch_ele2eta_lowpt_ptlt1 -> GetXaxis()->SetTitle("ele2 #eta");
  h_trueBestMatch_ele2eta_lowpt_ptlt1 -> SetTitle("");
  h_trueBestMatch_ele2eta_lowpt_ptlt1d5 -> SetLineColor(2);
  h_trueBestMatch_ele2eta_lowpt_ptlt1d5 -> SetLineWidth(2);
  h_trueBestMatch_ele2eta_lowpt_ptlt1d5 -> GetXaxis()->SetTitle("ele2 #eta");
  h_trueBestMatch_ele2eta_lowpt_ptlt1d5 -> SetTitle("");
  h_trueBestMatch_ele2eta_lowpt_ptlt2 -> SetLineColor(2);
  h_trueBestMatch_ele2eta_lowpt_ptlt2 -> SetLineWidth(2);
  h_trueBestMatch_ele2eta_lowpt_ptlt2 -> GetXaxis()->SetTitle("ele2 #eta");
  h_trueBestMatch_ele2eta_lowpt_ptlt2 -> SetTitle("");
  h_trueBestMatch_ele2eta_lowpt_idlt0 -> SetLineColor(2);
  h_trueBestMatch_ele2eta_lowpt_idlt0 -> SetLineWidth(2);
  h_trueBestMatch_ele2eta_lowpt_idlt0 -> GetXaxis()->SetTitle("ele2 #eta");
  h_trueBestMatch_ele2eta_lowpt_idlt0 -> SetTitle("");
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
  h_trueBestMatch_maxDrRecoGen   -> SetLineColor(2);
  h_trueBestMatch_maxDrRecoGen   -> SetLineWidth(2);
  h_trueBestMatch_maxDrRecoGen   -> GetXaxis()->SetTitle("max #DeltaR (gen, reco ele)");
  h_trueBestMatch_maxDrRecoGen   -> SetTitle("");
  h_trueBestMatch_minDrRecoGen   -> SetLineColor(2);
  h_trueBestMatch_minDrRecoGen   -> SetLineWidth(2);
  h_trueBestMatch_minDrRecoGen   -> GetXaxis()->SetTitle("min #DeltaR (gen, reco ele)");
  h_trueBestMatch_minDrRecoGen   -> SetTitle("");
  h_trueBestMatch_drRecoGenTrack -> SetLineColor(2);
  h_trueBestMatch_drRecoGenTrack -> SetLineWidth(2);
  h_trueBestMatch_drRecoGenTrack -> GetXaxis()->SetTitle("#DeltaR (gen, reco K)");
  h_trueBestMatch_drRecoGenTrack -> SetTitle("");
  //
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
  h_comb_eleptsum -> SetLineColor(4);
  h_comb_eleptsum -> SetLineWidth(2);
  h_comb_eleptsum -> GetXaxis()->SetTitle("ele1 pT + ele2 pT");
  h_comb_eleptsum -> SetTitle("");
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
  h_comb_ele2eta_lowpt_ptlt1 -> SetLineColor(4);
  h_comb_ele2eta_lowpt_ptlt1 -> SetLineWidth(2);
  h_comb_ele2eta_lowpt_ptlt1 -> GetXaxis()->SetTitle("ele2 #eta");
  h_comb_ele2eta_lowpt_ptlt1 -> SetTitle("");
  h_comb_ele2eta_lowpt_ptlt1d5 -> SetLineColor(4);
  h_comb_ele2eta_lowpt_ptlt1d5 -> SetLineWidth(2);
  h_comb_ele2eta_lowpt_ptlt1d5 -> GetXaxis()->SetTitle("ele2 #eta");
  h_comb_ele2eta_lowpt_ptlt1d5 -> SetTitle("");
  h_comb_ele2eta_lowpt_ptlt2 -> SetLineColor(4);
  h_comb_ele2eta_lowpt_ptlt2 -> SetLineWidth(2);
  h_comb_ele2eta_lowpt_ptlt2 -> GetXaxis()->SetTitle("ele2 #eta");
  h_comb_ele2eta_lowpt_ptlt2 -> SetTitle("");
  h_comb_ele2eta_lowpt_idlt0 -> SetLineColor(4);
  h_comb_ele2eta_lowpt_idlt0 -> SetLineWidth(2);
  h_comb_ele2eta_lowpt_idlt0 -> GetXaxis()->SetTitle("ele2 #eta");
  h_comb_ele2eta_lowpt_idlt0 -> SetTitle("");
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
  h_comb_maxDrRecoGen   -> SetLineColor(4);
  h_comb_maxDrRecoGen   -> SetLineWidth(2);
  h_comb_maxDrRecoGen   -> GetXaxis()->SetTitle("max #DeltaR (gen, reco ele)");
  h_comb_maxDrRecoGen   -> SetTitle("");
  h_comb_minDrRecoGen   -> SetLineColor(4);
  h_comb_minDrRecoGen   -> SetLineWidth(2);
  h_comb_minDrRecoGen   -> GetXaxis()->SetTitle("min #DeltaR (gen, reco ele)");
  h_comb_minDrRecoGen   -> SetTitle("");
  h_comb_drRecoGenTrack -> SetLineColor(4);
  h_comb_drRecoGenTrack -> SetLineWidth(2);
  h_comb_drRecoGenTrack -> GetXaxis()->SetTitle("#DeltaR (gen, reco K)");
  h_comb_drRecoGenTrack -> SetTitle("");
  //
  TLegend leg (0.55,0.7,0.95,0.9);
  leg.SetFillColor(0);
  leg.SetFillStyle(0);
  leg.SetBorderSize(0);
  leg.AddEntry(h_trueBestMatch_svProb,"Best match to MC-truth");
  leg.AddEntry(h_comb_svProb,"Combinatorics");
  //
  TCanvas c10("c10","",1);
  h_comb_svProb->DrawNormalized();
  h_trueBestMatch_svProb->DrawNormalized("same");
  leg.Draw("same");  
  c10.SaveAs("TrueVsComb_svprob.png");
  //
  TCanvas c11("c11","",1);
  h_comb_xysig->DrawNormalized();
  h_trueBestMatch_xysig->DrawNormalized("same");
  leg.Draw("same");  
  c11.SetLogy();
  c11.SaveAs("TrueVsComb_dxy.png");
  //
  TCanvas c12("c12","",1);
  h_trueBestMatch_cos2d->DrawNormalized();
  h_comb_cos2d->DrawNormalized("same");
  leg.Draw("same");  
  c12.SetLogy();
  c12.SaveAs("TrueVsComb_cos2d.png");
  //
  TCanvas c13("c13","",1);
  h_comb_eleptsum->DrawNormalized();
  h_trueBestMatch_eleptsum->DrawNormalized("same");
  leg.Draw("same");  
  c13.SetLogy();
  c13.SaveAs("TrueVsComb_eleptsum.png");
  //
  TCanvas c14("c14","",1);
  h_comb_kpt->DrawNormalized();
  h_trueBestMatch_kpt->DrawNormalized("same");
  leg.Draw("same");  
  c14.SetLogy();
  c14.SaveAs("TrueVsComb_kpt.png");
  //
  TCanvas c14a("c14a","",1);
  h_trueBestMatch_keta->DrawNormalized();
  h_comb_keta->DrawNormalized("same");
  leg.Draw("same");  
  c14a.SaveAs("TrueVsComb_keta.png");
  //
  TCanvas c14b("c14b","",1);
  h_trueBestMatch_ele1eta->DrawNormalized();
  h_comb_ele1eta->DrawNormalized("same");
  leg.Draw("same");  
  c14b.SaveAs("TrueVsComb_ele1eta.png");
  //
  TCanvas c14c("c14c","",1);
  h_trueBestMatch_ele2eta->DrawNormalized();
  h_comb_ele2eta->DrawNormalized("same");
  leg.Draw("same");  
  c14c.SaveAs("TrueVsComb_ele2eta.png");
  //
  TCanvas c14cc("c14cc","",1);
  h_trueBestMatch_ele2eta_lowpt_ptlt1->DrawNormalized();
  h_comb_ele2eta_lowpt_ptlt1->DrawNormalized("same");
  leg.Draw("same");  
  c14cc.SaveAs("TrueVsComb_ele2eta_lowpt_ptlt1.png");
  //
  TCanvas c14ccc("c14ccc","",1);
  h_trueBestMatch_ele2eta_lowpt_ptlt1d5->DrawNormalized();
  h_comb_ele2eta_lowpt_ptlt1d5->DrawNormalized("same");
  leg.Draw("same");  
  c14ccc.SaveAs("TrueVsComb_ele2eta_lowpt_ptlt1d5.png");
  //
  TCanvas c14cccc("c14cccc","",1);
  h_trueBestMatch_ele2eta_lowpt_ptlt2->DrawNormalized();
  h_comb_ele2eta_lowpt_ptlt2->DrawNormalized("same");
  leg.Draw("same");  
  c14cccc.SaveAs("TrueVsComb_ele2eta_lowpt_ptlt2.png");
  //
  TCanvas c14ccccc("c14ccccc","",1);
  h_trueBestMatch_ele2eta_lowpt_idlt0->DrawNormalized();
  h_comb_ele2eta_lowpt_idlt0->DrawNormalized("same");
  leg.Draw("same");  
  c14ccccc.SaveAs("TrueVsComb_ele2eta_lowpt_idlt0.png");
  //
  TCanvas c14d("c14d","",1);
  h_comb_ele1pt->DrawNormalized();
  h_trueBestMatch_ele1pt->DrawNormalized("same");
  leg.Draw("same");  
  c14d.SetLogy();
  c14d.SaveAs("TrueVsComb_ele1pt.png");
  //
  TCanvas c14e("c14e","",1);
  h_comb_ele2pt->DrawNormalized();
  h_trueBestMatch_ele2pt->DrawNormalized("same");
  leg.Draw("same");  
  c14e.SetLogy();
  c14e.SaveAs("TrueVsComb_ele2pt.png");
  //
  TCanvas c14f("c14f","",1);
  h_trueBestMatch_pfmva1->DrawNormalized();
  h_comb_pfmva1->DrawNormalized("same");
  leg.Draw("same");  
  c14f.SaveAs("TrueVsComb_ele1pfmvaId.png");
  //
  TCanvas c14fbis("c14fbis","",1);
  h_trueBestMatch_pfmva1->DrawNormalized();
  h_comb_pfmva1_badele1->DrawNormalized("same");
  leg.Draw("same");  
  c14fbis.SaveAs("TrueVsComb_ele1pfmvaId_ele1bad.png");
  //
  TCanvas c14g("c14g","",1);
  h_trueBestMatch_pfmva2->DrawNormalized();
  h_comb_pfmva2->DrawNormalized("same");
  leg.Draw("same");  
  c14g.SaveAs("TrueVsComb_ele2pfmvaId.png");
  //
  TCanvas c14gbis("c14gbis","",1);
  h_trueBestMatch_pfmva2->DrawNormalized();
  h_comb_pfmva2_badele2->DrawNormalized("same");
  leg.Draw("same");  
  c14gbis.SaveAs("TrueVsComb_ele2pfmvaId_ele2bad.png");
  //
  TCanvas c14h("c14h","",1);
  h_trueBestMatch_lptmva1->DrawNormalized();
  h_comb_lptmva1->DrawNormalized("same");
  leg.Draw("same");  
  c14h.SaveAs("TrueVsComb_ele1lptmvaId.png");
  //
  TCanvas c14hbis("c14hbis","",1);
  h_trueBestMatch_lptmva1->DrawNormalized();
  h_comb_lptmva1_badele1->DrawNormalized("same");
  leg.Draw("same");  
  c14hbis.SaveAs("TrueVsComb_ele1lptmvaId_ele1bad.png");
  //
  TCanvas c14i("c14i","",1);
  h_trueBestMatch_lptmva2->DrawNormalized();
  h_comb_lptmva2->DrawNormalized("same");
  leg.Draw("same");  
  c14i.SaveAs("TrueVsComb_ele2lptmvaId.png");
  //
  TCanvas c14ibis("c14ibis","",1);
  h_trueBestMatch_lptmva2->DrawNormalized();
  h_comb_lptmva2_badele2->DrawNormalized("same");
  leg.Draw("same");  
  c14ibis.SaveAs("TrueVsComb_ele2lptmvaId_ele2bad.png");
  //
  TCanvas c141("c141","",1);
  hz_comb_kpt->DrawNormalized();
  hz_trueBestMatch_kpt->DrawNormalized("same");
  leg.Draw("same");  
  c141.SaveAs("TrueVsComb_kpt_zoom.png");
  //
  TCanvas c141d("c141d","",1);
  hz_comb_ele1pt->DrawNormalized();
  hz_trueBestMatch_ele1pt->DrawNormalized("same");
  leg.Draw("same");  
  c141d.SaveAs("TrueVsComb_ele1pt_zoom.png");
  //
  TCanvas c141e("c141e","",1);
  hz_comb_ele2pt->DrawNormalized();
  hz_trueBestMatch_ele2pt->DrawNormalized("same");
  leg.Draw("same");  
  c141e.SaveAs("TrueVsComb_ele2pt_zoom.png");
  //
  TCanvas c141bis("c141bis","",1);
  hz_comb_kpt_badk->DrawNormalized();
  hz_trueBestMatch_kpt->DrawNormalized("same");
  leg.Draw("same");  
  c141bis.SaveAs("TrueVsComb_kpt_badKonly_zoom.png");
  //
  TCanvas c141dbis("c141dbis","",1);
  hz_comb_ele1pt_badele1->DrawNormalized();
  hz_trueBestMatch_ele1pt->DrawNormalized("same");
  leg.Draw("same");  
  c141dbis.SaveAs("TrueVsComb_ele1pt_badele1only_zoom.png");
  //
  TCanvas c141ebis("c141ebis","",1);
  hz_comb_ele2pt_badele2->DrawNormalized();
  hz_trueBestMatch_ele2pt->DrawNormalized("same");
  leg.Draw("same");  
  c141ebis.SaveAs("TrueVsComb_ele2pt_badele2only_zoom.png");
  //
  TCanvas c141f("c141f","",1);
  hz_comb_minpt->DrawNormalized();
  hz_trueBestMatch_minpt->DrawNormalized("same");
  leg.Draw("same");  
  c141f.SaveAs("TrueVsComb_minpt_zoom.png");
  //
  TCanvas c15("c15","",1);
  h_trueBestMatch_maxDrRecoGen->DrawNormalized();
  h_comb_maxDrRecoGen->DrawNormalized("same");
  leg.Draw("same");  
  c15.SetLogy();
  c15.SaveAs("TrueVsComb_maxDrRecoGen.png");
  //
  TCanvas c16("c16","",1);
  h_trueBestMatch_minDrRecoGen->DrawNormalized();
  h_comb_minDrRecoGen->DrawNormalized("same");
  leg.Draw("same");  
  c16.SetLogy();
  c16.SaveAs("TrueVsComb_minDrRecoGen.png");
  //
  TCanvas c17("c17","",1);
  h_trueBestMatch_drRecoGenTrack->DrawNormalized();
  h_comb_drRecoGenTrack->DrawNormalized("same");
  leg.Draw("same");  
  c17.SetLogy();
  c17.SaveAs("TrueVsComb_drRecoGenTrack.png");
  //
  TCanvas c18("c18","",1);
  hz_trueBestMatch_ele2Vskpt->Draw("colz");
  c18.SaveAs("BestMatch_ele2vsKpt.png");
  //
  TCanvas c18b("c18b","",1);
  hz_comb_ele2Vskpt->Draw("colz");
  c18b.SaveAs("Comb_ele2vsKpt.png");
  //
  //
  // Best SV prob
  h_svProbM_notok_ele1pt  -> SetLineColor(4);
  h_svProbM_notok_ele1pt  -> SetLineWidth(2);
  h_svProbM_notok_ele1pt  -> GetXaxis()->SetTitle("ele1 pT");
  h_svProbM_notok_ele1pt  -> SetTitle("");
  h_svProbM_notok_ele2pt  -> SetLineColor(4);
  h_svProbM_notok_ele2pt  -> SetLineWidth(2);
  h_svProbM_notok_ele2pt  -> GetXaxis()->SetTitle("ele2 pT");
  h_svProbM_notok_ele2pt  -> SetTitle("");
  h_svProbM_notok_kpt     -> SetLineColor(4);
  h_svProbM_notok_kpt     -> SetLineWidth(2);
  h_svProbM_notok_kpt     -> GetXaxis()->SetTitle("k pT");
  h_svProbM_notok_kpt     -> SetTitle("");
  h_svProbM_notok_ele1eta -> SetLineColor(4);
  h_svProbM_notok_ele1eta -> SetLineWidth(2);
  h_svProbM_notok_ele1eta -> GetXaxis()->SetTitle("ele1 #eta");
  h_svProbM_notok_ele1eta -> SetTitle("");
  h_svProbM_notok_ele2eta -> SetLineColor(4);
  h_svProbM_notok_ele2eta -> SetLineWidth(2);
  h_svProbM_notok_ele2eta -> GetXaxis()->SetTitle("ele2 #eta");
  h_svProbM_notok_ele2eta -> SetTitle("");
  h_svProbM_notok_keta    -> SetLineColor(4);
  h_svProbM_notok_keta    -> SetLineWidth(2);
  h_svProbM_notok_keta    -> GetXaxis()->SetTitle("k #eta");
  h_svProbM_notok_keta    -> SetTitle("");
  hz_svProbM_notok_ele1pt -> SetLineColor(4);
  hz_svProbM_notok_ele1pt -> SetLineWidth(2);
  hz_svProbM_notok_ele1pt -> GetXaxis()->SetTitle("ele1 pT");
  hz_svProbM_notok_ele1pt -> SetTitle("");
  hz_svProbM_notok_ele2pt -> SetLineColor(4);
  hz_svProbM_notok_ele2pt -> SetLineWidth(2);
  hz_svProbM_notok_ele2pt -> GetXaxis()->SetTitle("ele2 pT");
  hz_svProbM_notok_ele2pt -> SetTitle("");
  hz_svProbM_notok_kpt    -> SetLineColor(4);
  hz_svProbM_notok_kpt    -> SetLineWidth(2);
  hz_svProbM_notok_kpt    -> GetXaxis()->SetTitle("k pT");
  hz_svProbM_notok_kpt    -> SetTitle("");
  //
  h_svProbM_ok_ele1pt  -> SetLineColor(2);
  h_svProbM_ok_ele1pt  -> SetLineWidth(2);
  h_svProbM_ok_ele1pt  -> GetXaxis()->SetTitle("ele1 pT");
  h_svProbM_ok_ele1pt  -> SetTitle("");
  h_svProbM_ok_ele2pt  -> SetLineColor(2);
  h_svProbM_ok_ele2pt  -> SetLineWidth(2);
  h_svProbM_ok_ele2pt  -> GetXaxis()->SetTitle("ele2 pT");
  h_svProbM_ok_ele2pt  -> SetTitle("");
  h_svProbM_ok_kpt     -> SetLineColor(2);
  h_svProbM_ok_kpt     -> SetLineWidth(2);
  h_svProbM_ok_kpt     -> GetXaxis()->SetTitle("k pT");
  h_svProbM_ok_kpt     -> SetTitle("");
  h_svProbM_ok_ele1eta -> SetLineColor(2);
  h_svProbM_ok_ele1eta -> SetLineWidth(2);
  h_svProbM_ok_ele1eta -> GetXaxis()->SetTitle("ele1 #eta");
  h_svProbM_ok_ele1eta -> SetTitle("");
  h_svProbM_ok_ele2eta -> SetLineColor(2);
  h_svProbM_ok_ele2eta -> SetLineWidth(2);
  h_svProbM_ok_ele2eta -> GetXaxis()->SetTitle("ele2 #eta");
  h_svProbM_ok_ele2eta -> SetTitle("");
  h_svProbM_ok_keta    -> SetLineColor(2);
  h_svProbM_ok_keta    -> SetLineWidth(2);
  h_svProbM_ok_keta    -> GetXaxis()->SetTitle("k #eta");
  h_svProbM_ok_keta    -> SetTitle("");
  hz_svProbM_ok_ele1pt -> SetLineColor(2);
  hz_svProbM_ok_ele1pt -> SetLineWidth(2);
  hz_svProbM_ok_ele1pt -> GetXaxis()->SetTitle("ele1 pT");
  hz_svProbM_ok_ele1pt -> SetTitle("");
  hz_svProbM_ok_ele2pt -> SetLineColor(2);
  hz_svProbM_ok_ele2pt -> SetLineWidth(2);
  hz_svProbM_ok_ele2pt -> GetXaxis()->SetTitle("ele2 pT");
  hz_svProbM_ok_ele2pt -> SetTitle("");
  hz_svProbM_ok_kpt    -> SetLineColor(2);
  hz_svProbM_ok_kpt    -> SetLineWidth(2);
  hz_svProbM_ok_kpt    -> GetXaxis()->SetTitle("k pT");
  hz_svProbM_ok_kpt    -> SetTitle("");
  //
  //
  // Best XYsig
  h_xySigM_notok_ele1pt  -> SetLineColor(4);
  h_xySigM_notok_ele1pt  -> SetLineWidth(2);
  h_xySigM_notok_ele1pt  -> GetXaxis()->SetTitle("ele1 pT");
  h_xySigM_notok_ele1pt  -> SetTitle("");
  h_xySigM_notok_ele2pt  -> SetLineColor(4);
  h_xySigM_notok_ele2pt  -> SetLineWidth(2);
  h_xySigM_notok_ele2pt  -> GetXaxis()->SetTitle("ele2 pT");
  h_xySigM_notok_ele2pt  -> SetTitle("");
  h_xySigM_notok_kpt     -> SetLineColor(4);
  h_xySigM_notok_kpt     -> SetLineWidth(2);
  h_xySigM_notok_kpt     -> GetXaxis()->SetTitle("k pT");
  h_xySigM_notok_kpt     -> SetTitle("");
  h_xySigM_notok_pfmva1  -> SetLineColor(4);
  h_xySigM_notok_pfmva1  -> SetLineWidth(2);
  h_xySigM_notok_pfmva1  -> GetXaxis()->SetTitle("ele1 PF MVA");
  h_xySigM_notok_pfmva1  -> SetTitle("");
  h_xySigM_notok_pfmva2  -> SetLineColor(4);
  h_xySigM_notok_pfmva2  -> SetLineWidth(2);
  h_xySigM_notok_pfmva2  -> GetXaxis()->SetTitle("ele2 PF MVA");
  h_xySigM_notok_pfmva2  -> SetTitle("");
  h_xySigM_notok_lptmva1 -> SetLineColor(4);
  h_xySigM_notok_lptmva1 -> SetLineWidth(2);
  h_xySigM_notok_lptmva1 -> GetXaxis()->SetTitle("ele1 Low pT MVA");
  h_xySigM_notok_lptmva1 -> SetTitle("");
  h_xySigM_notok_lptmva2 -> SetLineColor(4);
  h_xySigM_notok_lptmva2 -> SetLineWidth(2);
  h_xySigM_notok_lptmva2 -> GetXaxis()->SetTitle("ele2 Low pT MVA");
  h_xySigM_notok_lptmva2 -> SetTitle("");
  h_xySigM_notok_costhetaSK -> SetLineColor(4);
  h_xySigM_notok_costhetaSK -> SetLineWidth(2);
  h_xySigM_notok_costhetaSK -> GetXaxis()->SetTitle("cos(Theta*) K");
  h_xySigM_notok_costhetaSK -> SetTitle("");
  h_xySigM_notok_costhetaSKCS -> SetLineColor(4);
  h_xySigM_notok_costhetaSKCS -> SetLineWidth(2);
  h_xySigM_notok_costhetaSKCS -> GetXaxis()->SetTitle("cos(Theta*) K");
  h_xySigM_notok_costhetaSKCS -> SetTitle("");
  h_xySigM_notok_costhetaL -> SetLineColor(4);
  h_xySigM_notok_costhetaL -> SetLineWidth(2);
  h_xySigM_notok_costhetaL -> GetXaxis()->SetTitle("cos(Theta) ele1");
  h_xySigM_notok_costhetaL -> SetTitle("");
  h_xySigM_notok_ele1eta -> SetLineColor(4);
  h_xySigM_notok_ele1eta -> SetLineWidth(2);
  h_xySigM_notok_ele1eta -> GetXaxis()->SetTitle("ele1 #eta");
  h_xySigM_notok_ele1eta -> SetTitle("");
  h_xySigM_notok_ele2eta -> SetLineColor(4);
  h_xySigM_notok_ele2eta -> SetLineWidth(2);
  h_xySigM_notok_ele2eta -> GetXaxis()->SetTitle("ele2 #eta");
  h_xySigM_notok_ele2eta -> SetTitle("");
  h_xySigM_notok_keta    -> SetLineColor(4);
  h_xySigM_notok_keta    -> SetLineWidth(2);
  h_xySigM_notok_keta    -> GetXaxis()->SetTitle("k #eta");
  h_xySigM_notok_keta    -> SetTitle("");
  hz_xySigM_notok_ele1pt -> SetLineColor(4);
  hz_xySigM_notok_ele1pt -> SetLineWidth(2);
  hz_xySigM_notok_ele1pt -> GetXaxis()->SetTitle("ele1 pT");
  hz_xySigM_notok_ele1pt -> SetTitle("");
  hz_xySigM_notok_ele2pt -> SetLineColor(4);
  hz_xySigM_notok_ele2pt -> SetLineWidth(2);
  hz_xySigM_notok_ele2pt -> GetXaxis()->SetTitle("ele2 pT");
  hz_xySigM_notok_ele2pt -> SetTitle("");
  hz_xySigM_notok_kpt    -> SetLineColor(4);
  hz_xySigM_notok_kpt    -> SetLineWidth(2);
  hz_xySigM_notok_kpt    -> GetXaxis()->SetTitle("k pT");
  hz_xySigM_notok_kpt    -> SetTitle("");
  hz_xySigM_notok_minpt    -> SetLineColor(4);
  hz_xySigM_notok_minpt    -> SetLineWidth(2);
  hz_xySigM_notok_minpt    -> GetXaxis()->SetTitle("min (ele2 pT, k pT)");
  hz_xySigM_notok_minpt    -> SetTitle("");
  hz_xySigM_notok_ele2Vskpt -> GetXaxis()->SetTitle("k pT");
  hz_xySigM_notok_ele2Vskpt -> GetYaxis()->SetTitle("ele2 pT");
  hz_xySigM_notok_ele2Vskpt -> SetTitle("");
  //
  h_xySigM_ok_ele1pt  -> SetLineColor(2);
  h_xySigM_ok_ele1pt  -> SetLineWidth(2);
  h_xySigM_ok_ele1pt  -> GetXaxis()->SetTitle("ele1 pT");
  h_xySigM_ok_ele1pt  -> SetTitle("");
  h_xySigM_ok_ele2pt  -> SetLineColor(2);
  h_xySigM_ok_ele2pt  -> SetLineWidth(2);
  h_xySigM_ok_ele2pt  -> GetXaxis()->SetTitle("ele2 pT");
  h_xySigM_ok_ele2pt  -> SetTitle("");
  h_xySigM_ok_kpt     -> SetLineColor(2);
  h_xySigM_ok_kpt     -> SetLineWidth(2);
  h_xySigM_ok_kpt     -> GetXaxis()->SetTitle("k pT");
  h_xySigM_ok_kpt     -> SetTitle("");
  h_xySigM_ok_pfmva1  -> SetLineColor(2);
  h_xySigM_ok_pfmva1  -> SetLineWidth(2);
  h_xySigM_ok_pfmva1  -> GetXaxis()->SetTitle("ele1 PF MVA");
  h_xySigM_ok_pfmva1  -> SetTitle("");
  h_xySigM_ok_pfmva2  -> SetLineColor(2);
  h_xySigM_ok_pfmva2  -> SetLineWidth(2);
  h_xySigM_ok_pfmva2  -> GetXaxis()->SetTitle("ele2 PF MVA");
  h_xySigM_ok_pfmva2  -> SetTitle("");
  h_xySigM_ok_lptmva1 -> SetLineColor(2);
  h_xySigM_ok_lptmva1 -> SetLineWidth(2);
  h_xySigM_ok_lptmva1 -> GetXaxis()->SetTitle("ele1 Low pT MVA");
  h_xySigM_ok_lptmva1 -> SetTitle("");
  h_xySigM_ok_lptmva2 -> SetLineColor(2);
  h_xySigM_ok_lptmva2 -> SetLineWidth(2);
  h_xySigM_ok_lptmva2 -> GetXaxis()->SetTitle("ele2 Low pT MVA");
  h_xySigM_ok_lptmva2 -> SetTitle("");
  h_xySigM_ok_costhetaSK_gen -> SetLineColor(2);
  h_xySigM_ok_costhetaSK_gen -> SetLineWidth(2);
  h_xySigM_ok_costhetaSK_gen -> GetXaxis()->SetTitle("cos(Theta*) K");
  h_xySigM_ok_costhetaSK_gen -> SetTitle("");
  h_xySigM_ok_costhetaSK -> SetLineColor(2);
  h_xySigM_ok_costhetaSK -> SetLineWidth(2);
  h_xySigM_ok_costhetaSK -> GetXaxis()->SetTitle("cos(Theta*) K");
  h_xySigM_ok_costhetaSK -> SetTitle("");
  h_xySigM_ok_costhetaSKCS -> SetLineColor(2);
  h_xySigM_ok_costhetaSKCS -> SetLineWidth(2);
  h_xySigM_ok_costhetaSKCS -> GetXaxis()->SetTitle("cos(Theta*) K");
  h_xySigM_ok_costhetaSKCS -> SetTitle("");
  h_xySigM_ok_costhetaL -> SetLineColor(2);
  h_xySigM_ok_costhetaL -> SetLineWidth(2);
  h_xySigM_ok_costhetaL -> GetXaxis()->SetTitle("cos(Theta) ele1");
  h_xySigM_ok_costhetaL -> SetTitle("");
  h_xySigM_ok_ele1eta -> SetLineColor(2);
  h_xySigM_ok_ele1eta -> SetLineWidth(2);
  h_xySigM_ok_ele1eta -> GetXaxis()->SetTitle("ele1 #eta");
  h_xySigM_ok_ele1eta -> SetTitle("");
  h_xySigM_ok_ele2eta -> SetLineColor(2);
  h_xySigM_ok_ele2eta -> SetLineWidth(2);
  h_xySigM_ok_ele2eta -> GetXaxis()->SetTitle("ele2 #eta");
  h_xySigM_ok_ele2eta -> SetTitle("");
  h_xySigM_ok_keta    -> SetLineColor(2);
  h_xySigM_ok_keta    -> SetLineWidth(2);
  h_xySigM_ok_keta    -> GetXaxis()->SetTitle("k #eta");
  h_xySigM_ok_keta    -> SetTitle("");
  hz_xySigM_ok_ele1pt -> SetLineColor(2);
  hz_xySigM_ok_ele1pt -> SetLineWidth(2);
  hz_xySigM_ok_ele1pt -> GetXaxis()->SetTitle("ele1 pT");
  hz_xySigM_ok_ele1pt -> SetTitle("");
  hz_xySigM_ok_ele2pt -> SetLineColor(2);
  hz_xySigM_ok_ele2pt -> SetLineWidth(2);
  hz_xySigM_ok_ele2pt -> GetXaxis()->SetTitle("ele2 pT");
  hz_xySigM_ok_ele2pt -> SetTitle("");
  hz_xySigM_ok_kpt    -> SetLineColor(2);
  hz_xySigM_ok_kpt    -> SetLineWidth(2);
  hz_xySigM_ok_kpt    -> GetXaxis()->SetTitle("k pT");
  hz_xySigM_ok_kpt    -> SetTitle("");
  hz_xySigM_ok_minpt  -> SetLineColor(2);
  hz_xySigM_ok_minpt  -> SetLineWidth(2);
  hz_xySigM_ok_minpt  -> GetXaxis()->SetTitle("min(ele2, k pT)");
  hz_xySigM_ok_minpt  -> SetTitle("");
  hz_xySigM_ok_ele2Vskpt -> GetXaxis()->SetTitle("k pT");
  hz_xySigM_ok_ele2Vskpt -> GetYaxis()->SetTitle("ele2 pT");
  hz_xySigM_ok_ele2Vskpt -> SetTitle("");
  //
  //
  // Best cos2D
  h_cos2DM_notok_ele1pt  -> SetLineColor(4);
  h_cos2DM_notok_ele1pt  -> SetLineWidth(2);
  h_cos2DM_notok_ele1pt  -> GetXaxis()->SetTitle("ele1 pT");
  h_cos2DM_notok_ele1pt  -> SetTitle("");
  h_cos2DM_notok_ele2pt  -> SetLineColor(4);
  h_cos2DM_notok_ele2pt  -> SetLineWidth(2);
  h_cos2DM_notok_ele2pt  -> GetXaxis()->SetTitle("ele2 pT");
  h_cos2DM_notok_ele2pt  -> SetTitle("");
  h_cos2DM_notok_kpt     -> SetLineColor(4);
  h_cos2DM_notok_kpt     -> SetLineWidth(2);
  h_cos2DM_notok_kpt     -> GetXaxis()->SetTitle("k pT");
  h_cos2DM_notok_kpt     -> SetTitle("");
  h_cos2DM_notok_ele1eta -> SetLineColor(4);
  h_cos2DM_notok_ele1eta -> SetLineWidth(2);
  h_cos2DM_notok_ele1eta -> GetXaxis()->SetTitle("ele1 #eta");
  h_cos2DM_notok_ele1eta -> SetTitle("");
  h_cos2DM_notok_ele2eta -> SetLineColor(4);
  h_cos2DM_notok_ele2eta -> SetLineWidth(2);
  h_cos2DM_notok_ele2eta -> GetXaxis()->SetTitle("ele2 #eta");
  h_cos2DM_notok_ele2eta -> SetTitle("");
  h_cos2DM_notok_keta    -> SetLineColor(4);
  h_cos2DM_notok_keta    -> SetLineWidth(2);
  h_cos2DM_notok_keta    -> GetXaxis()->SetTitle("k #eta");
  h_cos2DM_notok_keta    -> SetTitle("");
  hz_cos2DM_notok_ele1pt -> SetLineColor(4);
  hz_cos2DM_notok_ele1pt -> SetLineWidth(2);
  hz_cos2DM_notok_ele1pt -> GetXaxis()->SetTitle("ele1 pT");
  hz_cos2DM_notok_ele1pt -> SetTitle("");
  hz_cos2DM_notok_ele2pt -> SetLineColor(4);
  hz_cos2DM_notok_ele2pt -> SetLineWidth(2);
  hz_cos2DM_notok_ele2pt -> GetXaxis()->SetTitle("ele2 pT");
  hz_cos2DM_notok_ele2pt -> SetTitle("");
  hz_cos2DM_notok_kpt    -> SetLineColor(4);
  hz_cos2DM_notok_kpt    -> SetLineWidth(2);
  hz_cos2DM_notok_kpt    -> GetXaxis()->SetTitle("k pT");
  hz_cos2DM_notok_kpt    -> SetTitle("");
  //
  h_cos2DM_ok_ele1pt  -> SetLineColor(2);
  h_cos2DM_ok_ele1pt  -> SetLineWidth(2);
  h_cos2DM_ok_ele1pt  -> GetXaxis()->SetTitle("ele1 pT");
  h_cos2DM_ok_ele1pt  -> SetTitle("");
  h_cos2DM_ok_ele2pt  -> SetLineColor(2);
  h_cos2DM_ok_ele2pt  -> SetLineWidth(2);
  h_cos2DM_ok_ele2pt  -> GetXaxis()->SetTitle("ele2 pT");
  h_cos2DM_ok_ele2pt  -> SetTitle("");
  h_cos2DM_ok_kpt     -> SetLineColor(2);
  h_cos2DM_ok_kpt     -> SetLineWidth(2);
  h_cos2DM_ok_kpt     -> GetXaxis()->SetTitle("k pT");
  h_cos2DM_ok_kpt     -> SetTitle("");
  h_cos2DM_ok_ele1eta -> SetLineColor(2);
  h_cos2DM_ok_ele1eta -> SetLineWidth(2);
  h_cos2DM_ok_ele1eta -> GetXaxis()->SetTitle("ele1 #eta");
  h_cos2DM_ok_ele1eta -> SetTitle("");
  h_cos2DM_ok_ele2eta -> SetLineColor(2);
  h_cos2DM_ok_ele2eta -> SetLineWidth(2);
  h_cos2DM_ok_ele2eta -> GetXaxis()->SetTitle("ele2 #eta");
  h_cos2DM_ok_ele2eta -> SetTitle("");
  h_cos2DM_ok_keta    -> SetLineColor(2);
  h_cos2DM_ok_keta    -> SetLineWidth(2);
  h_cos2DM_ok_keta    -> GetXaxis()->SetTitle("k #eta");
  h_cos2DM_ok_keta    -> SetTitle("");
  hz_cos2DM_ok_ele1pt -> SetLineColor(2);
  hz_cos2DM_ok_ele1pt -> SetLineWidth(2);
  hz_cos2DM_ok_ele1pt -> GetXaxis()->SetTitle("ele1 pT");
  hz_cos2DM_ok_ele1pt -> SetTitle("");
  hz_cos2DM_ok_ele2pt -> SetLineColor(2);
  hz_cos2DM_ok_ele2pt -> SetLineWidth(2);
  hz_cos2DM_ok_ele2pt -> GetXaxis()->SetTitle("ele2 pT");
  hz_cos2DM_ok_ele2pt -> SetTitle("");
  hz_cos2DM_ok_kpt    -> SetLineColor(2);
  hz_cos2DM_ok_kpt    -> SetLineWidth(2);
  hz_cos2DM_ok_kpt    -> GetXaxis()->SetTitle("k pT");
  hz_cos2DM_ok_kpt    -> SetTitle("");

  //
  TLegend leg2 (0.55,0.7,0.95,0.9);
  leg2.SetFillColor(0);
  leg2.SetFillStyle(0);
  leg2.SetBorderSize(0);
  leg2.AddEntry(h_svProbM_ok_ele1pt,   "Matched to MC-truth");
  leg2.AddEntry(h_svProbM_notok_ele1pt,"Combinatorics");
  TLegend leg3 (0.15,0.7,0.55,0.9);
  leg3.SetFillColor(0);
  leg3.SetFillStyle(0);
  leg3.SetBorderSize(0);
  leg3.AddEntry(h_svProbM_ok_ele1pt,   "Matched to MC-truth");
  leg3.AddEntry(h_svProbM_notok_ele1pt,"Combinatorics");
  //
  TCanvas c21("c21","",1);
  h_svProbM_ok_ele1pt->DrawNormalized();
  h_svProbM_notok_ele1pt->DrawNormalized("same");
  leg2.Draw("same");  
  c21.SaveAs("BestSvProb_goodVsBad_ele1pt.png");
  //
  TCanvas c22("c22","",1);
  h_svProbM_notok_ele2pt->DrawNormalized();
  h_svProbM_ok_ele2pt->DrawNormalized("same");
  leg2.Draw("same");  
  c22.SaveAs("BestSvProb_goodVsBad_ele2pt.png");
  //
  TCanvas c23("c23","",1);
  h_svProbM_notok_kpt->DrawNormalized();
  h_svProbM_ok_kpt->DrawNormalized("same");
  leg2.Draw("same");  
  c23.SaveAs("BestSvProb_goodVsBad_kpt.png");
  //
  TCanvas c24("c24","",1);
  h_svProbM_ok_ele1eta->DrawNormalized();
  h_svProbM_notok_ele1eta->DrawNormalized("same");
  leg2.Draw("same");  
  c24.SaveAs("BestSvProb_goodVsBad_ele1eta.png");
  //
  TCanvas c25("c25","",1);
  h_svProbM_ok_ele2eta->DrawNormalized();
  h_svProbM_notok_ele2eta->DrawNormalized("same");
  leg2.Draw("same");  
  c25.SaveAs("BestSvProb_goodVsBad_ele2eta.png");
  //
  TCanvas c26("c26","",1);
  h_svProbM_ok_keta->DrawNormalized();
  h_svProbM_notok_keta->DrawNormalized("same");
  leg2.Draw("same");  
  c26.SaveAs("BestSvProb_goodVsBad_keta.png");
  //
  TCanvas c21b("c21b","",1);
  hz_svProbM_ok_ele1pt->DrawNormalized();
  hz_svProbM_notok_ele1pt->DrawNormalized("same");
  leg2.Draw("same");  
  c21b.SaveAs("BestSvProb_goodVsBad_ele1ptZ.png");
  //
  TCanvas c22b("c22b","",1);
  hz_svProbM_notok_ele2pt->DrawNormalized();
  hz_svProbM_ok_ele2pt->DrawNormalized("same");
  leg2.Draw("same");  
  c22b.SaveAs("BestSvProb_goodVsBad_ele2ptZ.png");
  //
  TCanvas c23b("c23b","",1);
  hz_svProbM_notok_kpt->DrawNormalized();
  hz_svProbM_ok_kpt->DrawNormalized("same");
  leg2.Draw("same");  
  c23b.SaveAs("BestSvProb_goodVsBad_kptZ.png");
  //
  //
  //
  TCanvas c31("c31","",1);
  h_xySigM_notok_ele1pt->DrawNormalized();
  h_xySigM_ok_ele1pt->DrawNormalized("same");
  leg2.Draw("same");  
  c31.SaveAs("BestXYsig_goodVsBad_ele1pt.png");
  //
  TCanvas c32("c32","",1);
  h_xySigM_notok_ele2pt->DrawNormalized();
  h_xySigM_ok_ele2pt->DrawNormalized("same");
  leg2.Draw("same");  
  c32.SaveAs("BestXYsig_goodVsBad_ele2pt.png");
  //
  TCanvas c33("c33","",1);
  h_xySigM_notok_kpt->DrawNormalized();
  h_xySigM_ok_kpt->DrawNormalized("same");
  leg2.Draw("same");  
  c33.SaveAs("BestXYsig_goodVsBad_kpt.png");
  //
  TCanvas c34("c34","",1);
  h_xySigM_ok_ele1eta->DrawNormalized();
  h_xySigM_notok_ele1eta->DrawNormalized("same");
  leg2.Draw("same");  
  c34.SaveAs("BestXYsig_goodVsBad_ele1eta.png");
  //
  TCanvas c35("c35","",1);
  h_xySigM_ok_ele2eta->DrawNormalized();
  h_xySigM_notok_ele2eta->DrawNormalized("same");
  leg2.Draw("same");  
  c35.SaveAs("BestXYsig_goodVsBad_ele2eta.png");
  //
  TCanvas c36("c36","",1);
  h_xySigM_ok_keta->DrawNormalized();
  h_xySigM_notok_keta->DrawNormalized("same");
  leg2.Draw("same");  
  c36.SaveAs("BestXYsig_goodVsBad_keta.png");
  //
  TCanvas c37("c37","",1);
  h_xySigM_notok_pfmva1->DrawNormalized();
  h_xySigM_ok_pfmva1->DrawNormalized("same");
  leg3.Draw("same");  
  c37.SaveAs("BestXYsig_goodVsBad_pfmva1.png");
  //
  TCanvas c38("c38","",1);
  h_xySigM_notok_pfmva2->DrawNormalized();
  h_xySigM_ok_pfmva2->DrawNormalized("same");
  leg3.Draw("same");  
  c38.SaveAs("BestXYsig_goodVsBad_pfmva2.png");
  //
  TCanvas c39("c39","",1);
  h_xySigM_ok_lptmva1->DrawNormalized();
  h_xySigM_notok_lptmva1->DrawNormalized("same");
  leg3.Draw("same");  
  c39.SaveAs("BestXYsig_goodVsBad_lptmva1.png");
  //
  TCanvas c40("c40","",1);
  h_xySigM_ok_lptmva2->DrawNormalized();
  h_xySigM_notok_lptmva2->DrawNormalized("same");
  leg3.Draw("same");  
  c40.SaveAs("BestXYsig_goodVsBad_lptmva2.png");
  //
  TCanvas c41a("c41a","",1);
  h_xySigM_notok_costhetaSK -> DrawNormalized();
  h_xySigM_ok_costhetaSK -> DrawNormalized("same");
  leg3.Draw("same");
  c41a.SaveAs("BestXYsig_goodVsBad_costhetaSK.png");
  //
  TCanvas c41aa("c41aa","",1);
  h_xySigM_ok_costhetaSK_gen -> DrawNormalized();
  c41aa.SaveAs("BestXYsig_costhetaSK_genLevel.png");
  //
  TCanvas c41aaa("c41aaa","",1);
  h_xySigM_notok_costhetaL -> DrawNormalized();
  h_xySigM_ok_costhetaL    -> DrawNormalized("same");
  leg3.Draw("same");
  c41aaa.SaveAs("BestXYsig_goodVsBad_costhetaL.png");
  //
  TCanvas c41aaaa("c41aaaa","",1);
  h_xySigM_notok_costhetaSKCS -> DrawNormalized();
  h_xySigM_ok_costhetaSKCS -> DrawNormalized("same");
  leg3.Draw("same");
  c41aaaa.SaveAs("BestXYsig_goodVsBad_costhetaSKCS.png");
  //
  //
  TCanvas c31b("c31b","",1);
  hz_xySigM_ok_ele1pt->DrawNormalized();
  hz_xySigM_notok_ele1pt->DrawNormalized("same");
  leg2.Draw("same");  
  c31b.SaveAs("BestXYsig_goodVsBad_ele1ptZ.png");
  //
  TCanvas c32b("c32b","",1);
  hz_xySigM_notok_ele2pt->DrawNormalized();
  hz_xySigM_ok_ele2pt->DrawNormalized("same");
  leg2.Draw("same");  
  c32b.SaveAs("BestXYsig_goodVsBad_ele2ptZ.png");
  //
  TCanvas c33b("c33b","",1);
  hz_xySigM_notok_kpt->DrawNormalized();
  hz_xySigM_ok_kpt->DrawNormalized("same");
  leg2.Draw("same");  
  c33b.SaveAs("BestXYsig_goodVsBad_kptZ.png");
  //
  TCanvas c34b("c34b","",1);
  hz_xySigM_notok_minpt->DrawNormalized();
  hz_xySigM_ok_minpt->DrawNormalized("same");
  leg2.Draw("same");  
  c34b.SaveAs("BestXYsig_goodVsBad_minptZ.png");
  //
  TCanvas c30c("c30c","",1);
  hz_xySigM_ok_ele2Vskpt -> Draw("colz");
  c30c.SaveAs("BestXYsig_Good_ele2VsKptZ.png");
  //
  TCanvas c30c2("c30c2","",1);
  hz_xySigM_notok_ele2Vskpt -> Draw("colz");
  c30c2.SaveAs("BestXYsig_Bad_ele2VsKptZ.png");
  //
  //
  TCanvas c41("c41","",1);
  h_cos2DM_ok_ele1pt->DrawNormalized();
  h_cos2DM_notok_ele1pt->DrawNormalized("same");
  leg2.Draw("same");  
  c41.SaveAs("BestCos2D_goodVsBad_ele1pt.png");
  //
  TCanvas c42("c42","",1);
  h_cos2DM_notok_ele2pt->DrawNormalized();
  h_cos2DM_ok_ele2pt->DrawNormalized("same");
  leg2.Draw("same");  
  c42.SaveAs("BestCos2D_goodVsBad_ele2pt.png");
  //
  TCanvas c43("c43","",1);
  h_cos2DM_notok_kpt->DrawNormalized();
  h_cos2DM_ok_kpt->DrawNormalized("same");
  leg2.Draw("same");  
  c43.SaveAs("BestCos2D_goodVsBad_kpt.png");
  //
  TCanvas c44("c44","",1);
  h_cos2DM_ok_ele1eta->DrawNormalized();
  h_cos2DM_notok_ele1eta->DrawNormalized("same");
  leg2.Draw("same");  
  c44.SaveAs("BestCos2D_goodVsBad_ele1eta.png");
  //
  TCanvas c45("c45","",1);
  h_cos2DM_ok_ele2eta->DrawNormalized();
  h_cos2DM_notok_ele2eta->DrawNormalized("same");
  leg2.Draw("same");  
  c45.SaveAs("BestCos2D_goodVsBad_ele2eta.png");
  //
  TCanvas c46("c46","",1);
  h_cos2DM_ok_keta->DrawNormalized();
  h_cos2DM_notok_keta->DrawNormalized("same");
  leg2.Draw("same");  
  c46.SaveAs("BestCos2D_goodVsBad_keta.png");
  //
  TCanvas c41b("c41b","",1);
  hz_cos2DM_ok_ele1pt->DrawNormalized();
  hz_cos2DM_notok_ele1pt->DrawNormalized("same");
  leg2.Draw("same");  
  c41b.SaveAs("BestCos2D_goodVsBad_ele1ptZ.png");
  //
  TCanvas c42b("c42b","",1);
  hz_cos2DM_notok_ele2pt->DrawNormalized();
  hz_cos2DM_ok_ele2pt->DrawNormalized("same");
  leg2.Draw("same");  
  c42b.SaveAs("BestCos2D_goodVsBad_ele2ptZ.png");
  //
  TCanvas c43b("c43b","",1);
  hz_cos2DM_notok_kpt->DrawNormalized();
  hz_cos2DM_ok_kpt->DrawNormalized("same");
  leg2.Draw("same");  
  c43b.SaveAs("BestCos2D_goodVsBad_kptZ.png");
  //
}
