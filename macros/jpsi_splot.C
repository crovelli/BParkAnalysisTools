// Based on ROOT/tutorials/roostats/rs301_splot.C

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooExponential.h"
#include "RooAddPdf.h"
#include "RooAddition.h"
#include "RooAbsPdf.h"
#include "RooFit.h"
#include "RooFitResult.h"
#include "RooWorkspace.h"
#include "RooDoubleCB.h"
#include "RooStats/SPlot.h"

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TH1.h"

#include <iostream>

// use this order for safety on library loading
using namespace RooFit;
using namespace RooStats;
using namespace std;

// see below for implementation
void AddModel(RooWorkspace*, float, float);
void DoSPlot(RooWorkspace*);
void MakePlots(RooWorkspace*);
void MakeHistos(RooWorkspace*);
void getDataSet(const char *, RooWorkspace*, float, float);

void jpsi_splot()
{
  // set range of observable
  Float_t lowRange  = 2.; 
  Float_t highRange = 4.;  

  // Create a new workspace to manage the project.
  RooWorkspace* wspace = new RooWorkspace("myWS");
  
  // add the signal and background models to the workspace.
  // Inside this function you will find a description our model.
  AddModel(wspace, lowRange, highRange);
  
  // add dataset from converted root tree
  // getDataSet("/afs/cern.ch/work/c/crovelli/BPhys/BParkAnalysisTools/macros/Formatted_BuToKJpsi_Toee_BParkNANO_mc_2020May16_ext_probeLowPt.root", wspace, lowRange, highRange);
  getDataSet("Formatted_ParkingBPH1_Run2018D_ALL_probeLowPt.root", wspace, lowRange, highRange);
  
  // inspect the workspace if you wish
  wspace->Print();
  
  // make a new dataset with sWeights added for every event.
  DoSPlot(wspace);

  // Make some plots showing the discriminating variable and
  // the control variable after unfolding.
  MakePlots(wspace);
  
  // Save variables in histos
  MakeHistos(wspace);

  // cleanup
  delete wspace;

}

// Signal and background fit models
void AddModel(RooWorkspace* ws, float lowRange, float highRange){
  
  // make a RooRealVar for the observables
  RooRealVar pair_mass("pair_mass", "M_{inv}", lowRange, highRange,"GeV");
  // 
  RooRealVar probeMvaId("probeMvaId", "probeMvaId", -10., 10., "");
  RooRealVar probePt("probePt", "probePt", 0., 15., "");
  RooRealVar probeEta("probeEta", "probeEta", -2.4, 2.4, "");
  RooRealVar probeFBrem("probeFBrem", "probeFBrem", 0., 1., "");
  RooRealVar probeDxySig("probeDxySig", "probeDxySig", -50., 50., "");
  RooRealVar probeDzSig("probeDzSig", "probeDzSig", -50., 50., "");

  // --------------------------------------
  // signal model
  std::cout << "make JPsi model" << std::endl;
  //
  RooRealVar m0("m0", "JPsi Mass", 3.0969, 3.09, 3.10);                 // Init MC: 3.09, 3.07, 3.11;    Init data: 3.0969, 3.09, 3.10
  RooRealVar sigma("sigma", "sigma",  0.05, 0.03, 0.1);                 // Init MC: 0.05, 0.01, 0.1;     Init data: 0.05, 0.03, 0.1
  RooRealVar alphaL("alphaL", "alpha left",  0.6, 0.5, 0.7);            // Init MC: 0.5,  0.1,  0.9;     Init data: 0.6, 0.5, 0.7
  RooRealVar alphaR("alphaR", "alpha right", 1.2, 1.1, 1.3);            // Init MC: 1.2,  0.8,  1.6;     Init data: 1.2, 1.1, 1.3  
  RooRealVar nL("nL", "N left",  3.6, 3.56, 3.64);                      // Init MC: 3.5, 2, 5;           Init data: 3.6, 3.56, 3.64
  RooRealVar nR("nR", "N right", 1.85, 1.80, 1.95);                     // Init MC: 2, 1, 3;             Init data: 1.85, 1.80, 1.95  
  RooDoubleCB jpsiModel("jpsiModel", "JPsi Model", pair_mass, m0, sigma, alphaL, alphaR, nL, nR);

  // Successful fit to MC only 
  // 1  alpha        9.99999e-01   1.24767e-01   5.00000e-01  -4.71072e+00
  // 2  alphaL       5.06119e-01   1.40426e-02   7.66633e-04   1.52980e-02
  // 3  alphaR       1.13909e+00   4.14844e-02   1.98050e-03  -1.52874e-01
  // 4  bkgYield     4.18975e+02   1.60700e+02   3.37642e-04  -1.39173e+00
  // 5  jpsiYield    5.01047e+04   2.74905e+02   2.00192e-04   7.33677e-01
  // 6  m0           3.09223e+00   6.70656e-04   2.13518e-03   1.11760e-01
  // 7  nL           3.64837e+00   1.17983e-01   6.26960e-04   9.90746e-02
  // 8  nR           2.23064e+00   1.54501e-01   3.03268e-03   2.32739e-01
  // 9  sigma        4.83608e-02   1.13788e-03   6.75603e-04  -1.48078e-01
  //
  //
  // Successful fit to data
  //1  alpha        4.53409e-01   4.62169e-03   1.40754e-02  -1.55934e-01
  //2  alphaL       5.95870e-01   2.54297e-02   1.89544e-01  -4.13119e-02
  //3  alphaR       1.10001e+00   1.31155e-01   5.00000e-01  -1.58137e+00
  //4  bkgYield     1.70276e+05   5.59852e+02   1.75545e-03  -3.27354e-01
  //5  jpsiYield    1.69462e+04   4.01390e+02   2.68122e-03  -1.21134e+00
  //6  m0           3.09558e+00   5.04847e-04   2.53231e-01  -3.25763e+00
  //7  nL           3.56000e+00   6.22620e-03   5.00000e-01  -1.56681e+00
  //8  nR           1.84336e+00   3.70475e-02   5.00000e-01  -4.35559e-01
  //9  sigma        4.84558e-02   1.88821e-03   3.44988e-02   3.63393e+00

  // we know JPsi mass
  // m0.setConstant();
  
  // --------------------------------------
  // background model
  std::cout << "make background model" << std::endl;
  RooRealVar alpha("alpha", "Decay const for background mass spectrum", 0.5, 0.2, 0.8,"1/GeV");    // Init MC: 0.6, 0, 1;  Init data: 0.5, 0.2, 0.8
  RooExponential bkgModel("bkgModel", "bkg Mass Model", pair_mass, alpha);
  
  // --------------------------------------
  // combined model
  
  // These variables represent the number of JPsi or background events to be fitted
  RooRealVar jpsiYield("jpsiYield","fitted yield for JPsi",     500 ,1.,500000) ;      // Init MC: 500 ,100.,60000; Init Data: 500 ,1000.,500000
  //RooRealVar jpsiYield("jpsiYield","fitted yield for JPsi",     500 ,1000.,500000) ;      // Init MC: 500 ,100.,60000; Init Data: 500 ,1000.,500000
  RooRealVar bkgYield("bkgYield","fitted yield for background", 500 ,1000.,500000) ;      // Init MC: 500 ,100.,40000; Init Data: 500 ,1000.,500000
  
  // now make the combined model
  std::cout << "make full model" << std::endl; 
  RooAddPdf model("model","jpsi+background models",
		  RooArgList(jpsiModel, bkgModel),
		  RooArgList(jpsiYield,bkgYield)); 

  std::cout << "import model" << std::endl;
  ws->import(model);
}

// Add s-weights to the dataset
void DoSPlot(RooWorkspace* ws){
  
  std::cout << "Calculate sWeights" << std::endl;
  
  // get what we need out of the workspace to do the fit
  RooAbsPdf* model = ws->pdf("model");
  RooRealVar* jpsiYield = ws->var("jpsiYield");
  RooRealVar* bkgYield = ws->var("bkgYield");
  RooDataSet* data = (RooDataSet*) ws->data("data");

  // fit the model to the data.
  model->fitTo(*data, Extended() );

  //TCanvas* cdata = new TCanvas("cdata","data fit", 1);
  //RooRealVar* pair_mass = ws->var("pair_mass");
  //RooPlot* frame = pair_mass->frame() ;
  //data->plotOn(frame) ;
  //model->plotOn(frame) ;
  //model->plotOn(frame,Components(*jpsiModel),LineStyle(kDashed), LineColor(kRed)) ;
  //model->plotOn(frame4,Components(*bkgModel),LineStyle(kDashed),LineColor(kGreen)) ;
  //frame->SetTitle("Fit of model to discriminating variable");
  //frame->Draw() ;
  //cdata->SaveAs("Fit.png");

  // The sPlot technique requires that we fix the parameters
  // of the model that are not yields after doing the fit.
  // This *could* be done with the lines below, however this is taken care of
  // by the RooStats::SPlot constructor (or more precisely the AddSWeight method).
  RooRealVar* m0  = ws->var("m0");
  RooRealVar* sigma = ws->var("sigma");
  RooRealVar* alphaL = ws->var("alphaL");
  RooRealVar* alphaR = ws->var("alphaR");
  RooRealVar* nL = ws->var("nL");
  RooRealVar* nR = ws->var("nR");
  RooRealVar* alpha  = ws->var("alpha");
  m0->setConstant();   
  sigma->setConstant();   
  alphaL->setConstant();   
  alphaR->setConstant();   
  nL->setConstant();   
  nR->setConstant();   
  alpha->setConstant();   

  RooMsgService::instance().setSilentMode(false);

  // Now we use the SPlot class to add SWeights to our data set
  // based on our model and our yield variables
  RooStats::SPlot* sData = new RooStats::SPlot("sData","An SPlot",
					       *data, model, RooArgList(*jpsiYield,*bkgYield) );


  // Check that our weights have the desired properties
  std::cout << "Check SWeights:" << std::endl;

  std::cout << std::endl <<  "Yield of JPsi is "
	    << jpsiYield->getVal() << ".  From sWeights it is "
	    << sData->GetYieldFromSWeight("jpsiYield") << std::endl;
  
  std::cout << "Yield of background is "
	    << bkgYield->getVal() << ".  From sWeights it is "
	    << sData->GetYieldFromSWeight("bkgYield") << std::endl
	    << std::endl;
  
  for(Int_t i=0; i < 10; i++) {
      std::cout << "jpsi Weight = " << sData->GetSWeight(i,"jpsiYield")
		<< ", bkg Weight = " << sData->GetSWeight(i,"bkgYield")
		<< ", Total Weight = " << sData->GetSumOfEventSWeight(i)
		<< std::endl;
    }
  
  std::cout << std::endl;

  // import this new dataset with sWeights
  std::cout << "import new dataset with sWeights" << std::endl;
  ws->import(*data, Rename("dataWithSWeights"));
}

// Control plots
void MakePlots(RooWorkspace* ws){
  
  // Here we make plots of the discriminating variable (pair_mass) after the fit
  // and of the control variable (ID, or whatever) after unfolding with sPlot.
  std::cout << std::endl;
  std::cout << "make plots" << std::endl;

  // make our canvas
  TCanvas* cdata = new TCanvas("sPlot","sPlot demo", 400, 600);
  cdata->Divide(1,3);

  // get what we need out of the workspace
  RooAbsPdf* model = ws->pdf("model");
  RooAbsPdf* jpsiModel = ws->pdf("jpsiModel");
  RooAbsPdf* bkgModel = ws->pdf("bkgModel");
  
  RooRealVar* probeMvaId = ws->var("probeMvaId");
  RooRealVar* pair_mass = ws->var("pair_mass");

  // note, we get the dataset with sWeights
  RooDataSet* data = (RooDataSet*) ws->data("dataWithSWeights");
  
  // this shouldn't be necessary, need to fix something with workspace
  // do this to set parameters back to their fitted values.
  model->fitTo(*data, Extended() );

  //plot pair_mass for data with full model and individual components overlaid
  cdata->cd(1);
  RooPlot* frame = pair_mass->frame() ;
  data->plotOn(frame ) ;
  model->plotOn(frame) ;
  model->plotOn(frame,Components(*jpsiModel),LineStyle(kDashed), LineColor(kRed)) ;
  model->plotOn(frame,Components(*bkgModel),LineStyle(kDashed),LineColor(kGreen)) ;
  frame->SetTitle("Fit of model to discriminating variable");
  frame->Draw() ;


  // Now use the sWeights to show our variable distribution for JPsi and background.
  //
  // Plot our variable for JPsi component.
  // Do this by plotting all events weighted by the sWeight for the JPsi component.
  // The SPlot class adds a new variable that has the name of the corresponding
  // yield + "_sw".
   cdata->cd(2);
   
  // create weighted data set
  RooDataSet * dataw_jpsi = new RooDataSet(data->GetName(),data->GetTitle(),data,*data->get(),0,"jpsiYield_sw") ;
  
  RooPlot* frame2 = probeMvaId->frame() ;
  dataw_jpsi->plotOn(frame2, DataError(RooAbsData::SumW2) ) ;
  frame2->SetTitle("ID distribution for JPsi");
  frame2->Draw() ;

  // Plot interesting variables for background
  cdata->cd(3);
  RooDataSet * dataw_bkg = new RooDataSet(data->GetName(),data->GetTitle(),data,*data->get(),0,"bkgYield_sw") ;
  RooPlot* frame3 = probeMvaId->frame() ;
  dataw_bkg->plotOn(frame3,DataError(RooAbsData::SumW2) ) ;
  frame3->SetTitle("ID distribution for background");
  frame3->Draw() ;
  
  cdata->SaveAs("SPlot.png");


  // Fit variable
  TCanvas* cdata2 = new TCanvas("cdata2","data fit", 1);
  RooPlot* frame4 = pair_mass->frame() ;
  data->plotOn(frame4 ) ;
  model->plotOn(frame4) ;
  model->plotOn(frame4,Components(*jpsiModel),LineStyle(kDashed), LineColor(kRed)) ;
  model->plotOn(frame4,Components(*bkgModel),LineStyle(kDashed),LineColor(kGreen)) ;
  frame4->SetTitle("Fit of model to discriminating variable");
  frame4->Draw() ;
  cdata2->SaveAs("Fit.png");
}

// Control plots
void MakeHistos(RooWorkspace* ws){
  
  gStyle->SetOptStat(0);

  std::cout << std::endl;
  std::cout << "save histos" << std::endl;

  RooRealVar* probeMvaId = ws->var("probeMvaId");
  RooRealVar* probePt  = ws->var("probePt");
  RooRealVar* probeEta = ws->var("probeEta");
  RooRealVar* probeFBrem = ws->var("probeFBrem");
  RooRealVar* probeDxySig = ws->var("probeDxySig");
  RooRealVar* probeDzSig = ws->var("probeDzSig");
  
  // note, we get the dataset with sWeights
  RooDataSet* data = (RooDataSet*) ws->data("dataWithSWeights");
  
  // create weighted data set
  RooDataSet * dataw_jpsi = new RooDataSet(data->GetName(),data->GetTitle(),data,*data->get(),0,"jpsiYield_sw") ;
  RooDataSet * dataw_bkg  = new RooDataSet(data->GetName(),data->GetTitle(),data,*data->get(),0,"bkgYield_sw") ;
  
  // convert to TH1
  TH1 *h1_probeMvaId_jpsi  = dataw_jpsi->createHistogram("h1_probeMvaId_jpsi",*probeMvaId,Binning(60));
  TH1 *h1_probeMvaId_bkg   = dataw_bkg->createHistogram("h1_probeMvaId_bkg",*probeMvaId,Binning(60));
  TH1 *h1_probePt_jpsi     = dataw_jpsi->createHistogram("h1_probePt_jpsi",*probePt,Binning(60));
  TH1 *h1_probePt_bkg      = dataw_bkg->createHistogram("h1_probePt_bkg",*probePt,Binning(60));
  TH1 *h1_probeEta_jpsi    = dataw_jpsi->createHistogram("h1_probeEta_jpsi",*probeEta,Binning(40));
  TH1 *h1_probeEta_bkg     = dataw_bkg->createHistogram("h1_probeEta_bkg",*probeEta,Binning(40));
  TH1 *h1_probeFBrem_jpsi  = dataw_jpsi->createHistogram("h1_probeFBrem_jpsi",*probeFBrem,Binning(50));
  TH1 *h1_probeFBrem_bkg   = dataw_bkg->createHistogram("h1_probeFBrem_bkg",*probeFBrem,Binning(50));
  TH1 *h1_probeDxySig_jpsi = dataw_jpsi->createHistogram("h1_probeDxySig_jpsi",*probeDxySig,Binning(50));
  TH1 *h1_probeDxySig_bkg  = dataw_bkg->createHistogram("h1_probeDxySig_bkg",*probeDxySig,Binning(50));
  TH1 *h1_probeDzSig_jpsi  = dataw_jpsi->createHistogram("h1_probeDzSig_jpsi",*probeDzSig,Binning(50));
  TH1 *h1_probeDzSig_bkg   = dataw_bkg->createHistogram("h1_probeDzSig_bkg",*probeDzSig,Binning(50));

  TFile myFileSPlots("myFileSPlots.root","RECREATE");
  myFileSPlots.cd();
  h1_probeMvaId_jpsi->Write();
  h1_probeMvaId_bkg->Write();
  h1_probePt_jpsi->Write();
  h1_probePt_bkg->Write();
  h1_probeEta_jpsi->Write();
  h1_probeEta_bkg->Write();
  h1_probeFBrem_jpsi->Write();
  h1_probeFBrem_bkg->Write();
  h1_probeDxySig_jpsi->Write();
  h1_probeDxySig_bkg->Write();
  h1_probeDzSig_jpsi->Write();
  h1_probeDzSig_bkg->Write();

  TCanvas* ch1 = new TCanvas("ch1","ch1", 1);
  h1_probeMvaId_jpsi->SetLineWidth(2);
  h1_probeMvaId_bkg->SetLineWidth(2);
  h1_probeMvaId_jpsi->SetLineColor(6);
  h1_probeMvaId_bkg->SetLineColor(4);
  h1_probeMvaId_jpsi->SetTitle("");
  h1_probeMvaId_bkg->SetTitle("");
  h1_probeMvaId_jpsi->DrawNormalized("hist");
  h1_probeMvaId_bkg->DrawNormalized("samehist");
  ch1->SaveAs("probeMvaIdH.png");

  TCanvas* ch2 = new TCanvas("ch2","ch2", 1);
  h1_probePt_jpsi->SetLineWidth(2);
  h1_probePt_bkg->SetLineWidth(2);
  h1_probePt_jpsi->SetLineColor(6);
  h1_probePt_bkg->SetLineColor(4);
  h1_probePt_jpsi->SetTitle("");
  h1_probePt_bkg->SetTitle("");
  h1_probePt_bkg->DrawNormalized("hist");
  h1_probePt_jpsi->DrawNormalized("samehist");
  ch2->SaveAs("probePt.png");

  TCanvas* ch3 = new TCanvas("ch3","ch3", 1);
  h1_probeEta_jpsi->SetLineWidth(2);
  h1_probeEta_bkg->SetLineWidth(2);
  h1_probeEta_jpsi->SetLineColor(6);
  h1_probeEta_bkg->SetLineColor(4);
  h1_probeEta_jpsi->SetTitle("");
  h1_probeEta_bkg->SetTitle("");
  h1_probeEta_jpsi->DrawNormalized("hist");
  h1_probeEta_bkg->DrawNormalized("samehist");
  ch3->SaveAs("probeEta.png");

  TCanvas* ch4 = new TCanvas("ch4","ch4", 1);
  h1_probeFBrem_jpsi->SetLineWidth(2);
  h1_probeFBrem_bkg->SetLineWidth(2);
  h1_probeFBrem_jpsi->SetLineColor(6);
  h1_probeFBrem_bkg->SetLineColor(4);
  h1_probeFBrem_jpsi->SetTitle("");
  h1_probeFBrem_bkg->SetTitle("");
  h1_probeFBrem_bkg->DrawNormalized("hist");
  h1_probeFBrem_jpsi->DrawNormalized("samehist");
  ch4->SaveAs("probeFBrem.png");

  TCanvas* ch5 = new TCanvas("ch5","ch5", 1);
  h1_probeDxySig_jpsi->SetLineWidth(2);
  h1_probeDxySig_bkg->SetLineWidth(2);
  h1_probeDxySig_jpsi->SetLineColor(6);
  h1_probeDxySig_bkg->SetLineColor(4);
  h1_probeDxySig_jpsi->SetTitle("");
  h1_probeDxySig_bkg->SetTitle("");
  h1_probeDxySig_bkg->DrawNormalized("hist");
  h1_probeDxySig_jpsi->DrawNormalized("samehist");
  ch5->SaveAs("probeDxySig.png");

  TCanvas* ch6 = new TCanvas("ch6","ch6", 1);
  h1_probeDzSig_jpsi->SetLineWidth(2);
  h1_probeDzSig_bkg->SetLineWidth(2);
  h1_probeDzSig_jpsi->SetLineColor(6);
  h1_probeDzSig_bkg->SetLineColor(4);
  h1_probeDzSig_jpsi->SetTitle("");
  h1_probeDzSig_bkg->SetTitle("");
  h1_probeDzSig_bkg->DrawNormalized("hist");
  h1_probeDzSig_jpsi->DrawNormalized("samehist");
  ch6->SaveAs("probeDzSig.png");
}

// Convert ROOT tree in RooDataset
void getDataSet(const char *rootfile, RooWorkspace *ws, float lowRange, float highRange) {    
  
  cout << "roofitting file " << rootfile << endl;
  
  // fit variables
  RooRealVar pair_mass("pair_mass", "M_{inv}", lowRange, highRange,"GeV");
  // discriminating variables
  RooRealVar probeMvaId("probeMvaId", "probeMvaId", -10., 10., "");
  RooRealVar probePt("probePt", "probePt", 0., 20., "");
  RooRealVar probeEta("probeEta", "probeEta", -2.5, 2.5, "");
  RooRealVar probeFBrem("probeFBrem", "probeFBrem", 0., 1., "");
  RooRealVar probeDxySig("probeDxySig", "probeDxySig", -50., 50., "");
  RooRealVar probeDzSig("probeDzSig", "probeDzSig", -50., 50., "");

  RooArgSet setall(pair_mass,probeMvaId,probePt,probeEta,probeFBrem,probeDxySig,probeDzSig);

  TFile *file = TFile::Open(rootfile);
  TTree *tree = (TTree*)file->Get("tnpAna/fitter_tree");

  RooDataSet *data = new RooDataSet("data","data",tree,setall,0); 

  // Barrel:
  // pt: 0.5-1.0 GeV 
  // data = (RooDataSet*)data->reduce("probePt>0.5 && probePt<1.0 && probeEta<1.5 && probeEta>-1.5");    
  // pt: 1.0-1.5 GeV 
  // data = (RooDataSet*)data->reduce("probePt>1.0 && probePt<1.5 && probeEta<1.5 && probeEta>-1.5");    
  // pt: 1.5-2.0 GeV 
  // data = (RooDataSet*)data->reduce("probePt>1.5 && probePt<2.0 && probeEta<1.5 && probeEta>-1.5");    
  // pt: 2.0-5.0 GeV 
  // data = (RooDataSet*)data->reduce("probePt>2.0 && probePt<5.0 && probeEta<1.5 && probeEta>-1.5");    
  // pt: >5.0 GeV 
  // data = (RooDataSet*)data->reduce("probePt>5.0 && probeEta<1.5 && probeEta>-1.5");    
  //
  // Endcap:
  // pt: 0.5-1.0 GeV 
  // data = (RooDataSet*)data->reduce("probePt>0.5 && probePt<1.0 && (probeEta<-1.5 || probeEta>1.5)");    
  // pt: 1.0-1.5 GeV 
  // data = (RooDataSet*)data->reduce("probePt>1.0 && probePt<1.5 && (probeEta<-1.5 || probeEta>1.5)");    
  // pt: 1.5-2.0 GeV 
  // data = (RooDataSet*)data->reduce("probePt>1.5 && probePt<2.0 && (probeEta<-1.5 || probeEta>1.5)");    
  // pt: 2.0-5.0 GeV 
  // data = (RooDataSet*)data->reduce("probePt>2.0 && probePt<5.0 && (probeEta<-1.5 || probeEta>1.5)");    
  // pt: >5.0 GeV 
  data = (RooDataSet*)data->reduce("probePt>5.0 && (probeEta<-1.5 || probeEta>1.5)");    

  data->Print();

  ws->import(*data);
}
