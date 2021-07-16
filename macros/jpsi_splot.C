// Based on ROOT/tutorials/roostats/rs301_splot.C

// To be run on data tnp formatted ntuples 
// - Fit mee invariant mass as a reference distribution
// - Extract electron-related distributions for signal and background with the sPlots technique

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
  Float_t lowRange  = 2.6;   
  Float_t highRange = 3.4;   

  // Create a new workspace to manage the project.
  RooWorkspace* wspace = new RooWorkspace("myWS");
  
  // add the signal and background models to the workspace.
  // Inside this function you will find a description our model.
  AddModel(wspace, lowRange, highRange);
  
  // add dataset from converted root tree
  getDataSet("/eos/cms/store/user/crovelli/LowPtEle/TnpData/March21noRegression/Formatted_March21_NoRegr_ParkingALL_tightCuts.root", wspace, lowRange, highRange);
  
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

  // --------------------------------------
  // signal model
  std::cout << "make JPsi model" << std::endl;


  // ------------------------------------------------
  //
  // Inclusive allEB+EE
  // loose selection
  /*
  RooRealVar m0("m0", "JPsi Mass", 3.089, 3.0885, 3.0895);                  
  RooRealVar alphaL("alphaL", "alpha left",  0.5, 0.3, 0.6);          
  RooRealVar alphaR("alphaR", "alpha right", 1., 0.999, 1.001);         
  RooRealVar sigma("sigma", "sigma",  0.04, 0.03, 0.05);  
  RooRealVar nL("nL", "N left",  3.50, 3.48, 3.52);        
  RooRealVar nR("nR", "N right", 1.60, 1.58, 1.62);       
  RooRealVar alpha("alpha", "Decay const for background mass spectrum", 0.0, -0.2, 0.2,"1/GeV");      
  RooRealVar jpsiYield("jpsiYield","fitted yield for JPsi",      2000 , 1., 500000) ;              
  RooRealVar bkgYield("bkgYield","fitted yield for background", 10000 , 1., 5000000) ;             
  */
  //
  // tight cuts
  RooRealVar m0("m0", "JPsi Mass", 3.0890, 3.0889, 3.0891);                  
  RooRealVar alphaL("alphaL", "alpha left",  0.3, 0.1, 0.6);          
  RooRealVar alphaR("alphaR", "alpha right", 1., 0.995, 1.005);         
  RooRealVar sigma("sigma", "sigma",  0.04, 0.03, 0.05);  
  RooRealVar nL("nL", "N left",  3.40, 2.8, 3.52);        
  RooRealVar nR("nR", "N right", 1.60, 0.8, 1.62);       
  RooRealVar alpha("alpha", "Decay const for background mass spectrum", 0.1, 0., 0.2,"1/GeV");      
  RooRealVar jpsiYield("jpsiYield","fitted yield for JPsi",      2000 , 1., 500000) ;              
  RooRealVar bkgYield("bkgYield","fitted yield for background", 10000 , 1., 5000000) ;             



  // ------------------------------------------------
  // 
  // EB: 1-1.5 
  // loose selection, not the best fit....
  /*  
  RooRealVar m0("m0", "JPsi Mass", 3.089, 3.0885, 3.0895);                  
  RooRealVar alphaL("alphaL", "alpha left",  0.3, 0.2, 0.5);          
  RooRealVar alphaR("alphaR", "alpha right", 1.0, 0.99, 1.01);         
  RooRealVar sigma("sigma", "sigma",  0.04, 0.03, 0.06);  
  RooRealVar nL("nL", "N left",  3.3, 3.2, 3.60);  
  RooRealVar nR("nR", "N right", 1.60, 1.50, 1.70);
  RooRealVar alpha("alpha", "Decay const for background mass spectrum", 0.3, 0., 0.6,"1/GeV");      
  RooRealVar jpsiYield("jpsiYield","fitted yield for JPsi",      2000 , 1., 500000) ;              
  RooRealVar bkgYield("bkgYield","fitted yield for background", 10000 , 1., 5000000) ;             
  m0.setConstant();    
  */
  // tight cuts
  /*
  RooRealVar m0("m0", "JPsi Mass", 3.0895, 3.0895, 3.0895);                  
  RooRealVar alphaL("alphaL", "alpha left",  0.3, 0.2, 0.5);          
  RooRealVar alphaR("alphaR", "alpha right", 1.0, 0.99, 1.01);         
  RooRealVar sigma("sigma", "sigma",  0.04, 0.03, 0.06);  
  RooRealVar nL("nL", "N left",  3.3, 3.2, 3.60);  
  RooRealVar nR("nR", "N right", 1.60, 1.50, 1.70);
  RooRealVar alpha("alpha", "Decay const for background mass spectrum", 0.3, 0., 0.6,"1/GeV");      
  RooRealVar jpsiYield("jpsiYield","fitted yield for JPsi",      2000 , 1., 500000) ;              
  RooRealVar bkgYield("bkgYield","fitted yield for background", 10000 , 1., 5000000) ;             
  m0.setConstant();    
  */


  // ------------------------------------------------
  //
  // EB: 1.5-2 
  // loose selection
  /*
  RooRealVar m0("m0", "JPsi Mass", 3.089, 3.0885, 3.0895);                  
  RooRealVar alphaL("alphaL", "alpha left",  0.5, 0.3, 0.6);          
  RooRealVar alphaR("alphaR", "alpha right", 1., 0.999, 1.001);         
  RooRealVar sigma("sigma", "sigma",  0.04, 0.03, 0.05);  
  RooRealVar nL("nL", "N left",  3.50, 3.48, 3.52);        
  RooRealVar nR("nR", "N right", 1.60, 1.58, 1.62);       
  RooRealVar alpha("alpha", "Decay const for background mass spectrum", 0.0, -0.2, 0.2,"1/GeV");      
  RooRealVar jpsiYield("jpsiYield","fitted yield for JPsi",      2000 , 1., 500000) ;              
  RooRealVar bkgYield("bkgYield","fitted yield for background", 10000 , 1., 5000000) ;             
  */
  //
  // tight cuts
  /*
  RooRealVar m0("m0", "JPsi Mass", 3.0885, 3.0885, 3.0885);                  
  RooRealVar alphaL("alphaL", "alpha left",  0.3, 0.1, 0.6);          
  RooRealVar alphaR("alphaR", "alpha right", 1., 0.999, 1.001);         
  RooRealVar sigma("sigma", "sigma",  0.04, 0.03, 0.05);  
  RooRealVar nL("nL", "N left",  3.40, 3.1, 3.52);        
  RooRealVar nR("nR", "N right", 1.60, 1.4, 1.62);       
  RooRealVar alpha("alpha", "Decay const for background mass spectrum", 0.0, -0.2, 0.2,"1/GeV");      
  RooRealVar jpsiYield("jpsiYield","fitted yield for JPsi",      2000 , 1., 500000) ;              
  RooRealVar bkgYield("bkgYield","fitted yield for background", 10000 , 1., 5000000) ;             
  m0.setConstant();    
  */


  // ----------------------------------------
  //
  // EB: 2-5   
  // loose selection
  /*
  RooRealVar m0("m0", "JPsi Mass", 3.090, 3.090, 3.090);  
  RooRealVar alphaL("alphaL", "alpha left",  0.5, 0.3, 0.6);          
  RooRealVar alphaR("alphaR", "alpha right", 1.3, 1.1, 1.5); 
  RooRealVar sigma("sigma", "sigma",  0.04, 0.03, 0.05);  
  RooRealVar nL("nL", "N left",  3.6, 3.56, 3.64);        
  RooRealVar nR("nR", "N right", 1.85, 1.75, 1.95);       
  RooRealVar alpha("alpha", "Decay const for background mass spectrum", 0.0, -0.1, 0.1,"1/GeV");      
  RooRealVar jpsiYield("jpsiYield","fitted yield for JPsi",      2000 , 1., 500000) ;              
  RooRealVar bkgYield("bkgYield","fitted yield for background", 10000 , 1., 5000000) ;             
  m0.setConstant();  
  */
  //
  //tight selection
  /*
  RooRealVar m0("m0", "JPsi Mass", 3.089, 3.0885, 3.0895);                  
  RooRealVar alphaL("alphaL", "alpha left",  0.5, 0.3, 0.7);          
  RooRealVar alphaR("alphaR", "alpha right", 1.3, 1.1, 1.5); 
  RooRealVar sigma("sigma", "sigma",  0.04, 0.03, 0.05);  
  RooRealVar nL("nL", "N left",  3.6, 3.40, 3.64);        
  RooRealVar nR("nR", "N right", 1.85, 1.75, 1.95);       
  RooRealVar alpha("alpha", "Decay const for background mass spectrum", 0.0, -0.15, 0.15,"1/GeV");      
  RooRealVar jpsiYield("jpsiYield","fitted yield for JPsi",      2000 , 1., 500000) ;              
  RooRealVar bkgYield("bkgYield","fitted yield for background", 10000 , 1., 5000000) ;             
  */

  // ------------------------------------------------
  //
  // EB: >5   
  /*
  // loose cuts 
  RooRealVar m0("m0", "JPsi Mass", 3.091, 3.091, 3.091);                  
  RooRealVar alphaL("alphaL", "alpha left",  0.6, 0.5, 0.8);
  RooRealVar alphaR("alphaR", "alpha right", 1.4, 1.2, 1.6);  
  RooRealVar sigma("sigma", "sigma",  0.05, 0.045, 0.055);  
  RooRealVar nL("nL", "N left",  3.6, 3.56, 3.64);        
  RooRealVar nR("nR", "N right", 1.85, 1.75, 1.95);       
  RooRealVar alpha("alpha", "Decay const for background mass spectrum",  -0.1, -0.5, 0.,"1/GeV"); 
  RooRealVar jpsiYield("jpsiYield","fitted yield for JPsi",      2000 , 1., 500000) ;              
  RooRealVar bkgYield("bkgYield","fitted yield for background", 10000 , 1., 5000000) ;             
  */
  //
  /*
  // tight cuts 
  RooRealVar m0("m0", "JPsi Mass", 3.089, 3.0885, 3.0895);                  
  RooRealVar alphaL("alphaL", "alpha left",  0.6, 0.5, 0.8);
  RooRealVar alphaR("alphaR", "alpha right", 1.4, 1.2, 1.6);  
  RooRealVar sigma("sigma", "sigma",  0.05, 0.045, 0.055);  
  RooRealVar nL("nL", "N left",  3.6, 3.5, 3.64);        
  RooRealVar nR("nR", "N right", 1.85, 1.7, 1.95);       
  RooRealVar alpha("alpha", "Decay const for background mass spectrum",  -0.1, -0.8, 0.,"1/GeV"); 
  RooRealVar jpsiYield("jpsiYield","fitted yield for JPsi",      2000 , 1., 500000) ;              
  RooRealVar bkgYield("bkgYield","fitted yield for background", 10000 , 1., 5000000) ;             
  */


  // ----------------------------------------------------
  // 
  //  EE: 2-5 
  // loose selection
  /*
  RooRealVar m0("m0", "JPsi Mass", 3.089, 3.0885, 3.0895);                  
  RooRealVar alphaL("alphaL", "alpha left",  0.3, 0.2, 0.4);          
  RooRealVar alphaR("alphaR", "alpha right", 0.9, 0.75, 1.0);         
  RooRealVar sigma("sigma", "sigma",  0.05, 0.045, 0.06);  
  RooRealVar nL("nL", "N left",  3.6, 3.56, 3.64);        
  RooRealVar nR("nR", "N right", 1.85, 1.75, 1.95);       
  RooRealVar alpha("alpha", "Decay const for background mass spectrum", 0.0, -0.2, 0.8,"1/GeV");      
  RooRealVar jpsiYield("jpsiYield","fitted yield for JPsi",      2000 , 1., 500000) ;              
  RooRealVar bkgYield("bkgYield","fitted yield for background", 10000 , 1., 5000000) ;             
  m0.setConstant();           
  */
  //
  // tight selection
  // not best fit
  /*
  RooRealVar m0("m0", "JPsi Mass", 3.0880, 3.0880, 3.0880);                  
  RooRealVar alphaL("alphaL", "alpha left",  0.3, 0.2, 0.4);          
  RooRealVar alphaR("alphaR", "alpha right", 0.9, 0.3, 1.0);         
  RooRealVar sigma("sigma", "sigma",  0.05, 0.045, 0.06);  
  RooRealVar nL("nL", "N left",  3.6, 0.5, 4.);        
  RooRealVar nR("nR", "N right", 1.85, 0.5, 1.95);       
  RooRealVar alpha("alpha", "Decay const for background mass spectrum", 0.0, -0.2, 0.8,"1/GeV");      
  RooRealVar jpsiYield("jpsiYield","fitted yield for JPsi",      2000 , 1., 500000) ;              
  RooRealVar bkgYield("bkgYield","fitted yield for background", 10000 , 1., 5000000) ;             
  m0.setConstant(); 
  */


  // ----------------------------------------------------
  // 
  //  EE: pt>5      
  // loose selection
  /*
  RooRealVar m0("m0", "JPsi Mass", 3.089, 3.0885, 3.0895);                  
  RooRealVar alphaL("alphaL", "alpha left",  0.5, 0.3, 0.6);          
  RooRealVar alphaR("alphaR", "alpha right", 0.9, 0.75, 1.0);         
  RooRealVar sigma("sigma", "sigma",  0.06, 0.05, 0.07);  
  RooRealVar nL("nL", "N left",  3.6, 3.56, 3.64);        
  RooRealVar nR("nR", "N right", 1.85, 1.75, 1.95);       
  RooRealVar alpha("alpha", "Decay const for background mass spectrum", 0.0, -0.2, 0.8,"1/GeV");      
  RooRealVar jpsiYield("jpsiYield","fitted yield for JPsi",      2000 , 1., 500000) ;              
  RooRealVar bkgYield("bkgYield","fitted yield for background", 10000 , 1., 5000000) ;             
  m0.setConstant();
  */
  //
  // tight selection
  /*
  RooRealVar m0("m0", "JPsi Mass", 3.0880, 3.0870, 3.0895);                  
  RooRealVar alphaL("alphaL", "alpha left",  0.3, 0.2, 0.4);          
  RooRealVar alphaR("alphaR", "alpha right", 0.9, 0.3, 1.0);         
  RooRealVar sigma("sigma", "sigma",  0.05, 0.045, 0.06);  
  RooRealVar nL("nL", "N left",  3.6, 0.5, 4.);        
  RooRealVar nR("nR", "N right", 1.85, 0.5, 1.95);       
  RooRealVar alpha("alpha", "Decay const for background mass spectrum", 0.0, -0.2, 0.8,"1/GeV");      
  RooRealVar jpsiYield("jpsiYield","fitted yield for JPsi",      2000 , 1., 500000) ;              
  RooRealVar bkgYield("bkgYield","fitted yield for background", 10000 , 1., 5000000) ;             
  */


  RooDoubleCB jpsiModel("jpsiModel", "JPsi Model", pair_mass, m0, sigma, alphaL, alphaR, nL, nR);

  
  // --------------------------------------
  // background model
  std::cout << "make background model" << std::endl;

  RooExponential bkgModel("bkgModel", "bkg Mass Model", pair_mass, alpha);
  
  // --------------------------------------
  // combined model
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

  // Plot pair_mass for data with full model and individual components overlaid
  cdata->cd(1);
  RooPlot* frame = pair_mass->frame() ;
  data->plotOn(frame ) ;
  model->plotOn(frame) ;
  model->plotOn(frame,Components(*jpsiModel),LineStyle(kDashed), LineColor(kRed)) ;
  model->plotOn(frame,Components(*bkgModel),LineStyle(kDashed),LineColor(kGreen)) ;
  frame->SetTitle("Fit of model to discriminating variable");
  frame->Draw() ;

  // ------------------------------------------------------------
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
  data->plotOn(frame4, Binning(100) ) ;
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

  RooRealVar* hlt_9ip6 = ws->var("hlt_9ip6");

  RooRealVar* probePt  = ws->var("probePt");
  RooRealVar* probeEta = ws->var("probeEta");

  RooRealVar* probeDxySig = ws->var("probeDxySig");
  RooRealVar* probeDzTrg = ws->var("probeDzTrg"); 
  RooRealVar* probeIso04Rel = ws->var("probeIso04Rel"); 

  // note, we get the dataset with sWeights
  RooDataSet* data = (RooDataSet*) ws->data("dataWithSWeights");

  // create weighted data set
  RooDataSet * dataw_jpsi = new RooDataSet(data->GetName(),data->GetTitle(),data,*data->get(),0,"jpsiYield_sw") ;
  RooDataSet * dataw_bkg  = new RooDataSet(data->GetName(),data->GetTitle(),data,*data->get(),0,"bkgYield_sw") ;
  
  // convert to TH1
  TH1 *h1_probeMvaId_jpsi  = dataw_jpsi->createHistogram("h1_probeMvaId_jpsi",*probeMvaId,Binning(40));         
  TH1 *h1_probeMvaId_bkg   = dataw_bkg->createHistogram("h1_probeMvaId_bkg",*probeMvaId,Binning(40));           

  TH1 *h1_probeDxySig_jpsi = dataw_jpsi->createHistogram("h1_probeDxySig_jpsi",*probeDxySig,Binning(120,-30.,30.)); 
  TH1 *h1_probeDxySig_bkg  = dataw_bkg->createHistogram("h1_probeDxySig_bkg",*probeDxySig,Binning(120,-30.,30.)); 

  TH1 *h1_probeDzTrg_jpsi = dataw_jpsi->createHistogram("h1_probeDzTrg_jpsi",*probeDzTrg,Binning(100,-1.,1.)); 
  TH1 *h1_probeDzTrg_bkg  = dataw_bkg->createHistogram("h1_probeDzTrg_bkg",*probeDzTrg,Binning(100,-1.,1.)); 
  
  TH1 *h1_probeIso04Rel_jpsi = dataw_jpsi->createHistogram("h1_probeIso04Rel_jpsi",*probeIso04Rel,Binning(50,0.,50.)); 
  TH1 *h1_probeIso04Rel_bkg  = dataw_bkg->createHistogram("h1_probeIso04Rel_bkg",*probeIso04Rel,Binning(50,0.,50.)); 

  TH1 *h1_probePt_jpsi = dataw_jpsi->createHistogram("h1_probePt_jpsi",*probePt,Binning(90, 0.,30.));
  TH1 *h1_probePt_bkg  = dataw_bkg->createHistogram("h1_probePt_bkg",*probePt,Binning(90,0.,30.));

  TH1 *h1_probeEta_jpsi = dataw_jpsi->createHistogram("h1_probeEta_jpsi",*probeEta,Binning(40));
  TH1 *h1_probeEta_bkg  = dataw_bkg->createHistogram("h1_probeEta_bkg",*probeEta,Binning(40));

  TFile myFileSPlots("myFileSPlots.root","RECREATE");
  myFileSPlots.cd();

  h1_probeDxySig_jpsi->Write();        
  h1_probeDxySig_bkg->Write();        
  h1_probeDzTrg_jpsi->Write();        
  h1_probeDzTrg_bkg->Write();        
  h1_probeIso04Rel_jpsi->Write();        
  h1_probeIso04Rel_bkg->Write();        

  h1_probeMvaId_jpsi->Write();       
  h1_probeMvaId_bkg->Write();        
  h1_probePt_jpsi->Write();
  h1_probePt_bkg->Write();
  h1_probeEta_jpsi->Write();
  h1_probeEta_bkg->Write();

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
  h1_probePt_jpsi->DrawNormalized("hist");
  h1_probePt_bkg->DrawNormalized("samehist");
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
  h1_probeDxySig_jpsi->SetLineWidth(2);
  h1_probeDxySig_bkg->SetLineWidth(2);
  h1_probeDxySig_jpsi->SetLineColor(6);
  h1_probeDxySig_bkg->SetLineColor(4);
  h1_probeDxySig_jpsi->SetTitle("");
  h1_probeDxySig_bkg->SetTitle("");
  h1_probeDxySig_bkg->DrawNormalized("hist");
  h1_probeDxySig_jpsi->DrawNormalized("samehist");
  ch4->SaveAs("probeDxySig.png");

  TCanvas* ch5 = new TCanvas("ch5","ch5", 1);
  h1_probeDzTrg_jpsi->SetLineWidth(2);
  h1_probeDzTrg_bkg->SetLineWidth(2);
  h1_probeDzTrg_jpsi->SetLineColor(6);
  h1_probeDzTrg_bkg->SetLineColor(4);
  h1_probeDzTrg_jpsi->SetTitle("");
  h1_probeDzTrg_bkg->SetTitle("");
  h1_probeDzTrg_bkg->DrawNormalized("hist");
  h1_probeDzTrg_jpsi->DrawNormalized("samehist");
  ch5->SaveAs("probeDzTrg.png");

  TCanvas* ch6 = new TCanvas("ch6","ch6", 1);
  h1_probeIso04Rel_jpsi->SetLineWidth(2);
  h1_probeIso04Rel_bkg->SetLineWidth(2);
  h1_probeIso04Rel_jpsi->SetLineColor(6);
  h1_probeIso04Rel_bkg->SetLineColor(4);
  h1_probeIso04Rel_jpsi->SetTitle("");
  h1_probeIso04Rel_bkg->SetTitle("");
  h1_probeIso04Rel_bkg->DrawNormalized("hist");
  h1_probeIso04Rel_jpsi->DrawNormalized("samehist");
  ch6->SaveAs("probeDzTrg.png");
}

// Convert ROOT tree in RooDataset
void getDataSet(const char *rootfile, RooWorkspace *ws, float lowRange, float highRange) {    
  
  cout << "roofitting file " << rootfile << endl;
  
  // fit variables
  RooRealVar pair_mass("pair_mass", "M_{inv}", lowRange, highRange,"GeV");
  // 
  // trigger
  RooRealVar hlt_9ip6("hlt_9ip6", "hlt_9ip6", -0.5, 1.5, "");
  // 
  // discriminating variables
  RooRealVar probeMvaId("probeMvaId", "probeMvaId", -4, 16, "");           
  RooRealVar probePt("probePt", "probePt", 0., 1000., "");
  RooRealVar probeEta("probeEta", "probeEta", -2.4, 2.4, "");

  RooRealVar probeDxySig("probeDxySig", "probeDxySig", -200., 200., ""); 
  RooRealVar probeDzTrg("probeDzTrg", "probeDzTrg", -10., 10., ""); 
  RooRealVar probeIso04Rel("probeIso04Rel", "probeIso04Rel", -1., 5000., "");

  RooArgSet setall(pair_mass,hlt_9ip6,probeMvaId,probePt,probeEta,probeDxySig,probeDzTrg,probeIso04Rel);

  TFile *file = TFile::Open(rootfile);
  TTree *tree = (TTree*)file->Get("tnpAna/fitter_tree");

  RooDataSet *data = new RooDataSet("data","data",tree,setall,0); 

  // Inclusive
  data = (RooDataSet*)data->reduce("hlt_9ip6==1 && ( (probePt>1.0 && probeEta>-1.5 && probeEta<1.5) || (probePt>2.0 && (probeEta<-1.5 || probeEta>1.5)) )");   

  // EB / EE
  // data = (RooDataSet*)data->reduce("hlt_9ip6==1 && probePt>1.0 && probeEta>-1.5 && probeEta<1.5");          
  // data = (RooDataSet*)data->reduce("hlt_9ip6==1 && probePt>2.0 && (probeEta<-1.5 || probeEta>1.5)");    

  // Barrel:
  // pt: 1-1.5 GeV    
  // data = (RooDataSet*)data->reduce("hlt_9ip6==1 && probePt>1.0 && probePt<1.5 && probeEta<1.5 && probeEta>-1.5");      
  // pT: 1.5-2 GeV
  // data = (RooDataSet*)data->reduce("hlt_9ip6==1 && probePt>1.5 && probePt<2.0 && probeEta<1.5 && probeEta>-1.5");    
  // pt: 2.0-5.0 GeV 
  // data = (RooDataSet*)data->reduce("hlt_9ip6==1 && probePt>2.0 && probePt<5.0 && probeEta<1.5 && probeEta>-1.5");    
  // pt: >5.0 GeV 
  // data = (RooDataSet*)data->reduce("hlt_9ip6==1 && probePt>5.0 && probeEta<1.5 && probeEta>-1.5");    
  //
  // Endcap:
  // pt: 1.0-2.0 GeV - test
  // data = (RooDataSet*)data->reduce("hlt_9ip6==1 && probePt>1.0 && probePt<2.0 && (probeEta<-1.5 || probeEta>1.5)");    
  // pt: 2.0-5.0 GeV 
  // data = (RooDataSet*)data->reduce("hlt_9ip6==1 && probePt>2.0 && probePt<5.0 && (probeEta<-1.5 || probeEta>1.5)");    
  // pt: >5.0 GeV 
  // data = (RooDataSet*)data->reduce("hlt_9ip6==1 && probePt>5.0 && (probeEta<-1.5 || probeEta>1.5)");            

  data->Print();

  ws->import(*data);
}
