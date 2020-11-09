#define weightsNanoVsTnp_cxx

// ROOT includes 
#include <TROOT.h>
#include <TStyle.h>
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <iostream>

using namespace std;

// To compute weights nani->tnp for signal distributions
// and to make control distributions
// before and after applying weights

void weightsNanoVsTnp()
{
  // Files
  TFile fileTnP("myFileMcAfterTnp.root");          // using TnP selection applied to MC (prepareInputsFromMcWithTnP)
  TFile fileNoTnP("myFileFromNani.root");          // without TnP selection, directy applied to nanpAODs (prepareInputsFromNaniInMc)

  // Histos: pT
  TH1F *ptSignalMcFromTnP   = (TH1F*)fileTnP.Get("probePtSignalMc");
  TH1F *ptSignalMcFromTnPWW = (TH1F*)fileTnP.Get("probePtSignalMcWW");
  TH1F *ptSignalMcNoTnP     = (TH1F*)fileNoTnP.Get("ptSignalMc");
  ptSignalMcFromTnP   -> Sumw2();   
  ptSignalMcFromTnPWW -> Sumw2();   
  ptSignalMcNoTnP     -> Sumw2(); 
  ptSignalMcFromTnP   -> Scale(1./ptSignalMcFromTnP->Integral());   
  ptSignalMcFromTnPWW -> Scale(1./ptSignalMcFromTnPWW->Integral());   
  ptSignalMcNoTnP     -> Scale(1./ptSignalMcNoTnP->Integral());   

  // Histos: eta
  TH1F *etaSignalMcFromTnP   = (TH1F*)fileTnP.Get("probeEtaSignalMc");
  TH1F *etaSignalMcFromTnPWW = (TH1F*)fileTnP.Get("probeEtaSignalMcWW");
  TH1F *etaSignalMcNoTnP     = (TH1F*)fileNoTnP.Get("etaSignalMc");
  etaSignalMcFromTnP   -> Sumw2();   
  etaSignalMcFromTnPWW -> Sumw2();   
  etaSignalMcNoTnP     -> Sumw2(); 
  etaSignalMcFromTnP   -> Scale(1./etaSignalMcFromTnP->Integral());   
  etaSignalMcFromTnPWW -> Scale(1./etaSignalMcFromTnPWW->Integral());   
  etaSignalMcNoTnP     -> Scale(1./etaSignalMcNoTnP->Integral());   
  //
  TH1F *etaSignalMcFromTnP0   = (TH1F*)fileTnP.Get("probeEtaSignalMc0");
  TH1F *etaSignalMcFromTnP0WW = (TH1F*)fileTnP.Get("probeEtaSignalMcWW0");
  TH1F *etaSignalMcNoTnP0     = (TH1F*)fileNoTnP.Get("etaSignalMc0");
  etaSignalMcFromTnP0   -> Sumw2();   
  etaSignalMcFromTnP0WW -> Sumw2();   
  etaSignalMcNoTnP0     -> Sumw2(); 
  etaSignalMcFromTnP0   -> Scale(1./etaSignalMcFromTnP0->Integral());   
  etaSignalMcFromTnP0WW -> Scale(1./etaSignalMcFromTnP0WW->Integral());   
  etaSignalMcNoTnP0     -> Scale(1./etaSignalMcNoTnP0->Integral());   
  //
  TH1F *etaSignalMcFromTnP1   = (TH1F*)fileTnP.Get("probeEtaSignalMc1");
  TH1F *etaSignalMcFromTnP1WW = (TH1F*)fileTnP.Get("probeEtaSignalMcWW1");
  TH1F *etaSignalMcNoTnP1     = (TH1F*)fileNoTnP.Get("etaSignalMc1");
  etaSignalMcFromTnP1   -> Sumw2();   
  etaSignalMcFromTnP1WW -> Sumw2();   
  etaSignalMcNoTnP1     -> Sumw2(); 
  etaSignalMcFromTnP1   -> Scale(1./etaSignalMcFromTnP1->Integral());   
  etaSignalMcFromTnP1WW -> Scale(1./etaSignalMcFromTnP1WW->Integral());   
  etaSignalMcNoTnP1     -> Scale(1./etaSignalMcNoTnP1->Integral());   
  //
  TH1F *etaSignalMcFromTnP2   = (TH1F*)fileTnP.Get("probeEtaSignalMc2");
  TH1F *etaSignalMcFromTnP2WW = (TH1F*)fileTnP.Get("probeEtaSignalMcWW2");
  TH1F *etaSignalMcNoTnP2     = (TH1F*)fileNoTnP.Get("etaSignalMc2");
  etaSignalMcFromTnP2   -> Sumw2();   
  etaSignalMcFromTnP2WW -> Sumw2();   
  etaSignalMcNoTnP2     -> Sumw2(); 
  etaSignalMcFromTnP2   -> Scale(1./etaSignalMcFromTnP2->Integral());   
  etaSignalMcFromTnP2WW -> Scale(1./etaSignalMcFromTnP2WW->Integral());   
  etaSignalMcNoTnP2     -> Scale(1./etaSignalMcNoTnP2->Integral());   
  //
  TH1F *etaSignalMcFromTnP3   = (TH1F*)fileTnP.Get("probeEtaSignalMc3");
  TH1F *etaSignalMcFromTnP3WW = (TH1F*)fileTnP.Get("probeEtaSignalMcWW3");
  TH1F *etaSignalMcNoTnP3     = (TH1F*)fileNoTnP.Get("etaSignalMc3");
  etaSignalMcFromTnP3   -> Sumw2();   
  etaSignalMcFromTnP3WW -> Sumw2();   
  etaSignalMcNoTnP3     -> Sumw2(); 
  etaSignalMcFromTnP3   -> Scale(1./etaSignalMcFromTnP3->Integral());   
  etaSignalMcFromTnP3WW -> Scale(1./etaSignalMcFromTnP3WW->Integral());   
  etaSignalMcNoTnP3     -> Scale(1./etaSignalMcNoTnP3->Integral());   
  //
  TH1F *etaSignalMcFromTnP4   = (TH1F*)fileTnP.Get("probeEtaSignalMc4");
  TH1F *etaSignalMcFromTnP4WW = (TH1F*)fileTnP.Get("probeEtaSignalMcWW4");
  TH1F *etaSignalMcNoTnP4     = (TH1F*)fileNoTnP.Get("etaSignalMc4");
  etaSignalMcFromTnP4   -> Sumw2();   
  etaSignalMcFromTnP4WW -> Sumw2();   
  etaSignalMcNoTnP4     -> Sumw2(); 
  etaSignalMcFromTnP4   -> Scale(1./etaSignalMcFromTnP4->Integral());   
  etaSignalMcFromTnP4WW -> Scale(1./etaSignalMcFromTnP4WW->Integral());   
  etaSignalMcNoTnP4     -> Scale(1./etaSignalMcNoTnP4->Integral());   
  
  // Histos: mva
  TH1F *mvaSignalMcFromTnP   = (TH1F*)fileTnP.Get("probeMvaSignalMc");
  TH1F *mvaSignalMcFromTnPWW = (TH1F*)fileTnP.Get("probeMvaSignalMcWW");
  TH1F *mvaSignalMcNoTnP     = (TH1F*)fileNoTnP.Get("mvaSignalMc");
  mvaSignalMcFromTnP   -> Sumw2();   
  mvaSignalMcFromTnPWW -> Sumw2();   
  mvaSignalMcNoTnP     -> Sumw2(); 
  mvaSignalMcFromTnP   -> Scale(1./mvaSignalMcFromTnP->Integral());   
  mvaSignalMcFromTnPWW -> Scale(1./mvaSignalMcFromTnPWW->Integral());   
  mvaSignalMcNoTnP     -> Scale(1./mvaSignalMcNoTnP->Integral());   

  // Histos: mva in small bins
  TH1F *mvaSignalEBMc0FromTnP = (TH1F*)fileTnP.Get("mvaSignalEBMc0");
  TH1F *mvaSignalEBMc1FromTnP = (TH1F*)fileTnP.Get("mvaSignalEBMc1");
  TH1F *mvaSignalEBMc2FromTnP = (TH1F*)fileTnP.Get("mvaSignalEBMc2");
  TH1F *mvaSignalEBMc3FromTnP = (TH1F*)fileTnP.Get("mvaSignalEBMc3");
  TH1F *mvaSignalEBMc4FromTnP = (TH1F*)fileTnP.Get("mvaSignalEBMc4");
  TH1F *mvaSignalEEMc0FromTnP = (TH1F*)fileTnP.Get("mvaSignalEEMc0");
  TH1F *mvaSignalEEMc1FromTnP = (TH1F*)fileTnP.Get("mvaSignalEEMc1");
  TH1F *mvaSignalEEMc2FromTnP = (TH1F*)fileTnP.Get("mvaSignalEEMc2");
  TH1F *mvaSignalEEMc3FromTnP = (TH1F*)fileTnP.Get("mvaSignalEEMc3");
  TH1F *mvaSignalEEMc4FromTnP = (TH1F*)fileTnP.Get("mvaSignalEEMc4");
  mvaSignalEBMc0FromTnP->Sumw2();
  mvaSignalEBMc1FromTnP->Sumw2();
  mvaSignalEBMc2FromTnP->Sumw2();
  mvaSignalEBMc3FromTnP->Sumw2();
  mvaSignalEBMc4FromTnP->Sumw2();
  mvaSignalEEMc0FromTnP->Sumw2();
  mvaSignalEEMc1FromTnP->Sumw2();
  mvaSignalEEMc2FromTnP->Sumw2();
  mvaSignalEEMc3FromTnP->Sumw2();
  mvaSignalEEMc4FromTnP->Sumw2();
  mvaSignalEBMc0FromTnP->Scale(1./mvaSignalEBMc0FromTnP->Integral());
  mvaSignalEBMc1FromTnP->Scale(1./mvaSignalEBMc1FromTnP->Integral());
  mvaSignalEBMc2FromTnP->Scale(1./mvaSignalEBMc2FromTnP->Integral());
  mvaSignalEBMc3FromTnP->Scale(1./mvaSignalEBMc3FromTnP->Integral());
  mvaSignalEBMc4FromTnP->Scale(1./mvaSignalEBMc4FromTnP->Integral());
  mvaSignalEEMc0FromTnP->Scale(1./mvaSignalEEMc0FromTnP->Integral());
  mvaSignalEEMc1FromTnP->Scale(1./mvaSignalEEMc1FromTnP->Integral());
  mvaSignalEEMc2FromTnP->Scale(1./mvaSignalEEMc2FromTnP->Integral());
  mvaSignalEEMc3FromTnP->Scale(1./mvaSignalEEMc3FromTnP->Integral());
  mvaSignalEEMc4FromTnP->Scale(1./mvaSignalEEMc4FromTnP->Integral());
  //
  TH1F *mvaSignalEBMc0NoTnP = (TH1F*)fileNoTnP.Get("mvaSignalEBMc0");
  TH1F *mvaSignalEBMc1NoTnP = (TH1F*)fileNoTnP.Get("mvaSignalEBMc1");
  TH1F *mvaSignalEBMc2NoTnP = (TH1F*)fileNoTnP.Get("mvaSignalEBMc2");
  TH1F *mvaSignalEBMc3NoTnP = (TH1F*)fileNoTnP.Get("mvaSignalEBMc3");
  TH1F *mvaSignalEBMc4NoTnP = (TH1F*)fileNoTnP.Get("mvaSignalEBMc4");
  TH1F *mvaSignalEEMc0NoTnP = (TH1F*)fileNoTnP.Get("mvaSignalEEMc0");
  TH1F *mvaSignalEEMc1NoTnP = (TH1F*)fileNoTnP.Get("mvaSignalEEMc1");
  TH1F *mvaSignalEEMc2NoTnP = (TH1F*)fileNoTnP.Get("mvaSignalEEMc2");
  TH1F *mvaSignalEEMc3NoTnP = (TH1F*)fileNoTnP.Get("mvaSignalEEMc3");
  TH1F *mvaSignalEEMc4NoTnP = (TH1F*)fileNoTnP.Get("mvaSignalEEMc4");
  mvaSignalEBMc0NoTnP->Sumw2();
  mvaSignalEBMc1NoTnP->Sumw2();
  mvaSignalEBMc2NoTnP->Sumw2();
  mvaSignalEBMc3NoTnP->Sumw2();
  mvaSignalEBMc4NoTnP->Sumw2();
  mvaSignalEEMc0NoTnP->Sumw2();
  mvaSignalEEMc1NoTnP->Sumw2();
  mvaSignalEEMc2NoTnP->Sumw2();
  mvaSignalEEMc3NoTnP->Sumw2();
  mvaSignalEEMc4NoTnP->Sumw2();
  mvaSignalEBMc0NoTnP->Scale(1./mvaSignalEBMc0NoTnP->Integral());
  mvaSignalEBMc1NoTnP->Scale(1./mvaSignalEBMc1NoTnP->Integral());
  mvaSignalEBMc2NoTnP->Scale(1./mvaSignalEBMc2NoTnP->Integral());
  mvaSignalEBMc3NoTnP->Scale(1./mvaSignalEBMc3NoTnP->Integral());
  mvaSignalEBMc4NoTnP->Scale(1./mvaSignalEBMc4NoTnP->Integral());
  mvaSignalEEMc0NoTnP->Scale(1./mvaSignalEEMc0NoTnP->Integral());
  mvaSignalEEMc1NoTnP->Scale(1./mvaSignalEEMc1NoTnP->Integral());
  mvaSignalEEMc2NoTnP->Scale(1./mvaSignalEEMc2NoTnP->Integral());
  mvaSignalEEMc3NoTnP->Scale(1./mvaSignalEEMc3NoTnP->Integral());
  mvaSignalEEMc4NoTnP->Scale(1./mvaSignalEEMc4NoTnP->Integral());
  //
  TH1F *mvaSignalEBMc0WWFromTnP = (TH1F*)fileTnP.Get("mvaSignalEBMc0WW");
  TH1F *mvaSignalEBMc1WWFromTnP = (TH1F*)fileTnP.Get("mvaSignalEBMc1WW");
  TH1F *mvaSignalEBMc2WWFromTnP = (TH1F*)fileTnP.Get("mvaSignalEBMc2WW");
  TH1F *mvaSignalEBMc3WWFromTnP = (TH1F*)fileTnP.Get("mvaSignalEBMc3WW");
  TH1F *mvaSignalEBMc4WWFromTnP = (TH1F*)fileTnP.Get("mvaSignalEBMc4WW");
  TH1F *mvaSignalEEMc0WWFromTnP = (TH1F*)fileTnP.Get("mvaSignalEEMc0WW");
  TH1F *mvaSignalEEMc1WWFromTnP = (TH1F*)fileTnP.Get("mvaSignalEEMc1WW");
  TH1F *mvaSignalEEMc2WWFromTnP = (TH1F*)fileTnP.Get("mvaSignalEEMc2WW");
  TH1F *mvaSignalEEMc3WWFromTnP = (TH1F*)fileTnP.Get("mvaSignalEEMc3WW");
  TH1F *mvaSignalEEMc4WWFromTnP = (TH1F*)fileTnP.Get("mvaSignalEEMc4WW");
  mvaSignalEBMc0WWFromTnP->Sumw2();
  mvaSignalEBMc1WWFromTnP->Sumw2();
  mvaSignalEBMc2WWFromTnP->Sumw2();
  mvaSignalEBMc3WWFromTnP->Sumw2();
  mvaSignalEBMc4WWFromTnP->Sumw2();
  mvaSignalEEMc0WWFromTnP->Sumw2();
  mvaSignalEEMc1WWFromTnP->Sumw2();
  mvaSignalEEMc2WWFromTnP->Sumw2();
  mvaSignalEEMc3WWFromTnP->Sumw2();
  mvaSignalEEMc4WWFromTnP->Sumw2();
  mvaSignalEBMc0WWFromTnP->Scale(1./mvaSignalEBMc0WWFromTnP->Integral());
  mvaSignalEBMc1WWFromTnP->Scale(1./mvaSignalEBMc1WWFromTnP->Integral());
  mvaSignalEBMc2WWFromTnP->Scale(1./mvaSignalEBMc2WWFromTnP->Integral());
  mvaSignalEBMc3WWFromTnP->Scale(1./mvaSignalEBMc3WWFromTnP->Integral());
  mvaSignalEBMc4WWFromTnP->Scale(1./mvaSignalEBMc4WWFromTnP->Integral());
  mvaSignalEEMc0WWFromTnP->Scale(1./mvaSignalEEMc0WWFromTnP->Integral());
  mvaSignalEEMc1WWFromTnP->Scale(1./mvaSignalEEMc1WWFromTnP->Integral());
  mvaSignalEEMc2WWFromTnP->Scale(1./mvaSignalEEMc2WWFromTnP->Integral());
  mvaSignalEEMc3WWFromTnP->Scale(1./mvaSignalEEMc3WWFromTnP->Integral());
  mvaSignalEEMc4WWFromTnP->Scale(1./mvaSignalEEMc4WWFromTnP->Integral());

  // 2Dim histos: pT vs eta
  TH2F *ptVsEtaSignalMcFromTnP = (TH2F*)fileTnP.Get("probePtVsEtaSignalMc");
  TH2F *ptVsEtaSignalMcNoTnP   = (TH2F*)fileNoTnP.Get("ptVsEtaSignalMc");
  //ptVsEtaSignalMcFromTnP -> Sumw2();   
  //ptVsEtaSignalMcNoTnP   -> Sumw2();   
  //ptVsEtaSignalMcFromTnP -> Scale(1./ptVsEtaSignalMcFromTnP->Integral());   
  //ptVsEtaSignalMcNoTnP   -> Scale(1./ptVsEtaSignalMcNoTnP->Integral());   
  //
  TH2F *ptVsEtaSignalMcFromTnPWW = (TH2F*)fileTnP.Get("probePtVsEtaSignalMcWW");
  //ptVsEtaSignalMcFromTnPWW -> Sumw2();   
  //ptVsEtaSignalMcFromTnPWW -> Scale(1./ptVsEtaSignalMcFromTnPWW->Integral());   

  // --------------------------------------------------------
  // Compute weights
  TH2F *ptVsEtaSignalWeights = (TH2F*)ptVsEtaSignalMcNoTnP->Clone("ptVsEtaSignalWeights");
  ptVsEtaSignalWeights->Sumw2(); 
  ptVsEtaSignalWeights->Divide(ptVsEtaSignalMcFromTnP);
  ptVsEtaSignalWeights->SetTitle("ptVsEtaSignalWeights");
  ptVsEtaSignalWeights->SetName("ptVsEtaSignalWeights");
  cout << "signal pre: min = " << ptVsEtaSignalWeights->GetMinimum() << ", max = " << ptVsEtaSignalWeights->GetMaximum() << endl;
  for (int iBinEta=0; iBinEta<(ptVsEtaSignalWeights->GetNbinsX()); iBinEta++) {
    for (int iBinPt=0; iBinPt<(ptVsEtaSignalWeights->GetNbinsY()); iBinPt++) {
      int iBinEtaP = iBinEta+1;
      int iBinPtP  = iBinPt+1;
      if ((ptVsEtaSignalWeights->GetBinContent(iBinEtaP,iBinPtP))>2) ptVsEtaSignalWeights->SetBinContent(iBinEtaP,iBinPtP,2);
    }}
  cout << "signal post: min = " << ptVsEtaSignalWeights->GetMinimum() << ", max = " << ptVsEtaSignalWeights->GetMaximum() << endl;
  ptVsEtaSignalWeights->Scale(1./ptVsEtaSignalWeights->Integral()); 
  cout << "signal scaled: min = " << ptVsEtaSignalWeights->GetMinimum() << ", max = " << ptVsEtaSignalWeights->GetMaximum() << endl;

  // --------------------------------------------------------
  // Cosmetics
  ptSignalMcFromTnP   -> SetLineWidth(2);  
  ptSignalMcFromTnPWW -> SetLineWidth(2);  
  ptSignalMcNoTnP     -> SetLineWidth(2);  
  ptSignalMcFromTnP   -> SetLineColor(2);  
  ptSignalMcFromTnPWW -> SetLineColor(3);  
  ptSignalMcNoTnP     -> SetLineColor(4);  
  ptSignalMcFromTnP   -> SetTitle("");
  ptSignalMcFromTnPWW -> SetTitle("");
  ptSignalMcNoTnP     -> SetTitle("");
  ptSignalMcFromTnP   -> GetXaxis()->SetTitle("pT [GeV]");
  ptSignalMcFromTnPWW -> GetXaxis()->SetTitle("pT [GeV]");
  ptSignalMcNoTnP     -> GetXaxis()->SetTitle("pT [GeV]");
  //
  etaSignalMcFromTnP   -> SetLineWidth(2);  
  etaSignalMcFromTnPWW -> SetLineWidth(2);  
  etaSignalMcNoTnP     -> SetLineWidth(2);  
  etaSignalMcFromTnP   -> SetLineColor(2);  
  etaSignalMcFromTnPWW -> SetLineColor(3);  
  etaSignalMcNoTnP     -> SetLineColor(4);  
  etaSignalMcFromTnP   -> SetTitle("");
  etaSignalMcFromTnPWW -> SetTitle("");
  etaSignalMcNoTnP     -> SetTitle("");
  etaSignalMcFromTnP   -> GetXaxis()->SetTitle("#eta");
  etaSignalMcFromTnPWW -> GetXaxis()->SetTitle("#eta");
  etaSignalMcNoTnP     -> GetXaxis()->SetTitle("#eta");
  //
  etaSignalMcFromTnP0   -> SetLineWidth(2);  
  etaSignalMcFromTnP0WW -> SetLineWidth(2);  
  etaSignalMcNoTnP0     -> SetLineWidth(2);  
  etaSignalMcFromTnP0   -> SetLineColor(2);  
  etaSignalMcFromTnP0WW -> SetLineColor(3);  
  etaSignalMcNoTnP0     -> SetLineColor(4);  
  etaSignalMcFromTnP0   -> SetTitle("");
  etaSignalMcFromTnP0WW -> SetTitle("");
  etaSignalMcNoTnP0     -> SetTitle("");
  etaSignalMcFromTnP0   -> GetXaxis()->SetTitle("#eta");
  etaSignalMcFromTnP0WW -> GetXaxis()->SetTitle("#eta");
  etaSignalMcNoTnP0     -> GetXaxis()->SetTitle("#eta");
  //
  etaSignalMcFromTnP1   -> SetLineWidth(2);  
  etaSignalMcFromTnP1WW -> SetLineWidth(2);  
  etaSignalMcNoTnP1     -> SetLineWidth(2);  
  etaSignalMcFromTnP1   -> SetLineColor(2);  
  etaSignalMcFromTnP1WW -> SetLineColor(3);  
  etaSignalMcNoTnP1     -> SetLineColor(4);  
  etaSignalMcFromTnP1   -> SetTitle("");
  etaSignalMcFromTnP1WW -> SetTitle("");
  etaSignalMcNoTnP1     -> SetTitle("");
  etaSignalMcFromTnP1   -> GetXaxis()->SetTitle("#eta");
  etaSignalMcFromTnP1WW -> GetXaxis()->SetTitle("#eta");
  etaSignalMcNoTnP1     -> GetXaxis()->SetTitle("#eta");
  //
  etaSignalMcFromTnP2   -> SetLineWidth(2);  
  etaSignalMcFromTnP2WW -> SetLineWidth(2);  
  etaSignalMcNoTnP2     -> SetLineWidth(2);  
  etaSignalMcFromTnP2   -> SetLineColor(2);  
  etaSignalMcFromTnP2WW -> SetLineColor(3);  
  etaSignalMcNoTnP2     -> SetLineColor(4);  
  etaSignalMcFromTnP2   -> SetTitle("");
  etaSignalMcFromTnP2WW -> SetTitle("");
  etaSignalMcNoTnP2     -> SetTitle("");
  etaSignalMcFromTnP2   -> GetXaxis()->SetTitle("#eta");
  etaSignalMcFromTnP2WW -> GetXaxis()->SetTitle("#eta");
  etaSignalMcNoTnP2     -> GetXaxis()->SetTitle("#eta");
  //
  etaSignalMcFromTnP3   -> SetLineWidth(2);  
  etaSignalMcFromTnP3WW -> SetLineWidth(2);  
  etaSignalMcNoTnP3     -> SetLineWidth(2);  
  etaSignalMcFromTnP3   -> SetLineColor(2);  
  etaSignalMcFromTnP3WW -> SetLineColor(3);  
  etaSignalMcNoTnP3     -> SetLineColor(4);  
  etaSignalMcFromTnP3   -> SetTitle("");
  etaSignalMcFromTnP3WW -> SetTitle("");
  etaSignalMcNoTnP3     -> SetTitle("");
  etaSignalMcFromTnP3   -> GetXaxis()->SetTitle("#eta");
  etaSignalMcFromTnP3WW -> GetXaxis()->SetTitle("#eta");
  etaSignalMcNoTnP3     -> GetXaxis()->SetTitle("#eta");
  //
  etaSignalMcFromTnP4   -> SetLineWidth(2);  
  etaSignalMcFromTnP4WW -> SetLineWidth(2);  
  etaSignalMcNoTnP4     -> SetLineWidth(2);  
  etaSignalMcFromTnP4   -> SetLineColor(2);  
  etaSignalMcFromTnP4WW -> SetLineColor(3);  
  etaSignalMcNoTnP4     -> SetLineColor(4);  
  etaSignalMcFromTnP4   -> SetTitle("");
  etaSignalMcFromTnP4WW -> SetTitle("");
  etaSignalMcNoTnP4     -> SetTitle("");
  etaSignalMcFromTnP4   -> GetXaxis()->SetTitle("#eta");
  etaSignalMcFromTnP4WW -> GetXaxis()->SetTitle("#eta");
  etaSignalMcNoTnP4     -> GetXaxis()->SetTitle("#eta");
  //
  mvaSignalMcFromTnP   -> SetLineWidth(2);  
  mvaSignalMcFromTnPWW -> SetLineWidth(2);  
  mvaSignalMcNoTnP     -> SetLineWidth(2);  
  mvaSignalMcFromTnP   -> SetLineColor(2);  
  mvaSignalMcFromTnPWW -> SetLineColor(3);  
  mvaSignalMcNoTnP     -> SetLineColor(4);  
  mvaSignalMcFromTnP   -> SetTitle("");
  mvaSignalMcFromTnPWW -> SetTitle("");
  mvaSignalMcNoTnP     -> SetTitle("");
  mvaSignalMcFromTnP   -> GetXaxis()->SetTitle("ID BDT output");
  mvaSignalMcFromTnPWW -> GetXaxis()->SetTitle("ID BDT output");
  mvaSignalMcNoTnP     -> GetXaxis()->SetTitle("ID BDT output");
  //
  mvaSignalEBMc0FromTnP->SetLineWidth(2);
  mvaSignalEBMc0FromTnP->SetLineColor(2);
  mvaSignalEBMc0FromTnP->SetTitle("");
  mvaSignalEBMc0FromTnP->GetXaxis()->SetTitle("ID BDT output");
  mvaSignalEBMc1FromTnP->SetLineWidth(2);
  mvaSignalEBMc1FromTnP->SetLineColor(2);
  mvaSignalEBMc1FromTnP->SetTitle("");
  mvaSignalEBMc1FromTnP->GetXaxis()->SetTitle("ID BDT output");
  mvaSignalEBMc2FromTnP->SetLineWidth(2);
  mvaSignalEBMc2FromTnP->SetLineColor(2);
  mvaSignalEBMc2FromTnP->SetTitle("");
  mvaSignalEBMc2FromTnP->GetXaxis()->SetTitle("ID BDT output");
  mvaSignalEBMc3FromTnP->SetLineWidth(2);
  mvaSignalEBMc3FromTnP->SetLineColor(2);
  mvaSignalEBMc3FromTnP->SetTitle("");
  mvaSignalEBMc3FromTnP->GetXaxis()->SetTitle("ID BDT output");
  mvaSignalEBMc4FromTnP->SetLineWidth(2);
  mvaSignalEBMc4FromTnP->SetLineColor(2);
  mvaSignalEBMc4FromTnP->SetTitle("");
  mvaSignalEBMc4FromTnP->GetXaxis()->SetTitle("ID BDT output");
  mvaSignalEEMc0FromTnP->SetLineWidth(2);
  mvaSignalEEMc0FromTnP->SetLineColor(2);
  mvaSignalEEMc0FromTnP->SetTitle("");
  mvaSignalEEMc0FromTnP->GetXaxis()->SetTitle("ID BDT output");
  mvaSignalEEMc1FromTnP->SetLineWidth(2);
  mvaSignalEEMc1FromTnP->SetLineColor(2);
  mvaSignalEEMc1FromTnP->SetTitle("");
  mvaSignalEEMc1FromTnP->GetXaxis()->SetTitle("ID BDT output");
  mvaSignalEEMc2FromTnP->SetLineWidth(2);
  mvaSignalEEMc2FromTnP->SetLineColor(2);
  mvaSignalEEMc2FromTnP->SetTitle("");
  mvaSignalEEMc2FromTnP->GetXaxis()->SetTitle("ID BDT output");
  mvaSignalEEMc3FromTnP->SetLineWidth(2);
  mvaSignalEEMc3FromTnP->SetLineColor(2);
  mvaSignalEEMc3FromTnP->SetTitle("");
  mvaSignalEEMc3FromTnP->GetXaxis()->SetTitle("ID BDT output");
  mvaSignalEEMc4FromTnP->SetLineWidth(2);
  mvaSignalEEMc4FromTnP->SetLineColor(2);
  mvaSignalEEMc4FromTnP->SetTitle("");
  mvaSignalEEMc4FromTnP->GetXaxis()->SetTitle("ID BDT output");
  //
  mvaSignalEBMc0WWFromTnP->SetLineWidth(2);
  mvaSignalEBMc0WWFromTnP->SetLineColor(3);
  mvaSignalEBMc0WWFromTnP->SetTitle("");
  mvaSignalEBMc0WWFromTnP->GetXaxis()->SetTitle("ID BDT output");
  mvaSignalEBMc1WWFromTnP->SetLineWidth(2);
  mvaSignalEBMc1WWFromTnP->SetLineColor(3);
  mvaSignalEBMc1WWFromTnP->SetTitle("");
  mvaSignalEBMc1WWFromTnP->GetXaxis()->SetTitle("ID BDT output");
  mvaSignalEBMc2WWFromTnP->SetLineWidth(2);
  mvaSignalEBMc2WWFromTnP->SetLineColor(3);
  mvaSignalEBMc2WWFromTnP->SetTitle("");
  mvaSignalEBMc2WWFromTnP->GetXaxis()->SetTitle("ID BDT output");
  mvaSignalEBMc3WWFromTnP->SetLineWidth(2);
  mvaSignalEBMc3WWFromTnP->SetLineColor(3);
  mvaSignalEBMc3WWFromTnP->SetTitle("");
  mvaSignalEBMc3WWFromTnP->GetXaxis()->SetTitle("ID BDT output");
  mvaSignalEBMc4WWFromTnP->SetLineWidth(2);
  mvaSignalEBMc4WWFromTnP->SetLineColor(3);
  mvaSignalEBMc4WWFromTnP->SetTitle("");
  mvaSignalEBMc4WWFromTnP->GetXaxis()->SetTitle("ID BDT output");
  mvaSignalEEMc0WWFromTnP->SetLineWidth(2);
  mvaSignalEEMc0WWFromTnP->SetLineColor(3);
  mvaSignalEEMc0WWFromTnP->SetTitle("");
  mvaSignalEEMc0WWFromTnP->GetXaxis()->SetTitle("ID BDT output");
  mvaSignalEEMc1WWFromTnP->SetLineWidth(2);
  mvaSignalEEMc1WWFromTnP->SetLineColor(3);
  mvaSignalEEMc1WWFromTnP->SetTitle("");
  mvaSignalEEMc1WWFromTnP->GetXaxis()->SetTitle("ID BDT output");
  mvaSignalEEMc2WWFromTnP->SetLineWidth(2);
  mvaSignalEEMc2WWFromTnP->SetLineColor(3);
  mvaSignalEEMc2WWFromTnP->SetTitle("");
  mvaSignalEEMc2WWFromTnP->GetXaxis()->SetTitle("ID BDT output");
  mvaSignalEEMc3WWFromTnP->SetLineWidth(2);
  mvaSignalEEMc3WWFromTnP->SetLineColor(3);
  mvaSignalEEMc3WWFromTnP->SetTitle("");
  mvaSignalEEMc3WWFromTnP->GetXaxis()->SetTitle("ID BDT output");
  mvaSignalEEMc4WWFromTnP->SetLineWidth(2);
  mvaSignalEEMc4WWFromTnP->SetLineColor(3);
  mvaSignalEEMc4WWFromTnP->SetTitle("");
  mvaSignalEEMc4WWFromTnP->GetXaxis()->SetTitle("ID BDT output");
  //
  mvaSignalEBMc0NoTnP->SetLineWidth(2);
  mvaSignalEBMc0NoTnP->SetLineColor(4);
  mvaSignalEBMc0NoTnP->SetTitle("");
  mvaSignalEBMc0NoTnP->GetXaxis()->SetTitle("ID BDT output");
  mvaSignalEBMc1NoTnP->SetLineWidth(2);
  mvaSignalEBMc1NoTnP->SetLineColor(4);
  mvaSignalEBMc1NoTnP->SetTitle("");
  mvaSignalEBMc1NoTnP->GetXaxis()->SetTitle("ID BDT output");
  mvaSignalEBMc2NoTnP->SetLineWidth(2);
  mvaSignalEBMc2NoTnP->SetLineColor(4);
  mvaSignalEBMc2NoTnP->SetTitle("");
  mvaSignalEBMc2NoTnP->GetXaxis()->SetTitle("ID BDT output");
  mvaSignalEBMc3NoTnP->SetLineWidth(2);
  mvaSignalEBMc3NoTnP->SetLineColor(4);
  mvaSignalEBMc3NoTnP->SetTitle("");
  mvaSignalEBMc3NoTnP->GetXaxis()->SetTitle("ID BDT output");
  mvaSignalEBMc4NoTnP->SetLineWidth(2);
  mvaSignalEBMc4NoTnP->SetLineColor(4);
  mvaSignalEBMc4NoTnP->SetTitle("");
  mvaSignalEBMc4NoTnP->GetXaxis()->SetTitle("ID BDT output");
  mvaSignalEEMc0NoTnP->SetLineWidth(2);
  mvaSignalEEMc0NoTnP->SetLineColor(4);
  mvaSignalEEMc0NoTnP->SetTitle("");
  mvaSignalEEMc0NoTnP->GetXaxis()->SetTitle("ID BDT output");
  mvaSignalEEMc1NoTnP->SetLineWidth(2);
  mvaSignalEEMc1NoTnP->SetLineColor(4);
  mvaSignalEEMc1NoTnP->SetTitle("");
  mvaSignalEEMc1NoTnP->GetXaxis()->SetTitle("ID BDT output");
  mvaSignalEEMc2NoTnP->SetLineWidth(2);
  mvaSignalEEMc2NoTnP->SetLineColor(4);
  mvaSignalEEMc2NoTnP->SetTitle("");
  mvaSignalEEMc2NoTnP->GetXaxis()->SetTitle("ID BDT output");
  mvaSignalEEMc3NoTnP->SetLineWidth(2);
  mvaSignalEEMc3NoTnP->SetLineColor(4);
  mvaSignalEEMc3NoTnP->SetTitle("");
  mvaSignalEEMc3NoTnP->GetXaxis()->SetTitle("ID BDT output");
  mvaSignalEEMc4NoTnP->SetLineWidth(2);
  mvaSignalEEMc4NoTnP->SetLineColor(4);
  mvaSignalEEMc4NoTnP->SetTitle("");
  mvaSignalEEMc4NoTnP->GetXaxis()->SetTitle("ID BDT output");


  // -------------------------------------------------------------
  // Plots
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);

  TLegend *leg;
  leg = new TLegend(0.45,0.55,0.75,0.80);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.05);
  leg->SetFillColor(0);
  leg->AddEntry(ptSignalMcFromTnP, "With tnp", "lp");
  leg->AddEntry(ptSignalMcNoTnP,   "No tnp", "lp");
  //
  TLegend *legww;
  legww = new TLegend(0.45,0.55,0.75,0.80);
  legww->SetFillStyle(0);
  legww->SetBorderSize(0);
  legww->SetTextSize(0.05);
  legww->SetFillColor(0);
  legww->AddEntry(ptSignalMcFromTnPWW, "With tnp, weighted", "lp");
  legww->AddEntry(ptSignalMcNoTnP,   "No tnp", "lp");

  TLegend *leg2;
  leg2 = new TLegend(0.10,0.65,0.40,0.90);
  leg2->SetFillStyle(0);
  leg2->SetBorderSize(0);
  leg2->SetTextSize(0.05);
  leg2->SetFillColor(0);
  leg2->AddEntry(ptSignalMcFromTnP, "With tnp", "lp");
  leg2->AddEntry(ptSignalMcNoTnP,   "No tnp", "lp");
  //
  TLegend *leg2ww;
  leg2ww = new TLegend(0.10,0.65,0.40,0.90);
  leg2ww->SetFillStyle(0);
  leg2ww->SetBorderSize(0);
  leg2ww->SetTextSize(0.05);
  leg2ww->SetFillColor(0);
  leg2ww->AddEntry(ptSignalMcFromTnPWW, "With tnp, weighted", "lp");
  leg2ww->AddEntry(ptSignalMcNoTnP,   "No tnp", "lp");

  TLegend *leg3;
  leg3 = new TLegend(0.35,0.10,0.65,0.35);
  leg3->SetFillStyle(0);
  leg3->SetBorderSize(0);
  leg3->SetTextSize(0.05);
  leg3->SetFillColor(0);
  leg3->AddEntry(ptSignalMcFromTnP, "With tnp", "lp");
  leg3->AddEntry(ptSignalMcNoTnP,   "No tnp", "lp");
  //
  TLegend *leg3ww;
  leg3ww = new TLegend(0.35,0.10,0.65,0.35);
  leg3ww->SetFillStyle(0);
  leg3ww->SetBorderSize(0);
  leg3ww->SetTextSize(0.05);
  leg3ww->SetFillColor(0);
  leg3ww->AddEntry(ptSignalMcFromTnPWW, "With tnp, weighted", "lp");
  leg3ww->AddEntry(ptSignalMcNoTnP,   "No tnp", "lp");

  TCanvas cs1("cs1","cs1",1);
  ptSignalMcNoTnP->DrawNormalized("hist");
  ptSignalMcFromTnP->DrawNormalized("samehist");
  leg->Draw();
  cs1.SaveAs("ptSignal_naniVsTnp.png");

  TCanvas cs2("cs2","cs2",1);
  ptSignalMcNoTnP->DrawNormalized("hist");
  ptSignalMcFromTnPWW->DrawNormalized("samehist");
  legww->Draw();
  cs2.SaveAs("ptSignal_naniVsTnp_withWeight.png");
  
  // -----------

  TCanvas cs3("cs3","cs3",1);
  etaSignalMcFromTnP->DrawNormalized("hist");
  etaSignalMcNoTnP->DrawNormalized("samehist");
  leg3->Draw();
  cs3.SaveAs("etaSignal_naniVsTnp.png");

  TCanvas cs3A("cs3A","cs3A",1);
  etaSignalMcFromTnP0->DrawNormalized("hist");
  etaSignalMcNoTnP0->DrawNormalized("samehist");
  leg3->Draw();
  cs3A.SaveAs("etaSignal_naniVsTnp0.png");

  TCanvas cs3B("cs3B","cs3B",1);
  etaSignalMcFromTnP1->DrawNormalized("hist");
  etaSignalMcNoTnP1->DrawNormalized("samehist");
  leg3->Draw();
  cs3B.SaveAs("etaSignal_naniVsTnp1.png");

  TCanvas cs3C("cs3C","cs3C",1);
  etaSignalMcFromTnP2->DrawNormalized("hist");
  etaSignalMcNoTnP2->DrawNormalized("samehist");
  leg3->Draw();
  cs3C.SaveAs("etaSignal_naniVsTnp2.png");

  TCanvas cs3D("cs3D","cs3D",1);
  etaSignalMcFromTnP3->DrawNormalized("hist");
  etaSignalMcNoTnP3->DrawNormalized("samehist");
  leg3->Draw();
  cs3D.SaveAs("etaSignal_naniVsTnp3.png");

  TCanvas cs3E("cs3E","cs3E",1);
  etaSignalMcFromTnP4->DrawNormalized("hist");
  etaSignalMcNoTnP4->DrawNormalized("samehist");
  leg3->Draw();
  cs3E.SaveAs("etaSignal_naniVsTnp4.png");

  TCanvas cs4("cs4","cs4",1);
  etaSignalMcFromTnPWW->DrawNormalized("hist");
  etaSignalMcNoTnP->DrawNormalized("samehist");
  leg3ww->Draw();
  cs4.SaveAs("etaSignal_naniVsTnp_withWeight.png");

  TCanvas cs4A("cs4A","cs4A",1);
  etaSignalMcFromTnP0WW->DrawNormalized("hist");
  etaSignalMcNoTnP0->DrawNormalized("samehist");
  leg3ww->Draw();
  cs4A.SaveAs("etaSignal_naniVsTnp0_withWeight.png");

  TCanvas cs4B("cs4B","cs4B",1);
  etaSignalMcFromTnP1WW->DrawNormalized("hist");
  etaSignalMcNoTnP1->DrawNormalized("samehist");
  leg3ww->Draw();
  cs4B.SaveAs("etaSignal_naniVsTnp1_withWeight.png");

  TCanvas cs4C("cs4C","cs4C",1);
  etaSignalMcFromTnP2WW->DrawNormalized("hist");
  etaSignalMcNoTnP2->DrawNormalized("samehist");
  leg3ww->Draw();
  cs4C.SaveAs("etaSignal_naniVsTnp2_withWeight.png");

  TCanvas cs4D("cs4D","cs4D",1);
  etaSignalMcFromTnP3WW->DrawNormalized("hist");
  etaSignalMcNoTnP3->DrawNormalized("samehist");
  leg3ww->Draw();
  cs4D.SaveAs("etaSignal_naniVsTnp3_withWeight.png");

  TCanvas cs4E("cs4E","cs4E",1);
  etaSignalMcFromTnP4WW->DrawNormalized("hist");
  etaSignalMcNoTnP4->DrawNormalized("samehist");
  leg3ww->Draw();
  cs4E.SaveAs("etaSignal_naniVsTnp4_withWeight.png");

  TCanvas cs5("cs5","cs5",1);
  mvaSignalMcNoTnP->DrawNormalized("hist");
  mvaSignalMcFromTnP->DrawNormalized("samehist");
  leg2->Draw();
  cs5.SaveAs("mvaSignal_naniVsTnp.png");

  TCanvas cs6("cs6","cs6",1);
  mvaSignalMcNoTnP->DrawNormalized("hist");
  mvaSignalMcFromTnPWW->DrawNormalized("samehist");
  leg2ww->Draw();
  cs6.SaveAs("mvaSignal_naniVsTnp_withWeight.png");

  TCanvas csw("csw","cs",1);
  ptVsEtaSignalWeights->Draw("colz");
  csw.SaveAs("ptVsEtaSignalWeights.png"); 

  // -----------------------------------------------
  TCanvas caeb("caeb","caeb",1);
  caeb.Divide(3,2);
  caeb.cd(1); 
  mvaSignalEBMc0NoTnP->DrawNormalized("hist");
  mvaSignalEBMc0FromTnP->DrawNormalized("samehist");
  caeb.cd(2); 
  mvaSignalEBMc1NoTnP->DrawNormalized("hist");
  mvaSignalEBMc1FromTnP->DrawNormalized("samehist");
  caeb.cd(3); 
  mvaSignalEBMc2NoTnP->DrawNormalized("hist");
  mvaSignalEBMc2FromTnP->DrawNormalized("samehist");
  caeb.cd(4); 
  mvaSignalEBMc3NoTnP->DrawNormalized("hist");
  mvaSignalEBMc3FromTnP->DrawNormalized("samehist");
  caeb.cd(5); 
  mvaSignalEBMc4NoTnP->DrawNormalized("hist");
  mvaSignalEBMc4FromTnP->DrawNormalized("samehist");
  caeb.SaveAs("mvaSignalEB_naniVsTnp.png");
  //
  TCanvas caee("caee","caee",1);
  caee.Divide(3,2);
  caee.cd(1); 
  mvaSignalEEMc0NoTnP->DrawNormalized("hist");
  mvaSignalEEMc0FromTnP->DrawNormalized("samehist");
  caee.cd(2); 
  mvaSignalEEMc1NoTnP->DrawNormalized("hist");
  mvaSignalEEMc1FromTnP->DrawNormalized("samehist");
  caee.cd(3); 
  mvaSignalEEMc2NoTnP->DrawNormalized("hist");
  mvaSignalEEMc2FromTnP->DrawNormalized("samehist");
  caee.cd(4); 
  mvaSignalEEMc3NoTnP->DrawNormalized("hist");
  mvaSignalEEMc3FromTnP->DrawNormalized("samehist");
  caee.cd(5); 
  mvaSignalEEMc4NoTnP->DrawNormalized("hist");
  mvaSignalEEMc4FromTnP->DrawNormalized("samehist");
  caee.SaveAs("mvaSignalEE_naniVsTnp.png");
  //
  TCanvas caebww("caebww","caebww",1);
  caebww.Divide(3,2);
  caebww.cd(1); 
  mvaSignalEBMc0NoTnP->DrawNormalized("hist");
  mvaSignalEBMc0WWFromTnP->DrawNormalized("samehist");
  caebww.cd(2); 
  mvaSignalEBMc1NoTnP->DrawNormalized("hist");
  mvaSignalEBMc1WWFromTnP->DrawNormalized("samehist");
  caebww.cd(3); 
  mvaSignalEBMc2NoTnP->DrawNormalized("hist");
  mvaSignalEBMc2WWFromTnP->DrawNormalized("samehist");
  caebww.cd(4); 
  mvaSignalEBMc3NoTnP->DrawNormalized("hist");
  mvaSignalEBMc3WWFromTnP->DrawNormalized("samehist");
  caebww.cd(5); 
  mvaSignalEBMc4NoTnP->DrawNormalized("hist");
  mvaSignalEBMc4WWFromTnP->DrawNormalized("samehist");
  caebww.SaveAs("mvaSignalEB_naniVsTnp_withWeights.png");
  //
  TCanvas caeeww("caeeww","caeeww",1);
  caeeww.Divide(3,2);
  caeeww.cd(1); 
  mvaSignalEEMc0NoTnP->DrawNormalized("hist");
  mvaSignalEEMc0WWFromTnP->DrawNormalized("samehist");
  caeeww.cd(2); 
  mvaSignalEEMc1NoTnP->DrawNormalized("hist");
  mvaSignalEEMc1WWFromTnP->DrawNormalized("samehist");
  caeeww.cd(3); 
  mvaSignalEEMc2NoTnP->DrawNormalized("hist");
  mvaSignalEEMc2WWFromTnP->DrawNormalized("samehist");
  caeeww.cd(4); 
  mvaSignalEEMc3NoTnP->DrawNormalized("hist");
  mvaSignalEEMc3WWFromTnP->DrawNormalized("samehist");
  caeeww.cd(5); 
  mvaSignalEEMc4NoTnP->DrawNormalized("hist");
  mvaSignalEEMc4WWFromTnP->DrawNormalized("samehist");
  caeeww.SaveAs("mvaSignalEE_naniVsTnp_withWeights.png");


  // Check weights 2dim
  TH2F *check2dSignal = (TH2F*)ptVsEtaSignalMcFromTnPWW->Clone("check2dSignal");  
  check2dSignal->Sumw2();
  check2dSignal->Divide(ptVsEtaSignalMcNoTnP);

  // ----------------------------------------------  
  TFile myFile("weightFile_tnpVsNani.root","RECREATE");
  ptVsEtaSignalWeights->Write();
  check2dSignal->Write();
  myFile.Close();
}
