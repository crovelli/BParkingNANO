#define testRoc_cxx
#include "testRoc.h"

// ROOT includes 
#include "TVector3.h"
#include <TROOT.h>
#include <TStyle.h>
#include <TH1.h>
#include <TH2.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <iostream>

// Utilities
#include "./FiguresOfMeritEvaluator.cc"

using namespace std;

void testRoc::Loop()
{
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  
  // Histos
  TH1F *mvaSignal = new TH1F("mvaSignal","mvaSignal",100, -10., 10.);
  TH1F *nnSignal  = new TH1F("nnSignal", "nnSignal", 100,  0.,  1.5);
  TH1F *mvaBack   = new TH1F("mvaBack",  "mvaBack",  100, -10., 10.);
  TH1F *nnBack    = new TH1F("nnBack",   "nnBack",   100,  0.,  1.5);
  TH1F *dR        = new TH1F("dR",       "dR",       100,  0.,  10.);
  mvaSignal->Sumw2();
  nnSignal->Sumw2();
  mvaBack->Sumw2();
  nnBack->Sumw2();
  dR->Sumw2(); 

  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    // Look for mc truth e-/e+
    int mcPos = -999;
    int mcEle = -999;
    for (int iGen=0; iGen<nGenPart; iGen++) {
      if ( GenPart_pdgId[iGen]==11 && GenPart_status[iGen]==1 )  mcPos = iGen;
      if ( GenPart_pdgId[iGen]==-11 && GenPart_status[iGen]==1 ) mcEle = iGen;
    } 
    TVector3 mcEleTV3(0,0,0);
    if (mcEle>=0) mcEleTV3.SetPtEtaPhi(GenPart_pt[mcEle], GenPart_eta[mcEle], GenPart_phi[mcEle]);
    TVector3 mcPosTV3(0,0,0);
    if (mcPos>=0) mcPosTV3.SetPtEtaPhi(GenPart_pt[mcPos], GenPart_eta[mcPos], GenPart_phi[mcPos]);

    // Loop over electrons
    int trueEle = -1;
    for (int iEle=0; iEle<nElectron; iEle++) {

      if (fabs(Electron_eta[iEle])>2.4) continue;
      if (Electron_pt[iEle]<0.5)        continue;
      if (Electron_isPF[iEle]==1)       continue;
      if (Electron_isLowPt[iEle]==0)    continue;  

      TVector3 recoEleTV3(0,0,0); 
      recoEleTV3.SetPtEtaPhi(Electron_pt[iEle], Electron_eta[iEle], Electron_phi[iEle]);
      if (mcEle>=0) dR->Fill( mcEleTV3.DeltaR(recoEleTV3) );
      if (mcPos>=0) dR->Fill( mcPosTV3.DeltaR(recoEleTV3) );

      float minDR = 999;
      if (mcEle>=0 && (mcEleTV3.DeltaR(recoEleTV3))<minDR) minDR = mcEleTV3.DeltaR(recoEleTV3);
      if (mcPos>=0 && (mcPosTV3.DeltaR(recoEleTV3))<minDR) minDR = mcPosTV3.DeltaR(recoEleTV3);
      
      if (minDR<0.03) trueEle=1;
      if (minDR>0.1)  trueEle=0;
      if (trueEle==1) { 
	mvaSignal->Fill(Electron_mvaId[iEle]);
	nnSignal->Fill(Electron_nnId[iEle]);
      } else if (trueEle==0) {       
	mvaBack->Fill(Electron_mvaId[iEle]);
	nnBack->Fill(Electron_nnId[iEle]);
      }
    } // Loop over electrons

  } // Loop over entries


  // Cosmetics
  mvaSignal->SetLineWidth(2);
  nnSignal ->SetLineWidth(2);
  mvaBack  ->SetLineWidth(2);
  nnBack   ->SetLineWidth(2);
  mvaSignal->SetLineColor(1);
  nnSignal ->SetLineColor(1);
  mvaBack  ->SetLineColor(2);
  nnBack   ->SetLineColor(2);

  // Plots
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);

  // Outputs
  TLegend *leg;
  leg = new TLegend(0.15,0.55,0.50,0.80);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.05);
  leg->SetFillColor(0);
  leg->AddEntry(mvaSignal, "Signal", "lp");
  leg->AddEntry(mvaBack, "Fakes", "lp");

  TCanvas cmva("cmva","cmva",1);
  mvaSignal->DrawNormalized();
  mvaBack->DrawNormalized("same");
  leg->Draw();
  cmva.SaveAs("outputBDT.png");

  TCanvas cnn("cnn","cnn",1);
  nnSignal->DrawNormalized();
  nnBack->DrawNormalized("same");
  leg->Draw();
  cnn.SaveAs("outputNN.png");

  // FOMs
  FiguresOfMeritEvaluator myRocBDT;
  myRocBDT.addSignal("BDT", mvaSignal);
  myRocBDT.addBackgrounds(mvaBack);
  myRocBDT.setCutDirection(">");
  TGraph *myGraphBDT= myRocBDT.getFOM("BDT",2);
  myGraphBDT->SetTitle("BDT");
  myGraphBDT->SetName("BDT");

  FiguresOfMeritEvaluator myRocNN;
  myRocNN.addSignal("NN", nnSignal);
  myRocNN.addBackgrounds(nnBack);
  myRocNN.setCutDirection(">");
  TGraph *myGraphNN= myRocNN.getFOM("NN",2);
  myGraphNN->SetTitle("NN");
  myGraphNN->SetName("NN");

  TLegend *leg2;
  leg2 = new TLegend(0.15,0.55,0.50,0.80);
  leg2->SetFillStyle(0);
  leg2->SetBorderSize(0);
  leg2->SetTextSize(0.05);
  leg2->SetFillColor(0);
  leg2->AddEntry(myGraphBDT, "BDT", "lp");
  leg2->AddEntry(myGraphNN, "NN", "lp");

  TCanvas croc("roc","",1);
  croc.SetGrid();
  TH2F *myH = new TH2F("myH","",100, 0., 1., 100, 0.,1.);
  myH->GetXaxis()->SetTitle("Mistag Rate");
  myH->GetYaxis()->SetTitle("Efficiency");
  myH->Draw();
  myGraphBDT->SetMarkerColor(4);
  myGraphBDT->SetMarkerStyle(20);
  myGraphBDT->SetMarkerSize(1);
  myGraphBDT->Draw("sameP");
  myGraphNN->SetMarkerColor(2);
  myGraphNN->SetMarkerStyle(21);
  myGraphNN->SetMarkerSize(1);
  myGraphNN->Draw("sameP");
  leg2->Draw();
  croc.SaveAs("roc.png");  
}
