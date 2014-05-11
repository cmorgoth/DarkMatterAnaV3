#include <iostream>
#include "math.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TString.h"
#include <fstream>
#include "TEfficiency.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include <vector>

int main(){

  gROOT->Reset();
  
  TChain* out = new TChain("outTree");
  out->Add("/media/data/cmorgoth/Data/DMData/DataFinal/HTMHT_Parked/HTMHTParked_ILV_RunB_Tot.root");
  out->Add("/media/data/cmorgoth/Data/DMData/DataFinal/HTMHT_Parked/HTMHTParked_ILV_RunC_Tot.root");
  out->Add("/media/data/cmorgoth/Data/DMData/DataFinal/HTMHT_Parked/HTMHTParked_ILV_RunD_Tot.root");

  double mr[4], rsq[4], Jet_PT[20], Jet_Eta[20], Jet_Phi[20], Npassed_In,
    metX[4], metY[4], metCorrX[4], metCorrY[4];
  double pTHem1, pTHem2, etaHem1, etaHem2, phiHem1, phiHem2;
  int btag, box, N_Jets;
  
  out->SetBranchStatus("*", 0);
  out->SetBranchStatus("MR", 1);
  out->SetBranchStatus("RSQ",1);
  out->SetBranchStatus("nBtag", 1);
  out->SetBranchStatus("BOX_NUM",1);
  out->SetBranchStatus("N_Jets",1);
  out->SetBranchStatus("Jet_PT",1);
  out->SetBranchStatus("Jet_Phi",1);
  out->SetBranchStatus("Jet_Eta",1);
  out->SetBranchStatus("pTHem1",1);
  out->SetBranchStatus("pTHem2",1);
  out->SetBranchStatus("etaHem1",1);
  out->SetBranchStatus("etaHem2",1);
  out->SetBranchStatus("phiHem1",1);
  out->SetBranchStatus("phiHem2",1);
  out->SetBranchStatus("metX",1);
  out->SetBranchStatus("metY",1);
  out->SetBranchStatus("metCorrX",1);
  out->SetBranchStatus("metCorrY",1);
  out->SetBranchAddress("MR", mr);
  out->SetBranchAddress("RSQ", rsq);
  out->SetBranchAddress("nBtag", &btag);
  out->SetBranchAddress("BOX_NUM", &box);
  out->SetBranchAddress("N_Jets", &N_Jets);
  out->SetBranchAddress("Jet_PT", Jet_PT);
  out->SetBranchAddress("Jet_Phi", Jet_Phi);
  out->SetBranchAddress("Jet_Eta", Jet_Eta);
  out->SetBranchAddress("pTHem1", &pTHem1);
  out->SetBranchAddress("pTHem2", &pTHem2);
  out->SetBranchAddress("etaHem1", &etaHem1);
  out->SetBranchAddress("etaHem2", &etaHem2);
  out->SetBranchAddress("phiHem1", &phiHem1);
  out->SetBranchAddress("phiHem2", &phiHem2);
  out->SetBranchAddress("metX", metX);
  out->SetBranchAddress("metY", metY);
  out->SetBranchAddress("metCorrX", metCorrX);
  out->SetBranchAddress("metCorrY", metCorrY);
  
  int N_out = out->GetEntries();
  double N_passed[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
  for(int j = 0; j < N_out; j++){
    out->GetEntry(j);
    double hlt_w = 1.0;
    
    TLorentzVector j1;
    TLorentzVector j2;
    j1.SetPtEtaPhiE(pTHem1, etaHem1, phiHem1, pTHem1*cosh(etaHem1));//Hemisphere
    j2.SetPtEtaPhiE(pTHem2, etaHem2, phiHem2, pTHem2*cosh(etaHem2));//Hemisphere
    
    TLorentzVector j1_b;//Normal Jet
    TLorentzVector j2_b;//Normal Jet
    j1_b.SetPtEtaPhiE(Jet_PT[0], Jet_Eta[0], Jet_Phi[0], Jet_PT[0]*cosh(Jet_Eta[0]));//Leading Jet
    j2_b.SetPtEtaPhiE(Jet_PT[1], Jet_Eta[1], Jet_Phi[1], Jet_PT[1]*cosh(Jet_Eta[1]));//SubLeading Jet
    
    
    double Dphi = j1.DeltaPhi(j2);
    double Dphi_mono = j1_b.DeltaPhi(j2_b);
    double MET = sqrt(metX[2]*metX[2]+metY[2]*metY[2]);
    
    if(mr[2] > 200.0 && rsq[2] > 0.5 && btag == 0 && box == 0 && fabs(Dphi) < 2.5)N_passed[0] += hlt_w;
    if(mr[2] > 200.0 && rsq[2] > 0.5 && btag == 0 && box == 0 && fabs(Dphi) < 2.5 && N_Jets == 2)N_passed[1] += hlt_w;
    if(mr[2] > 200.0 && rsq[2] > 0.5 && btag == 0 && box == 0 && fabs(Dphi) < 2.5 && N_Jets > 2)N_passed[2] += hlt_w;
    if(mr[2] > 200.0 && rsq[2] > 0.5 && btag == 0 && box == 0 && fabs(Dphi) < 2.5 && N_Jets == 2 
       && Jet_PT[0] > 110)N_passed[3] += hlt_w;
    if(mr[2] > 200.0 && rsq[2] > 0.5 && btag == 0 && box == 0 && fabs(Dphi) < 2.5 && N_Jets == 2 
       && Jet_PT[0] > 110 && MET > 400) N_passed[4] += hlt_w;
  }
  
  std::cout << N_passed[0] << " " <<
    N_passed[1] << " " <<
    N_passed[2] <<" " <<
    N_passed[3] <<" " <<
    N_passed[4] << std::endl;
  
  return 0;
}
