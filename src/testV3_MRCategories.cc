#include <iostream>
#include <sstream>
#include <fstream>
#include "TCanvas.h"
#include "TLatex.h"
#include "TLegend.h"
#include "DM_WJetsHTBinsV3.hh"
#include "DM_ZJetsNuNuV3.hh"
#include "DM_DY_HTBinsV3.hh"
#include "DM_TT_LSLHV3.hh"
#include "DM_DataV3.hh"
#include "DM_METPlotsV3.hh"
#include "DM_BaseV3.hh"
#include "THStack.h"
#include "TString.h"
#include "DM_StackPlotsV3.hh"
#include "DM_1DRatioV3.hh"

using namespace std;

//5x5 v2
//const float BaseDM::RSQ_BinArr[] = {0.5, 0.6, 0.725, 0.85, 1.1, 2.50};                       
//const float BaseDM::MR_BinArr[] = {200., 300., 400., 600., 900., 3500.}; 
//4x4
const float BaseDM::RSQ_BinArr[] = {0.5, 0.6, 0.725, 0.85, 2.50};
const float BaseDM::MR_BinArr[] = {200., 300., 400., 600., 3500.};

int main(){
  
  int bL, bM, bT;
  bL = bM = 1;
  bT = 1;

   
  TCanvas* C1 = new TCanvas("C1", "C1", 1024, 1024);
  C1->cd();
  std::ofstream ofs("TEXFiles/OneTightBtagCorr.tex", std::ofstream::out);
  
  
  DY* dy = new DY( 2 );
  dy->SetBtagCut(bL,bM,bT);
  dy->PrintEvents();
  std::vector<TH1F*> DYjets = dy->Plot_MRCat();
  for(int i = 0; i < 12; i++){
    std::cout << i << " " << (DYjets[i])->Integral() << std::endl;
  }
  
  ///////////////////////////////
  ////////// tt + jets//////////
  //////////////////////////////
  
  TTJets* TT = new TTJets(2);
  TT->SetBtagCut(bL,bM,bT);
  TT->PrintEvents();
  
  std::vector<TH1F*> TTjets = TT->Plot_MRCat();
  for(int i = 0; i < 12; i++){
  std::cout << i << " " << (TTjets[i])->Integral() << std::endl;
  }
 
  ///////////////////////////////////
  ////////////Z(nunu)+Jets///////////
  //////////////////////////////////
  ZJetsNuNu* Z = new ZJetsNuNu( 2 );
  Z->SetBtagCut(bL,bM,bT);
  Z->PrintEvents();

  std::vector<TH1F*> Zjets = Z->Plot_MRCat();
  for(int i = 0; i < 12; i++){
  std::cout << i << " " << (Zjets[i])->Integral() << std::endl;
  }
  
  //////////////////////////
  /////////W+Jets//////////
  /////////////////////////
  WJetsHTBins* W = new WJetsHTBins( 2 );
  W->SetBtagCut(bL,bM,bT);
  W->PrintEvents();

  std::vector<TH1F*> Wjets = W->Plot_MRCat();
  for(int i = 0; i < 12; i++){
  std::cout << i << " " << (Wjets[i])->Integral() << std::endl;
  }

  //////////////////////////
  ///////Data//////////////
  ////////////////////////
  Data* D = new Data("", 2 );
  D->SetBtagCut(bL,bM,bT);
  D->PrintEvents();
  
  std::vector<TH1F*> Data = D->Plot_MRCat();
  for(int i = 0; i < 12; i++){
    std::cout << i << " " << (Data[i])->Integral() << std::endl;
  }
  
  
  TFile* f1 = new TFile("ROOTFiles/OneTightBtagCorr_MRcategories.root", "RECREATE");
  for(int i = 0; i < 12; i++){
    (DYjets[i])->Write();
    (Zjets[i])->Write();
    (Wjets[i])->Write();
    (TTjets[i])->Write();
    (Data[i])->Write();
  }
  
  

  f1->Close();
 
  return 0;
  
}  






