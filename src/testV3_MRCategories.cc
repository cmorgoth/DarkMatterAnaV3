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
  //std::ofstream ofs("TEXFiles/VetoBtag_Loose_May_2014_NewTriger_AN.tex", std::ofstream::out);
  //std::ofstream ofs("TEXFiles/OneBtag_Loose_May_2014_NewTriger_AN_4B_Mcut_Z.tex", std::ofstream::out);
  //std::ofstream ofs("TEXFiles/TwoBtag_Tight_May_2014_NewTriger_AN_4B_McutV3.tex", std::ofstreout);
  std::ofstream ofs("TEXFiles/OneLBtag_FixBtag_PFNoPu_Sep2014_INCLUSIVE.tex", std::ofstream::out);
  //std::ofstream ofs("TEXFiles/OneBtag_Loose_May_2014_NewTriger_NoMU_NoMassCut_NB_V3.tex", std::ofstream::out);
  
  //////////////////////////
  ///////Data//////////////
  ////////////////////////                                                                  
  Data* D = new Data("", 2 );
  D->SetBtagCut(bL,bM,bT);
  //D->PrintEvents();

  std::vector<TH1F*> Data = D->Plot_MRCat();
  for(int i = 0; i < 12; i++){
    std::cout << i << " " << (Data[i])->Integral() << std::endl;
  }

  ////////////////////////
  ////////DY/////////////
  //////////////////////

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
  //TT->PrintEvents();
  
  std::vector<TH1F*> TTjets = TT->Plot_MRCat();
  for(int i = 0; i < 12; i++){
  std::cout << i << " " << (TTjets[i])->Integral() << std::endl;
  }
 
  ///////////////////////////////////
  ////////////Z(nunu)+Jets///////////
  //////////////////////////////////
  ZJetsNuNu* Z = new ZJetsNuNu( 2 );
  Z->SetBtagCut(bL,bM,bT);
  //Z->PrintEvents();

  std::vector<TH1F*> Zjets = Z->Plot_MRCat();
  for(int i = 0; i < 12; i++){
  std::cout << i << " " << (Zjets[i])->Integral() << std::endl;
  }
  
  //////////////////////////
  /////////W+Jets//////////
  /////////////////////////
  WJetsHTBins* W = new WJetsHTBins( 2 );
  W->SetBtagCut(bL,bM,bT);
  //W->PrintEvents();

  std::vector<TH1F*> Wjets = W->Plot_MRCat();
  for(int i = 0; i < 12; i++){
  std::cout << i << " " << (Wjets[i])->Integral() << std::endl;
  }
  
  double INT_DY[12];
  double ERR_DY[12];
  double INT_Z[12];
  double ERR_Z[12];
  double INT_W[12];
  double ERR_W[12];
  double INT_TT[12];
  double ERR_TT[12];
  double INT_TOT[12];
  double ERR_TOT[12];
  double INT_Data[12];
  double ERR_Data[12];
  
  for(int i = 0; i < 12; i++){
    INT_DY[i] = DYjets[i]->IntegralAndError(1, DYjets[i]->GetNbinsX(), ERR_DY[i], "");
    INT_Z[i] = Zjets[i]->IntegralAndError(1, Zjets[i]->GetNbinsX(), ERR_Z[i], "");
    INT_W[i] = Wjets[i]->IntegralAndError(1, Wjets[i]->GetNbinsX(), ERR_W[i], "");
    INT_TT[i] = TTjets[i]->IntegralAndError(1, TTjets[i]->GetNbinsX(), ERR_TT[i], "");
    INT_TOT[i] = INT_DY[i] + INT_Z[i] + INT_W[i] + INT_TT[i];
    ERR_TOT[i] = sqrt(ERR_DY[i]*ERR_DY[i] + ERR_Z[i]*ERR_Z[i] + ERR_W[i]*ERR_W[i] + ERR_TT[i]*ERR_TT[i]);
    INT_Data[i] = Data[i]->IntegralAndError(1, Data[i]->GetNbinsX(), ERR_Data[i], "");
  }
  
  for(int k = 0; k < 3; k++){
    ofs << "-------" << k << " Mu---------------\n";

    ofs << "VL & " << INT_Z[4*k+0] << "\\pm" << ERR_Z[4*k+0] << "&"
	<< INT_W[4*k+0] << "\\pm" << ERR_W[4*k+0] << "&"
	<< INT_DY[4*k+0] << "\\pm" << ERR_DY[4*k+0] << "&"
	<< INT_TT[4*k+0] << "\\pm" << ERR_TT[4*k+0] << "&"
	<< INT_TOT[4*k+0] << "\\pm" << ERR_TOT[4*k+0] << "&"
	<< "--" << "\\pm" << "--" << "&"
	<< INT_Data[4*k+0] << "\\pm" << ERR_Data[4*k+0] << "&\n";
    
    ofs << "L & " << INT_Z[4*k+1] << "\\pm" << ERR_Z[4*k+1] << "&"
	<< INT_W[4*k+1] << "\\pm" << ERR_W[4*k+1] << "&"
	<< INT_DY[4*k+1] << "\\pm" << ERR_DY[4*k+1] << "&"
	<< INT_TT[4*k+1] << "\\pm" << ERR_TT[4*k+1] << "&"
	<< INT_TOT[4*k+1] << "\\pm" << ERR_TOT[4*k+1] << "&"
	<< "--" << "\\pm" << "--" << "&"
	<< INT_Data[4*k+1] << "\\pm" << ERR_Data[4*k+1] << "&\n";
    
    ofs << "H & " << INT_Z[4*k+2] << "\\pm" << ERR_Z[4*k+2] << "&"
	<< INT_W[4*k+2] << "\\pm" << ERR_W[4*k+2] << "&"
	<< INT_DY[4*k+2] << "\\pm" << ERR_DY[4*k+2] << "&"
	<< INT_TT[4*k+2] << "\\pm" << ERR_TT[4*k+2] << "&"
	<< INT_TOT[4*k+2] << "\\pm" << ERR_TOT[4*k+2] << "&"
	<< "--" << "\\pm" << "--" << "&"
	<< INT_Data[4*k+2] << "\\pm" << ERR_Data[4*k+2] << "&\n";
    
    ofs << "VL & " << INT_Z[4*k+3] << "\\pm" << ERR_Z[4*k+3] << "&"
	<< INT_W[4*k+3] << "\\pm" << ERR_W[4*k+3] << "&"
	<< INT_DY[4*k+3] << "\\pm" << ERR_DY[4*k+3] << "&"
	<< INT_TT[4*k+3] << "\\pm" << ERR_TT[4*k+3] << "&"
	<< INT_TOT[4*k+3] << "\\pm" << ERR_TOT[4*k+3] << "&"
	<< "--" << "\\pm" << "--" << "&"
	<< INT_Data[4*k+3] << "\\pm" << ERR_Data[4*k+3] << "&\n";
  }
  ofs.close();
  
  //TFile* f1 = new TFile("ROOTFiles/VetoBtag_Loose_May_2014_NewTriger_AN.root", "RECREATE");
  //TFile* f1 = new TFile("ROOTFiles/TwoBtag_Tight_May_2014_NewTriger_AN_4B.root", "RECREATE");
  //TFile* f1 = new TFile("ROOTFiles/OneBtag_Loose_May_2014_NewTriger_AN_4B_Mcut_Z.root", "RECREATE");
  TFile* f1 = new TFile("ROOTFiles/OneLBtag_FixBtag_PFNoPu_Sep2014_INCLUSIVE.root", "RECREATE");
  //TFile* f1 = new TFile("ROOTFiles/OneBtag_Loose_May_2014_NewTriger_NoMU_NoMassCut_NB_V3.root", "RECREATE");
  for(int i = 0; i < 15; i++){
    (DYjets[i])->Write();
    (Zjets[i])->Write();
    (Wjets[i])->Write();
    (TTjets[i])->Write();
    (Data[i])->Write();
  }
  
  

  f1->Close();
 
  return 0;
  
}  






