#ifndef DM_BaseV3_HH
#define DM_BaseV3_HH 1

#include <iostream>
#include <map>
#include <vector>
#include <string>
#include <sstream>
#include <math.h>
#include <iomanip>

#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include "TChain.h"
#include "TString.h"
#include "TFile.h"
#include "TEfficiency.h"

class BaseDM{
  
public:
  
  TChain* T;
  TChain* effT;
  
  //static const int MR_Bins = 5;
  //static const int RSQ_Bins = 5;
  static const int MR_Bins = 4;
  static const int RSQ_Bins = 4;
  
  //static const float Lumi = 19.6;//fb^{-1}
  static const float Lumi = 18.836;
  
  //Choose btagIndex Accordingly
  //0->Veto Btag(Loose), 1-> Btag(Loose) >=Cut0, 2-> BtagMed >= cut1, 3->BtagTight == cut2,4->BtagTight >= cut2
  static const int btagIndex = 0;
  
  static const float RSQ_BinArr[RSQ_Bins+1];
  static const float MR_BinArr[MR_Bins+1];
  //float RSQ_BinArr[] = {0.5, 0.6, 0.725, 0.85, 2.50};
  //float MR_BinArr[] = {200., 300., 400., 600., 3500.};
  
  BaseDM();
  BaseDM(int , TString);
  BaseDM(const char* );
  
  ~BaseDM();
  
  virtual TH1F PlotMR_1Box();
  virtual TH1F PlotMR_0Box();
  
  virtual TH1F PlotRSQ_1Box();
  virtual TH1F PlotRSQ_0Box();
  
  virtual TH2F PlotRSQ_vs_MR_0Box();
  virtual TH2F PlotRSQ_vs_MR_1Box();

  virtual bool pfJetPassCSVM(double );
  virtual int pfJetPassCSVM(double*, int);

  virtual std::vector<TH2F*> Plot_2DRazor();
  //virtual std::vector<TH2F*> Plot_2D_Kine();
  virtual std::vector<TH1F*> Plot_1DRazor();
  virtual std::vector<TH1F*> PlotMETmag();
  virtual std::vector<TH1F*> PlotKine();
  virtual std::vector<TH1F*> Plot_MRCat();//MR Categories
  
  virtual bool PrintEvents();
  virtual bool SetStatus();
  virtual bool SetMetStatus();
  virtual bool SetStatusKine();
  
  virtual bool SetBtagCut(int a, int b, int c){nBtagCut[0]=a; nBtagCut[1]=b; nBtagCut[2]=c;};
  virtual double HLTscale(double, double);
  virtual double HLTscaleEle(double, double);
  virtual double GetTriggerEff(double , TString);

private:
  
  //TChain* T;
  //TChain* effT;
  
  TH2F* hlt;
  TH2F* hlt_ele;
  TEfficiency* eff;
  TEfficiency* eff_ele;
  
  TH1F* to_0mu;
  TH1F* to_12mu;
  
  //int metIndex;
  //float MRMin;
  //float RSQMin;
  //TString pName;
    
protected:

  int metIndex;
  float MRMin;
  float RSQMin;
  TString pName;
  
  bool fBtag[5];
  TString BtagBranch;
  int nBtagCut[3];
  
};

#endif
