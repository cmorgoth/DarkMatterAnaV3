#ifndef DM_DATA_V3_HH
#define DM_DATA_V3_HH 1

#include "TH1F.h"
#include "DM_BaseV3.hh"
#include "TH2F.h"
#include "TTree.h"
#include "TFile.h"
#include "TChain.h"
#include <math.h>

class Data: public BaseDM{
  
public:
  
  //static const int btagIndex = 0;//0->Veto Btag(Loose), 1-> Btag(Loose) >=1, 2-> BtagTight >=1
  
  Data();
  Data(const char* FileName, int MetIndex);
  //Data(char const*, int );
  
  ~Data();

  
  TH1F PlotMR_1Box();
  TH1F PlotMR_0Box();
  
  TH1F PlotRSQ_1Box();
  TH1F PlotRSQ_0Box();
  
  TH2F PlotRSQ_vs_MR_0Box();
  TH2F PlotRSQ_vs_MR_1Box();
  
  bool pfJetPassCSVM(double );
  int pfJetPassCSVM(double*, int);
  
  std::vector<TH2F*> Plot_2DRazor();
  std::vector<TH1F*> Plot_1DRazor();
  std::vector<TH1F*> PlotMETmag();
  std::vector<TH1F*> PlotKine();
  std::vector<TH1F*> Plot_MRCat();//MR Categories

  bool PrintEvents();
  bool SetStatus();
  bool SetStatusKine();
  bool SetMetStatus();
    
};


#endif
