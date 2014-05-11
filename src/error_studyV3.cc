#include <iostream>
#include <sstream>
#include "TCanvas.h"
#include "TLatex.h"
#include "TLegend.h"
#include "DM_WJetsHTBinsV3.hh"
#include "DM_ZJetsNuNuV3.hh"
#include "DM_DY_HTBinsV3.hh"
#include "DM_TT_LSLHV3.hh"
#include "DM_METPlotsV3.hh"
#include "DM_DataV3.hh"
#include "THStack.h"
#include "TString.h"
#include "DM_StackPlotsV3.hh"
#include "DM_KinePlotsV3.hh"
#include "DM_RatioPlotsV3.hh"

using namespace std;

//5x5 v2
//const float BaseDM::RSQ_BinArr[] = {0.5, 0.6, 0.725, 0.85, 1.1, 2.50};
//const float BaseDM::MR_BinArr[] = {200., 300., 400., 600., 900., 3500.};
//4x4
const float BaseDM::RSQ_BinArr[] = {0.5, 0.6, 0.725, 0.85, 2.50};
const float BaseDM::MR_BinArr[] = {200., 300., 400., 600., 3500.};

int main(){
  
  //CreateStackPlots();
  //CreateRatioPlots();
  //Create2DPlots();
  //CreateMetPlots();
  CreateKinePlots();
  
  return 0;
  
}  






