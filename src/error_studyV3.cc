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
#include "DM_DataV2.hh"
#include "THStack.h"
#include "TString.h"
#include "DM_StackPlotsV3.hh"


using namespace std;

//Oct29-2013 and before
/*
const float ZJetsNuNu::RSQ_BinArr[] = {0.5, 0.65, 0.8, 1.0, 2.50};
const float ZJetsNuNu::MR_BinArr[] = {200., 400., 600., 800., 3500.};

const float WJetsHTBins::RSQ_BinArr[] = {0.5, 0.65, 0.8, 1.0,  2.50};
const float WJetsHTBins::MR_BinArr[] = {200., 400., 600., 800., 3500.};

const float DY::RSQ_BinArr[] = {0.5, 0.65, 0.8, 1.0, 2.50};
const float DY::MR_BinArr[] = {200., 400., 600., 800., 3500.};

const float TTJets::RSQ_BinArr[] = {0.5, 0.65, 0.8, 1.0, 2.50};
const float TTJets::MR_BinArr[] = {200., 400., 600., 800., 3500.};

const float BaseDM::RSQ_BinArr[] = {0.5, 0.65, 0.8, 1.0, 2.50};
const float BaseDM::MR_BinArr[] = {200., 400., 600., 800., 3500.};
*/

const float BaseDM::RSQ_BinArr[] = {0.5, 0.6, 0.725, 0.85, 1.0, 2.50};
const float BaseDM::MR_BinArr[] = {200., 300., 400., 600., 800., 3500.};

int main(){
  
  CreateStackPlots();
  //CreateRatioPlots();
  //Create2DPlots();
  CreateMetPlots();
  
  return 0;
  
}  






