#ifndef DM_1DRATIOS_V3
#define DM_1DRATIOS_V3 1

#include <iostream>
#include <sstream>
#include "TCanvas.h"
#include "TLatex.h"
#include "TH1F.h"
#include "TLegend.h"
#include "TPad.h"
#include "TString.h"
#include "DM_BaseV3.hh"
#include "THStack.h"

int RatioPlots( TH1F*, TH1F*, TString, TString, TString, TString );
int RatioPlotsV2( THStack*, TH1F*, TH1F*, TString, TString, TString, TString, TLegend* );

#endif
