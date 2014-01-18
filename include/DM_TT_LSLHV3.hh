#ifndef DM_TTJetsV3_HH
#define DM_TTJetsV3_HH 1

#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "DM_BaseV3.hh"
#include "TEfficiency.h"
#include <vector>
#include <math.h>

class TTJets: public BaseDM{
  
public:
  
  TTJets();
  TTJets(int );
  
  ~TTJets();
  
};


#endif
