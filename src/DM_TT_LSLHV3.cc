#include "DM_TT_LSLHV3.hh"
#include <iostream>
#include <iomanip>

TTJets::TTJets(){
  
};

TTJets::TTJets(int MetIndex ): BaseDM( MetIndex, "TT" ){
  
  T = new TChain("outTree");
  T->Add("/media/data/cmorgoth/Data/DMData/TTJets/BtagCorrMC/TTJetsFullyLept_pu_mu_LooseAndTightBtag_sf_ISR.root");
  T->Add("/media/data/cmorgoth/Data/DMData/TTJets/BtagCorrMC/TTJetsSemiLept_pu_mu_LooseAndTightBtag_sf_ISR.root");
  T->Add("/media/data/cmorgoth/Data/DMData/TTJets/BtagCorrMC/TTJetsHad_pu_mu_LooseAndTightBtag_sf_ISR.root");
  
  effT = new TChain("effTree");
  effT->Add("/media/data/cmorgoth/Data/DMData/TTJets/BtagCorrMC/TTJetsFullyLept_pu_mu_LooseAndTightBtag_sf_ISR.root");
  effT->Add("/media/data/cmorgoth/Data/DMData/TTJets/BtagCorrMC/TTJetsSemiLept_pu_mu_LooseAndTightBtag_sf_ISR.root");
  effT->Add("/media/data/cmorgoth/Data/DMData/TTJets/BtagCorrMC/TTJetsHad_pu_mu_LooseAndTightBtag_sf_ISR.root");
  
};

TTJets::~TTJets(){
  delete T;
};

