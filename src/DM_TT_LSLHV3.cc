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
  //T->Add("/media/data2/Documents/cmorgoth/Data/TTbarBtagCorrected/TightCorrectedFinalNtuples/TTJetsFullyLept_pu_mu_sf_ISR_TightCorr.root");
  //T->Add("/media/data2/Documents/cmorgoth/Data/TTbarBtagCorrected/TightCorrectedFinalNtuples/TTJetsSemiLept_pu_mu_sf_ISR_TightCorr.root");
  //T->Add("/media/data2/Documents/cmorgoth/Data/TTbarBtagCorrected/TightCorrectedFinalNtuples/TTJetsHad_pu_mu_sf_ISR_TightCorr.root");

  effT = new TChain("effTree");
  effT->Add("/media/data/cmorgoth/Data/DMData/TTJets/BtagCorrMC/TTJetsFullyLept_pu_mu_LooseAndTightBtag_sf_ISR.root");
  effT->Add("/media/data/cmorgoth/Data/DMData/TTJets/BtagCorrMC/TTJetsSemiLept_pu_mu_LooseAndTightBtag_sf_ISR.root");
  effT->Add("/media/data/cmorgoth/Data/DMData/TTJets/BtagCorrMC/TTJetsHad_pu_mu_LooseAndTightBtag_sf_ISR.root");
  //effT->Add("/media/data2/Documents/cmorgoth/Data/TTbarBtagCorrected/TightCorrectedFinalNtuples/TTJetsFullyLept_pu_mu_sf_ISR_TightCorr.root");
  //effT->Add("/media/data2/Documents/cmorgoth/Data/TTbarBtagCorrected/TightCorrectedFinalNtuples/TTJetsSemiLept_pu_mu_sf_ISR_TightCorr.root");
  //effT->Add("/media/data2/Documents/cmorgoth/Data/TTbarBtagCorrected/TightCorrectedFinalNtuples/TTJetsHad_pu_mu_sf_ISR_TightCorr.root");
  
};

TTJets::~TTJets(){
  delete T;
};

