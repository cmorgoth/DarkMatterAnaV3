#include "DM_TT_LSLHV3.hh"
#include <iostream>
#include <iomanip>

TTJets::TTJets(){
  
};

TTJets::TTJets(int MetIndex ): BaseDM( MetIndex, "TT" ){
  
  T = new TChain("outTree");
  //T->Add("/mnt/hadoop/store/user/cmorgoth/FINALFINAL_8TeV_SM_MC_MAY2014/TTJets/TTJetsFullyLept_May2014_NopfMuonID_TLbtag.root");
  //T->Add("/mnt/hadoop/store/user/cmorgoth/FINALFINAL_8TeV_SM_MC_MAY2014/TTJets/TTJetsSemLept_May2014_NopfMuonID_TLbtag.root");
  //T->Add("/mnt/hadoop/store/user/cmorgoth/FINALFINAL_8TeV_SM_MC_MAY2014/TTJets/TTJetsHam_May2014_NopfMuonID_TLbtag.root");
  
  //Aug2014
  T->Add("/wntmp/cmorgoth/TTJetsFullyLept_PFNoPU_Fixed_Lbtag.root");
  T->Add("/wntmp/cmorgoth/TTJetsSemiLept_PFNoPU_Fixed_Lbtag.root");
  T->Add("/wntmp/cmorgoth/TTJetsHad_PFNoPU_Fixed_Lbtag.root");

  effT = new TChain("effTree");
  //effT->Add("/mnt/hadoop/store/user/cmorgoth/FINALFINAL_8TeV_SM_MC_MAY2014/TTJets/TTJetsFullyLept_May2014_NopfMuonID_TLbtag.root");
  //effT->Add("/mnt/hadoop/store/user/cmorgoth/FINALFINAL_8TeV_SM_MC_MAY2014/TTJets/TTJetsSemLept_May2014_NopfMuonID_TLbtag.root");
  //effT->Add("/mnt/hadoop/store/user/cmorgoth/FINALFINAL_8TeV_SM_MC_MAY2014/TTJets/TTJetsHam_May2014_NopfMuonID_TLbtag.root");
  
  //Aug2014                                                                             
  effT->Add("/wntmp/cmorgoth/TTJetsFullyLept_PFNoPU_Fixed_Lbtag.root");
  effT->Add("/wntmp/cmorgoth/TTJetsSemiLept_PFNoPU_Fixed_Lbtag.root");
  effT->Add("/wntmp/cmorgoth/TTJetsHad_PFNoPU_Fixed_Lbtag.root");
  
};

TTJets::~TTJets(){
  delete T;
};

