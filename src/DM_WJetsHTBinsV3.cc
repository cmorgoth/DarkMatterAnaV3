#include "DM_WJetsHTBinsV3.hh"
#include <iostream>
#include <iomanip>

WJetsHTBins::WJetsHTBins(){
  
};

WJetsHTBins::WJetsHTBins(int MetIndex ): BaseDM( MetIndex, "W" ){
  
  T = new TChain("outTree");
  //T->Add("/mnt/hadoop/store/user/cmorgoth/FINALFINAL_8TeV_SM_MC_MAY2014/WJets/WJets_150_HT_200_May2014_NopfMuonID_Lbtag.root");
  //T->Add("/mnt/hadoop/store/user/cmorgoth/FINALFINAL_8TeV_SM_MC_MAY2014/WJets/WJets_200_HT_250_May2014_NopfMuonID_Lbtag.root");
  //T->Add("/mnt/hadoop/store/user/cmorgoth/FINALFINAL_8TeV_SM_MC_MAY2014/WJets/WJets_250_HT_300_May2014_NopfMuonID_Lbtag.root");
  //T->Add("/mnt/hadoop/store/user/cmorgoth/FINALFINAL_8TeV_SM_MC_MAY2014/WJets/WJets_300_HT_400_May2014_NopfMuonID_Lbtag.root");
  //T->Add("/mnt/hadoop/store/user/cmorgoth/FINALFINAL_8TeV_SM_MC_MAY2014/WJets/WJets_400_HT_Inf_May2014_NopfMuonID_Lbtag.root");
  
  //Aug2014
  T->Add("/wntmp/cmorgoth/WJetsToLNu_150_HT_200_PFNoPU_Fixed_Lbtag.root");
  T->Add("/wntmp/cmorgoth/WJetsToLNu_200_HT_250_PFNoPU_Fixed_Lbtag.root");
  T->Add("/wntmp/cmorgoth/WJetsToLNu_250_HT_300_PFNoPU_Fixed_Lbtag.root");
  T->Add("/wntmp/cmorgoth/WJetsToLNu_300_HT_400_PFNoPU_Fixed_Lbtag.root");
  T->Add("/wntmp/cmorgoth/WJetsToLNu_400_HT_Inf_PFNoPU_Fixed_Lbtag.root");
  
  effT = new TChain("effTree");
  //effT->Add("/mnt/hadoop/store/user/cmorgoth/FINALFINAL_8TeV_SM_MC_MAY2014/WJets/WJets_150_HT_200_May2014_NopfMuonID_Lbtag.root");
  //effT->Add("/mnt/hadoop/store/user/cmorgoth/FINALFINAL_8TeV_SM_MC_MAY2014/WJets/WJets_200_HT_250_May2014_NopfMuonID_Lbtag.root");
  //effT->Add("/mnt/hadoop/store/user/cmorgoth/FINALFINAL_8TeV_SM_MC_MAY2014/WJets/WJets_250_HT_300_May2014_NopfMuonID_Lbtag.root");
  //effT->Add("/mnt/hadoop/store/user/cmorgoth/FINALFINAL_8TeV_SM_MC_MAY2014/WJets/WJets_300_HT_400_May2014_NopfMuonID_Lbtag.root");
  //effT->Add("/mnt/hadoop/store/user/cmorgoth/FINALFINAL_8TeV_SM_MC_MAY2014/WJets/WJets_400_HT_Inf_May2014_NopfMuonID_Lbtag.root");

  
  //Aug2014                                                                     
  effT->Add("/wntmp/cmorgoth/WJetsToLNu_150_HT_200_PFNoPU_Fixed_Lbtag.root");
  effT->Add("/wntmp/cmorgoth/WJetsToLNu_200_HT_250_PFNoPU_Fixed_Lbtag.root");
  effT->Add("/wntmp/cmorgoth/WJetsToLNu_250_HT_300_PFNoPU_Fixed_Lbtag.root");
  effT->Add("/wntmp/cmorgoth/WJetsToLNu_300_HT_400_PFNoPU_Fixed_Lbtag.root");
  effT->Add("/wntmp/cmorgoth/WJetsToLNu_400_HT_Inf_PFNoPU_Fixed_Lbtag.root");
  
};

WJetsHTBins::~WJetsHTBins(){
  delete T;
};

