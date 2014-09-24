#include "DM_ZJetsNuNuV3.hh"
#include <iostream>
#include <iomanip>

ZJetsNuNu::ZJetsNuNu(){
  
};

ZJetsNuNu::ZJetsNuNu(int MetIndex ): BaseDM( MetIndex, "Z" ){
  
  T = new TChain("outTree");
  //T->Add("/mnt/hadoop/store/user/cmorgoth/FINALFINAL_8TeV_SM_MC_MAY2014/ZJets/ZJetsToNuNu_50_HT_100_May2014_NopfMuonID_Lbtag.root");
  //T->Add("/mnt/hadoop/store/user/cmorgoth/FINALFINAL_8TeV_SM_MC_MAY2014/ZJets/ZJetsToNuNu_100_HT_200_May2014_NopfMuonID_Lbtag.root");
  //T->Add("/mnt/hadoop/store/user/cmorgoth/FINALFINAL_8TeV_SM_MC_MAY2014/ZJets/ZJetsToNuNu_200_HT_400_May2014_NopfMuonID_Lbtag.root");
  //T->Add("/mnt/hadoop/store/user/cmorgoth/FINALFINAL_8TeV_SM_MC_MAY2014/ZJets/ZJetsToNuNu_400_HT_Inf_May2014_NopfMuonID_Lbtag.root");
  
  //Aug2014
  T->Add("/wntmp/cmorgoth/ZJetsToNuNu_50_HT_100_PFNoPU_Fixed_Lbtag.root");
  T->Add("/wntmp/cmorgoth/ZJetsToNuNu_100_HT_200_PFNoPU_Fixed_Lbtag.root");
  T->Add("/wntmp/cmorgoth/ZJetsToNuNu_200_HT_400_PFNoPU_Fixed_Lbtag.root");
  T->Add("/wntmp/cmorgoth/ZJetsToNuNu_400_HT_Inf_PFNoPU_Fixed_Lbtag.root");

  effT = new TChain("effTree");
  //effT->Add("/mnt/hadoop/store/user/cmorgoth/FINALFINAL_8TeV_SM_MC_MAY2014/ZJets/ZJetsToNuNu_50_HT_100_May2014_NopfMuonID_Lbtag.root");
  //effT->Add("/mnt/hadoop/store/user/cmorgoth/FINALFINAL_8TeV_SM_MC_MAY2014/ZJets/ZJetsToNuNu_100_HT_200_May2014_NopfMuonID_Lbtag.root");
  //effT->Add("/mnt/hadoop/store/user/cmorgoth/FINALFINAL_8TeV_SM_MC_MAY2014/ZJets/ZJetsToNuNu_200_HT_400_May2014_NopfMuonID_Lbtag.root");
  //effT->Add("/mnt/hadoop/store/user/cmorgoth/FINALFINAL_8TeV_SM_MC_MAY2014/ZJets/ZJetsToNuNu_400_HT_Inf_May2014_NopfMuonID_Lbtag.root");

  //Aug2014                                                                     
  effT->Add("/wntmp/cmorgoth/ZJetsToNuNu_50_HT_100_PFNoPU_Fixed_Lbtag.root");
  effT->Add("/wntmp/cmorgoth/ZJetsToNuNu_100_HT_200_PFNoPU_Fixed_Lbtag.root");
  effT->Add("/wntmp/cmorgoth/ZJetsToNuNu_200_HT_400_PFNoPU_Fixed_Lbtag.root");
  effT->Add("/wntmp/cmorgoth/ZJetsToNuNu_400_HT_Inf_PFNoPU_Fixed_Lbtag.root");
  
};

ZJetsNuNu::~ZJetsNuNu(){
  delete T;
};

