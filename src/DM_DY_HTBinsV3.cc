#include "DM_DY_HTBinsV3.hh"
#include <iostream>
#include <iomanip>

DY::DY(){
  
};

DY::DY(int MetIndex ): BaseDM( MetIndex, "DY" ){
  
  T = new TChain("outTree");
  //T->Add("/mnt/hadoop/store/user/cmorgoth/FINALFINAL_8TeV_SM_MC_MAY2014/DYJets/DYJetsHT200To400_JEC_FIX_May2014_NopfMuonID_Lbtag.root");
  //T->Add("/mnt/hadoop/store/user/cmorgoth/FINALFINAL_8TeV_SM_MC_MAY2014/DYJets/DYJetsHT400_JEC_FIX_May2014_NopfMuonID_Lbtag.root");
  
  //Aug2014
  T->Add("/wntmp/cmorgoth/DYJetsHT200To400_PFNoPU_Fixed_Lbtag.root");
  T->Add("/wntmp/cmorgoth/DYJetsHT400_PFNoPU_Fixed_Lbtag.root");
  
  effT = new TChain("effTree");
  //effT->Add("/mnt/hadoop/store/user/cmorgoth/FINALFINAL_8TeV_SM_MC_MAY2014/DYJets/DYJetsHT200To400_JEC_FIX_May2014_NopfMuonID_Lbtag.root");
  //effT->Add("/mnt/hadoop/store/user/cmorgoth/FINALFINAL_8TeV_SM_MC_MAY2014/DYJets/DYJetsHT400_JEC_FIX_May2014_NopfMuonID_Lbtag.root");
  
  //Aug2014
  effT->Add("/wntmp/cmorgoth/DYJetsHT200To400_PFNoPU_Fixed_Lbtag.root");
  effT->Add("/wntmp/cmorgoth/DYJetsHT400_PFNoPU_Fixed_Lbtag.root");
  
};

DY::~DY(){
  delete T;
};

