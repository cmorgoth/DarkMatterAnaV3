#include "DM_DY_HTBinsV3.hh"
#include <iostream>
#include <iomanip>

DY::DY(){
  
};

DY::DY(int MetIndex ): BaseDM( MetIndex, "DY" ){
  
  T = new TChain("outTree");
  T->Add("/media/data/cmorgoth/Data/DMData/Pt100Data/DY_pt100_sf.root");
  //T->Add("/media/data/cmorgoth/Data/DMData/Pt100Data/Znunu_ptZ100_sf_xsec207p1.root");
  //T->Add("/media/data/cmorgoth/Data/DMData/DYJets/BtagCorrMC/DYJetsHT200To400_PU_MU_LooseBtag_sf.root");
  //T->Add("/media/data/cmorgoth/Data/DMData/DYJets/BtagCorrMC/DYJetsHT400_PU_MU_LooseBtag_sf.root");
  
  effT = new TChain("effTree");
  effT->Add("/media/data/cmorgoth/Data/DMData/Pt100Data/DY_pt100_sf.root");
  //effT->Add("/media/data/cmorgoth/Data/DMData/Pt100Data/Znunu_ptZ100_sf_xsec207p1.root");
  //effT->Add("/media/data/cmorgoth/Data/DMData/DYJets/BtagCorrMC/DYJetsHT200To400_PU_MU_LooseBtag_sf.root");
  //effT->Add("/media/data/cmorgoth/Data/DMData/DYJets/BtagCorrMC/DYJetsHT400_PU_MU_LooseBtag_sf.root");
  
};

DY::~DY(){
  delete T;
};

