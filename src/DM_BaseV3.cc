#include "DM_BaseV3.hh"
#include <iostream>
#include "TLorentzVector.h"

BaseDM::BaseDM(){ };

BaseDM::BaseDM(int metIndex, TString pName){
  
  MRMin = 200.0;//nominal 150.0
  RSQMin = 0.5;//nominal 0.5
  this->metIndex = metIndex;
  this->pName = pName;
  
  if(pName != "Data"){
    TFile* file = new TFile("/mnt/hadoop/store/user/cmorgoth/TriggerTurnOn_DM_Analysis_8TeV/hlt_eff_HTMHT_12mu_BOX_0bTag_Final.root");
    eff = (TEfficiency*)file->Get("Eff2d");
    
    TFile* file1 = new TFile("/mnt/hadoop/store/user/cmorgoth/TriggerTurnOn_DM_Analysis_8TeV/hlt_eff_HTMHT_0mu_BOX_0bTag_Final.root");
    eff_ele = (TEfficiency*)file1->Get("Eff2d");
    
    TFile* file2 = new TFile("/mnt/hadoop/store/user/cmorgoth/TriggerTurnOn_DM_Analysis_8TeV/TEST_0MU.root");
    to_0mu = (TH1F*)file2->Get("MRn1_0p5");
    //to_0mu = (TH1F*)file2->Get("R2n1");
    TFile* file3 = new TFile("/mnt/hadoop/store/user/cmorgoth/TriggerTurnOn_DM_Analysis_8TeV/TEST_OneTwoMU.root");
    to_12mu = (TH1F*)file3->Get("MRn1_0p5");
    //to_12mu = (TH1F*)file3->Get("R2n1"); 
  }else{
    std::cout << "======== IS DATA!==========" << std::endl;
  }
  //////////////////////////////////////////////////////////////// 
  ////////// Including Btag capability///////////////////////////                               
  ///////////////////////////////////////////////////////////////
  
  std::cout << "====Btag Index====> " << btagIndex << std::endl;
  if( btagIndex == 0 || btagIndex == 1 ){
    this->BtagBranch = "nBtag";
  }else{
    this->BtagBranch = "nBtagTight";
  }
  
  std::cout << "----------------Branch Name:  " << BtagBranch << std::endl;
  std::cout << "------------------BranchIndex: " << this->btagIndex << std::endl;

};


BaseDM::BaseDM(const char* FileName ){
  
  //F = new TFile(FileName);
  //T0 = (TTree*)F->Get("outTree");
  
};

BaseDM::~BaseDM(){
  delete T;
};


bool BaseDM::PrintEvents(){

  double NtotGen = 0, Nt_PV = 0, Nt_2J = 0, Nt_0b = 0, Nt_LepVeto = 0, N_In, N_PV, N_2J, N_0b, N_LepVeto;

  std::cout.precision(11);
  effT->SetBranchAddress("Npassed_In", &N_In);
  effT->SetBranchAddress("Npassed_PV", &N_PV);
  effT->SetBranchAddress("Npassed_2Jet", &N_2J);
  effT->SetBranchAddress("Npassed_0btag", &N_0b);
  effT->SetBranchAddress("Npassed_LepVeto", &N_LepVeto);
  
  
  for (int i = 0; i < effT->GetEntries(); i++){
    effT->GetEntry(i);
    NtotGen += N_In;
    Nt_PV += N_PV;
    Nt_2J += N_2J;
    Nt_0b += N_0b;
    Nt_LepVeto += N_LepVeto;
  }
  
  std::cout << "========================================" << std::endl;
  std::cout << "================ " << this->pName << " ==============" << std::endl;
  std::cout << "========================================" << std::endl;
  
  std::cout << this->pName << "Nt_In: " << NtotGen << std::endl;
  std::cout << this->pName << "Nt_PV: " << Nt_PV << std::endl;
  std::cout << this->pName << "Nt_2J: " << Nt_2J << std::endl;
  std::cout << this->pName << "Nt_0b: " << Nt_0b << std::endl;
  std::cout << this->pName << "Nt_LepVeto: " << Nt_LepVeto << std::endl;
  
  //After Selection cuts
  double RSQ[4], MR[4];
  int BOX, nBtag;
  
  double Nt_MR_RSQ_cut0BTag = 0.0, Nt_2muBox0BTag = 0.0,  Nt_1muBox0BTag = 0.0, Nt_0muBox0BTag = 0.0,\
    Nt_MR_RSQ_cut = 0.0, Nt_2muBox = 0.0,  Nt_1muBox = 0.0, Nt_0muBox = 0.0, Nt_Btags = 0.0, sf_w;

  SetStatus();
  T->SetBranchAddress("RSQ", RSQ);
  T->SetBranchAddress("MR", MR);
  T->SetBranchAddress("BOX_NUM", &BOX);
  T->SetBranchAddress("nBtag", &nBtag);
  T->SetBranchAddress("sf_w", &sf_w);
  
  for(int i = 0; i < T->GetEntries(); i++){
    T->GetEntry(i);
    if( RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin ){
      Nt_MR_RSQ_cut += sf_w;
      
      if( BOX == 0 )Nt_0muBox += sf_w;
      if( BOX == 1 )Nt_1muBox += sf_w;
      if( BOX == 2 )Nt_2muBox += sf_w;
      
      if( nBtag == 0 ){
        
        Nt_MR_RSQ_cut0BTag += sf_w;
        
        if( BOX == 0 )Nt_0muBox0BTag += sf_w;
        if( BOX == 1 )Nt_1muBox0BTag += sf_w;
        if( BOX == 2 )Nt_2muBox0BTag += sf_w;
        
      }
      
    }
    
  }
  T->SetBranchStatus("*", 0);//Disable all branches, trying to gain some performance

  
  std::cout << this->pName << "weighted Nt_MR_RSQ_cut: " << Nt_MR_RSQ_cut << std::endl;
  std::cout << this->pName << "weighted Nt_0muBox: " << Nt_0muBox << std::endl;
  std::cout << this->pName << "weighted Nt_1muBox: " << Nt_1muBox << std::endl;
  std::cout << this->pName << "weighted Nt_2muBox: " << Nt_2muBox << std::endl;

  std::cout << this->pName << "weighted Nt_MR_RSQ_cut0Btag: " << Nt_MR_RSQ_cut0BTag << std::endl;
  std::cout << this->pName << "weighted Nt_0muBox0Btag: " << Nt_0muBox0BTag << std::endl;
  std::cout << this->pName << "weighted Nt_1muBox0Btag: " << Nt_1muBox0BTag << std::endl;
  std::cout << this->pName << "weighted Nt_2muBox0Btag: " << Nt_2muBox0BTag << std::endl;

  std::cout << this->pName << "Btag EVENTS: " << Nt_MR_RSQ_cut - Nt_MR_RSQ_cut0BTag << std::endl;
  std::cout << "========================================" << std::endl;
  std::cout << "========================================" << "\n\n" << std::endl;
  
  return true;
  
  
};

TH1F BaseDM::PlotMR_1Box(){
  
  double RSQ[4], MR[4], CSV[30];
  double hltWeight = .0, sf_w;
  int BOX, N_Jets, nBtag[2];

  TH1F* MR1 = new TH1F("MR1", "MR1BOX", MR_Bins, MR_BinArr);

  SetStatus();
  T->SetBranchAddress("RSQ", RSQ);
  T->SetBranchAddress("MR", MR);
  T->SetBranchAddress("BOX_NUM", &BOX);
  T->SetBranchAddress("nBtag", &nBtag[0]);
  T->SetBranchAddress("nBtagTight", &nBtag[1]);
  T->SetBranchAddress("N_Jets", &N_Jets);
  T->SetBranchAddress("CSV", CSV);
  T->SetBranchAddress("sf_w", &sf_w);
    
  for(int i = 0; i < T->GetEntries(); i++){
    T->GetEntry(i);
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= nBtagCut[0]);
    fBtag[3] = (nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    fBtag[4] = (nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1]);

    if( BOX==1 && RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin && fBtag[btagIndex] ){
      hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
      if( hltWeight == 0.0 )hltWeight = 1.0;
      MR1->Fill(MR[metIndex], sf_w*hltWeight);
    }
  }
  T->SetBranchStatus("*", 0);
    
  return *MR1;
};

TH1F BaseDM::PlotMR_0Box(){
  
  double RSQ[4], MR[4], CSV[30];
  double hltWeight = .0, sf_w;
  int BOX, N_Jets, nBtag[2];

  TH1F* MR0 = new TH1F("MR0", "MR0BOX", MR_Bins, MR_BinArr);

  SetStatus();
  T->SetBranchAddress("RSQ", RSQ);
  T->SetBranchAddress("MR", MR);
  T->SetBranchAddress("BOX_NUM", &BOX);
  T->SetBranchAddress("nBtag", &nBtag[0]);
  T->SetBranchAddress("nBtagTight", &nBtag[1]);
  T->SetBranchAddress("N_Jets", &N_Jets);
  T->SetBranchAddress("CSV", CSV);
  T->SetBranchAddress("sf_w", &sf_w);
  
  for(int i = 0; i < T->GetEntries(); i++){
    T->GetEntry(i);
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= nBtagCut[0]);
    fBtag[3] = (nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1]);
    
    if( BOX == 0 && RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin && fBtag[btagIndex]){
      hltWeight = HLTscaleEle( MR[metIndex], RSQ[metIndex]);
      if( hltWeight == 0.0 )hltWeight = 1.0;
      MR0->Fill(MR[metIndex], sf_w*hltWeight);
    }
  }
  T->SetBranchStatus("*", 0);
    
  return *MR0;
};


TH1F  BaseDM::PlotRSQ_1Box(){
  
  double RSQ[4], MR[4], CSV[30];
  double hltWeight = 0.0, sf_w;
  int BOX, N_Jets, nBtag[2];

  TH1F* RSQ1 = new TH1F("RSQ1", "RSQ1BOX", RSQ_Bins, RSQ_BinArr);
  
  SetStatus();
  T->SetBranchAddress("RSQ", RSQ);
  T->SetBranchAddress("MR", MR);
  T->SetBranchAddress("BOX_NUM", &BOX);
  T->SetBranchAddress("nBtag", &nBtag[0]);
  T->SetBranchAddress("nBtagTight", &nBtag[1]);
  T->SetBranchAddress("N_Jets", &N_Jets);
  T->SetBranchAddress("CSV", CSV);
  T->SetBranchAddress("sf_w", &sf_w);
  
  for(int i = 0; i < T->GetEntries(); i++){
    T->GetEntry(i);
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= nBtagCut[0]);
    fBtag[3] = (nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1]);
    
    if( BOX == 1 && RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin && fBtag[btagIndex]){
      hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
      if( hltWeight == 0.0 )hltWeight = 1.0;
      RSQ1->Fill(RSQ[metIndex], sf_w*hltWeight);
    }
  }
  T->SetBranchStatus("*", 0);
    
  return *RSQ1;
    
};


TH1F  BaseDM::PlotRSQ_0Box(){
  
  double RSQ[4], MR[4], CSV[30];
  double hltWeight = 0.0, sf_w;
  int BOX, N_Jets, nBtag[2];
  TH1F* RSQ0 = new TH1F("RSQ0", "RSQ0BOX", RSQ_Bins, RSQ_BinArr);
  
  SetStatus();
  T->SetBranchAddress("RSQ", RSQ);
  T->SetBranchAddress("MR", MR);
  T->SetBranchAddress("BOX_NUM", &BOX);
  T->SetBranchAddress("nBtag", &nBtag[0]);
  T->SetBranchAddress("nBtagTight", &nBtag[1]);
  T->SetBranchAddress("N_Jets", &N_Jets);
  T->SetBranchAddress("CSV", CSV);
  T->SetBranchAddress("sf_w", &sf_w);
  
  for(int i = 0; i < T->GetEntries(); i++){
    T->GetEntry(i);
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= nBtagCut[0]);
    fBtag[3] = (nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1]);

    if( BOX == 0 && RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin && fBtag[btagIndex]){
      hltWeight = HLTscaleEle( MR[metIndex], RSQ[metIndex]);
      if( hltWeight == 0.0 )hltWeight = 1.0;
      RSQ0->Fill(RSQ[metIndex], sf_w*hltWeight);
    }
  }
  T->SetBranchStatus("*", 0);

  return *RSQ0;
  
};

TH2F BaseDM::PlotRSQ_vs_MR_0Box(){
  
  double RSQ[4], MR[4], CSV[30];
  double hltWeight = 0.0, sf_w;
  int BOX, N_Jets, nBtag[2];
  
  TH2F* RSQ_MR_0BOX = new TH2F("RSQ_MR_0BOX_W", "RSQ_VS_MR_0BOX_W", MR_Bins, MR_BinArr, RSQ_Bins, RSQ_BinArr);
  RSQ_MR_0BOX->Sumw2();
  
  SetStatus();
  T->SetBranchAddress("RSQ", RSQ);
  T->SetBranchAddress("MR", MR);
  T->SetBranchAddress("BOX_NUM", &BOX);
  T->SetBranchAddress("nBtag", &nBtag[0]);
  T->SetBranchAddress("nBtagTight", &nBtag[1]);
  T->SetBranchAddress("N_Jets", &N_Jets);
  T->SetBranchAddress("CSV", CSV);
  T->SetBranchAddress("sf_w", &sf_w);
  
  for(int i = 0; i < T->GetEntries(); i++){
    T->GetEntry(i);
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= nBtagCut[0]);
    fBtag[3] = (nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1]);
    
    if( BOX == 0 && RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin && fBtag[btagIndex]){
      hltWeight = HLTscaleEle( MR[metIndex], RSQ[metIndex]);
      if( hltWeight == 0.0 )hltWeight = 1.0;
      RSQ_MR_0BOX->Fill(MR[metIndex], RSQ[metIndex], sf_w*hltWeight);
    }
  }
  T->SetBranchStatus("*", 0);

  return *RSQ_MR_0BOX;
  
};

TH2F BaseDM::PlotRSQ_vs_MR_1Box(){
  
  double RSQ[4], MR[4], CSV[30];
  double hltWeight = 0.0, sf_w;
  int BOX, N_Jets, nBtag[2];
  
  TH2F* RSQ_MR_1BOX = new TH2F("RSQ_MR_1BOX_W", "RSQ_VS_MR_1BOX_W", MR_Bins, MR_BinArr, RSQ_Bins, RSQ_BinArr);
  RSQ_MR_1BOX->Sumw2();
  
  SetStatus();
  T->SetBranchAddress("RSQ", RSQ);
  T->SetBranchAddress("MR", MR);
  T->SetBranchAddress("BOX_NUM", &BOX);
  T->SetBranchAddress("nBtag", &nBtag[0]);
  T->SetBranchAddress("nBtagTight", &nBtag[1]);
  T->SetBranchAddress("N_Jets", &N_Jets);
  T->SetBranchAddress("CSV", CSV);
  T->SetBranchAddress("sf_w", &sf_w);
  
  for(int i = 0; i < T->GetEntries(); i++){
    T->GetEntry(i);
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= nBtagCut[0]);
    fBtag[3] = (nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1]);
    
    if( BOX == 1 && RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin && fBtag[btagIndex]){
      hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
      if( hltWeight == 0.0 )hltWeight = 1.0;
      RSQ_MR_1BOX->Fill(MR[metIndex], RSQ[metIndex], sf_w*hltWeight);
    }
  }
  T->SetBranchStatus("*", 0);
  
  return *RSQ_MR_1BOX;
  
};


std::vector<TH1F*> BaseDM::PlotMETmag(){
  double metX[4], metcorrX[4], metY[4], metcorrY[4], ht, RSQ[4], MR[4], run/*, ls, evNum*/;
  double mht[3], CSV[30], pu_w, mu_w, sf_w, ISR_w;
  int BOX, nBtag[3], N_Jets;
  double Jet_PT[20], Jet_Eta[20], Jet_Phi[20];

  double Mu_Px[2], Mu_Py[2], Mu_Pz[2], Mu_E[2];
  double pTHem1, pTHem2, etaHem1, etaHem2, phiHem1, phiHem2;
 
  std::vector< TH1F* > metvec;
  TH1F* MET[12];
  TString name;
  for(int l = 0; l < 3; l++ ){
    for( int m = 0; m < 3; m++ ){
      name = TString(Form("BaseDM_METmag_Box%d_plotType%d",l,m));
      MET[3*l + m] = new TH1F( name, name, 20, 0, 1500 );
    }
  }
  
  MET[9] = new TH1F( "NJETS0_W", "NJETS 0 BOX", 9, 1, 10);
  MET[10] = new TH1F( "NJETS1_W", "NJETS 1 BOX", 9, 1, 10);
  MET[11] = new TH1F( "NJETS2_W", "NJETS 2 BOX", 9, 1, 10);
  
  SetMetStatus();
  T->SetBranchStatus("N_Jets",1);
  T->SetBranchStatus("Jet_PT",1);
  T->SetBranchStatus("Jet_Phi",1);
  T->SetBranchStatus("Jet_Eta",1);
  
  T->SetBranchAddress("pu_w", &pu_w);
  T->SetBranchAddress("mu_w", &mu_w);
  T->SetBranchAddress("RSQ", RSQ);
  T->SetBranchAddress("MR", MR);
  T->SetBranchAddress("BOX_NUM", &BOX);
  //T->SetBranchAddress("nBtag", &nBtag[0]);
  T->SetBranchAddress("nBtagLCorr", &nBtag[0]);
  if(this->pName == "TT_LL"){
    T->SetBranchAddress("nBtagTCorr", &nBtag[2]);
  }else{
    T->SetBranchAddress("nBtagTight", &nBtag[2]);
  }
  T->SetBranchAddress("nBtagMed", &nBtag[1]);
  
  T->SetBranchAddress("N_Jets", &N_Jets);
  T->SetBranchAddress("CSV", CSV);
  T->SetBranchAddress("ht", &ht);
  T->SetBranchAddress("mht", &mht[0]);
  T->SetBranchAddress("metX", metX);
  T->SetBranchAddress("metCorrX", metcorrX);
  T->SetBranchAddress("metY", metY);
  T->SetBranchAddress("metCorrY", metcorrY);
  T->SetBranchAddress("N_Jets", &N_Jets);
  T->SetBranchAddress("Jet_PT", Jet_PT);
  T->SetBranchAddress("Jet_Eta", Jet_Eta);
  T->SetBranchAddress("Jet_Phi", Jet_Phi);
  
  T->SetBranchAddress("pTHem1", &pTHem1);
  T->SetBranchAddress("pTHem2", &pTHem2);
  T->SetBranchAddress("etaHem1", &etaHem1);
  T->SetBranchAddress("etaHem2", &etaHem2);
  T->SetBranchAddress("phiHem1", &phiHem1);
  T->SetBranchAddress("phiHem2", &phiHem2);
  T->SetBranchAddress("Mu_Px", Mu_Px);
  T->SetBranchAddress("Mu_Py", Mu_Py);
  T->SetBranchAddress("Mu_Pz", Mu_Pz);
  T->SetBranchAddress("Mu_E", Mu_E);
  T->SetBranchAddress("ISR_w", &sf_w);
  T->SetBranchAddress("ISR", &ISR_w);

  float metmag = .0;
  float metmagcorr = .0;
  double hltWeight = 1.;
  //double e_sf = 0.2719076637;//5.3 fb-1
  
  double e_sf = 1.0;
  if(this->pName == "Z" || this->pName == "DY"){
    e_sf = 1.1973;
  }else if(this->pName == "W"){
    e_sf = 1.2401;
  }else if(this->pName == "TT"){
    e_sf = 1.776;
  }
  
  for(int i = 0; i < T->GetEntries(); i++){
    T->GetEntry(i);
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = (nBtag[0] >= nBtagCut[0]);//Loose >= nBtagCut[0]    
    fBtag[2] = (nBtag[1] >= nBtagCut[1]);//Med >= nBtagCut[1]                               
    fBtag[3] = (nBtag[2] == nBtagCut[2]);//Tight == nBtagCut[2], exclusive
    fBtag[4] = (nBtag[2] >= nBtagCut[2]);//Tight >= nBtagCut[2]
    
    TLorentzVector j1;
    TLorentzVector j2, double_mu;

    j1.SetPtEtaPhiE(pTHem1, etaHem1, phiHem1, pTHem1*cosh(etaHem1));//Hemisphere
    j2.SetPtEtaPhiE(pTHem2, etaHem2, phiHem2, pTHem2*cosh(etaHem2));//Hemisphere
    double Dphi = j1.DeltaPhi(j2);

    metmag = sqrt(metX[metIndex]*metX[metIndex]+metY[metIndex]*metY[metIndex]);
    metmagcorr = sqrt(metcorrX[metIndex]*metcorrX[metIndex]+metcorrY[metIndex]*metcorrY[metIndex]);
    
    if( RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin  && fBtag[btagIndex]  && fabs(Dphi) < 2.5){
      //double met = sqrt(metX[2]*metX[2] + metY[2]*metY[2]);
      //if(!(Jet_PT[0] > 110.0 && fabs(Jet_Eta[0]) < 2.4 && N_Jets == 2 && met > 550.0))continue;
      bool EtaVeto = false;
      for(int k = 0; k < N_Jets; k++){
        if(fabs(Jet_Eta[k]) > 2.4){
          EtaVeto = true;
          break;//Veto events with at least 1 jet |eta|>2.4      
	}
      } 
      if(EtaVeto)continue;

      hltWeight = GetTriggerEff(MR[metIndex], "NoMu");
      if( hltWeight == 0.0 )hltWeight = 1.0;
      if( BOX == 0 ){
	MET[0]->Fill(metmag, sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
        MET[1]->Fill(metmagcorr, sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
        MET[2]->Fill(metmagcorr-metmag, sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
        MET[9]->Fill(N_Jets,sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
      }else if( BOX == 1 ){
	MET[3]->Fill(metmag, sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
        MET[4]->Fill(metmagcorr, sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
        MET[5]->Fill(metmagcorr-metmag, sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
        MET[10]->Fill(N_Jets, sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
      }else if( BOX == 2 ){
	double_mu.SetPxPyPzE(Mu_Px[0]+Mu_Px[1], Mu_Py[0]+Mu_Py[1], Mu_Pz[0]+Mu_Pz[1], Mu_E[0]+Mu_E[1]);
	
	if(double_mu.M() > 80.0 && double_mu.M() < 100.0){
	//if((double_mu.M() > 10.0 && double_mu.M() < 75.0 )||(double_mu.M() >105.0 && double_mu.M() < 200.0)){
	  MET[6]->Fill(metmag, sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
	  MET[7]->Fill(metmagcorr, sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
	  MET[8]->Fill(metmagcorr-metmag, sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
	  MET[11]->Fill(N_Jets, sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
	}
	
      }
    }
  }
  T->SetBranchStatus("*", 0);
  
  for(int j = 0; j < 12; j++){
    metvec.push_back(MET[j]);
  }

  return metvec;
};

std::vector<TH2F*> BaseDM::Plot_2DRazor(){
  double RSQ[4], MR[4], CSV[30];
  double pTHem1, pTHem2, etaHem1, etaHem2, phiHem1, phiHem2, pu_w, mu_w, sf_w, ISR_w;
  double Mu_Px[2], Mu_Py[2], Mu_Pz[2], Mu_E[2];
  int BOX, N_Jets, nBtag[2];
  
  std::vector< TH2F* > Razor2DVec;
  TH2F* Razor2D[3];
  TString name, name1;
  double hltWeight;
  for(int l = 0; l < 3; l++ ){
    name = TString(Form("MR_2D_TT_%dmu_Box",l));
    name1 = TString(Form("R2_2D_TT_%dmu_Box",l));
    Razor2D[l] = new TH2F( name, name, MR_Bins, MR_BinArr,  RSQ_Bins, RSQ_BinArr);
    Razor2D[l]->Sumw2();
  }

  SetStatus();
  T->SetBranchAddress("pu_w", &pu_w);
  T->SetBranchAddress("mu_w", &mu_w);
  T->SetBranchAddress("RSQ", RSQ);
  T->SetBranchAddress("MR", MR);
  T->SetBranchAddress("BOX_NUM", &BOX);
  //T->SetBranchAddress("nBtag", &nBtag[0]);
  T->SetBranchAddress("nBtagLCorr", &nBtag[0]);
  if(this->pName == "TT_L"){
    T->SetBranchAddress("nBtagTCorr", &nBtag[1]);
  }else{
    T->SetBranchAddress("nBtagTight", &nBtag[1]);
  }
  T->SetBranchAddress("N_Jets", &N_Jets);
  T->SetBranchAddress("CSV", CSV);
  T->SetBranchAddress("pTHem1", &pTHem1);
  T->SetBranchAddress("pTHem2", &pTHem2);
  T->SetBranchAddress("etaHem1", &etaHem1);
  T->SetBranchAddress("etaHem2", &etaHem2);
  T->SetBranchAddress("phiHem1", &phiHem1);
  T->SetBranchAddress("phiHem2", &phiHem2);
  T->SetBranchAddress("Mu_Px", Mu_Px);
  T->SetBranchAddress("Mu_Py", Mu_Py);
  T->SetBranchAddress("Mu_Pz", Mu_Pz);
  T->SetBranchAddress("Mu_E", Mu_E);
  T->SetBranchAddress("sf_w", &sf_w);
  T->SetBranchAddress("ISR", &ISR_w);

  //double e_sf = 0.2719076637;
  double e_sf = 1.0;//18.42 fb-1
  for(int i = 0; i < T->GetEntries(); i++){
    T->GetEntry(i);
    TLorentzVector j1;
    TLorentzVector j2, double_mu;
    
    j1.SetPtEtaPhiE(pTHem1, etaHem1, phiHem1, pTHem1*cosh(etaHem1));//Hemisphere
    j2.SetPtEtaPhiE(pTHem2, etaHem2, phiHem2, pTHem2*cosh(etaHem2));//Hemisphere
    double Dphi = j1.DeltaPhi(j2);
    
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= nBtagCut[0]);
    fBtag[3] = (nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    //fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1]);
    fBtag[4] = (nBtag[1] >= nBtagCut[2]);

    if( RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin  && fBtag[btagIndex] && fabs(Dphi) < 2.5){
      if( BOX == 0){
        hltWeight = HLTscaleEle( MR[metIndex], RSQ[metIndex]);
        if( hltWeight == 0.0 )hltWeight = 1.0;
        Razor2D[0]->Fill(MR[metIndex], RSQ[metIndex], sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
      }else if( BOX == 1 ){
        hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
        if( hltWeight == 0.0 )hltWeight = 1.0;
        Razor2D[1]->Fill(MR[metIndex], RSQ[metIndex], sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
      }else if( BOX == 2 ){
	double_mu.SetPxPyPzE(Mu_Px[0]+Mu_Px[1], Mu_Py[0]+Mu_Py[1], Mu_Pz[0]+Mu_Pz[1], Mu_E[0]+Mu_E[1]);
	
	if(double_mu.M() > 75.0 && double_mu.M() < 110.0){
	  hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
	  if( hltWeight == 0.0 )hltWeight = 1.0;
	  Razor2D[2]->Fill(MR[metIndex], RSQ[metIndex], sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
	}
	
      }
    }
  }
  T->SetBranchStatus("*", 0);

  for(int j = 0; j < 3; j++){
    Razor2DVec.push_back(Razor2D[j]);
  }
  
  return Razor2DVec;

};

bool BaseDM::pfJetPassCSVM(double btagOutput){
  if(btagOutput < 0.679)   return false;
  return true;
};

int BaseDM::pfJetPassCSVM(double* CSVM, int N_Jets){
  int nMBtag = 0;
  for(int i = 0; i < N_Jets; i++)if(CSVM[i] >= 0.679)nMBtag++;
  return nMBtag;
};

std::vector<TH1F*> BaseDM::Plot_1DRazor(){
  double RSQ[4], MR[4], CSV[30];
  double pTHem1, pTHem2, etaHem1, etaHem2, phiHem1, phiHem2, pu_w, mu_w, sf_w, ISR_w;
  double Mu_Px[2], Mu_Py[2], Mu_Pz[2], Mu_E[2];

  double metX[4], metY[4];
  double Jet_PT[20], Jet_Eta[20], Jet_Phi[20];

  int BOX, N_Jets, nBtag[3];
  
  std::vector< TH1F* > Razor1DVec;
  TH1F* Razor1D[12];
  TString name, name1;
  double hltWeight;
  for(int l = 0; l < 6; l++){
    name = TString(Form("MR_1D_%dmu_Box",l));
    name1 = TString(Form("R2_1D_%dmu_Box",l));
    //Razor1D[2*l] = new TH1F(name,name, MR_Bins, MR_BinArr);
    Razor1D[2*l] = new TH1F(name,name, 20, 200, 1700.0);
    //Razor1D[2*l+1] = new TH1F(name1,name1, RSQ_Bins, RSQ_BinArr);
    Razor1D[2*l+1] = new TH1F(name1,name1, 20, 0.5, 2.0); 
    if(l < 3){
      Razor1D[2*l]->Sumw2();
      Razor1D[2*l+1]->Sumw2();
    }
  }

  SetStatus();

  T->SetBranchStatus("N_Jets",1);
  T->SetBranchStatus("Jet_PT",1);
  T->SetBranchStatus("Jet_Phi",1);
  T->SetBranchStatus("Jet_Eta",1);
  T->SetBranchStatus("metX",1);
  T->SetBranchStatus("metY",1);

  T->SetBranchAddress("pu_w", &pu_w);
  T->SetBranchAddress("mu_w", &mu_w);
  T->SetBranchAddress("RSQ", RSQ);
  T->SetBranchAddress("MR", MR);
  T->SetBranchAddress("BOX_NUM", &BOX);
  //T->SetBranchAddress("nBtag", &nBtag[0]);
  T->SetBranchAddress("nBtagLCorr", &nBtag[0]);
  if(this->pName == "TT_LL"){
    T->SetBranchAddress("nBtagTCorr", &nBtag[2]);
  }else{
    T->SetBranchAddress("nBtagTight", &nBtag[2]);
  }
  T->SetBranchAddress("nBtagMed", &nBtag[1]);
  T->SetBranchAddress("N_Jets", &N_Jets);
  T->SetBranchAddress("Jet_PT", Jet_PT);
  T->SetBranchAddress("Jet_Phi", Jet_Phi);
  T->SetBranchAddress("Jet_Eta", Jet_Eta);
  T->SetBranchAddress("CSV", CSV);
  T->SetBranchAddress("pTHem1", &pTHem1);
  T->SetBranchAddress("pTHem2", &pTHem2);
  T->SetBranchAddress("etaHem1", &etaHem1);
  T->SetBranchAddress("etaHem2", &etaHem2);
  T->SetBranchAddress("phiHem1", &phiHem1);
  T->SetBranchAddress("phiHem2", &phiHem2);
  T->SetBranchAddress("metX", metX);
  T->SetBranchAddress("metY", metY);
  T->SetBranchAddress("Mu_Px", Mu_Px);
  T->SetBranchAddress("Mu_Py", Mu_Py);
  T->SetBranchAddress("Mu_Pz", Mu_Pz);
  T->SetBranchAddress("Mu_E", Mu_E);
  //T->SetBranchAddress("sf_w", &sf_w);
  T->SetBranchAddress("ISR_w", &sf_w);
  T->SetBranchAddress("ISR", &ISR_w);
  
  double e_sf = 1.0;
  if(this->pName == "Z" || this->pName == "DY"){
    e_sf = 1.1973;
  }else if(this->pName == "W"){
    e_sf = 1.2401;
  }else if(this->pName == "TT"){
    e_sf = 1.776;
  }
  
  for(int i = 0; i < T->GetEntries(); i++){
    T->GetEntry(i);
    TLorentzVector j1, double_mu;
    TLorentzVector j2;
    
    j1.SetPtEtaPhiE(pTHem1, etaHem1, phiHem1, pTHem1*cosh(etaHem1));//Hemisphere   
    j2.SetPtEtaPhiE(pTHem2, etaHem2, phiHem2, pTHem2*cosh(etaHem2));//Hemisphere
    
    double Dphi = j1.DeltaPhi(j2);
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = (nBtag[0] >= nBtagCut[0]);//Loose >= nBtagCut[0]
    fBtag[2] = (nBtag[1] >= nBtagCut[1]);//Med >= nBtagCut[1]
    fBtag[3] = (nBtag[2] == nBtagCut[2]);//Tight == nBtagCut[2], exclusive
    fBtag[4] = (nBtag[2] >= nBtagCut[2]);//Tight >= nBtagCut[2]
    
    if( RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin  && fBtag[btagIndex] && fabs(Dphi) < 2.5){
      //double MET = sqrt(metX[2]*metX[2] + metY[2]*metY[2]);
      //if(!(Jet_PT[0] > 110.0 && fabs(Jet_Eta[0]) < 2.4 && N_Jets == 2 && MET > 550.0))continue;
      bool EtaVeto = false;
      for(int k = 0; k < N_Jets; k++){
        if(fabs(Jet_Eta[k]) > 2.4){
          EtaVeto = true;
          break;//Veto events with at least 1 jet |eta|>2.4                                       
	}
      } 
      if(EtaVeto)continue;
      
      if(BOX == 0){
        //hltWeight = HLTscaleEle( MR[metIndex], RSQ[metIndex]);
	hltWeight = GetTriggerEff(MR[metIndex], "NoMu");
	if(MR[metIndex]> 1700.0)MR[metIndex] = 1699.0;
	if(RSQ[metIndex]> 2.0)RSQ[metIndex] = 1.99;
        if( hltWeight == 0.0 )hltWeight = 1.0;
        Razor1D[0]->Fill(MR[metIndex], sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
        Razor1D[1]->Fill(RSQ[metIndex], sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
        Razor1D[6]->Fill(MR[metIndex], sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
        Razor1D[7]->Fill(RSQ[metIndex], sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
      }else if(BOX == 1){
        //hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
	hltWeight = GetTriggerEff(MR[metIndex], "NoMu");
        if( hltWeight == 0.0 )hltWeight = 1.0;
	if(MR[metIndex]> 1700.0)MR[metIndex] = 1699.0;
	if(RSQ[metIndex]> 2.0)RSQ[metIndex] = 1.99;
        Razor1D[2]->Fill(MR[metIndex], sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
        Razor1D[3]->Fill(RSQ[metIndex], sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
        Razor1D[8]->Fill(MR[metIndex], sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
        Razor1D[9]->Fill(RSQ[metIndex], sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
      }else if(BOX == 2){
	double_mu.SetPxPyPzE(Mu_Px[0]+Mu_Px[1], Mu_Py[0]+Mu_Py[1], Mu_Pz[0]+Mu_Pz[1], Mu_E[0]+Mu_E[1]);
	
	//if(double_mu.M() > 75.0 && double_mu.M() < 110.0){	
	if(double_mu.M() > 80.0 && double_mu.M() < 100.0){
	//if((double_mu.M() > 10.0 && double_mu.M() < 75.0 )||(double_mu.M() >105.0 && double_mu.M() < 200.0)){
	  //hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
	  hltWeight = GetTriggerEff(MR[metIndex], "NoMu");
	  if(MR[metIndex]> 1700.0)MR[metIndex] = 1699.0;
	  if(RSQ[metIndex]> 2.0)RSQ[metIndex] = 1.99;
	  if( hltWeight == 0.0 )hltWeight = 1.0;
	  Razor1D[4]->Fill(MR[metIndex], sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
	  Razor1D[5]->Fill(RSQ[metIndex], sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
	  Razor1D[10]->Fill(MR[metIndex], sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
	  Razor1D[11]->Fill(RSQ[metIndex], sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
	}
	
      }
    }

  }

  for(int j = 0; j < 12; j++){
    Razor1DVec.push_back(Razor1D[j]);
  }

  return Razor1DVec;
  
};


std::vector<TH1F*> BaseDM::Plot_MRCat(){
 double RSQ[4], MR[4], CSV[30];
 double pTHem1, pTHem2, etaHem1, etaHem2, phiHem1, phiHem2, pu_w, mu_w, sf_w, ISR_w;
 double Mu_Px[2], Mu_Py[2], Mu_Pz[2], Mu_E[2];

 double metX[4], metY[4];
 double Jet_PT[20], Jet_Eta[20], Jet_Phi[20];

  int BOX, N_Jets, nBtag[3];
  int c1_bins_1 = 11;
  float c1B_1[] = {0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.9, 0.95, 1.0, 1.2};
  int c2_bins_1 = 6;
  float c2B_1[] = {0.50, 0.575, 0.65, 0.75, 0.85, .950, 1.2};
  int c3_bins_1 = 6;
  float c3B_1[] = {0.50, 0.575, 0.65, 0.75, 0.85, .950, 1.2};
  int c4_bins_1 = 4;
  float c4B_1[] = {0.50, 0.60, 0.70, .950, 1.20};
  //v2
  /*
  int c1_bins_2 = 7;
  float c1B_2[] = {0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 1.2};
  int c2_bins_2 = 4;
  float c2B_2[] = {0.50, 0.575, 0.65, 0.75, 1.2};
  int c3_bins_2 = 4;
  float c3B_2[] = {0.50, 0.575, 0.65, 0.75, 1.2};
  int c4_bins_2 = 4;
  float c4B_2[] = {0.50, 0.575, 0.65, 0.75, 1.2};
  */
  //v3
  
  int c1_bins_2 = 4;
  //float c1B_2[] = {0.50, 0.57, 0.65, 0.73, 0.81, 0.91, 1.2};
  //float c1B_2[] = {0.50, 0.6, 0.70, 0.85, 1.2};    
  //float c1B_2[] = {0.50, 0.6, 0.75, 0.95, 1.2};
  //float c1B_2[] = {0.50, 0.625, 0.80, 1.05, 1.2};//(last)
  float c1B_2[] = {0.50, 0.6, 0.75, .9, 1.2}; 
  //float c1B_2[] = {0.50, 0.57, 0.65, 0.73, 0.81, 0.91, 1.01, 1.2};
  int c2_bins_2 = 5;
  float c2B_2[] = {0.50, 0.57, 0.65, 0.73, 0.81, 1.2};
  int c3_bins_2 = 5;
  float c3B_2[] = {0.50, 0.57, 0.65, 0.73, 0.81, 1.2};
  int c4_bins_2 = 5;
  float c4B_2[] = {0.50, 0.57, 0.65, 0.73, 0.81, 1.2};
  
  
  int c1_bins, c2_bins, c3_bins, c4_bins;
  float* c1B;
  float* c2B;
  float* c3B;
  float* c4B;
  std::cout << "bTag Index: " << btagIndex << std::endl;
  if(btagIndex == 0){
    c1_bins = 11;
    c1B = c1B_1;
    c2_bins = 6;
    c2B = c2B_1;
    c3_bins = 6;
    c3B = c3B_1;
    c4_bins = 4;
    c4B = c4B_1;
  }else{
    c1_bins = 4;
    c1B = c1B_2;
    c2_bins = 4;
    c2B = c2B_2;
    c3_bins = 4;
    c3B = c3B_2;
    c4_bins = 4;
    c4B = c4B_2;
  }
  
  /*
    c1_bins = 11;
    c1B = c1B_3;
    c2_bins = 6;
    c2B = c2B_3;
    c3_bins = 6;
    c3B = c3B_3;
    c4_bins = 4;
    c4B = c4B_3;
  */
  
  std::vector< TH1F* > Razor1DVec;
  TH1F* Razor1D[12];
  TH1F* Razor1DInc[3];
  TString name1, name2, name3, name4;
  double hltWeight;
  for(int l = 0; l < 3; l++){
    //name1 = TString(Form(this->pName+"_cat1_1D_%dmu_Box",l));
    name2 = TString(Form(this->pName+"_cat2_1D_%dmu_Box",l));
    name3 = TString(Form(this->pName+"_cat3_1D_%dmu_Box",l));
    name4 = TString(Form(this->pName+"_cat4_1D_%dmu_Box",l));
    
    name1 = TString(Form(this->pName+"_INC_1D_%dmu_Box",l));
    Razor1DInc[l] = new TH1F(name1, name1, c1_bins, c1B);
    Razor1DInc[l]->Sumw2();
    
    name1 = TString(Form(this->pName+"_cat1_1D_%dmu_Box",l));
    Razor1D[4*l] = new TH1F(name1,name1, c1_bins, c1B);
    Razor1D[4*l+1] = new TH1F(name2,name2, c2_bins, c2B);
    Razor1D[4*l+2] = new TH1F(name3,name3, c3_bins, c3B);
    Razor1D[4*l+3] = new TH1F(name4,name4, c4_bins, c4B);
    
    if(l < 3){
      Razor1D[4*l]->Sumw2();
      Razor1D[4*l+1]->Sumw2();
      Razor1D[4*l+2]->Sumw2();
      Razor1D[4*l+3]->Sumw2();
    }
  }

  SetStatus();
  T->SetBranchStatus("N_Jets",1);
  T->SetBranchStatus("Jet_PT",1);
  T->SetBranchStatus("Jet_Phi",1);
  T->SetBranchStatus("Jet_Eta",1);
  T->SetBranchStatus("metX",1);
  T->SetBranchStatus("metY",1);
  
  T->SetBranchAddress("pu_w", &pu_w);
  T->SetBranchAddress("mu_w", &mu_w);
  T->SetBranchAddress("RSQ", RSQ);
  T->SetBranchAddress("MR", MR);
  T->SetBranchAddress("BOX_NUM", &BOX);
  T->SetBranchAddress("nBtagLCorr", &nBtag[0]);
  //T->SetBranchAddress("nBtag", &nBtag[0]);
  if(this->pName == "TT_LL"){
    T->SetBranchAddress("nBtagTCorr", &nBtag[2]);
  }else{
    T->SetBranchAddress("nBtagTight", &nBtag[2]);
  }
  T->SetBranchAddress("nBtagMed", &nBtag[1]);
  T->SetBranchAddress("N_Jets", &N_Jets);
  T->SetBranchAddress("Jet_PT", Jet_PT);
  T->SetBranchAddress("Jet_Phi", Jet_Phi);
  T->SetBranchAddress("Jet_Eta", Jet_Eta);
  T->SetBranchAddress("CSV", CSV);
  T->SetBranchAddress("pTHem1", &pTHem1);
  T->SetBranchAddress("pTHem2", &pTHem2);
  T->SetBranchAddress("etaHem1", &etaHem1);
  T->SetBranchAddress("etaHem2", &etaHem2);
  T->SetBranchAddress("phiHem1", &phiHem1);
  T->SetBranchAddress("phiHem2", &phiHem2);
  T->SetBranchAddress("metX", metX);
  T->SetBranchAddress("metY", metY);
  T->SetBranchAddress("Mu_Px", Mu_Px);
  T->SetBranchAddress("Mu_Py", Mu_Py);
  T->SetBranchAddress("Mu_Pz", Mu_Pz);
  T->SetBranchAddress("Mu_E", Mu_E);
  T->SetBranchAddress("ISR_w", &sf_w);//Corrected number of Initial events
  T->SetBranchAddress("ISR", &ISR_w);//individual weight for ISR
  //T->SetBranchAddress("sf_w", &sf_w);
  //ISR_w = 1.0;
  
  double e_sf = 1.00;//already 18.5 fb^-1
  //std::cout << "Loose: " << nBtagCut[0] << " Med: " << nBtagCut[1] << " Tight: " << nBtagCut[2] << std::endl;
  for(int i = 0; i < T->GetEntries(); i++){
    T->GetEntry(i);
    TLorentzVector j1, double_mu;
    TLorentzVector j2;

    j1.SetPtEtaPhiE(pTHem1, etaHem1, phiHem1, pTHem1*cosh(etaHem1));//Hemisphere   
    j2.SetPtEtaPhiE(pTHem2, etaHem2, phiHem2, pTHem2*cosh(etaHem2));//Hemisphere
    
    double Dphi = j1.DeltaPhi(j2);
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = (nBtag[0] >= nBtagCut[0]);//Loose >= nBtagCut[0]
    fBtag[2] = (nBtag[1] >= nBtagCut[1]);//Med >= nBtagCut[1]
    fBtag[3] = (nBtag[2] == nBtagCut[2]);//Tight == nBtagCut[2], exclusive
    fBtag[4] = (nBtag[2] >= nBtagCut[2]);//Tight >= nBtagCut[2]
    
    if( RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin  && fBtag[btagIndex] && fabs(Dphi) < 2.5){
      //double MET = sqrt(metX[2]*metX[2] + metY[2]*metY[2]);
      //if(!(Jet_PT[0] > 110.0 && fabs(Jet_Eta[0]) < 2.4 && N_Jets == 2 && MET > 550.0))continue;
      bool EtaVeto = false;
      for(int k = 0; k < N_Jets; k++){
        if(fabs(Jet_Eta[k]) > 2.4){
          EtaVeto = true;
          break;//Veto events with at least 1 jet |eta|>2.4                      
	}
      } 
      if(EtaVeto)continue;
      
      if(RSQ[metIndex] > 1.20)RSQ[metIndex] = 1.1;
      if(BOX == 0){
        //hltWeight = HLTscaleEle( MR[metIndex], RSQ[metIndex]);
	hltWeight = GetTriggerEff(MR[metIndex], "NoMu");
	if(i%100000 == 0)std::cout << "HLT No Mu: " << hltWeight << std::endl; 
        if( hltWeight == 0.0 )hltWeight = 1.0;
	Razor1DInc[0]->Fill(RSQ[metIndex], sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
	if(MR[metIndex] > 200.0 && MR[metIndex] <= 300.0 ){
	  Razor1D[0]->Fill(RSQ[metIndex], sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
	}else if(MR[metIndex] > 300.0 && MR[metIndex] <= 400.0){
	  Razor1D[1]->Fill(RSQ[metIndex], sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
	}else if(MR[metIndex] > 400.0 && MR[metIndex] <= 600.0){
	  Razor1D[2]->Fill(RSQ[metIndex], sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
	}else if(MR[metIndex] > 600.0 && MR[metIndex] <= 3500.0){
	  Razor1D[3]->Fill(RSQ[metIndex], sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
	}
      }else if(BOX == 1){
        //hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
	//hltWeight = HLTscaleEle( MR[metIndex], RSQ[metIndex]);
	//hltWeight = GetTriggerEff(MR[metIndex], "Mu");
	hltWeight = GetTriggerEff(MR[metIndex], "NoMu");
	if(i%100000 == 0)std::cout << "HLT Mu: " << hltWeight << std::endl;
        if( hltWeight == 0.0 )hltWeight = 1.0;
	Razor1DInc[1]->Fill(RSQ[metIndex], sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
	if(MR[metIndex] > 200.0 && MR[metIndex] <= 300.0 ){
	  Razor1D[4]->Fill(RSQ[metIndex], sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
	}else if(MR[metIndex] > 300.0 && MR[metIndex] <= 400.0){
	  Razor1D[5]->Fill(RSQ[metIndex], sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
	}else if(MR[metIndex] > 400.0 && MR[metIndex] <= 600.0){
	  Razor1D[6]->Fill(RSQ[metIndex], sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
	}else if(MR[metIndex] > 600.0 && MR[metIndex] <= 3500.0){
	  Razor1D[7]->Fill(RSQ[metIndex], sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
	}
      }else if(BOX == 2){
	double_mu.SetPxPyPzE(Mu_Px[0]+Mu_Px[1], Mu_Py[0]+Mu_Py[1], Mu_Pz[0]+Mu_Pz[1], Mu_E[0]+Mu_E[1]);
	if(double_mu.M() > 80.0 && double_mu.M() < 100.0){  
	  //if((double_mu.M() > 10.0 && double_mu.M() < 75.0 )||(double_mu.M() > 105.0 && double_mu.M() < 200.0)){
	  //hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
	  //hltWeight = HLTscaleEle( MR[metIndex], RSQ[metIndex]);
	  //hltWeight = GetTriggerEff(MR[metIndex], "Mu");
	  hltWeight = GetTriggerEff(MR[metIndex], "NoMu"); 
	  if( hltWeight == 0.0 )hltWeight = 1.0;
	  Razor1DInc[2]->Fill(RSQ[metIndex], sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
	  if(MR[metIndex] > 200.0 && MR[metIndex] <= 300.0 ){
	    Razor1D[8]->Fill(RSQ[metIndex], sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
	  }else if(MR[metIndex] > 300.0 && MR[metIndex] <= 400.0){
	    Razor1D[9]->Fill(RSQ[metIndex], sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
	  }else if(MR[metIndex] > 400.0 && MR[metIndex] <= 600.0){
	    Razor1D[10]->Fill(RSQ[metIndex], sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
	  }else if(MR[metIndex] > 600.0 && MR[metIndex] <= 3500.0){
	    Razor1D[11]->Fill(RSQ[metIndex], sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
	  }
	}
	
      }
    }

  }

  for(int j = 0; j < 12; j++){
    Razor1DVec.push_back(Razor1D[j]);
  }
  for(int j = 0; j < 3; j++){
    Razor1DVec.push_back(Razor1DInc[j]);
  }

  return Razor1DVec;
  
}

std::vector<TH1F*> BaseDM::PlotKine(){
  double RSQ[4], MR[4], CSV[30];
  double Jet_PT[30], Jet_Eta[30], Jet_Phi[30], Mu_Px[2], Mu_Py[2], Mu_Pz[2], Mu_E[2];
  double pTHem1, pTHem2, etaHem1, etaHem2, phiHem1, phiHem2, pu_w, mu_w, sf_w, ISR_w;
  int BOX, N_Jets, nBtag[3];

  double metX[4], metY[4];
  
  std::vector< TH1F* > Razor1DKVec;
  TH1F* Razor1DK[68];
  TString name[6];
  double hltWeight;
  
  for(int l = 0; l < 6; l++){
    name[0] = TString(Form("PT_1D_J1_%dmu_BoxMC",l));
    name[1] = TString(Form("ETA_1D_J1_%dmu_BoxMC",l));
    name[2] = TString(Form("PHI_1D_J1_%dmu_BoxMC",l));
    name[3] = TString(Form("PT_1D_J2_%dmu_BoxMC",l));
    name[4] = TString(Form("ETA_1D_J2_%dmu_BoxMC",l));
    name[5] = TString(Form("PHI_1D_J2_%dmu_BoxMC",l));
    
    Razor1DK[6*l] = new TH1F(name[0],name[0], 20, 80, 1500);
    Razor1DK[6*l+1] = new TH1F(name[1],name[1], 20, -3.0, 3.0);
    Razor1DK[6*l+2] = new TH1F(name[2],name[2], 20, -TMath::Pi(), TMath::Pi());
    Razor1DK[6*l+3] = new TH1F(name[3],name[3], 20, 80, 1500);
    Razor1DK[6*l+4] = new TH1F(name[4],name[4], 20, -3.0, 3.0);
    Razor1DK[6*l+5] = new TH1F(name[5],name[5], 20, -TMath::Pi(), TMath::Pi());
    
    if(l < 3){
      for(int k = 0; k < 6; k++)Razor1DK[6*l+k]->Sumw2();
    }
  }
  
  for(int l = 0; l < 2; l++){
    int box_mu = 2;
    name[0] = TString(Form("PT_1D_m1_%dmu_BoxMC_v%d",box_mu, l));
    name[1] = TString(Form("ETA_1D_m1_%dmu_BoxMC_v%d",box_mu, l));
    name[2] = TString(Form("PHI_1D_m1_%dmu_BoxMC_v%d",box_mu, l));
    name[3] = TString(Form("PT_1D_m2_%dmu_BoxMC_v%d",box_mu, l));
    name[4] = TString(Form("ETA_1D_m2_%dmu_BoxMC_v%d",box_mu, l));
    name[5] = TString(Form("PHI_1D_m2_%dmu_BoxMC_v%d",box_mu, l));
    
    Razor1DK[36+9*l] = new TH1F(name[0],name[0], 20, 0.0, 1500);
    Razor1DK[36+9*l+1] = new TH1F(name[1],name[1], 20, -3.0, 3.0);
    Razor1DK[36+9*l+2] = new TH1F(name[2],name[2], 20, -TMath::Pi(), TMath::Pi());
    Razor1DK[36+9*l+3] = new TH1F(name[3],name[3], 20, 0.0, 1500);
    Razor1DK[36+9*l+4] = new TH1F(name[4],name[4], 20, -3.0, 3.0);
    Razor1DK[36+9*l+5] = new TH1F(name[5],name[5], 20, -TMath::Pi(), TMath::Pi());
    
    box_mu = 1;
    name[0] = TString(Form("PT_1D_m1_%dmu_BoxMC_v%d",box_mu, l));
    name[1] = TString(Form("ETA_1D_m1_%dmu_BoxMC_v%d",box_mu, l));
    name[2] = TString(Form("PHI_1D_m1_%dmu_BoxMC_v%d",box_mu, l));
    
    Razor1DK[36+9*l+6] = new TH1F(name[0],name[0], 20, 0.0, 1500);
    Razor1DK[36+9*l+7] = new TH1F(name[1],name[1], 20, -3.0, 3.0);
    Razor1DK[36+9*l+8] = new TH1F(name[2],name[2], 20, -TMath::Pi(), TMath::Pi());
    
    if(l < 1){
      for(int k = 0; k < 9; k++)Razor1DK[36+9*l+k]->Sumw2();
    }
  }
  
  //Razor1DK[54] = new TH1F("diMuMass", "diMuMass", 10, 75.0, 110.0);
  //Razor1DK[55] = new TH1F("diMuMass_err", "diMuMass_err", 10, 75.0, 110.0);
  
  Razor1DK[54] = new TH1F("diMuMass", "diMuMass", 20, 80.0, 100.0);
  Razor1DK[55] = new TH1F("diMuMass_err", "diMuMass_err", 20, 80.0, 100.0);
  
  //Razor1DK[54] = new TH1F("diMuMass", "diMuMass", 50, 0.0, 200.0);
  //Razor1DK[55] = new TH1F("diMuMass_err", "diMuMass_err", 50, 0.0, 200.0);
  Razor1DK[55]->Sumw2();
  
  //HT and DeltaPhi
  for(int k = 0; k < 6; k++){
    TString n(Form("HT_%dmu_%d",k%3,k));
    Razor1DK[56+k] = new TH1F(n, n, 20, 0.0, 1400.0);
    n = Form("Dphi_%dmu_%d",k%3,k);
    Razor1DK[62+k] = new TH1F(n, n, 20,  -TMath::Pi(), TMath::Pi());
  }
  
  SetStatusKine();
  T->SetBranchStatus("metX",1);
  T->SetBranchStatus("metY",1);
  
  T->SetBranchAddress("pu_w", &pu_w);
  T->SetBranchAddress("mu_w", &mu_w);
  T->SetBranchAddress("RSQ", RSQ);
  T->SetBranchAddress("MR", MR);
  T->SetBranchAddress("BOX_NUM", &BOX);
  //T->SetBranchAddress("nBtag", &nBtag[0]);
  T->SetBranchAddress("nBtagLCorr", &nBtag[0]);
  if(this->pName == "TT_LL"){
    T->SetBranchAddress("nBtagTCorr", &nBtag[2]);
  }else{
    T->SetBranchAddress("nBtagTight", &nBtag[2]);
  }
  T->SetBranchAddress("nBtagMed", &nBtag[1]);
  T->SetBranchAddress("N_Jets", &N_Jets);
  T->SetBranchAddress("Jet_PT", Jet_PT);
  T->SetBranchAddress("Jet_Eta", Jet_Eta);
  T->SetBranchAddress("Jet_Phi", Jet_Phi);
  T->SetBranchAddress("CSV", CSV);
  T->SetBranchAddress("pTHem1", &pTHem1);
  T->SetBranchAddress("pTHem2", &pTHem2);
  T->SetBranchAddress("etaHem1", &etaHem1);
  T->SetBranchAddress("etaHem2", &etaHem2);
  T->SetBranchAddress("phiHem1", &phiHem1);
  T->SetBranchAddress("phiHem2", &phiHem2);
  T->SetBranchAddress("metX", metX);
  T->SetBranchAddress("metY", metY);
  T->SetBranchAddress("Mu_Px", Mu_Px);
  T->SetBranchAddress("Mu_Py", Mu_Py);
  T->SetBranchAddress("Mu_Pz", Mu_Pz);
  T->SetBranchAddress("Mu_E", Mu_E);
  //T->SetBranchAddress("sf_w", &sf_w);
  T->SetBranchAddress("ISR_w", &sf_w);
  T->SetBranchAddress("ISR", &ISR_w);
  
  double e_sf = 1.0;//18.42 fb-1
  if(this->pName == "Z" || this->pName == "DY"){
    e_sf = 1.1973;
  }else if(this->pName == "W"){
    e_sf = 1.2401;
  }else if(this->pName == "TT"){
    e_sf = 1.776;
  }
  
  for(int i = 0; i < T->GetEntries(); i++){
    T->GetEntry(i);
    TLorentzVector j1, j2;//Megajet
    TLorentzVector mu1, mu2, double_mu;
    
    j1.SetPtEtaPhiE(pTHem1, etaHem1, phiHem1, pTHem1*cosh(etaHem1));//Hemisphere1
    j2.SetPtEtaPhiE(pTHem2, etaHem2, phiHem2, pTHem2*cosh(etaHem2));//Hemisphere2
    
    double Dphi = j1.DeltaPhi(j2);
    float HT = 0.0;
    for(int iJ = 0; iJ < N_Jets; iJ++){
      HT += Jet_PT[iJ];
    }
    
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = (nBtag[0] >= nBtagCut[0]);//Loose >= nBtagCut[0]                           
    fBtag[2] = (nBtag[1] >= nBtagCut[1]);//Med >= nBtagCut[1]
    fBtag[3] = (nBtag[2] == nBtagCut[2]);//Tight == nBtagCut[2], exclusive
    fBtag[4] = (nBtag[2] >= nBtagCut[2]);//Tight >= nBtagCut[2]
    
    if( RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin  && fBtag[btagIndex] && fabs(Dphi) < 2.5){
      //double MET = sqrt(metX[2]*metX[2] + metY[2]*metY[2]);
      //if(!(Jet_PT[0] > 110.0 && fabs(Jet_Eta[0]) < 2.4 && N_Jets == 2 && MET > 550.0))continue;
      bool EtaVeto = false;
      for(int k = 0; k < N_Jets; k++){
	if(fabs(Jet_Eta[k]) > 2.4){
	  EtaVeto = true;
	  break;//Veto events with at least 1 jet |eta|>2.4           
	}
      }
      if(EtaVeto)continue;
      
      if(BOX == 0){
        hltWeight = GetTriggerEff(MR[metIndex], "NoMu");
        if( hltWeight == 0.0 )hltWeight = 1.0;
        Razor1DK[0]->Fill(Jet_PT[0], sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
        Razor1DK[1]->Fill(Jet_Eta[0], sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
        Razor1DK[2]->Fill(Jet_Phi[0], sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
        Razor1DK[3]->Fill(Jet_PT[1], sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
	Razor1DK[4]->Fill(Jet_Eta[1], sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
	Razor1DK[5]->Fill(Jet_Phi[1], sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
	
	Razor1DK[18]->Fill(Jet_PT[0], sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
        Razor1DK[18+1]->Fill(Jet_Eta[0], sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
        Razor1DK[18+2]->Fill(Jet_Phi[0], sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
        Razor1DK[18+3]->Fill(Jet_PT[1], sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
        Razor1DK[18+4]->Fill(Jet_Eta[1], sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
        Razor1DK[18+5]->Fill(Jet_Phi[1], sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
	
	//HT and Dphi
	Razor1DK[56]->Fill(HT, sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
	Razor1DK[56+3]->Fill(HT, sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);//No sumw2
	Razor1DK[62]->Fill(Dphi, sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
	Razor1DK[62+3]->Fill(Dphi, sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);//No sumw2
	
      }else if(BOX == 1){
        hltWeight = GetTriggerEff(MR[metIndex], "NoMu");
        if( hltWeight == 0.0 )hltWeight = 1.0;
        //Jet Plots
	Razor1DK[6]->Fill(Jet_PT[0], sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
        Razor1DK[7]->Fill(Jet_Eta[0], sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
        Razor1DK[8]->Fill(Jet_Phi[0], sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
        Razor1DK[9]->Fill(Jet_PT[1], sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
        Razor1DK[10]->Fill(Jet_Eta[1], sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
        Razor1DK[11]->Fill(Jet_Phi[1], sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
	
	Razor1DK[18+6]->Fill(Jet_PT[0], sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
        Razor1DK[18+7]->Fill(Jet_Eta[0], sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
        Razor1DK[18+8]->Fill(Jet_Phi[0], sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
        Razor1DK[18+9]->Fill(Jet_PT[1], sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
        Razor1DK[18+10]->Fill(Jet_Eta[1], sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
        Razor1DK[18+11]->Fill(Jet_Phi[1], sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
	
	//Mu Plots
	mu1.SetPxPyPzE(Mu_Px[0], Mu_Py[0], Mu_Pz[0], Mu_E[0]);
	
	Razor1DK[36+6]->Fill(mu1.Pt(), sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
	Razor1DK[36+7]->Fill(mu1.Eta(), sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
	Razor1DK[36+8]->Fill(mu1.Phi(), sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
	
	Razor1DK[36+15]->Fill(mu1.Pt(), sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
	Razor1DK[36+16]->Fill(mu1.Eta(), sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
	Razor1DK[36+17]->Fill(mu1.Phi(), sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
	
	//HT and Dphi
	Razor1DK[57]->Fill(HT, sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
	Razor1DK[57+3]->Fill(HT, sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);//No sumw2
	Razor1DK[63]->Fill(Dphi, sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
	Razor1DK[63+3]->Fill(Dphi, sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);//No sumw2
	
      }else if(BOX == 2){
	double_mu.SetPxPyPzE(Mu_Px[0]+Mu_Px[1], Mu_Py[0]+Mu_Py[1], Mu_Pz[0]+Mu_Pz[1], Mu_E[0]+Mu_E[1]);
	
	if(double_mu.M() > 80.0 && double_mu.M() < 100.0){
	//if((double_mu.M() > 10.0 && double_mu.M() < 75.0 )||(double_mu.M() >105.0 && double_mu.M() < 200.0)){
	  hltWeight = GetTriggerEff(MR[metIndex], "NoMu");
	  if( hltWeight == 0.0 )hltWeight = 1.0;
	  Razor1DK[54]->Fill(double_mu.M(), sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
	  Razor1DK[55]->Fill(double_mu.M(), sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
	  
	  Razor1DK[12]->Fill(Jet_PT[0], sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
	  Razor1DK[13]->Fill(Jet_Eta[0], sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
	  Razor1DK[14]->Fill(Jet_Phi[0], sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
	  Razor1DK[15]->Fill(Jet_PT[1], sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
	  Razor1DK[16]->Fill(Jet_Eta[1], sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
	  Razor1DK[17]->Fill(Jet_Phi[1], sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
	  
	  Razor1DK[18+12]->Fill(Jet_PT[0], sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
	  Razor1DK[18+13]->Fill(Jet_Eta[0], sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
	  Razor1DK[18+14]->Fill(Jet_Phi[0], sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
	  Razor1DK[18+15]->Fill(Jet_PT[1], sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
	  Razor1DK[18+16]->Fill(Jet_Eta[1], sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
	  Razor1DK[18+17]->Fill(Jet_Phi[1], sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
	  
	  //Mu Plots
	  mu1.SetPxPyPzE(Mu_Px[0], Mu_Py[0], Mu_Pz[0], Mu_E[0]);
	  mu2.SetPxPyPzE(Mu_Px[1], Mu_Py[1], Mu_Pz[1], Mu_E[1]);
	  
	  Razor1DK[36]->Fill(mu1.Pt(), sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
	  Razor1DK[36+1]->Fill(mu1.Eta(), sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
	  Razor1DK[36+2]->Fill(mu1.Phi(), sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
	  Razor1DK[36+3]->Fill(mu2.Pt(), sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
	  Razor1DK[36+4]->Fill(mu2.Eta(), sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
	  Razor1DK[36+5]->Fill(mu2.Phi(), sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
	  
	  Razor1DK[36+9]->Fill(mu1.Pt(), sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
	  Razor1DK[36+10]->Fill(mu1.Eta(), sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
	  Razor1DK[36+11]->Fill(mu1.Phi(), sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
	  Razor1DK[36+12]->Fill(mu2.Pt(), sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
	  Razor1DK[36+13]->Fill(mu2.Eta(), sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
	  Razor1DK[36+14]->Fill(mu2.Phi(), sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);

	  //HT and Dphi
	  Razor1DK[58]->Fill(HT, sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
	  Razor1DK[58+3]->Fill(HT, sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);//No sumw2
	  Razor1DK[64]->Fill(Dphi, sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);
	  Razor1DK[64+3]->Fill(Dphi, sf_w*hltWeight*pu_w*mu_w*e_sf*ISR_w);//No sumw2
	}
	
      }
    }

  }
  
  for(int j = 0; j < 68; j++){
    Razor1DKVec.push_back(Razor1DK[j]);
  }
  
  return Razor1DKVec;

};

bool BaseDM::SetStatus(){
  T->SetBranchStatus("*",0); //disable all branches
  T->SetBranchStatus("pu_w",1);
  T->SetBranchStatus("mu_w",1);
  T->SetBranchStatus("mu_w_up", 1);
  T->SetBranchStatus("mu_w_down",1);
  T->SetBranchStatus("RSQ",1);
  T->SetBranchStatus("RSQ_up",1);
  T->SetBranchStatus("RSQ_down",1);
  T->SetBranchStatus("MR",1);
  T->SetBranchStatus("MR_up",1);
  T->SetBranchStatus("MR_down",1);
  T->SetBranchStatus("BOX_NUM",1);
  T->SetBranchStatus("nBtag",1);
  T->SetBranchStatus("nBtagLCorr", 1);
  T->SetBranchStatus("nBtagLCorrUp", 1);
  T->SetBranchStatus("nBtagLCorrDown", 1);
  if(this->pName == "TT_L"){
    T->SetBranchStatus("nBtagTCorr", 1);
    T->SetBranchStatus("nBtagTCorrUp", 1);
    T->SetBranchStatus("nBtagTCorrDown", 1);
  }else{
    T->SetBranchStatus("nBtagTight", 1);
  }
  T->SetBranchStatus("N_Jets",1);
  T->SetBranchStatus("CSV",1);
  
  T->SetBranchStatus("pTHem1", 1);
  T->SetBranchStatus("pTHem1_up", 1);
  T->SetBranchStatus("pTHem1_down", 1);
  T->SetBranchStatus("pTHem2", 1);
  T->SetBranchStatus("pTHem2_up", 1);
  T->SetBranchStatus("pTHem2_down", 1);
  
  T->SetBranchStatus("etaHem1", 1);
  T->SetBranchStatus("etaHem1_up", 1);
  T->SetBranchStatus("etaHem1_down", 1);
  T->SetBranchStatus("etaHem2", 1);
  T->SetBranchStatus("etaHem2_up", 1);
  T->SetBranchStatus("etaHem2_down", 1);
  
  T->SetBranchStatus("phiHem1", 1);
  T->SetBranchStatus("phiHem1_up", 1);
  T->SetBranchStatus("phiHem1_down", 1);
  T->SetBranchStatus("phiHem2", 1);
  T->SetBranchStatus("phiHem2_up", 1);
  T->SetBranchStatus("phiHem2_down", 1);
  
  T->SetBranchStatus("Mu_Px",1);
  T->SetBranchStatus("Mu_Py",1);
  T->SetBranchStatus("Mu_Pz",1);
  T->SetBranchStatus("Mu_E",1);
  T->SetBranchStatus("sf_w", 1);
  T->SetBranchStatus("ISR", 1);
  T->SetBranchStatus("ISR_up", 1);
  T->SetBranchStatus("ISR_down", 1);
  T->SetBranchStatus("ISR_w", 1);
  T->SetBranchStatus("ISR_w_up", 1);
  T->SetBranchStatus("ISR_w_down", 1);
};

bool BaseDM::SetStatusKine(){
  T->SetBranchStatus("*",0); //disable all branches
  T->SetBranchStatus("pu_w",1);
  T->SetBranchStatus("mu_w",1);
  T->SetBranchStatus("mu_w_up", 1);
  T->SetBranchStatus("mu_w_down",1);
  T->SetBranchStatus("RSQ",1);
  T->SetBranchStatus("MR",1);
  T->SetBranchStatus("BOX_NUM",1);
  T->SetBranchStatus("nBtag",1);
  T->SetBranchStatus("nBtagLCorr", 1);
  T->SetBranchStatus("nBtagLCorrUp", 1);
  T->SetBranchStatus("nBtagLCorrDown", 1);
  if(this->pName == "TT_L"){
    T->SetBranchStatus("nBtagTCorr", 1);
    T->SetBranchStatus("nBtagTCorrUp", 1);
    T->SetBranchStatus("nBtagTCorrDown", 1);
  }else{
    T->SetBranchStatus("nBtagTight", 1);
  }
  T->SetBranchStatus("N_Jets",1);
  T->SetBranchStatus("Jet_PT",1);
  T->SetBranchStatus("Jet_Eta",1);
  T->SetBranchStatus("Jet_Phi",1);
  T->SetBranchStatus("CSV",1);
  T->SetBranchStatus("pTHem1", 1);
  T->SetBranchStatus("pTHem2", 1);
  T->SetBranchStatus("etaHem1", 1);
  T->SetBranchStatus("etaHem2", 1);
  T->SetBranchStatus("phiHem1", 1);
  T->SetBranchStatus("phiHem2", 1);
  T->SetBranchStatus("Mu_Px",1);
  T->SetBranchStatus("Mu_Py",1);
  T->SetBranchStatus("Mu_Pz",1);
  T->SetBranchStatus("Mu_E",1);
  T->SetBranchStatus("sf_w", 1);
  T->SetBranchStatus("ISR", 1);
  T->SetBranchStatus("ISR_up", 1);
  T->SetBranchStatus("ISR_down", 1);
};

bool BaseDM::SetMetStatus(){
  T->SetBranchStatus("*",0); //disable all branches                               
  T->SetBranchStatus("pu_w",1);
  T->SetBranchStatus("mu_w",1);
  T->SetBranchStatus("mu_w_up", 1);
  T->SetBranchStatus("mu_w_down",1);
  T->SetBranchStatus("RSQ",1);
  T->SetBranchStatus("MR",1);
  T->SetBranchStatus("BOX_NUM",1);
  T->SetBranchStatus("nBtag",1);
  T->SetBranchStatus("nBtagLCorr", 1);
  T->SetBranchStatus("nBtagLCorrUp", 1);
  T->SetBranchStatus("nBtagLCorrDown", 1);
  if(this->pName == "TT_L"){
    T->SetBranchStatus("nBtagTCorr", 1);
    T->SetBranchStatus("nBtagTCorrUp", 1);
    T->SetBranchStatus("nBtagTCorrDown", 1);
  }else{
    T->SetBranchStatus("nBtagTight", 1);
  }
  T->SetBranchStatus("N_Jets",1);
  T->SetBranchStatus("CSV",1);
  T->SetBranchStatus("pTHem1", 1);
  T->SetBranchStatus("pTHem2", 1);
  T->SetBranchStatus("etaHem1", 1);
  T->SetBranchStatus("etaHem2", 1);
  T->SetBranchStatus("phiHem1", 1);
  T->SetBranchStatus("phiHem2", 1);
  T->SetBranchStatus("ht",1);
  T->SetBranchStatus("metX",1);
  T->SetBranchStatus("metY",1);
  T->SetBranchStatus("metCorrX",1);
  T->SetBranchStatus("metCorrY",1);
  T->SetBranchStatus("N_Jets",1);
  T->SetBranchStatus("idMC", 1);
  T->SetBranchStatus("statusMC", 1);
  T->SetBranchStatus("Mu_Px",1);
  T->SetBranchStatus("Mu_Py",1);
  T->SetBranchStatus("Mu_Pz",1);
  T->SetBranchStatus("Mu_E",1);
  T->SetBranchStatus("sf_w", 1);
  T->SetBranchStatus("ISR", 1);
  T->SetBranchStatus("ISR_up", 1);
  T->SetBranchStatus("ISR_down", 1);
};

double BaseDM::HLTscale(double MR, double R2){

  int MRbin = -1;
  int R2bin = -1;

  const double R2A[] = {0.3, 0.4, 0.5, 0.6, 2.5};
  const double MRA[] = {200., 300., 400., 3500.};

  int Nbins = 3;
  int NbinsR2 = 4;

  for(int j = 0; j <= NbinsR2; j++){
    if( R2 > R2A[j]){
      if(R2 < R2A[j + 1]){
        R2bin = j+1;
        break;
      }
    }
  }
  
  for(int j = 0; j <= Nbins; j++){
    if( MR > MRA[j]){
      if(MR < MRA[j + 1]){
        MRbin = j+1;
        break;
      }
    }
  }
  
  return eff->GetEfficiency(eff->GetGlobalBin(MRbin, R2bin , 0));

};

double BaseDM::HLTscaleEle(double MR, double R2){

  int MRbin = -1;
  int R2bin = -1;

  const double R2A[] = {0.3, 0.4, 0.5, 0.6, 2.5};
  const double MRA[] = {200., 300., 400., 3500.};

  int Nbins = 3;
  int NbinsR2 = 4;

  for(int j = 0; j <= NbinsR2; j++){
    if( R2 > R2A[j]){
      if(R2 < R2A[j + 1]){
        R2bin = j+1;
        break;
      }
    }
  }

  for(int j = 0; j <= Nbins; j++){
    if( MR > MRA[j]){
      if(MR < MRA[j + 1]){
        MRbin = j+1;
        break;
      }
    }
  }

  return eff_ele->GetEfficiency(eff_ele->GetGlobalBin(MRbin, R2bin , 0));
  
};

double BaseDM::GetTriggerEff(double M_R, TString type = " "){
  int bin = -99;
  if(M_R >= 200.0 && M_R < 250.0){
    bin = 1;
  }else if(M_R >= 250.0 && M_R < 300.0){
    bin = 2;
  }else if(M_R >= 300.0){
    bin = 3;
  }
  if(type == "NoMu"){
    return to_0mu->GetBinContent(bin);
    //return to_0mu->GetBinContent(4);
  }else if(type == "Mu"){
    return to_12mu->GetBinContent(bin);
    //return to_12mu->GetBinContent(4);
  }else{
    std::cout << "INCORRECT TYPE FOR TRIGGER TURN ON; EXITING!!!!!" << std::endl;
    return -999999;
  }
}
