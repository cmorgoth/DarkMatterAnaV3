#include "DM_DataV3.hh"
#include <iostream>
#include <iomanip>
#include "TLorentzVector.h"
#include "TROOT.h"

Data::Data(){ /*empty*/};

Data::Data(const char* FileName, int MetIndex ):BaseDM(MetIndex, "Data"){
  gROOT->Reset();
  std::cout << "====DeBUG DATA====" << std::endl;
  T = new TChain("outTree");
  //T->Add("/media/data/cmorgoth/Data/DMData/DataFinal/Run2012A_06Aug2012v1_total.root");
  //T->Add("/media/data/cmorgoth/Data/DMData/DataFinal/Run2012A_13Jul2012_total.root");
  //T->Add("/media/data/cmorgoth/Data/DMData/DataFinal/NewTrigger/RunAB_NewTrigger.root");
  T->Add("/media/data/cmorgoth/Data/DMData/DataFinal/HTMHT_Parked/HTMHTParked_ILV_RunB_Tot.root");
  T->Add("/media/data/cmorgoth/Data/DMData/DataFinal/HTMHT_Parked/HTMHTParked_ILV_RunC_Tot.root");
  T->Add("/media/data/cmorgoth/Data/DMData/DataFinal/HTMHT_Parked/HTMHTParked_ILV_RunD_Tot.root");
  
  effT = new TChain("effTree");
  //effT->Add("/media/data/cmorgoth/Data/DMData/DataFinal/Run2012A_06Aug2012v1_total.root");
  //effT->Add("/media/data/cmorgoth/Data/DMData/DataFinal/Run2012A_13Jul2012_total.root");
  //effT->Add("/media/data/cmorgoth/Data/DMData/DataFinal/NewTrigger/RunAB_NewTrigger.root");
  effT->Add("/media/data/cmorgoth/Data/DMData/DataFinal/HTMHT_Parked/HTMHTParked_ILV_RunB_Tot.root");
  effT->Add("/media/data/cmorgoth/Data/DMData/DataFinal/HTMHT_Parked/HTMHTParked_ILV_RunC_Tot.root");
  effT->Add("/media/data/cmorgoth/Data/DMData/DataFinal/HTMHT_Parked/HTMHTParked_ILV_RunD_Tot.root");
  
  std::cout << "ENTRIES: " <<  T->GetEntries() << std::endl;
  std::cout << "====Adding File====> " << FileName << std::endl;
  std::cout << "====Btag Index====> " << btagIndex << std::endl;
};

Data::~Data(){
  delete T;
};


bool Data::PrintEvents(){

  double NtotGen = 0, Nt_PV = 0, Nt_2J = 0, Nt_0b = 0, Nt_LepVeto = 0, N_In, N_PV, N_2J, N_0b, N_LepVeto;

  
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
  std::cout << "================ Data ==============" << std::endl;
  std::cout << "========================================" << std::endl;
  
  std::cout << "Data Nt_In: " << NtotGen << std::endl;
  std::cout << "Data Nt_PV: " << Nt_PV << std::endl;
  std::cout << "Data Nt_2J: " << Nt_2J << std::endl;
  std::cout << "Data Nt_0b: " << Nt_0b << std::endl;
  std::cout << "Data Nt_LepVeto: " << Nt_LepVeto << std::endl;
  
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
  
  for(int i = 0; i < T->GetEntries(); i++){
    T->GetEntry(i);
    if( RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin ){
      Nt_MR_RSQ_cut++;
      
      if( BOX == 0 )Nt_0muBox++;
      if( BOX == 1 )Nt_1muBox++;
      if( BOX == 2 )Nt_2muBox++;
      
      if( nBtag == 0 ){
	
	Nt_MR_RSQ_cut0BTag++;
	
	if( BOX == 0 )Nt_0muBox0BTag++;
	if( BOX == 1 )Nt_1muBox0BTag++;
	if( BOX == 2 )Nt_2muBox0BTag++;
	
      }
      
    }
    
  }
  T->SetBranchStatus("*", 0);//Disable all branches, trying to gain some performance

  
  std::cout << "Data weighted Nt_MR_RSQ_cut: " << Nt_MR_RSQ_cut << std::endl;
  std::cout << "Data weighted Nt_0muBox: " << Nt_0muBox << std::endl;
  std::cout << "Data weighted Nt_1muBox: " << Nt_1muBox << std::endl;
  std::cout << "Data weighted Nt_2muBox: " << Nt_2muBox << std::endl;

  std::cout << "Data weighted Nt_MR_RSQ_cut0Btag: " << Nt_MR_RSQ_cut0BTag << std::endl;
  std::cout << "Data weighted Nt_0muBox0Btag: " << Nt_0muBox0BTag << std::endl;
  std::cout << "Data weighted Nt_1muBox0Btag: " << Nt_1muBox0BTag << std::endl;
  std::cout << "Data weighted Nt_2muBox0Btag: " << Nt_2muBox0BTag << std::endl;

  std::cout << "Data Btag EVENTS: " << Nt_MR_RSQ_cut - Nt_MR_RSQ_cut0BTag << std::endl;
  std::cout << "========================================" << std::endl;
  std::cout << "========================================" << "\n\n" << std::endl;
  
  return true;
  
  
};

TH1F Data::PlotMR_1Box(){
  
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
  
  for(int i = 0; i < T->GetEntries(); i++){
    T->GetEntry(i);
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= nBtagCut[0]);
    fBtag[3] = (nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    fBtag[4] = (nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1]);

    if( BOX==1 && RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin && fBtag[btagIndex] ){
      MR1->Fill(MR[metIndex]);
    }
  }
  T->SetBranchStatus("*", 0);
    
  return *MR1;
};

TH1F Data::PlotMR_0Box(){
  
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
  
  for(int i = 0; i < T->GetEntries(); i++){
    T->GetEntry(i);
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= nBtagCut[0]);
    fBtag[3] = (nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1]);
    
    if( BOX == 0 && RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin && fBtag[btagIndex]){
      MR0->Fill(MR[metIndex]);
    }
  }
  T->SetBranchStatus("*", 0);
    
  return *MR0;
};


TH1F  Data::PlotRSQ_1Box(){
  
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
  
  for(int i = 0; i < T->GetEntries(); i++){
    T->GetEntry(i);
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= nBtagCut[0]);
    fBtag[3] = (nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1]);
    
    if( BOX == 1 && RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin && fBtag[btagIndex]){
      RSQ1->Fill(RSQ[metIndex]);
    }
  }
  T->SetBranchStatus("*", 0);
    
  return *RSQ1;
    
};


TH1F  Data::PlotRSQ_0Box(){
  
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
  
  for(int i = 0; i < T->GetEntries(); i++){
    T->GetEntry(i);
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= nBtagCut[0]);
    fBtag[3] = (nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1]);

    if( BOX == 0 && RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin && fBtag[btagIndex]){
      RSQ0->Fill(RSQ[metIndex]);
    }
  }
  T->SetBranchStatus("*", 0);

  return *RSQ0;
  
};

TH2F Data::PlotRSQ_vs_MR_0Box(){
  
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
  
  for(int i = 0; i < T->GetEntries(); i++){
    T->GetEntry(i);
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= nBtagCut[0]);
    fBtag[3] = (nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1]);
    
    if( BOX == 0 && RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin && fBtag[btagIndex]){
      RSQ_MR_0BOX->Fill(MR[metIndex], RSQ[metIndex]);
    }
  }
  T->SetBranchStatus("*", 0);

  return *RSQ_MR_0BOX;
  
};

TH2F Data::PlotRSQ_vs_MR_1Box(){
  
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
  
  for(int i = 0; i < T->GetEntries(); i++){
    T->GetEntry(i);
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= nBtagCut[0]);
    fBtag[3] = (nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1]);
    
    if( BOX == 1 && RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin && fBtag[btagIndex]){
      RSQ_MR_1BOX->Fill(MR[metIndex], RSQ[metIndex]);
    }
  }
  T->SetBranchStatus("*", 0);
  
  return *RSQ_MR_1BOX;
  
};


std::vector<TH1F*> Data::PlotMETmag(){
  double metX[4], metcorrX[4], metY[4], metcorrY[4], ht, RSQ[4], MR[4], run/*, ls, evNum*/;
  double mht[3], CSV[30], pu_w, mu_w, sf_w;
  int BOX, nBtag[2], N_Jets;
  double pTHem1, pTHem2, etaHem1, etaHem2, phiHem1, phiHem2;
  double Mu_Px[2], Mu_Py[2], Mu_Pz[2], Mu_E[2];
 
  std::vector< TH1F* > metvec;
  TH1F* MET[12];
  TString name;
  for(int l = 0; l < 3; l++ ){
    for( int m = 0; m < 3; m++ ){
      name = TString(Form("Data_METmag_Box%d_plotType%d",l,m));
      MET[3*l + m] = new TH1F( name, name, 20, 0, 1000 );
    }
  }
  
  MET[9] = new TH1F( "NJETS0_W", "NJETS 0 BOX", 9, 1, 10);
  MET[10] = new TH1F( "NJETS1_W", "NJETS 1 BOX", 9, 1, 10);
  MET[11] = new TH1F( "NJETS2_W", "NJETS 2 BOX", 9, 1, 10);
  
  SetMetStatus();
  T->SetBranchAddress("RSQ", RSQ);
  T->SetBranchAddress("MR", MR);
  T->SetBranchAddress("BOX_NUM", &BOX);
  T->SetBranchAddress("nBtag", &nBtag[0]);
  T->SetBranchAddress("nBtagTight", &nBtag[1]);
  T->SetBranchAddress("N_Jets", &N_Jets);
  T->SetBranchAddress("CSV", CSV);
  T->SetBranchAddress("ht", &ht);
  T->SetBranchAddress("mht", &mht[0]);
  T->SetBranchAddress("metX", metX);
  T->SetBranchAddress("metCorrX", metcorrX);
  T->SetBranchAddress("metY", metY);
  T->SetBranchAddress("metCorrY", metcorrY);
  T->SetBranchAddress("N_Jets", &N_Jets);
  
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
  
  float metmag = .0;
  float metmagcorr = .0;
  double wt = 1.;
  for(int i = 0; i < T->GetEntries(); i++){
    T->GetEntry(i);
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= nBtagCut[0]);
    fBtag[3] = (nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    //fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1]);
    fBtag[4] = (nBtag[1] >= nBtagCut[2]);
    
    TLorentzVector j1;
    TLorentzVector j2, double_mu;

    j1.SetPtEtaPhiE(pTHem1, etaHem1, phiHem1, pTHem1*cosh(etaHem1));//Hemisphere
    j2.SetPtEtaPhiE(pTHem2, etaHem2, phiHem2, pTHem2*cosh(etaHem2));//Hemisphere
    double Dphi = j1.DeltaPhi(j2);

    metmag = sqrt(metX[metIndex]*metX[metIndex]+metY[metIndex]*metY[metIndex]);
    metmagcorr = sqrt(metcorrX[metIndex]*metcorrX[metIndex]+metcorrY[metIndex]*metcorrY[metIndex]);
    if( RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin  && fBtag[btagIndex]  && fabs(Dphi) < 2.5){
      if( BOX == 0 ){
	MET[0]->Fill(metmag);
	MET[1]->Fill(metmagcorr);
	MET[2]->Fill(metmagcorr-metmag);
	MET[9]->Fill(N_Jets);
      }else if( BOX == 1 ){
	MET[3]->Fill(metmag);
	MET[4]->Fill(metmagcorr);
	MET[5]->Fill(metmagcorr-metmag);
	MET[10]->Fill(N_Jets);
      }else if( BOX == 2 ){
	double_mu.SetPxPyPzE(Mu_Px[0]+Mu_Px[1], Mu_Py[0]+Mu_Py[1], Mu_Pz[0]+Mu_Pz[1], Mu_E[0]+Mu_E[1]);
	
	if(double_mu.M() > 75.0 && double_mu.M() < 110.0){
	  MET[6]->Fill(metmag);
	  MET[7]->Fill(metmagcorr);
	  MET[8]->Fill(metmagcorr-metmag);
	  MET[11]->Fill(N_Jets);
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

std::vector<TH2F*> Data::Plot_2DRazor(){
  double RSQ[4], MR[4], CSV[30];
  double pTHem1, pTHem2, etaHem1, etaHem2, phiHem1, phiHem2, pu_w, mu_w, sf_w;
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
  T->SetBranchAddress("RSQ", RSQ);
  T->SetBranchAddress("MR", MR);
  T->SetBranchAddress("BOX_NUM", &BOX);
  T->SetBranchAddress("nBtag", &nBtag[0]);
  T->SetBranchAddress("nBtagTight", &nBtag[1]);
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
	Razor2D[0]->Fill(MR[metIndex], RSQ[metIndex]);
      }else if( BOX == 1 ){
	Razor2D[1]->Fill(MR[metIndex], RSQ[metIndex]);
      }else if( BOX == 2 ){
	double_mu.SetPxPyPzE(Mu_Px[0]+Mu_Px[1], Mu_Py[0]+Mu_Py[1], Mu_Pz[0]+Mu_Pz[1], Mu_E[0]+Mu_E[1]);
	
	if(double_mu.M() > 75.0 && double_mu.M() < 110.0){
	  Razor2D[2]->Fill(MR[metIndex], RSQ[metIndex]);
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

bool Data::pfJetPassCSVM(double btagOutput){
  if(btagOutput < 0.679)   return false;
  return true;
};

int Data::pfJetPassCSVM(double* CSVM, int N_Jets){
  int nMBtag = 0;
  for(int i = 0; i < N_Jets; i++)if(CSVM[i] >= 0.679)nMBtag++;
  return nMBtag;
};

std::vector<TH1F*> Data::Plot_1DRazor(){
  double RSQ[4], MR[4], CSV[30];
  double pTHem1, pTHem2, etaHem1, etaHem2, phiHem1, phiHem2, pu_w, mu_w, sf_w;
  double Mu_Px[2], Mu_Py[2], Mu_Pz[2], Mu_E[2];
  int BOX, N_Jets, nBtag[2];

  std::vector< TH1F* > Razor1DVec;
  TH1F* Razor1D[12];
  TString name, name1;
  
  for(int l = 0; l < 6; l++){
    name = TString(Form("MR_1D_Data_%dmu_Box",l));                                      
    name1 = TString(Form("R2_1D_Data_%dmu_Box",l));
    Razor1D[2*l] = new TH1F(name,name, MR_Bins, MR_BinArr);
    Razor1D[2*l+1] = new TH1F(name1,name1, RSQ_Bins, RSQ_BinArr);
    if(l < 3){
      Razor1D[2*l]->Sumw2();
      Razor1D[2*l+1]->Sumw2();
    }
  }

  SetStatus();
  T->SetBranchAddress("RSQ", RSQ);
  T->SetBranchAddress("MR", MR);
  T->SetBranchAddress("BOX_NUM", &BOX);
  T->SetBranchAddress("nBtag", &nBtag[0]);
  T->SetBranchAddress("nBtagTight", &nBtag[1]);
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
      if(BOX == 0){
	Razor1D[0]->Fill(MR[metIndex]);
	Razor1D[1]->Fill(RSQ[metIndex]);
	Razor1D[6]->Fill(MR[metIndex]);
        Razor1D[7]->Fill(RSQ[metIndex]);
      }else if(BOX == 1){
	Razor1D[2]->Fill(MR[metIndex]);
        Razor1D[3]->Fill(RSQ[metIndex]);
        Razor1D[8]->Fill(MR[metIndex]);
        Razor1D[9]->Fill(RSQ[metIndex]);
      }else if(BOX == 2){
	double_mu.SetPxPyPzE(Mu_Px[0]+Mu_Px[1], Mu_Py[0]+Mu_Py[1], Mu_Pz[0]+Mu_Pz[1], Mu_E[0]+Mu_E[1]);
	
	if(double_mu.M() > 75.0 && double_mu.M() < 110.0){
	  Razor1D[4]->Fill(MR[metIndex]);
	  Razor1D[5]->Fill(RSQ[metIndex]);
	  Razor1D[10]->Fill(MR[metIndex]);
	  Razor1D[11]->Fill(RSQ[metIndex]);
	}
	
      }
    }

  }
  
  for(int j = 0; j < 12; j++){
    Razor1DVec.push_back(Razor1D[j]);
  }
  
  return Razor1DVec;
  
};

std::vector<TH1F*> Data::Plot_MRCat(){
 double RSQ[4], MR[4], CSV[30];
  double pTHem1, pTHem2, etaHem1, etaHem2, phiHem1, phiHem2;
  double Mu_Px[2], Mu_Py[2], Mu_Pz[2], Mu_E[2];
  
  int BOX, N_Jets, nBtag[2];
  int c1_bins = 11;
  float c1B[] = {0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.9, 0.95, 1.0, 2.5};
  int c2_bins = 6;
  float c2B[] = {0.50, 0.575, 0.65, 0.75, 0.85, .950, 2.5};
  int c3_bins = 6;
  float c3B[] = {0.50, 0.575, 0.65, 0.75, 0.85, .950, 2.5};
  int c4_bins = 4;
  float c4B[] = {0.50, 0.60, 0.70, .950, 2.50};
    
  std::vector< TH1F* > Razor1DVec;
  TH1F* Razor1D[12];
  TString name1, name2, name3, name4;
  for(int l = 0; l < 3; l++){
    name1 = TString(Form(this->pName+"_cat1_1D_%dmu_Box",l));
    name2 = TString(Form(this->pName+"_cat2_1D_%dmu_Box",l));
    name3 = TString(Form(this->pName+"_cat3_1D_%dmu_Box",l));
    name4 = TString(Form(this->pName+"_cat4_1D_%dmu_Box",l));
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
  
  T->SetBranchAddress("RSQ", RSQ);
  T->SetBranchAddress("MR", MR);
  T->SetBranchAddress("BOX_NUM", &BOX);
  T->SetBranchAddress("nBtag", &nBtag[0]);
  //T->SetBranchAddress("nBtagLCorr", &nBtag[0]);
  T->SetBranchAddress("nBtagTight", &nBtag[1]);
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
  
  
  for(int i = 0; i < T->GetEntries(); i++){
    T->GetEntry(i);
    TLorentzVector j1, double_mu;
    TLorentzVector j2;
    
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
      if(BOX == 0){
	if(MR[metIndex] > 200.0 && MR[metIndex] <= 300.0 ){
	  Razor1D[0]->Fill(RSQ[metIndex]);
	}else if(MR[metIndex] > 300.0 && MR[metIndex] <= 400.0){
	  Razor1D[1]->Fill(RSQ[metIndex]);
	}else if(MR[metIndex] > 400.0 && MR[metIndex] <= 600.0){
	  Razor1D[2]->Fill(RSQ[metIndex]);
	}else if(MR[metIndex] > 600.0 && MR[metIndex] <= 3500.0){
	  Razor1D[3]->Fill(RSQ[metIndex]);
	}
      }else if(BOX == 1){
	if(MR[metIndex] > 200.0 && MR[metIndex] <= 300.0 ){
	  Razor1D[4]->Fill(RSQ[metIndex]);
	}else if(MR[metIndex] > 300.0 && MR[metIndex] <= 400.0){
	  Razor1D[5]->Fill(RSQ[metIndex]);
	}else if(MR[metIndex] > 400.0 && MR[metIndex] <= 600.0){
	  Razor1D[6]->Fill(RSQ[metIndex]);
	}else if(MR[metIndex] > 600.0 && MR[metIndex] <= 3500.0){
	  Razor1D[7]->Fill(RSQ[metIndex]);
	}
      }else if(BOX == 2){
	double_mu.SetPxPyPzE(Mu_Px[0]+Mu_Px[1], Mu_Py[0]+Mu_Py[1], Mu_Pz[0]+Mu_Pz[1], Mu_E[0]+Mu_E[1]);
	//std::cout << "-----" << double_mu.M() << std::endl;
	if(double_mu.M() > 75.0 && double_mu.M() < 110.0){	
	  if(MR[metIndex] > 200.0 && MR[metIndex] <= 300.0 ){
	    Razor1D[8]->Fill(RSQ[metIndex]);
	  }else if(MR[metIndex] > 300.0 && MR[metIndex] <= 400.0){
	    Razor1D[9]->Fill(RSQ[metIndex]);
	  }else if(MR[metIndex] > 400.0 && MR[metIndex] <= 600.0){
	    Razor1D[10]->Fill(RSQ[metIndex]);
	  }else if(MR[metIndex] > 600.0 && MR[metIndex] <= 3500.0){
	    Razor1D[11]->Fill(RSQ[metIndex]);
	  }
	}
	
      }
    }

  }

  for(int j = 0; j < 12; j++){
    Razor1DVec.push_back(Razor1D[j]);
  }

  return Razor1DVec;
  
}

std::vector<TH1F*> Data::PlotKine(){
  double RSQ[4], MR[4], CSV[30];
  double Jet_PT[30], Jet_Eta[30], Jet_Phi[30], Mu_Px[2], Mu_Py[2], Mu_Pz[2], Mu_E[2];
  double pTHem1, pTHem2, etaHem1, etaHem2, phiHem1, phiHem2;
  int BOX, N_Jets, nBtag[2];

  std::vector< TH1F* > Razor1DKVec;
  TH1F* Razor1DK[56];
  TString name[6];
  double hltWeight;
  
  for(int l = 0; l < 6; l++){
    name[0] = TString(Form("PT_1D_J1_Data_%dmu_Box",l));
    name[1] = TString(Form("ETA_1D_J1_Data_%dmu_Box",l));
    name[2] = TString(Form("PHI_1D_J1_Data_%dmu_Box",l));
    name[3] = TString(Form("PT_1D_J2_Data_%dmu_Box",l));
    name[4] = TString(Form("ETA_1D_J2_Data_%dmu_Box",l));
    name[5] = TString(Form("PHI_1D_J2_Data_%dmu_Box",l));

    Razor1DK[6*l] = new TH1F(name[0],name[0], 20, 80, 1000);
    Razor1DK[6*l+1] = new TH1F(name[1],name[1], 20, -3.0, 3.0);
    Razor1DK[6*l+2] = new TH1F(name[2],name[2], 20, -TMath::Pi(), TMath::Pi());
    Razor1DK[6*l+3] = new TH1F(name[3],name[3], 20, 80, 1000);
    Razor1DK[6*l+4] = new TH1F(name[4],name[4], 20, -3.0, 3.0);
    Razor1DK[6*l+5] = new TH1F(name[5],name[5], 20, -TMath::Pi(), TMath::Pi());

    if(l < 3){
      for(int k = 0; k < 6; k++)Razor1DK[6*l+k]->Sumw2();
    }
  }

  for(int l = 0; l < 2; l++){
    int box_mu = 2;
    name[0] = TString(Form("PT_1D_m1_Data_%dmu_Box_v%d",box_mu, l));
    name[1] = TString(Form("ETA_1D_m1_Data_%dmu_Box_v%d",box_mu, l));
    name[2] = TString(Form("PHI_1D_m1_Data_%dmu_Box_v%d",box_mu, l));
    name[3] = TString(Form("PT_1D_m2_Data_%dmu_Box_v%d",box_mu, l));
    name[4] = TString(Form("ETA_1D_m2_Data_%dmu_Box_v%d",box_mu, l));
    name[5] = TString(Form("PHI_1D_m2_Data_%dmu_Box_v%d",box_mu, l));

    Razor1DK[36+9*l] = new TH1F(name[0],name[0], 20, 0.0, 1000);
    Razor1DK[36+9*l+1] = new TH1F(name[1],name[1], 20, -3.0, 3.0);
    Razor1DK[36+9*l+2] = new TH1F(name[2],name[2], 20, -TMath::Pi(), TMath::Pi());
    Razor1DK[36+9*l+3] = new TH1F(name[3],name[3], 20, 0.0, 1000);
    Razor1DK[36+9*l+4] = new TH1F(name[4],name[4], 20, -3.0, 3.0);
    Razor1DK[36+9*l+5] = new TH1F(name[5],name[5], 20, -TMath::Pi(), TMath::Pi());

    box_mu = 1;
    name[0] = TString(Form("PT_1D_m1_Data_%dmu_Box_v%d",box_mu, l));
    name[1] = TString(Form("ETA_1D_m1_Data_%dmu_Box_v%d",box_mu, l));
    name[2] = TString(Form("PHI_1D_m1_Data_%dmu_Box_v%d",box_mu, l));

    Razor1DK[36+9*l+6] = new TH1F(name[0],name[0], 20, 0.0, 1000);
    Razor1DK[36+9*l+7] = new TH1F(name[1],name[1], 20, -3.0, 3.0);
    Razor1DK[36+9*l+8] = new TH1F(name[2],name[2], 20, -TMath::Pi(), TMath::Pi());

    if(l < 1){
      for(int k = 0; k < 9; k++)Razor1DK[36+9*l+k]->Sumw2();
    }
  }
  
  Razor1DK[54] = new TH1F("diMuMass", "diMuMass", 10, 75.0, 110.0);
  Razor1DK[55] = new TH1F("diMuMass_err", "diMuMass_err", 10, 75.0, 110.0);
  Razor1DK[55]->Sumw2();
  
  SetStatusKine();
  
  T->SetBranchAddress("RSQ", RSQ);
  T->SetBranchAddress("MR", MR);
  T->SetBranchAddress("BOX_NUM", &BOX);
  T->SetBranchAddress("nBtag", &nBtag[0]);
  T->SetBranchAddress("nBtagTight", &nBtag[1]);
  T->SetBranchAddress("N_Jets", &N_Jets);
  T->SetBranchAddress("Jet_PT", &Jet_PT);
  T->SetBranchAddress("Jet_Eta", &Jet_Eta);
  T->SetBranchAddress("Jet_Phi", &Jet_Phi);
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
  
  for(int i = 0; i < T->GetEntries(); i++){
    T->GetEntry(i);
    TLorentzVector j1, j2;//Megajet
    TLorentzVector mu1, mu2, double_mu;

    j1.SetPtEtaPhiE(pTHem1, etaHem1, phiHem1, pTHem1*cosh(etaHem1));//Hemisphere1                          
    j2.SetPtEtaPhiE(pTHem2, etaHem2, phiHem2, pTHem2*cosh(etaHem2));//Hemisphere2

    double Dphi = j1.DeltaPhi(j2);
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= nBtagCut[0]);
    fBtag[3] = (nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    //fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1]);
    fBtag[4] = (nBtag[1] >= nBtagCut[2]);
    
    if( RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin  && fBtag[btagIndex] && fabs(Dphi) < 2.5){
      if(BOX == 0){
	Razor1DK[0]->Fill(Jet_PT[0]);
        Razor1DK[1]->Fill(Jet_Eta[0]);
        Razor1DK[2]->Fill(Jet_Phi[0]);
        Razor1DK[3]->Fill(Jet_PT[1]);
	Razor1DK[4]->Fill(Jet_Eta[1]);
	Razor1DK[5]->Fill(Jet_Phi[1]);
	
	Razor1DK[18]->Fill(Jet_PT[0]);
        Razor1DK[18+1]->Fill(Jet_Eta[0]);
        Razor1DK[18+2]->Fill(Jet_Phi[0]);
        Razor1DK[18+3]->Fill(Jet_PT[1]);
        Razor1DK[18+4]->Fill(Jet_Eta[1]);
        Razor1DK[18+5]->Fill(Jet_Phi[1]);
      }else if(BOX == 1){
	//Jet Plots
	Razor1DK[6]->Fill(Jet_PT[0]);
        Razor1DK[7]->Fill(Jet_Eta[0]);
        Razor1DK[8]->Fill(Jet_Phi[0]);
        Razor1DK[9]->Fill(Jet_PT[1]);
        Razor1DK[10]->Fill(Jet_Eta[1]);
        Razor1DK[11]->Fill(Jet_Phi[1]);
	//Replicas
	Razor1DK[18+6]->Fill(Jet_PT[0]);
        Razor1DK[18+7]->Fill(Jet_Eta[0]);
        Razor1DK[18+8]->Fill(Jet_Phi[0]);
        Razor1DK[18+9]->Fill(Jet_PT[1]);
        Razor1DK[18+10]->Fill(Jet_Eta[1]);
        Razor1DK[18+11]->Fill(Jet_Phi[1]);
	
	//Mu Plots
	mu1.SetPxPyPzE(Mu_Px[0], Mu_Py[0], Mu_Pz[0], Mu_E[0]);
	
	Razor1DK[36+6]->Fill(mu1.Pt());
	Razor1DK[36+7]->Fill(mu1.Eta());
	Razor1DK[36+8]->Fill(mu1.Phi());
	//Replicas
	Razor1DK[36+15]->Fill(mu1.Pt());
	Razor1DK[36+16]->Fill(mu1.Eta());
	Razor1DK[36+17]->Fill(mu1.Phi());
	
      }else if(BOX == 2){
	double_mu.SetPxPyPzE(Mu_Px[0]+Mu_Px[1], Mu_Py[0]+Mu_Py[1], Mu_Pz[0]+Mu_Pz[1], Mu_E[0]+Mu_E[1]);
	
	if(double_mu.M() > 75.0 && double_mu.M() < 110.0){
	  Razor1DK[54]->Fill(double_mu.M());
	  Razor1DK[55]->Fill(double_mu.M());
	  
	  Razor1DK[12]->Fill(Jet_PT[0]);
	  Razor1DK[13]->Fill(Jet_Eta[0]);
	  Razor1DK[14]->Fill(Jet_Phi[0]);
	  Razor1DK[15]->Fill(Jet_PT[1]);
	  Razor1DK[16]->Fill(Jet_Eta[1]);
	  Razor1DK[17]->Fill(Jet_Phi[1]);
	  //Replicas
	  Razor1DK[18+12]->Fill(Jet_PT[0]);
	  Razor1DK[18+13]->Fill(Jet_Eta[0]);
	  Razor1DK[18+14]->Fill(Jet_Phi[0]);
	  Razor1DK[18+15]->Fill(Jet_PT[1]);
	  Razor1DK[18+16]->Fill(Jet_Eta[1]);
	  Razor1DK[18+17]->Fill(Jet_Phi[1]);
	  
	  //Mu Plots              
	  mu1.SetPxPyPzE(Mu_Px[0], Mu_Py[0], Mu_Pz[0], Mu_E[0]);
	  mu2.SetPxPyPzE(Mu_Px[1], Mu_Py[1], Mu_Pz[1], Mu_E[1]);
	  
	  Razor1DK[36]->Fill(mu1.Pt());
	  Razor1DK[36+1]->Fill(mu1.Eta());
	  Razor1DK[36+2]->Fill(mu1.Phi());
	  Razor1DK[36+3]->Fill(mu2.Pt());
	  Razor1DK[36+4]->Fill(mu2.Eta());
	  Razor1DK[36+5]->Fill(mu2.Phi());
	  //Replicas
	  Razor1DK[36+9]->Fill(mu1.Pt());
	  Razor1DK[36+10]->Fill(mu1.Eta());
	  Razor1DK[36+11]->Fill(mu1.Phi());
	  Razor1DK[36+12]->Fill(mu2.Pt());
	  Razor1DK[36+13]->Fill(mu2.Eta());
	  Razor1DK[36+14]->Fill(mu2.Phi());
	}
	
      }
    }

  }
  
  for(int j = 0; j < 56; j++){
    Razor1DKVec.push_back(Razor1DK[j]);
  }
  
  return Razor1DKVec;

};
  

bool Data::SetStatus(){

  T->SetBranchStatus("*",0); //disable all branches
  
  T->SetBranchStatus("BOX_NUM", 1);
  T->SetBranchStatus("RSQ", 1);
  T->SetBranchStatus("MR", 1);
  T->SetBranchStatus("pTHem1", 1);
  T->SetBranchStatus("etaHem1", 1);
  T->SetBranchStatus("phiHem1", 1);
  T->SetBranchStatus("pTHem2", 1);
  T->SetBranchStatus("etaHem2", 1);
  T->SetBranchStatus("phiHem2", 1);
  T->SetBranchStatus("N_Jets", 1);
  T->SetBranchStatus("CSV", 1);
  T->SetBranchStatus("nBtag", 1);
  T->SetBranchStatus("nBtagTight", 1);
  T->SetBranchStatus("Mu_Px",1);
  T->SetBranchStatus("Mu_Py",1);
  T->SetBranchStatus("Mu_Pz",1);
  T->SetBranchStatus("Mu_E",1);
  
  return true;
  
};

bool Data::SetStatusKine(){
  T->SetBranchStatus("*",0); //disable all branches
  T->SetBranchStatus("RSQ",1);
  T->SetBranchStatus("MR",1);
  T->SetBranchStatus("BOX_NUM",1);
  T->SetBranchStatus("nBtag",1);
  T->SetBranchStatus("nBtagTight",1);
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
};

bool Data::SetMetStatus(){
  T->SetBranchStatus("*",0); //disable all branches                                                                                                     
  T->SetBranchStatus("RSQ",1);
  T->SetBranchStatus("MR",1);
  T->SetBranchStatus("BOX_NUM",1);
  T->SetBranchStatus("nBtag",1);
  T->SetBranchStatus("nBtagTight",1);
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
  T->SetBranchStatus("Mu_Px",1);
  T->SetBranchStatus("Mu_Py",1);
  T->SetBranchStatus("Mu_Pz",1);
  T->SetBranchStatus("Mu_E",1);
  //T->SetBranchStatus("idMC", 1);
  //T->SetBranchStatus("statusMC", 1);

  return true;
  
};
