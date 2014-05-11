#include "DM_KinePlotsV3.hh"
#include "DM_1DRatioV3.hh"

void CreateKinePlots(){
  
  int bL, bM, bT;
  bL = bM = 0;
  bT = 0;
  
  std::cout << "bTag Loose: " << bL << " bTag Med: " << bM << " bTag Tight: " << bT << std::endl;
  WJetsHTBins* W = new WJetsHTBins( 2 );
  ///////////////////////////////////////////
  ///////////////WJETS 1D HISTOS/////////////
  ///////////////////////////////////////////
  W->SetBtagCut(bL,bM,bT);
  W->PrintEvents();
  
  std::vector<TH1F*> Wjets = W->PlotKine();
  for(int i = 0; i < 54; i++)std::cout << i << " " << (Wjets[i])->Integral() << std::endl;

  TLegend* leg;
  std::cout << "debug 0" << std::endl;
  /////////////////////////////////////////////
  ///////////////ZJetsNuNu ANA/////////////////
  /////////////////////////////////////////////
  
  ZJetsNuNu* Z = new ZJetsNuNu( 2 );
  Z->SetBtagCut(bL,bM,bT);
  Z->PrintEvents();

  std::vector<TH1F*> Zjets = Z->PlotKine();
  for(int i = 0; i < 6; i++)std::cout << i << " " << (Zjets[i])->Integral() << std::endl;
  
  //////////////////////////////////////////////////////
  ///////////////// Drell-Yan///////////////////////////
  /////////////////////////////////////////////////////

  DY* dy = new DY( 2 );
  dy->SetBtagCut(bL,bM,bT);
  dy->PrintEvents();
  
  std::vector<TH1F*> dy_jets = dy->PlotKine();
  for(int i = 0; i < 12; i++)std::cout << i << " " << (dy_jets[i])->Integral() << std::endl;
  
  //////////////////////////////////////////////////
  /////////////////////////////////////////////////
  /////////////////TTbar + Jets////////////////////
  ////////////////////////////////////////////////
  ///////////////////////////////////////////////
  
  
  TTJets* TT = new TTJets(2);
  TT->SetBtagCut(bL,bM,bT);
  TT->PrintEvents();

  std::vector<TH1F*> TTjets = TT->PlotKine();
  for(int i = 0; i < 12; i++)std::cout << i << " " << (TTjets[i])->Integral() << std::endl;  
  
  ////////////////////////////////////////////////////////////////
  //////////////////// DATA//////////////////////////////////////
  ///////////////////////////////////////////////////////////////
  
  const char* data_file = "/media/data/cmorgoth/Data/DMData/FullHTMHTRereco/HTMHT_ABCD_FullLumi20128TeV.root";
  
  Data* data = new Data(data_file, 2);
  data->SetBtagCut(bL,bM,bT);
  data->PrintEvents();
  
  std::vector<TH1F*> data_histos = data->PlotKine();
  for(int i = 0; i < 6; i++)std::cout << i << " " << (data_histos[i])->Integral() << std::endl;
  
  TCanvas* C1 = new TCanvas("C1", "C1", 1024, 1024);
  C1->cd();
  
  ////////////////////////////////////////////////////////
  /////////////////RSQ O muon Box Stack plots////////////
  //////////////////////////////////////////////////////
  
  // THStack* stack1 = new THStack("stack1", "");
  
  for(int k = 0; k < 68; k++){
    data_histos[k]->SetMarkerStyle(20);
    data_histos[k]->SetLineColor(1);
    data_histos[k]->SetMarkerSize(1.5);
    data_histos[k]->SetTitle("");
    data_histos[k]->SetStats(0);
    
    TTjets[k]->SetFillColor(kPink+9);
    TTjets[k]->SetTitle("");
    TTjets[k]->SetStats(0);
    dy_jets[k]->SetFillColor(kViolet+9);
    dy_jets[k]->SetTitle("");
    dy_jets[k]->SetStats(0);
    Wjets[k]->SetFillColor(kYellow-4);
    Wjets[k]->SetTitle("");
    Wjets[k]->SetStats(0);
    
    Zjets[k]->SetFillColor(kSpring+4);
    Zjets[k]->SetTitle("");
    Zjets[k]->SetStats(0);
  }
  
  TH1F* aux2;
  THStack* stack1;
  
  for(int m = 0; m < 18; m++){
    stack1 = new THStack("stack1", "");
    leg = new TLegend(0.7,0.7,0.9,0.92);
    
    leg->AddEntry(Wjets[18+m],"W + jets","f");
    leg->AddEntry(Zjets[18+m],"Z(#nu#bar{#nu}) + jets","f");
    leg->AddEntry(TTjets[18+m],"t #bar{t} + jets","f");
    leg->AddEntry(dy_jets[18+m],"Z/#gamma^{*}(ll) + jets","f");
    leg->AddEntry(data_histos[m],"Data","lep");
    
    stack1->Add(dy_jets[18+m]);//DY
    stack1->Add(TTjets[18+m]);//TTbar
    stack1->Add(Wjets[18+m]);//Wjets
    stack1->Add(Zjets[18+m]);//ZJets
    
    TH1F* aux2 = new TH1F( *dy_jets[m] );
    aux2->Sumw2();
    aux2->Add(TTjets[m], 1);
    aux2->Add(Zjets[m], 1);
    aux2->Add(Wjets[m], 1);
    
    TString type;
    TString name;
    if(m == 0 || m == 3 || m == 6 || m == 9 || m == 12 || m == 15){
      type = "PT_J";
      switch(m){
      case 0:
	name = "Kine/PT_J1_BOX0_Stack";
	break;
      case 3:
	name = "Kine/PT_J2_BOX0_Stack";
	break;
      case 6:
	name = "Kine/PT_J1_BOX1_Stack";
	break;
      case 9:
	name = "Kine/PT_J2_BOX1_Stack";
	break;
      case 12:
	name = "Kine/PT_J1_BOX2_Stack";
	break;
      case 15:
	name = "Kine/PT_J2_BOX2_Stack";
	break;
      default:
	std::cout << "value of x unknown" << std::endl;
      }
    }else if(m == 1 || m == 4 || m == 7 || m == 10 || m == 13 || m == 16){
      type = "Eta_J";
      switch(m){
      case 1:
	name = "Kine/Eta_J1_BOX0_Stack";
	break;
      case 4:
	name = "Kine/Eta_J2_BOX0_Stack";
	break;
      case 7:
	name = "Kine/Eta_J1_BOX1_Stack";
	break;
      case 10:
	name = "Kine/Eta_J2_BOX1_Stack";
	break;
      case 13:
	name = "Kine/Eta_J1_BOX2_Stack";
	break;
      case 16:
	name = "Kine/Eta_J2_BOX2_Stack";
	break;
      default:
	std::cout << "value of x unknown" << std::endl;
      }
    }else if(m == 2 || m == 5 || m == 8 || m == 11 || m == 14 || m == 17){
      type = "Phi_J";
      switch(m){
      case 2:
	name = "Kine/Phi_J1_BOX0_Stack";
	break;
      case 5:
	name = "Kine/Phi_J2_BOX0_Stack";
	break;
      case 8:
	name = "Kine/Phi_J1_BOX1_Stack";
	break;
      case 11:
	name = "Kine/Phi_J2_BOX1_Stack";
	break;
      case 14:
	name = "Kine/Phi_J1_BOX2_Stack";
	break;
      case 17:
	name = "Kine/Phi_J2_BOX2_Stack";
	break;
      default:
	std::cout << "value of x unknown" << std::endl;
      }
    }
    

    RatioPlotsV2(stack1, data_histos[m] , aux2, "MC 0 #mu BOX", "Data 0 #mu BOX", name, type, leg);
    delete leg;
    delete aux2;
    delete stack1;
  }
  
  for(int m = 36; m < 45; m++){
    
    stack1 = new THStack("stack1", "");
    leg = new TLegend(0.7,0.7,0.9,0.92);
    
    leg->AddEntry(Wjets[m+9],"W + jets","f");
    leg->AddEntry(Zjets[m+9],"Z(#nu#bar{#nu}) + jets","f");
    leg->AddEntry(TTjets[m+9],"t #bar{t} + jets","f");
    leg->AddEntry(dy_jets[m+9],"Z/#gamma^{*}(ll) + jets","f");
    leg->AddEntry(data_histos[m],"Data","lep");
    
    stack1->Add(dy_jets[m+9]);//DY
    stack1->Add(TTjets[m+9]);//TTbar
    stack1->Add(Wjets[m+9]);//Wjets
    stack1->Add(Zjets[m+9]);//ZJets
    
    TH1F* aux2 = new TH1F( *dy_jets[m] );
    aux2->Sumw2();
    aux2->Add(TTjets[m], 1);
    aux2->Add(Zjets[m], 1);
    aux2->Add(Wjets[m], 1);
    
    TString type;
    TString name;
    if(m == 36 || m == 39 || m == 42){
      type = "PT_mu";
      switch(m){
      case 36:
	name = "Kine/PT_mu1_BOX2_Stack";
	break;
      case 39:
	name = "Kine/PT_mu2_BOX2_Stack";
	break;
      case 42:
	name = "Kine/PT_mu1_BOX1_Stack";
	break;
      default:
	std::cout << "value of x unknown" << std::endl;
      }
    }else if(m == 37 || m == 40 || m == 43){
      type = "Eta_mu";
      switch(m){
      case 37:
	name = "Kine/Eta_mu1_BOX2_Stack";
	break;
      case 40:
	name = "Kine/Eta_mu2_BOX2_Stack";
	break;
      case 43:
	name = "Kine/Eta_mu1_BOX1_Stack";
	break;
      default:
	std::cout << "value of x unknown" << std::endl;
      }
    }else if(m == 38 || m == 41 || m == 44){
      type = "Phi_mu";
      switch(m){
      case 38:
	name = "Kine/Phi_mu1_BOX2_Stack";
	break;
      case 41:
	name = "Kine/Phi_mu2_BOX2_Stack";
	break;
      case 44:
	name = "Kine/Phi_mu1_BOX1_Stack";
	break;
      default:
	std::cout << "value of x unknown" << std::endl;
      }
    }
    
    RatioPlotsV2(stack1, data_histos[m] , aux2, "MC 0 #mu BOX", "Data 0 #mu BOX", name, type, leg);
    delete leg;
    delete aux2;
    delete stack1;
  }


  //Double Mu Invariant Mass
  stack1 = new THStack("stack1", "");
  leg = new TLegend(0.7,0.7,0.9,0.92);
  
  leg->AddEntry(Wjets[54],"W + jets","f");
  leg->AddEntry(Zjets[54],"Z(#nu#bar{#nu}) + jets","f");
  leg->AddEntry(TTjets[54],"t #bar{t} + jets","f");
  leg->AddEntry(dy_jets[54],"Z/#gamma^{*}(ll) + jets","f");
  leg->AddEntry(data_histos[55],"Data","lep");
  
  stack1->Add(dy_jets[54]);//DY
  stack1->Add(TTjets[54]);//TTbar
  stack1->Add(Wjets[54]);//Wjets
  stack1->Add(Zjets[54]);//ZJets
  
  aux2 = new TH1F( *dy_jets[55] );
  aux2->Sumw2();
  aux2->Add(TTjets[55], 1);
  aux2->Add(Zjets[55], 1);
  aux2->Add(Wjets[55], 1);

  RatioPlotsV2(stack1, data_histos[55] , aux2, "MC 2 #mu BOX", "Data 2 #mu BOX", "Kine/DoubleMuInvMass_70_100", "Mass", leg);
  delete leg;
  delete aux2;
  delete stack1;

  for(int m = 56; m < 59; m++){
    stack1 = new THStack("stack1", "");
    leg = new TLegend(0.7,0.7,0.9,0.92);
    
    leg->AddEntry(Wjets[m+3],"W + jets","f");
    leg->AddEntry(Zjets[m+3],"Z(#nu#bar{#nu}) + jets","f");
    leg->AddEntry(TTjets[m+3],"t #bar{t} + jets","f");
    leg->AddEntry(dy_jets[m+3],"Z/#gamma^{*}(ll) + jets","f");
    leg->AddEntry(data_histos[m],"Data","lep");
    
    stack1->Add(dy_jets[m+3]);//DY
    stack1->Add(TTjets[m+3]);//TTbar
    stack1->Add(Wjets[m+3]);//Wjets
    stack1->Add(Zjets[m+3]);//ZJets
    
    TH1F* aux2 = new TH1F( *Wjets[m] );
    //aux2->Sumw2();
    aux2->Add(TTjets[m], 1);
    aux2->Add(Zjets[m], 1);
    aux2->Add(dy_jets[m], 1);
    
    TString type("HT");
    TString name;
    TString lbl;
    switch(m){
    case 56:
      name = "Kine/HT_BOX0_Stack";
      lbl = "0#mu";
      break;
    case 57:
      name = "Kine/HT_BOX1_Stack";
      lbl = "1#mu";
      break;
    case 58:
      name = "Kine/HT_BOX2_Stack";
      lbl = "2#mu";
      break;
    default:
      std::cout << "value of x unknown" << std::endl;
      break;
    }
    
    RatioPlotsV2(stack1, data_histos[m] , aux2, "MC "+lbl, "Data "+lbl, name, type, leg);
    delete leg;
    delete aux2;
    delete stack1;
  }

  for(int m = 62; m < 65; m++){
    stack1 = new THStack("stack1", "");
    leg = new TLegend(0.7,0.7,0.9,0.92);
    
    leg->AddEntry(Wjets[m+3],"W + jets","f");
    leg->AddEntry(Zjets[m+3],"Z(#nu#bar{#nu}) + jets","f");
    leg->AddEntry(TTjets[m+3],"t #bar{t} + jets","f");
    leg->AddEntry(dy_jets[m+3],"Z/#gamma^{*}(ll) + jets","f");
    leg->AddEntry(data_histos[m],"Data","lep");
    
    stack1->Add(dy_jets[m+3]);//DY
    stack1->Add(TTjets[m+3]);//TTbar
    stack1->Add(Wjets[m+3]);//Wjets
    stack1->Add(Zjets[m+3]);//ZJets
    
    TH1F* aux2 = new TH1F( *Wjets[m] );
    //aux2->Sumw2();
    aux2->Add(TTjets[m], 1);
    aux2->Add(Zjets[m], 1);
    aux2->Add(dy_jets[m], 1);
    
    TString type("Dphi");
    TString name;
    TString lbl;
    switch(m){
    case 62:
      name = "Kine/Dphi_BOX0_Stack";
      lbl = "0#mu";
      break;
    case 63:
      name = "Kine/Dphi_BOX1_Stack";
      lbl = "1#mu";
      break;
    case 64:
      name = "Kine/Dphi_BOX2_Stack";
      lbl = "2#mu";
      break;
    default:
      std::cout << "value of x unknown" << std::endl;
      break;
    }
    
    RatioPlotsV2(stack1, data_histos[m] , aux2, "MC "+lbl, "Data "+lbl, name, type, leg);
    delete leg;
    delete aux2;
    delete stack1;
  }

}
