#include "DM_1DRatioV3.hh"
#include "TMath.h"

int RatioPlots(TH1F* h1, TH1F* h2, TString h1Name = "h1Name", TString h2Name = "h2Name", TString fname = "default_name", TString type = "defaultType"){
	
  TCanvas* C = new TCanvas("C", "C	", 400, 500);
  C->cd();  
  TLegend* leg;
  
  TH1F*  RATIO;
  TString label;
  if(type == "MR"){
    RATIO = new TH1F("RATIO", fname + "_" + type , BaseDM::MR_Bins, BaseDM::MR_BinArr);
    label = "M_{R}";
    h1->GetXaxis()->SetRangeUser(200,3500);
    h2->GetXaxis()->SetRangeUser(200,3500);
    RATIO->GetXaxis()->SetRangeUser(200,3500);

  }else if(type == "RSQ" ){
    RATIO = new TH1F("RATIO", fname + "_" + type , BaseDM::RSQ_Bins, BaseDM::RSQ_BinArr);
    label = "R^{2}";
    h1->GetXaxis()->SetRangeUser(0.5, 2.5);
    h2->GetXaxis()->SetRangeUser(0.5, 2.5);
    RATIO->GetXaxis()->SetRangeUser(0.5, 2.5);

  }else{
    delete RATIO;
    delete C;
    std::cout << "Unknown Type, please use: MR or RSQ" << std::endl;
    return -1;
  }
  std::cout << "=====================Dividing Histograms=====================" << std::endl;
  RATIO->Divide(h1, h2, 1, 1, "");
  RATIO->GetYaxis()->SetRangeUser(.0, 2.05);

  h1->SetLineColor(2);
  h1->SetMarkerSize(1.);
  h1->SetMarkerColor(2);
  h1->SetMarkerStyle(20);
  h1->SetFillColor(2);
  
  h2->SetLineColor(4);
  h2->SetMarkerSize(1.);
  h2->SetMarkerColor(4);
  h2->SetMarkerStyle(20);
  h2->SetFillColor(4);
  
  h1->SetStats(0);
  h1->SetTitle("");
  h2->SetTitle("");
  h2->SetStats(0);
  h2->SetXTitle( type );
  h1->SetXTitle( type );
  
  //RATIO->GetXaxis()->SetLabelFont(63); //font in pixels
  //RATIO->GetXaxis()->SetLabelSize(16); //in pixels
  //RATIO->GetYaxis()->SetLabelFont(63); //font in pixels
  //RATIO->GetYaxis()->SetLabelSize(16); //in pixels

  //RATIO->GetYaxis()->SetTitleOffset(0.24);
  //RATIO->GetXaxis()->SetTitleOffset(0.29);

  //RATIO->GetYaxis()->SetTitleSize(0.14);
  //RATIO->GetXaxis()->SetTitleSize(0.12);

  std::cout << "GET X Title Size:  " << RATIO->GetYaxis()->GetTitleSize() << std::endl;
   
  TPad *pad1 = new TPad("pad1","pad1",0,0.25,1,1);
  pad1->SetBottomMargin(0);
  pad1->Draw();
  pad1->cd();
  h1->DrawCopy();
  h2->Draw("same");
  C->cd();
  
  leg = new TLegend(0.55, 0.7, 0.89, 0.9);//(xmin, ymin, xmax, ymax)
  leg->AddEntry(h1, label + " " + h1Name ,"f");
  leg->AddEntry(h2, label + " " + h2Name ,"f");
  //leg->SetHeader("");
  leg->SetTextSize(.022);
  leg->SetFillColor(0);
  leg->Draw();
  pad1->SetLogy();
  C->Update();
  
  TPad *pad2 = new TPad("pad2","pad2",0,0.0,1,0.25);
  pad2->SetTopMargin(0);
  pad2->SetTopMargin(0.008);
  pad2->SetBottomMargin(0.25);
  pad2->SetGridy();
  pad2->Draw();
  pad2->cd();

  RATIO->SetLineColor(4);
  RATIO->SetStats(0);
  RATIO->SetTitle("");
  //RATIO->GetXaxis()->SetLabelSize(0.05);
  //RATIO->GetYaxis()->SetLabelSize(0.05);
  //RATIO->GetYaxis()->SetTitleOffset(0.35);
  //RATIO->GetXaxis()->SetTitleOffset(0.4);
  //RATIO->GetYaxis()->SetTitleSize(0.11);
  //RATIO->GetXaxis()->SetTitleSize(0.11);
  RATIO->GetXaxis()->SetLabelSize(0.1);
  RATIO->GetYaxis()->SetLabelSize(0.08);
  RATIO->GetYaxis()->SetTitleOffset(0.35);
  RATIO->GetXaxis()->SetTitleOffset(0.88);
  RATIO->GetYaxis()->SetTitleSize(0.11);
  RATIO->GetXaxis()->SetTitleSize(0.11);
  RATIO->SetXTitle( label );
  RATIO->SetYTitle("Ratio");
  RATIO->SetLineColor(4);
  RATIO->SetMarkerSize(.7);
  RATIO->SetMarkerColor(4);
  RATIO->SetMarkerStyle(20);
  RATIO->SetFillColor(4);

  /*
  pad2->SetLogy();
  RATIO->SetLineColor(4);
  RATIO->SetStats(0);
  RATIO->SetTitle("");
  RATIO->SetLabelSize(0.2);
  RATIO->GetXaxis()->SetLabelSize(0.046);
  RATIO->GetXaxis()->SetTitleOffset(.9);
  RATIO->SetXTitle( label );
  RATIO->SetYTitle("Ratio");
  RATIO->SetLineColor(4);
  RATIO->SetMarkerSize(1.3);
  RATIO->SetMarkerColor(4);
  RATIO->SetMarkerStyle(20);
  RATIO->SetFillColor(4);
  */
 
  RATIO->Draw();
  C->cd();
  C->SaveAs(fname + ".C");
  C->SaveAs(fname + ".root");
  C->SaveAs(fname + ".pdf");
  C->SaveAs(fname + ".png");
  
  delete leg;
  delete C;
  delete RATIO;
  
  return 0;
  
};


int RatioPlotsV2(THStack* s, TH1F* h2, TH1F* h1, TString h1Name = "h1Name", TString h2Name = "h2Name", TString fname = "default_name", TString type = "defaultType", TLegend* le = NULL){

  TCanvas* C = new TCanvas("C", "C      ", 400, 500);
  C->cd();

  TH1F*  RATIO;
  TString label;
  if(type == "MR"){
    RATIO = new TH1F("RATIO", fname + "_" + type , BaseDM::MR_Bins, BaseDM::MR_BinArr);
    label = "M_{R}";
    h1->GetXaxis()->SetRangeUser(200.,3500.);
    h2->GetXaxis()->SetRangeUser(200.,3500.);
    RATIO->GetXaxis()->SetRangeUser(200.,3500.);
    RATIO->GetYaxis()->SetRangeUser(.0, 2.0);
    s->SetMaximum(500000.);
  }else if(type == "RSQ" ){
    RATIO = new TH1F("RATIO", fname + "_" + type , BaseDM::RSQ_Bins, BaseDM::RSQ_BinArr);
    label = "R^{2}";
    h1->GetXaxis()->SetRangeUser(0.5, 2.50);
    h2->GetXaxis()->SetRangeUser(0.5, 2.50);
    RATIO->GetXaxis()->SetRangeUser(0.5, 2.50);
    RATIO->GetYaxis()->SetRangeUser(.0, 2.0);
    s->SetMaximum(100000.);
  }else if(type == "MET"){
    RATIO = new TH1F("RATIO", fname + "_" + type , 20, 0, 1000);
    label = "#slash{E}_{T}  GeV";
    s->SetMaximum(10000.);
  }else if(type == "NJETS"){
    RATIO = new TH1F("RATIO", fname + "_" + type , 9, 1, 10);
    label = "Jet Multiplicity";
    s->SetMaximum(100000.);
  }else if(type=="Mass"){
    RATIO = new TH1F("RATIO", fname + "_" + type , 10, 70.0, 110.0);
    label = "M_{#mu#mu}";
    s->SetMaximum(1e4);
  }else if(type=="Angle"){
    std::cout << "========ANgle=======" << std::endl;
    RATIO = new TH1F("RATIO", fname + "_" + type , 15, .0, 2*3.1416);
    label = "#Delta#theta";
  }else if(type=="PT_J"){
    std::cout << "========PT_J=======" << std::endl;
    RATIO = new TH1F("RATIO", fname + "_" + type , 20, 80, 1000);
    label = "P^{J}_{T}";
    s->SetMaximum(3e4);
  }else if(type=="Eta_J"){
    std::cout << "========eta=======" << std::endl;
    RATIO = new TH1F("RATIO", fname + "_" + type , 20, -3.0, 3.0);
    label = "#eta";
    std::cout << "========eta1=======" << std::endl;
    s->SetMaximum(5e4);
  }else if(type=="Phi_J"){
    std::cout << "========Phi=======" << std::endl;
    RATIO = new TH1F("RATIO", fname + "_" + type , 20, -TMath::Pi(), TMath::Pi());
    label = "#phi";
    s->SetMaximum(5e4);
  }else if(type=="PT_mu"){
    std::cout << "========PT_mu=======" << std::endl;
    RATIO = new TH1F("RATIO", fname + "_" + type , 20, 0, 1000);
    label = "P^{#mu}_{T}";
    s->SetMaximum(3e4);
  }else if(type=="Eta_mu"){
    std::cout << "========eta_mu=======" << std::endl;
    RATIO = new TH1F("RATIO", fname + "_" + type , 20, -3.0, 3.0);
    label = "#eta_{#mu}";
    std::cout << "========eta1=======" << std::endl;
    s->SetMaximum(5e4);
  }else if(type=="Phi_mu"){
    std::cout << "========Phi_mu=======" << std::endl;
    RATIO = new TH1F("RATIO", fname + "_" + type , 20, -TMath::Pi(), TMath::Pi());
    label = "#phi_{#mu}";
    s->SetMaximum(5e4);
  }else if(type == "HT"){
     std::cout << "========HT=======" << std::endl;
     RATIO = new TH1F("RATIO", fname + "_" + type , 20, 0.0, 1400.0);
    label = "HT";
    s->SetMaximum(5e4);
  }else if(type == "Dphi"){
     std::cout << "========Dphi=======" << std::endl;
     RATIO = new TH1F("RATIO", fname + "_" + type , 20, -TMath::Pi(), TMath::Pi());
    label = "#Delta#phi(J_{1},J_{2})";
    s->SetMaximum(9e4);
  }else{
    delete RATIO;
    delete C;
    std::cout << "Unknown Type, please use: MR or RSQ" << std::endl;
    return -1;
  }
  
  std::cout << "=====================Dividing Histograms=====================" << std::endl;
  std::cout << h1->GetBinContent(1) << " " << h1->GetBinContent(2) << " " << h1->GetBinContent(3) << " " << h1->GetBinContent(4) << std::endl;
  std::cout << h2->GetBinContent(1) << " " << h2->GetBinContent(2) << " " << h2->GetBinContent(3) << " " << h2->GetBinContent(4) << std::endl;
  double r1 =  h2->GetBinContent(1)/h1->GetBinContent(1);
  double r2 =  h2->GetBinContent(2)/h1->GetBinContent(2);
  double r3 =  h2->GetBinContent(3)/h1->GetBinContent(3);
  double r4 =  h2->GetBinContent(4)/h1->GetBinContent(4);
  
  std::cout << r1 << " " << r2 << " " << r3 << " " << r4 << std::endl;
  RATIO->Divide(h2, h1, 1.0, 1.0, "");
  RATIO->GetYaxis()->SetRangeUser(.0, 3.05);
  h2->SetMarkerSize(.7);
  h2->SetStats(0);

  std::cout << RATIO->GetBinContent(1) << " " << RATIO->GetBinContent(2) << " " << RATIO->GetBinContent(3) << " " << RATIO->GetBinContent(4) << std::endl;
  
  std::cout << "========deb-1=======" << std::endl;
  s->SetMinimum(1.);
  
  TPad *pad1 = new TPad("pad1","pad1",0,0.25,1,1);
  pad1->SetBottomMargin(0.0);
  pad1->Draw();
  pad1->cd();
  s->Draw();
  if(type == "MET"){
    s->GetYaxis()->SetTitle("Events/20 GeV");
  }else{
    s->GetYaxis()->SetTitle("Events");
  }
  s->GetYaxis()->SetTitleOffset(1.25);
  gPad->Modified();
  h2->SetStats(0);
  h2->Draw("same");
  C->cd();
  
  le->SetFillColor(0);
  le->SetBorderSize(0);
  le->SetTextSize(0.02);
  le->Draw();
  pad1->SetLogy();
  C->Update();
  TLatex *t;
  TString n1 = "CMS Preliminary";
  TString n2 = "#sqrt{s} = 8 TeV";
  //TString n3 = "#int L dt = 19.6 fb^{-1}";
  //TString n3 = "#int L dt = 5.3 fb^{-1}";
  TString n3 = "#int L dt = 18.4 fb^{-1}";

  if(type=="Mass"){
    t = new TLatex();
    t->SetNDC();
    t->SetTextAlign(22);
    t->SetTextSize(0.03);
    t->DrawLatex(0.45,0.87, n1);
    t->DrawLatex(0.45,0.83, n2);
    t->DrawLatex(0.45,0.77, n3);
  }else if(type == "Angle"){
    TLatex *t = new TLatex();
    t->SetNDC();
    t->SetTextAlign(22);
    t->SetTextSize(0.03);
    t->DrawLatex(0.45,0.87, n1);
    t->DrawLatex(0.45,0.83, n2);
    t->DrawLatex(0.45,0.77, n3);
  }else if(type=="PT_J"){
    TLatex *t = new TLatex();
    t->SetNDC();
    t->SetTextAlign(22);
    t->SetTextSize(0.03);
    t->DrawLatex(0.45,0.87, n1);
    t->DrawLatex(0.45,0.83, n2);
    t->DrawLatex(0.45,0.77, n3);
  }else if(type=="Eta_J"){
    TLatex *t = new TLatex();
    t->SetNDC();
    t->SetTextAlign(22);
    t->SetTextSize(0.03);
    t->DrawLatex(0.45,0.87, n1);
    t->DrawLatex(0.45,0.83, n2);
    t->DrawLatex(0.45,0.77, n3);
  }else{
    TLatex *t = new TLatex();
    t->SetNDC();
    t->SetTextAlign(22);
    t->SetTextSize(0.03);
    t->DrawLatex(0.5,0.87, n1);
    t->DrawLatex(0.5,0.83, n2);
    t->DrawLatex(0.5,0.77, n3);    
  }
  
  TPad *pad2 = new TPad("pad2","pad2",0,0.0,1,0.25);
  pad2->SetTopMargin(0.008);
  pad2->SetBottomMargin(0.25);
  pad2->SetGridy();
  pad2->Draw();
  pad2->cd();
  RATIO->SetLineColor(4);
  RATIO->SetStats(0);
  RATIO->SetTitle("");
  RATIO->GetXaxis()->SetLabelSize(0.1); 
  RATIO->GetYaxis()->SetLabelSize(0.1); 
  RATIO->GetYaxis()->SetTitleOffset(0.35);
  RATIO->GetXaxis()->SetTitleOffset(0.88);
  RATIO->GetYaxis()->SetTitleSize(0.11);
  RATIO->GetXaxis()->SetTitleSize(0.11);
  RATIO->SetXTitle( label );
  RATIO->SetYTitle("Data/MC");
  RATIO->SetLineColor(4);
  RATIO->SetMarkerSize(.7);
  RATIO->SetMarkerColor(4);
  RATIO->SetMarkerStyle(20);
  RATIO->SetFillColor(4);
  RATIO->Draw();
  C->cd();
  C->SaveAs(fname + ".C");
  C->SaveAs(fname + ".root");
  C->SaveAs(fname + ".pdf");
  C->SaveAs(fname + ".png");

  delete C;
  delete RATIO;

  return 0;

};
