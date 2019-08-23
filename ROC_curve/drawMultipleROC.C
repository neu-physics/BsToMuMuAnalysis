#include "TFile.h"
#include "string"
#include "TGraph.h"
#include "TCanvas.h"

TGraph* returnROC(string topDir, string subdet = "BTL")
{

  TString s_sig = ("zmumu_" + topDir+ "/muonIsolation_output_zmumu_" + topDir + ".root").c_str();
  TString s_bkg = ("ttbar_" + topDir+ "/muonIsolation_output_ttbar_" + topDir + ".root").c_str();
  
  TFile *f_sig = new TFile(s_sig); 
  TFile *f_bkg = new TFile(s_bkg);
  
  //std::cout << "signal has " << f_sig->IsOpen() << " entries" << std::endl;
  //std::cout << "background has " << f_bkg->IsOpen() << " entries" << std::endl;
  
  TH1F *h_sig = (TH1F*)f_sig->Get( ("MuonIsolationAnalyzer/h_muon_pfCandIso03_" + subdet).c_str() );
  TH1F *h_bkg = (TH1F*)f_bkg->Get( ("MuonIsolationAnalyzer/h_muon_pfCandIso03_" + subdet).c_str() );
  TH1F *h_ttbar_ptCand = (TH1F*)f_bkg->Get("MuonIsolationAnalyzer/h_muon_pT" );
  TH1F *h_z_numCand = (TH1F*)f_sig->Get("MuonIsolationAnalyzer/h_muon_cutflow" );
  TH1F *h_t_numCand = (TH1F*)f_bkg->Get("MuonIsolationAnalyzer/h_muon_cutflow" );
  
  int nbins = h_sig->GetNbinsX();
  float sig_integral = h_sig->Integral(1,nbins);
  float bkg_integral = h_bkg->Integral(1,nbins);
  
  std::vector<float> sigPoints(nbins);
  std::vector<float> bkgPoints(nbins);
  int firstBin = 1;
  
  for ( int i = firstBin; i < nbins; ++i ) {
    float sig_slice_integral = h_sig->Integral(firstBin,i+1);
    float bkg_slice_integral = h_bkg->Integral(firstBin,i+1);
    sigPoints.push_back(sig_slice_integral/sig_integral);
    bkgPoints.push_back(bkg_slice_integral/bkg_integral);
    //std::cout << "sig:" << sig_slice_integral << std::endl;
    //std::cout << "bkg:" << bkg_slice_integral << std::endl;
  }
  TGraph *g = new TGraph(sigPoints.size(),&sigPoints[0],&bkgPoints[0]);
 
  return g;
}

void drawThreeCurves( string topDir,  string subdet = "BTL")
{

    TGraph *g_0PU   = returnROC( ("0PU_" + topDir).c_str(), subdet );
    TGraph *g_200PU = returnROC( ("200PU_" + topDir).c_str(), subdet );
    TGraph *g_TDR   = returnROC( ("MTDTDR_200PU_" + topDir).c_str(), subdet );

    TCanvas *c1 = new TCanvas("c1","c1",600,600);
    c1->cd();
    c1->SetGrid();
    g_0PU->SetTitle("prompt muon ROC curve(BTL);prompt eff;non-prompt eff");
    //g->GetXaxis()->SetRangeUser(0.8,1);
    //g->GetYaxis()->SetRangeUser(0,0.35);
    g_0PU->SetLineWidth(3);
    g_0PU->Draw("AL");
    g_0PU->GetXaxis()->SetRangeUser(.85, 1.01);
    g_0PU->GetYaxis()->SetRangeUser(0, .1);

    g_200PU->SetLineColor(kRed);
    g_200PU->SetLineWidth(3);
    g_200PU->Draw("L same");

    g_TDR->SetLineColor(kBlue);
    g_TDR->SetLineWidth(3);
    g_TDR->Draw("L same");

    c1->Print( ("muonIsolationROC_multi_" + topDir + "_" + subdet + ".png").c_str() );

}

void drawMultipleROC(string topDir="", string subdet = "BTL")
{
    gStyle->SetOptStat(0);
    gStyle->SetLegendBorderSize(0);

    drawThreeCurves(topDir, "BTL");
    drawThreeCurves(topDir, "ETL");

    /*
    TGraph *g_0PU   = returnROC( ("0PU_" + topDir).c_str(), subdet );
    TGraph *g_200PU = returnROC( ("200PU_" + topDir).c_str(), subdet );
    TGraph *g_TDR   = returnROC( ("MTDTDR_200PU_" + topDir).c_str(), subdet );

    
    TCanvas *c1 = new TCanvas("c1","c1",600,600);
    c1->cd();
    c1->SetGrid();
    g_0PU->SetTitle("prompt muon ROC curve(BTL);prompt eff;non-prompt eff");
    //g->GetXaxis()->SetRangeUser(0.8,1);
    //g->GetYaxis()->SetRangeUser(0,0.35);
    g_0PU->SetLineWidth(3);
    g_0PU->Draw("AL");
    g_0PU->GetXaxis()->SetRangeUser(.85, 1.01);
    g_0PU->GetYaxis()->SetRangeUser(0, .1);

    g_200PU->SetLineColor(kRed);
    g_200PU->SetLineWidth(3);
    g_200PU->Draw("L same");

    g_TDR->SetLineColor(kBlue);
    g_TDR->SetLineWidth(3);
    g_TDR->Draw("L same");

    c1->Print( ("muonIsolationROC_multi_" + topDir + ".png").c_str() );
    */
}
