#include "TFile.h"
#include "string"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"

TGraph* returnROC(string topDir, string subdet, string pfCandIso)
{
  TString s_sig = ("zmumu_" + topDir+ "/muonIsolation_output_zmumu_" + topDir + ".root").c_str();
  TString s_bkg = ("ttbar_" + topDir+ "/muonIsolation_output_ttbar_" + topDir + ".root").c_str();
  
  TFile *f_sig = new TFile(s_sig); 
  TFile *f_bkg = new TFile(s_bkg);
  
  //std::cout << "signal has " << f_sig->IsOpen() << " entries" << std::endl;
  //std::cout << "background has " << f_bkg->IsOpen() << " entries" << std::endl;
  
  TH1F *h_sig = (TH1F*)f_sig->Get( ("MuonIsolationAnalyzer/h_muon_" + pfCandIso + "_" + subdet).c_str() );
  TH1F *h_bkg = (TH1F*)f_bkg->Get( ("MuonIsolationAnalyzer/h_muon_" + pfCandIso + "_" + subdet).c_str() );
  
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

void drawThreeCurves( string topDir,  string subdet = "BTL", string pfCandIso = "pfCandIso03")
{
  
  TGraph *g_0PU   = returnROC( ("0PU_" + topDir).c_str(), subdet, pfCandIso );
  TGraph *g_200PU = returnROC( ("200PU_" + topDir).c_str(), subdet, pfCandIso );
  //TGraph *g_TDR   = returnROC( ("MTDTDR_200PU_" + topDir).c_str(), subdet );
  
  TCanvas *c1 = new TCanvas("c1","c1",600,600);
  c1->cd();
  c1->SetGrid();
  g_0PU->SetTitle( ("prompt muon ROC curve (" + subdet + ");prompt eff;non-prompt eff").c_str() );
  //g->GetXaxis()->SetRangeUser(0.8,1);
  //g->GetYaxis()->SetRangeUser(0,0.35);
  g_0PU->SetLineWidth(3);
  g_0PU->Draw("AL");
  g_0PU->GetXaxis()->SetRangeUser(.85, 1.01);
  g_0PU->GetYaxis()->SetRangeUser(0, .1);
  
  g_200PU->SetLineColor(kRed);
  g_200PU->SetLineWidth(3);
  g_200PU->Draw("L same");
  
  //g_TDR->SetLineColor(kBlue);
  //g_TDR->SetLineWidth(3);
  //g_TDR->Draw("L same");

  TLegend *leg = new TLegend(0.2, 0.8, .35, .9);
  leg->AddEntry(g_0PU, "No MTD, No PU", "l");
  leg->AddEntry(g_200PU, "No MTD, 200PU", "l");
  //leg->AddEntry(g_200PU_dt30, "MTD, dz < 1 mm, dt < 3#sigma (#sigma = 30ps)", "l");
  leg->Draw("same");
  
  c1->Print( ("muonIsolationROC_multi_" + topDir + "_" + pfCandIso + "_" + subdet + ".png").c_str() );
  
}

void drawThreeCurveWithTiming( string topDir,  string subdet = "BTL", string useDxy = "")
{
  
  TGraph *g_0PU_nom    = returnROC( ("MTDTDR_0PU_" + topDir).c_str(), subdet, ("pfCandIso03" + useDxy).c_str() );
  TGraph *g_200PU_nom  = returnROC( ("MTDTDR_200PU_" + topDir).c_str(), subdet, ("pfCandIso03" + useDxy).c_str() );
  TGraph *g_200PU_dt30 = returnROC( ("MTDTDR_200PU_" + topDir).c_str(), subdet, ("pfCandIso03_dt30" + useDxy).c_str() );
  //TGraph *g_TDR   = returnROC( ("MTDTDR_200PU_" + topDir).c_str(), subdet );
  
  TCanvas *c1 = new TCanvas("c1","c1",600,600);
  c1->cd();
  c1->SetGrid();
  g_0PU_nom->SetTitle( ("prompt muon ROC curve (" + subdet + ");prompt eff;non-prompt eff").c_str() );

  g_0PU_nom->SetLineWidth(3);
  g_0PU_nom->SetLineStyle(9);
  g_0PU_nom->Draw("AL");
  g_0PU_nom->GetXaxis()->SetRangeUser(.85, 1.01);
  g_0PU_nom->GetYaxis()->SetRangeUser(0, .1);
  
  g_200PU_nom->SetLineColor(kBlue);
  g_200PU_nom->SetLineWidth(3);
  g_200PU_nom->Draw("L same");

  g_200PU_dt30->SetLineColor(kRed);
  g_200PU_dt30->SetLineWidth(3);
  g_200PU_dt30->Draw("L same");
  
  TLegend *leg = new TLegend(0.15, 0.65, .4, .85);
  leg->AddEntry(g_0PU_nom, "No MTD, No PU", "l");
  leg->AddEntry(g_200PU_nom, "No MTD, 200PU", "l");
  leg->AddEntry(g_200PU_dt30, "MTD, dz < 1 mm, dt < 3#sigma (#sigma = 30ps)", "l");
  leg->SetTextSize(0.025);
  leg->Draw("same");
  
  c1->Print( ("muonIsolationROC_TDR" + useDxy + "_" + topDir + "_" + subdet + ".png").c_str() );
  
}

void drawMultipleROC(string topDir="", string subdet = "BTL")
{
    gStyle->SetOptStat(0);
    gStyle->SetLegendBorderSize(0);
    
    drawThreeCurves(topDir, "BTL", "pfCandIso03");
    drawThreeCurves(topDir, "ETL", "pfCandIso03");
    drawThreeCurves(topDir, "BTL", "pfCandIso03_noDxy");
    drawThreeCurves(topDir, "ETL", "pfCandIso03_noDxy");

    drawThreeCurves(topDir, "BTL", "pfCandIso03_dt30");
    drawThreeCurves(topDir, "ETL", "pfCandIso03_dt30");
    drawThreeCurves(topDir, "BTL", "pfCandIso03_dt30_noDxy");
    drawThreeCurves(topDir, "ETL", "pfCandIso03_dt30_noDxy");
    
    drawThreeCurveWithTiming( topDir,  "BTL");
    drawThreeCurveWithTiming( topDir,  "ETL");
    drawThreeCurveWithTiming( topDir,  "BTL", "_noDxy");
    drawThreeCurveWithTiming( topDir,  "ETL", "_noDxy");

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
