#include "TFile.h"
#include "string"
#include "TGraph.h"
#include "TCanvas.h"

void drawROC(string topDir="")
{
    gStyle->SetOptStat(0);
    gStyle->SetLegendBorderSize(0);
    //TString s_sig = "01-08-19_0PU/muonIsolation_output_zmumu.root";
    //TString s_bkg = "01-08-19_0PU/muonIsolation_output_ttbar.root";
    TString s_sig = ("zmumu_" + topDir+ "/muonIsolation_output_zmumu_" + topDir + ".root").c_str();
    TString s_bkg = ("ttbar_" + topDir+ "/muonIsolation_output_ttbar_" + topDir + ".root").c_str();

    TFile *f_sig = new TFile(s_sig); 
    TFile *f_bkg = new TFile(s_bkg);

    //std::cout << "signal has " << f_sig->IsOpen() << " entries" << std::endl;
    //std::cout << "background has " << f_bkg->IsOpen() << " entries" << std::endl;

    TH1F *h_sig = (TH1F*)f_sig->Get("MuonIsolationAnalyzer/h_muon_pfCandIso03_BTL");
    TH1F *h_bkg = (TH1F*)f_bkg->Get("MuonIsolationAnalyzer/h_muon_pfCandIso03_BTL");
    TH1F *h_ttbar_ptCand = (TH1F*)f_bkg->Get("MuonIsolationAnalyzer/h_muon_pT");
    TH1F *h_z_numCand = (TH1F*)f_sig->Get("MuonIsolationAnalyzer/h_muon_cutflow");
    TH1F *h_t_numCand = (TH1F*)f_bkg->Get("MuonIsolationAnalyzer/h_muon_cutflow");

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
    TCanvas *c1 = new TCanvas("c1","c1",600,600);
    c1->cd();
    c1->SetGrid();
    g->SetTitle("prompt muon ROC curve(BTL);prompt eff;non-prompt eff");
    //g->GetXaxis()->SetRangeUser(0.8,1);
    //g->GetYaxis()->SetRangeUser(0,0.35);
    g->SetLineWidth(3);
    g->Draw("AL");
    g->GetXaxis()->SetRangeUser(.84, 1.01);
    g->GetYaxis()->SetRangeUser(0, .5);

    c1->Print( ("muonIsolationROC_" + topDir + ".png").c_str() );

    /*h_sig->SetLineColor(kRed);
    h_bkg->SetLineColor(kBlack);
    h_sig->GetYaxis()->SetRangeUser(0.1,1000);
    h_sig->Draw();
    h_bkg->Draw("same");
    c1->SetLogy();*/

    /*TCanvas *c2 = new TCanvas("c2","c2",600,600);
    c2->cd();
    h_ttbar_ptCand->GetXaxis()->SetRangeUser(0,100);
    h_ttbar_ptCand->SetTitle("muon_pT;pT;#events");
    //h_ttbar_ptCand->Draw();
    h_z_numCand->SetTitle("muon cut flow;;# muons");
    h_z_numCand->SetLineColor(kBlue);
    h_t_numCand->SetLineColor(kRed);
    //h_z_numCand->GetXaxis()->SetRangeUser(0,20);
    h_t_numCand->Draw();
    h_z_numCand->Draw("same");
    
    TLegend* leg = new TLegend(0.75,0.75,0.9,0.9);
    //gStyle->SetLegendBorderSize(0);
    leg->AddEntry(h_z_numCand,"z->mumu");
    leg->AddEntry(h_t_numCand,"ttbar");
    leg->Draw();
    */
}
