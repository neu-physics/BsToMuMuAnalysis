void drawROC()
{
    TString s_sig = "muonIsolation_output_zmumu.root";
    TString s_bkg = "muonIsolation_output_ttbar.root";

    TFile *f_sig = new TFile(s_sig); 
    TFile *f_bkg = new TFile(s_bkg);

    TH1F *h_sig = (TH1F*)f_sig->Get("MuonIsolationAnalyzer/h_muon_pfCandIso03_BTL");
    TH1F *h_bkg = (TH1F*)f_bkg->Get("MuonIsolationAnalyzer/h_muon_pfCandIso03_BTL");
    TH1F *h_ttbar_ptCand = (TH1F*)f_bkg->Get("MuonIsolationAnalyzer/h_muon_pT");

    int nbins = h_sig->GetNbinsX();
    float sig_integral = h_sig->Integral(1,nbins);
    float bkg_integral = h_bkg->Integral(1,nbins);
    
    std::vector<float> sigPoints(nbins);
    std::vector<float> bkgPoints(nbins);
    for ( int i = 0; i < nbins; ++i ) {
    float sig_slice_integral = h_sig->Integral(1,i+1);
    float bkg_slice_integral = h_bkg->Integral(1,i+1);
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
    g->GetXaxis()->SetRangeUser(0.8,1);
    g->GetYaxis()->SetRangeUser(0,0.35);
    g->Draw();
    /*h_sig->SetLineColor(kRed);
    h_bkg->SetLineColor(kBlack);
    h_sig->GetYaxis()->SetRangeUser(0.1,1000);
    h_sig->Draw();
    h_bkg->Draw("same");
    c1->SetLogy();*/
    TCanvas *c2 = new TCanvas("c2","c2",600,800);
    c2->cd();
    h_ttbar_ptCand->GetXaxis()->SetRangeUser(0,100);
    h_ttbar_ptCand->SetTitle("muon_pT;pT;#events");
    h_ttbar_ptCand->Draw();
}
