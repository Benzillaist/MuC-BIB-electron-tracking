void actsckf_eff() {
  auto fileName = "actsseedckf-BIBTest.root";
  
  TFile *myFile = new TFile(fileName);

  TH1 *reco_all_pt = (TH1*)myFile->Get("MyTrackPerf/all/reco_pt");
  TH1 *reco_real_pt = (TH1*)myFile->Get("MyTrackPerf/real/reco_pt");

  //draws an efficiency graph for transverse momentum
  TEfficiency* reco_real_eff = new TEfficiency("reco_real_eff", "Efficiency of transverse momentum reconstruction;Transverse momentum (GeV);Efficiency", 20, 0, 2000);
  TFile* eff_File = new TFile("effFile.root", "recreate");

  //checks to make sure if the two histograms are compatable to make an efficiency histogram
  if(TEfficiency::CheckConsistency(*reco_real_pt, *reco_all_pt)) {
    reco_real_eff->SetPassedHistogram(*reco_real_pt, "f");
    reco_real_eff->SetTotalHistogram(*reco_all_pt, "f");
    eff_File->Write();
  }

  //creates a multigraph for the transverse momentum efficiencies
  TMultiGraph *mg = new TMultiGraph();
  TCanvas *c = new TCanvas();

  mg->Add(reco_real_eff->CreateGraph());
  
  mg->SetTitle("Transverse momentum reconstruction efficiency");
  mg->GetXaxis()->SetTitle("Transverse momentum (GeV)");
  mg->GetYaxis()->SetTitle("Efficiency");

  mg->Draw("aZ");

  c->SaveAs("actsckf_pt_eff.png");
  c->Close();
}
