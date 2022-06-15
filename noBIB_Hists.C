void noBIB_Hists()
{
  gStyle->SetOptTitle(0);

  auto fileName = "ntuple.root";
  auto treeName = "MyLCTuple";

  TFile *myFile = new TFile("ntuple.root");
  TTree *myTree = (TTree*)myFile->Get("MyLCTuple");

  TCanvas *c = new TCanvas();

  //Creating reconstructed transverse momentum
  c->SetTitle("Electron reconstructed transverse momentum");
  myTree->Draw("sqrt(rcmox*rcmox + rcmoy*rcmoy)>>reco_pt(20, 0, 2000)", "abs(rctyp) == 11", "");
  TH1 *reco_pt = (TH1*)gDirectory->Get("reco_pt");
  reco_pt->SetTitle("Electron reconstructed transverse momentum (no BIB)");
  reco_pt->GetXaxis()->SetTitle("Transverse momentum (GeV)");
  reco_pt->GetYaxis()->SetTitle("Particle count");
  c->SaveAs("reco_pt.png");

  c->Clear();

  //Creating real particle pt graph
  c->SetTitle("Electron transverse momentum");
  myTree->Draw("sqrt(mcmox*mcmox + mcmoy*mcmoy)>>real_pt(20, 0, 2000)", "abs(mcpdg) == 11 && mcgst", "");
  TH1 *real_pt = (TH1*)gDirectory->Get("real_pt");
  real_pt->SetTitle("Electron real transverse momentum (no BIB)");
  real_pt->GetXaxis()->SetTitle("Transverse momentum (GeV)");
  real_pt->GetYaxis()->SetTitle("Particle count");
  c->SaveAs("real_pt.png");

  c->Clear();

  //Creating reco azimuth graph
  c->SetTitle("Electron reconstructed azimuth");
  myTree->Draw("atan(rcmoz / sqrt(rcmox*rcmox + rcmoy*rcmoy))>>reco_azimuth(20, -1.6, 1.6)", "abs(rctyp) == 11", "");
  TH1 *reco_azimuth = (TH1*)gDirectory->Get("reco_azimuth");
  reco_azimuth->SetTitle("Electron reconstructed transverse azimuth (no BIB)");
  reco_azimuth->GetXaxis()->SetTitle("Azimuth (Rads)");
  reco_azimuth->GetYaxis()->SetTitle("Particle count");
  c->SaveAs("reco_azimuth.png");

  c->Clear();
  
  //Creating real particle azimuth graph
  c->SetTitle("Electron reconstructed azimuth");
  reco_azimuth->SetTitle("Electron reconstructed transverse azimuth (no BIB)");
  myTree->Draw("atan(mcmoz / sqrt(mcmox*mcmox + mcmoy*mcmoy))>>real_azimuth(20, -1.6, 1.6)", "abs(mcpdg) == 11 && mcgst", "");
  TH1 *real_azimuth = (TH1*)gDirectory->Get("real_azimuth");
  real_azimuth->SetTitle("Electron real transverse azimuth (no BIB)");
  real_azimuth->GetXaxis()->SetTitle("Azimuth (Rads)");
  real_azimuth->GetYaxis()->SetTitle("Particle count");
  c->SaveAs("real_azimuth.png");

  c->Clear();

  THStack *pt_Stack = new THStack("hs", "Stack of transverse momenta");

  //sets the color scheme of the graph
  gStyle->SetPalette(kRust);

  //draws the histograms using the same color scheme (the "PLC" argument)
  c->SetTitle("Combined real and reconstructed transverse momenta of electrons");
  reco_pt->Draw("PLC");
  real_pt->Draw("PLC SAME");

  //create the legend
  TLegend *legend = new TLegend(0.1, 0.1, 0.48, 0.3);
  legend->SetHeader("Transverse momenta", "C");
  legend->AddEntry("reco_pt", "Reconstructed");
  legend->AddEntry("real_pt", "Real");
  legend->Draw();

  c->SaveAs("comb_pt.png");
  c->Clear();

  //opens the file to be read
  auto openFile = TFile::Open("ntuple.root");

  //defines the histogram
  c->SetTitle("Real reconstructed particles");
  TH1F *truthPt_Hist = new TH1F("h1", "MyLCTuple", 20, 0, 2000);
  TH1F *truthAzimuth_Hist = new TH1F("h1", "MyLCTuple", 20, -1.6, 1.6);

  //creates a reader that will traverse the events of the simulation
  TTreeReader myReader("MyLCTuple", openFile);

  //arrays that hold the values of branches
  TTreeReaderArray<int> r2f_RA(myReader, "r2f");
  TTreeReaderArray<Float_t> mcmox_RA(myReader, "mcmox");
  TTreeReaderArray<Float_t> mcmoy_RA(myReader, "mcmoy");
  TTreeReaderArray<Float_t> mcmoz_RA(myReader, "mcmoz");
  TTreeReaderArray<int> rctyp_RA(myReader, "rctyp");

  //temporary variable that will reduce the number of get operatons
  int r2fTemp;
  Float_t mcmoxTemp, mcmoyTemp;

  //loops over each event
  while(myReader.Next()) {
    
    //loops over each array
    for(int i = 0; i < r2f_RA.GetSize(); i++) {
      
      //finds which particle is real out of the reconstructed particles
      r2fTemp = r2f_RA.At(i);
      
      //
      if(abs(rctyp_RA.At(r2fTemp)) == 11) {
	mcmoxTemp = mcmox_RA.At(r2fTemp);
	mcmoyTemp = mcmoy_RA.At(r2fTemp);

	//calculates the transverse momentum and adds it to the histogram
	truthPt_Hist->Fill(sqrt((mcmoxTemp*mcmoxTemp) + (mcmoyTemp*mcmoyTemp)));
	truthAzimuth_Hist->Fill(atan(mcmoz_RA.At(r2fTemp) / sqrt((mcmoxTemp*mcmoxTemp) + (mcmoyTemp*mcmoyTemp))));
      }
    }
  }

  truthPt_Hist->SetLineWidth(5);

  truthPt_Hist->Draw("PLC");
  real_pt->Draw("PLC SAME");
  //reco_pt->Draw("PLC SAME");

  legend = new TLegend(0.1, 0.1, 0.48, 0.3);
  legend->SetHeader("Transverse momenta", "C");
  legend->AddEntry("truthPt_Hist", "Real reconstructed");
  legend->AddEntry("real_pt", "All real");
  //legend->AddEntry("reco_pt", "Reconstructed");
  legend->Draw();

  truthPt_Hist->GetXaxis()->SetTitle("Transverse momentum (MeV)");
  truthPt_Hist->GetYaxis()->SetTitle("Count");

  c->SaveAs("recoLink_pt.png");

  real_azimuth->Draw("PLC");
  truthAzimuth_Hist->Draw("PLC SAME");

  legend = new TLegend(0.1, 0.1, 0.48, 0.3);
  legend->SetHeader("Azimuth", "C");
  legend->AddEntry("truthAzimuth_Hist", "Real reconstructed");
  legend->AddEntry("real_azimuth", "All real");
  legend->Draw();

  truthAzimuth_Hist->GetXaxis()->SetTitle("Azimuth (Rads)");
  truthAzimuth_Hist->GetYaxis()->SetTitle("Count");

  c->SaveAs("recoLink_Azimuth.png");
  
  TEfficiency* pt_Eff = 0;
  TFile* eff_File = new TFile("effFile.root", "recreate");
  
  if(TEfficiency::CheckConsistency(*truthPt_Hist, *real_pt)) {
    pt_Eff = new TEfficiency(*truthPt_Hist, *real_pt);
    eff_File->Write();
  }

  c->Clear();
  
  pt_Eff->Draw();

  c->SaveAs("recoEff_pt.png");

  TEfficiency* azimuth_Eff = 0;
  eff_File = new TFile("effFile.root", "recreate");

  if(TEfficiency::CheckConsistency(*truthAzimuth_Hist, *real_azimuth)) {
    azimuth_Eff = new TEfficiency(*truthAzimuth_Hist, *real_azimuth);
    eff_File->Write();
  }

  c->Clear();

  azimuth_Eff->Draw();

  c->SaveAs("recoEff_azimuth.png");
}
