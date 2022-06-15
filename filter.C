void filter()
{
  gStyle->SetOptTitle(0);
  
  auto fileName = "ntuple.root";
  auto treeName = "MyLCTuple";

  TFile *myFile = new TFile("ntuple.root");
  TTree *myTree = (TTree*)myFile->Get("MyLCTuple");

  TCanvas *c = new TCanvas();
  
  //opens the file to be read
  auto openFile = TFile::Open("ntuple.root");

  //defines the histograms
  c->SetTitle("Real reconstructed particles");
  TH1F *realPassed_pt = new TH1F("h1", "MyLCTuple", 20, 0, 2000);
  TH1F *realAll_pt = new TH1F("h1", "MyLCTuple", 20, 0, 2000);
  TH1F *truthPt_Hist = new TH1F("h1", "MyLCTuple", 20, 0, 2000);
  TH1F *realPassed_azimuth = new TH1F("h1", "MyLCTuple", 20, -1.6, 1.6);
  TH1F *realAll_azimuth = new TH1F("h1", "MyLCTuple", 20, -1.6, 1.6);
  TH1F *truthAzimuth_Hist = new TH1F("h1", "MyLCTuple", 20, -1.6, 1.6);

  //creates a reader that will traverse the events of the simulation
  TTreeReader myReader("MyLCTuple", openFile);

  //arrays that hold the values of branches
  TTreeReaderArray<int> r2f_RA(myReader, "r2f");
  TTreeReaderArray<int> r2t_RA(myReader, "r2t");
  TTreeReaderArray<Float_t> mcmox_RA(myReader, "mcmox");
  TTreeReaderArray<Float_t> mcmoy_RA(myReader, "mcmoy");
  TTreeReaderArray<Float_t> mcmoz_RA(myReader, "mcmoz");
  TTreeReaderArray<Float_t> rcmox_RA(myReader, "rcmox");
  TTreeReaderArray<Float_t> rcmoy_RA(myReader, "rcmoy");
  TTreeReaderArray<Float_t> rcmoz_RA(myReader, "rcmoz");
  TTreeReaderArray<int> rctyp_RA(myReader, "rctyp");
  TTreeReaderArray<int> mcpdg_RA(myReader, "mcpdg");

  //temporary variable that will reduce the number of get operatons
  int r2fTemp, r2tTemp;
  Float_t mcmoxTemp, mcmoyTemp, mcmozTemp, rcmoxTemp, rcmoyTemp, rcmozTemp;

  //loops over each event
  while(myReader.Next()) {
    
    //loops over each array
      
    //finds which particle is real out of the reconstructed particles
    
    mcmoxTemp = mcmox_RA.At(0);
    mcmoyTemp = mcmoy_RA.At(0);
    mcmozTemp = mcmoz_RA.At(0);
	
    //calculates the transverse momentum and adds it to the histogram
    if(r2t_RA.At(0) == 0 && abs(rctyp_RA.At(0)) == 11) {
      realPassed_pt->Fill(sqrt((mcmoxTemp*mcmoxTemp) + (mcmoyTemp*mcmoyTemp)));
      realPassed_azimuth->Fill(atan(mcmozTemp / sqrt((mcmoxTemp*mcmoxTemp) + (mcmoyTemp*mcmoyTemp))));
    }
	
    realAll_pt->Fill(sqrt((mcmoxTemp*mcmoxTemp) + (mcmoyTemp*mcmoyTemp)));
    realAll_azimuth->Fill(atan(mcmozTemp / sqrt((mcmoxTemp*mcmoxTemp) + (mcmoyTemp*mcmoyTemp))));
    //truthPt_Hist->Fill(sqrt((mcmoxTemp*mcmoxTemp) + (mcmoyTemp*mcmoyTemp)));

    //calculates the polar angle and adds it to the histogram
    //truthAzimuth_Hist->Fill(atan(mcmozTemp / sqrt((mcmoxTemp*mcmoxTemp) + (mcmoyTemp*mcmoyTemp))));

    //realAll_pt->Fill(sqrt((mcmoxTemp*mcmoxTemp) + (mcmoyTemp*mcmoyTemp)));
    //realAll_azimuth->Fill()
  }

  gStyle->SetPalette(kRust);
  
  realAll_pt->Draw("PLC");
  realPassed_pt->Draw("PLC SAME");
  //reco_pt->Draw("PLC SAME");

  TLegend *legend = new TLegend(0.1, 0.1, 0.48, 0.3);
  legend->SetHeader("Transverse momenta", "C");
  legend->AddEntry("realAll_pt", "All real");
  legend->AddEntry("realPassed_pt", "Real reconstructed");
  //legend->AddEntry("reco_pt", "Reconstructed");
  legend->Draw();

  realAll_pt->GetXaxis()->SetTitle("Transverse momentum (GeV)");
  realAll_pt->GetYaxis()->SetTitle("Count");

  c->SaveAs("compLink_pt.png");

  realAll_azimuth->Draw("PLC");
  realPassed_azimuth->Draw("PLC SAME");

  legend = new TLegend(0.1, 0.1, 0.48, 0.3);
  legend->SetHeader("Azimuth", "C");
  legend->AddEntry("realAll_azimuth", "All real");
  legend->AddEntry("realPassed_azimuth", "Real reconstructed");
  legend->Draw();

  realAll_azimuth->GetXaxis()->SetTitle("Azimuth (Rads)");
  realAll_azimuth->GetYaxis()->SetTitle("Count");

  c->SaveAs("compLink_Azimuth.png");
  
  
  TEfficiency* pt_Eff = 0;
  TFile* eff_File = new TFile("effFile.root", "recreate");
  
  if(TEfficiency::CheckConsistency(*realPassed_pt, *realAll_pt)) {
    pt_Eff = new TEfficiency(*realPassed_pt, *realAll_pt);
    eff_File->Write();
  }

  c->Clear();
  
  pt_Eff->Draw();

  c->SaveAs("Eff_pt.png");

  TEfficiency* azimuth_Eff = 0;
  eff_File = new TFile("effFile.root", "recreate");

  if(TEfficiency::CheckConsistency(*realPassed_azimuth, *realAll_azimuth)) {
    azimuth_Eff = new TEfficiency(*realPassed_azimuth, *realAll_azimuth);
    eff_File->Write();
  }

  c->Clear();

  azimuth_Eff->Draw();

  c->SaveAs("Eff_azimuth.png");

  c->Clear();
  //pt_Eff->SetBins(20, 0, 200);
  //pt_Eff->Draw();
  
}
