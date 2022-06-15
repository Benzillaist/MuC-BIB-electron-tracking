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
  TH1F *realPassed_pt = new TH1F("rP_pt", "MyLCTuple", 20, 0, 2000);
  TH1F *realAll_pt = new TH1F("rA_pt", "MyLCTuple", 20, 0, 2000);
  TH1F *truthPt_Hist = new TH1F("a_pt", "MyLCTuple", 20, 0, 2000);
  TH1F *realPassed_azimuth = new TH1F("rP_a", "MyLCTuple", 20, -1.6, 1.6);
  TH1F *realAll_azimuth = new TH1F("rA_a", "MyLCTuple", 20, -1.6, 1.6);
  TH1F *truthAzimuth_Hist = new TH1F("a_a", "MyLCTuple", 20, -1.6, 1.6);

  TH1D *realWeights_Hist = new TH1D("rW_H", "MyLCTuple", 20, 0, 1);

  //creates a reader that will traverse the events of the simulation
  TTreeReader myReader("MyLCTuple", openFile);

  //arrays that hold the values of branches
  TTreeReaderArray<int> r2f_RA(myReader, "r2f");
  TTreeReaderArray<int> r2t_RA(myReader, "r2t");
  TTreeReaderArray<Float_t> r2w_RA(myReader, "r2w");
  TTreeReaderArray<Float_t> mcmox_RA(myReader, "mcmox");
  TTreeReaderArray<Float_t> mcmoy_RA(myReader, "mcmoy");
  TTreeReaderArray<Float_t> mcmoz_RA(myReader, "mcmoz");
  TTreeReaderArray<Float_t> rcmox_RA(myReader, "rcmox");
  TTreeReaderArray<Float_t> rcmoy_RA(myReader, "rcmoy");
  TTreeReaderArray<Float_t> rcmoz_RA(myReader, "rcmoz");
  TTreeReaderArray<int> rctyp_RA(myReader, "rctyp");
  TTreeReaderArray<int> mcpdg_RA(myReader, "mcpdg");
  TTreeReaderArray<int> mcgst_RA(myReader, "mcgst");

  //temporary variable that will reduce the number of get operatons
  int r2fTemp, r2tTemp;
  Float_t mcmoxTemp, mcmoyTemp, mcmozTemp, rcmoxTemp, rcmoyTemp, rcmozTemp, r2wTemp;
  Float_t r2wMax = 0;
  int r2wMax_Index = -1;
  int count = 0;

  //loops over each event
  while(myReader.Next()) {
    
    //loops over each array
      
    //finds which particle is real out of the reconstructed particles
	
    //calculates the transverse momentum and adds it to the histogram
    for(int i = 0; i < r2f_RA.GetSize(); i++) {
      r2tTemp = r2t_RA.At(i);
      r2fTemp = r2f_RA.At(i);
      r2wTemp = r2w_RA.At(i);
      
      if(abs(mcpdg_RA.At(r2tTemp)) == 11 && abs(rctyp_RA.At(r2fTemp)) == 11 && mcgst_RA.At(r2tTemp)) {
	if(r2wTemp > r2wMax) {
	  r2wMax = r2wTemp;
	  r2wMax_Index = i;
	}
	realWeights_Hist->Fill(r2wTemp);
      }
	
    //truthPt_Hist->Fill(sqrt((mcmoxTemp*mcmoxTemp) + (mcmoyTemp*mcmoyTemp)));

    //calculates the polar angle and adds it to the histogram
    //truthAzimuth_Hist->Fill(atan(mcmozTemp / sqrt((mcmoxTemp*mcmoxTemp) + (mcmoyTemp*mcmoyTemp))));

    //realAll_pt->Fill(sqrt((mcmoxTemp*mcmoxTemp) + (mcmoyTemp*mcmoyTemp)));
    //realAll_azimuth->Fill()
    }

    if(r2wMax_Index != -1) {
      r2tTemp = r2t_RA.At(r2wMax_Index);
      r2fTemp = r2f_RA.At(r2wMax_Index);
      r2wTemp = r2w_RA.At(r2wMax_Index);
      
      mcmoxTemp = mcmox_RA.At(r2tTemp);
      mcmoyTemp = mcmoy_RA.At(r2tTemp);
      mcmozTemp = mcmoz_RA.At(r2tTemp);
    
      cout << "count: " << count << " instance: " << r2wMax_Index << " r2f: " << r2fTemp << " r2t: " << r2tTemp << " rctyp: " << rctyp_RA.At(r2fTemp) <<  " mcpdg: " << mcpdg_RA.At(r2tTemp) << " r2w: " << r2wTemp << endl;
      realPassed_pt->Fill(sqrt((mcmoxTemp*mcmoxTemp) + (mcmoyTemp*mcmoyTemp)));
      realPassed_azimuth->Fill(atan(mcmozTemp / sqrt((mcmoxTemp*mcmoxTemp) + (mcmoyTemp*mcmoyTemp))));
    }
    
    for(int i = 0; i < mcmox_RA.GetSize(); i++) {
      if(mcgst_RA.At(i)) {
	mcmoxTemp = mcmox_RA.At(i);
	mcmoyTemp = mcmoy_RA.At(i);
	mcmozTemp = mcmoz_RA.At(i);
	realAll_pt->Fill(sqrt((mcmoxTemp*mcmoxTemp) + (mcmoyTemp*mcmoyTemp)));
	realAll_azimuth->Fill(atan(mcmozTemp / sqrt((mcmoxTemp*mcmoxTemp) + (mcmoyTemp*mcmoyTemp))));
      }
    }
    count++;
    r2wMax = 0;
    r2wMax_Index = -1;
  }

  gStyle->SetPalette(kRust);
  
  realAll_pt->Draw("PLC");
  realPassed_pt->Draw("PLC SAME");

  TLegend *legend = new TLegend(0.1, 0.1, 0.48, 0.3);
  legend->SetHeader("Transverse momenta", "C");
  legend->AddEntry("realAll_pt", "All real");
  legend->AddEntry("realPassed_pt", "Real reconstructed");
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

  c->Clear();
  
  realWeights_Hist->Draw();
  realWeights_Hist->GetXaxis()->SetTitle("Weight");
  realWeights_Hist->GetYaxis()->SetTitle("Count");
  
  c->SaveAs("realWeights.png");

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
