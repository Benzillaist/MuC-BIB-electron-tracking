void filter()
{ 
  auto fileName = "ntuple.root";
  auto treeName = "MyLCTuple";

  TFile *myFile = new TFile("ntuple.root");
  TTree *myTree = (TTree*)myFile->Get("MyLCTuple");

  TCanvas *c = new TCanvas();
  
  //opens the file to be read
  auto openFile = TFile::Open("ntuple.root");

  //defines the histograms
  c->SetTitle("Real reconstructed particles");
  TH1F *realPassed_pt = new TH1F("rP_pt", "Linked electrons", 20, 0, 2000);
  TH1F *realAll_pt = new TH1F("rA_pt", "All electrons", 20, 0, 2000);
  TH1F *realPassed_PA = new TH1F("rP_PA", "Linked electrons", 20, 0, 1.6);
  TH1F *realAll_PA = new TH1F("rA_PA", "All electrons", 20, 0, 1.6);
  TH1F *realPassed_azimuth = new TH1F("rP_a", "Linked electrons", 20, -1.6, 1.6);
  TH1F *realAll_azimuth = new TH1F("rA_a", "All electrons", 20, -1.6, 1.6);

  TH1D *realWeights_Hist = new TH1D("realWeights", "MyLCTuple", 26, 0, 1.3);
  TH1F *diffPt_Hist = new TH1F("diffPt", "MyLCTuple", 40, -500, 1500);
  TH1F *deltaR_Hist = new TH1F("deltaR", "MyLCTuple", 25, 0, 0.0025);

  TH1F *electronReco_pt = new TH1F("elRec_pt", "Electron pt", 20, 0, 2000);
  TH1F *photonReco_pt = new TH1F("phRec_pt", "Photon pt", 20, 0, 2000);
  TH1F *neutronReco_pt = new TH1F("nuRec_pt", "Neutron pt", 20, 0, 2000);
  TH1F *electronReco_PA = new TH1F("elRec_PA", "Electron PA", 20, 0, 1.6);
  TH1F *photonReco_PA = new TH1F("phRec_PA", "Photon PA", 20, 0, 1.6);
  TH1F *neutronReco_PA = new TH1F("nuRec_PA", "Neutron PA", 20, 0, 1.6);
  TH1F *electronReco_azimuth = new TH1F("elRec_a", "Electron PA", 20, -1.6, 1.6);
  TH1F *photonReco_azimuth = new TH1F("phRec_a", "Photon PA", 20, -1.6, 1.6);
  TH1F *neutronReco_azimuth = new TH1F("nuRec_a", "Neutron PA", 20, -1.6, 1.6);

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

  //temporary variables that will reduce the number of get operatons
  int r2fTemp, r2tTemp;
  Float_t mcmoxTemp, mcmoyTemp, mcmozTemp, rcmoxTemp, rcmoyTemp, rcmozTemp, r2wTemp, mPtTemp, mPATemp, mATemp, rPtTemp, rPATemp, rATemp, deltaRTemp, dPATemp, dATemp;
  Float_t r2wMax = 0;
  int r2wMax_Index = -1;
  int count = 0;

  /////////////////////
  //HISTOGRAM FILLING//
  /////////////////////
  
  //loops over each event
  while(myReader.Next()) {

    //searches for the largest weight of the relationship that links a reconstructed particle to a particle that we know all the details of
    for(int i = 0; i < r2f_RA.GetSize(); i++) {
      r2tTemp = r2t_RA.At(i);
      r2fTemp = r2f_RA.At(i);
      r2wTemp = r2w_RA.At(i);

      //checks to see if those particles are electrons and if the paticle is an originial generating particle
      if(abs(mcpdg_RA.At(r2tTemp)) == 11 && abs(rctyp_RA.At(r2fTemp)) == 11 && mcgst_RA.At(r2tTemp)) {
	if(r2wTemp > r2wMax && r2wTemp > 0.5) {
	  r2wMax = r2wTemp;
	  r2wMax_Index = i;
	}
	//adds the weight to a histogram
	realWeights_Hist->Fill(r2wTemp);
      }
	
    }

    //if a relationship was found between those particles, it is added to histograms
    if(r2wMax_Index != -1) {
      r2tTemp = r2t_RA.At(r2wMax_Index);
      r2fTemp = r2f_RA.At(r2wMax_Index);
      r2wTemp = r2w_RA.At(r2wMax_Index);
      
      mcmoxTemp = mcmox_RA.At(r2tTemp);
      mcmoyTemp = mcmoy_RA.At(r2tTemp);
      mcmozTemp = mcmoz_RA.At(r2tTemp);
      rcmoxTemp = rcmox_RA.At(r2fTemp);
      rcmoyTemp = rcmoy_RA.At(r2fTemp);
      rcmozTemp = rcmoz_RA.At(r2fTemp);

      mPtTemp = sqrt((mcmoxTemp*mcmoxTemp) + (mcmoyTemp*mcmoyTemp));
      mPATemp = atan2(mPtTemp, mcmozTemp);
      mATemp = atan2(mcmoxTemp, mcmoyTemp);

      rPtTemp = sqrt((rcmoxTemp*rcmoxTemp) + (rcmoyTemp*rcmoyTemp));
      rPATemp = atan2(rPtTemp, rcmozTemp);
      rATemp = atan2(rcmoxTemp, rcmoyTemp);

      dPATemp = mPATemp - rPATemp;
      dATemp = mATemp - rATemp;

      realPassed_pt->Fill(mPtTemp);
      realPassed_PA->Fill(mPATemp);
      realPassed_azimuth->Fill(mATemp);

      diffPt_Hist->Fill(mPtTemp - rPtTemp);
      deltaR_Hist->Fill(sqrt((dPATemp*dPATemp) + (dATemp*dATemp)));
    }

    //loops over all particles to find the generating particle
    for(int i = 0; i < mcmox_RA.GetSize(); i++) {
      if(mcgst_RA.At(i)) {
	mcmoxTemp = mcmox_RA.At(i);
	mcmoyTemp = mcmoy_RA.At(i);
	mcmozTemp = mcmoz_RA.At(i);

	mPtTemp = sqrt((mcmoxTemp*mcmoxTemp) + (mcmoyTemp*mcmoyTemp));
	mPATemp = atan2(mPtTemp, mcmozTemp);
	mATemp = atan2(mcmoxTemp, mcmoyTemp);
	
	realAll_pt->Fill(mPtTemp);
	realAll_PA->Fill(mPATemp);
	realAll_azimuth->Fill(mATemp);
	break; //this might cause issues in the future, double check to make sure that there are not two or more generating particles
      }
    }

    for(int i = 0; i < rcmox_RA.GetSize(); i++) {
      rcmoxTemp = rcmox_RA.At(i);
      rcmoyTemp = rcmoy_RA.At(i);
      rcmozTemp = rcmoz_RA.At(i);

      rPtTemp = sqrt((rcmoxTemp*rcmoxTemp) + (rcmoyTemp*rcmoyTemp));
      rPATemp = atan2(rPtTemp, rcmozTemp);
      rATemp = atan2(rcmoxTemp, rcmoyTemp);
      
      if(abs(rctyp_RA.At(i)) == 11) {
	electronReco_pt->Fill(rPtTemp);
	electronReco_PA->Fill(rPATemp);
	electronReco_azimuth->Fill(rATemp);
      } else if(abs(rctyp_RA.At(i)) == 22) {
	photonReco_pt->Fill(rPtTemp);
	photonReco_PA->Fill(rPATemp);
	photonReco_azimuth->Fill(rATemp);
      } else if(abs(rctyp_RA.At(i)) == 2112) {
	neutronReco_pt->Fill(rPtTemp);
	neutronReco_PA->Fill(rPATemp);
	neutronReco_azimuth->Fill(rATemp);
      }
    }
    count++;
    r2wMax = 0;
    r2wMax_Index = -1;
  }

  gStyle->SetPalette(kRust);
  
  //////////////////
  //COMBINED HISTS//
  //////////////////

  //draws the pt histograms for electrons, photons, and neutrons on top of each other
  THStack *hs = new THStack("hs", "Transverse momenta of electrons, photons, and neutrons;Transverse momentum (GeV);Count");

  c->SetLogy();
  
  electronReco_pt->SetLineColor(kBlue);
  photonReco_pt->SetLineColor(kRed);
  neutronReco_pt->SetLineColor(kBlack);

  hs->Add(electronReco_pt);
  hs->Add(photonReco_pt);
  hs->Add(neutronReco_pt);
  
  hs->Draw("nostack");

  gPad->SetGrid(1, 0);
  gPad->BuildLegend(0.75, 0.75, 0.9, 0.9, "");
  
  c->SaveAs("filterOutput/allPart_pt.png");
  c->Close();
  c = new TCanvas();

  //draws the polar angle histograms for electrons, photons, and neutrons on top of each other
  hs = new THStack("hs", "Polar angle of electrons, photons, and neutrons;Polar angle (Rads);Count");

  c->SetLogy();

  electronReco_PA->SetLineColor(kBlue);
  photonReco_PA->SetLineColor(kRed);
  neutronReco_PA->SetLineColor(kBlack);

  hs->Add(electronReco_PA);
  hs->Add(photonReco_PA);
  hs->Add(neutronReco_PA);
    
  hs->Draw("nostack");

  gPad->SetGrid(1, 0);
  gPad->BuildLegend(0.7, 0.1, 0.9, 0.3, "");
  
  c->SaveAs("filterOutput/allPart_PA.png");
  c->Close();
  c = new TCanvas();

  //draws the azimuth histograms for electrons, photons, and neutrons on top of each other
  hs = new THStack("hs", "Azimuth of electrons, photons, and neutrons;Azimuth (Rads);Count");

  c->SetLogy();

  electronReco_azimuth->SetLineColor(kBlue);
  photonReco_azimuth->SetLineColor(kRed);
  neutronReco_azimuth->SetLineColor(kBlack);

  hs->Add(electronReco_azimuth);
  hs->Add(photonReco_azimuth);
  hs->Add(neutronReco_azimuth);

  hs->Draw("nostack");

  gPad->SetGrid(1, 0);
  gPad->BuildLegend(0.4, 0.1, 0.6, 0.3, "");

  c->SaveAs("filterOutput/allPart_azimuth.png");
  c->Close();
  c = new TCanvas();
  
  //draws the two pt histograms on each other
  hs = new THStack("hs", "pt of linked generating electrons and all generating electrons (no BIB);Transverse momentum (GeV);Count");

  realAll_pt->SetLineColor(kBlue);
  realAll_pt->SetLineColor(kRed);

  hs->Add(realAll_pt);
  hs->Add(realPassed_pt);

  hs->Draw("nostack");

  gPad->SetGrid(1, 0);
  gPad->BuildLegend(0.1, 0.1, 0.45, 0.3, "");
  
  c->SaveAs("filterOutput/compLink_pt.png");
  c->Close();
  c = new TCanvas();

  //rinse and repeat for the polar angle instead of the transverse momentum
  hs = new THStack("hs", "Polar angle of linked generating electrons and all generating electrons (no BIB);Polar angle (Rads);Count");
  
  realAll_PA->SetLineColor(kBlue);
  realAll_PA->SetLineColor(kRed);

  hs->Add(realAll_PA);
  hs->Add(realPassed_PA);

  hs->Draw("nostack");

  gPad->SetGrid(1, 0);
  gPad->BuildLegend(0.1, 0.1, 0.3, 0.25, "");

  c->SaveAs("filterOutput/compLink_PA.png");
  c->Close();
  c = new TCanvas();

  //rinse and repeat for the azimuth instead of the transverse momentum
  hs = new THStack("hs", "Azimuth of linked generating electrons and all generating electrons (no BIB);Azimuth (Rads);Count");

  realAll_azimuth->SetLineColor(kBlue);
  realAll_azimuth->SetLineColor(kRed);

  hs->Add(realAll_azimuth);
  hs->Add(realPassed_azimuth);

  hs->Draw("nostack");

  gPad->SetGrid(1, 0);
  gPad->BuildLegend(0.1, 0.1, 0.3, 0.25, "");

  c->SaveAs("filterOutput/compLink_azimuth.png");
  c->Close();
  c = new TCanvas();

  ///////////////
  //SOLO GRAPHS//
  ///////////////
  
  //draws the histogram that shows the distrobution of the relation weights
  realWeights_Hist->Draw();
  realWeights_Hist->SetTitle("Weights of linked particles");
  realWeights_Hist->GetXaxis()->SetTitle("Weight");
  realWeights_Hist->GetYaxis()->SetTitle("Count");
  
  c->SaveAs("filterOutput/recoLinkWeights.png");
  c->Close();
  c = new TCanvas();

  //draws a hist that hold the differences in electron pt
  diffPt_Hist->Draw();
  diffPt_Hist->SetTitle("Difference in count between all and reconstructed particles");
  diffPt_Hist->GetXaxis()->SetTitle("Transverse momentum (GeV)");
  diffPt_Hist->GetYaxis()->SetTitle("Count");

  c->SaveAs("filterOutput/recoPt_diff.png");
  c->Close();
  c = new TCanvas();

  //draws a histogram for delta r (distance between polar angle and lambda on a plot between the real and reco data)
  c->SetLogy();
    
  deltaR_Hist->Draw();
  deltaR_Hist->SetTitle("Delta r between the real and reco data points");
  deltaR_Hist->GetXaxis()->SetTitle("delta R (Rads)");
  deltaR_Hist->GetYaxis()->SetTitle("Count");

  c->SaveAs("filterOutput/deltaR.png");
  c->Close();
  c = new TCanvas();

  ////////////////
  //EFFICIENCIES//
  ////////////////
  
  //draws an efficiency graph for transverse momentum
  TEfficiency* pt_Eff = new TEfficiency("pt_Eff", "Efficiency of transverse momentum reconstruction;Transverse momentum (GeV);Efficiency", 20, 0, 2000);
  TFile* eff_File = new TFile("effFile.root", "recreate");
  
  if(TEfficiency::CheckConsistency(*realPassed_pt, *realAll_pt)) {
    pt_Eff->SetPassedHistogram(*realPassed_pt, "f");
    pt_Eff->SetTotalHistogram(*realAll_pt, "f");
    eff_File->Write();
  }

  c->Clear();
  
  pt_Eff->Draw();

  c->SaveAs("filterOutput/Eff_pt.png");
  c->Close();
  c = new TCanvas();

  //draws an efficiency graph for the polar angle
  TEfficiency* PA_Eff = new TEfficiency("PA_Eff", "Efficiency of polar angle reconstruction;Polar angle (Rads);Efficiency", 20, 0, 1.6);
  eff_File = new TFile("effFile.root", "recreate");

  if(TEfficiency::CheckConsistency(*realPassed_PA, *realAll_PA)) {
    PA_Eff->SetPassedHistogram(*realPassed_PA, "f");
    PA_Eff->SetTotalHistogram(*realAll_PA, "f");
    eff_File->Write();
  }

  c->Clear();

  PA_Eff->Draw();

  c->SaveAs("filterOutput/Eff_PA.png");
  c->Close();
  c = new TCanvas();

  //draws an efficiency graph for the azimuth
  TEfficiency* azimuth_Eff = new TEfficiency("azimuth_Eff", "Efficiency of azimuth reconstruction;Azimuth (Rads);Efficiency", 20, -1.6, 1.6);
  eff_File = new TFile("effFile.root", "recreate");

  if(TEfficiency::CheckConsistency(*realPassed_azimuth, *realAll_azimuth)) {
    azimuth_Eff->SetPassedHistogram(*realPassed_azimuth, "f");
    azimuth_Eff->SetTotalHistogram(*realAll_azimuth, "f");
    eff_File->Write();
  }

  c->Clear();

  azimuth_Eff->Draw();

  c->SaveAs("filterOutput/Eff_azimuth.png");

  c->Clear();
  c->Close();
}
