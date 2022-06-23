void filter()
{ 
  auto fileName = "ntuple.root";
  auto treeName = "MyLCTuple";

  TFile *myFile = new TFile("ntuple.root");
  TTree *myTree = (TTree*)myFile->Get("MyLCTuple");

  TCanvas *c = new TCanvas();
  
  //opens the file to be read
  auto openFile = TFile::Open("ntuple.root");

  /////////////////////////
  //Histogram definitions//
  /////////////////////////

  //histograms for attributes of electrons
  //transverse momentum
  TH1F *realPassed_pt = new TH1F("rP_pt", "Linked electrons", 20, 0, 2000); //transverse momentum of electrons that are linked to reconstructed electrons
  TH1F *realAllPassed_pt = new TH1F("AP_pt", "Linked particles", 20, 0, 2000); //transverse momentum All truth electrons
  TH1F *realAll_pt = new TH1F("rA_pt", "All electrons", 20, 0, 2000); //Links the most likely particles together, regardless of reconstructed type
  //Polar angle
  TH1F *realPassed_PA = new TH1F("rP_PA", "Linked electrons", 20, 0, 3.2); //see above
  TH1F *realAllPassed_PA = new TH1F("AP_PA", "Linked particles", 20, 0, 3.2);
  TH1F *realAll_PA = new TH1F("rA_PA", "All electrons", 20, 0, 3.2);
  //Azimuth
  TH1F *realPassed_azimuth = new TH1F("rP_a", "Linked electrons", 20, -3.2, 3.2); //see above
  TH1F *realAllPassed_azimuth = new TH1F("AP_a", "Linked particles", 20, -3.2, 3.2);
  TH1F *realAll_azimuth = new TH1F("rA_a", "All electrons", 20, -3.2, 3.2);

  //histograms for weights of particles
  TH1D *realElWeights_Hist = new TH1D("realElWeights", "Reco link electron weights", 26, 0, 1.3); //link weights of reconstructed electrons linked to truth electrons
  TH1D *realAllWeights_Hist = new TH1D("realAllWeights", "Reco link all weights", 26, 0, 1.3); // link weights of particles linked to truth electrons

  //histograms for storing differences between values
  TH1F *diffPt_Hist = new TH1F("diffPt", "MyLCTuple", 40, -500, 1500); //difference between reconstructed linked particles and the truth particles in transverse momentum
  TH1F *deltaR_Hist = new TH1F("deltaR", "MyLCTuple", 25, 0, 0.0025); //delta R (delta polar angle and azimuth combined)
  TH1F *accEPSum_pt = new TH1F("aEP_pt", "Summed electrons and protons", 20, 0, 2000); //accuracy of transverse momentum reconstruction when the reconstructed particle is made up of all reconstructed electrons and photons
  TH1F *accEPSum_PA = new TH1F("aEP_PA", "Summed electrons and protons", 20, 0, 3.2); //accuracy of polar angle reconstructed when the reconstruction particle is made up of all reconstructed electrons and photons
  TH1F *accEPSum_azimuth = new TH1F("aEP_a", "Summed electrons and protons", 20, -3.2, 3.2); //accuracy of azimuth reconstruction when the reconstructed particle is made up of all reconstructed electrons and photons
  TH1F *accBestMatch_pt = new TH1F("aBM_pt", "Transverse momentum accuracy of best match link", 20, 0, 2000); //accuracy of transverse momentum reconstruction when the reconstructed particle is the best match
  TH1F *accBestMatch_PA = new TH1F("aMB_PA", "Polar angle accuracy of best match link", 20, 0, 3.2); //accuracy of polar angle reconstructed when the reconstruction particle is the best match
  TH1F *accBestMatch_azimuth = new TH1F("aMB_a", "Azimuth accuracy of best match link", 20, -3.2, 3.2); //accuracy of azimuth reconstruction when the reconstructed particle is made up of the best match

  //histograms for storing attributes unique to certain reconstructed types of particles
  //transverse momentum
  TH1F *electronReco_pt = new TH1F("elRec_pt", "Electron pt", 20, 0, 2000);
  TH1F *photonReco_pt = new TH1F("phRec_pt", "Photon pt", 20, 0, 2000);
  TH1F *neutronReco_pt = new TH1F("nuRec_pt", "Neutron pt", 20, 0, 2000);

  //polar angle
  TH1F *electronReco_PA = new TH1F("elRec_PA", "Electron PA", 20, 0, 3.2);
  TH1F *photonReco_PA = new TH1F("phRec_PA", "Photon PA", 20, 0, 3.2);
  TH1F *neutronReco_PA = new TH1F("nuRec_PA", "Neutron PA", 20, 0, 3.2);

  //azimuth
  TH1F *electronReco_azimuth = new TH1F("elRec_a", "Electron azimuth", 20, -3.2, 3.2);
  TH1F *photonReco_azimuth = new TH1F("phRec_a", "Photon azimuth", 20, -3.2, 3.2);
  TH1F *neutronReco_azimuth = new TH1F("nuRec_a", "Neutron azimuth", 20, -3.2, 3.2);

  //creates a reader that will traverse the events of the simulation
  TTreeReader myReader("MyLCTuple", openFile);

  //arrays that hold the values of branches
  TTreeReaderArray<int> r2f_RA(myReader, "r2f"); //link reco particle index number
  TTreeReaderArray<int> r2t_RA(myReader, "r2t"); //link truth partcie index number
  TTreeReaderArray<Float_t> r2w_RA(myReader, "r2w"); //link weight
  TTreeReaderArray<Float_t> mcmox_RA(myReader, "mcmox"); //truth x momentum
  TTreeReaderArray<Float_t> mcmoy_RA(myReader, "mcmoy"); //truth y momentum
  TTreeReaderArray<Float_t> mcmoz_RA(myReader, "mcmoz"); //truth z momentum
  TTreeReaderArray<Float_t> rcmox_RA(myReader, "rcmox"); //reco x momentum
  TTreeReaderArray<Float_t> rcmoy_RA(myReader, "rcmoy"); //reco y momentum
  TTreeReaderArray<Float_t> rcmoz_RA(myReader, "rcmoz"); //reco z momentum
  TTreeReaderArray<int> rctyp_RA(myReader, "rctyp"); //reconstructed type (11 = electron, 22 = photon, 2112 = neutron)
  TTreeReaderArray<int> mcpdg_RA(myReader, "mcpdg"); //truth type (see list above)
  TTreeReaderArray<int> mcgst_RA(myReader, "mcgst"); //is this particle a generating particle or not

  //temporary variables that will reduce the number of get operatons
  int r2fTemp, r2tTemp;
  Float_t mcmoxTemp, mcmoyTemp, mcmozTemp, rcmoxTemp, rcmoyTemp, rcmozTemp, r2wTemp, mPtTemp, mPATemp, mATemp, rPtTemp, rPATemp, rATemp, deltaRTemp, dPATemp, dATemp, rcmoxSumTemp, rcmoySumTemp, rcmozSumTemp, BMptTemp, BMPATemp, BMATemp;
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
	if(r2wTemp > r2wMax) {
	  r2wMax = r2wTemp;
	  r2wMax_Index = i;
	}
	//adds the weight to a histogram
	realElWeights_Hist->Fill(r2wTemp);
      }
    }

    //if a relationship was found between those particles, it is added to histograms
    if(r2wMax_Index != -1) {

      //find the particles that the particle linker refers to
      r2tTemp = r2t_RA.At(r2wMax_Index);
      r2fTemp = r2f_RA.At(r2wMax_Index);
      r2wTemp = r2w_RA.At(r2wMax_Index);

      //finds the momenta of the reco and truth particles
      mcmoxTemp = mcmox_RA.At(r2tTemp);
      mcmoyTemp = mcmoy_RA.At(r2tTemp);
      mcmozTemp = mcmoz_RA.At(r2tTemp);
      rcmoxTemp = rcmox_RA.At(r2fTemp);
      rcmoyTemp = rcmoy_RA.At(r2fTemp);
      rcmozTemp = rcmoz_RA.At(r2fTemp);

      //calculates the attributes (transverse momentum, polar angle, and azimuth)
      mPtTemp = sqrt((mcmoxTemp*mcmoxTemp) + (mcmoyTemp*mcmoyTemp));
      mPATemp = atan2(mPtTemp, mcmozTemp);
      mATemp = atan2(mcmoxTemp, mcmoyTemp);
      rPtTemp = sqrt((rcmoxTemp*rcmoxTemp) + (rcmoyTemp*rcmoyTemp));
      rPATemp = atan2(rPtTemp, rcmozTemp);
      rATemp = atan2(rcmoxTemp, rcmoyTemp);

      //finds the differences between truth and reconstructed attributes
      dPATemp = mPATemp - rPATemp;
      dATemp = mATemp - rATemp;

      //fills in histograms with attributes
      realPassed_pt->Fill(mPtTemp);
      realPassed_PA->Fill(mPATemp);
      realPassed_azimuth->Fill(mATemp);

      //fills in histograms with differences
      diffPt_Hist->Fill(mPtTemp - rPtTemp);
      deltaR_Hist->Fill(sqrt((dPATemp*dPATemp) + (dATemp*dATemp)));
      }

    //if no link was found between a reconstructed and truth electron, I expand the search to reconstructed photons as well
    if(r2wMax_Index == -1) {
      for(int i = 0; i < r2f_RA.GetSize(); i++) {
	//find the two particles in which the link occurs between
	r2tTemp = r2t_RA.At(i);
	r2fTemp = r2f_RA.At(i);
	r2wTemp = r2w_RA.At(i);

	//checks to see if the truth particle is an electron, if the reconstructed particle is a photon, and if the truth particle is an originial generating particle
	if(abs(mcpdg_RA.At(r2tTemp)) == 11 && abs(rctyp_RA.At(r2fTemp)) == 22 && mcgst_RA.At(r2tTemp)) {
	  if(r2wTemp > r2wMax) {
	    r2wMax = r2wTemp;
	    r2wMax_Index = i;
	  }
	  //adds the weight to a histogram
	  realAllWeights_Hist->Fill(r2wTemp);
	}
      }
    }

    //repeats the previous steps that adds data to histograms but instead adds it to the histograms that holds links between photons and electrons as well
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

      BMptTemp = sqrt((mcmoxTemp*mcmoxTemp) + (mcmoyTemp*mcmoyTemp));
      BMPATemp = atan2(mPtTemp, mcmozTemp);
      BMATemp = atan2(mcmoxTemp, mcmoyTemp);

      realAllPassed_pt->Fill(BMptTemp);
      realAllPassed_PA->Fill(BMPATemp);
      realAllPassed_azimuth->Fill(BMATemp);
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

	//adds generating particle to histogram
	realAll_pt->Fill(mPtTemp);
	realAll_PA->Fill(mPATemp);
	realAll_azimuth->Fill(mATemp);
	break; //this might cause issues in the future, double check to make sure that there are not two or more generating particles
      }
    }

    //resets data sum variables
    rcmoxSumTemp = 0;
    rcmoySumTemp = 0;
    rcmozSumTemp = 0;

    //loops over all reconstructed particles
    for(int i = 0; i < rcmox_RA.GetSize(); i++) {
      rcmoxTemp = rcmox_RA.At(i);
      rcmoyTemp = rcmoy_RA.At(i);
      rcmozTemp = rcmoz_RA.At(i);

      //if that particle is an electron or photon, add it the sum variables (this accounts for particle mis-identification)
      if(abs(rctyp_RA.At(i)) == 11 || abs(rctyp_RA.At(i)) == 22) {
	rcmoxSumTemp += rcmoxTemp;
	rcmoySumTemp += rcmoyTemp;
	rcmozSumTemp += rcmozTemp;
      }

      //calculates the 
      rPtTemp = sqrt((rcmoxTemp*rcmoxTemp) + (rcmoyTemp*rcmoyTemp));
      rPATemp = atan2(rPtTemp, rcmozTemp);
      rATemp = atan2(rcmoxTemp, rcmoyTemp);

      //if the particle is a electron, add it to one set of histograms
      if(abs(rctyp_RA.At(i)) == 11) {
	electronReco_pt->Fill(rPtTemp);
	electronReco_PA->Fill(rPATemp);
	electronReco_azimuth->Fill(rATemp);

	//if the particle is a photon, add it to another set of histograms
      } else if(abs(rctyp_RA.At(i)) == 22) {
	photonReco_pt->Fill(rPtTemp);
	photonReco_PA->Fill(rPATemp);
	photonReco_azimuth->Fill(rATemp);

	//if the particle is a neutron, add it to a third set of histograms
      } else if(abs(rctyp_RA.At(i)) == 2112) {
	neutronReco_pt->Fill(rPtTemp);
	neutronReco_PA->Fill(rPATemp);
	neutronReco_azimuth->Fill(rATemp);
      }
    }
    count++;
    r2wMax = 0;
    r2wMax_Index = -1;

    //finds the difference between the truth and reconstructed particles (sum of photons and electrons)
    accEPSum_pt->Fill(mPtTemp - sqrt((rcmoxSumTemp*rcmoxSumTemp) + (rcmoySumTemp*rcmoySumTemp)));
    accEPSum_PA->Fill(mPATemp - atan2(sqrt((rcmoxSumTemp*rcmoxSumTemp) + (rcmoySumTemp*rcmoySumTemp)), rcmozSumTemp));
    accEPSum_azimuth->Fill(mATemp - atan2(rcmoxSumTemp, rcmoySumTemp));

    //finds the difference between the truth and best match reconstructed particles
    accBestMatch_pt->Fill(mPtTemp - BMptTemp);
    accBestMatch_PA->Fill(mPATemp - BMPATemp);
    accBestMatch_azimuth->Fill(mATemp - BMATemp);
  }

  gStyle->SetPalette(kRust);
  
  //////////////////
  //COMBINED HISTS//
  //////////////////
  
  //=====================================//
  //all types of reco particles overlayed//
  //=====================================//

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

  //==============================//
  //most likely electron reco link//
  //==============================//
  
  //draws the two pt histograms on each other
  hs = new THStack("hs", "pt of linked generating electrons and all generating electrons (no BIB);Transverse momentum (GeV);Count");

  realAll_pt->SetLineColor(kBlue);
  realPassed_pt->SetLineColor(kRed);

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
  realPassed_PA->SetLineColor(kRed);

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
  realPassed_azimuth->SetLineColor(kRed);

  hs->Add(realAll_azimuth);
  hs->Add(realPassed_azimuth);

  hs->Draw("nostack");

  gPad->SetGrid(1, 0);
  gPad->BuildLegend(0.1, 0.1, 0.3, 0.25, "");

  c->SaveAs("filterOutput/compLink_azimuth.png");
  c->Close();
  c = new TCanvas();
  
  //====================================================//
  //most likely linked particles regardless of reco type//
  //====================================================//

  //draws the two pt histograms on each other
  hs = new THStack("hs", "pt of linked particles and all generating electrons (no BIB);Transverse momentum (GeV);Count");

  realAll_pt->SetLineColor(kBlue);
  realAllPassed_pt->SetLineColor(kRed);

  hs->Add(realAll_pt);
  hs->Add(realAllPassed_pt);

  hs->Draw("nostack");

  gPad->SetGrid(1, 0);
  gPad->BuildLegend(0.1, 0.1, 0.45, 0.3, "");

  c->SaveAs("filterOutput/compAllLink_pt.png");
  c->Close();
  c = new TCanvas();

  //rinse and repeat for the polar angle instead of the transverse momentum
  hs = new THStack("hs", "Polar angle of linked particles and all generating electrons (no BIB);Polar angle (Rads);Count");

  realAll_PA->SetLineColor(kBlue);
  realAllPassed_PA->SetLineColor(kRed);

  hs->Add(realAll_PA);
  hs->Add(realAllPassed_PA);

  hs->Draw("nostack");

  gPad->SetGrid(1, 0);
  gPad->BuildLegend(0.1, 0.1, 0.3, 0.25, "");

  c->SaveAs("filterOutput/compAllLink_PA.png");
  c->Close();
  c = new TCanvas();

  //rinse and repeat for the azimuth instead of the transverse momentum
  hs = new THStack("hs", "Azimuth of linked particles and all generating electrons (no BIB);Azimuth (Rads);Count");

  realAll_azimuth->SetLineColor(kBlue);
  realAllPassed_azimuth->SetLineColor(kRed);

  hs->Add(realAll_azimuth);
  hs->Add(realAllPassed_azimuth);

  hs->Draw("nostack");

  gPad->SetGrid(1, 0);
  gPad->BuildLegend(0.1, 0.1, 0.3, 0.25, "");

  c->SaveAs("filterOutput/compAllLink_azimuth.png");
  c->Close();
  c = new TCanvas();

  
  ///////////////
  //SOLO GRAPHS//
  ///////////////
  
  //draws the histogram that shows the distrobution of the linked electron weights
  realElWeights_Hist->Draw();
  realElWeights_Hist->SetTitle("Weights of linked electrons");
  realElWeights_Hist->GetXaxis()->SetTitle("Weight");
  realElWeights_Hist->GetYaxis()->SetTitle("Count");
  
  c->SaveAs("filterOutput/recoLinkElWeights.png");
  c->Close();
  c = new TCanvas();

  //draws the histogram that shows the distrobution of the relation weights for all particles
  realAllWeights_Hist->Draw();
  realAllWeights_Hist->SetTitle("Weights of linked particles");
  realAllWeights_Hist->GetXaxis()->SetTitle("Weight");
  realAllWeights_Hist->GetYaxis()->SetTitle("Count");

  c->SaveAs("filterOutput/recoLinkAllWeights.png");
  c->Close();
  c = new TCanvas();

  //draws a hist that hold the differences in electron pt
  diffPt_Hist->Draw();
  diffPt_Hist->SetTitle("Difference in pt between all and reconstructed particles");
  diffPt_Hist->GetXaxis()->SetTitle("Transverse momentum (GeV)");
  diffPt_Hist->GetYaxis()->SetTitle("Count");

  c->SaveAs("filterOutput/recoPt_diff.png");
  c->Close();
  c = new TCanvas();

  //=============================================//
  //summed data from linked electrons and photons//
  //=============================================//

  //draws a hist that hold the differences in electron and photon sum pt
  accEPSum_pt->Draw();
  accEPSum_pt->SetTitle("Difference in pt between all and sum of electrons and photons");
  accEPSum_pt->GetXaxis()->SetTitle("Transverse momentum (GeV)");
  accEPSum_pt->GetYaxis()->SetTitle("Count");

  c->SaveAs("filterOutput/accEPSumLink_pt.png");
  c->Close();
  c = new TCanvas();
  
  //draws a hist that hold the differences in electron and photon sum polar angle
  accEPSum_PA->Draw();
  accEPSum_PA->SetTitle("Difference in polar angle between all and sum of electrons and photons");
  accEPSum_PA->GetXaxis()->SetTitle("Polar angle (Rads)");
  accEPSum_PA->GetYaxis()->SetTitle("Count");

  c->SaveAs("filterOutput/accEPSumLink_PA.png");
  c->Close();
  c = new TCanvas();
  
  //draws a hist that hold the differences in electron and photon sum azimuth
  accEPSum_azimuth->Draw();
  accEPSum_azimuth->SetTitle("Difference in azimuth between all and sum of electrons and photons");
  accEPSum_azimuth->GetXaxis()->SetTitle("Azimmuth (Rads)");
  accEPSum_azimuth->GetYaxis()->SetTitle("Count");

  c->SaveAs("filterOutput/accEPSumLink_azimuth.png");
  c->Close();
  c = new TCanvas();

  ////////////////////////////
  //Best match accuracy data//
  ////////////////////////////

  //draws a hist that hold the differences in pt in between the truth particle and the best match link regardless of reco type
  accBestMatch_pt->Draw();
  accBestMatch_pt->SetTitle("Difference in pt between truth and best match link");
  accBestMatch_pt->GetXaxis()->SetTitle("Transverse momentum (GeV)");
  accBestMatch_pt->GetYaxis()->SetTitle("Count");

  c->SaveAs("filterOutput/accBestMatchLink_pt.png");
  c->Close();
  c = new TCanvas();

  //draws a hist that hold the differences in polar angle in between the truth particle and the best match link regardless of reco type
  accBestMatch_PA->Draw();
  accBestMatch_PA->SetTitle("Difference in polar angle between truth and best match link");
  accBestMatch_PA->GetXaxis()->SetTitle("Polar angle (Rads)");
  accBestMatch_PA->GetYaxis()->SetTitle("Count");

  c->SaveAs("filterOutput/accBestMatchLink_PA.png");
  c->Close();
  c = new TCanvas();

  //draws a hist that hold the differences in azimuth in between the truth particle and the best match link regardless of reco type
  accBestMatch_azimuth->Draw();
  accBestMatch_azimuth->SetTitle("Difference in azimuth between truth and best match link");
  accBestMatch_azimuth->GetXaxis()->SetTitle("Azimmuth (Rads)");
  accBestMatch_azimuth->GetYaxis()->SetTitle("Count");

  c->SaveAs("filterOutput/accBestMatchLink_azimuth.png");
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

  //checks to make sure if the two histograms are compatable to make an efficiency histogram
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
  
  //checks to make sure if the two histograms are compatable to make an efficiency histogram
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

  //checks to make sure if the two histograms are compatable to make an efficiency histogram
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

  //draws an efficiency graph for transverse momentum
  TEfficiency* ptAll_Eff = new TEfficiency("pt_Eff", "Efficiency of transverse momentum reconstruction for all particles;Transverse momentum (GeV);Efficiency", 20, 0, 2000);
  eff_File = new TFile("effFile.root", "recreate");
  
  //checks to make sure if the two histograms are compatable to make an efficiency histogram
  if(TEfficiency::CheckConsistency(*realAllPassed_pt, *realAll_pt)) {
    ptAll_Eff->SetPassedHistogram(*realAllPassed_pt, "f");
    ptAll_Eff->SetTotalHistogram(*realAll_pt, "f");
    eff_File->Write();
  }

  c->Clear();

  ptAll_Eff->Draw();

  c->SaveAs("filterOutput/EffAll_pt.png");
  c->Close();
  c = new TCanvas();

  //draws an efficiency graph for the polar angle
  TEfficiency* PAAll_Eff = new TEfficiency("PA_Eff", "Efficiency of polar angle reconstruction for all particles;Polar angle (Rads);Efficiency", 20, 0, 1.6);
  eff_File = new TFile("effFile.root", "recreate");

  //checks to make sure if the two histograms are compatable to make an efficiency histogram
  if(TEfficiency::CheckConsistency(*realAllPassed_PA, *realAll_PA)) {
    PAAll_Eff->SetPassedHistogram(*realAllPassed_PA, "f");
    PAAll_Eff->SetTotalHistogram(*realAll_PA, "f");
    eff_File->Write();
  }

  c->Clear();

  PAAll_Eff->Draw();

  c->SaveAs("filterOutput/EffAll_PA.png");
  c->Close();
  c = new TCanvas();

  //draws an efficiency graph for the azimuth for 
  TEfficiency* azimuthAll_Eff = new TEfficiency("azimuth_Eff", "Efficiency of azimuth reconstruction for all particles;Azimuth (Rads);Eff\
iciency", 20, -1.6, 1.6);
  eff_File = new TFile("effFile.root", "recreate");

  //checks to make sure if the two histograms are compatable to make an efficiency histogram
  if(TEfficiency::CheckConsistency(*realAllPassed_azimuth, *realAll_azimuth)) {
    azimuthAll_Eff->SetPassedHistogram(*realAllPassed_azimuth, "f");
    azimuthAll_Eff->SetTotalHistogram(*realAll_azimuth, "f");
    eff_File->Write();
  }

  c->Clear();

  azimuthAll_Eff->Draw();

  c->SaveAs("filterOutput/EffAll_azimuth.png");
  c->Close();
}
