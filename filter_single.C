void filter_single()
{
  /////////////////////////////
  //GENERAL ANALYSIS SETTINGS//
  /////////////////////////////
  
  auto fileName = "10k-ACTS-redo.root";
  auto treeName = "MyLCTuple";

  TString saveDir = "10k-ACTS-redo";

  /////////////////////////////

  //TFile *myFile = new TFile("ntuple.root");
  TFile *myFile = new TFile(fileName);
  TTree *myTree = (TTree*)myFile->Get(treeName);

  TCanvas *c = new TCanvas();
  
  //opens the file to be read
  auto openFile = TFile::Open(fileName);

  /////////////////////////
  //Histogram definitions//
  /////////////////////////

  //histograms for attributes of electrons
  //transverse momentum
  TH1F *realPassed_pt = new TH1F("rP_pt", "Linked electrons", 20, 0, 2000); //transverse momentum of electrons that are linked to reconstructed electrons
  TH1F *realAllPassed_pt = new TH1F("AP_pt", "Best match", 20, 0, 2000); //best match transverse momentum from linked electrons and photons
  TH1F *realEPSum_pt = new TH1F("EPS_pt", "Summed particles", 20, 0, 2000); //summed transverse momentum from linked electrons and photons
  TH1F *recoTrackLink_pt = new TH1F("recTck_pt", "Reco track", 40, -2000, 2000); //transverse momentum from individual track reconstruction segments
  TH1F *realAll_pt = new TH1F("rA_pt", "All electrons", 20, 0, 2000); //transverse momentum of the truth particles
  TH1F *realPion_pt = new TH1F("rPi_pt", "Only pions", 20, 0, 2000); //transverse momentum of linked pions
  //Polar angle
  TH1F *realPassed_PA = new TH1F("rP_PA", "Linked electrons", 20, 0, 3.2); //see above
  TH1F *realAllPassed_PA = new TH1F("AP_PA", "Best match", 20, 0, 3.2);
  TH1F *realEPSum_PA = new TH1F("EPS_PA", "Summed particles", 20, 0, 3.2);
  TH1F *realAll_PA = new TH1F("rA_PA", "All electrons", 20, 0, 3.2);
  TH1F *realPion_PA = new TH1F("rPi_PA", "Only pions", 20, 0, 3.2);
  //Azimuth
  TH1F *realPassed_azimuth = new TH1F("rP_a", "Linked electrons", 20, -3.2, 3.2); //see above
  TH1F *realAllPassed_azimuth = new TH1F("AP_a", "Best match", 20, -3.2, 3.2);
  TH1F *realEPSum_azimuth = new TH1F("EPS_a", "Summed particles", 20, -3.2, 3.2);
  TH1F *realAll_azimuth = new TH1F("rA_a", "All electrons", 20, -3.2, 3.2);
  TH1F *realPion_azimuth = new TH1F("rPi_a", "Only pions", 20, -3.2, 3.2);

  //histograms for weights of particles
  TH1D *realElWeights_Hist = new TH1D("realElWeights", "Reco link electron weights", 26, 0, 1.3); //link weights of reconstructed electrons linked to truth electrons
  TH1D *realAllWeights_Hist = new TH1D("realAllWeights", "Reco link all weights", 26, 0, 1.3); // link weights of particles linked to truth electrons

  //histograms for storing differences between values
  TH1F *diffPt_Hist = new TH1F("diffPt", "Simple electron linking", 40, -2000, 2000); //difference between reconstructed linked particles and the truth particles in transverse momentum
  TH1F *diffPA_Hist = new TH1F("diffPA", "Simple electron linking", 20, -3.2, 3.2); //difference between reconstructed linked particles and the truth particles in transverse momentum
  TH1F *diffA_Hist = new TH1F("diffA", "Simple electron linking", 20, -3.2, 3.2); //difference between reconstructed linked particles and the truth particles in transverse momentum
  TH1F *deltaR_Hist = new TH1F("deltaR", "MyLCTuple", 25, 0, 0.0025); //delta R (delta polar angle and azimuth combined)
  TH1F *resEPSum_pt = new TH1F("rEP_pt", "Summed electrons and protons", 40, -2000, 2000); //resolution of transverse momentum reconstruction when the reconstructed particle is made up of all reconstructed electrons and photons
  TH1F *resEPSum_PA = new TH1F("rEP_PA", "Summed electrons and protons", 20, -3.2, 3.2); //resolution of polar angle reconstructed when the reconstruction particle is made up of all reconstructed electrons and photons
  TH1F *resEPSum_azimuth = new TH1F("rEP_a", "Summed electrons and protons", 20, -3.2, 3.2); //resolution of azimuth reconstruction when the reconstructed particle is made up of all reconstructed electrons and photons
  TH1F *resBestMatch_pt = new TH1F("rBM_pt", "Best match link", 40, -2000, 2000); //resolution of transverse momentum reconstruction when the reconstructed particle is the best match
  TH1F *resBestMatch_PA = new TH1F("rMB_PA", "Best match link", 20, -3.2, 3.2); //resolution of polar angle reconstructed when the reconstruction particle is the best match
  TH1F *resBestMatch_azimuth = new TH1F("rMB_a", "Best match link", 20, -3.2, 3.2); //resolution of azimuth reconstruction when the reconstructed particle is made up of the best match
  TH1F *resPhotonReco_pt = new TH1F("rPR_pt", "Best match photon link", 40, -2000, 2000); //resolution of transverse momentum reconstruction when the reconstructed particle is the best photon match
  TH1F *resPhotonReco_PA = new TH1F("rPR_PA", "Best match photon link", 20, -3.2, 3.2); //resolution of transverse momentum reconstruction when the reconstructed particle is the best photon match
  TH1F *resPhotonReco_azimuth = new TH1F("rPR_A", "Best match photon link", 20, -3.2, 3.2); //resolution of transverse momentum reconstruction when the reconstructed particle is the best photon match
  TH1F *resPionReco_pt = new TH1F("rPiR_pt", "Best match pion link", 40, -2000, 2000); //resolution of transverse momentum reconstruction when the reconstructed particle is the best photon match
  TH1F *resPionReco_PA = new TH1F("rPiR_PA", "Best match pion link", 20, -3.2, 3.2); //resolution of transverse momentum reconstruction when the reconstructed particle is the best photon match
  TH1F *resPionReco_azimuth = new TH1F("rPiR_A", "Best match pion link", 20, -3.2, 3.2); //resolution of transverse momentum reconstruction when the reconstructed particle is the best photon match
  TH1F *resTrackLink_pt = new TH1F("rTck_pt", "Reco track", 40, -2000, 2000); //transverse momentum from individual track reconstruction segments

  //histograms for storing relative differences between values
  TH1F *relResEPSum_pt = new TH1F("rrEP_pt", "Summed electrons and protons", 40, -2, 2); //relative resolution of transverse momentum reconstruction when the reconstructed particle is made up of all reconstructed electrons and photons
  TH1F *relResEPSum_PA = new TH1F("rrEP_PA", "Summed electrons and protons", 20, -0.05, 0.05); //relative resolution of polar angle reconstruction when the reconstructed particle is made up of all reconstructed electrons and photons
  TH1F *relResEPSum_azimuth = new TH1F("rrEP_A", "Summed electrons and protons", 20, -0.2, 0.2); //relative resolution of azimuth reconstruction when the reconstructed particle is made up of all reconstructed electrons and photons

  //histograms for storing attributes unique to certain reconstructed types of particles
  //transverse momentum
  TH1F *electronReco_pt = new TH1F("elRec_pt", "Electron pt", 20, 0, 2000);
  TH1F *photonReco_pt = new TH1F("phRec_pt", "Photon pt", 20, 0, 2000);
  TH1F *neutronReco_pt = new TH1F("nuRec_pt", "Neutron pt", 20, 0, 2000);
  TH1F *pionReco_pt = new TH1F("piRec_pt", "Pion pt", 20, 0, 2000);

  //polar angle
  TH1F *electronReco_PA = new TH1F("elRec_PA", "Electron PA", 20, 0, 3.2);
  TH1F *photonReco_PA = new TH1F("phRec_PA", "Photon PA", 20, 0, 3.2);
  TH1F *neutronReco_PA = new TH1F("nuRec_PA", "Neutron PA", 20, 0, 3.2);
  TH1F *pionReco_PA = new TH1F("piRec_PA", "Pion PA", 20, 0, 3.2);

  //azimuth
  TH1F *electronReco_azimuth = new TH1F("elRec_a", "Electron azimuth", 20, -3.2, 3.2);
  TH1F *photonReco_azimuth = new TH1F("phRec_a", "Photon azimuth", 20, -3.2, 3.2);
  TH1F *neutronReco_azimuth = new TH1F("nuRec_a", "Neutron azimuth", 20, -3.2, 3.2);
  TH1F *pionReco_azimuth = new TH1F("piRec_a", "Pion azimuth", 20, -3.2, 3.2);

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
  TTreeReaderArray<Float_t> tsome_RA(myReader, "tsome"); //reco track curviture
  TTreeReaderArray<int> rctyp_RA(myReader, "rctyp"); //reconstructed type (11 = electron, 22 = photon, 2112 = neutron)
  TTreeReaderArray<int> mcpdg_RA(myReader, "mcpdg"); //truth type (see list above)
  TTreeReaderArray<int> mcgst_RA(myReader, "mcgst"); //is this particle a generating particle or not

  //temporary variables that will reduce the number of get operatons
  int r2fTemp, r2tTemp;
  bool electronMatch = false;
  bool photonMatch = false;
  bool pionMatch = false;
  Float_t r2wTemp, dPATemp, dATemp;
  TVector3 mcmoTemp, rcmoTemp, rcmoSumTemp, BMTemp, BMPTemp;
  Float_t r2wMax = 0;
  int r2wMax_Index = -1;
  int count = 0;

  /////////////////////
  //HISTOGRAM FILLING//
  /////////////////////
  
  //loops over each event
  while(myReader.Next()) {

    //loops over all particles to find the generating particle
    for(int i = 0; i < mcmox_RA.GetSize(); i++) {
      if(mcgst_RA.At(i)) {

	mcmoTemp.SetXYZ(mcmox_RA.At(i), mcmoy_RA.At(i), mcmoz_RA.At(i));

	//adds generating particle to histogram
	realAll_pt->Fill(mcmoTemp.Perp());
	realAll_PA->Fill(mcmoTemp.Theta());
	realAll_azimuth->Fill(mcmoTemp.Phi());
	break; //this might cause issues in the future, double check to make sure that there are not two or more generating particles
      }
    }

    //resets data sum variables
    rcmoTemp.SetXYZ(0, 0, 0);
    BMTemp.SetXYZ(0, 0, 0);
    rcmoSumTemp.SetXYZ(0, 0, 0);

    //loops over all reconstructed particles
    for(int i = 0; i < rcmox_RA.GetSize(); i++) {
      rcmoTemp.SetXYZ(rcmox_RA.At(i), rcmoy_RA.At(i), rcmoz_RA.At(i));

      //if that particle is an electron or photon, add it the sum variables (this accounts for particle mis-identification)
      if(abs(rctyp_RA.At(i)) == 11 || abs(rctyp_RA.At(i)) == 22 || abs(rctyp_RA.At(i)) == 211) {
	for(int j = 0; j < r2f_RA.GetSize(); j++) {
	  if(r2f_RA.At(j) == i && mcgst_RA.At(r2t_RA.At(j))) {
	    rcmoSumTemp.SetXYZ(rcmoSumTemp.X() + rcmoTemp.X(), rcmoSumTemp.Y() + rcmoTemp.Y(), rcmoSumTemp.Z() + rcmoTemp.Z());
	    break;
	  }
	}
      }

      //if the particle is a elect ron, add it to one set of histograms
      if(abs(rctyp_RA.At(i)) == 11) {
	electronReco_pt->Fill(rcmoTemp.Perp());
	electronReco_PA->Fill(rcmoTemp.Theta());
	electronReco_azimuth->Fill(rcmoTemp.Phi());

	//if the particle is a photon, add it to another set of histograms
      } else if(abs(rctyp_RA.At(i)) == 22) {
	photonReco_pt->Fill(rcmoTemp.Perp());
	photonReco_PA->Fill(rcmoTemp.Theta());
	photonReco_azimuth->Fill(rcmoTemp.Phi());

	//if the particle is a neutron, add it to a third set of histograms
      } else if(abs(rctyp_RA.At(i)) == 2112) {
	neutronReco_pt->Fill(rcmoTemp.Perp());
	neutronReco_PA->Fill(rcmoTemp.Theta());
	neutronReco_azimuth->Fill(rcmoTemp.Phi());

	//if the particle is a pion, add it to a fourth set of histograms
      } else if(abs(rctyp_RA.At(i)) == 211) {
	pionReco_pt->Fill(rcmoTemp.Perp());
	pionReco_PA->Fill(rcmoTemp.Theta());
	pionReco_azimuth->Fill(rcmoTemp.Phi());
      }
    }

    //searches for the largest weight of the relationship that links a reconstructed particle to a particle that we know all the details of
    
    for(int i = 0; i < r2f_RA.GetSize(); i++) {
      r2tTemp = r2t_RA.At(i);
      r2fTemp = r2f_RA.At(i);
      r2wTemp = r2w_RA.At(i);

      //checks to see if those particles are electrons and if the paticle is an originial generating particle
      if(abs(mcpdg_RA.At(r2tTemp)) == 11 && abs(rctyp_RA.At(r2fTemp)) == 11 && mcgst_RA.At(r2tTemp)) {
	if(r2wTemp > r2wMax) {
	  electronMatch = true;
	  r2wMax = r2wTemp;
	  r2wMax_Index = i;
	}
      }
    }
    
    if(r2wMax_Index != -1) {
      //adds the weight to a histogram
      realElWeights_Hist->Fill(r2wMax);
    }

    //if a relationship was found between those particles, it is added to histograms
    if(r2wMax_Index != -1) {

      //find the particles that the particle linker refers to
      r2tTemp = r2t_RA.At(r2wMax_Index);
      r2fTemp = r2f_RA.At(r2wMax_Index);
      r2wTemp = r2w_RA.At(r2wMax_Index);

      //finds the reco momenta of the particle
      rcmoTemp.SetXYZ(rcmox_RA.At(r2fTemp), rcmoy_RA.At(r2fTemp), rcmoz_RA.At(r2fTemp));
      BMTemp.SetXYZ(rcmox_RA.At(r2fTemp), rcmoy_RA.At(r2fTemp), rcmoz_RA.At(r2fTemp));
      
      //finds the differences between truth and reconstructed attributes
      dPATemp = mcmoTemp.Theta() - rcmoTemp.Theta();
      dATemp = mcmoTemp.Phi() - rcmoTemp.Phi();

      //fills in histograms with attributes
      realPassed_pt->Fill(mcmoTemp.Perp());
      realPassed_PA->Fill(mcmoTemp.Theta());
      realPassed_azimuth->Fill(mcmoTemp.Phi());
      
      //fills in histograms with differences
      diffPt_Hist->Fill(mcmoTemp.Perp() - rcmoTemp.Perp());
      diffPA_Hist->Fill(mcmoTemp.Theta() - rcmoTemp.Theta());
      diffA_Hist->Fill(mcmoTemp.Phi() - rcmoTemp.Phi());
      
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
	if(abs(mcpdg_RA.At(r2tTemp)) == 11 && (abs(rctyp_RA.At(r2fTemp)) == 22 || abs(rctyp_RA.At(r2fTemp)) == 211) && mcgst_RA.At(r2tTemp)) {
	  if(r2wTemp > r2wMax) {
	    r2wMax = r2wTemp;
	    r2wMax_Index = i;
	  }
	}
      }
      if(r2wMax_Index != -1) {	
	if(rctyp_RA.At(r2wMax_Index) == 22) {
	  photonMatch = true;
	} else if(rctyp_RA.At(r2wMax_Index) == 211) {
	  pionMatch = true;
	}
      }
    }

    //repeats the previous steps that adds data to histograms but instead adds it to the histograms that holds links between photons and electrons as well
    if(photonMatch || electronMatch || pionMatch) {
      r2tTemp = r2t_RA.At(r2wMax_Index);
      r2fTemp = r2f_RA.At(r2wMax_Index);
      r2wTemp = r2w_RA.At(r2wMax_Index);

      //adds the weight of the link to a histogram
      realAllWeights_Hist->Fill(r2wMax);

      rcmoTemp.SetXYZ(rcmox_RA.At(r2fTemp), rcmoy_RA.At(r2fTemp), rcmoz_RA.At(r2fTemp));

      realAllPassed_pt->Fill(mcmoTemp.Perp());
      realAllPassed_PA->Fill(mcmoTemp.Theta());
      realAllPassed_azimuth->Fill(mcmoTemp.Phi());

      //finds the difference between the truth and best match reconstructed particles
      resBestMatch_pt->Fill(mcmoTemp.Perp() - rcmoTemp.Perp());
      resBestMatch_PA->Fill(mcmoTemp.Theta() - rcmoTemp.Theta());
      resBestMatch_azimuth->Fill(mcmoTemp.Phi() - rcmoTemp.Phi());

      if(photonMatch) {
	resPhotonReco_pt->Fill(mcmoTemp.Perp() - rcmoTemp.Perp());
	resPhotonReco_PA->Fill(mcmoTemp.Theta() - rcmoTemp.Theta());
	resPhotonReco_azimuth->Fill(mcmoTemp.Phi() - rcmoTemp.Phi());

	if(tsome_RA.GetSize() > 0){
	  recoTrackLink_pt->Fill(((0.3 * 3.57) / tsome_RA.At(0)) / 1000);
	  resTrackLink_pt->Fill(mcmoTemp.Perp() - abs(((0.3 * 3.57) / tsome_RA.At(0)) / 1000));
	}
      }
      if(pionMatch) {
	resPionReco_pt->Fill(mcmoTemp.Perp() - rcmoTemp.Perp());
	resPionReco_PA->Fill(mcmoTemp.Theta() - rcmoTemp.Theta());
	resPionReco_azimuth->Fill(mcmoTemp.Phi() - rcmoTemp.Phi());
      }
    }
    
    count++;
    r2wMax = 0;
    r2wMax_Index = -1;
    electronMatch = false;
    photonMatch = false;
    pionMatch = false;

    //finds the difference between the truth and reconstructed particles (sum of photons and electrons)
    if(rcmoSumTemp.X() != 0 || rcmoSumTemp.Y() != 0 || rcmoSumTemp.Z() != 0){

      //fills in the electron and proton sum histograms
      realEPSum_pt->Fill(rcmoSumTemp.Perp());
      realEPSum_PA->Fill(rcmoSumTemp.Theta());
      realEPSum_azimuth->Fill(rcmoSumTemp.Phi());
      
      //fills in the electron and proton sum resolution histograms
      resEPSum_pt->Fill(mcmoTemp.Perp() - rcmoSumTemp.Perp());
      resEPSum_PA->Fill(mcmoTemp.Theta() - rcmoSumTemp.Theta());
      resEPSum_azimuth->Fill(mcmoTemp.Phi() - rcmoSumTemp.Phi());

      //fills in the electron and proton sum relative resolution histograms
      relResEPSum_pt->Fill((mcmoTemp.Perp() - rcmoSumTemp.Perp()) / mcmoTemp.Perp());
      relResEPSum_PA->Fill((mcmoTemp.Theta() - rcmoSumTemp.Theta()) / mcmoTemp.Theta());
      relResEPSum_azimuth->Fill((mcmoTemp.Phi() - rcmoSumTemp.Phi()) / mcmoTemp.Phi());
    }    
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
  pionReco_pt->SetLineColor(kGreen);

  hs->Add(electronReco_pt);
  hs->Add(photonReco_pt);
  hs->Add(neutronReco_pt);
  hs->Add(pionReco_pt);
  
  hs->Draw("nostack");

  gPad->SetGrid(1, 0);
  gPad->BuildLegend(0.75, 0.75, 0.9, 0.9, "");
  
  c->SaveAs(saveDir + "/allPart_pt.png");
  c->Close();
  c = new TCanvas();

  //draws the polar angle histograms for electrons, photons, and neutrons on top of each other
  hs = new THStack("hs", "Polar angle of electrons, photons, and neutrons;Polar angle (Rads);Count");

  c->SetLogy();

  electronReco_PA->SetLineColor(kBlue);
  photonReco_PA->SetLineColor(kRed);
  neutronReco_PA->SetLineColor(kBlack);
  pionReco_PA->SetLineColor(kGreen);

  hs->Add(electronReco_PA);
  hs->Add(photonReco_PA);
  hs->Add(neutronReco_PA);
  hs->Add(pionReco_PA);
  
  hs->Draw("nostack");

  gPad->SetGrid(1, 0);
  gPad->BuildLegend(0.7, 0.1, 0.9, 0.3, "");
  
  c->SaveAs(saveDir + "/allPart_PA.png");
  c->Close();
  c = new TCanvas();

  //draws the azimuth histograms for electrons, photons, and neutrons on top of each other
  hs = new THStack("hs", "Azimuth of electrons, photons, and neutrons;Azimuth (Rads);Count");

  c->SetLogy();

  electronReco_azimuth->SetLineColor(kBlue);
  photonReco_azimuth->SetLineColor(kRed);
  neutronReco_azimuth->SetLineColor(kBlack);
  pionReco_azimuth->SetLineColor(kGreen);

  hs->Add(electronReco_azimuth);
  hs->Add(photonReco_azimuth);
  hs->Add(neutronReco_azimuth);
  hs->Add(pionReco_azimuth);

  hs->Draw("nostack");

  gPad->SetGrid(1, 0);
  gPad->BuildLegend(0.4, 0.1, 0.6, 0.3, "");

  c->SaveAs(saveDir + "/allPart_azimuth.png");
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
  
  c->SaveAs(saveDir + "/compLink_pt.png");
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

  c->SaveAs(saveDir + "/compLink_PA.png");
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

  c->SaveAs(saveDir + "/compLink_azimuth.png");
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

  c->SaveAs(saveDir + "/compAllLink_pt.png");
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

  c->SaveAs(saveDir + "/compAllLink_PA.png");
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

  c->SaveAs(saveDir + "/compAllLink_azimuth.png");
  c->Close();
  c = new TCanvas();

  //=======================================//
  //stacked reconstruction resolution hists//
  //=======================================//

  //draws the resolution of the transverse momentum of different linking methods
  hs = new THStack("hs", "Resolution of different linking methods (no BIB);Transverse momentum (GeV);Count");

  diffPt_Hist->SetLineColor(kBlue);
  resEPSum_pt->SetLineColor(kRed);
  resBestMatch_pt->SetLineColor(kBlack);

  hs->Add(diffPt_Hist);
  hs->Add(resEPSum_pt);
  hs->Add(resBestMatch_pt);

  hs->Draw("nostack");

  gPad->SetGrid(1, 0);
  gPad->BuildLegend(0.1, 0.7, 0.4, 0.9, "");

  c->SaveAs(saveDir + "/resolution/compResLink_pt.png");
  c->Close();
  c = new TCanvas();

  //draws the resolution of the polar angle of different linking methods
  hs = new THStack("hs", "Resolution of different linking methods (no BIB);Polar angle (GeV);Count");

  diffPA_Hist->SetLineColor(kBlue);
  resEPSum_PA->SetLineColor(kRed);
  resBestMatch_PA->SetLineColor(kBlack);

  hs->Add(diffPA_Hist);
  hs->Add(resEPSum_PA);
  hs->Add(resBestMatch_PA);

  hs->Draw("nostack");

  gPad->SetGrid(1, 0);
  gPad->BuildLegend(0.1, 0.7, 0.4, 0.9, "");

  c->SaveAs(saveDir + "/resolution/compResLink_PA.png");
  c->Close();
  c = new TCanvas();

  //draws the resolution of the azimuth of different linking methods
  hs = new THStack("hs", "Resolution of different linking methods (no BIB);Azimuth (GeV);Count");

  diffA_Hist->SetLineColor(kBlue);
  resEPSum_azimuth->SetLineColor(kRed);
  resBestMatch_azimuth->SetLineColor(kBlack);

  hs->Add(diffA_Hist);
  hs->Add(resEPSum_azimuth);
  hs->Add(resBestMatch_azimuth);

  hs->Draw("nostack");

  gPad->SetGrid(1, 0);
  gPad->BuildLegend(0.1, 0.7, 0.4, 0.9, "");

  c->SaveAs(saveDir + "/resolution/compResLink_azimuth.png");
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
  
  c->SaveAs(saveDir + "/weight/recoLinkElWeights.png");
  c->Close();
  c = new TCanvas();

  //draws the histogram that shows the distrobution of the relation weights for all particles
  realAllWeights_Hist->Draw();
  realAllWeights_Hist->SetTitle("Weights of linked particles");
  realAllWeights_Hist->GetXaxis()->SetTitle("Weight");
  realAllWeights_Hist->GetYaxis()->SetTitle("Count");

  c->SaveAs(saveDir + "/weight/recoLinkAllWeights.png");
  c->Close();
  c = new TCanvas();

  //draws a hist that hold the reco track link data for best match photons
  recoTrackLink_pt->Draw();
  recoTrackLink_pt->SetTitle("Reconstructed track link transverse momentum for best match photons");
  recoTrackLink_pt->GetXaxis()->SetTitle("Truth - reconstructed (Gev)");
  recoTrackLink_pt->GetYaxis()->SetTitle("Count");

  c->SaveAs(saveDir + "/recoTrackLink_pt.png");
  c->Close();
  c = new TCanvas();

  //============================//
  //Resolution of electron linking//
  //============================//

  //draws a hist that hold the differences in electron pt
  diffPt_Hist->Draw();
  diffPt_Hist->SetTitle("Difference in pt between all and reconstructed particles");
  diffPt_Hist->GetXaxis()->SetTitle("Transverse momentum (GeV)");
  diffPt_Hist->GetYaxis()->SetTitle("Count");

  c->SaveAs(saveDir + "/resolution/recoPt_diff.png");
  c->Close();
  c = new TCanvas();

  //draws a hist that hold the differences in electron polar angle
  diffPA_Hist->Draw();
  diffPA_Hist->SetTitle("Difference in polar angle between all and reconstructed particles");
  diffPA_Hist->GetXaxis()->SetTitle("Polar angle (Rads)");
  diffPA_Hist->GetYaxis()->SetTitle("Count");

  c->SaveAs(saveDir + "/resolution/recoPA_diff.png");
  c->Close();
  c = new TCanvas();

  //draws a hist that hold the differences in electron azimuth
  diffA_Hist->Draw();
  diffA_Hist->SetTitle("Difference in azimuth between all and reconstructed particles");
  diffA_Hist->GetXaxis()->SetTitle("Azimuth (Rads)");
  diffA_Hist->GetYaxis()->SetTitle("Count");

  c->SaveAs(saveDir + "/resolution/recoAzimuth_diff.png");
  c->Close();
  c = new TCanvas();

  //=============================================//
  //summed data from linked electrons and photons//
  //=============================================//

  //draws a hist that hold the electron and photon sum pt resolution
  resEPSum_pt->Draw();
  resEPSum_pt->SetTitle("All truth particle and sum of electrons and photons pt resolution");
  resEPSum_pt->GetXaxis()->SetTitle("Truth - reconstructed p_{T} (GeV)");
  resEPSum_pt->GetYaxis()->SetTitle("Count");

  c->SaveAs(saveDir + "/resolution/resEPSumLink_pt.png");
  c->Close();
  c = new TCanvas();
  
  //draws a hist that hold the electron and photon sum polar angle resolution
  resEPSum_PA->Draw();
  resEPSum_PA->SetTitle("All truth particle and sum of electrons and photons polar angle resolution");
  resEPSum_PA->GetXaxis()->SetTitle("Truth - reconstructed polar angle (Rads)");
  resEPSum_PA->GetYaxis()->SetTitle("Count");

  c->SaveAs(saveDir + "/resolution/resEPSumLink_PA.png");
  c->Close();
  c = new TCanvas();
  
  //draws a hist that hold the electron and photon sum azimuth resolution
  resEPSum_azimuth->Draw();
  resEPSum_azimuth->SetTitle("All truth particle and sum of electrons and photons azimuth resolution");
  resEPSum_azimuth->GetXaxis()->SetTitle("Truth - reconstructed azimuth (Rads)");
  resEPSum_azimuth->GetYaxis()->SetTitle("Count");

  c->SaveAs(saveDir + "/resolution/resEPSumLink_azimuth.png");
  c->Close();
  c = new TCanvas();

  //==========================//
  //Best match resolution data//
  //==========================//

  //draws a hist that hold the truth particle and the best match link regardless of reco type pt resolution
  resBestMatch_pt->Draw();
  resBestMatch_pt->SetTitle("Truth and best match link pt resolution");
  resBestMatch_pt->GetXaxis()->SetTitle("Truth - reconstructed p_{T} (GeV)");
  resBestMatch_pt->GetYaxis()->SetTitle("Count");

  c->SaveAs(saveDir + "/resolution/resBestMatchLink_pt.png");
  c->Close();
  c = new TCanvas();

  //draws a hist that hold the truth particle and the best match link regardless of reco type polar angle resolution
  resBestMatch_PA->Draw();
  resBestMatch_PA->SetTitle("Truth and best match link polar angle resolution");
  resBestMatch_PA->GetXaxis()->SetTitle("Truth - reconstructed polar angle (Rads)");
  resBestMatch_PA->GetYaxis()->SetTitle("Count");

  c->SaveAs(saveDir + "/resolution/resBestMatchLink_PA.png");
  c->Close();
  c = new TCanvas();

  //draws a hist that hold the truth particle and the best match link regardless of reco type azimuth resolution
  resBestMatch_azimuth->Draw();
  resBestMatch_azimuth->SetTitle("Truth and best match link azimuth resolution");
  resBestMatch_azimuth->GetXaxis()->SetTitle("Truth - reconstructed azimuth (Rads)");
  resBestMatch_azimuth->GetYaxis()->SetTitle("Count");

  c->SaveAs(saveDir + "/resolution/resBestMatchLink_azimuth.png");
  c->Close();
  c = new TCanvas();
  
  //=============================================//
  //Best match resolution data for photon matches//
  //=============================================//

  //draws a hist that hold the truth particle and the best match link pt resolution for photons
  resPhotonReco_pt->Draw();
  resPhotonReco_pt->SetTitle("Truth and best match link pt resolution for photons");
  resPhotonReco_pt->GetXaxis()->SetTitle("Truth - reconstructed p_{T} (GeV)");
  resPhotonReco_pt->GetYaxis()->SetTitle("Count");

  c->SaveAs(saveDir + "/resolution/resBMPhotonLink_pt.png");
  c->Close();
  c = new TCanvas();

  //draws a hist that hold the truth particle and the best match link polar angle resolution for photons
  resPhotonReco_PA->Draw();
  resPhotonReco_PA->SetTitle("Truth and best match link polar angle resolution for photons");
  resPhotonReco_PA->GetXaxis()->SetTitle("Truth - reconstructed polar angle (Rads)");
  resPhotonReco_PA->GetYaxis()->SetTitle("Count");

  c->SaveAs(saveDir + "/resolution/resBMPhotonLinkLink_PA.png");
  c->Close();
  c = new TCanvas();

  //draws a hist that hold the truth particle and the best match link azimuth resolution for photons
  resPhotonReco_azimuth->Draw();
  resPhotonReco_azimuth->SetTitle("Truth and best match link azimuth resolution for photons");
  resPhotonReco_azimuth->GetXaxis()->SetTitle("Truth - reconstructed azimuth (Rads)");
  resPhotonReco_azimuth->GetYaxis()->SetTitle("Count");

  c->SaveAs(saveDir + "/resolution/resBMPhotonLink_azimuth.png");
  c->Close();
  c = new TCanvas();

  //=============================================//
  //Best match resolution data for pion matches//
  //=============================================//

  //draws a hist that hold the truth particle and the best match link pt resolution for pions
  resPionReco_pt->Draw();
  resPionReco_pt->SetTitle("Truth and best match link pt resolution for pions");
  resPionReco_pt->GetXaxis()->SetTitle("Truth - reconstructed p_{T} (GeV)");
  resPionReco_pt->GetYaxis()->SetTitle("Count");

  c->SaveAs(saveDir + "/resolution/resBMPionLink_pt.png");
  c->Close();
  c = new TCanvas();

  //draws a hist that hold the truth particle and the best match link polar angle resolution for pions
  resPionReco_PA->Draw();
  resPionReco_PA->SetTitle("Truth and best match link polar angle resolution for pions");
  resPionReco_PA->GetXaxis()->SetTitle("Truth - reconstructed polar angle (Rads)");
  resPionReco_PA->GetYaxis()->SetTitle("Count");

  c->SaveAs(saveDir + "/resolution/resBMPionLinkLink_PA.png");
  c->Close();
  c = new TCanvas();

  //draws a hist that hold the truth particle and the best match link azimuth resolution for pions
  resPionReco_azimuth->Draw();
  resPionReco_azimuth->SetTitle("Truth and best match link azimuth resolution for pions");
  resPionReco_azimuth->GetXaxis()->SetTitle("Truth - reconstructed azimuth (Rads)");
  resPionReco_azimuth->GetYaxis()->SetTitle("Count");

  c->SaveAs(saveDir + "/resolution/resBMPionLink_azimuth.png");
  c->Close();
  c = new TCanvas();
  
  //======================================================//
  //Resolution of track link reconstruction for BM photons//
  //======================================================//

  //draws a hist that hold the truth particle and the best match photon track link pt  resolution
  resTrackLink_pt->Draw();
  resTrackLink_pt->SetTitle("Truth and best match link photon track transverse momentum resolution");
  resTrackLink_pt->GetXaxis()->SetTitle("Truth - reconstructed track p_{T} (GeV)");
  resTrackLink_pt->GetYaxis()->SetTitle("Count");

  c->SaveAs(saveDir + "/resolution/resTrackLink_pt.png");
  c->Close();
  c = new TCanvas();
  
  //================================================//
  //relative resolution summed protons and electrons//
  //================================================//

  //draws a hist that hold the truth particle and the best match link regardless of reco type pt relative resolution
  relResEPSum_pt->Draw();
  relResEPSum_pt->SetTitle("Truth and best match link pt relative resolution");
  relResEPSum_pt->GetXaxis()->SetTitle("(Truth - reconstructed) / truth p_{T}");
  relResEPSum_pt->GetYaxis()->SetTitle("Count");

  c->SaveAs(saveDir + "/resolution/relResEPSumLink_pt.png");
  c->Close();
  c = new TCanvas();

  //draws a hist that hold the truth particle and the best match link regardless of reco type polar angle relative resolution
  relResEPSum_PA->Draw();
  relResEPSum_PA->SetTitle("Truth and best match link polar angle relative resolution");
  relResEPSum_PA->GetXaxis()->SetTitle("(Truth - reconstructed) / truth polar angle");
  relResEPSum_PA->GetYaxis()->SetTitle("Count");

  c->SaveAs(saveDir + "/resolution/relResEPSumLink_PA.png");
  c->Close();
  c = new TCanvas();

  //draws a hist that hold the truth particle and the best match link regardless of reco type azimuth relative resolution
  relResEPSum_azimuth->Draw();
  relResEPSum_azimuth->SetTitle("Truth and best match link azimuth relative resolution");
  relResEPSum_azimuth->GetXaxis()->SetTitle("(Truth - reconstructed) / truth azimuth");
  relResEPSum_azimuth->GetYaxis()->SetTitle("Count");

  c->SaveAs(saveDir + "/resolution/relResEPSumLin_azimuth.png");
  c->Close();
  c = new TCanvas();

  ////////
  //MISC//
  ////////

  //draws a histogram for delta r (distance between polar angle and lambda on a plot between the real and reco data)
  c->SetLogy();
    
  deltaR_Hist->Draw();
  deltaR_Hist->SetTitle("Delta r between the real and reco data points");
  deltaR_Hist->GetXaxis()->SetTitle("delta R (Rads)");
  deltaR_Hist->GetYaxis()->SetTitle("Count");

  c->SaveAs(saveDir + "/deltaR.png");
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

  c->SaveAs(saveDir + "/efficiency/Eff_pt.png");
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

  c->SaveAs(saveDir + "/efficiency/Eff_PA.png");
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

  c->SaveAs(saveDir + "/efficiency/Eff_azimuth.png");

  c->Clear();
  c->Close();
  c = new TCanvas();

  //draws an efficiency graph for transverse momentum for best match particles
  TEfficiency* ptAll_Eff = new TEfficiency("ptBM_Eff", "Efficiency of transverse momentum reconstruction for all particles;Transverse momentum (GeV);Efficiency", 20, 0, 2000);
  eff_File = new TFile("effFile.root", "recreate");
  
  //checks to make sure if the two histograms are compatable to make an efficiency histogram
  if(TEfficiency::CheckConsistency(*realAllPassed_pt, *realAll_pt)) {
    ptAll_Eff->SetPassedHistogram(*realAllPassed_pt, "f");
    ptAll_Eff->SetTotalHistogram(*realAll_pt, "f");
    eff_File->Write();
  }

  c->Clear();

  ptAll_Eff->Draw();

  c->SaveAs(saveDir + "/efficiency/EffAll_pt.png");
  c->Close();
  c = new TCanvas();

  //draws an efficiency graph for the polar angle for best match particles
  TEfficiency* PAAll_Eff = new TEfficiency("PABM_Eff", "Efficiency of polar angle reconstruction for all particles;Polar angle (Rads);Efficiency", 20, 0, 1.6);
  eff_File = new TFile("effFile.root", "recreate");

  //checks to make sure if the two histograms are compatable to make an efficiency histogram
  if(TEfficiency::CheckConsistency(*realAllPassed_PA, *realAll_PA)) {
    PAAll_Eff->SetPassedHistogram(*realAllPassed_PA, "f");
    PAAll_Eff->SetTotalHistogram(*realAll_PA, "f");
    eff_File->Write();
  }

  c->Clear();

  PAAll_Eff->Draw();

  c->SaveAs(saveDir + "/efficiency/EffAll_PA.png");
  c->Close();
  c = new TCanvas();

  //draws an efficiency graph for the azimuth for best match particles
  TEfficiency* azimuthAll_Eff = new TEfficiency("azimuthBM_Eff", "Efficiency of azimuth reconstruction for all particles;Azimuth (Rads);Eff\
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

  c->SaveAs(saveDir + "/efficiency/EffAll_azimuth.png");
  c->Close();

  //setting new line colors so the plots are differentable when put on the same graph
  pt_Eff->SetLineColor(4);
  PA_Eff->SetLineColor(4);
  azimuth_Eff->SetLineColor(4);
  
  ptAll_Eff->SetLineColor(2);
  PAAll_Eff->SetLineColor(2);
  azimuthAll_Eff->SetLineColor(2);

  pt_Eff->SetTitle("Just electrons");
  PA_Eff->SetTitle("Just electrons");
  azimuth_Eff->SetTitle("Just electrons");

  ptAll_Eff->SetTitle("Electrons and photons");
  PAAll_Eff->SetTitle("Electrons and photons");
  azimuthAll_Eff->SetTitle("Electrons and photons");

  //creates a multigraph for the transverse momentum efficiencies
  TMultiGraph *mg = new TMultiGraph();
  c = new TCanvas();
  
  mg->Add(pt_Eff->CreateGraph());
  mg->Add(ptAll_Eff->CreateGraph());

  mg->SetTitle("Transverse momentum reconstruction efficiencies");
  mg->GetXaxis()->SetTitle("Transverse momentum (GeV)");
  mg->GetYaxis()->SetTitle("Efficiency");

  mg->Draw("aZ");

  gPad->SetGrid(1, 0);
  gPad->BuildLegend(0.1, 0.1, 0.3, 0.25, "");

  c->SaveAs(saveDir + "/efficiency/multiBMEff_pt.png");
  c->Close();

  //creates a multigraph for the polar angle efficiencies
  mg = new TMultiGraph();
  c = new TCanvas();

  mg->Add(PA_Eff->CreateGraph());
  mg->Add(PAAll_Eff->CreateGraph());

  mg->SetTitle("Polar angle reconstruction efficiencies");
  mg->GetXaxis()->SetTitle("Polar angle (Rads)");
  mg->GetYaxis()->SetTitle("Efficiency");
  
  mg->Draw("aZ");

  gPad->SetGrid(1, 0);
  gPad->BuildLegend(0.1, 0.1, 0.3, 0.25, "");

  c->SaveAs(saveDir + "/efficiency/multiBMEff_PA.png");
  c->Close();

  //creates a multigraph for the azimuth efficiencies
  mg = new TMultiGraph();
  c = new TCanvas();

  mg->Add(azimuth_Eff->CreateGraph());
  mg->Add(azimuthAll_Eff->CreateGraph());

  mg->SetTitle("Azimuth reconstruction efficiencies");
  mg->GetXaxis()->SetTitle("Azimuth (Rads)");
  mg->GetYaxis()->SetTitle("Efficiency");
  
  mg->Draw("aZ");

  gPad->SetGrid(1, 0);
  gPad->BuildLegend(0.1, 0.1, 0.3, 0.25, "");

  c->SaveAs(saveDir + "/efficiency/multiBMEff_azimuth.png");
  c->Close();

}
