void filter_compare()
{
  /////////////////////////////
  //GENERAL ANALYSIS SETTINGS//
  /////////////////////////////
  
  auto fileName1 = "10k-ACTS-Filter6.root";
  auto treeName1 = "MyLCTuple";
  TString filePrefix1 = "Filter6";

  auto fileName2 = "10k-ACTS-Filter4.root";
  auto treeName2 = "MyLCTuple";
  TString filePrefix2 = "Filter4";

  TString saveDir = "10k-Compare-Filter6+4";

  /////////////////////////////

  TFile *myFile1 = new TFile(fileName1);
  TTree *myTree1 = (TTree*)myFile1->Get(treeName1);

  TFile *myFile2 = new TFile(fileName2);
  TTree *myTree2 = (TTree*)myFile2->Get(treeName2);

  TCanvas *c = new TCanvas();
  
  //opens the file to be read
  auto openFile1 = TFile::Open(fileName1);

  auto openFile2 = TFile::Open(fileName2);

  /////////////////////////
  //Histogram definitions//
  /////////////////////////

  //FILE 1
  //histograms for attributes of electrons
  //transverse momentum
  TH1F *realPassed_pt1 = new TH1F("rP_pt1", filePrefix1 + ":Linked electrons", 40, 0, 2000); //transverse momentum of electrons that are linked to reconstructed electrons
  TH1F *realAllPassed_pt1 = new TH1F("AP_pt1", filePrefix1 + ":Best match", 40, 0, 2000); //best match transverse momentum from linked electrons and photons
  TH1F *realEPSum_pt1 = new TH1F("EPS_pt1", filePrefix1 + ":Summed particles", 40, 0, 2000); //summed transverse momentum from linked electrons and photons
  TH1F *recoTrackLink_pt1 = new TH1F("recTck_pt1", filePrefix1 + ":Reco track", 40, -2000, 2000); //transverse momentum from individual track reconstruction segments
  TH1F *realAll_pt1 = new TH1F("rA_pt1", filePrefix1 + ":All electrons", 40, 0, 2000); //transverse momentum of the truth particles
  //Polar angle
  TH1F *realPassed_PA1 = new TH1F("rP_PA1", filePrefix1 + ":Linked electrons", 40, 0, 3.2); //see above
  TH1F *realAllPassed_PA1 = new TH1F("AP_PA1", filePrefix1 + ":Best match", 40, 0, 3.2);
  TH1F *realEPSum_PA1 = new TH1F("EPS_PA1", filePrefix1 + ":Summed particles", 40, 0, 3.2);
  TH1F *realAll_PA1 = new TH1F("rA_PA1", filePrefix1 + ":All electrons", 40, 0, 3.2);
  //Azimuth
  TH1F *realPassed_azimuth1 = new TH1F("rP_a1", filePrefix1 + ":Linked electrons", 40, -3.2, 3.2); //see above
  TH1F *realAllPassed_azimuth1 = new TH1F("AP_a1", filePrefix1 + ":Best match", 40, -3.2, 3.2);
  TH1F *realEPSum_azimuth1 = new TH1F("EPS_a1", filePrefix1 + ":Summed particles", 40, -3.2, 3.2);
  TH1F *realAll_azimuth1 = new TH1F("rA_a1", filePrefix1 + ":All electrons", 40, -3.2, 3.2);

  //histograms for weights of particles
  TH1D *recoWeight_El1 = new TH1D("realElWeights1", filePrefix1 + ":Reco link electron weights", 26, 0, 1.3); //link weights of reconstructed electrons linked to truth electrons
  TH1D *recoWeight_BM1 = new TH1D("realAllWeights1", filePrefix1 + ":Reco link all weights", 26, 0, 1.3); // link weights of particles linked to truth electrons

  //histogam for the number of hits of best match particles
  TH1D *numHits_BM1 = new TH1D("nHitBM1", filePrefix1 + ":Number of BM hits", 30, 0, 30); //number of hits that the track encounters

  //histograms for storing differences between values
  TH1F *relResOG_pt1 = new TH1F("diffPt1", filePrefix1 + ":Simple electron linking", 40, -1, 1); //difference between reconstructed linked particles and the truth particles in transverse momentum
  TH1F *relResOG_PA1 = new TH1F("diffPA1", filePrefix1 + ":Simple electron linking", 40, -1, 1); //difference between reconstructed linked particles and the truth particles in transverse momentum
  TH1F *relResOG_azimuth1 = new TH1F("diffA1", filePrefix1 + ":Simple electron linking", 40, -1, 1); //difference between reconstructed linked particles and the truth particles in transverse momentum
  TH1F *deltaR_Hist1 = new TH1F("deltaR1", filePrefix1 + ":MyLCTuple", 25, 0, 0.0025); //delta R (delta polar angle and azimuth combined)
  TH1F *relResEPSum_pt1 = new TH1F("rEP_pt1", filePrefix1 + ":Electron + photon sum", 40, -1, 1); //resolution of transverse momentum reconstruction when the reconstructed particle is made up of all reconstructed electrons and photons
  TH1F *relResEPSum_PA1 = new TH1F("rEP_PA1", filePrefix1 + ":Electron + photon sum", 40, -1, 1); //resolution of polar angle reconstructed when the reconstruction particle is made up of all reconstructed electrons and photons
  TH1F *relResEPSum_azimuth1 = new TH1F("rEP_a1", filePrefix1 + ":Electron + photon sum", 40, -1, 1); //resolution of azimuth reconstruction when the reconstructed particle is made up of all reconstructed electrons and photons
  TH1F *relResBestMatch_pt1 = new TH1F("rBM_pt1", filePrefix1 + ":Best match link", 40, -1, 1); //resolution of transverse momentum reconstruction when the reconstructed particle is the best match
  TH1F *relResBestMatch_PA1 = new TH1F("rMB_PA1", filePrefix1 + ":Best match link", 40, -1, 1); //resolution of polar angle reconstructed when the reconstruction particle is the best match
  TH1F *relResBestMatch_azimuth1 = new TH1F("rMB_a1", filePrefix1 + ":Best match link", 40, -2, 2); //resolution of azimuth reconstruction when the reconstructed particle is made up of the best match
  TH1F *relResPhotonReco_pt1 = new TH1F("rPR_pt1", filePrefix1 + ":Best match photon link", 40, -1, 1); //resolution of transverse momentum reconstruction when the reconstructed particle is the best photon match
  TH1F *relResPhotonReco_PA1 = new TH1F("rPR_PA1", filePrefix1 + ":Best match photon link", 40, -1, 1); //resolution of transverse momentum reconstruction when the reconstructed particle is the best photon match
  TH1F *relResPhotonReco_azimuth1 = new TH1F("rPR_a1", filePrefix1 + ":Best match photon link", 40, -1, 1); //resolution of transverse momentum reconstruction when the reconstructed particle is the best photon match
  TH1F *relResPionReco_pt1 = new TH1F("rPiR_pt1", filePrefix1 + " Pion p_{T}", 40, -1, 1); //resolution of transverse momentum reconstruction when the reconstructed particle is the best pion match
  TH1F *relResPionReco_PA1 = new TH1F("rPiR_PA1", filePrefix1 + " Pion p_{T}", 40, -1, 1); //resolution of transverse momentum reconstruction when the reconstructed particle is the best pion match
  TH1F *relResPionReco_azimuth1 = new TH1F("rPiR_a1", filePrefix1 + " Pion p_{T}", 40, -1, 1); //resolution of transverse momentum reconstruction when the reconstructed particle is the best pion match
  TH1F *relResTrackLink_pt1 = new TH1F("rTck_pt1", filePrefix1 + ":Reco track", 40, -1, 1); //transverse momentum from individual track reconstruction segments

  //histograms for storing attributes unique to certain reconstructed types of particles
  //transverse momentum
  TH1F *electronReco_pt1 = new TH1F("elRec_pt1", filePrefix1 + ":Electron pt", 30, 0, 3000);
  TH1F *photonReco_pt1 = new TH1F("phRec_pt1", filePrefix1 + ":Photon pt", 30, 0, 3000);
  TH1F *neutronReco_pt1 = new TH1F("nuRec_pt1", filePrefix1 + ":Neutron pt", 30, 0, 3000);
  TH1F *pionReco_pt1 = new TH1F("piRec_pt1", filePrefix1 + ":Pion pt", 30, 0, 3000);

  //polar angle
  TH1F *electronReco_PA1 = new TH1F("elRec_PA1", filePrefix1 + ":Electron PA", 40, 0, 3.2);
  TH1F *photonReco_PA1 = new TH1F("phRec_PA1", filePrefix1 + ":Photon PA", 40, 0, 3.2);
  TH1F *neutronReco_PA1 = new TH1F("nuRec_PA1", filePrefix1 + ":Neutron PA", 40, 0, 3.2);
  TH1F *pionReco_PA1 = new TH1F("piRec_PA1", filePrefix1 + ":Pion PA", 40, 0, 3.2);

  //azimuth
  TH1F *electronReco_azimuth1 = new TH1F("elRec_a1", filePrefix1 + ":Electron azimuth", 40, -3.2, 3.2);
  TH1F *photonReco_azimuth1 = new TH1F("phRec_a1", filePrefix1 + ":Photon azimuth", 40, -3.2, 3.2);
  TH1F *neutronReco_azimuth1 = new TH1F("nuRec_a1", filePrefix1 + ":Neutron azimuth", 40, -3.2, 3.2);
  TH1F *pionReco_azimuth1 = new TH1F("piRec_a1", filePrefix1 + ":Pion azimuth", 40, -3.2, 3.2);
  
  //FILE 2
  //histograms for attributes of electrons
  //transverse momentum
  TH1F *realPassed_pt2 = new TH1F("rP_pt2", filePrefix2 + ":Linked electrons", 40, 0, 2000); //transverse momentum of electrons that are linked to reconstructed electrons
  TH1F *realAllPassed_pt2 = new TH1F("AP_pt2", filePrefix2 + ":Best match", 40, 0, 2000); //best match transverse momentum from linked electrons and photons
  TH1F *realEPSum_pt2 = new TH1F("EPS_pt2", filePrefix2 + ":Summed particles", 40, 0, 2000); //summed transverse momentum from linked electrons and photons
  TH1F *recoTrackLink_pt2 = new TH1F("recTck_pt2", filePrefix2 + ":Reco track", 40, -2000, 2000); //transverse momentum from individual track reconstruction segments
  TH1F *realAll_pt2 = new TH1F("rA_pt2", filePrefix2 + ":All electrons", 40, 0, 2000); //transverse momentum of the truth particles
  //Polar angle
  TH1F *realPassed_PA2 = new TH1F("rP_PA2", filePrefix2 + ":Linked electrons", 40, 0, 3.2); //see above
  TH1F *realAllPassed_PA2 = new TH1F("AP_PA2", filePrefix2 + ":Best match", 40, 0, 3.2);
  TH1F *realEPSum_PA2 = new TH1F("EPS_PA2", filePrefix2 + ":Summed particles", 40, 0, 3.2);
  TH1F *realAll_PA2 = new TH1F("rA_PA2", filePrefix2 + ":All electrons", 40, 0, 3.2);
  //Azimuth
  TH1F *realPassed_azimuth2 = new TH1F("rP_a2", filePrefix2 + ":Linked electrons", 40, -3.2, 3.2); //see above
  TH1F *realAllPassed_azimuth2 = new TH1F("AP_a2", filePrefix2 + ":Best match", 40, -3.2, 3.2);
  TH1F *realEPSum_azimuth2 = new TH1F("EPS_a2", filePrefix2 + ":Summed particles", 40, -3.2, 3.2);
  TH1F *realAll_azimuth2 = new TH1F("rA_a2", filePrefix2 + ":All electrons", 40, -3.2, 3.2);

  //histograms for weights of particles
  TH1D *recoWeight_El2 = new TH1D("realElWeights2", filePrefix2 + ":Reco link electron weights", 26, 0, 1.3); //link weights of reconstructed electrons linked to truth electrons
  TH1D *recoWeight_BM2 = new TH1D("realAllWeights2", filePrefix2 + ":Reco link all weights", 26, 0, 1.3); // link weights of particles linked to truth electrons

  //histogam for the number of hits of best match particles
  TH1D *numHits_BM2 = new TH1D("nHitBM2", filePrefix2 + ":Number of BM hits", 30, 0, 30); //number of hits that the track encounters
    
  //histograms for storing differences between values
  TH1F *relResOG_pt2 = new TH1F("diffPt2", filePrefix2 + ":Simple electron linking", 40, -1, 1); //difference between reconstructed linked particles and the truth particles in transverse momentum
  TH1F *relResOG_PA2 = new TH1F("diffPA2", filePrefix2 + ":Simple electron linking", 40, -1, 1); //difference between reconstructed linked particles and the truth particles in transverse momentum
  TH1F *relResOG_azimuth2 = new TH1F("diffA2", filePrefix2 + ":Simple electron linking", 40, -1, 1); //difference between reconstructed linked particles and the truth particles in transverse momentum
  TH1F *deltaR_Hist2 = new TH1F("deltaR2", filePrefix2 + ":MyLCTuple", 25, 0, 0.0025); //delta R (delta polar angle and azimuth combined)
  TH1F *relResEPSum_pt2 = new TH1F("rEP_pt2", filePrefix2 + ":Electron + photon sum", 40, -1, 1); //resolution of transverse momentum reconstruction when the reconstructed particle is made up of all reconstructed electrons and photons
  TH1F *relResEPSum_PA2 = new TH1F("rEP_PA2", filePrefix2 + ":Electron + photon sum", 40, -1, 1); //resolution of polar angle reconstructed when the reconstruction particle is made up of all reconstructed electrons and photons
  TH1F *relResEPSum_azimuth2 = new TH1F("rEP_a2", filePrefix2 + ":Electron + photon sum", 40, -1, 1); //resolution of azimuth reconstruction when the reconstructed particle is made up of all reconstructed electrons and photons
  TH1F *relResBestMatch_pt2 = new TH1F("rBM_pt2", filePrefix2 + ":Best match link", 40, -1, 1); //resolution of transverse momentum reconstruction when the reconstructed particle is the best match
  TH1F *relResBestMatch_PA2 = new TH1F("rMB_PA2", filePrefix2 + ":Best match link", 40, -1, 1); //resolution of polar angle reconstructed when the reconstruction particle is the best match
  TH1F *relResBestMatch_azimuth2 = new TH1F("rMB_a2", filePrefix2 + ":Best match link", 40, -1, 1); //resolution of azimuth reconstruction when the reconstructed particle is made up of the best match
  TH1F *relResPhotonReco_pt2 = new TH1F("rPR_pt2", filePrefix2 + ":Best match photon link", 40, -1, 1); //resolution of transverse momentum reconstruction when the reconstructed particle is the best photon match
  TH1F *relResPhotonReco_PA2 = new TH1F("rPR_PA2", filePrefix2 + ":Best match photon link", 40, -1, 1); //resolution of transverse momentum reconstruction when the reconstructed particle is the best photon match
  TH1F *relResPhotonReco_azimuth2 = new TH1F("rPR_a2", filePrefix2 + ":Best match photon link", 40, -1, 1); //resolution of transverse momentum reconstruction when the reconstructed particle is the best photon match
  TH1F *relResPionReco_pt2 = new TH1F("rPiR_pt2", filePrefix2 + " Pion p_{T}", 40, -1, 1); //resolution of transverse momentum reconstruction when the reconstructed particle is the best pion match
  TH1F *relResPionReco_PA2 = new TH1F("rPiR_PA2", filePrefix2 + " Pion p_{T}", 40, -1, 1); //resolution of transverse momentum reconstruction when the reconstructed particle is the best pion match
  TH1F *relResPionReco_azimuth2 = new TH1F("rPiR_a2", filePrefix2 + " Pion p_{T}", 40, -1, 1); //resolution of transverse momentum reconstruction when the reconstructed particle is the best pion match
  TH1F *relResTrackLink_pt2 = new TH1F("rTck_pt2", filePrefix2 + ":Reco track", 40, -1, 1); //transverse momentum from individual track reconstruction segments

  //histograms for storing attributes unique to certain reconstructed types of particles
  //transverse momentum
  TH1F *electronReco_pt2 = new TH1F("elRec_pt2", filePrefix2 + ":Electron pt", 30, 0, 3000);
  TH1F *photonReco_pt2 = new TH1F("phRec_pt2", filePrefix2 + ":Photon pt", 30, 0, 3000);
  TH1F *neutronReco_pt2 = new TH1F("nuRec_pt2", filePrefix2 + ":Neutron pt", 30, 0, 3000);
  TH1F *pionReco_pt2 = new TH1F("piRec_pt2", filePrefix2 + ":Pion pt", 30, 0, 3000);

  //polar angle
  TH1F *electronReco_PA2 = new TH1F("elRec_PA2", filePrefix2 + ":Electron PA", 40, 0, 3.2);
  TH1F *photonReco_PA2 = new TH1F("phRec_PA2", filePrefix2 + ":Photon PA", 40, 0, 3.2);
  TH1F *neutronReco_PA2 = new TH1F("nuRec_PA2", filePrefix2 + ":Neutron PA", 40, 0, 3.2);
  TH1F *pionReco_PA2 = new TH1F("piRec_PA2", filePrefix2 + ":Pion PA", 40, 0, 3.2);

  //azimuth
  TH1F *electronReco_azimuth2 = new TH1F("elRec_a2", filePrefix2 + ":Electron azimuth", 40, -3.2, 3.2);
  TH1F *photonReco_azimuth2 = new TH1F("phRec_a2", filePrefix2 + ":Photon azimuth", 40, -3.2, 3.2);
  TH1F *neutronReco_azimuth2 = new TH1F("nuRec_a2", filePrefix2 + ":Neutron azimuth", 40, -3.2, 3.2);
  TH1F *pionReco_azimuth2 = new TH1F("piRec_a2", filePrefix2 + ":Pion azimuth", 40, -3.2, 3.2);

  //creates a reader that will traverse the events of the simulation
  TTreeReader myReader1("MyLCTuple", openFile1);
  TTreeReader myReader2("MyLCTuple", openFile2);

  //arrays that hold the values of branchesTTreeReader myReader1("MyLCTuple", openFile);
  TTreeReaderArray<int> r2f_RA1(myReader1, "r2f"); //link reco particle index number
  TTreeReaderArray<int> r2t_RA1(myReader1, "r2t"); //link truth partcie index number
  TTreeReaderArray<Float_t> r2w_RA1(myReader1, "r2w"); //link weight
  TTreeReaderArray<int> trthn_RA1(myReader1, "trthn"); //number of track hits
  TTreeReaderArray<Float_t> mcmox_RA1(myReader1, "mcmox"); //truth x momentum
  TTreeReaderArray<Float_t> mcmoy_RA1(myReader1, "mcmoy"); //truth y momentum
  TTreeReaderArray<Float_t> mcmoz_RA1(myReader1, "mcmoz"); //truth z momentum
  TTreeReaderArray<Float_t> rcmox_RA1(myReader1, "rcmox"); //reco x momentum
  TTreeReaderArray<Float_t> rcmoy_RA1(myReader1, "rcmoy"); //reco y momentum
  TTreeReaderArray<Float_t> rcmoz_RA1(myReader1, "rcmoz"); //reco z momentum
  TTreeReaderArray<Float_t> tsome_RA1(myReader1, "tsome"); //reco track curviture
  TTreeReaderArray<int> rctyp_RA1(myReader1, "rctyp"); //reconstructed type (11 = electron, 22 = photon, 2112 = neutron)
  TTreeReaderArray<int> mcpdg_RA1(myReader1, "mcpdg"); //truth type (see list above)
  TTreeReaderArray<int> mcgst_RA1(myReader1, "mcgst"); //is this particle a generating particle or not

  TTreeReaderArray<int> r2f_RA2(myReader2, "r2f"); //link reco particle index number
  TTreeReaderArray<int> r2t_RA2(myReader2, "r2t"); //link truth partcie index number
  TTreeReaderArray<Float_t> r2w_RA2(myReader2, "r2w"); //link weight
  TTreeReaderArray<int> trthn_RA2(myReader2, "trthn"); //number of track hits
  TTreeReaderArray<Float_t> mcmox_RA2(myReader2, "mcmox"); //truth x momentum
  TTreeReaderArray<Float_t> mcmoy_RA2(myReader2, "mcmoy"); //truth y momentum
  TTreeReaderArray<Float_t> mcmoz_RA2(myReader2, "mcmoz"); //truth z momentum
  TTreeReaderArray<Float_t> rcmox_RA2(myReader2, "rcmox"); //reco x momentum
  TTreeReaderArray<Float_t> rcmoy_RA2(myReader2, "rcmoy"); //reco y momentum
  TTreeReaderArray<Float_t> rcmoz_RA2(myReader2, "rcmoz"); //reco z momentum
  TTreeReaderArray<Float_t> tsome_RA2(myReader2, "tsome"); //reco track curviture
  TTreeReaderArray<int> rctyp_RA2(myReader2, "rctyp"); //reconstructed type (22 = electron, 22 = photon, 2112 = neutron)
  TTreeReaderArray<int> mcpdg_RA2(myReader2, "mcpdg"); //truth type (see list above)
  TTreeReaderArray<int> mcgst_RA2(myReader2, "mcgst"); //is this particle a generating particle or not

  //temporary variables that will reduce the number of get operatons
  int r2fTemp1, r2tTemp1, trthnTemp1;
  bool electronMatch1 = false;
  bool photonMatch1 = false;
  bool pionMatch1 = false;
  Float_t r2wTemp1, dPATemp1, dATemp1;
  TVector3 mcmoTemp1, rcmoTemp1, rcmoSumTemp1, BMTemp1, BMPTemp1;
  Float_t r2wMax1 = 0;
  int r2wMax_Index1 = -1;
  int count1 = 0;

  int r2fTemp2, r2tTemp2, trthnTemp2;
  bool electronMatch2 = false;
  bool photonMatch2 = false;
  bool pionMatch2 = false;
  Float_t r2wTemp2, dPATemp2, dATemp2;
  TVector3 mcmoTemp2, rcmoTemp2, rcmoSumTemp2, BMTemp2, BMPTemp2;
  Float_t r2wMax2 = 0;
  int r2wMax_Index2 = -1;
  int count2 = 0;

  /////////////////////
  //HISTOGRAM FILLING//
  /////////////////////
  
  //loops over each event
  while(myReader1.Next()) {
    myReader2.Next();

    if(mcmox_RA1.At(0) != mcmox_RA2.At(0)) {
      break;
    }

    //loops over all particles to find the generating particle
    for(int i = 0; i < mcmox_RA1.GetSize(); i++) {
      if(mcgst_RA1.At(i)) {

	mcmoTemp1.SetXYZ(mcmox_RA1.At(i), mcmoy_RA1.At(i), mcmoz_RA1.At(i));

	//adds generating particle to histogram
	realAll_pt1->Fill(mcmoTemp1.Perp());
	realAll_PA1->Fill(mcmoTemp1.Theta());
	realAll_azimuth1->Fill(mcmoTemp1.Phi());
	break; //this might cause issues in the future, double check to make sure that there are not two or more generating particles
      }
    }

    //loops over all particles to find the generating particle
    for(int i = 0; i < mcmox_RA2.GetSize(); i++) {
      if(mcgst_RA2.At(i)) {

	mcmoTemp2.SetXYZ(mcmox_RA2.At(i), mcmoy_RA2.At(i), mcmoz_RA2.At(i));

	//adds generating particle to histogram
	realAll_pt2->Fill(mcmoTemp2.Perp());
	realAll_PA2->Fill(mcmoTemp2.Theta());
	realAll_azimuth2->Fill(mcmoTemp2.Phi());
	break; //this might cause issues in the future, double check to make sure that there are not two or more generating particles
      }
    }

    //resets data sum variables
    rcmoTemp1.SetXYZ(0, 0, 0);
    BMTemp1.SetXYZ(0, 0, 0);
    rcmoSumTemp1.SetXYZ(0, 0, 0);

    //loops over all reconstructed particles
    for(int i = 0; i < rcmox_RA1.GetSize(); i++) {
      rcmoTemp1.SetXYZ(rcmox_RA1.At(i), rcmoy_RA1.At(i), rcmoz_RA1.At(i));

      //if that particle is an electron or photon, add it the sum variables (this accounts for particle mis-identification)
      if(abs(rctyp_RA1.At(i)) == 11 || abs(rctyp_RA1.At(i)) == 22 || abs(rctyp_RA1.At(i)) == 211) {
	for(int j = 0; j < r2f_RA1.GetSize(); j++) {
	  if(r2f_RA1.At(j) == i && mcgst_RA1.At(r2t_RA1.At(j))) {
	    rcmoSumTemp1.SetXYZ(rcmoSumTemp1.X() + rcmoTemp1.X(), rcmoSumTemp1.Y() + rcmoTemp1.Y(), rcmoSumTemp1.Z() + rcmoTemp1.Z());
	    break;
	  }
	}
      }

      //if the particle is a elect ron, add it to one set of histograms
      if(abs(rctyp_RA1.At(i)) == 11) {
	electronReco_pt1->Fill(rcmoTemp1.Perp());
	electronReco_PA1->Fill(rcmoTemp1.Theta());
	electronReco_azimuth1->Fill(rcmoTemp1.Phi());

	//if the particle is a photon, add it to another set of histograms
      } else if(abs(rctyp_RA1.At(i)) == 22) {
	photonReco_pt1->Fill(rcmoTemp1.Perp());
	photonReco_PA1->Fill(rcmoTemp1.Theta());
	photonReco_azimuth1->Fill(rcmoTemp1.Phi());

	//if the particle is a neutron, add it to a third set of histograms
      } else if(abs(rctyp_RA1.At(i)) == 2112) {
	neutronReco_pt1->Fill(rcmoTemp1.Perp());
	neutronReco_PA1->Fill(rcmoTemp1.Theta());
	neutronReco_azimuth1->Fill(rcmoTemp1.Phi());

	//if the particle is a pion, add it to a fourth set of histograms
      } else if(abs(rctyp_RA1.At(i)) == 211) {
	pionReco_pt1->Fill(rcmoTemp1.Perp());
	pionReco_PA1->Fill(rcmoTemp1.Theta());
	pionReco_azimuth1->Fill(rcmoTemp1.Phi());
      }
    }

    //resets data sum variables
    rcmoTemp2.SetXYZ(0, 0, 0);
    BMTemp2.SetXYZ(0, 0, 0);
    rcmoSumTemp2.SetXYZ(0, 0, 0);

    //loops over all reconstructed particles
    for(int i = 0; i < rcmox_RA2.GetSize(); i++) {
      rcmoTemp2.SetXYZ(rcmox_RA2.At(i), rcmoy_RA2.At(i), rcmoz_RA2.At(i));

      //if that particle is an electron or photon, add it the sum variables (this accounts for particle mis-identification)
      if(abs(rctyp_RA2.At(i)) == 11 || abs(rctyp_RA2.At(i)) == 22 || abs(rctyp_RA2.At(i)) == 211) {
	for(int j = 0; j < r2f_RA2.GetSize(); j++) {
	  if(r2f_RA2.At(j) == i && mcgst_RA2.At(r2t_RA2.At(j))) {
	    rcmoSumTemp2.SetXYZ(rcmoSumTemp2.X() + rcmoTemp2.X(), rcmoSumTemp2.Y() + rcmoTemp2.Y(), rcmoSumTemp2.Z()+ rcmoTemp2.Z());
	    break;
	  }
	}
      }

      //if the particle is a elect ron, add it to one set of histograms
      if(abs(rctyp_RA2.At(i)) == 11) {
	electronReco_pt2->Fill(rcmoTemp2.Perp());
	electronReco_PA2->Fill(rcmoTemp2.Theta());
	electronReco_azimuth2->Fill(rcmoTemp2.Phi());

	//if the particle is a photon, add it to another set of histograms
      } else if(abs(rctyp_RA2.At(i)) == 22) {
	photonReco_pt2->Fill(rcmoTemp2.Perp());
	photonReco_PA2->Fill(rcmoTemp2.Theta());
	photonReco_azimuth2->Fill(rcmoTemp2.Phi());

	//if the particle is a neutron, add it to a third set of histograms
      } else if(abs(rctyp_RA2.At(i)) == 2112) {
	neutronReco_pt2->Fill(rcmoTemp2.Perp());
	neutronReco_PA2->Fill(rcmoTemp2.Theta());
	neutronReco_azimuth2->Fill(rcmoTemp2.Phi());

	//if the particle is a pion, add it to a fourth set of histograms
      } else if(abs(rctyp_RA2.At(i)) == 211) {
	pionReco_pt2->Fill(rcmoTemp2.Perp());
	pionReco_PA2->Fill(rcmoTemp2.Theta());
	pionReco_azimuth2->Fill(rcmoTemp2.Phi());
      }
    }

    //searches for the largest weight of the relationship that links a reconstructed particle to a particle that we know all the details of
    
    for(int i = 0; i < r2f_RA1.GetSize(); i++) {
      r2tTemp1 = r2t_RA1.At(i);
      r2fTemp1 = r2f_RA1.At(i);
      r2wTemp1 = r2w_RA1.At(i);

      //checks to see if those particles are electrons and if the paticle is an originial generating particle
      if(abs(mcpdg_RA1.At(r2tTemp1)) == 11 && abs(rctyp_RA1.At(r2fTemp1)) == 11 && mcgst_RA1.At(r2tTemp1)) {
	if(r2wTemp1 > r2wMax1) {
	  electronMatch1 = true;
	  r2wMax1 = r2wTemp1;
	  r2wMax_Index1 = i;
	}
      }
    }

    for(int i = 0; i < r2f_RA2.GetSize(); i++) {
      r2tTemp2 = r2t_RA2.At(i);
      r2fTemp2 = r2f_RA2.At(i);
      r2wTemp2 = r2w_RA2.At(i);

      //checks to see if those particles are electrons and if the paticle is an originial generating particle
      if(abs(mcpdg_RA2.At(r2tTemp2)) == 11 && abs(rctyp_RA2.At(r2fTemp2)) == 11 && mcgst_RA2.At(r2tTemp2)) {
	if(r2wTemp2 > r2wMax2) {
	  electronMatch2 = true;
	  r2wMax2 = r2wTemp2;
	  r2wMax_Index2 = i;
	}
      }
    }
    
    //if a relationship was found between those particles, it is added to histograms
    if(electronMatch1 && electronMatch2) {

      //adds the weight to a histogram
      recoWeight_El1->Fill(r2wMax1);

      //find the particles that the particle linker refers to
      r2tTemp1 = r2t_RA1.At(r2wMax_Index1);
      r2fTemp1 = r2f_RA1.At(r2wMax_Index1);
      r2wTemp1 = r2w_RA1.At(r2wMax_Index1);

      //finds the reco momenta of the particle
      rcmoTemp1.SetXYZ(rcmox_RA1.At(r2fTemp1), rcmoy_RA1.At(r2fTemp1), rcmoz_RA1.At(r2fTemp1));
      BMTemp1.SetXYZ(rcmox_RA1.At(r2fTemp1), rcmoy_RA1.At(r2fTemp1), rcmoz_RA1.At(r2fTemp1));
      
      //finds the differences between truth and reconstructed attributes
      dPATemp1 = mcmoTemp1.Theta() - rcmoTemp1.Theta();
      dATemp1 = mcmoTemp1.Phi() - rcmoTemp1.Phi();

      //fills in histograms with attributes
      realPassed_pt1->Fill(mcmoTemp1.Perp());
      realPassed_PA1->Fill(mcmoTemp1.Theta());
      realPassed_azimuth1->Fill(mcmoTemp1.Phi());
      
      //fills in histograms with differences
      relResOG_pt1->Fill((mcmoTemp1.Perp() - rcmoTemp1.Perp()) / mcmoTemp1.Perp());
      relResOG_PA1->Fill((mcmoTemp1.Theta() - rcmoTemp1.Theta()) / mcmoTemp1.Theta());
      relResOG_azimuth1->Fill((mcmoTemp1.Phi() - rcmoTemp1.Phi()) / mcmoTemp1.Phi());
      
      deltaR_Hist1->Fill(sqrt((dPATemp1*dPATemp1) + (dATemp1*dATemp1)));
    }

    if(electronMatch2) {
      //adds the weight to a histogram
      recoWeight_El2->Fill(r2wMax2);
      
      //find the particles that the particle linker refers to
      r2tTemp2 = r2t_RA2.At(r2wMax_Index2);
      r2fTemp2 = r2f_RA2.At(r2wMax_Index2);
      r2wTemp2 = r2w_RA2.At(r2wMax_Index2);

      //finds the reco momenta of the particle
      rcmoTemp2.SetXYZ(rcmox_RA2.At(r2fTemp2), rcmoy_RA2.At(r2fTemp2), rcmoz_RA2.At(r2fTemp2));
      BMTemp2.SetXYZ(rcmox_RA2.At(r2fTemp2), rcmoy_RA2.At(r2fTemp2), rcmoz_RA2.At(r2fTemp2));

      //finds the differences between truth and reconstructed attributes
      dPATemp2 = mcmoTemp2.Theta() - rcmoTemp2.Theta();
      dATemp2 = mcmoTemp2.Phi() - rcmoTemp2.Phi();

      //fills in histograms with attributes
      realPassed_pt2->Fill(mcmoTemp2.Perp());
      realPassed_PA2->Fill(mcmoTemp2.Theta());
      realPassed_azimuth2->Fill(mcmoTemp2.Phi());

      //fills in histograms with differences
      relResOG_pt2->Fill((mcmoTemp2.Perp() - rcmoTemp2.Perp()) / mcmoTemp2.Perp());
      relResOG_PA2->Fill((mcmoTemp2.Theta() - rcmoTemp2.Theta()) / mcmoTemp2.Theta());
      relResOG_azimuth2->Fill((mcmoTemp2.Phi() - rcmoTemp2.Phi()) / mcmoTemp2.Phi());

      deltaR_Hist2->Fill(sqrt((dPATemp2*dPATemp2) + (dATemp2*dATemp2)));
    }

    //if no link was found between a reconstructed and truth electron, I expand the search to reconstructed photons as well
    if(r2wMax_Index1 == -1) {
      for(int i = 0; i < r2f_RA1.GetSize(); i++) {
	//find the two particles in which the link occurs between
	r2tTemp1 = r2t_RA1.At(i);
	r2fTemp1 = r2f_RA1.At(i);
	r2wTemp1 = r2w_RA1.At(i);

	//checks to see if the truth particle is an electron, if the reconstructed particle is a photon, and if the truth particle is an originial generating particle
	if(abs(mcpdg_RA1.At(r2tTemp1)) == 11 && (abs(rctyp_RA1.At(r2fTemp1)) == 22 || abs(rctyp_RA1.At(r2fTemp1)) == 211) && mcgst_RA1.At(r2tTemp1)) {
	  if(r2wTemp1 > r2wMax1) {
	    r2wMax1 = r2wTemp1;
	    r2wMax_Index1 = i;
	    photonMatch1 = true;
	  }
	}
      }
      if(r2wMax_Index1 != -1) {
	if(rctyp_RA1.At(r2wMax_Index1) == 22) {
	  photonMatch1 = true;
	} else if(rctyp_RA1.At(r2wMax_Index1) == 211) {
	  pionMatch1 = true;
	}
      }
    }

    //if no link was found between a reconstructed and truth electron, I expand the search to reconstructed photons as well
    if(r2wMax_Index2 == -1) {
      for(int i = 0; i < r2f_RA2.GetSize(); i++) {
	//find the two particles in which the link occurs between
	r2tTemp2 = r2t_RA2.At(i);
	r2fTemp2 = r2f_RA2.At(i);
	r2wTemp2 = r2w_RA2.At(i);

	//checks to see if the truth particle is an electron, if the reconstructed particle is a photon, and if the truth particle is an originial generating particle
	if(abs(mcpdg_RA2.At(r2tTemp2)) == 11 && (abs(rctyp_RA2.At(r2fTemp2)) == 22 || abs(rctyp_RA1.At(r2fTemp1)) == 211) && mcgst_RA2.At(r2tTemp2)) {
	  if(r2wTemp2 > r2wMax2) {
	    r2wMax2 = r2wTemp2;
	    r2wMax_Index2 = i;
	    photonMatch2 = true;
	  }
	}
      }
      if(r2wMax_Index2 != -1) {
	if(rctyp_RA2.At(r2wMax_Index2) == 22) {
	  photonMatch2 = true;
	} else if(rctyp_RA2.At(r2wMax_Index2) == 211) {
	  pionMatch2 = true;
	}
      }
    }

    //repeats the previous steps that adds data to histograms but instead adds it to the histograms that holds links between photons and electrons as well
    if(photonMatch1 || electronMatch1 || pionMatch1) {
      r2tTemp1 = r2t_RA1.At(r2wMax_Index1);
      r2fTemp1 = r2f_RA1.At(r2wMax_Index1);
      r2wTemp1 = r2w_RA1.At(r2wMax_Index1);

      //adds the weight of the link to a histogram
      recoWeight_BM1->Fill(r2wMax1);

      //adds the number of hits of the link to a histogram
      if(r2wMax_Index1 == 0) {
	numHits_BM1->Fill(trthn_RA1.At(0));
      }

      rcmoTemp1.SetXYZ(rcmox_RA1.At(r2fTemp1), rcmoy_RA1.At(r2fTemp1), rcmoz_RA1.At(r2fTemp1));

      realAllPassed_pt1->Fill(mcmoTemp1.Perp());
      realAllPassed_PA1->Fill(mcmoTemp1.Theta());
      realAllPassed_azimuth1->Fill(mcmoTemp1.Phi());

      //finds the difference between the truth and best match reconstructed particles
      relResBestMatch_pt1->Fill((mcmoTemp1.Perp() - rcmoTemp1.Perp()) / mcmoTemp1.Perp());
      relResBestMatch_PA1->Fill((mcmoTemp1.Theta() - rcmoTemp1.Theta()) / mcmoTemp1.Theta());
      relResBestMatch_azimuth1->Fill((mcmoTemp1.Phi() - rcmoTemp1.Phi()) / mcmoTemp1.Phi());

      if(photonMatch1) {
	relResPhotonReco_pt1->Fill((mcmoTemp1.Perp() - rcmoTemp1.Perp()) / mcmoTemp1.Perp());
	relResPhotonReco_PA1->Fill((mcmoTemp1.Theta() - rcmoTemp1.Theta()) / mcmoTemp1.Theta());
	relResPhotonReco_azimuth1->Fill((mcmoTemp1.Phi() - rcmoTemp1.Phi()) / mcmoTemp1.Phi());

	if(tsome_RA1.GetSize() > 0){
	  recoTrackLink_pt1->Fill(((0.3 * 3.57) / tsome_RA1.At(0)) / 1000);
	  relResTrackLink_pt1->Fill((mcmoTemp1.Perp() - abs(((0.3 * 3.57) / tsome_RA1.At(0)) / 1000)) / mcmoTemp1.Perp());
	}
      }
      if(pionMatch1) {
	relResPionReco_pt1->Fill((mcmoTemp1.Perp() - rcmoTemp1.Perp()) / mcmoTemp1.Perp());
	relResPionReco_PA1->Fill((mcmoTemp1.Theta() - rcmoTemp1.Theta()) / mcmoTemp1.Theta());
	relResPionReco_azimuth1->Fill((mcmoTemp1.Phi() - rcmoTemp1.Phi()) / mcmoTemp1.Phi());
      }
    }

    if(photonMatch2 || electronMatch2 || pionMatch2) {
      r2tTemp2 = r2t_RA2.At(r2wMax_Index2);
      r2fTemp2 = r2f_RA2.At(r2wMax_Index2);
      r2wTemp2 = r2w_RA2.At(r2wMax_Index2);

      //adds the weight of the link to a histogram
      recoWeight_BM2->Fill(r2wMax2);

      //adds the number of hits of the link to a histogram
      if(r2wMax_Index2 == 0) {
	numHits_BM2->Fill(trthn_RA2.At(0));
      }

      rcmoTemp2.SetXYZ(rcmox_RA2.At(r2fTemp2), rcmoy_RA2.At(r2fTemp2), rcmoz_RA2.At(r2fTemp2));

      realAllPassed_pt2->Fill(mcmoTemp2.Perp());
      realAllPassed_PA2->Fill(mcmoTemp2.Theta());
      realAllPassed_azimuth2->Fill(mcmoTemp2.Phi());

      //finds the difference between the truth and best match reconstructed particles
      relResBestMatch_pt2->Fill((mcmoTemp2.Perp() - rcmoTemp2.Perp()) / mcmoTemp2.Perp());
      relResBestMatch_PA2->Fill((mcmoTemp2.Theta() - rcmoTemp2.Theta()) / mcmoTemp2.Theta());
      relResBestMatch_azimuth2->Fill((mcmoTemp2.Phi() - rcmoTemp2.Phi()) / mcmoTemp2.Phi());

      if(photonMatch2) {
	relResPhotonReco_pt2->Fill((mcmoTemp2.Perp() - rcmoTemp2.Perp()) / mcmoTemp2.Perp());
	relResPhotonReco_PA2->Fill((mcmoTemp2.Theta() - rcmoTemp2.Theta()) / mcmoTemp2.Theta());
	relResPhotonReco_azimuth2->Fill((mcmoTemp2.Phi() - rcmoTemp2.Phi()) / mcmoTemp2.Phi());

	if(tsome_RA2.GetSize() > 0){
	  recoTrackLink_pt2->Fill(((0.3 * 3.57) / tsome_RA2.At(0)) / 1000);
	  relResTrackLink_pt2->Fill((mcmoTemp2.Perp() - abs(((0.3 * 3.57) / tsome_RA2.At(0)) / 1000)) / mcmoTemp2.Perp());
	}
      }
      if(pionMatch2) {
	relResPionReco_pt2->Fill((mcmoTemp2.Perp() - rcmoTemp2.Perp()) / mcmoTemp2.Perp());
	relResPionReco_PA2->Fill((mcmoTemp2.Theta() - rcmoTemp2.Theta()) / mcmoTemp2.Theta());
	relResPionReco_azimuth2->Fill((mcmoTemp2.Phi() - rcmoTemp2.Phi()) / mcmoTemp2.Phi());
      }
    }
    
    count1++;
    r2wMax1 = 0;
    r2wMax_Index1 = -1;
    electronMatch1 = false;
    photonMatch1 = false;
    pionMatch1 = false;

    count2++;
    r2wMax2 = 0;
    r2wMax_Index2 = -1;
    electronMatch2 = false;
    photonMatch2 = false;
    pionMatch2 = false;

    //finds the difference between the truth and reconstructed particles (sum of photons and electrons)
    if(rcmoSumTemp1.X() != 0 || rcmoSumTemp1.Y() != 0 || rcmoSumTemp1.Z() != 0){

      //fills in the electron and photon sum histograms
      realEPSum_pt1->Fill(rcmoSumTemp1.Perp());
      realEPSum_PA1->Fill(rcmoSumTemp1.Theta());
      realEPSum_azimuth1->Fill(rcmoSumTemp1.Phi());
      
      //fills in the electron and photon sum resolution histograms
      relResEPSum_pt1->Fill((mcmoTemp1.Perp() - rcmoSumTemp1.Perp()) / mcmoTemp1.Perp());
      relResEPSum_PA1->Fill((mcmoTemp1.Theta() - rcmoSumTemp1.Theta()) / mcmoTemp1.Theta());
      relResEPSum_azimuth1->Fill((mcmoTemp1.Phi() - rcmoSumTemp1.Phi()) / mcmoTemp1.Phi());
    }

    if(rcmoSumTemp2.X() != 0 || rcmoSumTemp2.Y() != 0 || rcmoSumTemp2.Z() != 0) {
      
      //fills in the electron and photon sum histograms
      realEPSum_pt2->Fill(rcmoSumTemp2.Perp());
      realEPSum_PA2->Fill(rcmoSumTemp2.Theta());
      realEPSum_azimuth2->Fill(rcmoSumTemp2.Phi());

      //fills in the electron and photon sum resolution histograms
      relResEPSum_pt2->Fill((mcmoTemp2.Perp() - rcmoSumTemp2.Perp()) / mcmoTemp2.Perp());
      relResEPSum_PA2->Fill((mcmoTemp2.Theta() - rcmoSumTemp2.Theta()) / mcmoTemp2.Theta());
      relResEPSum_azimuth2->Fill((mcmoTemp2.Phi() - rcmoSumTemp2.Phi()) / mcmoTemp2.Phi());
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
  THStack *hs = new THStack("hs", "Transverse momenta of electrons, photons, neutrons, and pions;Transverse momentum (GeV);Count");

  c->SetLogy();
  
  electronReco_pt1->SetLineColor(kBlue);
  photonReco_pt1->SetLineColor(kBlack);
  neutronReco_pt1->SetLineColor(kRed);
  pionReco_pt1->SetLineColor(kGreen);

  electronReco_pt2->SetLineColor(kBlue);
  photonReco_pt2->SetLineColor(kBlack);
  neutronReco_pt2->SetLineColor(kRed);
  pionReco_pt2->SetLineColor(kGreen);

  electronReco_pt2->SetLineStyle(6);
  photonReco_pt2->SetLineStyle(6);
  neutronReco_pt2->SetLineStyle(6);
  pionReco_pt2->SetLineStyle(6);

  hs->Add(electronReco_pt1);
  hs->Add(photonReco_pt1);
  hs->Add(neutronReco_pt1);
  hs->Add(pionReco_pt1);

  hs->Add(electronReco_pt2);
  hs->Add(photonReco_pt2);
  hs->Add(neutronReco_pt2);
  hs->Add(pionReco_pt2);
  
  hs->Draw("nostack");

  gPad->SetGrid(1, 0);
  gPad->BuildLegend(0.6, 0.6, 0.9, 0.9, "");
  
  c->SaveAs(saveDir + "/reconstructed-types/allPart_pt.png");
  c->Close();
  c = new TCanvas();

  //draws the pt histograms for electrons on top of each other
  hs = new THStack("hs", "Transverse momenta of electrons;Transverse momentum (GeV);Count");

  c->SetLogy();

  hs->Add(electronReco_pt1);
  hs->Add(electronReco_pt2);

  hs->Draw("nostack");

  gPad->SetGrid(1, 0);
  gPad->BuildLegend(0.6, 0.6, 0.9, 0.9, "");

  c->SaveAs(saveDir + "/reconstructed-types/recoElectron_pt.png");
  c->Close();
  c = new TCanvas();

  //draws the pt histograms for photons on top of each other
  hs = new THStack("hs", "Transverse momenta of photons;Transverse momentum (GeV);Count");

  c->SetLogy();

  hs->Add(photonReco_pt1);
  hs->Add(photonReco_pt2);

  hs->Draw("nostack");

  gPad->SetGrid(1, 0);
  gPad->BuildLegend(0.6, 0.6, 0.9, 0.9, "");

  c->SaveAs(saveDir + "/reconstructed-types/recoPhoton_pt.png");
  c->Close();
  c = new TCanvas();

  //draws the pt histograms for neutrons on top of each other
  hs = new THStack("hs", "Transverse momenta of neutrons;Transverse momentum (GeV);Count");

  c->SetLogy();

  hs->Add(neutronReco_pt1);
  hs->Add(neutronReco_pt2);

  hs->Draw("nostack");

  gPad->SetGrid(1, 0);
  gPad->BuildLegend(0.6, 0.6, 0.9, 0.9, "");

  c->SaveAs(saveDir + "/reconstructed-types/recoNeutron_pt.png");
  c->Close();
  c = new TCanvas();

  //draws the pt histograms for pions on top of each other
  hs = new THStack("hs", "Transverse momenta of pions;Transverse momentum (GeV);Count");

  hs->Add(pionReco_pt1);
  hs->Add(pionReco_pt2);

  hs->Draw("nostack");

  gPad->SetGrid(1, 0);
  gPad->BuildLegend(0.6, 0.6, 0.9, 0.9, "");

  c->SaveAs(saveDir + "/reconstructed-types/recoPion_pt.png");
  c->Close();
  c = new TCanvas();

  //draws the polar angle histograms for electrons, photons, and neutrons on top of each other
  hs = new THStack("hs", "Polar angle of electrons, photons, neutrons, and pions;Polar angle (Rads);Count");

  c->SetLogy();

  electronReco_PA1->SetLineColor(kBlue);
  photonReco_PA1->SetLineColor(kBlack);
  neutronReco_PA1->SetLineColor(kRed);
  pionReco_PA1->SetLineColor(kGreen);

  electronReco_PA2->SetLineColor(kBlue);
  photonReco_PA2->SetLineColor(kBlack);
  neutronReco_PA2->SetLineColor(kRed);
  pionReco_PA2->SetLineColor(kGreen);

  electronReco_PA2->SetLineStyle(6);
  photonReco_PA2->SetLineStyle(6);
  neutronReco_PA2->SetLineStyle(6);
  pionReco_PA2->SetLineStyle(6);
  
  hs->Add(electronReco_PA1);
  hs->Add(photonReco_PA1);
  hs->Add(neutronReco_PA1);
  hs->Add(pionReco_PA1);

  hs->Add(electronReco_PA2);
  hs->Add(photonReco_PA2);
  hs->Add(neutronReco_PA2);
  hs->Add(pionReco_PA2);
  
  hs->Draw("nostack");

  gPad->SetGrid(1, 0);
  gPad->BuildLegend(0.15, 0.1, 0.4, 0.4, "");
  
  c->SaveAs(saveDir + "/reconstructed-types/allPart_PA.png");
  c->Close();
  c = new TCanvas();

  //draws the polar angle histograms for electrons on top of each other
  hs = new THStack("hs", "Polar angle of electrons;Polar angle (Rads);Count");

  c->SetLogy();

  hs->Add(electronReco_PA1);
  hs->Add(electronReco_PA2);

  hs->Draw("nostack");

  gPad->SetGrid(1, 0);
  gPad->BuildLegend(0.4, 0.1, 0.6, 0.3, "");

  c->SaveAs(saveDir + "/reconstructed-types/recoElectron_PA.png");
  c->Close();
  c = new TCanvas();

  //draws the polar angle histograms for photons on top of each other
  hs = new THStack("hs", "Polar angle of photons;Polar angle (Rads);Count");

  c->SetLogy();

  hs->Add(photonReco_PA1);
  hs->Add(photonReco_PA2);

  hs->Draw("nostack");

  gPad->SetGrid(1, 0);
  gPad->BuildLegend(0.4, 0.1, 0.6, 0.3, "");

  c->SaveAs(saveDir + "/reconstructed-types/recoPhoton_PA.png");
  c->Close();
  c = new TCanvas();

  //draws the polar angle histograms for neutrons on top of each other
  hs = new THStack("hs", "Polar angle of neutrons;Polar angle (Rads);Count");

  c->SetLogy();

  hs->Add(neutronReco_PA1);
  hs->Add(neutronReco_PA2);

  hs->Draw("nostack");

  gPad->SetGrid(1, 0);
  gPad->BuildLegend(0.4, 0.1, 0.6, 0.3, "");

  c->SaveAs(saveDir + "/reconstructed-types/recoNeutron_PA.png");
  c->Close();
  c = new TCanvas();

  //draws the polar angle histograms for pions on top of each other
  hs = new THStack("hs", "Polar angle of pions;Polar angle (Rads);Count");

  hs->Add(pionReco_PA1);
  hs->Add(pionReco_PA2);

  hs->Draw("nostack");

  gPad->SetGrid(1, 0);
  gPad->BuildLegend(0.4, 0.7, 0.6, 0.9, "");

  c->SaveAs(saveDir + "/reconstructed-types/recoPion_PA.png");
  c->Close();
  c = new TCanvas();

  //draws the azimuth histograms for electrons, photons, and neutrons on top of each other
  hs = new THStack("hs", "Azimuth of electrons, photons, neutrons, and pions;Azimuth (Rads);Count");

  c->SetLogy();

  electronReco_azimuth1->SetLineColor(kBlue);
  photonReco_azimuth1->SetLineColor(kBlack);
  neutronReco_azimuth1->SetLineColor(kRed);
  pionReco_azimuth1->SetLineColor(kGreen);

  electronReco_azimuth2->SetLineColor(kBlue);
  photonReco_azimuth2->SetLineColor(kBlack);
  neutronReco_azimuth2->SetLineColor(kRed);
  pionReco_azimuth2->SetLineColor(kGreen);

  electronReco_azimuth2->SetLineStyle(6);
  photonReco_azimuth2->SetLineStyle(6);
  neutronReco_azimuth2->SetLineStyle(6);
  pionReco_azimuth2->SetLineStyle(6);
  
  hs->Add(electronReco_azimuth1);
  hs->Add(photonReco_azimuth1);
  hs->Add(neutronReco_azimuth1);
  hs->Add(pionReco_azimuth1);

  hs->Add(electronReco_azimuth2);
  hs->Add(photonReco_azimuth2);
  hs->Add(neutronReco_azimuth2);
  hs->Add(pionReco_azimuth2);

  hs->Draw("nostack");

  gPad->SetGrid(1, 0);
  gPad->BuildLegend(0.35, 0.1, 0.65, 0.4, "");

  c->SaveAs(saveDir + "/reconstructed-types/allPart_azimuth.png");
  c->Close();
  c = new TCanvas();

  //draws the azimuth histograms for electrons on top of each other
  hs = new THStack("hs", "Azimuth of electrons;Azimuth (Rads);Count");

  c->SetLogy();

  hs->Add(electronReco_azimuth1);
  hs->Add(electronReco_azimuth2);

  hs->Draw("nostack");

  gPad->SetGrid(1, 0);
  gPad->BuildLegend(0.4, 0.1, 0.6, 0.3, "");

  c->SaveAs(saveDir + "/reconstructed-types/recoElectron_azimuth.png");
  c->Close();
  c = new TCanvas();

  //draws the azimuth histograms for photons on top of each other
  hs = new THStack("hs", "Azimuth of photons;Azimuth (Rads);Count");

  c->SetLogy();

  hs->Add(photonReco_azimuth1);
  hs->Add(photonReco_azimuth2);

  hs->Draw("nostack");

  gPad->SetGrid(1, 0);
  gPad->BuildLegend(0.4, 0.1, 0.6, 0.3, "");

  c->SaveAs(saveDir + "/reconstructed-types/recoPhoton_azimuth.png");
  c->Close();
  c = new TCanvas();

  //draws the azimuth histograms for neutrons on top of each other
  hs = new THStack("hs", "Azimuth of neutrons;Azimuth (Rads);Count");

  c->SetLogy();

  hs->Add(neutronReco_azimuth1);
  hs->Add(neutronReco_azimuth2);

  hs->Draw("nostack");

  gPad->SetGrid(1, 0);
  gPad->BuildLegend(0.4, 0.1, 0.6, 0.3, "");

  c->SaveAs(saveDir + "/reconstructed-types/recoNeutron_azimuth.png");
  c->Close();
  c = new TCanvas();

  //draws the azimuth histograms for pions on top of each other
  hs = new THStack("hs", "Azimuth of pions;Azimuth (Rads);Count");

  hs->Add(pionReco_azimuth1);
  hs->Add(pionReco_azimuth2);

  hs->Draw("nostack");

  gPad->SetGrid(1, 0);
  gPad->BuildLegend(0.4, 0.1, 0.6, 0.3, "");

  c->SaveAs(saveDir + "/reconstructed-types/recoPion_azimuth.png");
  c->Close();
  c = new TCanvas();

  //==============================//
  //most likely electron reco link//
  //==============================//
  
  //draws the two pt histograms on each other
  hs = new THStack("hs", "pt of linked generating electrons and all generating electrons;Transverse momentum (GeV);Count");

  realAll_pt1->SetLineColor(kBlue);
  realPassed_pt1->SetLineColor(kBlack);

  realAll_pt2->SetLineColor(kBlue);
  realPassed_pt2->SetLineColor(kBlack);

  realAll_pt2->SetLineStyle(6);
  realPassed_pt2->SetLineStyle(6);

  hs->Add(realAll_pt1);
  hs->Add(realPassed_pt1);

  hs->Add(realAll_pt2);
  hs->Add(realPassed_pt2);

  hs->Draw("nostack");

  gPad->SetGrid(1, 0);
  gPad->BuildLegend(0.6, 0.6, 0.9, 0.9, "");
  
  c->SaveAs(saveDir + "/compLink_pt.png");
  c->Close();
  c = new TCanvas();

  //rinse and repeat for the polar angle instead of the transverse momentum
  hs = new THStack("hs", "Polar angle of linked generating electrons and all generating electrons;Polar angle (Rads);Count");
  
  realAll_PA1->SetLineColor(kBlue);
  realPassed_PA1->SetLineColor(kBlack);

  realAll_PA2->SetLineColor(kBlue);
  realPassed_PA2->SetLineColor(kBlack);

  realAll_PA2->SetLineStyle(6);
  realPassed_PA2->SetLineStyle(6);

  hs->Add(realAll_PA1);
  hs->Add(realPassed_PA1);

  hs->Add(realAll_PA2);
  hs->Add(realPassed_PA2);

  hs->Draw("nostack");

  gPad->SetGrid(1, 0);
  gPad->BuildLegend(0.35, 0.1, 0.65, 0.4, "");

  c->SaveAs(saveDir + "/compLink_PA.png");
  c->Close();
  c = new TCanvas();

  //rinse and repeat for the azimuth instead of the transverse momentum
  hs = new THStack("hs", "Azimuth of linked generating electrons and all generating electrons;Azimuth (Rads);Count");

  realAll_azimuth1->SetLineColor(kBlue);
  realPassed_azimuth1->SetLineColor(kBlack);

  realAll_azimuth2->SetLineColor(kBlue);
  realPassed_azimuth2->SetLineColor(kBlack);

  realAll_azimuth2->SetLineStyle(6);
  realPassed_azimuth2->SetLineStyle(6);

  hs->Add(realAll_azimuth1);
  hs->Add(realPassed_azimuth1);

  hs->Add(realAll_azimuth2);
  hs->Add(realPassed_azimuth2);

  hs->Draw("nostack");

  gPad->SetGrid(1, 0);
  gPad->BuildLegend(0.35, 0.1, 0.65, 0.4, "");

  c->SaveAs(saveDir + "/compLink_azimuth.png");
  c->Close();
  c = new TCanvas();
  
  //====================================================//
  //most likely linked particles regardless of reco type//
  //====================================================//

  //draws the two pt histograms on each other
  hs = new THStack("hs", "pt of linked particles and all generating electrons;Transverse momentum (GeV);Count");

  realAll_pt1->SetLineColor(kBlue);
  realAllPassed_pt1->SetLineColor(kBlack);

  realAll_pt2->SetLineColor(kBlue);
  realAllPassed_pt2->SetLineColor(kBlack);

  realAll_pt2->SetLineStyle(6);
  realAllPassed_pt2->SetLineStyle(6);

  hs->Add(realAll_pt1);
  hs->Add(realAllPassed_pt1);

  hs->Add(realAll_pt2);
  hs->Add(realAllPassed_pt2);
  
  hs->Draw("nostack");

  gPad->SetGrid(1, 0);
  gPad->BuildLegend(0.1, 0.1, 0.4, 0.4, "");

  c->SaveAs(saveDir + "/compAllLink_pt.png");
  c->Close();
  c = new TCanvas();

  //rinse and repeat for the polar angle instead of the transverse momentum
  hs = new THStack("hs", "Polar angle of linked particles and all generating electrons;Polar angle (Rads);Count");

  realAll_PA1->SetLineColor(kBlue);
  realAllPassed_PA1->SetLineColor(kBlack);

  realAll_PA2->SetLineColor(kBlue);
  realAllPassed_PA2->SetLineColor(kBlack);

  realAll_PA2->SetLineStyle(6);
  realAllPassed_PA2->SetLineStyle(6);

  hs->Add(realAll_PA1);
  hs->Add(realAllPassed_PA1);

  hs->Add(realAll_PA2);
  hs->Add(realAllPassed_PA2);

  hs->Draw("nostack");

  gPad->SetGrid(1, 0);
  gPad->BuildLegend(0.35, 0.1, 0.65, 0.4, "");

  c->SaveAs(saveDir + "/compAllLink_PA.png");
  c->Close();
  c = new TCanvas();

  //rinse and repeat for the azimuth instead of the transverse momentum
  hs = new THStack("hs", "Azimuth of linked particles and all generating electrons;Azimuth (Rads);Count");

  realAll_azimuth1->SetLineColor(kBlue);
  realAllPassed_azimuth1->SetLineColor(kBlack);

  realAll_azimuth2->SetLineColor(kBlue);
  realAllPassed_azimuth2->SetLineColor(kBlack);

  realAll_azimuth2->SetLineStyle(6);
  realAllPassed_azimuth2->SetLineStyle(6);

  hs->Add(realAll_azimuth1);
  hs->Add(realAllPassed_azimuth1);

  hs->Add(realAll_azimuth2);
  hs->Add(realAllPassed_azimuth2);

  hs->Draw("nostack");

  gPad->SetGrid(1, 0);
  gPad->BuildLegend(0.35, 0.1, 0.65, 0.4, "");

  c->SaveAs(saveDir + "/compAllLink_azimuth.png");
  c->Close();
  c = new TCanvas();

  //=======================================//
  //stacked reconstruction resolution hists//
  //=======================================//

  //draws the resolution of the transverse momentum of different linking methods
  hs = new THStack("hs", "Relative resolution of different linking methods;(Truth - reconstructed) / Truth p_{T};Count");

  relResOG_pt1->SetLineColor(kBlue);
  relResEPSum_pt1->SetLineColor(kRed);
  relResBestMatch_pt1->SetLineColor(kBlack);

  relResOG_pt2->SetLineColor(kBlue);
  relResEPSum_pt2->SetLineColor(kRed);
  relResBestMatch_pt2->SetLineColor(kBlack);

  relResOG_pt2->SetLineStyle(6);
  relResEPSum_pt2->SetLineStyle(6);
  relResBestMatch_pt2->SetLineStyle(6);
  
  hs->Add(relResOG_pt1);
  hs->Add(relResEPSum_pt1);
  hs->Add(relResBestMatch_pt1);

  hs->Add(relResOG_pt2);
  hs->Add(relResEPSum_pt2);
  hs->Add(relResBestMatch_pt2);
  
  hs->Draw("nostack");

  gPad->SetGrid(1, 0);
  gPad->BuildLegend(0.6, 0.6, 0.9, 0.9, "");

  c->SaveAs(saveDir + "/resolution/compResLink_pt.png");
  c->Close();
  c = new TCanvas();

  //draws the resolution of the polar angle of different linking methods
  hs = new THStack("hs", "Relative resolution of different linking methods;(Truth - reconstructed) / Truth polar angle;Count");

  relResOG_PA1->SetLineColor(kBlue);
  relResEPSum_PA1->SetLineColor(kRed);
  relResBestMatch_PA1->SetLineColor(kBlack);

  relResOG_PA2->SetLineColor(kBlue);
  relResEPSum_PA2->SetLineColor(kRed);
  relResBestMatch_PA2->SetLineColor(kBlack);

  relResOG_PA2->SetLineStyle(6);
  relResEPSum_PA2->SetLineStyle(6);
  relResBestMatch_PA2->SetLineStyle(6);

  hs->Add(relResOG_PA1);
  hs->Add(relResEPSum_PA1);
  hs->Add(relResBestMatch_PA1);

  hs->Add(relResOG_PA2);
  hs->Add(relResEPSum_PA2);
  hs->Add(relResBestMatch_PA2);

  hs->Draw("nostack");

  gPad->SetGrid(1, 0);
  gPad->BuildLegend(0.6, 0.6, 0.9, 0.9, "");

  c->SaveAs(saveDir + "/resolution/compResLink_PA.png");
  c->Close();
  c = new TCanvas();

  //draws the resolution of the azimuth of different linking methods
  hs = new THStack("hs", "Relative resolution of different linking methods;(Truth - reconstructed) / Truth azimuth;Count");

  relResOG_azimuth1->SetLineColor(kBlue);
  relResEPSum_azimuth1->SetLineColor(kRed);
  relResBestMatch_azimuth1->SetLineColor(kBlack);

  relResOG_azimuth2->SetLineColor(kBlue);
  relResEPSum_azimuth2->SetLineColor(kRed);
  relResBestMatch_azimuth2->SetLineColor(kBlack);

  relResOG_azimuth2->SetLineStyle(6);
  relResEPSum_azimuth2->SetLineStyle(6);
  relResBestMatch_azimuth2->SetLineStyle(6);

  hs->Add(relResOG_azimuth1);
  hs->Add(relResEPSum_azimuth1);
  hs->Add(relResBestMatch_azimuth1);

  hs->Add(relResOG_azimuth2);
  hs->Add(relResEPSum_azimuth2);
  hs->Add(relResBestMatch_azimuth2);
  
  hs->Draw("nostack");

  gPad->SetGrid(1, 0);
  gPad->BuildLegend(0.6, 0.6, 0.9, 0.9, "");

  c->SaveAs(saveDir + "/resolution/compResLink_azimuth.png");
  c->Close();
  c = new TCanvas();
  
  ///////////////
  //SOLO GRAPHS//
  ///////////////
  
  //draws the histogram that shows the distrobution of the linked electron weights
  hs = new THStack("hs", "Weight of linked electrons;Weight;Count");

  recoWeight_El1->SetLineColor(kBlue);
  recoWeight_El2->SetLineColor(kBlack);

  hs->Add(recoWeight_El1);
  hs->Add(recoWeight_El2);

  hs->Draw("nostack");

  gPad->SetGrid(1, 0);
  gPad->BuildLegend(0.1, 0.7, 0.4, 0.9, "");
  
  c->SaveAs(saveDir + "/weight/recoLinkElWeights.png");
  c->Close();
  c = new TCanvas();

  //draws the histogram that shows the distrobution of the relation weights for all particles
  hs = new THStack("hs", "Weight of linked electrons and photons;Weight;Count");

  recoWeight_BM1->SetLineColor(kBlue);
  recoWeight_BM2->SetLineColor(kBlack);

  hs->Add(recoWeight_BM1);
  hs->Add(recoWeight_BM2);

  hs->Draw("nostack");

  gPad->SetGrid(1, 0);
  gPad->BuildLegend(0.1, 0.7, 0.4, 0.9, "");

  c->SaveAs(saveDir + "/weight/recoLinkBMWeights.png");
  c->Close();
  c = new TCanvas();

  //draws the number of hits of reconstructed tracks for BM matching
  hs = new THStack("hs", "Number of hits of reconstructed tracks;Number of tracks;Count");

  numHits_BM1->SetLineColor(kBlue);
  numHits_BM2->SetLineColor(kBlack);

  hs->Add(numHits_BM1);
  hs->Add(numHits_BM2);

  hs->Draw("nostack");

  gPad->SetGrid(1, 0);
  gPad->BuildLegend(0.1, 0.7, 0.4, 0.9, "");

  c->SaveAs(saveDir + "/numHitsBM.png");
  c->Close();
  c = new TCanvas();
  
  //============================//
  //Resolution of electron linking//
  //============================//

  //draws a hist that hold the differences in electron pt
  hs = new THStack("hs", "Relative resolution in pt between all and reconstructed particles;(Truth - reconstructed) / Truth transverse momentum;Count");

  relResOG_pt1->SetLineColor(kBlue);
  relResOG_pt2->SetLineColor(kBlack);

  hs->Add(relResOG_pt1);
  hs->Add(relResOG_pt2);

  relResOG_pt2->SetLineStyle(1);

  hs->Draw("nostack");

  gPad->SetGrid(1, 0);
  gPad->BuildLegend(0.6, 0.7, 0.9, 0.9, "");

  c->SaveAs(saveDir + "/resolution/recoPt_diff.png");
  c->Close();
  c = new TCanvas();

  //draws a hist that hold the differences in electron polar angle
  hs = new THStack("hs", "Relative resolution in polar angle between all and reconstructed particles;(Truth - reconstructed) / Truth polar angle;Count");

  relResOG_PA1->SetLineColor(kBlue);
  relResOG_PA2->SetLineColor(kBlack);

  hs->Add(relResOG_PA1);
  hs->Add(relResOG_PA2);

  relResOG_PA2->SetLineStyle(1);

  hs->Draw("nostack");

  gPad->SetGrid(1, 0);
  gPad->BuildLegend(0.6, 0.7, 0.9, 0.9, "");

  c->SaveAs(saveDir + "/resolution/recoPA_diff.png");
  c->Close();
  c = new TCanvas();
  
  //draws a hist that hold the differences in electron azimuth
  hs = new THStack("hs", "Relative resolution in pt between all and reconstructed particles;(Truth - reconstructed) / Truth azimuth;Count");

  relResOG_azimuth1->SetLineColor(kBlue);
  relResOG_azimuth2->SetLineColor(kBlack);

  hs->Add(relResOG_azimuth1);
  hs->Add(relResOG_azimuth2);

  relResOG_azimuth2->SetLineStyle(1);

  hs->Draw("nostack");

  gPad->SetGrid(1, 0);
  gPad->BuildLegend(0.6, 0.7, 0.9, 0.9, "");

  c->SaveAs(saveDir + "/resolution/recoAzimuth_diff.png");
  c->Close();
  c = new TCanvas();

  //=============================================//
  //summed data from linked electrons and photons//
  //=============================================//

  //draws a hist that hold the electron and photon sum pt resolution
  hs = new THStack("hs", "All truth particle and sum of electrons and photons pt relative resolution;(Truth - reconstructed) / Truth p_{T};Count");

  relResEPSum_pt1->SetLineColor(kBlue);
  relResEPSum_pt2->SetLineColor(kBlack);

  hs->Add(relResEPSum_pt1);
  hs->Add(relResEPSum_pt2);

  relResEPSum_pt2->SetLineStyle(1);

  hs->Draw("nostack");

  gPad->SetGrid(1, 0);
  gPad->BuildLegend(0.6, 0.7, 0.9, 0.9, "");

  c->SaveAs(saveDir + "/resolution/relResEPSumLink_pt.png");
  c->Close();
  c = new TCanvas();
  
  //draws a hist that hold the electron and photon sum polar angle resolution
  hs = new THStack("hs", "All truth particle and sum of electrons and photons polar angle relative resolution;(Truth - reconstructed) / Truth polar angle;Count");

  relResEPSum_PA1->SetLineColor(kBlue);
  relResEPSum_PA2->SetLineColor(kBlack);

  hs->Add(relResEPSum_PA1);
  hs->Add(relResEPSum_PA2);

  relResEPSum_PA2->SetLineStyle(1);

  hs->Draw("nostack");

  gPad->SetGrid(1, 0);
  gPad->BuildLegend(0.6, 0.7, 0.9, 0.9, "");

  c->SaveAs(saveDir + "/resolution/relResEPSumLink_PA.png");
  c->Close();
  c = new TCanvas();
  
  //draws a hist that hold the electron and photon sum azimuth resolution
  hs = new THStack("hs", "All truth particle and sum of electrons and photons azimuth relative resolution;(Truth - reconstructed) / Truth azimuth;Count");

  relResEPSum_azimuth1->SetLineColor(kBlue);
  relResEPSum_azimuth2->SetLineColor(kBlack);

  hs->Add(relResEPSum_azimuth1);
  hs->Add(relResEPSum_azimuth2);

  relResEPSum_azimuth1->SetLineStyle(1);

  hs->Draw("nostack");

  gPad->SetGrid(1, 0);
  gPad->BuildLegend(0.6, 0.7, 0.9, 0.9, "");

  c->SaveAs(saveDir + "/resolution/relResEPSumLink_azimuth.png");
  c->Close();
  c = new TCanvas();

  //==========================//
  //Best match resolution data//
  //==========================//

  //draws a hist that hold the truth particle and the best match link regardless of reco type pt resolution
  hs = new THStack("hs", "Truth and best match link p_{T} relative resolution;(Truth - reconstructed) / Truth p_{T};Count");

  relResBestMatch_pt1->SetLineColor(kBlue);
  relResBestMatch_pt2->SetLineColor(kBlack);

  hs->Add(relResBestMatch_pt1);
  hs->Add(relResBestMatch_pt2);

  relResBestMatch_pt2->SetLineStyle(1); 

  hs->Draw("nostack");

  gPad->SetGrid(1, 0);
  gPad->BuildLegend(0.6, 0.7, 0.9, 0.9, "");

  c->SaveAs(saveDir + "/resolution/relResBestMatchLink_pt.png");
  c->Close();
  c = new TCanvas();

  //draws a hist that hold the truth particle and the best match link regardless of reco type polar angle resolution
  hs = new THStack("hs", "Truth and best match link polar angle relative resolution;(Truth - reconstructed) / Truth polar angle;Count");

  relResBestMatch_PA1->SetLineColor(kBlue);
  relResBestMatch_PA2->SetLineColor(kBlack);

  hs->Add(relResBestMatch_PA1);
  hs->Add(relResBestMatch_PA2);

  relResBestMatch_PA2->SetLineStyle(1);

  hs->Draw("nostack");

  gPad->SetGrid(1, 0);
  gPad->BuildLegend(0.6, 0.7, 0.9, 0.9, "");

  c->SaveAs(saveDir + "/resolution/relResBestMatchLink_PA.png");
  c->Close();
  c = new TCanvas();
  
  //draws a hist that hold the truth particle and the best match link regardless of reco type azimuth resolution
  hs = new THStack("hs", "Truth and best match link azimuth relative resolution;(Truth - reconstructed) / Truth azimuth;Count");

  relResBestMatch_azimuth1->SetLineColor(kBlue);
  relResBestMatch_azimuth2->SetLineColor(kBlack);

  hs->Add(relResBestMatch_azimuth1);
  hs->Add(relResBestMatch_azimuth2);

  relResBestMatch_azimuth2->SetLineStyle(1);

  hs->Draw("nostack");

  gPad->SetGrid(1, 0);
  gPad->BuildLegend(0.6, 0.7, 0.9, 0.9, "");

  c->SaveAs(saveDir + "/resolution/relResBestMatchLink_azimuth.png");
  c->Close();
  c = new TCanvas();
  
  //=============================================//
  //Best match resolution data for photon matches//
  //=============================================//

  //draws a hist that hold the truth particle and the best match link pt resolution for photons
  hs = new THStack("hs", "Truth and best match link p_{T} relative resolution for photons;(Truth - reconstructed) / Truth p_{T};Count");

  relResPhotonReco_pt1->SetLineColor(kBlue);
  relResPhotonReco_pt2->SetLineColor(kBlack);

  hs->Add(relResPhotonReco_pt1);
  hs->Add(relResPhotonReco_pt2);

  relResPhotonReco_pt2->SetLineStyle(1);

  hs->Draw("nostack");

  gPad->SetGrid(1, 0);
  gPad->BuildLegend(0.6, 0.7, 0.9, 0.9, "");

  c->SaveAs(saveDir + "/resolution/resBMPhotonLink_pt.png");
  c->Close();
  c = new TCanvas();

  //draws a hist that hold the truth particle and the best match link polar angle resolution for photons
  hs = new THStack("hs", "Truth and best match link polar angle relative resolution for photons;(Truth - reconstructed) / Truth polar angle;Count");

  relResPhotonReco_PA1->SetLineColor(kBlue);
  relResPhotonReco_PA2->SetLineColor(kBlack);

  hs->Add(relResPhotonReco_PA1);
  hs->Add(relResPhotonReco_PA2);

  relResPhotonReco_PA2->SetLineStyle(1);

  hs->Draw("nostack");

  gPad->SetGrid(1, 0);
  gPad->BuildLegend(0.6, 0.7, 0.9, 0.9, "");

  c->SaveAs(saveDir + "/resolution/resBMPhotonLink_PA.png");
  c->Close();
  c = new TCanvas();

  //draws a hist that hold the truth particle and the best match link azimuth resolution for photons
  hs = new THStack("hs", "Truth and best match link azimuth relative resolution for photons;(Truth - reconstructed) / Truth azimuth;Count");

  relResPhotonReco_azimuth1->SetLineColor(kBlue);
  relResPhotonReco_azimuth2->SetLineColor(kBlack);

  hs->Add(relResPhotonReco_azimuth1);
  hs->Add(relResPhotonReco_azimuth2);

  relResPhotonReco_azimuth2->SetLineStyle(1);

  hs->Draw("nostack");

  gPad->SetGrid(1, 0);
  gPad->BuildLegend(0.6, 0.7, 0.9, 0.9, "");

  c->SaveAs(saveDir + "/resolution/resBMPhotonLink_azimuth.png");
  c->Close();
  c = new TCanvas();

  //=============================================//
  //Best match resolution data for photon matches//
  //=============================================//

  //draws a hist that hold the truth particle and the best match link pt resolution for pions
  hs = new THStack("hs", "Truth and best match link p_{T} relative resolution for pions;(Truth - reconstructed) / Truth p_{T};Count");

  relResPionReco_pt1->SetLineColor(kBlue);
  relResPionReco_pt2->SetLineColor(kBlack);

  hs->Add(relResPionReco_pt1);
  hs->Add(relResPionReco_pt2);

  relResPionReco_pt2->SetLineStyle(1);

  hs->Draw("nostack");

  gPad->SetGrid(1, 0);
  gPad->BuildLegend(0.1, 0.7, 0.4, 0.9, "");

  c->SaveAs(saveDir + "/resolution/resBMPionLink_pt.png");
  c->Close();
  c = new TCanvas();

  //draws a hist that hold the truth particle and the best match link polar angle resolution for pions
  hs = new THStack("hs", "Truth and best match link polar angle relative resolution for pions;(Truth - reconstructed) / Truth polar angle;Count");

  relResPionReco_PA1->SetLineColor(kBlue);
  relResPionReco_PA2->SetLineColor(kBlack);

  hs->Add(relResPionReco_PA1);
  hs->Add(relResPionReco_PA2);

  relResPionReco_PA2->SetLineStyle(1);

  hs->Draw("nostack");

  gPad->SetGrid(1, 0);
  gPad->BuildLegend(0.6, 0.7, 0.9, 0.9, "");

  c->SaveAs(saveDir + "/resolution/resBMPionLink_PA.png");
  c->Close();
  c = new TCanvas();

  //draws a hist that hold the truth particle and the best match link azimuth resolution for pions
  hs = new THStack("hs", "Truth and best match link azimuth relative resolution for pions;(Truth - reconstructed) / Truth azimuth;Count");

  relResPionReco_azimuth1->SetLineColor(kBlue);
  relResPionReco_azimuth2->SetLineColor(kBlack);

  hs->Add(relResPionReco_azimuth1);
  hs->Add(relResPionReco_azimuth2);

  relResPionReco_azimuth2->SetLineStyle(1);

  hs->Draw("nostack");

  gPad->SetGrid(1, 0);
  gPad->BuildLegend(0.6, 0.7, 0.9, 0.9, "");

  c->SaveAs(saveDir + "/resolution/resBMPionLink_azimuth.png");
  c->Close();
  c = new TCanvas();
  
  //======================================================//
  //Resolution of track link reconstruction for BM photons//
  //======================================================//

  //draws a hist that hold the truth particle and the best match photon track link pt  resolution
  hs = new THStack("hs", "Truth and best match link photon track transverse momentum relative resolution;(Truth - reconstructed) / Truth p_{T};Count");

  relResTrackLink_pt1->SetLineColor(kBlue);
  relResTrackLink_pt2->SetLineColor(kBlack);

  hs->Add(relResTrackLink_pt1);
  hs->Add(relResTrackLink_pt2);

  relResTrackLink_pt2->SetLineStyle(1);

  hs->Draw("nostack");

  gPad->SetGrid(1, 0);
  gPad->BuildLegend(0.1, 0.7, 0.4, 0.9, "");

  c->SaveAs(saveDir + "/resolution/relResTrackLink_pt.png");
  c->Close();
  c = new TCanvas();
  
  ////////
  //MISC//
  ////////

  //draws a histogram for delta r (distance between polar angle and lambda on a plot between the real and reco data)
  //draws a hist that hold the truth particle and the best match photon track link pt  resolution
  hs = new THStack("hs", "Delta r between the real and reco data points;delta R (Rads);Count");

  deltaR_Hist1->SetLineColor(kBlue);
  deltaR_Hist2->SetLineColor(kBlack);

  hs->Add(deltaR_Hist1);
  hs->Add(deltaR_Hist2);

  hs->Draw("nostack");

  gPad->SetGrid(1, 0);
  gPad->BuildLegend(0.7, 0.7, 0.9, 0.9, "");

  c->SaveAs(saveDir + "/deltaR.png");
  c->Close();
  c = new TCanvas();

  ////////////////
  //EFFICIENCIES//
  ////////////////
  
  //draws an efficiency graph for transverse momentum
  TEfficiency* pt_Eff1 = new TEfficiency("pt_Eff1", "Efficiency of transverse momentum reconstruction;Transverse momentum (GeV);Efficiency", 40, 0, 2000);
  TFile* eff_File = new TFile("effFile.root", "recreate");

  //checks to make sure if the two histograms are compatable to make an efficiency histogram
  if(TEfficiency::CheckConsistency(*realPassed_pt1, *realAll_pt1)) {
    pt_Eff1->SetPassedHistogram(*realPassed_pt1, "f");
    pt_Eff1->SetTotalHistogram(*realAll_pt1, "f");
    eff_File->Write();
  }

  TEfficiency* pt_Eff2 = new TEfficiency("pt_Eff2", "Efficiency of transverse momentum reconstruction;Transverse momentum (GeV);Efficiency", 40, 0, 2000);
  eff_File = new TFile("effFile.root", "recreate");
  if(TEfficiency::CheckConsistency(*realPassed_pt2, *realAll_pt2)) {
    pt_Eff2->SetPassedHistogram(*realPassed_pt2, "f");
    pt_Eff2->SetTotalHistogram(*realAll_pt2, "f");
    eff_File->Write();
  }

  //draws an efficiency graph for the polar angle
  TEfficiency* PA_Eff1 = new TEfficiency("PA_Eff1", "Efficiency of polar angle reconstruction;Polar angle (Rads);Efficiency", 40, 0, 1.6);
  eff_File = new TFile("effFile.root", "recreate");
  
  //checks to make sure if the two histograms are compatable to make an efficiency histogram
  if(TEfficiency::CheckConsistency(*realPassed_PA1, *realAll_PA1)) {
    PA_Eff1->SetPassedHistogram(*realPassed_PA1, "f");
    PA_Eff1->SetTotalHistogram(*realAll_PA1, "f");
    eff_File->Write();
  }

  TEfficiency* PA_Eff2 = new TEfficiency("PA_Eff2", "Efficiency of polar angle reconstruction;Polar angle (Rads);Efficiency", 40, 0, 2000);
  eff_File = new TFile("effFile.root", "recreate");
  if(TEfficiency::CheckConsistency(*realPassed_PA2, *realAll_PA2)) {
    PA_Eff2->SetPassedHistogram(*realPassed_PA2, "f");
    PA_Eff2->SetTotalHistogram(*realAll_PA2, "f");
    eff_File->Write();
  }
  
  //draws an efficiency graph for the azimuth
  TEfficiency* azimuth_Eff1 = new TEfficiency("azimuth_Eff", "Efficiency of azimuth reconstruction;Azimuth (Rads);Efficiency", 40, -1.6, 1.6);
  eff_File = new TFile("effFile.root", "recreate");

  //checks to make sure if the two histograms are compatable to make an efficiency histogram
  if(TEfficiency::CheckConsistency(*realPassed_azimuth1, *realAll_azimuth1)) {
    azimuth_Eff1->SetPassedHistogram(*realPassed_azimuth1, "f");
    azimuth_Eff1->SetTotalHistogram(*realAll_azimuth1, "f");
    eff_File->Write();
  }

  TEfficiency* azimuth_Eff2 = new TEfficiency("azimuth_Eff2", "Efficiency of azimuth reconstruction;Azimuth (GeV);Efficiency", 40, 0, 2000);
  eff_File = new TFile("effFile.root", "recreate");
  if(TEfficiency::CheckConsistency(*realPassed_azimuth2, *realAll_azimuth2)) {
    azimuth_Eff2->SetPassedHistogram(*realPassed_azimuth2, "f");
    azimuth_Eff2->SetTotalHistogram(*realAll_azimuth2, "f");
    eff_File->Write();
  }
  
  //draws an efficiency graph for transverse momentum for best match particles
  TEfficiency* ptAll_Eff1 = new TEfficiency("ptBM_Eff", "Efficiency of transverse momentum reconstruction for all particles;Transverse momentum (GeV);Efficiency", 40, 0, 2000);
  eff_File = new TFile("effFile.root", "recreate");
  
  //checks to make sure if the two histograms are compatable to make an efficiency histogram
  if(TEfficiency::CheckConsistency(*realAllPassed_pt1, *realAll_pt1)) {
    ptAll_Eff1->SetPassedHistogram(*realAllPassed_pt1, "f");
    ptAll_Eff1->SetTotalHistogram(*realAll_pt1, "f");
    eff_File->Write();
  }

  TEfficiency* ptAll_Eff2 = new TEfficiency("ptAll_Eff2", "Efficiency of transverse momentum reconstruction;Transverse momentum (GeV);Efficiency", 40, 0, 2000);
  eff_File = new TFile("effFile.root", "recreate");
  if(TEfficiency::CheckConsistency(*realAllPassed_pt2, *realAll_pt2)) {
    ptAll_Eff2->SetPassedHistogram(*realAllPassed_pt2, "f");
    ptAll_Eff2->SetTotalHistogram(*realAll_pt2, "f");
    eff_File->Write();
  }

  //draws an efficiency graph for the polar angle for best match particles
  TEfficiency* PAAll_Eff1 = new TEfficiency("PABM_Eff", "Efficiency of polar angle reconstruction for all particles;Polar angle (Rads);Efficiency", 40, 0, 1.6);
  eff_File = new TFile("effFile.root", "recreate");

  //checks to make sure if the two histograms are compatable to make an efficiency histogram
  if(TEfficiency::CheckConsistency(*realAllPassed_PA1, *realAll_PA1)) {
    PAAll_Eff1->SetPassedHistogram(*realAllPassed_PA1, "f");
    PAAll_Eff1->SetTotalHistogram(*realAll_PA1, "f");
    eff_File->Write();
  }

  TEfficiency* PAAll_Eff2 = new TEfficiency("PAAll_Eff2", "Efficiency of polar angle reconstruction;Polar angle (Rads);Efficiency", 40, 0, 2000);
  eff_File = new TFile("effFile.root", "recreate");
  if(TEfficiency::CheckConsistency(*realAllPassed_PA2, *realAll_PA2)) {
    PAAll_Eff2->SetPassedHistogram(*realAllPassed_PA2, "f");
    PAAll_Eff2->SetTotalHistogram(*realAll_PA2, "f");
    eff_File->Write();
  }

  //draws an efficiency graph for the azimuth for best match particles
  TEfficiency* azimuthAll_Eff1 = new TEfficiency("azimuthBM_Eff", "Efficiency of azimuth reconstruction for all particles;Azimuth (Rads);Efficiency", 40, -1.6, 1.6);
  eff_File = new TFile("effFile.root", "recreate");

  //checks to make sure if the two histograms are compatable to make an efficiency histogram
  if(TEfficiency::CheckConsistency(*realAllPassed_azimuth1, *realAll_azimuth1)) {
    azimuthAll_Eff1->SetPassedHistogram(*realAllPassed_azimuth1, "f");
    azimuthAll_Eff1->SetTotalHistogram(*realAll_azimuth1, "f");
    eff_File->Write();
  }

  TEfficiency* azimuthAll_Eff2 = new TEfficiency("azimuthAll_Eff2", "Efficiency of azimuth reconstruction;Azimuth (Rads);Efficiency", 40, 0, 2000);
  eff_File = new TFile("effFile.root", "recreate");
  if(TEfficiency::CheckConsistency(*realAllPassed_azimuth2, *realAll_azimuth2)) {
    azimuthAll_Eff2->SetPassedHistogram(*realAllPassed_azimuth2, "f");
    azimuthAll_Eff2->SetTotalHistogram(*realAll_azimuth2, "f");
    eff_File->Write();
  }

  //setting new line colors and styles so the plots are differentable when put on the same graph
  pt_Eff1->SetLineColor(4);
  PA_Eff1->SetLineColor(4);
  azimuth_Eff1->SetLineColor(4);
  
  ptAll_Eff1->SetLineColor(2);
  PAAll_Eff1->SetLineColor(2);
  azimuthAll_Eff1->SetLineColor(2);

  pt_Eff2->SetLineColor(4);
  PA_Eff2->SetLineColor(4);
  azimuth_Eff2->SetLineColor(4);

  ptAll_Eff2->SetLineColor(2);
  PAAll_Eff2->SetLineColor(2);
  azimuthAll_Eff2->SetLineColor(2);

  pt_Eff2->SetLineStyle(6);
  PA_Eff2->SetLineStyle(6);
  azimuth_Eff2->SetLineStyle(6);

  ptAll_Eff2->SetLineStyle(6);
  PAAll_Eff2->SetLineStyle(6);
  azimuthAll_Eff2->SetLineStyle(6);

  pt_Eff1->SetTitle(filePrefix1 + ":Just electrons");
  PA_Eff1->SetTitle(filePrefix1 + ":Just electrons");
  azimuth_Eff1->SetTitle(filePrefix1 + ":Just electrons");

  ptAll_Eff1->SetTitle(filePrefix1 + ":Electrons and photons");
  PAAll_Eff1->SetTitle(filePrefix1 + ":Electrons and photons");
  azimuthAll_Eff1->SetTitle(filePrefix1 + ":Electrons and photons");

  pt_Eff2->SetTitle(filePrefix2 + ":Just electrons");
  PA_Eff2->SetTitle(filePrefix2 + ":Just electrons");
  azimuth_Eff2->SetTitle(filePrefix2 + ":Just electrons");

  ptAll_Eff2->SetTitle(filePrefix2 + ":Electrons and photons");
  PAAll_Eff2->SetTitle(filePrefix2 + ":Electrons and photons");
  azimuthAll_Eff2->SetTitle(filePrefix2 + ":Electrons and photons");

  //creates a multigraph for the transverse momentum efficiencies
  TMultiGraph *mg = new TMultiGraph();
  c = new TCanvas();
  
  mg->Add(pt_Eff1->CreateGraph());
  mg->Add(ptAll_Eff1->CreateGraph());

  mg->Add(pt_Eff2->CreateGraph());
  mg->Add(ptAll_Eff2->CreateGraph());

  mg->SetTitle("Transverse momentum reconstruction efficiencies");
  mg->GetXaxis()->SetTitle("p_{T} (GeV)");
  mg->GetYaxis()->SetTitle("Efficiency");

  mg->Draw("aZ");

  gPad->SetGrid(1, 0);
  gPad->BuildLegend(0.1, 0.1, 0.4, 0.4, "");

  c->SaveAs(saveDir + "/efficiency/multiBMEff_pt.png");
  c->Close();

  //creates a multigraph for the transverse momentum efficiencies
  mg = new TMultiGraph();
  c = new TCanvas();

  mg->Add(pt_Eff2->CreateGraph());
  mg->Add(ptAll_Eff2->CreateGraph());

  mg->SetTitle("Transverse momentum reconstruction efficiencies");
  mg->GetXaxis()->SetTitle("Transverse momentum (GeV)");
  mg->GetYaxis()->SetTitle("Efficiency");

  mg->Draw("aZ");

  gPad->SetGrid(1, 0);
  gPad->BuildLegend(0.65, 0.3, 0.95, 0.6, "");

  c->SaveAs(saveDir + "/efficiency/multiBMEffConf_pt.png");
  c->Close();

  //creates a multigraph for the polar angle efficiencies
  mg = new TMultiGraph();
  c = new TCanvas();

  mg->Add(PA_Eff1->CreateGraph());
  mg->Add(PAAll_Eff1->CreateGraph());

  mg->Add(PA_Eff2->CreateGraph());
  mg->Add(PAAll_Eff2->CreateGraph());

  mg->SetTitle("Polar angle reconstruction efficiencies");
  mg->GetXaxis()->SetTitle("Polar angle (Rads)");
  mg->GetYaxis()->SetTitle("Efficiency");
  
  mg->Draw("aZ");

  gPad->SetGrid(1, 0);
  gPad->BuildLegend(0.35, 0.1, 0.65, 0.4, "");

  c->SaveAs(saveDir + "/efficiency/multiBMEff_PA.png");
  c->Close();

  //creates a multigraph for the azimuth efficiencies
  mg = new TMultiGraph();
  c = new TCanvas();

  mg->Add(azimuth_Eff1->CreateGraph());
  mg->Add(azimuthAll_Eff1->CreateGraph());

  mg->Add(azimuth_Eff2->CreateGraph());
  mg->Add(azimuthAll_Eff2->CreateGraph());

  mg->SetTitle("Azimuth reconstruction efficiencies");
  mg->GetXaxis()->SetTitle("Azimuth (Rads)");
  mg->GetYaxis()->SetTitle("Efficiency");
  
  mg->Draw("aZ");

  gPad->SetGrid(1, 0);
  gPad->BuildLegend(0.35, 0.1, 0.65, 0.4, "");

  c->SaveAs(saveDir + "/efficiency/multiBMEff_azimuth.png");
  c->Close();
}
