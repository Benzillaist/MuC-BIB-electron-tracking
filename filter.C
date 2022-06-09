void filter()
{
  gStyle->SetOptTitle(0);
  
  auto fileName = "ntuple.root";
  auto treeName = "MyLCTuple";

  TFile *myFile = new TFile("ntuple.root");
  TTree *myTree = (TTree*)myFile->Get("MyLCTuple");

  ROOT::RDataFrame d(treeName, fileName, {"rcmox"});

  TCanvas *c = new TCanvas();

  //Creating reconstructed transverse momentum
  c->SetTitle("Electron reconstructed transverse momentum");
  myTree->Draw("sqrt(rcmox*rcmox + rcmoy*rcmoy)>>reco_pt(20, 0, 2000)", "rctyp == 11", "");
  TH1 *reco_pt = (TH1*)gDirectory->Get("reco_pt");
  reco_pt->GetXaxis()->SetTitle("Transverse momentum");
  reco_pt->GetYaxis()->SetTitle("Particle count");
  c->SaveAs("reco_pt.png");

  c->Clear();

  //Creating real particle pt graph
  c->SetTitle("Electron transverse momentum");
  myTree->Draw("sqrt(mcmox*mcmox + mcmoy*mcmoy)>>real_pt(20, 0, 2000)", "rctyp == 11", "");
  TH1 *real_pt = (TH1*)gDirectory->Get("real_pt");
  real_pt->GetXaxis()->SetTitle("Transverse momentum");
  real_pt->GetYaxis()->SetTitle("Particle count");
  c->SaveAs("real_pt.png");

  c->Clear();

  //Creating reco azimuth graph
  c->SetTitle("Electron reconstructed azimuth");
  myTree->Draw("atan(rcmoz / sqrt(rcmox*rcmox + rcmoy*rcmoy))>>reco_azimuth(20, -1.6, 1.6)", "rctyp == 11", "");
  TH1 *reco_azimuth = (TH1*)gDirectory->Get("reco_azimuth");
  reco_azimuth->GetXaxis()->SetTitle("Azimuth (Rads)");
  reco_azimuth->GetYaxis()->SetTitle("Particle count");
  c->SaveAs("reco_azimuth.png");

  c->Clear();
  
  //Creating real particle azimuth graph
  c->SetTitle("Electron reconstructed azimuth");
  myTree->Draw("atan(mcmoz/ sqrt(mcmox*mcmox + mcmoy*mcmoy))>>real_azimuth(20, -1.6, 1.6)", "rctyp == 11", "");
  TH1 *real_azimuth = (TH1*)gDirectory->Get("real_azimuth");
  real_azimuth->GetXaxis()->SetTitle("Azimuth (Rads)");
  real_azimuth->GetYaxis()->SetTitle("Particle count");
  c->SaveAs("real_azimuth.png");

  c->Clear();

  auto pt_Stack = new THStack("hs", "Stack of transverse momenta");

  gStyle->SetPalette(kRust);

  reco_pt->Draw("PLC");
  real_pt->Draw("PLC SAME");

  auto legend = new TLegend(0.1, 0.1, 0.48, 0.3);
  legend->SetHeader("Transverse momenta", "C");
  legend->AddEntry("reco_pt", "Reconstructed");
  legend->AddEntry("real_pt", "Real");
  legend->Draw();

  TEfficiency* pt_Eff = 0;
  TFile* eff_File = new TFile("effFile.root", "recreate");
  
  if(TEfficiency::CheckConsistency(*reco_pt, *real_pt)) {
    pt_Eff = new TEfficiency(*reco_pt, *real_pt);
    eff_File->Write();
  } 
   
}
