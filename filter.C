void filter()
{
  gStyle->SetOptTitle(0);
  
  auto fileName = "ntuple.root";
  auto treeName = "MyLCTuple";

  TFile *myFile = new TFile("ntuple.root");
  TTree *myTree = (TTree*)myFile->Get("MyLCTuple");

  ROOT::RDataFrame d(treeName, fileName, {"rcmox"});

  TCanvas *c = new TCanvas("c", "Test Title");

  //Creating reco pt graph
  c->SetTitle("Electron reconstructed transverse momentum");
  myTree->Draw("sqrt(rcmox*rcmox + rcmoy*rcmoy)>>Reconstructed_pt(20, 0, 2000)", "rctyp == 11", "");
  c->SaveAs("reco_pt.png");

  c->Clear();

  //Creating real particle pt graph
  c->SetTitle("Electron transverse momentum");
  myTree->Draw("sqrt(mcmox*mcmox + mcmoy*mcmoy)>>Real_pt(20, 0, 2000)", "rctyp == 11", "");
  c->SaveAs("real_pt.png");

  c->Clear();

  //Creating reco azimuth graph
  c->SetTitle("Electron reconstructed azimuth");
  myTree->Draw("atan(sqrt(rcmox*rcmox + rcmoy*rcmoy) / rcmoz)>>Reco_azimuth(20, -1.6, 1.6)", "rctyp == 11", "");  
  c->SaveAs("reco_azimuth.png");

  c->Clear();
  
  //Creating real particle azimuth graph
  c->SetTitle("Electron reconstructed azimuth");
  myTree->Draw("atan(sqrt(mcmox*mcmox + mcmoy*mcmoy) / mcmoz)>>Real_azimuth(20, -1.6, 1.6)", "rctyp == 11", "");
  c->SaveAs("real_azimuth.png");
  
  //auto def_rcmpt = d.Define("rcmpt", [](Float_t rcmox, Float_t rcmoy) { return sqrt(rcmox[0]*rcmox[0] + rcmoy[0]*rcmoy[0]); }, {"rcmox", "rcmoy"});

}
