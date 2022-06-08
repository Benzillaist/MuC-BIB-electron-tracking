void filter()
{
  auto fileName = "ntuple.root";
  auto treeName = "MyLCTuple";

  auto myHist = new TH1F("h1", "ntuple", 6141, -4, 4);

  TFile *myFile = new TFile("ntuple.root");
  TTree *myTree = (TTree*)myFile->Get("MyLCTuple");

  myTree->Draw("sqrt(rcmox*rcmox + rcmoy*rcmoy)");
}
