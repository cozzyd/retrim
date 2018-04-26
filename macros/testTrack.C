void testTrack(int which=0, const char * save=0, const char * collfile = "data/f/30torr/coll/COLLISON_200.00keV.txt")
{ 
  gStyle->SetOptStat(0); 
  gSystem->Load("lib/libretrim.so"); 
  retrim::CollisionReader c(collfile); 
  retrim::TableReader * tables[2]; 
  tables[0] = new retrim::TableReader("data/srim/f_in_cf4_30torr.txt"); 
  tables[1] = new retrim::TableReader("data/srim/c_in_cf4_30torr.txt"); 

  retrim::TrackMaker t(&c,2,tables); 
  t.enableHistograms(2.5e-2); 
  t.makeTrack(which); 
  retrim::SimpleIonizationModel m; 
  m.setSeed(0); 
  cout << t.makeElectrons(&m) << endl; 

  TCanvas * ct = new TCanvas("track","track"); 
  ct->cd(); 
  t.draw("s"); 

  TCanvas * ca = new TCanvas(TString::Format("c%d",which), "Canvas", 800,800); 
  ca->Divide(2,2); 

  ca->cd(1); 
  t.draw("e"); 
  //gPad->GetView()->ShowAxis(); 
  //gPad->GetView()->Top(); 
  

  ca->cd(2); 
  gStyle->SetNumberContours(255); 
  TH2 * img = new TH2I (TString::Format("img%d",which),"Electrons Generated", 50,-2.5,2.5,50,-2.5,2.5); 
  img->GetXaxis()->SetTitle("x (mm)"); 
  img->GetYaxis()->SetTitle("y (mm)"); 
  t.fillElectrons(img); 
  img->Draw("colz"); 

  ca->cd(3); 
  img->ProjectionX()->Draw(); 
  ca->cd(4); 
  TH1 * total = t.getIonizationIon()->Clone("total"); 
  total->SetTitle("Total Ionization"); 
  total->Add(t.getIonizationRecoils()); 
  total->SetLineColor(3); 
  total->DrawCopy(); 
  t.getIonizationIon()->DrawCopy("same"); 
  t.getIonizationRecoils()->DrawCopy("same"); 
  
  if (save) ca->SaveAs(save); 

//  t.getStats()->Draw("(ElecLoss-predictedElecLoss)/ElecLoss"); 
}
