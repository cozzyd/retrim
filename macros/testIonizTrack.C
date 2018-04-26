
void testIonizTrack(int track = 1, TH1 * diff = 0, int energy=150)
{ 
  gSystem->Load("lib/libretrim.so"); 
  gStyle->SetOptStat(0); 
  retrim::CollisionReader c(TString::Format("data/f/30torr/many/%dkeV/COLLISON%d.txt",energy,track)); 
  retrim::TableReader * tables[2]; 

  tables[0] = new retrim::TableReader("data/srim/f_in_cf4_30torr.txt"); 
  tables[1] = new retrim::TableReader("data/srim/c_in_cf4_30torr.txt"); 

  
  retrim::TrackMaker t(&c,2, tables,0.97); 
  t.setVerbose(1); 
  t.enableHistograms(1e-1); 
//  t.enableLeftoverCut(1); 
  t.makeTrack(0); 

  TCanvas *c1 = new TCanvas; 
  TCanvas *c2 = new TCanvas; 
  c2->cd(); 
  t.draw(); 

  c1->Divide(4,1); 
  c1->cd(1); 

  TH1 * sum = t.getIonizationRecoils()->Clone("recosum"); 
  sum->SetTitle("Reco Ionization"); 
  t.getIonizationRecoils()->SetLineColor(2); 
  sum->Add(t.getIonizationIon()); 
  sum->SetLineColor(3); 
  sum->DrawCopy(); 
  t.getIonizationIon()->DrawCopy("same"); 
  t.getIonizationRecoils()->DrawCopy("same"); 

//  c1->cd(2); 
//  t.getStats()->Draw("ElecLoss:predictedElecLoss","predictedElecLoss < 0.1","colz"); 

  c1->cd(2); 
 
  TH1F * true_ion = 0; 
  TH1F * true_recoil = 0; 

  retrim::HistogramReader::readIONIZ(TString::Format("data/f/30torr/many/%dkeV/IONIZ%d.txt",energy,track),&true_ion, &true_recoil); 
  TH1F * true_sum = true_ion->Clone("truesum"); 
  true_sum->SetTitle("TRIM Ionization"); 
  true_sum->Add(true_recoil); 
  true_sum->SetLineColor(3); 
  true_sum->Draw(); 
  true_recoil->SetLineColor(2); 
  true_ion->Draw("same"); 
  true_recoil->Draw("same"); 


  c1->cd(3); 

  TH1 * phononsum = t.getPhononsRecoils()->Clone("recosumphonon"); 
  phononsum->SetTitle("Reco Phonons"); 
  t.getPhononsRecoils()->SetLineColor(2); 
  phononsum->Add(t.getPhononsIon()); 
  phononsum->SetLineColor(3); 
  phononsum->DrawCopy(); 
  t.getPhononsIon()->DrawCopy("same"); 
  t.getPhononsRecoils()->DrawCopy("same"); 

  c1->cd(4); 

  TH1F * phonon_true_ion = 0; 
  TH1F * phonon_true_recoil = 0; 

  retrim::HistogramReader::readPHONON(TString::Format("data/f/30torr/many/%dkeV/PHONON%d.txt",energy,track),&phonon_true_ion, &phonon_true_recoil); 
  TH1F * phonon_true_sum = phonon_true_ion->Clone("phonon_truesum"); 
  phonon_true_sum->SetTitle("TRIM Phonons"); 
  phonon_true_sum->Add(phonon_true_recoil); 
  phonon_true_sum->SetLineColor(3); 
  phonon_true_sum->Draw(); 
  phonon_true_recoil->SetLineColor(2); 
  phonon_true_ion->Draw("same"); 
  phonon_true_recoil->Draw("same"); 


  printf("\n\nIONIZATION:\n\nTRUTH: Total: %f\t Ions: %f\t Recoils: %f\nRECO: Total: %f\t Ions: %f\t Recoils %f:\n", true_sum->Integral(), true_ion->Integral(), true_recoil->Integral(), sum->Integral(), t.getIonizationIon()->Integral(), t.getIonizationRecoils()->Integral()); 
  printf("\n\nPHONONS:\n\nTRUTH: Total: %f\t Ions: %f\t Recoils: %f\nRECO: Total: %f\t Ions: %f\t Recoils %f:\n", phonon_true_sum->Integral(), phonon_true_ion->Integral(), phonon_true_recoil->Integral(), phononsum->Integral(), t.getPhononsIon()->Integral(), t.getPhononsRecoils()->Integral()); 

}
