void testIoniz()
{ 
  gSystem->Load("lib/libretrim.so"); 
//  retrim::CollisionReader c("data/coll/COLLISON_200keV.txt"); 
  retrim::CollisionReader c("data/c/60torr/coll/COLLISON_150.00keV.txt"); 
  retrim::TableReader * tables[3]; 
  tables[1] = new retrim::TableReader("data/srim/f_in_cf4_60torr.txt"); 
  tables[0] = new retrim::TableReader("data/srim/c_in_cf4_60torr.txt"); 

  retrim::TrackMaker t(&c,2,tables,0.98); 

 // t.setVerbose(1); 
  t.enableHistograms(6e-2); 
  for (int i = 0; i < c.nTracks(); i++)
  {
    t.makeTrack(i); 
  }

  TCanvas *c1 = new TCanvas; 

  c1->Divide(2,2); 
  c1->cd(1); 
  t.getIonizationIon()->Scale(1./t.nHistogram()); 
  t.getIonizationRecoils()->Scale(1./t.nHistogram()); 
  t.getIonizationRecoils()->SetLineColor(2); 
  t.getIonizationIon()->DrawCopy(); 
  t.getIonizationRecoils()->DrawCopy("same"); 

//  c1->cd(2); 
//  t.getStats()->Draw("ElecLoss:predictedElecLoss","predictedElecLoss < 0.1","colz"); 

  c1->cd(2); 

  TH1F * true_ion = 0; 
  TH1F * true_recoil = 0; 

  retrim::HistogramReader::readIONIZ("data/c/60torr/ioniz/IONIZ_150.00keV.txt",&true_ion, &true_recoil); 

  true_ion->Draw(); 
  true_recoil->SetLineColor(2); 
  true_recoil->Draw("same"); 

  c1->cd(3); 

  t.getPhononsRecoils()->Scale(1./t.nHistogram()); 
  t.getPhononsIon()->Scale(1./t.nHistogram()); 
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

  retrim::HistogramReader::readPHONON(TString::Format("data/c/60torr/phonon/PHONON150.00keV.txt"),&phonon_true_ion, &phonon_true_recoil); 
  TH1F * phonon_true_sum = phonon_true_ion->Clone("phonon_truesum"); 
  phonon_true_sum->SetTitle("TRIM Phonons"); 
  phonon_true_sum->Add(phonon_true_recoil); 
  phonon_true_sum->SetLineColor(3); 
  phonon_true_sum->Draw(); 
  phonon_true_recoil->SetLineColor(2); 
  phonon_true_ion->Draw("same"); 
  phonon_true_recoil->Draw("same"); 


//  printf("\n\nIONIZATION:\n\nTRUTH: Total: %f\t Ions: %f\t Recoils: %f\nRECO: Total: %f\t Ions: %f\t Recoils %f:\n", true_sum->Integral(), true_ion->Integral(), true_recoil->Integral(), sum->Integral(), t.getIonizationIon()->Integral(), t.getIonizationRecoils()->Integral()); 
//  printf("\n\nPHONONS:\n\nTRUTH: Total: %f\t Ions: %f\t Recoils: %f\nRECO: Total: %f\t Ions: %f\t Recoils %f:\n", phonon_true_sum->Integral(), phonon_true_ion->Integral(), phonon_true_recoil->Integral(), phononsum->Integral(), t.getPhononsIon()->Integral(), t.getPhononsRecoils()->Integral()); 


 


}
