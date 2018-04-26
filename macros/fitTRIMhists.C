void fitTRIMhists(const char * out, const char * dir, int max=553) 
{
  gSystem->Load("lib/libretrim.so"); 

  TFile of(out,"RECREATE"); 

  double diff, diff1; 
  double integral; 

  TTree * tree = new TTree("truth","truth"); 
  TH1F * ion_ioniz = 0; 
  TH1F * recoil_ioniz = 0; 
  TH1F * total_ioniz = 0; 

  retrim::HistogramReader::readIONIZ(TString::Format("%s/IONIZ0.txt",dir), &ion_ioniz, &recoil_ioniz); 

  total_ioniz = (TH1F*) ion_ioniz->Clone("totalIoniz"); 
  total_ioniz->Add(recoil_ioniz); 
  total_ioniz->SetTitle("total_ionization"); 


  int lastbin; 
  tree->Branch("diff",&diff);  
  tree->Branch("integral",&integral);  
  tree->Branch("lastbin",&lastbin);  
  tree->Branch("diff1",&diff1);  
  tree->Branch("total_ioniz",&total_ioniz);  
  tree->Branch("ion_ioniz",&ion_ioniz);  
  tree->Branch("recoil_ioniz",&recoil_ioniz);  


  for (int i = 0; i < max; i++) 
  {

    if (i > 0)
    {
      retrim::HistogramReader::readIONIZ(TString::Format("%s/IONIZ%d.txt",dir,i), &ion_ioniz, &recoil_ioniz); 
      total_ioniz->Reset(); 
      total_ioniz->Add(ion_ioniz); 
      total_ioniz->Add(recoil_ioniz); 
    }

    lastbin = total_ioniz->GetNbinsX(); 
    for (int b = lastbin; b > 0; b--)
    {
      if (total_ioniz->GetBinContent(b) > 0)
      {
        lastbin = b; 
        break; 
      }
    }

    diff = total_ioniz->Integral(1,lastbin/2) - total_ioniz->Integral((lastbin+1)/2+1,lastbin); 
    diff1 = total_ioniz->Integral(1,(lastbin+1)/2) - total_ioniz->Integral((lastbin+1)/2+1,lastbin); 
    integral = total_ioniz->Integral(); 


    tree->Fill(); 

  }


  tree->Write(); 
  of.Close(); 




}
