#include "include/CollisionReader.hh" 
#include "include/TableReader.hh" 
#include "include/HistogramReader.hh" 
#include "include/TrackMaker.hh" 
#include "include/IonizationModel.hh" 
#include <dirent.h>


/* You probably need to run .L macros/loader.C+ first */ 

void sven(const char * outfile = "f_in_20TorrSF6_740TorrHe.root", const char * dir = "data_he/f/20+740torr/")
{


  /** Setup output*/ 

  TFile out(outfile,"RECREATE"); 

  const std::vector<TLorentzVector> * x = 0;  
  const std::vector<TLorentzVector> * p = 0;  
  TH1F * trim_ioniz_recoil = 0; 
  TH1F * trim_ioniz_ion = 0; 
  TH1F * reco_ioniz_recoil = 0; 
  TH1F * reco_ioniz_ion = 0; 
  int e = 0; 
  int i= 0; 


  TTree * tree = new TTree("sf6_20torr_he_740torr","20 Torr SF6, 740 Torr He"); 

  tree->Branch("electron_x",&x); 
  tree->Branch("electron_p",&p); 
  tree->Branch("trim_ioniz_recoil",&trim_ioniz_recoil); 
  tree->Branch("trim_ioniz_ion",&trim_ioniz_ion); 
  tree->Branch("reco_ioniz_recoil",&reco_ioniz_recoil); 
  tree->Branch("reco_ioniz_ion",&reco_ioniz_ion); 
  tree->Branch("EkeV",&e); 
  tree->Branch("indexAtEnergy",&i); 

  const retrim::TableReader * tables[3]; 
  tables[0] = new retrim::TableReader("data_he/srim/F_in_SF6_20Torr_He_740Torr.txt"); 
  tables[1] = new retrim::TableReader("data_he/srim/S_in_SF6_20Torr_He_740Torr.txt"); 
  tables[2] = new retrim::TableReader("data_he/srim/He_in_SF6_20Torr_He_740Torr.txt"); 
  retrim::SimpleIonizationModel m; 



  TString colldir; 
  colldir.Form("%s/coll",dir); 
  DIR * d = opendir(colldir.Data()); 
  TString iondir; 
  iondir.Form("%s/ioniz",dir); 
  TString phondir; 
  phondir.Form("%s/phonon",dir); 

  while(dirent * ent = readdir(d))
  {

    if (ent->d_name[0]=='.') continue; 
    printf("%s\n",ent->d_name); 
    sscanf(ent->d_name,"COLLISON_%d.00keV.%d.txt",&e,&i); 
    printf("%d keV:%d",e,i); 

    retrim::CollisionReader c(TString::Format("%s/COLLISON_%d.00keV.%d.txt",colldir.Data(), e,i)); 
    retrim::HistogramReader::readIONIZ(TString::Format("%s/IONIZ_%d.00keV.%d.txt",iondir.Data(), e,i),&trim_ioniz_recoil, &trim_ioniz_ion); 
    retrim::TrackMaker t(&c,3,&tables[0],1); 

    double width = trim_ioniz_recoil->GetBinWidth(1); 
    t.enableHistograms(width); 
//    t.setVerbose(true); 
//    gSystem->SetFPEMask(kInvalid); 
    t.makeTrack(0); 
 //   gSystem->SetFPEMask(kDefault); 
    t.makeElectrons(&m); 
    x =  t.getElectronX(); 
    p =  t.getElectronP(); 
    reco_ioniz_ion = (TH1F*) t.getIonizationIon(); 
    reco_ioniz_recoil = (TH1F*) t.getIonizationRecoils(); 


    tree->Fill(); 

  }


  tree->Write();


}
