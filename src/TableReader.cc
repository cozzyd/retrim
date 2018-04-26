#include "TableReader.hh" 
#include <cstdio> 
#include "TTree.h"
#include <cstring>
#include "TFile.h" 
#include "TAxis.h" 

#define BUF_SIZE 256 

#define DEBUG 0 

retrim::TableReader::TableReader(const char * fname) 
{

  FILE* tablefile = fopen(fname,"r"); 

  if (!tablefile) 
  {
    fprintf(stderr, "Cannot open file %s\n",fname); 
    return; 
  }

  char buf[BUF_SIZE] = {0}; 

  //skip until ion info  lines
  //
  while(!strstr(buf,"Ion =")) fgets(buf,BUF_SIZE, tablefile); 

  //Get Ion info 
  if (DEBUG) printf(buf); 
  sscanf(buf," Ion = %80s [%d] , Mass = %f amu", ion_name, &ion_Z, &ion_mass); 
  
  //skip line
  fgets(buf,BUF_SIZE,tablefile); 
  if (DEBUG) printf(buf); 

  //get density 
  fgets(buf,BUF_SIZE,tablefile); 
  if (DEBUG) printf(buf); 
  sscanf(buf," Target Density = %e g/cm3 = %e atoms/cm3", &density, &number_density); 


  //get state
  fgets(buf,BUF_SIZE,tablefile); 
  if (strstr(buf,"GAS")) 
  {
     gas = true; 
     fgets(buf,BUF_SIZE,tablefile);  // read an extra line
  }
  else gas = false; 

  //skip 4 lines 
  for (int i = 0; i < 4; i++) fgets(buf,BUF_SIZE,tablefile); 

  //get composition
  while(buf[1] != '=')
  {
    if (DEBUG) printf(buf); 
    char atom[3]; 
    int Z; 
    float atomic_pct; 
    float mass_pct; 
    sscanf(buf, "%3s %d %f %f",atom,&Z,&atomic_pct,&mass_pct);  
    composition_atom.push_back(strdup(atom)); 
    composition_Z.push_back(Z); 
    composition_atomic_pct.push_back(atomic_pct); 
    composition_mass_pct.push_back(mass_pct); 
    fgets(buf,BUF_SIZE,tablefile); 
  }

  //get Bragg Correction 
  fgets(buf,BUF_SIZE,tablefile); 
  sscanf(buf," Bragg Correction = %f\%",&bragg_correction); 


  //skip 6 lines
  for (int i = 0; i < 6; i++) fgets(buf,BUF_SIZE,tablefile); 

  float E=0, Se=0, Sn=0, r=0, lgs=0, lts=0; 

  tree = new TTree("srim","SRIM Table"); 
  tree->SetDirectory(0); 

  tree->Branch("E",&E); 
  tree->Branch("Se",&Se); 
  tree->Branch("Sn",&Sn); 
  tree->Branch("range",&r); 
  tree->Branch("long_straggle",&lgs); 
  tree->Branch("lat_straggle",&lts); 

  //fill all zeros to help with future interpolation
  tree->Fill(); 


  //get conversions first... so we have to skip the table then come back to it
  long start_of_table = ftell(tablefile); 

  fgets(buf,BUF_SIZE,tablefile); 
  while (buf[0]!='-')
  {
    fgets(buf,BUF_SIZE,tablefile); 
  }
  
  while (!strstr(buf, "keV / micron"))
    fgets(buf,BUF_SIZE,tablefile); 

  float conv; 
  sscanf(buf,"%e",&conv); 

  // micron ->mm 
  conv *= 1e3; 

  //go back to start of table 
  fseek(tablefile,start_of_table,SEEK_SET); 

  //read table
  fgets(buf,BUF_SIZE,tablefile); 
  while (buf[0]!='-')
  {
    char E_unit[6]; 
    char r_unit[6]; 
    char lgs_unit[6]; 
    char lts_unit[6]; 

    sscanf(buf, "%f %5s %e %e %f %5s %f %5s %f %5s", &E, E_unit, &Se, &Sn, &r, &r_unit, &lgs, lgs_unit, &lts, lts_unit); 
    

    //do conversions to make Energy in keV, ranges in mm

    if (strcmp(E_unit,"keV"))
    {
      E *=  (!strcmp(E_unit,"eV")  ? 1e-3 : 
             !strcmp(E_unit,"MeV") ? 1e3 :
             !strcmp(E_unit,"GeV") ? 1e6 :
             0 //unrecognized unit
             ); 
    }

    if (strcmp(r_unit,"mm"))
    {
      r*=  (!strcmp(r_unit,"um") ? 1e-3 : 
            !strcmp(r_unit,"m") ? 1e3 :
            !strcmp(r_unit,"A") ? 1e-7 :
            0 //unrecognized unit
           ); 
    }
    
    if (strcmp(lgs_unit,"mm"))
    {
      lgs *=  (!strcmp(lgs_unit,"um") ? 1e-3 : 
               !strcmp(r_unit,"A") ? 1e-7 :
               0 //unrecognized unit
              ); 
    }
 
    if (strcmp(lts_unit,"mm"))
    {
      lts *=  (!strcmp(lts_unit,"um") ? 1e-3 : 
               !strcmp(lts_unit,"m") ? 1e3 :
               !strcmp(r_unit,"A") ? 1e-7 :
               0 //unrecognized unit
              ); 
    }
 
    //convert stopping power to keV / mm 
    Se *=conv; 
    Sn *=conv; 

    tree->Fill(); 
    fgets(buf,BUF_SIZE,tablefile); 
  }


  //done 
  fclose(tablefile); 

  //make graphs
  Se_energy = new TGraph(tree->Draw("E:Se","","goff"), tree->GetV1(), tree->GetV2()); 
  Se_energy->GetXaxis()->SetTitle("Ion Energy (keV)"); 
  Se_energy->GetYaxis()->SetTitle("Electronic Stopping Power (keV/mm)"); 
  Sn_energy = new TGraph(tree->Draw("E:Sn","","goff"), tree->GetV1(), tree->GetV2()); 
  Sn_energy->GetXaxis()->SetTitle("Ion Energy (keV)"); 
  Sn_energy->GetYaxis()->SetTitle("Nuclear Stopping Power (keV/mm)"); 
  range_energy = new TGraph(tree->Draw("E:range","","goff"), tree->GetV1(), tree->GetV2()); 
  range_energy->GetXaxis()->SetTitle("Ion Energy (keV)"); 
  range_energy->GetYaxis()->SetTitle("Projected Range (mm)"); 
  longStraggle_energy = new TGraph(tree->Draw("E:long_straggle","","goff"), tree->GetV1(), tree->GetV2()); 
  longStraggle_energy->GetXaxis()->SetTitle("Ion Energy (keV)"); 
  longStraggle_energy->GetYaxis()->SetTitle("Longitudinal Straggling (mm)"); 
  latStraggle_energy = new TGraph(tree->Draw("E:lat_straggle","","goff"), tree->GetV1(), tree->GetV2()); 
  latStraggle_energy->GetXaxis()->SetTitle("Ion Energy (keV)"); 
  latStraggle_energy->GetYaxis()->SetTitle("Latitudinal Straggling (mm)"); 
} 


retrim::TableReader::~TableReader()
{

  delete tree; 
  delete Se_energy; 
  delete Sn_energy; 
  delete range_energy; 
  delete longStraggle_energy; 
  delete latStraggle_energy; 
}



void retrim::TableReader::writeROOTFile(const char * fname, const char *option) const
{
  TFile f(fname, option); 
  TTree * copy = tree->CloneTree(); 
  copy->Write(); 
  f.Close(); 
}




