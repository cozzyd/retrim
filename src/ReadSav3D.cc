#include "ReadSav3D.hh" 
#include "TH2.h" 
#include <cstdio>


TH2D * retrim::readSav3D(const char * file, double binwidth,  const char * name) 
{
  FILE * fptr = fopen(file,"rb") ; 
  if (!fptr) 
  {
    fprintf(stderr,"Cannot open %s\n", file); 
    return 0; 
  }

  TH2D * hist = new TH2D(name,TString::Format("Parsed %s",file), 100,0,100*binwidth, 100,-50*binwidth, 50*binwidth);  

  unsigned short byte; 
  double val; 

  int ibinx = 1; 
  int ibiny = 1; 
  while (!feof(fptr)) 
  {
    fread(&byte,2,1,fptr); 
    if (byte==5) 
    {
      fread(&val,8,1,fptr); 
      hist->SetBinContent(ibinx,ibiny,val); 
    }


    ibinx++; 
    if (ibinx > 100) 
    {
      ibiny++; 
      ibinx = 1; 
    }

    if (ibiny > 100) 
      break; 
  }


  hist->Scale(binwidth*10); //SRIM Units are in eV/binwidth in angstrom, so this gives units of keV in output
  hist->GetXaxis()->SetTitle("depth (um)"); 
  hist->GetYaxis()->SetTitle("lateral (um)"); 
  hist->GetZaxis()->SetTitle("ioniz (keV)"); 

  fclose(fptr);

  return hist; 
}





