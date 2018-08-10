#include "HistogramReader.hh"
#include <cstdio> 


#include "TH1.h" 

#define BUF_SIZE 256 

static int next(char * buf, FILE * f) 
{
    if (!fgets(buf,BUF_SIZE, f)) 
    {
      fprintf(stderr,"ERROR OR EOF!!!\n");
      return 1; 
    }

    return 0; 

}


static int skipN(int n, char * buf, FILE * f)
{
  for (int i = 0; i <n; i++)
  {
    if (next(buf,f))
    {
      return 1;  
    }
  }

  return 0; 
}

static int skipUntil(const char * needle, char * buf, FILE * f)
{
  while(!strstr(buf,needle))
  {
    if (next(buf,f)) return 1; 
  }

  return 0; 
}


int retrim::HistogramReader::readPHONON(const char * file, TH1F ** ion, TH1F** recoil)
{

   FILE * f = fopen(file,"r"); 

   if (!f) 
   {
      fprintf(stderr, "Cannot open file %s\n",file); 
      return 1; 
   }


   char buf[BUF_SIZE] = {0}; 

   skipUntil("DEPTH", buf, f); 
   skipN(3,buf,f); 

   //get binning 
   
   float depth, si, sr; 

   sscanf(buf,"%e %e %e", &depth, &si,&sr); 
//   printf("%s\n",buf); 
//   printf("depth: %f", depth); 


   depth *= 1e-7; //A -> mm 
   si*= 1e4; // ev/A -> keV/mm

   if (*ion) 
   {
     *ion = new (*ion) TH1F("phonon_ion","phonons from ions", 100, 0, 100 * depth); 
   }
   else
   {
     *ion = new TH1F("phonon_ion","phonons from ions", 100, 0, 100 * depth); 
   }

   if (*recoil) 
   {
     *recoil = new (*recoil) TH1F("phonon_recoil","phonons from recoils", 100, 0, 100 * depth); 
   }
   else
   {
     *recoil = new TH1F("phonon_recoil","phonons from recoils", 100, 0, 100 * depth); 
   }

   (*ion)->SetBinContent(1,si); 
   (*recoil)->SetBinContent(1,sr); 
   for (int i = 2; i <= 100; i++)
   {
     if (next(buf,f)) return 1 ; 
//     printf("%s\n",buf); 
     sscanf(buf,"%e %e %e", &depth, &si, &sr); 
     (*ion)->SetBinContent(i,si*1e4); 
     (*recoil)->SetBinContent(i,sr*1e4); 
   }
   fclose(f);
}



int retrim::HistogramReader::readIONIZ(const char * file, TH1F ** ion, TH1F** recoil)
{

   FILE * f = fopen(file,"r"); 

   if (!f) 
   {
      fprintf(stderr, "Cannot open file %s\n",file); 
      return 1; 
   }


   char buf[BUF_SIZE] = {0}; 

   skipUntil("DEPTH", buf, f); 
   skipN(3,buf,f); 

   //get binning 
   
   float depth, si, sr; 

   sscanf(buf,"%e %e %e", &depth, &si,&sr); 
//   printf("%s\n",buf); 
//   printf("depth: %f", depth); 


   depth *= 1e-7; //A -> mm 
   si*= 1e4; // ev/A -> keV/mm

   if (*ion) 
   {
     *ion = new (*ion) TH1F("ioniz_ion","ionization by ions", 100, 0, 100 * depth); 
   }
   else
   {
     *ion = new TH1F("ioniz_ion","ionization by ions", 100, 0, 100 * depth); 
   }

   if (*recoil) 
   {
     *recoil = new (*recoil) TH1F("ioniz_recoil","ionization by recoils", 100, 0, 100 * depth); 
   }
   else
   {
     *recoil = new TH1F("ioniz_recoil","ionization by recoils", 100, 0, 100 * depth); 
   }

   (*ion)->SetDirectory(0);
   (*recoil)->SetDirectory(0);
   (*ion)->SetBinContent(1,si); 
   (*recoil)->SetBinContent(1,sr); 
   for (int i = 2; i <= 100; i++)
   {
     if (next(buf,f)) return 1 ; 
//     printf("%s\n",buf); 
     sscanf(buf,"%e %e %e", &depth, &si, &sr); 
     (*ion)->SetBinContent(i,si*1e4); 
     (*recoil)->SetBinContent(i,sr*1e4); 
   }

   fclose(f);


}






