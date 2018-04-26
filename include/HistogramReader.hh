#ifndef RETRIM_HISTOGRAM_READER_HH
#define RETRIM_HISTOGRAM_READER_HH

class TH1;
class TH1F;


namespace retrim
{

  /* functions to read TRIM output histograms 
   *
   */ 
  namespace HistogramReader
  {

    /* read IONIZ.TXT. 
     *
     *  if *(ion_ioniz) or *(recoil_ioniz), placement new will be used so they better be initialized if non-zero! 
     */
    
    int readIONIZ(const char * file, TH1F ** ion_ioniz, TH1F **  recoil_ioniz); 
    int readPHONON(const char * file, TH1F ** ion_ioniz, TH1F **  recoil_ioniz); 




  }
}

#endif
