#ifndef RETRIM_READ_SAV_3D_HH
#define RETRIM_READ_SAV_3D_HH
//Cosmin Deaconu
//cozzyd@gmail.com

class TH2D; 


/** This will parse TRIM's binary save format for 3D (really 2D) histograms.
 *
 *  This is useful because there isn't an option to save 3D (well, really 2D) histograms
 *  in batch mode, which makes it impossible to collect large numbers of them for individual ions. 
 *
 *  However, it is possible to tell TRIM to save after each ion, in which case it will create the 
 *  e.g. Ioniz-3d.sav which can be copied into another directory and renamed. 
 *
 *  This has been reverse-engineered specifically for the 2D ionization data, but should work for all forms of 
 *  the histograms.  
 *  
 */

namespace retrim
{
    TH2D  * readSav3D(const char * file, double binwidth_in_um=1, const char * name = "TRIMSav"); 
}



#endif
