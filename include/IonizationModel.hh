#ifndef RETRIM_IONIZATION_MODEL_HH
#define RETRIM_IONIZATION_MODEL_HH


#include <vector>
#include "TLorentzVector.h"
#include "TRandom3.h"

namespace retrim
{

  class IonizationModel 
  {

    public: 
      /** Create electrons on a segment between start and end with energy  loss dE. The positions / times and energies /momenta of the electrons are added to the vectors. The number of electrons created is returned. 
       *
       *
       **/ 

      virtual int makeElectrons(size_t n, const TVector3 * startpoints, const TVector3* endpoints,
                                  const double *t0, const double *t1, const double * dE, const double *Eat, const int * species, 
                                  std::vector<TLorentzVector> * x, std::vector<TLorentzVector> * p) const = 0; 
      virtual ~IonizationModel() {;}

  }; 


  class SimpleIonizationModel : public IonizationModel 
  {
    /** This ionization model just uses the work function and fano factor to compute the number of electrons created. 
     */ 

    public: 
      
      SimpleIonizationModel(double W_in_eV=34 , double fano_factor=0.18 ); 
      virtual int makeElectrons(size_t n, const TVector3 * startpoints, const TVector3 * endpoints, const double *t0, 
                                const double *t1, const double *dE, const double * Eat, const int * species, 
                                std::vector<TLorentzVector> * x, std::vector<TLorentzVector> * p) const; 
      void setSeed(unsigned seed = 0) { rand.SetSeed(seed); } 
    private: 
      double W; 
      double fano; 
      mutable TRandom3 rand; 
  }; 


}


#endif
