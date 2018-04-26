#ifndef RETRIM_TABLE_READER__HH
#define RETRIM_TABLE_READER__HH

class TTree; 
#include "TGraph.h" 
#include <vector> 

/* SRIM table reader */ 


namespace retrim
{
  class TableReader
  {
    public: 
      TableReader(const char * fname); 
      virtual ~TableReader(); 

      void writeROOTFile(const char * fname, const char * option = "RECREATE") const; 
      double getSe(double E) const { return Se_energy->Eval(E); } 
      double getSn(double E) const  { return Sn_energy->Eval(E); } 
      double getRange(double E) const { return range_energy->Eval(E); } 
      double getLongStraggle(double E) const { return longStraggle_energy->Eval(E); } 
      double getLatStraggle(double E) const { return latStraggle_energy->Eval(E); } 

      TTree* getTree(){ return tree; } 
      TGraph* getSeGraph(){ return Se_energy; } 
      TGraph* getSnGraph(){ return Sn_energy; } 
      TGraph* getRangeGraph(){ return range_energy; } 
      TGraph* getLongStraggleGraph(){ return longStraggle_energy; } 
      TGraph* getLatStraggleGraph(){ return latStraggle_energy; } 

      const char * getIonName() const { return ion_name; } 
      int ionZ() const { return ion_Z; } 

      // amu 
      float getIonMass() const {return ion_mass; } 

      // g/cm^3
      float getDensity() const { return density; } 

      // cm^-3
      float getNumberDensity() const { return number_density; } 

      bool isGas() const { return gas; } 

      float braggCorrection() const { return bragg_correction; } 


      unsigned getTargetNElems() const { return composition_atom.size(); } 
      const char * getTargetAtom(int i ) const { return composition_atom[i]; } 
      int getTargetZ(int i ) const { return composition_Z[i]; } 
      float getTargetAtomicPercentage(int i ) const { return composition_atomic_pct[i]; } 
      float getTargetMassPercentage(int i ) const { return composition_mass_pct[i]; } 



    private:
      TTree * tree; 
      TGraph *Se_energy; 
      TGraph *Sn_energy; 
      TGraph *range_energy; 
      TGraph *longStraggle_energy; 
      TGraph *latStraggle_energy; 

      char ion_name[80]; 
      int ion_Z; 
      float ion_mass; //amu
      float density; //g / cm^3
      float number_density; // cm^-3
      bool gas; 
      
      std::vector<char*> composition_atom; 
      std::vector<int> composition_Z; 
      std::vector<float> composition_atomic_pct; 
      std::vector<float> composition_mass_pct; 
      float bragg_correction; 


  }; 
}

#endif


