#ifndef RETRIM_ENERGY_LOSS_MODEL_HH
#define RETRIM_ENERGY_LOSS_MODEL_HH
class TRandom; 
#include <vector>

namespace retrim
{
  class TableReader; 


  struct EnergyLoss
  {
    float elec; 
    float nuke; 
  };

  class EnergyLossModel
  {

    public: 
      virtual EnergyLoss dE(int species, double startE, double distance) const = 0; 
      virtual void setNukeEnabled(bool set) { nuke_enabled = set; } 
      virtual double getRange(int species, double E) const = 0; 
      virtual float getIonMass(int species) const = 0; 
    protected: 
        bool nuke_enabled; 
  };

  class ConstantEnergyLossModel : public EnergyLossModel 
  {
    public: 

      ConstantEnergyLossModel(int ntables, const TableReader ** tables, double ratio = 1, double max_nuke = 0.025);
      EnergyLoss dE(int species, double startE, double distance) const; 
      void setEstimateDelta(double delta) { estimate_delta = delta; }
      virtual float getIonMass(int species) const; 
      virtual double getRange(int species, double E) const; 

    private: 
      double ratio; 
      const TableReader ** srim; 
      std::vector<int> Zs; 
      double max_nuke; 
      double estimate_delta; 

  };

  /*
  class RandomizingEnergyLossModel : public EnergyLossModel 
  {
    public: 

      RandomizingEnergyLossModel(int ntables, TableReader const ** tables, double ratio);
      EnergyLoss dE(int species, double startE, double distance) const; 


    private: 
      double ratio; 
      const TableReader ** tables; 
      std::vector<int> Zs; 
      mutable TRandom* rand; 

  };
  */

}


#endif
