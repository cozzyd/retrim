#include "EnergyLossModel.hh"
#include "TRandom3.h" 
#include <assert.h>
#include "TableReader.hh" 


#define MAX_Z 200 

retrim::ConstantEnergyLossModel::ConstantEnergyLossModel(int ntables, const TableReader ** tables, double ratio, double maxnuke)
  : Zs(MAX_Z), ratio (ratio), max_nuke(maxnuke), srim (tables)
{
  for (int i = 0; i < ntables; i++)
  {
    Zs[tables[i]->ionZ()] = i; 
  }

  estimate_delta = 1e-3; 
}

float retrim::ConstantEnergyLossModel::getIonMass(int species) const 
{
  return srim[Zs[species]]->getIonMass(); 
}

double retrim::ConstantEnergyLossModel::getRange(int species, double E) const 
{
  return srim[Zs[species]]->getRange(E); 
}


retrim::EnergyLoss retrim::ConstantEnergyLossModel::dE(int species, double Estart, double distance) const
{

  int index = Zs[species]; 
  EnergyLoss loss; 
  loss.elec = 0;
  loss.nuke = 0;

  if (Estart < 0) return loss;

  int nsegs = int(distance / estimate_delta + 1); 
  double delta = distance / nsegs; 

  double elecloss = 0; 
  double E = Estart; 

  double nukeloss = 0; 

  for (int i = 0; i < nsegs; i++)
  {
    double this_Eloss=  srim[index]->getSe(E) * delta * ratio; 
    double this_Nloss=  nuke_enabled ? srim[index]->getSn(E) * delta * ratio : 0; 

    if (nukeloss + this_Nloss > max_nuke) 
    {
      this_Nloss = max_nuke - nukeloss; 
    }
    nukeloss += this_Nloss; 
    elecloss +=  this_Eloss; 

    E= E - this_Eloss - this_Nloss; 
    if (E <= 0) break; 
  }

  if (nukeloss > max_nuke) 
  {
    nukeloss = max_nuke; 
  }

  loss.nuke = nukeloss; 
  loss.elec = elecloss; 
  assert(loss.nuke >=0); 
  assert(loss.elec >=0); 

  return loss;

}
