#include "IonizationModel.hh" 

#define ELECMASS 510.998910

#define AMU 931494 

static double getM(int species)
{
  switch (species)
  {
    case 1: 
      return AMU; 
    case 2: 
      return 4*AMU;
    case 6:
      return 12 *AMU; 
    case 9:
      return 19 *AMU; 

    default:
      return 0; 
  }

  return 0; 

}

retrim::SimpleIonizationModel::SimpleIonizationModel(double Wev, double fano_factor) 
{

  W = Wev * 1e-3; //eV -> keV 
  fano = fano_factor; 
  setSeed(1337); 
}

int retrim::SimpleIonizationModel::makeElectrons(size_t n, const TVector3 * start, const TVector3 * end, const double * t0, const double * t1, const  double *  dE,
                                                            const double *Eat, const int * species, std::vector<TLorentzVector> * x, std::vector<TLorentzVector> * p) const
{




  //build up cumsum 
  std::vector<double> cumsum; 

  cumsum.reserve(n); 
  double sum = 0; 
  for (size_t j = 0; j < n; j++)
  {
    sum += dE[j];
    cumsum.push_back(sum); 
  }


  double nEraw = rand.Gaus(sum / W , sqrt(sum * fano/W)); 
  int nE = int(nEraw) + (rand.Uniform() < (nEraw - int(nEraw))); 



  

 for (int i = 0; i < nE; i++)
 {
   
    double picker = rand.Uniform() * sum; 
    size_t j = 0; 

    while (cumsum[j] < picker) 
    {
      j++;
    }

    double M = getM(species[j]); 
    double Emax = 4 *Eat[j]* (M * ELECMASS) / pow(M+ELECMASS,2); 


    TVector3 vec = end[j]-start[j]; 
    double distance = vec.Mag(); 
    vec = vec.Unit(); 

    double t_total = t1[j] - t0[j]; 



    double position = rand.Uniform(); 
    double pmag = pow(Emax , rand.Uniform()); 
    double E = sqrt(pmag*pmag + ELECMASS*ELECMASS); 
    double px,py,pz; 
    rand.Sphere(px,py,pz,pmag); 


    x->push_back(TLorentzVector(start[j] + vec * position * distance, t0[j] + position * t_total)); 
    p->push_back(TLorentzVector(TVector3(px,py,pz), E)); 


 }

 return nE; 
}
