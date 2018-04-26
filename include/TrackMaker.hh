#ifndef DMTPC_MC_RETRIM_TRACK_MAKER_HH
#define DMTPC_MC_RETRIM_TRACK_MAKER_HH

#include "TVector3.h"
#include "TLorentzVector.h"
#include "TRandom3.h" 
#include <vector> 


class TH1; 
class TTree; 
namespace retrim
{
  class CollisionReader; 
  class IonizationModel; 
  class TableReader;

  class TrackMaker
  {

    public:

      /* you need to pass srim tables for all possible recoils as well. first srim table passed assumed to be that of primary. pressure ratio is ratio between target pressure and SRIM table pressure. */ 
      TrackMaker(const CollisionReader * trim_coll, int nsrim, const TableReader ** srim_tables, double pressure_ratio = 1); 
      virtual ~TrackMaker(); 
      void setSRIM(int nsrimtables, const TableReader ** srim_tables, double pressure_ratio = 1)
        {srim = srim_tables; nsrim = nsrimtables; ratio=pressure_ratio; }
      void setTRIM(const CollisionReader * trim_coll) { trim=trim_coll; }
      int makeTrack(int tracknum, const TVector3 * origin = 0, const TVector3 * direction = 0, double t0 = 0); 
      int enableHistograms(double binsize_in_mm, int nbins=100); 
      int ntracks() const; 
      double E() const; 
      void clearHistograms(); 
      TH1 * getIonizationIon() const { return ioniz_ion; }
      TH1 * getIonizationRecoils() const { return ioniz_recoils; }
      TH1 * getPhononsIon() const { return phonon_ion; }
      TH1 * getPhononsRecoils() const { return phonon_recoils; }
      TTree * getStats() { return stats; } 



      /* s = segments, e = electrons. Returns vector of allocated objects so they may be properly deleted. */ 
      std::vector<TObject*> draw(const char * opt = "se") const; 

      /* Fill the histogram with the number of electrons in each bin. */
      void fillElectrons(TH1 * hist) const; 
      struct particle;  

      unsigned getNSegments() const { return startpoints.size(); } 
      const TVector3 * getSegmentStart(int i) const { return & (startpoints[i]); }
      const TVector3 * getSegmentEnd(int i) const { return & (endpoints[i]); }
      double getSegmentIoniz(int i) const  { return  ELoss[i]; }
      double getSegmentNuke(int i) const  { return  nukeLoss[i]; }
      int getSegmentSpecies(int i ) const  { return species[i]; } 
      int getSegmentParticleNumber(int i ) const  { return particle_number[i]; } 

      void setVerbose(bool b) {verbose = b;}
      int makeElectrons(const IonizationModel * ioniz); 
      /* Get the ionization in the hypercube defined by the input coordinates. If any of the bounds equal each other, that dimension will be integrated over */ 
      double getIonization(double x0, double x1, double y0=0, double y1=0, double z0 = 0, double z1 = 0) const; 
      double getTotalIonization() const { return ionization; }
      int getNElectrons() const { return electron_x.size(); } 
      const std::vector<TLorentzVector> * getElectronX() const { return &electron_x; } 
      const std::vector<TLorentzVector> * getElectronP() const { return &electron_p; } 
      const TLorentzVector * getElectronX(int i) const { return & ( electron_x[i]); }
      const TLorentzVector * getElectronP(int i) const { return & ( electron_p[i]); }
      const std::vector<double> * getLeftover() const {return &leftover; } 
      double getMaxLeftover() const { return maxleftover; } 
      double getMinLeftover() const { return maxleftover; } 
      void setLeftoverCut(double min, double max) {  cut_maxleftover = max; cut_minleftover = min;} 
      void enableLeftoverCut(bool enable) { leftover_cut = enable; } 
      void enableDaughterRandomization(bool enable) { daughter_randomization = true; } 
      void setMaxIter(int n) {max_iter = n; }
      int getElectron(int i, double & x, double & y, double & z, double & t, double & E, double & px, double & py, double & pz); 
      int nHistogram() const { return nhist; } 



    private:
        int fillTree(particle *p, bool record); 

        int clearTrack(); 
        void transform(const TVector3 * origin, const TVector3 * dir); 
        void fillIonizHists(); 
        std::vector<TVector3> startpoints; 
        std::vector<TVector3> endpoints; 
        std::vector<double> starttimes; 
        std::vector<double> endtimes; // BEWARE THE RAPTURE
        std::vector<double> ELoss; 
        std::vector<double> nukeLoss; 
        std::vector<double> Eat; 
        std::vector<int> particle_number; 
        std::vector<int> species; 
        std::vector<bool> randomized; 

        std::vector<TLorentzVector> electron_x; 
        std::vector<TLorentzVector> electron_p; 

        std::vector<double> leftover; 
        int nhist; 

        TTree * stats; 
        double ElecLoss; 
        double predictedElecLoss; 
        double minleftover; 
        double maxleftover; 
        double cut_maxleftover; 
        double cut_minleftover; 
        double score; 
        int max_iter; 
        double ratio; 
        const TableReader ** srim; 
        const CollisionReader * trim; 
        int nsrim; 
        int particleCounter; 
        TRandom3 rng; 
        TH1 * ioniz_ion; 
        TH1 * ioniz_recoils; 
        TH1 * phonon_ion; 
        TH1 * phonon_recoils; 
        TH1 * e2r_ion;
        TH1 * e2r_recoils;
        bool verbose; 
        bool daughter_randomization; 
        bool leftover_cut; 
        double ionization; 

  }; 
}

#endif
