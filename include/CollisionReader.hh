#ifndef RETRIM_COLLISON_READER_HH
#define RETRIM_COLLISON_READER_HH

#include <vector>

class TTree; 

namespace retrim
{

  /* Collision reader for ion and recoil collisions!
   *
   * use IonCollisionReader for ion only */ 
  class CollisionReader
  {
    public: 

    struct secondary_collision
    {
      float x,y,z,Er; 
      short vac, repl,recoil,atom; 
    }; 

    struct collision
    {
      float x,y,z; 
      float Er; 
      float Ei; 
      float Se; 
      int atom_hit; 
      int i; 
      std::vector<secondary_collision> secondaries; 
    }; 

      CollisionReader (const char * collisions_file); 
      unsigned nTracks() const { return tracks.size(); } 
      unsigned nCollisions(int track) const { return tracks[track].size(); } 
      unsigned nSecondaries(int track, int coll) { return tracks[track][coll].secondaries.size(); } 
      const collision * getCollision(int track, int coll) const { return &(tracks[track][coll]); } 
      const secondary_collision * getSecondaryCollision(int track, int coll, int sec) const { return &(tracks[track][coll].secondaries[sec]); } 
      int drawTrack(int track, const char * option="l") const; 
      int drawCollision(int track, int coll, const char * option="m") const; 
      float getE() const { return E; } 
      void setLineColor(int c) { line_color = c;}
      void setMarkerColor(int c) { marker_color = c;}
      int writeBinary(const char * out) const; 

      double getDisplacementEnergy(int Z) const;
      double getLatticeEnergy(int Z) const;
      double getSurfaceEnergy(int Z) const;

    private:
      int readBinary(const char * collisions_file); 
      int readAscii(const char * collisions_file); 


      std::vector<std::vector<collision> > tracks;
      std::vector<float> displacementEnergies; 
      std::vector<float> latticeEnergies; 
      std::vector<float> surfaceEnergies; 
      std::vector<int> Zs; 
      float E; 
      int marker_color; 
      int line_color; 
  }; 
}


#endif
