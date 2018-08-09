#include "CollisionReader.hh" 
#include "TPolyLine3D.h" 
#include <assert.h>
#include "TPolyMarker3D.h" 
#include "PeriodicTable.hh"

#define BUF_SIZE 256 

const char magic[] = "^_^"; 

#define WR(x,f) (assert(fwrite(&(x),sizeof(x),1,f)==1))
#define RE(x,f) (assert(fread(&(x),sizeof(x),1,f)==1))

//UGLY HACK 
// replace extended ascii characters with mundane ones in order to avoid encoding issues in scanf 

static char * fixbuf(char * buf) 
{
  for (int i = 0; i < strlen(buf); i++) 
  {

    if (buf[i] == char(186)) buf[i] = '$'; 
    else if (buf[i] == char(179)) buf[i] = '|'; 
    else if (buf[i] == char(219)) buf[i] = '@'; 
  }
  return buf; 
}

double retrim::CollisionReader::getDisplacementEnergy(int Z) const
{
  for (unsigned i = 0; i < Zs.size(); i++)
  {
    if (Z == Zs[i]) 
      return displacementEnergies[i]; 
  }
  return 0; 
}

double retrim::CollisionReader::getLatticeEnergy(int Z) const
{
  for (unsigned i = 0; i < Zs.size(); i++)
  {
    if (Z == Zs[i]) 
      return latticeEnergies[i]; 
  }
  return 0; 
}

double retrim::CollisionReader::getSurfaceEnergy(int Z) const
{
  for (unsigned i = 0; i < Zs.size(); i++)
  {
    if (Z == Zs[i]) 
      return surfaceEnergies[i]; 
  }
  return 0; 
}

int checkBinary(const char * fname) 
{

  FILE * file = fopen(fname,"rb"); 

  if (!file)
  {
    fprintf(stderr, "Cannot open file %s\n",fname); 
    return -1; 
  }

  char buf[3]; 
  fread(buf,1,3,file); 

  fclose(file); 

  if (!memcmp(buf,magic,3)) 
  {
    return 1; 
  }
  else if  (!memcmp(buf,"===",3))
  {
    return 0; 
  }

  fprintf(stderr, "file %s neither ascii nor binary...\n",fname); 
  return -1; 

}


retrim::CollisionReader::CollisionReader(const char * fname) 
{

  marker_color = 2; 
  line_color = 4; 

  int binary = checkBinary(fname);

  if (binary < 0) return; 

  else if (binary) 
  {
    readBinary(fname); 
  }
  else
  {
    readAscii(fname); 
  }

}

int retrim::CollisionReader::readAscii(const char * fname)
{

  FILE* collfile = fopen(fname,"r"); 
  
  if (!collfile)
  {
    fprintf(stderr, "Cannot open file %s\n",fname); 
    return 1; 
  }

  char buf[BUF_SIZE]; 

  //get 9th line
  for (int i = 0; i < 9; i++) fgets(buf,BUF_SIZE, collfile); 

  // get ion energy 
  sscanf(fixbuf(buf),"$     Ion Energy =%f keV                          $",&E); 

  // get collision energies 
  while(buf[0]!='>')  fgets(buf,BUF_SIZE, collfile); 

while(buf[0] == '>') 
{
  char elem0[2]; 
  char elem1[2]; 
  float dispE; 
  float lattE; 
  float surfE; 
  sscanf(buf,">      %s  Displacement Energy of  %s = %f eV", elem0,elem1, &dispE); 
  fgets(buf,BUF_SIZE, collfile); 
  sscanf(buf,"  Latt.Binding Energy of  %s = %f eV", elem1, &lattE); 
  fgets(buf,BUF_SIZE, collfile); 
  sscanf(buf,"  SurfaceBind. Energy of  %s = %f eV", elem1, &surfE); 
  fgets(buf,BUF_SIZE, collfile); 

  int Z = PeriodicTable::getZ(elem0); 
//  printf("elem0: %s, elem1: %s, Z: %d\n", elem0, elem1, Z); 

  displacementEnergies.push_back(dispE*1e-3); 
  latticeEnergies.push_back(lattE*1e-3); 
  surfaceEnergies.push_back(surfE*1e-3); 
  Zs.push_back(Z); 
}

  //skip to first collision 

  std::vector<collision> track; 

  while(buf[0]!=char(179)) fgets(buf,BUF_SIZE, collfile); 

  
  //go until we run out of tracks
  while (true) 
  {
    collision c; 
 

    char atom_hit_buf[3]; 
    int nread = sscanf(fixbuf(buf),"|%d|%e|%e|%e|%e|%e| %02s |%e|",&(c.i), &(c.Ei), &(c.x), &(c.y), &(c.z), &(c.Se), atom_hit_buf, &(c.Er) );   

    c.atom_hit = PeriodicTable::getZ(atom_hit_buf); 

    //printf(buf); 

    //convert Er to keV 
    c.Er *= 1e-3; 

    //convert A to mm 
    c.x *= 1e-7; 
    c.y *= 1e-7; 
    c.z *= 1e-7; 

    //printf("read collision: ion=%d, E=%f, x=(%f,%f,%f), Er=%f, hit=%s\n",c.i,c.Ei,c.x,c.y,c.z,c.Er, c.atom_hit);

    bool formatting_bug = false; 
    //check for cascades
    if (strstr(buf,"Start of New Cascade"))
    {
      //skip 3 lines
      for (int i = 0; i< 3; i++) fgets(buf,BUF_SIZE, collfile); 
             
       while (buf[0]!='=' && buf[0] != char(179))
       {
         secondary_collision s; 

         sscanf(fixbuf(buf),"@ %hd %hd %e %e %e %e %hd %hd",&(s.recoil), &(s.atom), &(s.Er), &(s.x), &(s.y), &(s.z), &(s.vac), &(s.repl)); 
         //printf(buf); 

         //convert Er to keV 
         s.Er *= 1e-3; 

         //convert A to mm 
         s.x *= 1e-7; 
         s.y *= 1e-7; 
         s.z *= 1e-7; 

         //printf("read secondary collision: recoil=%d, E=%f, x=(%f,%f,%f),  hit=%d\n",s.recoil, s.Er, s.x,s.y,s.z,s.atom);

         if (!fgets(buf,BUF_SIZE, collfile)) break; 
         c.secondaries.push_back(s); 
       }
       //eat another line (summary line) if reached = 
       if (buf[0] == '=') fgets(buf,BUF_SIZE, collfile); 
       else formatting_bug = true; 
       track.push_back(c); 
    }
    else if (nread==8) track.push_back(c); 

    //reached eof 
    if (!formatting_bug && !fgets(buf,BUF_SIZE, collfile)) 
    {
      if (track.size() > 0) 
        tracks.push_back(track); 
      break; 
    }

    // check for summary 
    if (buf[0]=='=')
    {
      if (track.size())
        tracks.push_back(track); 
      track.clear(); 
      while(buf[1]!='-' && fgets(buf,BUF_SIZE, collfile)); 
    }

    
  }

  printf("Read %d tracks with E = %f\n",tracks.size(), E);  
  fclose(collfile); 
  return  0; 

}

int retrim::CollisionReader::drawTrack(int t, const char * opt) const
{
  // set the range of the view 
  
  unsigned npoints = nCollisions(t); 

  /*
  double max_x = tracks[t][npoints-1].x * 1.5; 

  max_x*=2; 
  double min_x = -0.5 * max_x; 
  double min_y = -0.5 * max_x; 
  double max_y = 0.5 * max_x; 
  double min_z = -0.5 * max_x; 
  double max_z = 0.5 * max_x; 
  
  view->SetRange(min_x,min_y,min_z,max_x,max_y,max_z); 
  view->SetSystem(1); 

  */


  if (strchr(opt,'m'))
  {
    TPolyMarker3D * line = new TPolyMarker3D(npoints); 
    for (unsigned i = 0; i < npoints; i++)
    {
      line->SetPoint(i,tracks[t][i].x,tracks[t][i].y,tracks[t][i].z); 
      drawCollision(t,i,"m"); 
    }

    line->Draw(); 
  }
  else if (strchr(opt,'l'))
  {
    TPolyLine3D * line = new TPolyLine3D(npoints); 
    for (unsigned i = 0; i < npoints; i++)
    {
      line->SetPoint(i,tracks[t][i].x,tracks[t][i].y,tracks[t][i].z); 
      drawCollision(t,i,"m"); 
    }

    line->SetLineWidth(2); 
    line->SetLineColor(line_color); 
    line->Draw(); 
  }

  return npoints; 

} 


int retrim::CollisionReader::drawCollision(int t, int i, const char * opt) const
{

    unsigned nsecondaries = tracks[t][i].secondaries.size(); 

    if (nsecondaries > 0)
    {
      if (strchr(opt,'m'))
      {

        for (unsigned j =0; j < nsecondaries; j++)
        {
          TPolyMarker3D * marker = new TPolyMarker3D(1); 
          marker->SetPoint(0, tracks[t][i].secondaries[j].x,  tracks[t][i].secondaries[j].y,  tracks[t][i].secondaries[j].z);  
          marker->SetMarkerColor(marker_color); 
          marker->SetMarkerStyle(7); 
          marker->SetName(TString::Format("collision_%d_%d_%d",t,i,j)); 
          //marker->SetMarkerSize(tracks[t][i].secondaries[j].Er/10); 
          marker->Draw(); 
        }

      }
      else
      {
        TPolyLine3D * marker = new TPolyLine3D(nsecondaries); 
        for (unsigned j =0; j < nsecondaries; j++)
        {
          marker->SetPoint(j, tracks[t][i].secondaries[j].x,  tracks[t][i].secondaries[j].y,  tracks[t][i].secondaries[j].z);  
        }

        marker->SetLineColor(2); 
        marker->Draw(); 
      }
    }

    return nsecondaries;
}



int retrim::CollisionReader::readBinary(const char * file) 
{
  FILE  * f = fopen(file,"rb"); 
  if (!f) 
  {
    fprintf(stderr,"Cannot open destination file %s\n",file); 
    return 1; 
  }

  char buf[3]; 

  fread(buf,1,3,f); 
  assert(!memcmp(buf,magic,3)); 

  RE(E,f); 
 // printf("E: %f\n",E); 

  unsigned nions; 
  RE(nions,f); 
//  printf("nions: %u\n",nions); 

  Zs.insert(Zs.begin(),nions,0); 
  displacementEnergies.insert(displacementEnergies.begin(),nions,0); 
  latticeEnergies.insert(latticeEnergies.begin(),nions,0); 
  surfaceEnergies.insert(surfaceEnergies.begin(),nions,0); 

  int *Z = &(Zs[0]); 
  float *d = &(displacementEnergies[0]); 
  float *l = &(latticeEnergies[0]); 
  float *s = &(surfaceEnergies[0]); 

  fread(Z, sizeof(*Z),nions,f);  
  fread(d, sizeof(*d),nions,f);  
  fread(l, sizeof(*l),nions,f);  
  fread(s, sizeof(*s),nions,f);  



  unsigned ntracks; 
  RE(ntracks,f); 

  tracks.reserve(ntracks); 
  for (unsigned i = 0; i < ntracks; i++)
  {
  
    std::vector<collision> colls; 
    unsigned ncoll; 
    RE(ncoll,f); 
    colls.reserve(ncoll); 
    

    for (unsigned j = 0; j < ncoll; j++)
    {
      collision c; 
      RE(c.x,f); 
      RE(c.y,f); 
      RE(c.z,f); 
      RE(c.Er,f); 
      RE(c.Ei,f); 
      RE(c.Se,f); 
      RE(c.atom_hit,f); 
      RE(c.i,f); 
 
      unsigned coll_nsec; 
      RE(coll_nsec, f); 

      c.secondaries.reserve(coll_nsec); 

      for (unsigned s = 0; s < coll_nsec; s++) 
      {
        secondary_collision sec; 
        RE(sec,f); 
        c.secondaries.push_back(sec); 
      }

      colls.push_back(c); 
    }
    tracks.push_back(colls); 
  }
  


  fclose(f); 
  printf("Read %d %f keV tracks\n",tracks.size(),E);  
}



int retrim::CollisionReader::writeBinary(const char * file)  const
{
  FILE  * f = fopen(file,"wb"); 
  if (!f) 
  {
    fprintf(stderr,"Cannot open destination file %s\n",file); 
    return 1; 
  }

  //write magic number
  fwrite(magic,3,1,f); 

  //write the energy 
  WR(E,f); 
  
  //write the number of atoms  
  unsigned nions = Zs.size(); 
  WR(nions,f); 


  //write the Zs, and energies for each ion
  
  const int *Z = &(Zs[0]); 
  const float *d = &(displacementEnergies[0]); 
  const float *l = &(latticeEnergies[0]); 
  const float *s = &(surfaceEnergies[0]); 


  fwrite(Z, sizeof(*Z),nions,f);  
  fwrite(d, sizeof(*d),nions,f);  
  fwrite(l, sizeof(*l),nions,f);  
  fwrite(s, sizeof(*s),nions,f);  

  //write the number of tracks
  
  unsigned ntracks = tracks.size(); 
  WR(ntracks,f);

  for (unsigned i = 0; i < ntracks; i++)
  {
    //write the number of collisions
    unsigned ncoll = tracks[i].size(); 
    WR(ncoll,f); 
    for (unsigned j = 0; j < ncoll; j++)
    {
      //write collision info 

      collision c = tracks[i][j]; 
      WR(c.x,f); 
      WR(c.y,f); 
      WR(c.z,f); 
      WR(c.Er,f); 
      WR(c.Ei,f); 
      WR(c.Se,f); 
      WR(c.atom_hit,f); 
      WR(c.i,f); 
      unsigned coll_nsec = tracks[i][j].secondaries.size(); 
      WR(coll_nsec,f); 

      for (unsigned s = 0; s < coll_nsec; s++)
      {
        secondary_collision sec = tracks[i][j].secondaries[s]; 
        WR(sec,f);  //pod, so we can do this
      }
    }
  }
  fclose(f); 
}




