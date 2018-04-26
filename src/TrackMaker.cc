#include "TrackMaker.hh"
#include "CollisionReader.hh" 
#include "TPolyLine3D.h"
#include "TPolyMarker3D.h" 
#include <list>
#include "IonizationModel.hh"
#include "TableReader.hh" 
#include "TH1.h" 
#include <set>
#include <assert.h> 
#include "TTree.h"

#include <cstdio>

/** THERE BE DRAGONS HERE !!!!
 *
 * */ 

const double estimate_delta = 1e-3; //10 microns 
const double elec_fudge = 0.9; 
const double nuke_fudge = 0.3; 

struct retrim::TrackMaker::particle 
{
           particle * parent; 
           std::list<particle*> daughters; 
           double startE; 
           TVector3 start; 
           TVector3 startdir; 
           short species; 
           short index; 
           double startT; 
}; 


//TODO, make relativistic so more applicable for higher energies
static double getT(double T, double dist, double mass_in_amu)
{
  double m = 931494.061*mass_in_amu; //amu -> keV
  double v =  sqrt(2*fabs(T)/m) * 300;  //c = 300 mm/ns
  if (v == 0) return 0; 
  return dist / v; 
}


static void printTree(FILE * f, retrim::TrackMaker::particle * p, int offset)
{

  for (int i = 0; i < offset;i++) 
    fprintf(f,"\t"); 

  fprintf(f,"Particle i=%d at (%f,%f,%f) with E=%f\n", p->index, p->start.x(), p->start.y(), p->start.z(), p->startE); 

  for(std::list<retrim::TrackMaker::particle*>::iterator it = p->daughters.begin(); it!=p->daughters.end(); it++)
  {
    printTree(f,*it, offset+1); 
  }
}




/*** TODO: Do a better job here... The assumptions here are very sketchy. */ 
const bool nuke_enabled = true;
const double max_nuke_abs = 0.028; 
const double max_nuke_frac = 1.0; 
const double nuke_suppress_factor = 0.8; 
static double estimateElectronicLoss(const retrim::TableReader * srim, double Estart, double distance, double ratio, double * nuke = 0) 
{
  if (Estart < 0) return 0; 
  int nsegs = int(distance / estimate_delta + 1); 
  double delta = distance / nsegs; 

  double loss = 0; 
  double E = Estart; 

  double nukeloss = 0; 
  double max_nuke  = TMath::Min(max_nuke_abs, max_nuke_frac * Estart); 

  for (int i = 0; i < nsegs; i++)
  {
    double this_Eloss=  srim->getSe(E) * delta * ratio; 
    double this_Nloss=  nuke_enabled ? srim->getSn(E) * delta * ratio * nuke_suppress_factor : 0; 

    if (nukeloss + this_Nloss > max_nuke) 
    {
      this_Nloss = max_nuke - nukeloss; 
    }

    double frac =1;  
    if (E < this_Nloss + this_Eloss)  frac = E / (this_Nloss + this_Eloss); 

    nukeloss += this_Nloss * frac; 
    loss +=  this_Eloss * frac; 

    E= E - this_Eloss - this_Nloss; 
    if (E <= 0)
    {
      break; 
    }
  }

  if (nukeloss > max_nuke) 
  {
    nukeloss = max_nuke; 
  }
  if (nuke)
  {
    *nuke += nukeloss;
  }

  return loss; 
} 

static int whichSRIM(int atom, int nsrim, const retrim::TableReader **srim) 
{

  for ( int i = 0; i < nsrim; i++)
  {
//    printf("atom=%d, ionZ=%d\n", atom, srim[i]->ionZ()); 
    if (atom == srim[i]->ionZ()) 
    {
      return i; 
    }
  }
  fprintf(stderr,"Uhoh, could not find %d\n",atom); 
  return -1; 
}

int retrim::TrackMaker::fillTree(particle * p, bool record)
{

  if (record) startpoints.push_back(p->start); 
  double t = p->startT; 
  if (record) starttimes.push_back(t); 
  TVector3 v = p->start; 
  TVector3 d = p->startdir; 
  int ip = ++particleCounter; 
  double E = p->startE - trim->getLatticeEnergy(p->species);
  int srim_index = whichSRIM(p->species,nsrim,srim); 

  double nuke = 0;
  for(std::list<particle*>::iterator it = p->daughters.begin(); it!=p->daughters.end(); it++)
  {
    particle * q = *it; 
    double dist = (q->start-v).Mag(); 
    q->startT = t + getT(E,dist, srim[srim_index]->getIonMass()); 
    t = q->startT; 
    d= q->start-v; 
    v = q->start; 
    if (record) endpoints.push_back(q->start); 
    if (record) endtimes.push_back(t); 
    double El = estimateElectronicLoss(srim[srim_index], E, dist, ratio,&nuke); 
    if (record) ELoss.push_back(El); 
    if (record) nukeLoss.push_back(nuke); 
    nuke = 0; 
    if(record) randomized.push_back(false); 
    if(record) species.push_back(p->species); 
    if (record) Eat.push_back(E); 
    E=E - El-q->startE; 
    if (record) particle_number.push_back(ip); 
    fillTree(q,record); 
    if (record) startpoints.push_back(q->start); 
    if (record) starttimes.push_back(t); 
  }

  
  if(record) Eat.push_back(E); 
  if(record) leftover.push_back(E); 

  if (E > maxleftover) 
  {
    maxleftover = E; 
  }
  if (E < minleftover) 
  {
    minleftover = E; 
  }

  score += fabs(E); 

  if (E < cut_minleftover) 
  {
    if (verbose) 
    {
      printf("i: %d E: %f\n",p->index, E); 
      printf("^---WARNING!! < minleftover\n"); 
    }
    score += fabs(E); 
  }

  if (E > cut_maxleftover) 
  {
    if (verbose) 
    {
      printf("i: %d E: %f\n",p->index, E); 
      printf("^---WARNING!! > maxleftover\n"); 
    }
    score += fabs(E); 
    if (leftover_cut) 
      E = cut_maxleftover; 
  }

  double range = 0; 
  if (E < trim->getDisplacementEnergy(p->species)) 
  {
    nuke+= trim->getDisplacementEnergy(p->species); 
    if(record) ELoss.push_back(0); 
    if(record) randomized.push_back(false); 
  }
  else 
  {
    range = srim[srim_index]->getRange(E); 
    if(record) randomized.push_back(true); 
    if(record) ELoss.push_back(estimateElectronicLoss(srim[srim_index],E,range,ratio,&nuke)); 
  }
  if(record) nukeLoss.push_back(nuke); 

  if(record) species.push_back(p->species); 
  if(record) particle_number.push_back(ip); 

  if(record) endpoints.push_back(v + range * d.Unit() );  
  if(record) endtimes.push_back(t + getT(E,range, srim[srim_index]->getIonMass())); 
}


int retrim::TrackMaker::enableHistograms(double binsize_in_mm, int nbins)
{

  if (ioniz_ion) delete ioniz_ion; 
  ioniz_ion = new TH1F("ioniz_ion", "Ionization from Ion",nbins, 0, nbins * binsize_in_mm); 
  ioniz_ion->GetXaxis()->SetTitle("depth (mm)"); 
  ioniz_ion->GetYaxis()->SetTitle("ionization (keV/mm)"); 
  ioniz_ion->SetLineColor(4); 
  ioniz_ion->SetDirectory(0); 

  if (ioniz_recoils) delete ioniz_recoils; 
  ioniz_recoils = new TH1F("ioniz_recoils", "Ionization from Recoils",nbins, 0, nbins * binsize_in_mm); 
  ioniz_recoils->GetXaxis()->SetTitle("depth (mm)"); 
  ioniz_recoils->GetYaxis()->SetTitle("ionization (keV/mm)"); 
  ioniz_recoils->SetLineColor(2); 
  ioniz_recoils->SetDirectory(0); 

  if (phonon_ion) delete phonon_ion; 
  phonon_ion = new TH1F("phonon_ion", "Phononsfrom Ion",nbins, 0, nbins * binsize_in_mm); 
  phonon_ion->GetXaxis()->SetTitle("depth (mm)"); 
  phonon_ion->GetYaxis()->SetTitle("phononsation (keV/mm)"); 
  phonon_ion->SetLineColor(4); 
  phonon_ion->SetDirectory(0); 

  if (phonon_recoils) delete phonon_recoils; 
  phonon_recoils = new TH1F("phonon_recoils", "Fhononsfrom Recoils",nbins, 0, nbins * binsize_in_mm); 
  phonon_recoils->GetXaxis()->SetTitle("depth (mm)"); 
  phonon_recoils->GetYaxis()->SetTitle("phononsation (keV/mm)"); 
  phonon_recoils->SetLineColor(2); 
  phonon_recoils->SetDirectory(0); 



  if (e2r_ion) delete e2r_ion; 
  if (e2r_recoils) delete e2r_recoils; 

}


void retrim::TrackMaker::clearHistograms() 
{
  if (ioniz_ion) ioniz_ion->Reset(); 
  if (ioniz_recoils) ioniz_recoils->Reset(); 
  if (phonon_ion) phonon_ion->Reset(); 
  if (phonon_recoils) phonon_recoils->Reset(); 
  nhist = 0; 
}

retrim::TrackMaker::TrackMaker(const CollisionReader * trim_coll, int nsrim, const TableReader ** srim_tables, double pressure_ratio ) 
{
  setTRIM(trim_coll); 
  setSRIM(nsrim, srim_tables, pressure_ratio); 

  stats = new TTree("trim_stats","TRIM Stast"); 
  stats->SetDirectory(0); 
  stats->Branch("ElecLoss",&ElecLoss); 
  stats->Branch("predictedElecLoss",&predictedElecLoss); 

  ionization = 0; 
  ioniz_ion = 0; 
  ioniz_recoils = 0; 
  phonon_ion = 0; 
  phonon_recoils = 0; 
  e2r_ion = 0; 
  e2r_recoils = 0; 
  verbose = false; 
  clearTrack(); 
  cut_maxleftover = 0.2; 
  cut_minleftover = -0.1; 
  max_iter = 100; 
  leftover_cut = true; 
  daughter_randomization = true; 
}

retrim::TrackMaker::~TrackMaker()
{
  delete stats; 
  if (ioniz_ion) ioniz_ion->Delete(); 
  if (ioniz_recoils)  ioniz_recoils->Delete(); 
  if (e2r_ion) delete e2r_ion; 
  if (e2r_recoils) delete e2r_recoils; 
}

bool between(double x, double x0, double x1) 
{
  return (x > x0 && x < x1) || (x > x1 && x < x0); 
}

double retrim::TrackMaker::getIonization(double x0, double x1, double y0, double y1, double z0, double z1) const
{
  double total = 0; 
  for (unsigned i = 0; i < startpoints.size(); i++)
  {

    bool sx_between = x0 == x1 || between(startpoints[i].x(),x0,x1);
    bool ex_between = x0 == x1 || between(endpoints[i].x(),x0,x1); 

    bool sy_between = y0 == y1 || between(startpoints[i].y(),y0,y1);
    bool ey_between = y0 == y1 || between(endpoints[i].y(),y0,y1); 

    bool sz_between = z0 == z1 || between(startpoints[i].z(),z0,z1);
    bool ez_between = z0 == z1 || between(endpoints[i].z(),z0,z1); 


    //all coordinates between bounds means full energy loss
    if (sx_between && ex_between && sy_between && ey_between && sz_between && ez_between)
    {
      total += ELoss[i]; 
    }

    //otherwise we must compute intersection of segment with bounds... 


  
  }
  //not implemented
  assert(0); 
}

int retrim::TrackMaker::clearTrack() 
{

  minleftover = 0; 
  maxleftover = 0; 

  ionization = 0; 
  particle_number.clear(); 
  startpoints.clear(); 
  endpoints.clear(); 
  starttimes.clear(); 
  endtimes.clear(); 
  ELoss.clear(); 
  nukeLoss.clear(); 
  randomized.clear(); 
  Eat.clear(); 
  species.clear(); 
  leftover.clear(); 
  particleCounter = 0; 
}

int retrim::TrackMaker::makeTrack(int itrack, const TVector3 * origin, const TVector3 * direction, double t0) 
{

  clearTrack(); 

  int ncoll = trim->nCollisions(itrack); 
  int primary_species = srim[0]->ionZ(); 

  double last_ion_time =t0;
  for (unsigned i = 0; i < ncoll; i++)
  {
    const CollisionReader::collision * c = trim->getCollision(itrack,i); 
    TVector3 endpoint(c->x, c->y, c->z); 
    endpoints.push_back(endpoint); 
    particle_number.push_back(0); 
    species.push_back(primary_species); 
    double nuke = 0; 
    if (i == 0)
    {
      startpoints.push_back(TVector3(0,0,0)); 
      starttimes.push_back(t0); 
      ElecLoss = trim->getE()- c->Ei;
      Eat.push_back(trim->getE()); 
      double dist = endpoints[endpoints.size()-1].Mag(); 
      endtimes.push_back(getT(trim->getE(),dist, srim[0]->getIonMass())); 
      predictedElecLoss = estimateElectronicLoss(srim[0], trim->getE(), dist , ratio,&nuke); 
//      ELoss.push_back(ElecLoss); 
      ELoss.push_back(predictedElecLoss); 
      nukeLoss.push_back(nuke); 
      randomized.push_back(false);
      stats->Fill(); 
      last_ion_time = endtimes[0]; 

    }
    else 
    {
      const CollisionReader::collision *  pc = trim->getCollision(itrack,i-1); 
      startpoints.push_back(TVector3(pc->x, pc->y, pc->z)); 
      starttimes.push_back(last_ion_time); 
      double Eireal = pc->Ei - pc->Er; 
      double Efreal = c->Ei; 
      ElecLoss = Eireal -Efreal; 
      Eat.push_back(Eireal); 
      double dist = (startpoints[startpoints.size()-1] - endpoints[endpoints.size()-1]).Mag(); 
      predictedElecLoss = estimateElectronicLoss(srim[0], Eireal,dist , ratio,&nuke); 
//      ELoss.push_back(ElecLoss); 
      ELoss.push_back(predictedElecLoss); 
      nukeLoss.push_back(nuke);
      randomized.push_back(false); 
      endtimes.push_back(last_ion_time + getT(Eireal, dist, srim[0]->getIonMass())); 
      last_ion_time = endtimes[endtimes.size()-1]; 
      stats->Fill(); 
    }

//    printf("collision is %d, Er = %f\n", i, c->Er); 

    assert(startpoints.size() ==  endpoints.size()); 

    //loop through secondaries 
    
    int nsec =  c->secondaries.size(); 
//    printf ("n secondaries: %d\n",nsec); 



    TVector3 lastdiff = endpoints[endpoints.size()-1] - startpoints[startpoints.size()-1]; 
    TVector3 lastprimarydiff = lastdiff; 

    //BUILD UP PARTICLE TREE 
    double tmp_maxleftover = maxleftover; 
    double tmp_minleftover = minleftover; 

    int niter = 0; 
    std::vector<particle> *best_particle_store = 0; 
    int best_iter = -1; 
    double bestscore = DBL_MAX; 

    while (true)
    {
      const CollisionReader::secondary_collision * s = &(c->secondaries[0]); 
      if (verbose) 
      {
        printf("iter: %d\n", niter); 
      }

      maxleftover = 0; 
      minleftover = 0; 
      std::vector<particle>* particle_store = new std::vector<particle>(nsec); 

      particle*  prime_recoil = &particle_store->at(0);  

      prime_recoil->parent = NULL; 
      prime_recoil->start = TVector3(s->x,s->y,s->z); 
      prime_recoil->startE = s->Er; 
      prime_recoil->startT = last_ion_time; 
      prime_recoil->species = s->atom; 
      prime_recoil->index = 0 ;
      prime_recoil->startdir = lastdiff; 

      particle * par = prime_recoil; 

      for (unsigned j = 1 ; j < nsec ; j++)
      {
         s = &(c->secondaries[j]); 

         particle * p = &particle_store->at(j); 
         p->startE = s->Er; 
         p->start = TVector3(s->x,s->y,s->z); 
         p->index = j; 
         p->startT = -1; //figure this out later
         p->species = s->atom; 

         //check feasibility of parent 
         while (par->parent!=NULL)
         {

           bool infeasible = false; 

           //check if between parent and first daughter?  //TODO check this logic and enable to protect against flagrant momentum violation
           
          /* if (!par->daughters.empty())
           {
             TVector3 par2me = p->start - par->start; 
             TVector3 me2d = par->daughters.front()->start - par->start; 

             if (par2me.Dot(me2d) > 0) 
               infeasible = true; 
           }
           */


           if (!infeasible)
           {
             double nuke =0;
             // check if energy compatible
             TVector3 par_pos = par->start; 
             double parE = par->startE - trim->getLatticeEnergy(par->species); 
             int par_srim= whichSRIM(par->species,nsrim,srim); 


             if (parE < trim->getDisplacementEnergy(par->species))
             {
                infeasible = true; 
             }

             else
             {
                 
               parE -= elec_fudge * estimateElectronicLoss(srim[par_srim], parE, (p->start - par_pos).Mag(), ratio,&nuke); 
               parE -= nuke_fudge * nuke; 
               parE -= p->startE; 
               par_pos = p->start; 


               for (std::list<particle*>::iterator it = par->daughters.begin(); it != par->daughters.end(); it++)
               {
                 if (parE < trim->getDisplacementEnergy(par->species)) 
                 {
                   infeasible = true; 
                   break;
                 }
                 particle * daughter = *it; 
                 parE -= elec_fudge * estimateElectronicLoss(srim[par_srim], parE, (daughter->start - par_pos).Mag(), ratio,&nuke); 
                 parE -= nuke_fudge * nuke; 
                 parE -= daughter->startE; 
                 par_pos = daughter->start; 
               }

               if (parE < 0) infeasible = true; 
             }

           }

           if (infeasible || (daughter_randomization && niter > 0 && rng.Uniform() < pow(p->startE / par->startE,4)) )
           {
             par = par->parent; 
           }
           else
           {
             break; 
           }
         }

         p->parent = par; 
         p->startdir = p->start - par->startdir; 
         par->daughters.push_front(p); 
         par = p; 

      }


      //compute score
      score = 0; 
      fillTree(prime_recoil,false); 

      if (verbose) 
      {
        printf("SCORE for iter %d is %f\n", niter,score); 
      }

      if ( score < bestscore)
      {
        if (best_particle_store) delete best_particle_store; 
        best_particle_store = particle_store; 
        bestscore = score; 
        best_iter = niter; 
      }
      else
      {
        delete particle_store; 
      }

      if (!daughter_randomization || ++niter > max_iter || ( minleftover > cut_minleftover && maxleftover < cut_maxleftover))
        break; 

    }

    particle * best_prime_recoil = &(best_particle_store->at(0)); 
    assert(best_prime_recoil->parent == NULL); 
    fillTree(best_prime_recoil,true); 

    if (verbose)
    {
      printf("Best iteration: %d\n",best_iter);  
      printTree(stdout, best_prime_recoil,0); 
    }

    if (best_prime_recoil) delete best_prime_recoil; 

    maxleftover = TMath::Max(maxleftover, tmp_maxleftover); 
    minleftover = TMath::Min(minleftover, tmp_minleftover); 


    //handle last segment
    if (i == ncoll - 1) 
    {
 //     printf("last segment!\n"); 
      startpoints.push_back(TVector3(c->x, c->y, c->z)); 
      starttimes.push_back(last_ion_time); 
      double Eactual = c->Ei - c->Er; 
      double howfar=0;
//      printf(" E->r = %f->%f\n",Eactual,howfar); 
      Eat.push_back(Eactual); 
      nuke = 0;
      if (Eactual < trim->getDisplacementEnergy(primary_species))
      {
        ELoss.push_back(0); 
        randomized.push_back(false); 
        nuke = trim->getDisplacementEnergy(primary_species); 
      }
       else
      { 
        howfar = srim[0]->getRange(Eactual); 
        ELoss.push_back(estimateElectronicLoss(srim[0], Eactual, howfar, ratio, &nuke)); 
        randomized.push_back(true); 
      }

      nukeLoss.push_back(nuke); 
      endtimes.push_back(last_ion_time + getT(Eactual, howfar, srim[0]->getIonMass())); 
      TVector3 last = startpoints[startpoints.size()-1] + lastprimarydiff.Unit() * howfar; 
      endpoints.push_back(last); 
      species.push_back(primary_species);
      particle_number.push_back(0); 
    }


  //  printf("startpoints: %d endpoints: %d \n", startpoints.size(), endpoints.size()); 
  }

  assert (startpoints.size() == endpoints.size()); 
  assert (startpoints.size() == nukeLoss.size()); 
  assert (startpoints.size() == starttimes.size()); 
  assert (startpoints.size() == endtimes.size()); 
  assert (startpoints.size() == ELoss.size()); 
  assert (startpoints.size() == Eat.size()); 
  assert (startpoints.size() == species.size()) ; 
  assert (startpoints.size() == particle_number.size()); 
  assert (startpoints.size() == randomized.size()); 


  fillIonizHists(); 

  transform(origin,direction); 

  for (unsigned i = 0; i < ELoss.size(); i++)
  {

    ionization += ELoss[i]; 
  }

}


void retrim::TrackMaker::fillIonizHists()
{

  if ((!ioniz_ion || !ioniz_recoils)) return; 

  nhist++; 

  for (unsigned i = 0; i < startpoints.size(); i++)
  {
    double x0 = startpoints[i].x() ; 
    double x1 = endpoints[i].x() ; 
    double E = ELoss[i]; 
    double N = nukeLoss[i]; 


    TH1 * ioniz = particle_number[i] == 0 ? ioniz_ion : ioniz_recoils;
    TH1 * phonon = particle_number[i] == 0 ? phonon_ion : phonon_recoils; 

    double dist_total = TMath::Abs(x1-x0); 
    double binsize = ioniz->GetBinWidth(1); 

    int startbin = ioniz->FindBin(x0);  
    int endbin = ioniz->FindBin(x1); 

    if (startbin == endbin)
    {
      ioniz->Fill(x0, E/binsize); 
      phonon->Fill(x0, N/binsize); 
    }

    else
    {
      for (int bin = startbin; bin <= endbin; bin++)
      {
        double frac = 0; 
        if (bin == startbin)
        {
          frac = x1 > x0 ? (ioniz->GetXaxis()->GetBinLowEdge(bin) + binsize - x0) / dist_total
                          : (x0 - ioniz->GetXaxis()->GetBinLowEdge(bin)) / dist_total; 
        }
        else if (bin == endbin)
        {

          frac = x0 > x1 ? (ioniz->GetXaxis()->GetBinLowEdge(bin) + binsize - x1) / dist_total
                          : (x1 - ioniz->GetXaxis()->GetBinLowEdge(bin)) / dist_total; 

        }
        else 
        {
          frac = binsize/dist_total; 
        }

//        printf("%g %g %g\n", frac,E,N); 

        ioniz->AddBinContent(bin,frac * E/binsize); 
        phonon->AddBinContent(bin,frac * N/binsize); 
      }
    }
  }


}

void retrim::TrackMaker::transform(const TVector3 * origin, const TVector3 * dir)
{

  TVector3  unit, orth, norm; 

  if (dir) 
  {
    unit = dir->Unit(); 
    orth = unit.Orthogonal(); 
    norm = unit.Cross(orth); 
  }

  for (unsigned i = 0; i < startpoints.size(); i++)
  {
    if(dir)
    {
      startpoints[i] = (startpoints[i].x() * unit + startpoints[i].y() * orth + startpoints[i].z() * norm); 
      endpoints[i] = (endpoints[i].x() * unit + endpoints[i].y() * orth + endpoints[i].z() * norm); 
    }
    if (origin)
    {
      startpoints[i] += *origin; 
      endpoints[i] += *origin; 
    }
  }
}


static int particle_colors[] = {kRed, kOrange, kGreen, kMagenta }; 

std::vector<TObject*> retrim::TrackMaker::draw(const char * opt) const
{

  std::map<int,int> species_color; 
  int color_index = kRed; 

  std::vector<TObject *> ret; 
  if (strchr(opt,'s'))
  {

    for (unsigned i = 0; i < startpoints.size(); i++)
    {
      TPolyLine3D * line = new TPolyLine3D (2); 
      line->SetPoint(0, startpoints[i].x(), startpoints[i].y(), startpoints[i].z()); 
      line->SetPoint(1, endpoints[i].x(), endpoints[i].y(), endpoints[i].z()); 

      if (!species_color.count(species[i]))
      {
        species_color[species[i]] = particle_colors[color_index++ % 4];
      }

      if (randomized[i])
      {
        line->SetLineColor(11); 
      }
      else if (particle_number[i] == 0) 
      {
        line->SetLineColor(4); 
      }
      else
      {
        line->SetLineColor(species_color[species[i]]); 
      }
      
//    line->SetLineWidth(int(ELoss[i]*10)+1) ; 
      line->Draw(); 
      ret.push_back(line); 
    }

  }

  if (strchr(opt,'e') && getNElectrons())
  {
    TPolyMarker3D * p = new TPolyMarker3D(getNElectrons()); 

    for (int i = 0; i < electron_x.size(); i++)
    {
      p->SetPoint(i, electron_x[i].X(), electron_x[i].Y(), electron_x[i].Z()); 
    }

    p->Draw(); 
    ret.push_back(p); 
  }

  return ret; 
}

void retrim::TrackMaker::fillElectrons(TH1 * fillme) const 
{
  for (int i = 0; i < getNElectrons(); i++)
  {
     double x = electron_x[i].X(); 
     double y = electron_x[i].Y(); 
     double z = electron_x[i].Z(); 
     int bin = fillme->FindBin(x,y,z); 
//     printf("%d %f %f %f\n",bin,x,y,z); 
     fillme->SetBinContent(bin,fillme->GetBinContent(bin)+1); 
  }
}

int retrim::TrackMaker::getElectron(int i, double & x, double & y, double & z, double & t, double & E, double & px, double & py, double & pz)
{
  if (i >= getNElectrons()) return 1; 
  
  x = electron_x[i].X(); 
  y = electron_x[i].Y(); 
  z = electron_x[i].Z(); 
  t = electron_x[i].T(); 

  px = electron_p[i].X(); 
  py = electron_p[i].Y(); 
  pz = electron_p[i].Z(); 
  E = electron_p[i].T(); 

  return 0; 
}
double retrim::TrackMaker::E() const
{
  return trim->getE(); 

}
int retrim::TrackMaker::ntracks() const
{

  return trim->nTracks(); 

}

int retrim::TrackMaker::makeElectrons(const IonizationModel * m)
{
  electron_x.clear(); 
  electron_p.clear(); 

  m->makeElectrons(startpoints.size(), &startpoints[0],&endpoints[0], &starttimes[0], &endtimes[0], &ELoss[0], 
      &Eat[0], &species[0], &electron_x, &electron_p); 

  return getNElectrons(); 
}

