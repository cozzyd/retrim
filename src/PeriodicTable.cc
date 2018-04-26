#include "PeriodicTable.hh" 
#include <cstring>


static const char * elements[] = 
{
  "H", //1
  "He", 
  "Li",  
  "Be",
  "B", 
  "C", //6 
  "N",
  "O", 
  "F", 
  "Ne", //10
  "Na", 
  "Mg", 
  "Al", 
  "Si", 
  "P",  //15
  "S", 
  "Cl", 
  "Ar", 
  "K", 
  "Ca" //20 

  //TODO: Finish this!!! 

};

static const int nimpl = 20; 


int retrim::PeriodicTable::getZ(const char * element) 
{
  for (int i = 0; i < nimpl; i++)
  {
    if (!strncmp(element,elements[i],strlen(elements[i]))) return i+1; 
  }

  return 0; 
}


const char * retrim::PeriodicTable::getSymbol(int Z)
{
  if (Z > 0 && Z < nimpl)
  {
    return elements[Z-1]; 
  }

  else return 0; 
}

