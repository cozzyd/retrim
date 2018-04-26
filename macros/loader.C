/* Loader for TLorentzVector... not sure why it's necessary */ 
#include "TLorentzVector.h"
#include <vector>
#ifdef __MAKECINT__
#pragma link C++ class std::vector<TLorentzVector>+;
#endif
