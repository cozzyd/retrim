/** \file LinkDef.h
 CINT LinkDef file for this package
*/
#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ nestedtypedefs;
#pragma link C++ nestedclasses;

#pragma link C++ namespace retrim;
#pragma link C++ namespace retrim::PeriodicTable;
#pragma link C++ namespace retrim::HistogramReader;
#pragma link C++ class retrim::TableReader;
#pragma link C++ class retrim::CollisionReader;
#pragma link C++ class retrim::TrackMaker;
#pragma link C++ class retrim::EnergyLossModel;
#pragma link C++ class retrim::ConstantEnergyLossModel;
#pragma link C++ class retrim::IonizationModel;
#pragma link C++ class retrim::SimpleIonizationModel;
#pragma link C++ class retrim::CollisionReader::collision;
#pragma link C++ class retrim::CollisionReader::secondary_collision;

#endif

