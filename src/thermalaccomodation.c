//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>


void ThermalAccomodation(real dt) {
  
  
  // Centered drag coefficient pre-factor
  FARGO_SAFE(ThermalAccomodation_Coeff()); //store dragcoeff in Qs
  //-------------------------------------------------------------------------------------
  Reset_field(DensStar);  
  MULTIFLUID(ThermalAccomodation_Sumrho(dt)); 
  // WriteMerging(DensStar, 0); //write DensStar to disk
  // exit(22);
  Reset_field(Slope);
  MULTIFLUID(ThermalAccomodation_Sumpressure(dt));
  
  MULTIFLUID(ThermalAccomodation_UpdateEnergy(dt)); //update velocities
}