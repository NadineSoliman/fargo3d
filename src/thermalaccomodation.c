//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>


void ThermalAccomodation(real dt) {
  
  //Module starts here
  //-------------------------------------------------------------------------------------
  MULTIFLUID(ThermalAccomodation_Coeff(dt));
  //-------------------------------------------------------------------------------------
#ifdef THERMALRELAXATION
   MULTIFLUID(if(Fluidtype==DUST) ThermalRelaxation(dt)); 
#endif  
  
  Reset_field(DensStar);  
  MULTIFLUID(if(Fluidtype==DUST) ThermalAccomodation_Sumrho(dt)); // Return "B" in in MDIRK method
  
  //-------------------------------------------------------------------------------
  // compute q1
  MULTIFLUID(ThermalAccomodation_ComputeQ(dt,0.)); //0 means we compute q1
  Reset_field(Slope);
  MULTIFLUID(if(Fluidtype==DUST) ThermalAccomodation_Sumpressure(dt)); // Return "A" in in MDIRK method
  //compute k1
  MULTIFLUID(ThermalAccomodation_ComputeK(dt,1)); 
  //-------------------------------------------------------------------------------

  //-------------------------------------------------------------------------------
  // compute q2
  MULTIFLUID(ThermalAccomodation_ComputeQ(dt,1.)); //1 means we compute q2
  Reset_field(Slope);
  MULTIFLUID( if(Fluidtype==DUST) ThermalAccomodation_Sumpressure(dt)); // Return "A" in in MDIRK method  
  //compute k2
  MULTIFLUID(ThermalAccomodation_ComputeK(dt,2)); 
  //-------------------------------------------------------------------------------
  
  //update energy
  MULTIFLUID(ThermalAccomodation_UpdateEnergy(dt)); 
  //Module ends here
  //-------------------------------------------------------------------------------------
  
}