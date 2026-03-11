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
  //-------------------------------------------------------------------------------
  Reset_field(DensStar);  
  MULTIFLUID( ThermalAccomodation_Sumrho(dt)); 

  Reset_field(Slope);
  MULTIFLUID(ThermalAccomodation_Sumpressure(dt)); 

  MULTIFLUID(ThermalAccomodation_UpdateEnergy(dt)); 
  //Module ends here
  //-------------------------------------------------------------------------------------
  }