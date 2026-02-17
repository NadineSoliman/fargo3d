//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>


void ThermalAccomodation(real dtl) {
  
  
 

    //Module starts here
    //-------------------------------------------------------------------------------------
    MULTIFLUID(ThermalAccomodation_Coeff(dtl));
    //-------------------------------------------------------------------------------------

    //-------------------------------------------------------------------------------
    Reset_field(DensStar);  
    MULTIFLUID(if(Fluidtype==DUST) ThermalAccomodation_Sumrho(dtl)); // Return "B" in in MDIRK method
  
    //-------------------------------------------------------------------------------
    // compute q1
    MULTIFLUID(ThermalAccomodation_ComputeQ(dtl,0.)); //0 means we compute q1
    Reset_field(Slope);
    MULTIFLUID(if(Fluidtype==DUST) ThermalAccomodation_Sumpressure(dtl)); // Return "A" in in MDIRK method
    //compute k1
    MULTIFLUID(ThermalAccomodation_ComputeK(dtl,1)); 
    //-------------------------------------------------------------------------------

    //-------------------------------------------------------------------------------
    // compute q2
    MULTIFLUID(ThermalAccomodation_ComputeQ(dtl,1.)); //1 means we compute q2
    Reset_field(Slope);
    MULTIFLUID( if(Fluidtype==DUST) ThermalAccomodation_Sumpressure(dtl)); // Return "A" in in MDIRK method  
    //compute k2
    MULTIFLUID(ThermalAccomodation_ComputeK(dtl,2)); 
    //-------------------------------------------------------------------------------
  
    //update energy
    MULTIFLUID(ThermalAccomodation_UpdateEnergy(dtl)); 
    //Module ends here
    //-------------------------------------------------------------------------------------
   


  }