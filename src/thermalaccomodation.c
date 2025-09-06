//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>


void ThermalAccomodation(real dt) {
  
  
  //Subcycling for the thermal accomodation
  int nsub;
  real dtl;
  real sumdt=0.;
  real dtlocal = MIN(dt,StepTimeCol); //if dt<StepTimeCol we do not subcycle
  int nsteps   = (int)(dt/dtlocal);  // nsteps=1 if dt<StepTimeCol
  dtl = max2(dt,StepTimeCol)/nsteps; //dtl is the subcycled timestep 

 
  for (nsub = 0; nsub < nsteps; nsub++){
    sumdt += dtl;
    if(sumdt>dt) dtl = dt - (sumdt-dtl); 

    //Module starts here
    //-------------------------------------------------------------------------------------
    MULTIFLUID(ThermalAccomodation_Coeff(dtl));
    //-------------------------------------------------------------------------------------
#ifdef THERMALRELAXATION
    MULTIFLUID(if(Fluidtype==DUST) ThermalRelaxation(dtl)); 
#endif  
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

  }