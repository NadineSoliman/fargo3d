//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>


// Function to compute S_N = sum(tan(i/N)) for i=0 to N (uniform)
real computesn(int n) {
    real sum = 0.0;
    real u;
    int i;
    for (i = 0; i < n; i++) {
      u = (real)i / ((n-1)+1e-19);
      sum += tan(u*u);
    }
    return sum;
}

void ThermalAccomodation(real dt) {
  
  
  //Subcycling for the thermal accomodation
  int nsub,m, nsteps;
  real dtl, dtmax, sn, amp, dtnsteps, u;
  
  dtnsteps = MIN(dt,StepTimeCol); 
  dtmax    = max2(dt,StepTimeCol); 
  nsteps   = (int)(dt/dtnsteps);  
  if(dt<StepTimeCol) NSUBTH=1;
  nsub     = nsteps/NSUBTH; //

  sn  = computesn(nsub);
  amp = (dtmax - nsub*StepTimeCol)/sn;
 
  real sumdt=0.;
  for (m = 0; m < nsub; m++){
    u = (real)m /(nsub-1+1e-19);
    dtl = amp*tan( u*u )+dtnsteps;
    sumdt += dtl;
    //if(sumdt>dt) dtl = dt - (sumdt-dtl); 

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