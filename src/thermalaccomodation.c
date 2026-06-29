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
  MULTIFLUID(ThermalRelaxation(dt)); 
#endif  
  //-------------------------------------------------------------------------------
  Reset_field(DensStar);  
  MULTIFLUID( ThermalAccomodation_Sumrho(dt)); 

  Reset_field(Slope);
  MULTIFLUID(ThermalAccomodation_Sumpressure(dt)); 

  MULTIFLUID(ThermalAccomodation_UpdateEnergy(dt)); 
  //Module ends here
  //-------------------------------------------------------------------------------------


  //-------------------------------------------------------------------------------------
  // COOLING TIME FOR SINGLE FLUID MODEL (well-mizex dust)

  if(NFLUIDS==1){
      int n;
      real stokes_plus[NDUST+1];
      real stokes[NDUST];
      real epsilons[NDUST];
      real smax = TSMAX;
      real smin = TSMIN;
      real ds = (log(smax) - log(smin)) / NDUST;

    //size distributioin
    for ( n = 0; n < NDUST+1; n++) stokes_plus[n] = smin * exp(ds * n);
    real slope = 4.0 - SQ;
    for ( n = 0; n < NDUST; n++) {
        if (slope != 0.0) {
            epsilons[n] = (pow(stokes_plus[n+1], slope) - pow(stokes_plus[n], slope)) * (EPSILON / (pow(smax, slope) - pow(smin, slope)));
        } else {
            epsilons[n] = log(stokes_plus[n+1] / stokes_plus[n]) * (EPSILON / log(smax / smin));
        }
        stokes[n] = sqrt(stokes_plus[n] * stokes_plus[n+1]);
    }

    Reset_field(Fluids[0]->Betarad);
    for ( n = 0; n < NDUST; n++) MULTIFLUID(if(Fluidtype==GAS) CoolingTime(dt, R0_CGS/stokes[n], RHOSOLID*R0_CGS*R0_CGS*R0_CGS/MSTAR_CGS,  epsilons[n], n));
  }  
      
  //-------------------------------------------------------------------------------------

  
}
