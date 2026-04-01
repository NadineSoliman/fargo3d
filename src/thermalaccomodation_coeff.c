//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>

void ThermalAccomodation_Coeff_cpu(real dt) {

//<USER_DEFINED>
  INPUT(Density);
  INPUT(Energy);
  OUTPUT(Alphacol);
//<\USER_DEFINED>

//<EXTERNAL>
  real* alpha      = Alphacol->field_cpu;
  real* dens_gas   = Fluids[0]->Density->field_cpu;
  real* energy_gas = Fluids[0]->Energy->field_cpu; 
  real* dens = Density->field_cpu;
  real* energy = Energy->field_cpu;
  int pitch  = Pitch_cpu;
  int stride = Stride_cpu;
  int size_x = Nx;
  int size_y = Ny+2*NGHY;
  int size_z = Nz+2*NGHZ;
  real invparticlesize = Coeffval[1];
  real rhosolid        = Coeffval[2];
  int fluidtype = Fluidtype;
  real cpdg=CPDG;
//<\EXTERNAL>

//<INTERNAL>
  int i;
  int j;
  int k;
  int ll;
  real tempgas;
  real cpdust;
  real cpgas;
  real omega;
  #ifdef CONSTANTTHERMALCOEFF
  real invthermaltime=invparticlesize;
#endif
//<\INTERNAL>

//<CONSTANT>
// real xmin(Nx+2*NGHX+1);
// real ymin(Ny+2*NGHY+1);
// real GAMMA(1);
//<\CONSTANT>
  
//<MAIN_LOOP>

  i = j = k = 0;

#ifdef Z
  for (k=0; k<size_z; k++) {
#endif
#ifdef Y
    for (j=0; j<size_y; j++) {
#endif
#ifdef X
      for (i=0; i<size_x; i++ ) {
#endif
//<#>
	ll = l;
	 omega     = sqrt(G*MSTAR/ymed(j)/ymed(j)/ymed(j));
   cpgas      = GAMMA*R_MU/(GAMMA-1.0);
   cpdust    = cpdg* cpgas;

  if (fluidtype == GAS) {
    alpha[ll] = 0.0; //gas
  }
  else {
#ifdef CONSTANTTHERMALCOEFF
      alpha[ll] = 1.0*invthermaltime;
#endif
#ifdef DUSTSIZE
  tempgas    =  (GAMMA-1.0)*energy_gas[ll]/(dens_gas[ll]*R_MU);
  alpha[ll]  = 0.75 *pow(KBOLTZ/MH, 1.5)* pow(tempgas, 0.5) * dens_gas[ll];
  alpha[ll] *= invparticlesize/rhosolid/cpdust;
#endif	
  }
/* Radiative relaxation rate for dust is computed inside thermalrelaxation.c (no Betarad field). */
//<\#>
#ifdef X
      }
#endif
#ifdef Y
    }
#endif
#ifdef Z
  }
#endif
//<\MAIN_LOOP>
}
