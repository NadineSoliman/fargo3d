//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>

void ThermalAccomodation_Coeff_cpu () {

//<USER_DEFINED>
#ifdef DUSTSIZE
  INPUT(Density);
  INPUT(Energy);
#endif
  OUTPUT(Qs);
//<\USER_DEFINED>

//<EXTERNAL>
  real* coeff      = Qs->field_cpu;
#ifdef DUSTSIZE
  real* dens_gas   = Fluids[0]->Density->field_cpu;
  real* energy_gas = Fluids[0]->Energy->field_cpu; 
#endif
  int pitch  = Pitch_cpu;
  int stride = Stride_cpu;
  int size_x = Nx;
  int size_y = Ny+2*NGHY;
  int size_z = Nz+2*NGHZ;
  real thermalcoeff = 	THERMALACCOMODATIONCOEFF;
//<\EXTERNAL>

//<INTERNAL>
  int i;
  int j;
  int k;
  int ll;
  real omega;
  real cs;
  real tempgas;
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

	cs = sqrt(GAMMA*(GAMMA-1)*energy_gas[ll]/dens_gas[ll]);
	
  tempgas   =  (GAMMA-1.0)*e[ll]/(rho[ll]*R_MU);

  coeff[ll] = thermalcoeff*KBOLTZ*MH..... tempgas;
	
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
