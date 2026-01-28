//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>

void RTD_Temperature_cpu () {

//<USER_DEFINED>
  INPUT(Energy);
  INPUT(Density);
  OUTPUT(Temperature);
//<\USER_DEFINED>


//<EXTERNAL>
  real* energy   = Energy->field_cpu;
  real* temp= Temperature->field_cpu;
  real* dens = Density->field_cpu;
  int pitch  = Pitch_cpu;
  int stride = Stride_cpu;
  int size_x = Nx+2*NGHX;
  int size_y = Ny+2*NGHY;
  int size_z = Nz+2*NGHZ;
  int fluidtype = Fluidtype;
  real cpdg = CPDG;
//<\EXTERNAL>

//<INTERNAL>
  int i;
  int j;
  int k;
  int ll;
  real cpgas;
  real cpdust;
//<\INTERNAL>

//<CONSTANT>
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
    cpgas  = GAMMA*R_MU/(GAMMA-1.0);
    cpdust = cpdg* cpgas;

	if(fluidtype==GAS)  temp[ll] = (GAMMA-1.0)*energy[ll]/(dens[ll]*R_MU);
    if(fluidtype==DUST) temp[ll] = energy[ll] / (dens[ll]*cpdust);
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
