//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>

void ThermalAccomodation_Sumrho_cpu (real dt) {
 
//<USER_DEFINED>
  INPUT(Density);
  INPUT(Qs);
  INPUT(DensStar);
  OUTPUT(DensStar);
//<\USER_DEFINED>

//<EXTERNAL>
  real* dens = Density->field_cpu;
  real* sumrho    = DensStar->field_cpu;
  real* pref = Qs->field_cpu;
  int pitch  = Pitch_cpu;
  int stride = Stride_cpu;
  int size_x = Nx;
  int size_y = Ny+2*NGHY;
  int size_z = Nz+2*NGHZ;
  int fluidtype = Fluidtype;
  real invparticlesize = Coeffval[1];
  real rhosolid        = Coeffval[2];
//<\EXTERNAL>

//<INTERNAL>
  int i,j,k;
  int ll;
  real alphak;
  real sk;
//<\INTERNAL>

//<CONSTANT>
// real xmin(Nx+2*NGHX+1);
// real ymin(Ny+2*NGHY+1);
//<\CONSTANT>


//<MAIN_LOOP>

  i = j =k =0;

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

#ifdef DUSTSIZE
	alphak  = pref[ll]*invparticlesize/rhosolid;

	sk      = (CP_DUST/CP_GAS) * dt*alphak/(1+dt*alphak);

	if (fluidtype == GAS)  sk = 1.0;
	sumrho[ll] += dens[ll]*sk;
#endif
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
