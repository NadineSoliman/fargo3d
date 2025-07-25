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
  real alphak;
  real sk;
  real cpgas;
  real cpdust;
  real temp;
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


  cpgas  = GAMMA*R_MU/(GAMMA-1.0);
  cpdust = cpdg* cpgas;
  alphak = 0.0;
#ifdef CONSTANTTHERMALCOEFF
        alphak = pref[ll]*invthermaltime;
#endif
#ifdef DUSTSIZE
  alphak = pref[ll]*invparticlesize/rhosolid/cpdust;
#endif
	sk      = (GAMMA * cpdust/cpgas) * dt*alphak/(1+dt*alphak);

	if (fluidtype == GAS)  sk = 1.0;
	sumrho[ll] += dens[ll]*sk;

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
