//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>

void mon_kinetic_cpu () {

//<USER_DEFINED>
  INPUT(Density);
  INPUT(Vz);
  INPUT(Vy);
  OUTPUT(Slope);
//<\USER_DEFINED>


//<EXTERNAL>
  real* dens = Density->field_cpu;
  real* vz = Vz->field_cpu;
  real* vy = Vy->field_cpu;
  real* interm = Slope->field_cpu;
  int pitch  = Pitch_cpu;
  int pitch2d = Pitch2D;
  int stride = Stride_cpu;
  int size_x = Nx+2*NGHX;
  int size_y = Ny+2*NGHY-1;
  int size_z = Nz+2*NGHZ;
//<\EXTERNAL>

//<INTERNAL>
  int i;
  int j;
  int k;
  int ll;
  real vyc;
  real vzc;
//<\INTERNAL>

//<CONSTANT>
// real Sxi(Nx);
// real Syk(Nz+2*NGHZ);
// real InvVj(Ny+2*NGHY);
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
  vyc = 0.5*(vy[ll]+vy[lyp]);
  vzc = 0.5*(vz[ll]+vz[lzp]);
  
	interm[ll] = 0.5*dens[ll]*(vyc*vyc + vzc*vzc)*Vol(i,j,k);
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
