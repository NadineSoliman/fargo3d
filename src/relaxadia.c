//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>

void RelaxAdia_cpu (real dt) {

//<USER_DEFINED>
  INPUT(Energy);
  INPUT(Density);
  INPUT2D(Energy0);
  INPUT2D(Density0);
  OUTPUT(Energy);
//<\USER_DEFINED>

//<EXTERNAL>
  real* e     = Energy->field_cpu;
  real* rho   = Density->field_cpu;
  real* e0     = Energy0->field_cpu;
  real* rho0   = Density0->field_cpu;
  int pitch2d = Pitch2D;
  int pitch  = Pitch_cpu;
  int stride = Stride_cpu;
  int size_x = XIP; 
  int size_y = Ny+2*NGHY-1;
  int size_z = Nz+2*NGHZ-1;
  real beta = EBETA;
//<\EXTERNAL>

//<INTERNAL>
  int i; //Variables reserved
  int j; //for the topology
  int k; //of the kernels
  int ll;
#ifdef X
  int llxp;
#endif
#ifdef Y
  int llyp;
#endif
#ifdef Z
  int llzp;
#endif
  real omega;
  real temp;
  real temp0;
  real tau;
  real soundspeed;
  real scaleheight;
  real trelax;
  //<\INTERNAL>
  
//<CONSTANT>
// real NU(1);
// real GAMMA(1);
// real ALPHA(1);
// real Sxj(Ny+2*NGHY);
// real Syj(Ny+2*NGHY);
// real Szj(Ny+2*NGHY);
// real Sxk(Nz+2*NGHZ);
// real Syk(Nz+2*NGHZ);
// real Szk(Nz+2*NGHZ);
// real ymin(Ny+2*NGHY+1);
// real zmin(Nz+2*NGHZ+1);
// real InvVj(Ny+2*NGHY);
//<\CONSTANT>


//<MAIN_LOOP>
  
  i = j = k = 0;
  
#ifdef Z
  for(k=0; k<size_z; k++) {
#endif
#ifdef Y
    for(j=0; j<size_y; j++) {
#endif
#ifdef X
      for(i=0; i<size_x; i++) {
#endif
//<#>

	ll = l;
#ifdef SHEARINGBOX
	omega = OMEGAFRAME;
#endif
#ifdef CYLINDRICAL
  omega = sqrt(G*MSTAR/ymed(j)/ymed(j)/ymed(j));
#endif
	//equilibrium (initial) temperature
	temp0 = (GAMMA-1.0)*e0[l2D]/(rho0[l2D]*R_MU);
	// 	

	trelax = beta/omega;
	temp   = ( (GAMMA-1.0)*e[ll]/(rho[ll]*R_MU)  + temp0*dt/trelax)/(1.+dt/trelax);

	e[ll]  = rho[ll]*temp*R_MU/(GAMMA-1);
	
//<\#>
#ifdef X
      }
#endif
#ifdef Y
    }
#endif
#ifdef Z
  }
  //    exit(33);

#endif
//<\MAIN_LOOP>
}
