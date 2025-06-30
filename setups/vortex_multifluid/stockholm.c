//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>

void StockholmBoundary_cpu(real dt) {

//<USER_DEFINED>
  INPUT(Density);
  INPUT2D(Density0);
  OUTPUT(Density);
#ifdef ADIABATIC
  INPUT(Energy);
  INPUT2D(Energy0);
  OUTPUT(Energy);
#endif
#ifdef X
  INPUT(Vx);
  INPUT2D(Vx0);
  OUTPUT(Vx);
#endif
#ifdef Y
  INPUT(Vy);
  INPUT2D(Vy0);
  OUTPUT(Vy);
#endif
#ifdef Z
  INPUT(Vz);
  INPUT2D(Vz0);
  OUTPUT(Vz);
#endif
//<\USER_DEFINED>

//<EXTERNAL>
  real* rho  = Density->field_cpu;
  real* rho0 = Density0->field_cpu;
#ifdef X
  real* vx  = Vx->field_cpu;
  real* vx0 = Vx0->field_cpu;
#endif
#ifdef Y
  real* vy  = Vy->field_cpu;
  real* vy0 = Vy0->field_cpu;
#endif
#ifdef Z
  real* vz  = Vz->field_cpu;
  real* vz0 = Vz0->field_cpu;
#endif
#ifdef ADIABATIC
  real* e    = Energy->field_cpu;
  real* e0   = Energy0->field_cpu;
#endif
  int pitch   = Pitch_cpu;
  int stride  = Stride_cpu;
  int size_x  = Nx+2*NGHX;
  int size_y  = Ny+2*NGHY;
  int size_z  = Nz+2*NGHZ;
  int pitch2d = Pitch2D;
  real y_min = YMIN;
  real y_max = YMAX;
  real z_min = ZMIN;
  real z_max = ZMAX;
  real dampingzone = DAMPINGZONE;
  real kbcol = KILLINGBCCOLATITUDE;
  real of    = OMEGAFRAME;
  real of0   = OMEGAFRAME0;
  real r0 = R0;
  real ds = TAUDAMP;
  real dy = H0;
  int periodic_z = PERIODICZ;
//<\EXTERNAL>

//<INTERNAL>
  int i;
  int j;
  int k;
  //  Similar to Benitez-Llambay et al. (2016), Eq. 7.
  real Y_inf = y_min + dy;
  real Y_sup = y_max - dy;
  real radius;
  real vx0_target;
  real rampy;
  real rampz;
  real rampzz;
  real rampi;
  real ramp;
  real tau;
  real taud;
//<\INTERNAL>

//<CONSTANT>
// real xmin(Nx+1);
// real ymin(Ny+2*NGHY+1);
// real zmin(Nz+2*NGHZ+1);
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
      for (i=0; i<size_x; i++) {
#endif
//<#>
	rampy = 0.0;
	rampz = 0.0;
	rampzz = 0.0;

	if(ymed(j) > Y_sup) {
	  rampy   = (ymed(j)-Y_sup)/(y_max-Y_sup);
	}
	if(ymed(j) < Y_inf) {
	  rampy   = (Y_inf-ymed(j))/(Y_inf-y_min);
	}
	rampy *= rampy;		/* Parabolic ramp as in De Val Borro et al (2006) */
	
	ramp = rampy;
	rampi= rampy;
	tau = ds;
	
	if(ramp>0.0) {
	  taud = tau/ramp;
	  vy[l] = (vy[l]*taud+vy0[l2D]*dt)/(dt+taud);
	  vx[l] = (vx[l]*taud+vx0[l2D]*dt)/(dt+taud);
	}
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
