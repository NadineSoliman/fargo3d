//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>

void ThermalAccomodation_ComputeQ_cpu (real dt, int option) {
 
//<USER_DEFINED>
  INPUT(Energy);
  INPUT(Density);
#ifdef THERMALRELAXATION
  INPUT(Betarad);
#endif
  INPUT(Rkk1);
  INPUT2D(Energy0);  
  INPUT2D(Density0);
  OUTPUT(Qvec)
//<\USER_DEFINED>

//<EXTERNAL>
  real* energy = Energy->field_cpu;
  real* dens = Density->field_cpu;
  real* dens0 = Density0->field_cpu; 
  real* energy0 = Energy0->field_cpu;
  real* q = Qvec->field_cpu;
#ifdef THERMALRELAXATION
  real* beta = Betarad->field_cpu;
#endif
  real* rkk1 = Rkk1->field_cpu;
  real* grk  = Gammark->field_cpu; 
  int pitch2d = Pitch2D;  
  int pitch  = Pitch_cpu;
  int stride = Stride_cpu;
  int size_x = Nx;
  int size_y = Ny+2*NGHY;
  int size_z = Nz+2*NGHZ;
  int fluidtype = Fluidtype;
  real cpdg=CPDG;
//<\EXTERNAL>

//<INTERNAL>
  int i;
  int j;
  int k;
  int ll;
  real cpdust;
  real cpgas;
  real temp;
  real temp0;
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

    if(fluidtype==GAS) {
        temp=(GAMMA-1.0)*energy[ll]/(dens[ll]*R_MU);
        temp0 = (GAMMA-1.0)*energy0[l2D] / (dens0[l2D]*R_MU);
    }
    else {
        temp = energy[ll] / (dens[ll]*(cpdust));
        temp0 = energy0[l2D] / (dens0[l2D]*cpdust);
    }

    q[ll] = temp  + option*(1.0-grk[ll])*dt*rkk1[ll];
#ifdef THERMALRELAXATION
    q[ll] += beta[ll]*dt*temp0;
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
