//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>


void ThermalAccomodation_Sumpressure_cpu(real dt) {

//<USER_DEFINED>
  INPUT(Density);
  INPUT(Alphacol);
  INPUT(Qvec);
  INPUT(Slope);
  OUTPUT(Slope);
//<\USER_DEFINED>

//<EXTERNAL>
  real* gasdens = Fluids[0]->Density->field_cpu;
  real* dens = Density->field_cpu;
  real* sumpressure = Slope->field_cpu;
  real* q = Qvec->field_cpu;
  real* alpha = Alphacol->field_cpu;
  real* grk    = Gammark->field_cpu; 
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
  real sk;
  real cpgas;
  real cpdust;
  real beta;
//<\INTERNAL>

//<CONSTANT>
// real xmin(Nx+2*NGHX+1);
// real ymin(Ny+2*NGHY+1);
// real GAMMA(1);
//<\CONSTANT>


//<MAIN_LOOP>

  i =j =k =0;
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

	sk   =  alpha[ll]/(1+grk[ll]*dt*alpha[ll]);
	sumpressure[ll] +=  (GAMMA * cpdust/cpgas) *dens[ll]/gasdens[ll]*sk*q[ll];

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
