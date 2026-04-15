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
  INPUT(Energy);
  INPUT(Alphacol);
  INPUT(Betarad);
  INPUT2D(Density0);
  INPUT2D(Energy0);
  INPUT(Slope);
  OUTPUT(Slope);
//<\USER_DEFINED>

//<EXTERNAL>
  real* gasdens = Fluids[0]->Density->field_cpu;
  real* dens = Density->field_cpu;
  real* energy = Energy->field_cpu;
  real* sumpressure = Slope->field_cpu;
  real* alpha = Alphacol->field_cpu;
  real* beta = Betarad->field_cpu;
  real* dens0 = Density0->field_cpu;
  real* energy0 = Energy0->field_cpu;
  
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
  real cpgas = GAMMA*R_MU/(GAMMA-1.0);
  real cpdust = cpdg* cpgas;
  real rhotemp;
  real temp0;
  real tempn;
  real dtl;
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
      
  // dtl = (exp(  alpha[ll] * dt)  - 1.0)/alpha[ll];
  dtl = dt;
  rhotemp = energy[ll]/cpdust;
  temp0 = energy0[l2D]/cpdust/dens0[l2D];
  tempn = energy[ll]/cpdust/dens[ll];
	sk   =  cpdg* dtl* alpha[ll]*(1 + temp0/tempn * beta[ll]*dtl)  /(1+dtl*(alpha[ll] + beta[ll]));

  if (fluidtype == GAS) {
    sk = 1.0;
    rhotemp = energy[ll]/cpgas;
  }
  sumpressure[ll] +=  rhotemp *sk;

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
