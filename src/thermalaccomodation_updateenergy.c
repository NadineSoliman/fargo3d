//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>


void ThermalAccomodation_UpdateEnergy_cpu(real dt) {

//<USER_DEFINED>
  INPUT(Alphacol);
  INPUT(Energy);
  INPUT(Density);
  INPUT(DensStar);
  INPUT(Slope);
  OUTPUT(Energy);
//<\USER_DEFINED>

//<EXTERNAL>
  real* energy = Energy->field_cpu;
  real* dens = Density->field_cpu;
  real* alpha = Alphacol->field_cpu;
  real* sumpressure = Slope->field_cpu;
  real* sumrho = DensStar->field_cpu;
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
  real tempgas;
  real tempdustn;
  real tempdust;
  real tempgasn;
  real dtl;
//<\INTERNAL>

//<CONSTANT>
// real xmin(Nx+2*NGHX+1);
// real ymin(Ny+2*NGHY+1);
// real GAMMA(1);
//<\CONSTANT>

//<MAIN_LOOP>

  i =j= k= 0;

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

  
	if (fluidtype == GAS)  {
    tempgas    = sumpressure[ll]/( sumrho[ll] );
    energy[ll] = (dens[ll]*R_MU) * tempgas/(GAMMA - 1.0); 
  }
	else{
    dtl = (exp(alpha[ll] * dt) - 1.0)/alpha[ll];
    tempdustn = energy[ll] / (dens[ll]*cpdust);
    tempdust = (dtl*alpha[ll])/(1.+ dtl*alpha[ll])*sumpressure[ll]/( sumrho[ll] ) + tempdustn/(1.+ dtl*alpha[ll]);
    energy[ll] = (dens[ll]*cpdust) * tempdust; 
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
