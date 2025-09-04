//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>


void ThermalAccomodation_UpdateEnergy_cpu(real dt) {

//<USER_DEFINED>
  INPUT(Rkk1);
  INPUT(Rkk2);
  INPUT(Energy);
  INPUT(Density);
  OUTPUT(Energy);
//<\USER_DEFINED>

//<EXTERNAL>
  real* energy = Energy->field_cpu;
  real* dens = Density->field_cpu;
  real* rkk1 = Rkk1->field_cpu;
  real* rkk2 = Rkk2->field_cpu;
  real* grk  = Gammark->field_cpu; 
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
    tempgasn   =  (GAMMA-1.0)*energy[ll]/(dens[ll]*R_MU);
    tempgas    = tempgasn + dt*(1.0-grk[ll])*rkk1[ll] + grk[ll]*dt*rkk2[ll];
    energy[ll] = (dens[ll]*R_MU) * tempgas/(GAMMA - 1.0); 
  }
	else{
    tempdustn = energy[ll] / (dens[ll]*cpdust);
    tempdust  = tempdustn + dt*(1.0-grk[ll])*rkk1[ll] + grk[ll]*dt*rkk2[ll];
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
