//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>


void ThermalRelaxation_cpu(real dt) {

//<USER_DEFINED>
  INPUT(Energy);
  INPUT(Density);
  OUTPUT(Energy);
//<\USER_DEFINED>

//<EXTERNAL>
  real* energy = Energy->field_cpu;
  real* dens = Density->field_cpu;
  int pitch  = Pitch_cpu;
  int stride = Stride_cpu;
  int size_x = Nx;
  int size_y = Ny+2*NGHY;
  int size_z = Nz+2*NGHZ;
  int fluidtype = Fluidtype;
  real cpdg=CPDG;
  real tgasup=TGASUP;
  real tgasmid=TGASMID;
  real tdustup=TDUSTUP;
  real tdustmid=TDUSTMID;
//<\EXTERNAL>

//<INTERNAL>
  int i;
  int j;
  int k;
  int ll;
  real temp;
  real cpdust;
  real cpgas;
  real trgas;
  real trdust;
  real tempgas;
  real tempdust;
  real omega;
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
	omega = sqrt(G*MSTAR/ymed(j)/ymed(j)/ymed(j));

	if (fluidtype == GAS)  {
    //Gas op. thin cooling time due to line cooling
    tempgas = (tgasmid + (tgasup-tgasmid)*pow(cos(zmed(k))/0.15,8))/ymed(j);
    trgas   = 0.1/omega;
    temp   = ( (GAMMA-1.0)*energy[ll]/(dens[ll]*R_MU)  + tempgas*dt/trgas)/(1.+dt/trgas);
    energy[ll] = (dens[ll]*R_MU) * temp/(GAMMA - 1.0); 
  }
	else{
    //Dust op. thin cooling time due to radiative cooling
    tempdust = (tdustmid + (tdustup-tdustmid)*pow(cos(zmed(k))/0.15,8))/ymed(j);
    trdust   = 0.01/omega;
    temp   = ( energy[ll] / (dens[ll]*(cpdust))  + tempdust*dt/trdust)/(1.+dt/trdust);
    energy[ll] = (dens[ll]*cpdust) * temp; 
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
