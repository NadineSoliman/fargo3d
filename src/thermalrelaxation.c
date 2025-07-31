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
  INPUT2D(Energy0);  
  INPUT2D(Density0);
  INPUT(Qs);
  OUTPUT(Energy);
  OUTPUT(Trad);
//<\USER_DEFINED>

//<EXTERNAL>
  real* energy = Energy->field_cpu;
  real* dens = Density->field_cpu;
  real* dens0 = Density0->field_cpu; 
  real* energy0 = Energy0->field_cpu;
  real* trad = Trad->field_cpu;
  int pitch2d = Pitch2D;  
  int pitch  = Pitch_cpu;
  int stride = Stride_cpu;
  int size_x = Nx;
  int size_y = Ny+2*NGHY;
  int size_z = Nz+2*NGHZ;
  int fluidtype = Fluidtype;
  real cpdg=CPDG;
  real invparticlesize = Coeffval[1];
  real rhosolid = Coeffval[2];
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
  real trdust_inv;
  real tempgas;
  real tempdust0;
  real tempdustn;
  real omega;
  real Q;
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

  //Dust op. thin cooling time due to radiative cooling
  tempdust0 = energy0[l2D] / (dens0[l2D]*cpdust);
  tempdustn = energy[ll] / (dens[ll]*cpdust);
  Q = 8 * 3.14159265359  * KBOLTZ * tempdustn/ invparticlesize / PLANCK / C0;
  if (Q >= 1.0) {
    trdust_inv = 12.0 * invparticlesize/ ( rhosolid * cpdust )* STEFANK * pow(tempdustn, 3.0);
  }
  else{
    trdust_inv = 3.14159265359 * 120.0 * STEFANK * KBOLTZ / (PLANCK * C0 * rhosolid * cpdust) * pow(tempdustn, 4.0);
  }
  temp   = ( tempdustn + tempdust0*dt*trdust_inv)/(1.+dt*trdust_inv);
  energy[ll] = dens[ll]* temp * cpdust; 
  trad[ll] = 1/trdust_inv;
  
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
