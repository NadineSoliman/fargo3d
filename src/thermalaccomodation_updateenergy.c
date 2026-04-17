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
  INPUT(Betarad);
  INPUT(Energy);
  INPUT2D(Energy0);
  INPUT2D(Density0);
  INPUT(Density);
  INPUT(DensStar);
  INPUT(Slope);
  OUTPUT(Energy);
  OUTPUT(Tcol);
//<\USER_DEFINED>

//<EXTERNAL>
  real* energy = Energy->field_cpu;
  real* energy0   = Energy0->field_cpu;
  real* dens = Density->field_cpu;
  real* dens0 = Density0->field_cpu;
  real* alpha = Alphacol->field_cpu;
  real* beta = Betarad->field_cpu;
  real* sumpressure = Slope->field_cpu;
  real* sumrho = DensStar->field_cpu;
  real* tcol = Tcol->field_cpu;
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
  real tempi;
  real tempdust;
  real temp0;
  real tempn;
  // real tempgasn;
  real energyn;
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

  cpgas  = R_MU/(GAMMA-1.0);
  cpdust = cpdg* cpgas;
  
	if (fluidtype == GAS)  {
    tempn =energy[ll]/dens[ll]/cpgas;
    tempgas = sumpressure[ll]/ sumrho[ll];
    energy[ll] = (dens[ll]* tempgas * cpgas);
    temp0 =energy0[l2D]/dens0[l2D]/cpgas;
    tcol[ll] = 1/(log((tempgas - temp0)/(tempn - temp0))/dt);

  }
	else{
    dtl = dt;//(1.0 - exp( -1.0* alpha[ll] * dt));
    tempdustn = energy[ll] / (dens[ll]*cpdust);
    temp0 = energy0[l2D]/dens0[l2D]/cpdust;
    tempdust = tempdustn  + (beta[ll] * dtl* temp0) + (alpha[ll] * dtl * sumpressure[ll]/sumrho[ll]);
    tempdust /= (1 + dtl * (alpha[ll] + beta[ll]));
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
