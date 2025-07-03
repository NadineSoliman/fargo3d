//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>


void ThermalAccomodation_UpdateEnergy_cpu(real dt) {

//<USER_DEFINED>
  INPUT(DensStar);
  INPUT(Slope);
  INPUT(Qs);
  INPUT(Energy);
  OUTPUT(Energy);
//<\USER_DEFINED>

//<EXTERNAL>
  real* sumpressure = Slope->field_cpu;
  real* sumrho = DensStar->field_cpu;
  real* pref = Qs->field_cpu;
  real* energy = =Energy->field_cpu;
  int pitch  = Pitch_cpu;
  int stride = Stride_cpu;
  int size_x = Nx;
  int size_y = Ny+2*NGHY;
  int size_z = Nz+2*NGHZ;
  int fluidtype = Fluidtype;
  real invstokesnumber = Coeffval[0];
  real invparticlesize = Coeffval[1];
  real rhosolid        = Coeffval[2];
//<\EXTERNAL>

//<INTERNAL>
  int i;
  int ll;
  real alphak;
  real sk;
  real omega;
  real temp;
//<\INTERNAL>

//<CONSTANT>
// real xmin(Nx+2*NGHX+1);
// real ymin(Ny+2*NGHY+1);
//<\CONSTANT>

//<MAIN_LOOP>

  i = 0;

#ifdef X
      for (i=0; i<size_x; i++ ) {
#endif
//<#>
	ll = l;

#ifdef SHEARINGBOX
	omega = OMEGAFRAME;
#endif
#ifdef CYLINDRICAL
  omega = sqrt(G*MSTAR/ymed(j)/ymed(j)/ymed(j));
#endif
#ifdef STOKESNUMBER
	alphak  = 0.5*(pref[ll]+pref[lm])*invstokesnumber;
#endif
#ifdef DUSTSIZE
	alphak  = pref[ll]*invparticlesize/rhosolid;
#endif
	sk     = dt*alphak/(1+dt*alphak);

	if (fluidtype == GAS)  {
    temp = (GAMMA-1.0)*energy[ll]/(rho[ll]*R_MU);
    temp=    sumpressure[ll]/( sumrho[ll] );
    energy[ll] = (rho[ll]*R_MU) * temp/(GAMMA - 1.0) 
  }
	else{
    temp = energy[ll] / (rho[ll]*CD);
    temp = sk*sumpressure[ll]/( sumrho[ll] ) + temp/(1.+ dt*alphak);
    energy[ll] = (rho[ll]*CD) * temp; // just doing e = m c_d T for the dust for now to be updated later
  }

//<\#>
#ifdef X
      }
#endif
//<\MAIN_LOOP>
}
