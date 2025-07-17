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
  INPUT(Density);
  OUTPUT(Energy);
//<\USER_DEFINED>

//<EXTERNAL>
  real* sumpressure = Slope->field_cpu;
  real* sumrho = DensStar->field_cpu;
  real* pref = Qs->field_cpu;
  real* energy = Energy->field_cpu;
  real* dens = Density->field_cpu;
  int pitch  = Pitch_cpu;
  int stride = Stride_cpu;
  int size_x = Nx;
  int size_y = Ny+2*NGHY;
  int size_z = Nz+2*NGHZ;
  int fluidtype = Fluidtype;
#ifdef DUSTSIZE
  real invparticlesize = Coeffval[1];
  real rhosolid        = Coeffval[2];
#endif
#ifdef CONSTANTTHERMALCOEFF
  real invthermaltime = Coeffval[1];
#endif
//<\EXTERNAL>

//<INTERNAL>
  int i,j,k;
  int ll;
  real alphak;
  real sk;
  real temp;
//<\INTERNAL>

//<CONSTANT>
// real xmin(Nx+2*NGHX+1);
// real ymin(Ny+2*NGHY+1);
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
#ifdef CONSTANTTHERMALCOEFF
        alphak = pref[ll]*invthermaltime;
#endif
#ifdef DUSTSIZE
  alphak = pref[ll]*invparticlesize/rhosolid/CP_DUST;
#endif

  sk     = dt*alphak/(1+dt*alphak);

	if (fluidtype == GAS)  {
    temp=    sumpressure[ll]/( sumrho[ll] );
    energy[ll] = (dens[ll]*R_MU) * temp/(GAMMA - 1.0); 
  }
	else{
    temp = energy[ll] / (dens[ll]*CP_DUST);
    temp = sk*sumpressure[ll]/( sumrho[ll] ) + temp/(1.+ dt*alphak);
    energy[ll] = (dens[ll]*CP_DUST) * temp; // just doing e = m c_d T for the dust for now to be updated later
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
