//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>

void RTD_StellarFlux_cpu () {

//<USER_DEFINED>
  INPUT(Tau);
  INPUT(Density);
  INPUT(Kappa);
  INPUT(Total_Density);
  INPUT(KappaP);
  OUTPUT(GammaRad);
  real rstar = RSTAR/R0_CGS*R0;
  real tstar = TSTAR*(G*MSTAR/R0/R_MU)/(G_CGS*MSTAR_CGS/R0_CGS/R_MU_CGS);
//<\USER_DEFINED>

//<EXTERNAL>
  real* tau        = Tau->field_cpu;
  real* gammarad   = GammaRad->field_cpu;
  real* dens       = Density->field_cpu;
  real* kappa      = Kappa->field_cpu;
  real* kappaP= KappaP->field_cpu;
  real* totaldens = Total_Density->field_cpu;
  real  lum        = STEFANK * rstar*rstar * tstar*tstar*tstar*tstar;
  int pitch        = Pitch_cpu;
  int stride       = Stride_cpu;
  int size_x       = Nx+2*NGHX;
  int size_y       = Ny+2*NGHY;
  int size_z       = Nz+2*NGHZ;
//<\EXTERNAL>


//<INTERNAL>
  int i;
  int j;
  int k;
  int ll;
//<\INTERNAL>

//<CONSTANT>
// real xmin(Nx+2*NGHX+1);
// real Syj(Ny+2*NGHY);
// real Syk(Nz+2*NGHZ);
// real ymin(Ny+2*NGHY+1);
// real zmin(Nz+2*NGHZ+1);
// real InvVj(Ny+2*NGHY);
// real Sxi(Nx+2*NGHX);
//<\CONSTANT>

  
//<MAIN_LOOP>

  i = j = k = 0;

#ifdef Z
  for (k=0; k<size_z; k++) {
#endif
#ifdef Y
    for (j=1; j<size_y; j++) {
#endif
#ifdef X
      for (i=0; i<size_x; i++ ) {
#endif
//<#>
	ll = l;
	
    //Attenuated flux due to integrated absorption (tau: sum of all dust species)
    gammarad[ll] = lum*SurfY(i,j,k)/(Vol(i,j,k)*ymin(j)*ymin(j))*(exp(-tau[lym])-exp(-tau[ll]));

    //Effective stellar energy flux on a given species
	gammarad[ll] *= kappa[ll]*dens[ll]/(kappaP[ll]*totaldens[ll]);
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
