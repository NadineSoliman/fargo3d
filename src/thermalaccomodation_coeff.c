//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>

void ThermalAccomodation_Coeff_cpu(real dt) {

//<USER_DEFINED>
  INPUT(Density);
  INPUT(Energy);
  OUTPUT(Alphacol);
#ifdef THERMALRELAXATION
  OUTPUT(Betarad);
#endif
  OUTPUT(Gammark);
//<\USER_DEFINED>

//<EXTERNAL>
  real* alpha      = Alphacol->field_cpu;
#ifdef THERMALRELAXATION
  real* beta       = Betarad->field_cpu;
#endif
  real* dens_gas   = Fluids[0]->Density->field_cpu;
  real* energy_gas = Fluids[0]->Energy->field_cpu; 
  real* dens = Density->field_cpu;
  real* energy = Energy->field_cpu;
  real* grk   = Gammark->field_cpu;
  int pitch  = Pitch_cpu;
  int stride = Stride_cpu;
  int size_x = Nx;
  int size_y = Ny+2*NGHY;
  int size_z = Nz+2*NGHZ;
  real invparticlesize = Coeffval[1];
  real rhosolid        = Coeffval[2];
  int fluidtype = Fluidtype;
  real cpdg=CPDG;
  real gammark2=GAMMARK2;
//<\EXTERNAL>

//<INTERNAL>
  int i;
  int j;
  int k;
  int ll;
  real tempgas;
  real tempdustn;
  real cpdust;
  real cpgas;
  real qlocal;  
  real omega;
//<\INTERNAL>

//<CONSTANT>
// real xmin(Nx+2*NGHX+1);
// real ymin(Ny+2*NGHY+1);
// real GAMMA(1);
//<\CONSTANT>
  
//<MAIN_LOOP>

  i = j = k = 0;

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
	 omega     = sqrt(G*MSTAR/ymed(j)/ymed(j)/ymed(j));
   cpgas      = GAMMA*R_MU/(GAMMA-1.0);
   cpdust    = cpdg* cpgas;

  if (fluidtype == GAS) {
    alpha[ll] = 0.0; //gas
    grk[ll]   = gammark2;
  }
  else {
#ifdef CONSTANTTHERMALCOEFF
      alpha[ll] = 1.0*invthermaltime;
#endif
#ifdef DUSTSIZE
  tempgas    =  (GAMMA-1.0)*energy_gas[ll]/(dens_gas[ll]*R_MU);
  alpha[ll]  = 0.75 *pow(KBOLTZ/MH, 1.5)* pow(tempgas, 0.5) * dens_gas[ll];
  alpha[ll] *= invparticlesize/rhosolid/cpdust;
  grk[ll]    = max2(grk[ll],gammark2);
  if(alpha[ll]*dt>1.0 && grk[ll] != 1.0) grk[ll]=0.5; 
#endif	
  }
#ifdef THERMALRELAXATION
  if (fluidtype == GAS) beta[ll] = 0.0;
  else {
   tempdustn = energy[ll] / (dens[ll]*cpdust);
   qlocal    = 8.0 * M_PI   * KBOLTZ * tempdustn/ invparticlesize / PLANCK / C0;
   if (qlocal >= 1.0) beta[ll] = 12.0 * invparticlesize/ ( rhosolid * cpdust )* STEFANK * pow(tempdustn, 3.0);
   else beta[ll] = M_PI  * 120.0 * STEFANK * KBOLTZ / (PLANCK * C0 * rhosolid * cpdust) * pow(tempdustn, 4.0);
  }
#endif
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
