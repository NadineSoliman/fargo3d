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
  OUTPUT(Betarad);
//<\USER_DEFINED>

//<EXTERNAL>
  real* energy = Energy->field_cpu;
  real* dens = Density->field_cpu;
  real* beta = Betarad->field_cpu;
  int size_x = Nx;
  int size_y = Ny + 2 * NGHY;
  int size_z = Nz + 2 * NGHZ;
  int pitch  = Pitch_cpu;
  int stride = Stride_cpu;
  real cpdg = CPDG;
  real invparticlesize = Coeffval[1];
  real rhosolid = Coeffval[2];
  int fluidtype = Fluidtype;
//<\EXTERNAL>

//<INTERNAL>
  int i;
  int j;
  int k;
  int ll;
  real cpdust;
  real cpgas;
  real qlocal;
  real tempdustn;
//<\INTERNAL>

//<MAIN_LOOP>

  i = j = k = 0;

#ifdef Z
  for (k = 0; k < size_z; k++) {
#endif
#ifdef Y
    for (j = 0; j < size_y; j++) {
#endif
#ifdef X
      for (i = 0; i < size_x; i++) {
#endif
//<#>
        ll = l;
        if (fluidtype == GAS) {
          beta[ll] = 0.0; 
        }
        else{
          cpgas = R_MU / (GAMMA - 1.0);
          cpdust = cpdg * cpgas;

          tempdustn = energy[ll] / (dens[ll] * cpdust);
          qlocal = 8.0 * M_PI * KBOLTZ * tempdustn / invparticlesize / PLANCK / C0;
          if (qlocal >= 1.0) {
            beta[ll] = 12.0 * invparticlesize / (rhosolid * cpdust) * STEFANK *
                  pow(tempdustn, 3.0);
          } else {
            beta[ll] = M_PI * 120.0 * STEFANK * KBOLTZ /
                  (PLANCK * C0 * rhosolid * cpdust) * pow(tempdustn, 4.0);
          }
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
