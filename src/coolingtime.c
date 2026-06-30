//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>

void CoolingTime_cpu(real dt, real invparticlesize, real rhosolid, real eps, int nd) {

//<USER_DEFINED>
  INPUT(Energy);
  INPUT(Density);
  INPUT(Betarad);
  OUTPUT(Betarad);
#ifdef GPU
  real* temptabpointer = TempTable_d;
  real* dsharppointer  = Dsharp_d;
#else
  real* temptabpointer = TempTable;
  real* dsharppointer  = Dsharp;
#endif
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
  real* temptab = temptabpointer;
  real* dsharp  = dsharppointer;
  real cpdg = CPDG;
  int ntab = NTABLE;
  real dlog_T = (log10(TMAXTAB) - log10(TMINTAB)) / (NTABLE - 1);
  real logT_min = log10(TMINTAB);
  real tmintab = TMINTAB;
  real tmaxtab = TMAXTAB;
//<\EXTERNAL>
  
//<INTERNAL>
  int i;
  int j;
  int k;
  int ll;
  int left;
  int right;
  int mid;
  int offset;
  real qlocal;
  real tempdustn;
  real cpgas;
  real cpdust;
  real tdust;
  real t0;
  real t1;
  real f0;
  real f1;
  real coolingfactor;
  real betar;
  real alphar;
  real tempgas;
  real inv_dlog_T;
  real logT;
//<\INTERNAL>

//<CONSTANT>
// real xmin(Nx+2*NGHX+1);
// real ymin(Ny+2*NGHY+1);
// real GAMMA(1);
//<\CONSTANT>

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
        cpgas      = R_MU/(GAMMA-1.0);
        cpdust     = cpdg* cpgas;

	tempgas = energy[ll]/(dens[ll]*cpgas);

	tempdustn = tempgas;

 
	//radiative cooling

#ifdef DSHARP
	
	  tdust = tempdustn *TUNITS;
    inv_dlog_T = 1.0 / dlog_T;
	  
    if (tdust <= tmintab) {
      coolingfactor = dsharp[nd * ntab + 0];
    } else if (tdust >= tmaxtab) {
      coolingfactor = dsharp[nd * ntab + ntab - 1];
    } 
    else {
      logT = log10(tdust);
      left = (int)((logT - logT_min) * inv_dlog_T);
      right = left + 1;
     
      // Linear interpolation
      t0 = temptab[left];
      t1 = temptab[right];
      f0 = dsharp[nd * ntab + left];
      f1 = dsharp[nd * ntab + right];

      coolingfactor = f0 + (f1 - f0) * (tdust - t0) / (t1 - t0);
	  }
	  
	  betar =  coolingfactor/pow(TUNITS,3.0);
	  betar *= 3.*STEFANK*invparticlesize/(rhosolid*cpdust);
	 
#else
          qlocal = 8.0 * M_PI * KBOLTZ * tempdustn / invparticlesize / PLANCK / C0;
          if (qlocal >= 1.0) {
            betar = 12.0 * invparticlesize / (rhosolid * cpdust) * STEFANK *
                  pow(tempdustn, 3.0)/10.0;
          } else {
            betar = M_PI * 120.0 * STEFANK * KBOLTZ /
                  (PLANCK * C0 * rhosolid * cpdust) * pow(tempdustn, 4.0)/10.0;
          }
#endif
     

      //Coolisional cooling
      alphar  = 0.75 *pow(KBOLTZ/MH, 1.5)* pow(tempgas, 0.5) * dens[ll];
      alphar  *= invparticlesize/(rhosolid*cpdust);

      //Total
      beta[ll] += cpdg*eps*alphar*betar/(alphar+betar);

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
