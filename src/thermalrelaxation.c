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
  real invparticlesize = Coeffval[1];
  real rhosolid = Coeffval[2];
  int fluidtype = Fluidtype;
  int ndust = FluidIndex-1;
  int ntab = NTABLE;
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
        if (fluidtype == GAS) {
          beta[ll] = 0.0; 
        }
        else{
          tempdustn = energy[ll] / (dens[ll] * cpdust);

#ifdef DSHARP
	  
	  tdust = tempdustn *TUNITS;

	  if (tdust <= temptab[0]) coolingfactor = dsharp[ndust*ntab + 0];
	  else if (tdust >= temptab[ntab - 1]) coolingfactor = dsharp[ndust*ntab + ntab - 1];
	  else{
	    //Binary search to find the correct temperature bin
	    left = 0;
	    right = ntab - 1;
	    
	    while (right - left > 1) {
	      mid = left + (right - left) / 2;
	      if (tdust >= temptab[mid]) {
	  	left = mid;
	      } else {
	  	right = mid;
	      }
	    }
	    
	    // 3. Linear interpolation
	    offset = ndust*ntab;
	    t0 = temptab[left];
	    t1 = temptab[right];
	    f0 = dsharp[offset+left];
	    f1 = dsharp[offset+right];
	    
	    coolingfactor = f0 + (f1 - f0) * (tdust - t0) / (t1 - t0);
	    
	  }
	  
	  beta[ll] =  coolingfactor/pow(TUNITS,3.0);
	  beta[ll] *= 3.*STEFANK*invparticlesize/(rhosolid*cpdust);
	 
#else
          qlocal = 8.0 * M_PI * KBOLTZ * tempdustn / invparticlesize / PLANCK / C0;
          if (qlocal >= 1.0) {
            beta[ll] = 12.0 * invparticlesize / (rhosolid * cpdust) * STEFANK *
                  pow(tempdustn, 3.0)/10.0;
          } else {
            beta[ll] = M_PI * 120.0 * STEFANK * KBOLTZ /
                  (PLANCK * C0 * rhosolid * cpdust) * pow(tempdustn, 4.0)/10.0;
          }
#endif
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
