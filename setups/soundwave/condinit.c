#include "fargo3d.h"
#include <complex.h>

real mode(real fr, real fi, real kmode, real x) {
  return ( fr*cos(kmode*x) - fi*sin(kmode*x) );
}

void _CondInit(int index) {

  int feedback = YES;
  int k,m;
  
  real* rho   = Density->field_cpu;
  real* cs    = Energy->field_cpu;
  real* v     = Vz->field_cpu;
  real kmode  = 2.*M_PI/(ZMAX-ZMIN);
 
  
  //Input parameters
  real Ts[NFLUIDS-1] = { 0.1 , 0.21544346900318834 , 0.46415888336127786 , 1.0 };
  real eps[NFLUIDS-1] = { 0.1 , 0.23333333333333334 , 0.3666666666666667 , 0.5 };
  double complex lambda =  0.9124135035415595 + -5.49379966793634 *I;
  
  double complex vdust[NFLUIDS-1];
  double complex rhodust[NFLUIDS-1];
  double complex delta_rhog = DENSGR;
  double complex vgas       = -I*lambda/kmode*delta_rhog/RHOG;

  for(m=0;m<NFLUIDS-1;m++){
    vdust[m]   =  -I*lambda/( kmode*(1.0-lambda*Ts[m]) ) *delta_rhog/RHOG;
    rhodust[m] =  eps[m]/(1.0-lambda*Ts[m])*delta_rhog;
    //printf("Dust %d: v=%.3f+%.3fi, rho=%.3f+%.3fi\n", m, creal(vdust[m]), cimag(vdust[m]), creal(rhodust[m]), cimag(rhodust[m]));
  }
  
  if(index > 0 )    Coeffval[0]   = 1.0/Ts[index-1];
  for (k = 0; k<Nz+2*NGHZ; k++) {
    
    if (index == 0) {
      rho[k] = RHOG + AMPLITUDE*mode( DENSGR, DENSGI, kmode, zmed(k) );
      v[k]   = AMPLITUDE*mode( creal(vgas), cimag(vgas), kmode, zmin(k) );
      cs[k]  = CS;
    }
    else {
      rho[k]  = eps[index-1]*RHOG;
      rho[k] += AMPLITUDE*mode( creal(rhodust[index-1]), cimag(rhodust[index-1]), kmode, zmed(k) );
      v[k]    = AMPLITUDE*mode( creal(vdust[index-1])  , cimag(vdust[index-1])  , kmode, zmin(k) );
      cs[k]   = 0.0;
    }
  }
}


void CondInit() {
  
  int id_gas = 0;
  int feedback = YES;
  //We first create the gaseous fluid and store it in the array Fluids[]
  Fluids[id_gas] = CreateFluid("gas",GAS);

  //We now select the fluid
  SelectFluid(id_gas);

  //and fill its fields
  _CondInit(0);

  //We repeat the process for the dust fluids
  char dust_name[MAXNAMELENGTH];
  int id_dust;

  for(id_dust = 1; id_dust<NFLUIDS; id_dust++) {
    sprintf(dust_name,"dust%d",id_dust); //We assign different names to the dust fluids

    Fluids[id_dust]  = CreateFluid(dust_name, DUST);
    SelectFluid(id_dust);
    _CondInit(id_dust);

  }

}
