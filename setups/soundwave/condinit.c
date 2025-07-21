#include "fargo3d.h"
#include <complex.h>

real mode(real fr, real fi, real kmode, real x) {
  return ( fr*cos(kmode*x) - fi*sin(kmode*x) );
}

void _CondInit(int index) {

  int feedback = YES;
  int k,m;
  
  real* rho   = Density->field_cpu;
  real* e    = Energy->field_cpu;
  real* v     = Vz->field_cpu;
  real kmode  = 2.*M_PI/(ZMAX-ZMIN);


  //Input parameters
  real Ts[NFLUIDS-1] = { 0.1 };//, 0.21544346900318834 , 0.46415888336127786 , 1.0 };
  real Th[NFLUIDS-1] = { 0.1};// , 0.21544346900318834 , 0.46415888336127786 , 1.0 };
  real eps[NFLUIDS-1] = { 0.1}; //, 0.23333333333333334 , 0.3666666666666667 , 0.5 };
  double complex lambda =  -0.18186145573318502 - 7.16601928447143877*I; 
  
 #ifdef CONSTANTTHERMALCOEFF
  if(index > 0 ) Coeffval[1] = 1.0/Th[index-1];
  #endif
  if(index > 0 )  Coeffval[0]   = 1.0/Ts[index-1];

  double complex vdust[NFLUIDS-1];
  double complex rhodust[NFLUIDS-1];
  double complex Tdust[NFLUIDS-1];
  double complex edust[NFLUIDS-1];
  double complex delta_rhog = DENSGR * RHOG;
  double complex vgas       = -I*lambda/(CS*kmode)*delta_rhog/RHOG * CS;

  double complex sum=0.0;
    for(m=0;m<NFLUIDS-1;m++){
      sum += eps[m] / (1.0 + lambda*Ts[m]);
  }
  double complex delta_Tg = - (1 - (((lambda * lambda) / (CS * kmode * CS * kmode)) *(1 + sum))) * delta_rhog / RHOG * TGAS;
  double complex delta_e = CP_GAS * (delta_rhog  * TGAS) + (delta_Tg * RHOG);

  for(m=0;m<NFLUIDS-1;m++){

    Tdust[m]   =  - 1/ (1 + lambda*Th[m]) * delta_Tg;
    vdust[m]   = vgas/ (1 + lambda *Ts[m]);
    rhodust[m] =  eps[m]/(1.0+lambda*Ts[m])*delta_rhog;
    edust[m]   = CP_DUST * (rhodust[m] * TGAS  + Tdust[m] * RHOG * eps[m]);

  }
  

  for (k = 0; k<Nz+2*NGHZ; k++) {
    
    if (index == 0) {
      rho[k] = RHOG + AMPLITUDE*mode( DENSGR, DENSGI, kmode, zmed(k) );
      v[k]   = AMPLITUDE*mode( creal(vgas), cimag(vgas), kmode, zmin(k) );
#ifdef ADIABATIC
      e[k]  = (RHOG * TGAS * CP_GAS);
      e[k] +=  AMPLITUDE*mode( creal(delta_e), cimag(delta_e), kmode, zmin(k) );
#else
      e[k]  = CS;
#endif
    }
    else {
      rho[k]  = eps[index-1]*RHOG;
      rho[k] += AMPLITUDE*mode( creal(rhodust[index-1]), cimag(rhodust[index-1]), kmode, zmed(k) );
      v[k]    = AMPLITUDE*mode( creal(vdust[index-1])  , cimag(vdust[index-1])  , kmode, zmin(k) );
#ifdef ADIABATIC
      e[k]   = TDUST * CP_DUST * eps[index-1] * RHOG;
      e[k]  += AMPLITUDE*mode( creal(edust[index-1]), cimag(edust[index-1]), kmode, zmed(k) );
#else
      e[k]   =0;
#endif
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