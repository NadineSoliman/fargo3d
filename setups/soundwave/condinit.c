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
  real CS = sqrt(GAMMA * R_MU * TGAS); // Speed of sound in the gas
  real cp_gas = GAMMA * R_MU / (GAMMA - 1.0); // Specific heat capacity of the gas
  real cp_dust = 0.01 * cp_gas; // Specific heat capacity of the dust
  
  //Input parameters
  real Ts[NFLUIDS-1] = { 0.1 };//, 0.21544346900318834 , 0.46415888336127786 , 1.0 };
  real Th[NFLUIDS-1] = { 0.1 };// , 0.21544346900318834 , 0.46415888336127786 , 1.0 };
  real eps[NFLUIDS-1] = { 2.24 }; //, 0.23333333333333334 , 0.3666666666666667 , 0.5 };
  double complex lambda=-0.60418653999843897 - 4.14725080837136062*I;//-0.16757341098197820 - 7.18794806155152788 * I; 
  

 #ifdef CONSTANTTHERMALCOEFF
  if(index > 0 ) Coeffval[1] = 1.0/Th[index-1];
  #endif
  if(index > 0 )  Coeffval[0]   = 1.0/Ts[index-1];

  double complex vdust[NFLUIDS-1];
  double complex rhodust[NFLUIDS-1];
  double complex Tdust[NFLUIDS-1];
  double complex edust[NFLUIDS-1];
  double complex delta_rhog = DENSGR * RHOG;
  double complex vgas       = I* lambda/(kmode)*delta_rhog/RHOG;

  printf("v_gas = %g + %gi\n", creal(vgas), cimag(vgas));
  double complex sum=0.0;
    for(m=0;m<NFLUIDS-1;m++){
      sum += eps[m] / (1.0 + lambda*Ts[m]);
  }
  double complex delta_Tg = - (1.0 - ((GAMMA * (-lambda * lambda) / (CS * kmode * CS * kmode)) *(1 + sum))) * delta_rhog / RHOG * TGAS;
  double complex delta_e = (cp_gas/GAMMA) * ((delta_rhog  * TGAS) + (delta_Tg * RHOG));

  printf("Background = %g \n", (RHOG * TGAS * (cp_gas /GAMMA)));
  printf("Mode = %g\n", mode( creal(delta_e), cimag(delta_e), kmode, 0 ));
  printf("Velocity = %g\n", mode( creal(vgas), cimag(vgas), kmode, 0 ));
  for(m=0;m<NFLUIDS-1;m++){

    Tdust[m]   =  1/ (1 + lambda*Th[m]) * delta_Tg;
    vdust[m]   = vgas/ (1 + lambda *Ts[m]);
    rhodust[m] =  eps[m]/(1.0 +lambda*Ts[m])*delta_rhog;
    edust[m]   = (cp_dust) * (rhodust[m] * TGAS  + Tdust[m] * RHOG * eps[m]);

  }
  

  for (k = 0; k<Nz+2*NGHZ; k++) {
    
    if (index == 0) {
      rho[k] = RHOG + AMPLITUDE*mode( DENSGR, DENSGI, kmode, zmed(k) );
      v[k]   = AMPLITUDE*mode( creal(vgas), cimag(vgas), kmode, zmin(k) );
#ifdef ADIABATIC
      e[k]  = (RHOG * TGAS * (cp_gas /GAMMA));
      e[k] +=  AMPLITUDE*mode( creal(delta_e), cimag(delta_e), kmode, zmed(k) );
      printf("e[k] = %f\n", e[k]);
#else
      e[k]  = CS;
#endif
    }
    else {
      rho[k]  = eps[index-1]*RHOG;
      rho[k] += AMPLITUDE*mode( creal(rhodust[index-1]), cimag(rhodust[index-1]), kmode, zmed(k) );
      v[k]    = AMPLITUDE*mode( creal(vdust[index-1])  , cimag(vdust[index-1])  , kmode, zmin(k) );
#ifdef ADIABATIC
      e[k]   = TDUST * (cp_dust)* eps[index-1] * RHOG;
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