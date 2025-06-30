#include "fargo3d.h"

void _CondInit(int id) {
  int i,j,k;

  real *vx   = Vx->field_cpu; 
  real *vy   = Vy->field_cpu; 
  real *rho  = Density->field_cpu;
  real *e    = Energy->field_cpu;

  OUTPUT(Density);
  OUTPUT(Energy);
  OUTPUT(Vx);
  OUTPUT(Vy);

  real stokes_plus[NFLUIDS];
  real stokes[NFLUIDS-1];
  real epsilons[NFLUIDS-1];
  real sq = SQ;
  real slope;
  // Stokes numbers                                                                      
  real smax = TSMAX;
  real smin = TSMIN;
  real ds   = (log(smax)-log(smin))/(NFLUIDS-1);
  for(int n=0;n<NFLUIDS;n++){
    stokes_plus[n] = smin*exp(ds*n);
  }

  //Dust size-distribution                                                                 
  slope = 4-sq;
  for(int n=0; n<NFLUIDS-1; n++){
    if(slope != 0.) {
      epsilons[n]  = pow(stokes_plus[n+1],slope) - pow(stokes_plus[n],slope) ;
      epsilons[n] *= EPSILON/(pow(smax, slope) - pow(smin,slope));
    }
    else{
      epsilons[n]  = log(stokes_plus[n+1]/stokes_plus[n]);
      epsilons[n] *= EPSILON/log(smax/smin);
    }
    stokes[n] = stokes_plus[n+1];
    if( NFLUIDS == 2) stokes[n] = TSMAX;
  }

#ifdef DRAGFORCE
  if(id > 0) {
#ifdef STOKESNUMBER
   Coeffval[0]   = 1.0/stokes[id-1];
#endif
#ifdef DUSTSIZE
    Coeffval[1]   = 1.0/(stokes[id-1]*R0/R0_CGS);    
    Coeffval[2]   = RHOSOLID/(MSTAR_CGS/(R0_CGS*R0_CGS*R0_CGS))*(MSTAR/(R0*R0*R0));
#endif
    if(CPU_Master) printf("Ts %1.16f \t eps %1.16f \n", stokes[id-1], epsilons[id-1]);
  }
#endif
  //Initial profiles
  real x,y;
  real fbump, fbumpp;
  real pressure;

  i = j = k = 0;
  
#ifdef Y
  for (j=0; j<Ny+2*NGHY; j++) {
#endif
  y = Ymed(j); // "y" is equivalent to radial (local) coordinate
#ifdef X
    for (i=0; i<Nx+2*NGHX; i++) {
#endif
    x = Xmed(i); // "x" is equivalent to azimuthal (local) coordinate
    fbump  = exp( -0.5*y*y/(GW*GW) ); // Gaussian with width w
    fbumpp = fbump*(-y/(GW*GW)); // derivative of Gaussian

if (Fluidtype == GAS) {
	
    rho[l] = SIGMA0*(1.0 + GA*fbump);
    pressure = CS*CS*rho[l];

    vx[l] = -SHEARPARAM*OMEGAFRAME*y + 0.5*CS*CS*GA*fbumpp/rho[l]/OMEGAFRAME; // vx ix azimuthal velocity
    vy[l] = 0.5*PERT*(sin(2.*M_PI*x/5./H0)-sin(4.*M_PI*x/4./H0))*fbump; // vy is the radial velocity
    
    #ifdef ISOTHERMAL
      e[l] = CS; // isothermal; energy
    #else
      e[l] = pressure/(GAMMA-1.0); 
    #endif
}
if (Fluidtype == DUST) {
	
    rho[l] = SIGMA0*(1.0 + GA*fbump)*epsilons[id-1];
    
    vx[l] = -SHEARPARAM*OMEGAFRAME*y; // vx ix azimuthal velocity
    vy[l] = 0.0;
    e[l] = 0.0;
}

#ifdef X
      }
#endif
#ifdef Y
    }
#endif
}

void CondInit() {
  
  int id_gas = 0;
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

