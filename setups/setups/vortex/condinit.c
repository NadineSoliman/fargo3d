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


    rho[l] = SIGMA0*(1.0 + GA*fbump);
    pressure = CS*CS*rho[l];

    vx[l] = -SHEARPARAM*OMEGAFRAME*y + 0.5*CS*CS*GA*fbumpp/rho[l]/OMEGAFRAME; // vx ix azimuthal velocity
    vy[l] = 0.5*PERT*(sin(2.*M_PI*x/5./H0)-sin(4.*M_PI*x/4./H0))*fbump; // vy is the radial velocity
    
    #ifdef ISOTHERMAL
      e[l] = CS; // isothermal; energy
    #else
      e[l] = pressure/(GAMMA-1.0); 
    #endif

#ifdef X
      }
#endif
#ifdef Y
    }
#endif
}

void CondInit() {
   Fluids[0] = CreateFluid("gas",GAS);
   SelectFluid(0);
   _CondInit(0);
}
