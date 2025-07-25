#include "fargo3d.h"

void _CondInit(int id) {
  int i,j,k;
  real *v1;
  real *v2;
  real *v3;
  real *e;
  real *rho;
  
  real omega;
  real r, r3;
  real soundspeed;
  rho = Density->field_cpu;
  e   = Energy->field_cpu;
  v1  = Vx->field_cpu;
  v2  = Vy->field_cpu;
  v3  = Vz->field_cpu;

  real stokes_plus[NFLUIDS];
  real stokes[NFLUIDS-1];
  real epsilons[NFLUIDS-1];
  real sq = SQ;
  real slope;
  // Stokes numbers                                                                      
  real smax = TSMAX;
  real smin = TSMIN;
  real cv    = 1./(GAMMA-1.0);
  real cdust=CPDG*GAMMA/(GAMMA-1.0); 
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

  for (k=0; k<Nz+2*NGHZ; k++) {
    for (j=0; j<Ny+2*NGHY; j++) {
      r = Ymed(j);
      r3 = r*r*r;
      omega = sqrt(G*MSTAR/(r3));
      for (i=0; i<Nx; i++) {
	  v2[l] = v3[l] = 0.0;
	  v1[l] = omega*r;

	  real xi = SIGMASLOPE+1.+FLARINGINDEX;
	  real beta = 1.-2*FLARINGINDEX;
	  real h = ASPECTRATIO*pow(r/R0,FLARINGINDEX);
	  if (FLARINGINDEX == 0.0) {
	   rho[l] = SIGMA0/sqrt(2.0*M_PI)/(R0*ASPECTRATIO)*pow(r/R0,-xi)* \
	    pow(sin(Zmed(k)),-beta-xi+1./(h*h));
	  } else {
	    rho[l] = SIGMA0/sqrt(2.0*M_PI)/(R0*ASPECTRATIO)*pow(r/R0,-xi)* \
	    pow(sin(Zmed(k)),-xi-beta)*					\
	    exp((1.-pow(sin(Zmed(k)),-2.*FLARINGINDEX))/2./FLARINGINDEX/(h*h));
	  }
	
    if(Fluidtype==DUST)       rho[l]  *= epsilons[id-1];

  #ifdef ISOTHERMAL
	  e[l] = h*sqrt(G*MSTAR/r);
    if(Fluidtype==DUST) e[l] = 0.;
  #else
	  e[l] = cv*rho[l]*h*h*G*MSTAR/r;
    if(Fluidtype==DUST) e[l]    *= 0.1*CPDG*GAMMA; //Dust energy assuming  Td=Tg
  #endif

  if(Fluidtype==GAS) v1[l] *= sqrt(pow(sin(Zmed(k)),-2.*FLARINGINDEX)-(beta+xi)*h*h);

  //Frame rotation
	v1[l] -= OMEGAFRAME*r*sin(Zmed(k));

  //Initial noise
  soundspeed  = h*sqrt(G*MSTAR/r);
  v2[l]+= soundspeed*NOISE*(drand48()-.5);
  v3[l]+= soundspeed*NOISE*(drand48()-.5);



      }
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
