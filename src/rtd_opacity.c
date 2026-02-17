//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>

void RTD_Opacity_cpu () {

//<USER_DEFINED>
  INPUT(Temperature);
  INPUT(Density);
  INPUT(Total_Density);
  INPUT(KappaP);
  OUTPUT(KappaP);
  OUTPUT(Kappa);
//<\USER_DEFINED>


//<EXTERNAL>
  real* kappa= Kappa->field_cpu;
  real* kappaP= KappaP->field_cpu;
  real* temp = Temperature->field_cpu;
  real* dens = Density->field_cpu;
  real* totaldens = Total_Density->field_cpu;
  int pitch  = Pitch_cpu;
  int stride = Stride_cpu;
  int size_x = Nx+2*NGHX;
  int size_y = Ny+2*NGHY;
  int size_z = Nz+2*NGHZ;
  real invparticlesize = Coeffval[1];
  real rhosolid        = RHOSOLID;
  int fluidtype = Fluidtype;
//<\EXTERNAL>

//<INTERNAL>
  int i;
  int j;
  int k;
  int ll;
  real qeff;
  real size;
  real kext;
  real tempu;
//<\INTERNAL>
  
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
    if(fluidtype==DUST){
        size = R0_CGS/R0/invparticlesize; //cm
        tempu=temp[ll]/(G*MSTAR/(R0*R_MU))*(G_CGS*MSTAR_CGS/(R0_CGS*R_MU_CGS)); //K

	      qeff    = MIN(1.0, tempu*size/6.0e-2 );
	      kext    = 3./4.*qeff/(size*rhosolid);
        // cross-section per unit mass of the dust species
        kappa[ll] = kext/(R0_CGS*R0_CGS/MSTAR_CGS)*R0*R0/MSTAR;
        // opacity
    }
    else{
        kappa[ll] = 1.0e-4/(R0_CGS*R0_CGS/MSTAR_CGS)*R0*R0/MSTAR;
    }

    // Mass-average total opacity
    kappaP[ll] += kappa[ll]*dens[ll]/totaldens[ll];
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
