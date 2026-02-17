//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>


void RTD_MatterRadiation_UpdateTemp_cpu(real dt) {

//<USER_DEFINED>
  INPUT(Density);
  INPUT(Energy);
  INPUT(DensStar);
  INPUT(Temperature);
  INPUT(Kappa);
  INPUT(GammaRad);
  INPUT(Energyrad);
  OUTPUT(DensStar);
//<\USER_DEFINED>

//<EXTERNAL>
  real* dens = Density->field_cpu;
  real* temp = Temperature->field_cpu;
  real* kappa = Kappa->field_cpu;
  real* erad = Energyrad->field_cpu;
  real* heat  = GammaRad->field_cpu;
  real* sum  = DensStar->field_cpu;
  real* energy = Energy->field_cpu;
  int pitch  = Pitch_cpu;
  int stride = Stride_cpu;
  int size_x = Nx;
  int size_y = Ny+2*NGHY;
  int size_z = Nz+2*NGHZ;
  int fluidtype = Fluidtype;
  real cpdg=CPDG;
//<\EXTERNAL>

//<INTERNAL>
  int i;
  int j;
  int k;
  int ll;
  real cpgas;
  real cpdust;
  real beta;
  real xi;
  real cv;
//<\INTERNAL>

//<CONSTANT>
// real xmin(Nx+2*NGHX+1);
// real ymin(Ny+2*NGHY+1);
// real GAMMA(1);
//<\CONSTANT>


//<MAIN_LOOP>

  i =j =k =0;
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
      
    cpgas  = GAMMA*R_MU/(GAMMA-1.0);
    cpdust = cpdg* cpgas;

    cv = cpgas/GAMMA;
    if(fluidtype==DUST) cv = cpdust;

    xi = 16.*dt*kappa[ll]*STEFANK*pow(temp[ll],3)/(cv);

	temp[ll] += dt* ( heat[ll] + dens[ll]*kappa[ll]*(C0*erad[ll]- 4.*STEFANK*pow(temp[ll],4)))/(dens[ll]*cv*(1+xi));

  //Update energy
  energy[ll] = cv*dens[ll]*temp[ll]; 
  
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