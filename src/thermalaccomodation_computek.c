//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>

void ThermalAccomodation_ComputeK(real dt, int option) {

  if (option == 1)
    _ThermalAccomodation_ComputeK(dt, Rkk1);
  if (option == 2)
    _ThermalAccomodation_ComputeK(dt, Rkk2);
  
}

void _ThermalAccomodation_ComputeK_cpu(real dt,Field *K) {

//<USER_DEFINED>
#ifdef THERMALRELAXATION
  INPUT(Betarad);
#endif
  INPUT(Alphacol);
  INPUT(Qvec);
  INPUT(DensStar);
  INPUT(Slope);
  OUTPUT(K);
//<\USER_DEFINED>

//<EXTERNAL>
  real* q = Qvec->field_cpu;
  real* qgas = Fluids[0]->Qvec->field_cpu;
#ifdef THERMALRELAXATION
  real* Betagas = Fluids[0]->Betarad->field_cpu;
  real* Beta = Betarad->field_cpu;
#endif
  real* alpha = Alphacol->field_cpu;
  real* rkk = K->field_cpu;
  real* Bsum = DensStar->field_cpu;
  real* Asum = Slope->field_cpu;
  real* grk    = Gammark->field_cpu; 
  int pitch2d = Pitch2D;  
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
  real k0;
  real beta;
  real betagas;
//<\INTERNAL>

//<CONSTANT>
// real xmin(Nx+2*NGHX+1);
// real ymin(Ny+2*NGHY+1);
// real GAMMA(1);
//<\CONSTANT>


//<MAIN_LOOP>

  i = j =k =0;

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
    

#ifdef THERMALRELAXATION
    betagas = Betagas[ll];
    beta = Beta[ll];
#else   
    betagas = 0.0;
    beta = 0.0;
#endif

    k0 = (Asum[ll] - qgas[ll]*(Bsum[ll]+betagas))/(1 + grk[ll]*dt*(Bsum[ll]+betagas));

    if(fluidtype==GAS) rkk[ll] = k0;
    if(fluidtype==DUST) {
        rkk[ll] =  1.0/(1+grk[ll]*dt*(alpha[ll]+beta));
        rkk[ll] *= (alpha[ll]*qgas[ll] - (alpha[ll]+beta)*q[ll] + alpha[ll]*grk[ll]*dt*k0);
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
