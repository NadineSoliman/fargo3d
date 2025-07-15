#include "fargo3d.h"

void Init() {
  int i,j,k;
  real* rho = Density->field_cpu;
  real* e = Energy->field_cpu;
  real* v1 = Vz->field_cpu;
  int dim1 = Nz + 2*NGHZ;
  Coeffval[0]   = 1.0;    
  Coeffval[1]   = 1.0;    
  Coeffval[2]   = 1.0;

  for (i = 0; i<dim1; i++) {
    rho[i] = 1.0;
    v1[i]  = 0.0;
    
    if(Fluidtype==GAS) e[i]   = 2.0/(GAMMA-1.0);
    
    if(Fluidtype==DUST) e[i]   = 1.0;
    
  }
}


void CondInit() {
   Fluids[0] = CreateFluid("gas",GAS);
   SelectFluid(0);
   Init();

    Fluids[1] = CreateFluid("dust",DUST);
   SelectFluid(1);
   Init();
}
