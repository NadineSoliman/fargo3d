#include "fargo3d.h"

void CflFluidsMin() {
  int i;
  real step = 1e30;
  real min;
  real stepcol = 1e30;
  real mincol;

  for (i=0;i<NFLUIDS;i++) {
    if (step > Min[i])
      step = Min[i];
  }
  
  for (i=0;i<NFLUIDS;i++) {
    if (stepcol > Mincol[i])
      stepcol = Mincol[i];
  }
#ifdef FLOAT
  MPI_Allreduce(&step, &StepTime, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
#else
  MPI_Allreduce(&step, &StepTime, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
#endif
  if(StepTime <= SMALLTIME) {
    masterprint("Error!!!--> Null dt\n");
    prs_exit(1);
  }

#ifdef FLOAT
  MPI_Allreduce(&stepcol, &StepTimeCol, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
#else
  MPI_Allreduce(&stepcol, &StepTimeCol, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
#endif

}
