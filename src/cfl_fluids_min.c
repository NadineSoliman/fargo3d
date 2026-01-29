#include "fargo3d.h"

void CflFluidsMin() {
  int i;
  real step = 1e30;
  real min;
  real stepcol = 1e30;
  real mincol;
  real steprt = 1e30;
  real minrt;

  for (i=0;i<NFLUIDS;i++) {
    if (step > Min[i])
      step = Min[i];
  }
  
  for (i=0;i<NFLUIDS;i++) {
    if (stepcol > Mincol[i])
      stepcol = Mincol[i];
  }
   
   minrt = reduction_full_MIN(Qs, NGHY, Ny+NGHY, NGHZ, Nz+NGHZ);
   steprt = minrt;

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


#ifdef FLOAT
  MPI_Allreduce(&steprt, &StepTimeRT, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
#else
  MPI_Allreduce(&steprt, &StepTimeRT, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
#endif


#ifdef STS
  StepTimeSTS = StepTimeRT*NSTS*NSTS/NUSTS;
  if(StepTimeSTS > StepTime) printf("STS error: Number of substeps NSTS must be reduced");
  else: StepTime = StepTimeSTS;
#endif
}
