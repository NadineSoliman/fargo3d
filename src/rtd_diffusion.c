//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>

void RTD_SolveDiffusion_cpu(real dt) {
  
//<USER_DEFINED>
  INPUT(Erad);
  INPUT(DiffCoef);
  OUTPUT(Erad);
//<\USER_DEFINED>

//<EXTERNAL>
  real* diff  = DiffCoef->field_cpu;
  real* erad  = Erad->field_cpu;
  real* aux  = Pressure->field_cpu;
  int pitch  = Pitch_cpu;
  int stride = Stride_cpu;
  int size_x = Nx+2*NGHX;
  int size_y = Ny+2*NGHY-1;
  int size_z = Nz+2*NGHZ-1;
//<\EXTERNAL>

//<INTERNAL>
  int i;
  int j;
  int k;
  int ll;
  real c;
  real update;
#ifdef X
  real d1;
  real d2;
  real cxp;
  real cxm;
  int llxm;
  int llxp;
#endif
#ifdef Y
  real cyp;
  real cym;
  real d3;
  real d4;
  int llyp;
  int llym;
#endif
#ifdef Z
  real czp;
  real czm;
  real d5;
  real d6;
  int llzp;
  int llzm;
#endif
  real dxmedm;
  real dxmedp;
  real dxmin;
//<\INTERNAL>

//<CONSTANT>
// real xmin(Nx+2*NGHX+1);
// real ymin(Ny+2*NGHY+1);
// real zmin(Nz+2*NGHZ+1);
// real InvDiffXmed(Nx+2*NGHX);
//<\CONSTANT>

//<MAIN_LOOP>

  i = j = k = 0;

#ifdef Z
  for (k=1; k<size_z; k++) {
#endif
#ifdef Y
    for (j=1; j<size_y; j++) {
#endif
#ifdef X
      for (i=0; i<size_x; i++ ) {
#endif
//<#>
#ifdef X
	ll = l;
	llxm = lxm;
	llxp = lxp;
#endif
#ifdef Y
	llyp = lyp;
	llym = lym;
#endif
#ifdef Z
	llzp = lzp;
	llzm = lzm;
#endif

	update = 0.0;
	c    = erad[ll]; //Cell centered
	
	// FLUX DIFFUSION ALONG X-DIRECTION
	dxmin  = xmin(i+1)-xmin(i);
	dxmedp = 1.0/InvDiffXmed(ixp);
	dxmedm = 1.0/InvDiffXmed(i);
	
	
#ifdef X
        d1   = 0.5*(diff[llxp]+diff[ll]); //face centered in X
        d2   = 0.5*(diff[llxm]+diff[ll]); //face centered in X
	      cxp  = erad[llxp];                                                   //Cell centered
        cxm  = erad[llxm];                                                   //Cell centered
	
	
#ifdef SPHERICAL
        update += 1.0/ymed(j)/ymed(j)/sin(zmed(k))/sin(zmed(k))/(dxmin)*( d1*(cxp-c)/(dxmin) - d2*(c-cxm)/(dxmin) );
#endif
#endif //X
	
	//FLUX DIFFUSION ALONG Y-DIRECTION
#ifdef Y
	d3   = 0.5*(diff[llyp]+diff[ll]);//face centered in Y
	d4   = 0.5*(diff[llym]+diff[ll]);  //face centered in Y
	cyp  = erad[llyp];                                //Cell centered
	cym  = erad[llym];                                //Cell centered
	
	
#ifdef SPHERICAL
	update += 1.0/ymed(j)/ymed(j)/(ymin(j+1)-ymin(j))*( ymin(j+1)*ymin(j+1)*d3*(cyp-c)/(ymed(j+1)-ymed(j  )) -
							    ymin(j  )*ymin(j  )*d4*(c-cym)/(ymed(j  )-ymed(j-1)) );
#endif
#endif //Y

        // FLUX DIFFUSION ALONG Z-DIRECTION
#ifdef Z
	d5   = 0.5*(diff[llzp]+diff[ll]); // face centered in Z
	d6   = 0.5*(diff[llzm]+diff[ll]);   // face centered in Z
	czp  = erad[llzp];                                 // Cell centered
	czm  = erad[llzm];                                 // Cell centered

#ifdef SPHERICAL
	update += 1.0/ymed(j)/ymed(j)/sin(zmed(k))/(zmin(k+1)-zmin(k))*( sin(zmin(k+1))*d5*(czp-c)/(zmed(k+1)-zmed(k  )) -
									 sin(zmin(k  ))*d6*(c-czm)/(zmed(k  )-zmed(k-1)) );
#endif
#endif // Z
      	aux[ll] = erad[ll] + dt*update;
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