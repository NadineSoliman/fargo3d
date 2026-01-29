//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>

//Radiative Transfer Dust module
void RTD_main(real dt) {


    real dtsts;
    int n;
        FARGO_SAFE(FillGhosts(ENERGYRAD));

    (WriteMerging(Erad,0));

    for(n=1;n<=NSTS;n++){

        dtsts = dt/((NUSTS-1)*cos((2.*n-1.)*M_PI*0.5/NSTS)+NUSTS+1);
        printf("%d %1.8f \n",n,dtsts);
        //Compute Temperature
        MULTIFLUID(RTD_Temperature());

        //Compute Opacity of each species
        Reset_field(KappaP);  
        MULTIFLUID(RTD_Opacity());

        //Compute Optical depth
        Reset_field(Tau);  
        MULTIFLUID(RTD_Tau());
        RTD_Scan(); //not suitable for domain decompsition in j-direction

        //Compute Stellar flux
        //MULTIFLUID(RTD_StellarFlux());

        //For debugging
        MULTIFLUID(WriteMerging(GammaRad,0));



        //--------------------------------------------------------------------------------------------------------------
        //Radiative Energy (1) Solution matter-radiation coupling
        //summation A
        Reset_field(DensStar); 
        MULTIFLUID(RTD_MatterRadiationSumA(dtsts));
        //summation B
        Reset_field(Slope);
        MULTIFLUID(RTD_MatterRadiationSumB(dtsts));

        //Update Radiative Energy
        FARGO_SAFE(RTD_MatterRadiation_UpdateErad());
        //Update Temperature and Energy
        MULTIFLUID(RTD_MatterRadiation_UpdateTemp(dtsts));

        FARGO_SAFE(FillGhosts(ENERGYRAD));
 
        (WriteMerging(Erad,111));

        //--------------------------------------------------------------------------------------------------------------
        //Radiative Energy (2) Diffusion with limmitter and reduced speed-of-light
        //Solve diffusion equation with super-time-stepping to update radiative energy
        FARGO_SAFE(RTD_DiffusionCoeff());
        FARGO_SAFE(RTD_SolveDiffusion(dtsts));
        FARGO_SAFE(copy_field(Erad,Pressure));
        FARGO_SAFE(FillGhosts(ENERGYRAD));

        MULTIFLUID(WriteMerging(Erad,999));

    }
}