//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>

//Radiative Transfer Dust module
void RTD_main(real dt) {

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
    MULTIFLUID(RTD_StellarFlux());

    //For debugging
    //MULTIFLUID(WriteMerging(GammaRad,0));

    //Radiative Energy (1) Solution matter-radiation coupling
    // Add stellar heating to dust energy equation (and viscous heating for gas)
    // Compute radiative cooling
    // Solve coupled quartic

    //Radiative Energy (2) Diffusion with limmitter and reduced speed-of-light
    //Solve diffusion equation with super-time-stepping to update radiative energy

}