#include "fargo3d.h"

// --- TEMPERATURE PROVIDER FUNCTIONS ---

/**
 * Returns the gas temperature based on the selected physics (Isothermal vs Gradient)
 */
static real get_radiation_temperature(real r, int k) {
#ifdef TEMPGRAD
    real h = ASPECTRATIO * pow(r / R0, FLARINGINDEX);
    real Tunit = (TCMB / 2.73);
    real Tatm = 80.0 * pow(r, -1) * Tunit;
    real Tmid = h * h * G * MSTAR / r;
    real zq = 3.0;
    real z = cos(Zmed(k)) / ASPECTRATIO;
    return (Tmid + (Tatm - Tmid) * pow(cos(M_PI / 2.0 * z / zq + M_PI / 2.0), 4));
#else
    real h = ASPECTRATIO * pow(r / R0, FLARINGINDEX);
    return h * h * G * MSTAR / r;
#endif
}

/**
 * Returns the specific temperature of a dust grain based on its size (Stokes number)
 * and the local radiation regime.
 */
static real get_dust_temperature(real teff, real stokes_val) {
#ifdef DIFFTEMP
    real x = 2.0 * M_PI * stokes_val;
    real lambda_irr =0.01;//PLANCK * C0 / (4.0 * (teff) * KBOLTZ * R0 / R0_CGS); // cm UV/optical
    real grain_size_cgs = stokes_val * R0 / R0_CGS;
    printf("l=%g, g = %g \n", lambda_irr, x);
    // Regime 1: Small grains
    if (x < lambda_irr) {
        return pow(PLANCK * C0 / (4.0 * (lambda_irr * R0 / R0_CGS) * KBOLTZ), 0.2) * pow(teff, 0.8);
    } 
    // Regime 2: Intermediate grains
    else if (grain_size_cgs < PLANCK * C0 / (4.0 * KBOLTZ * teff)) {
        return pow(PLANCK * C0 / (8.0 * M_PI * grain_size_cgs * KBOLTZ), 0.2) * pow(teff, 0.8);
    } 
#endif
    // Regime 3 / Default: Large grains or no DIFFTEMP
    return teff;
}

// --- SYNCHRONIZATION ---

static void sync_dust_energy_to_gas(void) {
    real cv = R_MU / (GAMMA - 1.0);
    real cdust = CPDG  *cv;
    real *rho_g = Fluids[0]->Density->field_cpu;
    real *e_g   = Fluids[0]->Energy->field_cpu;
    int i, j, k;

    for (k = 0; k < Nz + 2 * NGHZ; k++) {
        for (j = 0; j < Ny + 2 * NGHY; j++) {
            for (i = 0; i < Nx + 2 * NGHX; i++) {
                int ll = l;                 
                real Tgas = e_g[ll] / (rho_g[ll] * cv);
                for (int id = 1; id < NFLUIDS; id++) {
                    real *rho_d = Fluids[id]->Density->field_cpu;
                    real *e_d   = Fluids[id]->Energy->field_cpu;
                    e_d[ll] = rho_d[ll] * cdust * Tgas;
                }
            }
        }
    }
}

// --- CORE INITIALIZATION ---

void _CondInit(int id) {
    int i, j, k;
    real *v1 = Vx->field_cpu;
    real *v2 = Vy->field_cpu;
    real *v3 = Vz->field_cpu;
    real *rho = Density->field_cpu;
    real *e = Energy->field_cpu;
    real *gas_energy = Fluids[0]->Energy->field_cpu;
    real *gas_dens = Fluids[0]->Density->field_cpu;

    real omega, r, r3, soundspeed;
    real stokes_plus[NFLUIDS];
    real stokes[NFLUIDS - 1];
    real epsilons[NFLUIDS - 1];
    real smax = TSMAX;
    real smin = TSMIN;
    real cv = R_MU / (GAMMA - 1.0);
    real cdust = CPDG * cv;
    real denom = 0.0;
    real ds = (log(smax) - log(smin)) / (NFLUIDS - 1);

    // 1. Calculate size distribution and Stokes numbers
    for (int n = 0; n < NFLUIDS; n++) stokes_plus[n] = smin * exp(ds * n);

    real slope = 4.0 - SQ;
    for (int n = 0; n < NFLUIDS - 1; n++) {
        if (slope != 0.0) {
            epsilons[n] = (pow(stokes_plus[n+1], slope) - pow(stokes_plus[n], slope)) * (EPSILON / (pow(smax, slope) - pow(smin, slope)));
        } else {
            epsilons[n] = log(stokes_plus[n+1] / stokes_plus[n]) * (EPSILON / log(smax / smin));
        }
        stokes[n] = sqrt(stokes_plus[n] * stokes_plus[n+1]);
        if (NFLUIDS == 2) {
            stokes[n] = TSMAX;
            epsilons[n] = EPSILON;
        }
        denom += epsilons[n] / (stokes[n]);
    }

#ifdef DRAGFORCE
    if (id > 0) {
        #ifdef STOKESNUMBER
        Coeffval[0] = 1.0 / stokes[id-1];
        #endif
        #ifdef DUSTSIZE
        Coeffval[1] = 1.0 / (stokes[id-1] * R0 / R0_CGS);
        Coeffval[2] = RHOSOLID / (MSTAR_CGS / (R0_CGS * R0_CGS * R0_CGS)) * (MSTAR / (R0 * R0 * R0));
        #endif
        if (CPU_Master) printf("Fluid ID %d: Ts %f, eps %f\n", id, stokes[id-1], epsilons[id-1]);
    }
#endif

    for (k = 0; k < Nz + 2 * NGHZ; k++) {
        for (j = 0; j < Ny + 2 * NGHY; j++) {
            r = Ymed(j);
            r3 = r * r * r;
            omega = sqrt(G * MSTAR / r3);
            real h = ASPECTRATIO * pow(r / R0, FLARINGINDEX);
            real trad = get_radiation_temperature(r, k);

            for (i = 0; i < Nx; i++) {
                int ll = l; 

                // --- DENSITY ---
                real xi = SIGMASLOPE + 1. + FLARINGINDEX;
                real beta = 1. - 2 * FLARINGINDEX;
                if (FLARINGINDEX == 0.0) {
                    rho[ll] = SIGMA0 / sqrt(2.0 * M_PI) / (R0 * ASPECTRATIO) * pow(r / R0, -xi) *
                             pow(sin(Zmed(k)), -beta - xi + 1. / (h * h));
                } else {
                    rho[ll] = SIGMA0 / sqrt(2.0 * M_PI) / (R0 * ASPECTRATIO) * pow(r / R0, -xi) *
                             pow(sin(Zmed(k)), -xi - beta) *
                             exp((1. - pow(sin(Zmed(k)), -2. * FLARINGINDEX)) / 2. / FLARINGINDEX / (h * h));
                }
                if (Fluidtype == DUST) rho[ll] *= epsilons[id-1];

                // --- ENERGY ---
                #ifdef ISOTHERMAL
                    e[ll] = (Fluidtype == GAS) ? h * sqrt(G * MSTAR / r) : 0.0;
                #else
                    if (Fluidtype == DUST) {
                        #ifdef DIFFTEMP
                        real tdust = get_dust_temperature(trad, stokes[id-1]);
                        real t_contrib = (epsilons[id-1] / (stokes[id-1])) * tdust;
                        gas_energy[ll] += (gas_dens[ll]) * cv * (t_contrib / denom);
                        #else
                        real tdust = trad;
                        #endif
                        e[ll] = cdust * tdust * rho[ll];
                    } else {
                      #ifndef DIFFTEMP
                        e[ll] = cv * rho[ll] * trad;
                      #endif
                    }
                #endif

                // --- VELOCITY ---
                v2[ll] = v3[ll] = 0.0;
                v1[ll] = omega * r;
                if (Fluidtype == GAS) {
                    v1[ll] *= sqrt(pow(sin(Zmed(k)), -2. * FLARINGINDEX) - (beta + xi) * h * h);
                }
                v1[ll] -= OMEGAFRAME * r * sin(Zmed(k));

                soundspeed = h * sqrt(G * MSTAR / r);
                v2[ll] += soundspeed * NOISE * (drand48() - 0.5);
                v3[ll] += soundspeed * NOISE * (drand48() - 0.5);
            }
        }
    }
}

void CondInit() {
    // 1. Initialize Gas
    Fluids[0] = CreateFluid("gas", GAS);
    SelectFluid(0);
    _CondInit(0);

    // 2. Initialize Dust Fluids
    char dust_name[MAXNAMELENGTH];
    for (int id_dust = 1; id_dust < NFLUIDS; id_dust++) {
        sprintf(dust_name, "dust%d", id_dust);
        Fluids[id_dust] = CreateFluid(dust_name, DUST);
        SelectFluid(id_dust);
        _CondInit(id_dust);
    }

#ifdef SAMETEMP
    sync_dust_energy_to_gas();
#endif
}