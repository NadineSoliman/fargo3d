#include "fargo3d.h"
#define N_NU_DATA 210
#define N_S_DATA  200
#define H_CGS 6.62607015e-27
#define KB_CGS 1.380649e-16
#define C_CGS 2.99792458e10

// --- Global variables  ---
// Add this to your globals
int N_TABLE = 256;
real *Temp_Table;
real *FT_Table_Multi; // Size will be NFLUIDS-1 * N_TABLE
real *kappa_1d_slice;
real *nu_data;
real *s_data;
real *kappa_2d;



void ReadKappaASCII() {
    FILE *fp;

    // 1. Allocate memory on the Host (CPU)
    // kappa_2d is allocated as a flat 1D array for memory efficiency
    nu_data  = (real *)malloc(N_NU_DATA * sizeof(real));
    s_data   = (real *)malloc(N_S_DATA * sizeof(real));
    kappa_2d = (real *)malloc(N_NU_DATA * N_S_DATA * sizeof(real));

    if (nu_data == NULL || s_data == NULL || kappa_2d == NULL) {
        printf("Error: Memory allocation failed for opacity arrays.\n");
        exit(1);
    }

    // 2. Read Frequency (nu) Data
    fp = fopen("setups/vsi/nu_data.dat", "r");
    if (fp == NULL) { 
        printf("Error: Cannot find nu_data.dat\n"); 
        exit(1); 
    }
    for (int i = 0; i < N_NU_DATA; i++) {
        fscanf(fp, "%lf", &nu_data[i]);
    }
    fclose(fp);

    // 3. Read Size (s) Data
    fp = fopen("setups/vsi/s_data.dat", "r");
    if (fp == NULL) { 
        printf("Error: Cannot find s_data.dat\n"); 
        exit(1); 
    }
    for (int i = 0; i < N_S_DATA; i++) {
        fscanf(fp, "%lf", &s_data[i]);
    }
    fclose(fp);

    // 4. Read the 2D Opacity (kappa) Matrix
    fp = fopen("setups/vsi/kappa_ext.dat", "r");
    if (fp == NULL) { 
        printf("Error: Cannot find kappa_2d.dat\n"); 
        exit(1); 
    }
    
    
    for (int i = 0; i < N_NU_DATA * N_S_DATA; i++) {
        fscanf(fp, "%lf", &kappa_2d[i]);
    }
    fclose(fp);
    
    printf("Successfully loaded nu, s, and kappa arrays into memory.\n");
}




// Planck Function B_nu(T)
real planck_Bnu(real nu, real T) {
    real x = (H_CGS * nu) / (KB_CGS * T);
    if (x > 500.0) return 0.0; // Prevent floating-point overflow at low T
    return (2.0 * H_CGS * pow(nu, 3) / pow(C_CGS, 2)) / (exp(x) - 1.0);
}

// Derivative of Planck Function wrt Temperature
real planck_dBnu_dT(real nu, real T) {
    real x = (H_CGS * nu) / (KB_CGS * T);
    if (x > 500.0) return 0.0;
    real exp_x = exp(x);
    return planck_Bnu(nu, T) * (x / T) * (exp_x / (exp_x - 1.0));
}


void CreateKappaSlice(real target_size) {


    // SAFETY CHECK 1: Ensure global arrays actually exist in memory
    if (s_data == NULL || kappa_2d == NULL || nu_data == NULL) {
        printf("CRITICAL ERROR: Global opacity arrays are NULL.\n");
        printf("Make sure ReadKappaASCII() is called before building the table!\n");
        exit(1);
    }

    // SAFETY CHECK 2: Ensure target_size is a valid number
    if (isnan(target_size) || target_size <= 0.0) {
        printf("CRITICAL ERROR: Invalid target_size passed: %lf\n", target_size);
        exit(1);
    }

    // Allocate the 1D slice
    kappa_1d_slice = (real *)malloc(N_NU_DATA * sizeof(real));
    if (kappa_1d_slice == NULL) {
        printf("CRITICAL ERROR: Failed to allocate memory for kappa_1d_slice.\n");
        exit(1);
    }

    // Find bounding indices for target_size in s_data
    int left = 0, right = N_S_DATA - 1;
    
    // Strict boundary clamping
    if (target_size <= s_data[0]) { 
        left = 0; 
        right = 0; 
    } else if (target_size >= s_data[N_S_DATA - 1]) { 
        left = N_S_DATA - 1; 
        right = N_S_DATA - 1; 
    } else {
        // Binary search
        while (right - left > 1) {
            int mid = left + (right - left) / 2;
            if (target_size >= s_data[mid]) left = mid;
            else right = mid;
        }
    }

    // SAFETY CHECK 3: Ensure indices are strictly within bounds
    if (left < 0 || left >= N_S_DATA || right < 0 || right >= N_S_DATA) {
        printf("CRITICAL ERROR: Binary search yielded out-of-bounds indices! Left: %d, Right: %d\n", left, right);
        exit(1);
    }

    // Log-log interpolation weights for the size axis
    real log_s_target = log10(target_size);
    real log_s0 = log10(s_data[left]);
    real log_s1 = (right == left) ? log_s0 : log10(s_data[right]);
    real weight = (right == left) ? 0.0 : (log_s_target - log_s0) / (log_s1 - log_s0);

    // Apply weight to all frequencies
    for (int i = 0; i < N_NU_DATA; i++) {
        // SAFETY CHECK 4: Ensure flat index math doesn't overflow
        int index_left = i * N_S_DATA + left;
        int index_right = i * N_S_DATA + right;
        
        if (index_left >= N_NU_DATA * N_S_DATA || index_right >= N_NU_DATA * N_S_DATA) {
            printf("CRITICAL ERROR: 2D array index out of bounds at row %d!\n", i);
            exit(1);
        }

        real k0 = kappa_2d[index_left];
        real k1 = kappa_2d[index_right];
        
        if (right == left) {
            kappa_1d_slice[i] = k0;
        } else {
            real log_k_target = log10(k0) + weight * (log10(k1) - log10(k0));
            kappa_1d_slice[i] = pow(10.0, log_k_target);
        }
    }
}

// Log-log interpolation along the frequency axis
real InterpolateKappa1D(real target_nu) {
    if (target_nu <= nu_data[0]) return kappa_1d_slice[0];
    if (target_nu >= nu_data[N_NU_DATA - 1]) return kappa_1d_slice[N_NU_DATA - 1];
    
    int left = 0, right = N_NU_DATA - 1;
    while (right - left > 1) {
        int mid = left + (right - left) / 2;
        if (target_nu >= nu_data[mid]) left = mid;
        else right = mid;
    }
    
    real log_nu0 = log10(nu_data[left]);
    real log_nu1 = log10(nu_data[right]);
    real log_k0  = log10(kappa_1d_slice[left]);
    real log_k1  = log10(kappa_1d_slice[right]);
    
    real log_target = log10(target_nu);
    real log_k_target = log_k0 + (log_target - log_nu0) * (log_k1 - log_k0) / (log_nu1 - log_nu0);
    
    return pow(10.0, log_k_target);
}

void BuildMultiFluidCoolingTable(real *dust_sizes, int num_fluids, real rhosolid) {


    //Load the raw opacity data from ASCII files
    ReadKappaASCII();

    int N_DUST = num_fluids;
    
    
    // Setup T grid (shared by all fluids)
    real T_min = 8.0, T_max = 1700.0;
    real dlog_T = (log10(T_max) - log10(T_min)) / (N_TABLE - 1);
    for (int i = 0; i < N_TABLE; i++) {
        Temp_Table[i] = pow(10.0, log10(T_min) + i * dlog_T);
    }


    // Integration Parameters
    int N_INT = 1000; 
    real nu_min = nu_data[0];
    real nu_max = nu_data[N_NU_DATA - 1];
    real dlog_nu = (log10(nu_max) - log10(nu_min)) / (N_INT - 1);

    printf("Building cooling tables for %d dust fluids...\n", N_DUST);

    // --- LOOP OVER EACH DUST FLUID ---
    for (int n = 0; n < N_DUST; n++) {
        real current_size = dust_sizes[n];
        // 1. Slice the 2D kappa array for THIS specific fluid's size
        CreateKappaSlice(current_size); 

        // 2. Loop over Temperature
        for (int i = 0; i < N_TABLE; i++) {
            real T = Temp_Table[i];
            real sum_Bnu = 0.0, sum_dBnu_dT = 0.0;
            real sum_Bnu_Q = 0.0, sum_dBnu_dT_Q = 0.0;

            // 3. Integrate over Frequency
            for (int j = 0; j < N_INT; j++) {
                real log_nu = log10(nu_min) + j * dlog_nu;
                real nu = pow(10.0, log_nu);
                real d_nu = nu * dlog_nu * log(10.0);
                
                real Bnu = planck_Bnu(nu, T);
                real dBnu = planck_dBnu_dT(nu, T);
                
                real kappa = InterpolateKappa1D(nu); 

                real Q_ext_d = (4.0 / 3.0) * kappa * current_size * RHOSOLID; 

                sum_Bnu       += Bnu * d_nu;
                sum_dBnu_dT   += dBnu * d_nu;
                sum_Bnu_Q     += Bnu * Q_ext_d * d_nu;
                sum_dBnu_dT_Q += dBnu * Q_ext_d * d_nu;
            }

            real Qave = sum_Bnu_Q / sum_Bnu;
            real dQ_dT = (sum_dBnu_dT_Q / sum_Bnu) - (Qave * sum_dBnu_dT / sum_Bnu);
            
            int table_index = n * N_TABLE + i;
            Dsharp[table_index] = 4.0 * pow(T, 3) * Qave + pow(T, 4) * dQ_dT;

        }

        // Free the 1D slice for this fluid so we can recreate it for the next fluid
        free(kappa_1d_slice); 
    }

    // Cleanup massive raw opacity arrays
    free(nu_data); free(s_data); free(kappa_2d);

#ifdef GPU
  DevMemcpyH2D(Dsharp_d,Dsharp,sizeof(real)*(NFLUIDS-1)*256);
  DevMemcpyH2D(Temp_Table_d,Temp_Table,sizeof(real)*256);
#endif
}

real Interpolate_FT(real T, real *T_tab, real *FT_tab, int n_tab) {
    // 1. Boundary clamping
    if (T <= T_tab[0]) return FT_tab[0];
    if (T >= T_tab[n_tab - 1]) return FT_tab[n_tab - 1];
    
    // 2. Binary search to find the correct temperature bin
    int left = 0;
    int right = n_tab - 1;
    
    while (right - left > 1) {
        int mid = left + (right - left) / 2;
        if (T >= T_tab[mid]) {
            left = mid;
        } else {
            right = mid;
        }
    }
    
    // 3. Linear interpolation
    real t0 = T_tab[left];
    real t1 = T_tab[right];
    real f0 = FT_tab[left];
    real f1 = FT_tab[right];
    
    return f0 + (f1 - f0) * (T - t0) / (t1 - t0);
}