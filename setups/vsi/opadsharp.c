#include "fargo3d.h"
#define N_NU_DATA 210
#define N_S_DATA  200
#define H_CGS 6.62607015e-27
#define KB_CGS 1.380649e-16
#define C_CGS 2.99792458e10

// --- Global variables  ---
// Add this to your globals
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

int FindNearestIndex(real target, real *array, int n_elements) {
    // 1. Determine if the array is ascending or descending
    int is_ascending = (array[n_elements - 1] > array[0]);

    // 2. Boundary checks based on sort order
    if (is_ascending) {
        if (target <= array[0]) return 0;
        if (target >= array[n_elements - 1]) return n_elements - 1;
    } else {
        if (target >= array[0]) return 0;
        if (target <= array[n_elements - 1]) return n_elements - 1;
    }
    
    // 3. Binary search
    int left = 0, right = n_elements - 1;
    while (right - left > 1) {
        int mid = left + (right - left) / 2;
        
        if (is_ascending) {
            if (target >= array[mid]) left = mid;
            else right = mid;
        } else {
            // Flipped logic for descending arrays
            if (target <= array[mid]) left = mid; 
            else right = mid;
        }
    }
    
    // 4. Calculate absolute differences to find the truly closest index
    real diff_left  = (target > array[left])  ? (target - array[left])  : (array[left]  - target);
    real diff_right = (target > array[right]) ? (target - array[right]) : (array[right] - target);
    
    if (diff_left < diff_right) {
        return left;
    } else {
        return right;
    }
}

void CreateKappaSlice(real target_size) {
    if (s_data == NULL || kappa_2d == NULL ||  nu_data == NULL) {
        printf("CRITICAL ERROR: Global opacity arrays are NULL.\n");
        exit(1);
    }

    kappa_1d_slice = (real *)malloc(N_NU_DATA * sizeof(real));
    if (kappa_1d_slice == NULL) {
        printf("CRITICAL ERROR: Memory allocation failed for kappa_1d_slice.\n");
        exit(1);
    }

    // Find the nearest size index
    int nearest_s_idx = FindNearestIndex(target_size, s_data, N_S_DATA);
    
    printf("Fluid Size %1.4e mapped to table size %1.4e (Index %d)\n", 
           target_size, s_data[nearest_s_idx], nearest_s_idx);

    // Extract that exact column for all frequencies
    for (int i = 0; i < N_NU_DATA; i++) {
        kappa_1d_slice[i] = kappa_2d[i + N_NU_DATA*nearest_s_idx];
    }
}

void SaveKappaSliceToASCII(int fluid_index, real target_size) {
    // Generate a unique filename for each fluid
    char filename[256];
    sprintf(filename, "kappa_slice_fluid_%d.dat", fluid_index);
    
    FILE *fp = fopen(filename, "w");
    if (fp == NULL) {
        printf("Warning: Could not open %s to save kappa slice.\n", filename);
        return; // Non-fatal error, just return without crashing
    }
    
    // Write a header with the specific size for reference
    fprintf(fp, "# Frequency(Hz)    Kappa(cm^2/g)   | Size = %1.6e cm\n", target_size);
    
    // Loop through the frequency grid and save the interpolated slice
    for (int i = 0; i < N_NU_DATA; i++) {
        fprintf(fp, "%1.8e\t%1.8e\n", nu_data[i], kappa_1d_slice[i]);
    }
    
    fclose(fp);
    printf("Saved 1D opacity slice to %s\n", filename);
}


real GetKappaNearest1D(real target_nu) {
    // Find the closest frequency index
    int nearest_nu_idx = FindNearestIndex(target_nu, nu_data, N_NU_DATA);
  
    // Return the exact value from our 1D slice
    return kappa_1d_slice[nearest_nu_idx];
}

void BuildMultiFluidCoolingTable(real *dust_sizes, real rhosolid) {


    //Load the raw opacity data from ASCII files
    ReadKappaASCII();

    
    
    // Setup T grid (shared by all fluids)
    real dlog_T = (log10(TMAXTAB) - log10(TMINTAB)) / (NTABLE - 1);
    for (int i = 0; i < NTABLE; i++) {
        TempTable[i] = pow(10.0, log10(TMINTAB) + i * dlog_T);
    }


    // Integration Parameters
    int N_INT = 1000; 
    real nu_min = nu_data[0];
    real nu_max = nu_data[N_NU_DATA - 1];
    real dlog_nu = (log10(nu_max) - log10(nu_min)) / (N_INT - 1);

    printf("Building cooling tables for %d dust fluids...\n", NDUST);

    // --- LOOP OVER EACH DUST FLUID ---
    for (int f = 0; f < NDUST; f++) {
        real current_size = dust_sizes[f];
        
        // 1. Slice the 2D kappa array for THIS specific fluid's size
        CreateKappaSlice(current_size); 
	    SaveKappaSliceToASCII(f,current_size);
	
        // Calculate the native logarithmic step size from your data file
        // We use fabs() to ensure it's positive regardless of sort order
        real dlog_nu = fabs(log10(nu_data[0]) - log10(nu_data[N_NU_DATA - 1])) / (N_NU_DATA - 1);

        // 2. Loop over Temperature
        for (int i = 0; i < NTABLE; i++) {
            real T = TempTable[i];
            
            real sum_Bnu = 0.0, sum_dBnu_dT = 0.0;
            real sum_Bnu_Q = 0.0, sum_dBnu_dT_Q = 0.0;
            
            // 3. Integrate DIRECTLY over the 210 native frequency points!
            for (int j = 0; j < N_NU_DATA; j++) {
                real nu = nu_data[j]; // Grab true frequency
                real d_nu = nu * dlog_nu * log(10.0);
                
                real Bnu = planck_Bnu(nu, T);
                real dBnu = planck_dBnu_dT(nu, T);
                
                // NO INTERPOLATION NEEDED! We just grab the matching opacity.
                real kappa = kappa_1d_slice[j]; 
                
                real Q_ext_d = (4.0 / 3.0) * kappa * current_size * rhosolid;
                
                sum_Bnu       += Bnu * d_nu;
                sum_dBnu_dT   += dBnu * d_nu;
                sum_Bnu_Q     += Bnu * Q_ext_d * d_nu;
                sum_dBnu_dT_Q += dBnu * Q_ext_d * d_nu;
            }
            
            real Qave = sum_Bnu_Q / sum_Bnu;
            real dQ_dT = (sum_dBnu_dT_Q / sum_Bnu) - (Qave * sum_dBnu_dT / sum_Bnu);
            
  
            // 4. Store F(T)
            int table_index = f * NTABLE + i;
            Dsharp[table_index] = 4.0 * pow(T, 3) * Qave + pow(T, 4) * dQ_dT;
        }

    }


#ifdef GPU
  DevMemcpyH2D(Dsharp_d,Dsharp,sizeof(real)*NDUST*NTABLE);
  DevMemcpyH2D(TempTable_d,TempTable,sizeof(real)*NTABLE);
#endif
}

real Interpolate_FT(real T, real *T_tab, real *FT_tab, int n_tab, int ndust) {
    // 1. Boundary clamping
    if (T <= T_tab[0]) return FT_tab[ndust*n_tab + 0];
    if (T >= T_tab[n_tab - 1]) return FT_tab[ndust*n_tab + n_tab - 1];
    
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
    int offset = ndust*n_tab;
    real t0 = T_tab[left];
    real t1 = T_tab[right];
    real f0 = FT_tab[offset+left];
    real f1 = FT_tab[offset+right];
    
    return f0 + (f1 - f0) * (T - t0) / (t1 - t0);
}
