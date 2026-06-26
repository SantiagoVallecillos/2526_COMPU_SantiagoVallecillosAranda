#include <iostream>
#include <cmath>
#include <chrono>
#include <cstdlib>
#include <vector>
#include <omp.h>   // === CAMBIO OPENMP: Librería para paralelización ===
#include <ctime>   // === CAMBIO OPENMP: Para generar semillas de tiempo ===

#define N 128            // Tamaño de las matrices
#define TOTAL_TRIALS 100000000

using namespace std;

static int sord[N][N];
static int sdesord[N][N];
static int sord_next[N][N];
static int sdesord_next[N][N];

// === FORWARD DECLARATIONS ===
void export_density_profile(FILE* file, const int spins[N][N], int L, double T);
int delta_energy_swap_L(const int spins[128][128], int n1, int m1, int n2, int m2, int L);
int total_energy_L(const int spins[128][128], int L);
void select_neighbor_bc(int n, int m, int L, int& n2, int& m2);
void copy_lattice(const int src[128][128], int dst[128][128], int L);

// === CAMBIO OPENMP: Declaración de versiones thread-safe ===
void select_neighbor_bc_r(int n, int m, int L, int& n2, int& m2, unsigned int& seed);
void initialize_kawasaki_m0_r(int spins[128][128], double m0, int L, unsigned int& seed);


// ---------------------------------------------------------
// CÁLCULO DE ERRORES ESTADÍSTICOS (Montecarlo)
// ---------------------------------------------------------
void calcular_error_montecarlo(const std::vector<double>& medidas, double& media, double& error) 
{
    int N_exp = medidas.size();
    media = 0.0;
    
    for (int i = 0; i < N_exp; ++i) {
        media += medidas[i];
    }
    media /= N_exp;

    double varianza = 0.0;
    for (int i = 0; i < N_exp; ++i) {
        varianza += (medidas[i] - media) * (medidas[i] - media);
    }
    varianza /= N_exp;

    error = sqrt(varianza / N_exp);
}

// Devuelve un número aleatorio uniforme en [0, 1) - VERSIÓN GLOBAL
double random_double()
{
    return static_cast<double>(rand()) / (RAND_MAX + 1.0);
}

// === CAMBIO OPENMP: VERSIÓN THREAD-SAFE ===
inline double random_double_r(unsigned int &seed)
{
    return static_cast<double>(rand_r(&seed)) / (RAND_MAX + 1.0);
}

// Calcula el cambio de energía para el espín en la posición (n, m)
int delta_energy(const int spins[N][N], int n, int m)
{
    int up = spins[(n + 1) % N][m];
    int down = spins[(n - 1 + N) % N][m];
    int right = spins[n][(m + 1) % N];
    int left = spins[n][(m - 1 + N) % N];
    return 2 * spins[n][m] * (up + down + right + left);
}

// Calcula el cambio de energía al intercambiar dos espines
int delta_energy_swap_L(const int spins[128][128], int n1, int m1, int n2, int m2, int L)
{
    if (n1 == n2 && m1 == m2) return 0;

    int up1 = spins[(n1 + 1) % L][m1];
    int down1 = spins[(n1 - 1 + L) % L][m1];
    int right1 = spins[n1][(m1 + 1) % L];
    int left1 = spins[n1][(m1 - 1 + L) % L];
    int sum1 = up1 + down1 + right1 + left1;

    int up2 = spins[(n2 + 1) % L][m2];
    int down2 = spins[(n2 - 1 + L) % L][m2];
    int right2 = spins[n2][(m2 + 1) % L];
    int left2 = spins[n2][(m2 - 1 + L) % L];
    int sum2 = up2 + down2 + right2 + left2;

    return (spins[n1][m1] - spins[n2][m2]) * ((sum1 - spins[n2][m2]) - (sum2 - spins[n1][m1]));
}

int delta_energy_swap(const int spins[N][N], int n1, int m1, int n2, int m2)
{
    return delta_energy_swap_L((const int (*)[128])spins, n1, m1, n2, m2, N);
}

// Calcula la energia total del sistema 
int total_energy_L(const int spins[128][128], int L)
{
    int energy = 0;
    // === CAMBIO OPENMP: Paralelización con reducción para acelerar la suma ===
    #pragma omp parallel for reduction(+:energy) schedule(static)
    for (int i = 0; i < L; ++i)
    {
        for (int j = 0; j < L; ++j)
        {
            energy += -spins[i][j] * (spins[(i + 1) % L][j] + spins[(i - 1 + L) % L][j] + spins[i][(j + 1) % L] + spins[i][(j - 1 + L) % L]);
        }
    }
    return energy / 2;
}

int total_energy(const int spins[N][N])
{
    return total_energy_L((const int (*)[128])spins, N);
}

// Inicializa la configuración ordenada
void initialize_ordered(int spins[N][N])
{
    // === CAMBIO OPENMP: Aceleración de matrices ===
    #pragma omp parallel for collapse(2) schedule(static)
    for (int n = 0; n < N; ++n)
    {
        for (int m = 0; m < N; ++m)
        {
            spins[n][m] = 1;
        }
    }
}

// Inicializa la configuración desordenada
void initialize_disordered(int spins[N][N])
{
    #pragma omp parallel
    {
        unsigned int seed = (unsigned int)time(NULL) ^ (unsigned int)omp_get_thread_num();
        #pragma omp for collapse(2) schedule(static)
        for (int n = 0; n < N; ++n)
        {
            for (int m = 0; m < N; ++m)
            {
                spins[n][m] = (random_double_r(seed) <= 0.5) ? 1 : -1;
            }
        }
    }
}

// Inicializa la configuración desordenada con condiciones de frontera fijas
void initialize_disordered_bc(int spins[N][N])
{
    #pragma omp parallel
    {
        unsigned int seed = (unsigned int)time(NULL) ^ (unsigned int)omp_get_thread_num();
        #pragma omp for collapse(2) schedule(static)
        for (int n = 0; n < N; ++n)
        {
            for (int m = 0; m < N; ++m)
            {
                if (n == 0 || n == N - 1)
                    spins[n][m] = 1;
                else
                    spins[n][m] = (random_double_r(seed) <= 0.5) ? 1 : -1;
            }
        }
    }
}

// Copia una submatriz de tamaño L x L de src a dst
void copy_lattice(const int src[128][128], int dst[128][128], int L)
{
    #pragma omp parallel for collapse(2) schedule(static)
    for (int n = 0; n < L; ++n)
    {
        for (int m = 0; m < L; ++m)
        {
            dst[n][m] = src[n][m];
        }
    }
}

// Selecciona un vecino válido (VERSIÓN GLOBAL)
void select_neighbor_bc(int n, int m, int L, int& n2, int& m2)
{
    int candidates[4][2];
    int count = 0;

    if (n + 1 != L - 1) { candidates[count][0] = n + 1; candidates[count][1] = m; count++; }
    if (n - 1 != 0) { candidates[count][0] = n - 1; candidates[count][1] = m; count++; }
    candidates[count][0] = n; candidates[count][1] = (m + 1) % L; count++;
    candidates[count][0] = n; candidates[count][1] = (m - 1 + L) % L; count++;

    int choice = rand() % count;
    n2 = candidates[choice][0];
    m2 = candidates[choice][1];
}

// === CAMBIO OPENMP: Selecciona un vecino válido (VERSIÓN THREAD-SAFE) ===
void select_neighbor_bc_r(int n, int m, int L, int& n2, int& m2, unsigned int& seed)
{
    int candidates[4][2];
    int count = 0;

    if (n + 1 != L - 1) { candidates[count][0] = n + 1; candidates[count][1] = m; count++; }
    if (n - 1 != 0) { candidates[count][0] = n - 1; candidates[count][1] = m; count++; }
    candidates[count][0] = n; candidates[count][1] = (m + 1) % L; count++;
    candidates[count][0] = n; candidates[count][1] = (m - 1 + L) % L; count++;

    int choice = rand_r(&seed) % count;
    n2 = candidates[choice][0];
    m2 = candidates[choice][1];
}

// Inicialización m0 (VERSIÓN GLOBAL)
void initialize_kawasaki_m0(int spins[128][128], double m0, int L)
{
    for (int n = 0; n < L; ++n) {
        for (int m = 0; m < L; ++m) {
            if (n < 128 && m < 128) spins[n][m] = -1;
        }
    }
    double x = (1.0 + m0) / 2.0;
    int total_spins = L * L;
    int count_plus_needed = static_cast<int>(round(x * total_spins));

    int count_ones = 0;
    while (count_ones < count_plus_needed) {
        int rn = rand() % L;
        int rm = rand() % L;
        if (spins[rn][rm] == -1) {
            spins[rn][rm] = 1;
            count_ones++;
        }
    }
}

// === CAMBIO OPENMP: Inicialización m0 (VERSIÓN THREAD-SAFE) ===
void initialize_kawasaki_m0_r(int spins[128][128], double m0, int L, unsigned int& seed)
{
    for (int n = 0; n < L; ++n) {
        for (int m = 0; m < L; ++m) {
            if (n < 128 && m < 128) spins[n][m] = -1;
        }
    }
    double x = (1.0 + m0) / 2.0;
    int total_spins = L * L;
    int count_plus_needed = static_cast<int>(round(x * total_spins));

    int count_ones = 0;
    while (count_ones < count_plus_needed) {
        int rn = rand_r(&seed) % L;
        int rm = rand_r(&seed) % L;
        if (spins[rn][rm] == -1) {
            spins[rn][rm] = 1;
            count_ones++;
        }
    }
}

// Escribe la matriz de espines en el fichero con formato CSV
void write_lattice(FILE *file, const int spins[N][N])
{
    for (int n = 0; n < N; ++n)
    {
        for (int m = 0; m < N; ++m)
        {
            if (m < N - 1)
                fprintf(file, "%i,", spins[n][m]);
            else
                fprintf(file, "%i", spins[n][m]);
        }
        fprintf(file, "\n");
    }
    fprintf(file, "\n");
}

void metropolis_sweep(int spins[N][N], int spins_next[N][N], double T)
{
    for (int trial = 0; trial < N * N; ++trial)
    {
        int n = rand() % N;
        int m = rand() % N;
        int n2 = n, m2 = m;
        int vecino = rand() % 4;
        if (vecino == 0) n2 = (n + 1) % N;
        else if (vecino == 1) n2 = (n - 1 + N) % N;
        else if (vecino == 2) m2 = (m + 1) % N;
        else m2 = (m - 1 + N) % N;

        int dE = delta_energy_swap(spins, n, m, n2, m2);
        double p = min(1.0, exp(-dE / T));
        if (random_double() < p)
        {
            int temp = spins[n][m];
            spins[n][m] = spins[n2][m2];
            spins[n2][m2] = temp;
            spins_next[n][m] = spins[n][m];
            spins_next[n2][m2] = spins[n2][m2];
        }
        else
        {
            spins_next[n][m] = spins[n][m];
            spins_next[n2][m2] = spins[n2][m2];
        }
    }
}

void metropolis_sweep_bc(int spins[N][N], int spins_next[N][N], double T)
{
    for (int trial = 0; trial < N * N; ++trial)
    {
        int n = 1 + rand() % (N - 2); 
        int m = rand() % N;
        int n2, m2;
        select_neighbor_bc(n, m, N, n2, m2);

        int dE = delta_energy_swap(spins, n, m, n2, m2);
        double p = min(1.0, exp(-dE / T));
        if (random_double() < p)
        {
            int temp = spins[n][m];
            spins[n][m] = spins[n2][m2];
            spins[n2][m2] = temp;
            spins_next[n][m] = spins[n][m];
            spins_next[n2][m2] = spins[n2][m2];
        }
        else
        {
            spins_next[n][m] = spins[n][m];
            spins_next[n2][m2] = spins[n2][m2];
        }
    }
}

// ====================================================================================
// TERMODINÁMICA CON OPENMP PARA m0 = 0
// ====================================================================================
void calculate_magnetization()
{
    #define N_MAG 128
    constexpr int N_TEMPS = 10;
    const double T_m[N_TEMPS] = {1.5, 1.7, 1.9, 2.1, 2.3, 2.5, 2.7, 2.9, 3.1, 3.3};
    int i; // Única variable global para el bucle for externo

    FILE *mag32 = fopen("magn32.txt", "w");
    FILE *mag64 = fopen("magn64.txt", "w");
    FILE *mag128 = fopen("magn128.txt", "w");
    FILE *prof32 = fopen("perfil32.txt", "w");
    FILE *prof64 = fopen("perfil64.txt", "w");
    FILE *prof128 = fopen("perfil128.txt", "w");
    FILE *files_prof[3] = {prof32, prof64, prof128};

    if (!mag32 || !mag64 || !mag128)
    {
        std::cerr << "Error abriendo archivos de magnetización." << std::endl;
        return;
    }

    const int N_experimentos = 10;
    int sizes[3] = {32, 64, 128};
    FILE *files[3] = {mag32, mag64, mag128};

    for (i = 0; i < N_TEMPS; i++)
    {
        for (int s_idx = 0; s_idx < 3; ++s_idx) 
        {
            int L = sizes[s_idx];
            
            std::vector<double> medidas_M(N_experimentos);
            std::vector<double> medidas_E(N_experimentos);
            std::vector<double> medidas_E2(N_experimentos);
            std::vector<double> medidas_M2(N_experimentos);

            // === CAMBIO OPENMP: Paralelizamos la repetición de experimentos ===
            #pragma omp parallel for schedule(dynamic)
            for (int rep = 0; rep < N_experimentos; ++rep)
            {
                // Semilla única por hilo
                unsigned int seed = (unsigned int)time(NULL) ^ (unsigned int)rep ^ (unsigned int)omp_get_thread_num();
                
                // === CAMBIO OPENMP: Matriz local exclusiva para este hilo ===
                int sN_local[N_MAG][N_MAG];
                
                for (int nn = 0; nn < L; nn++) {
                    for (int mm = 0; mm < L; mm++) {
                        sN_local[nn][mm] = -1; 
                    }
                }
                
                int half_spins = (L * L) / 2;
                int count_ones = 0;
                while (count_ones < half_spins) 
                {
                    int rn = rand_r(&seed) % L;
                    int rm = rand_r(&seed) % L;
                    if (sN_local[rn][rm] == -1) 
                    {
                        sN_local[rn][rm] = 1;
                        count_ones++;
                    }
                }

                int sN_next[N_MAG][N_MAG];
                for (int nn = 0; nn < L; ++nn) {
                    for (int mm = 0; mm < L; ++mm) {
                        sN_next[nn][mm] = sN_local[nn][mm];
                    }
                }

                const long long total_mc_trials = TOTAL_TRIALS;
                int mc_steps = total_mc_trials / (L * L);
                int extra_trials = total_mc_trials % (L * L);

                double suma_E = 0.0;
                double suma_E2 = 0.0;
                int numero_medidas = 0; 
                
                double energia_actual = total_energy_L((const int (*)[128])sN_local, L);

                for (int step = 0; step < mc_steps; ++step)
                {
                    for (int trial = 0; trial < L * L; ++trial)
                    {
                        int n = 1 + rand_r(&seed) % (L - 2);
                        int m = rand_r(&seed) % L;
                        int n2, m2;
                        
                        select_neighbor_bc_r(n, m, L, n2, m2, seed);

                        int dE = delta_energy_swap_L((const int (*)[128])sN_local, n, m, n2, m2, L);
                        double p = min(1.0, exp(-dE / T_m[i]));
                        if (random_double_r(seed) < p)
                        {
                            int temp = sN_local[n][m];
                            sN_local[n][m] = sN_local[n2][m2];
                            sN_local[n2][m2] = temp;
                            sN_next[n][m] = sN_local[n][m];
                            sN_next[n2][m2] = sN_local[n2][m2];
                            energia_actual += dE;
                        }
                    }
                    suma_E += energia_actual;
                    suma_E2 += (energia_actual * energia_actual);
                    numero_medidas++;
                }

                if (extra_trials > 0)
                {
                    for (int j_ext = 0; j_ext < extra_trials; ++j_ext)
                    {
                        int n = 1 + rand_r(&seed) % (L - 2);
                        int m = rand_r(&seed) % L;
                        int n2, m2;
                        
                        select_neighbor_bc_r(n, m, L, n2, m2, seed);

                        int dE = delta_energy_swap_L((const int (*)[128])sN_local, n, m, n2, m2, L);
                        double p = min(1.0, exp(-dE / T_m[i]));
                        if (random_double_r(seed) < p)
                        {
                            int temp = sN_local[n][m];
                            sN_local[n][m] = sN_local[n2][m2];
                            sN_local[n2][m2] = temp;
                            sN_next[n][m] = sN_local[n][m];
                            sN_next[n2][m2] = sN_local[n2][m2];
                            energia_actual += dE;
                        }
                    }
                    suma_E += energia_actual;
                    suma_E2 += (energia_actual * energia_actual);
                    numero_medidas++;
                }

                double m_top = 0.0;
                double m_bottom = 0.0;

                for (int jj = 0; jj < L; ++jj)
                {
                    for (int kk = 0; kk < L; ++kk)
                    {
                        if (kk < L / 2) 
                            m_top += sN_local[jj][kk];
                        else 
                            m_bottom += sN_local[jj][kk];
                    }
                }
                
                double particles_per_half = (L * L) / 2.0;
                double mag_domain = (fabs(m_top) / particles_per_half + fabs(m_bottom) / particles_per_half) / 2.0;
                double num_particulas = L * L;
                
                // ESCRITURA SEGURA: Dado que "rep" es único para este hilo, es thread-safe
                medidas_M[rep] = mag_domain; 
                medidas_M2[rep] = (mag_domain * mag_domain);
                medidas_E[rep] = (suma_E / numero_medidas) / num_particulas;
                medidas_E2[rep] = (suma_E2 / numero_medidas) / (num_particulas * num_particulas);
                
                if (rep == N_experimentos - 1)
                {
                    // === CAMBIO OPENMP: Protegemos la escritura al disco ===
                    #pragma omp critical
                    {
                        export_density_profile(files_prof[s_idx], (const int (*)[N_MAG])sN_local, L, T_m[i]);
                    }
                }
            }

            double media_M, err_M, media_M2, err_M2, media_E, err_E, media_E2, err_E2;
            calcular_error_montecarlo(medidas_M, media_M, err_M);
            calcular_error_montecarlo(medidas_M2, media_M2, err_M2);
            calcular_error_montecarlo(medidas_E, media_E, err_E);
            calcular_error_montecarlo(medidas_E2, media_E2, err_E2);

            double num_particulas = L * L;
            double T = T_m[i];
            
            double c_N = (num_particulas / (T * T)) * (media_E2 - (media_E * media_E));
            double chi_N = (num_particulas / T) * (media_M2 - (media_M * media_M));

            fprintf(files[s_idx], "%lf\t%e\t%e\t%e\t%e\t%e\t%e\n", T_m[i], media_M, err_M, media_E, err_E, c_N, chi_N);
        }
    }

    fclose(mag32); fclose(mag64); fclose(mag128);
    fclose(prof32); fclose(prof64); fclose(prof128);
}


void export_global_spin_count(FILE* file, const int spins[128][128], int L, double T)
{
    int count_plus = 0;
    int count_minus = 0;
    for (int i = 0; i < L; ++i) {
        for (int j = 0; j < L; ++j) {
            if (spins[i][j] == 1) count_plus++;
            else if (spins[i][j] == -1) count_minus++;
        }
    }
    fprintf(file, "%lf\t%d\t%d\n", T, count_plus, count_minus);
}

void export_density_profile(FILE* file, const int spins[N][N], int L, double T)
{
    for (int m = 0; m < L; ++m)
    {
        int count_plus = 0;
        int count_minus = 0;
        for (int n = 0; n < L; ++n)
        {
            if (spins[n][m] == 1) count_plus++;
            else if (spins[n][m] == -1) count_minus++;
        }
        fprintf(file, "%lf\t%d\t%d\t%d\n", T, m, count_plus, count_minus);
    }
    fprintf(file, "\n"); 
}

// ====================================================================================
// TERMODINÁMICA CON OPENMP PARA m0 != 0
// ====================================================================================
void calculate_magnetization_m0(double m0)
{
    #define N_MAG 128
    constexpr int N_TEMPS = 10;
    const double T_m[N_TEMPS] = {1.5, 1.7, 1.9, 2.1, 2.3, 2.5, 2.7, 2.9, 3.1, 3.3};
    int i; 

    char name_mag32[100], name_mag64[100], name_mag128[100];
    char name_prof32[100], name_prof64[100], name_prof128[100];

    snprintf(name_mag32, sizeof(name_mag32), "magn32_m0_%.2f.txt", m0);
    snprintf(name_mag64, sizeof(name_mag64), "magn64_m0_%.2f.txt", m0);
    snprintf(name_mag128, sizeof(name_mag128), "magn128_m0_%.2f.txt", m0);

    snprintf(name_prof32, sizeof(name_prof32), "perfil32_m0_%.2f.txt", m0);
    snprintf(name_prof64, sizeof(name_prof64), "perfil64_m0_%.2f.txt", m0);
    snprintf(name_prof128, sizeof(name_prof128), "perfil128_m0_%.2f.txt", m0);

    FILE *mag32 = fopen(name_mag32, "w");
    FILE *mag64 = fopen(name_mag64, "w");
    FILE *mag128 = fopen(name_mag128, "w");
    FILE *prof32 = fopen(name_prof32, "w");
    FILE *prof64 = fopen(name_prof64, "w");
    FILE *prof128 = fopen(name_prof128, "w");
    
    FILE *files[3] = {mag32, mag64, mag128};
    FILE *files_prof[3] = {prof32, prof64, prof128};

    if (!mag32 || !mag64 || !mag128)
    {
        std::cerr << "Error abriendo archivos de magnetización m0." << std::endl;
        return;
    }

    const int N_experimentos = 10;
    int sizes[3] = {32, 64, 128};
    double x_frac = (1.0 + m0) / 2.0;

    for (i = 0; i < N_TEMPS; i++)
    {
        for (int s_idx = 0; s_idx < 3; ++s_idx) 
        {
            int L = sizes[s_idx];
            int L_split = static_cast<int>(round(x_frac * L)); 
            
            std::vector<double> medidas_M(N_experimentos);
            std::vector<double> medidas_E(N_experimentos);
            std::vector<double> medidas_E2(N_experimentos);
            std::vector<double> medidas_M2(N_experimentos);

            // === CAMBIO OPENMP: Bucle parelelizado y seguro ===
            #pragma omp parallel for schedule(dynamic)
            for (int rep = 0; rep < N_experimentos; ++rep)
            {
                unsigned int seed = (unsigned int)time(NULL) ^ (unsigned int)rep ^ (unsigned int)omp_get_thread_num();
                
                int sN_local[N_MAG][N_MAG];
                initialize_kawasaki_m0_r((int (*)[N_MAG])sN_local, m0, L, seed);

                int sN_next[N_MAG][N_MAG];
                for (int nn = 0; nn < L; ++nn) {
                    for (int mm = 0; mm < L; ++mm) {
                        sN_next[nn][mm] = sN_local[nn][mm];
                    }
                }

                const long long total_mc_trials = TOTAL_TRIALS;
                int mc_steps = total_mc_trials / (L * L);
                int extra_trials = total_mc_trials % (L * L);

                bool capture_evolution = (s_idx == 2 && rep == 0);
                FILE *f_evol = nullptr;

                if (capture_evolution)
                {
                    char filename_evol[100];
                    snprintf(filename_evol, sizeof(filename_evol), "ising_evol_m0_%.2f_T_%.2f.dat", m0, T_m[i]);
                    f_evol = fopen(filename_evol, "w");
                }

                double suma_E = 0.0;
                double suma_E2 = 0.0;
                int numero_medidas = 0; 
                
                double energia_actual = total_energy_L((const int (*)[128])sN_local, L);

                for (int step = 0; step < mc_steps; ++step)
                {
                    for (int trial = 0; trial < L * L; ++trial)
                    {
                        int n = 1 + rand_r(&seed) % (L - 2);
                        int m = rand_r(&seed) % L;
                        int n2, m2;
                        
                        select_neighbor_bc_r(n, m, L, n2, m2, seed);

                        int dE = delta_energy_swap_L((const int (*)[128])sN_local, n, m, n2, m2, L);
                        double p = min(1.0, exp(-dE / T_m[i]));
                        if (random_double_r(seed) < p)
                        {
                            int temp = sN_local[n][m];
                            sN_local[n][m] = sN_local[n2][m2];
                            sN_local[n2][m2] = temp;
                            sN_next[n][m] = sN_local[n][m];
                            sN_next[n2][m2] = sN_local[n2][m2];
                            
                            energia_actual += dE;
                        }
                    }

                    if (capture_evolution && f_evol)
                    {
                        if (step % std::max(1, mc_steps / 10) == 0)
                        {
                            // === CAMBIO OPENMP: Protegemos el fotograma ===
                            #pragma omp critical
                            {
                                write_lattice(f_evol, (const int (*)[N])sN_local);
                            }
                        }
                    }

                    suma_E += energia_actual;
                    suma_E2 += (energia_actual * energia_actual);
                    numero_medidas++;
                }

                if (extra_trials > 0)
                {
                    for (int j_ext = 0; j_ext < extra_trials; ++j_ext)
                    {
                        int n = 1 + rand_r(&seed) % (L - 2);
                        int m = rand_r(&seed) % L;
                        int n2, m2;
                        
                        select_neighbor_bc_r(n, m, L, n2, m2, seed);

                        int dE = delta_energy_swap_L((const int (*)[128])sN_local, n, m, n2, m2, L);
                        double p = min(1.0, exp(-dE / T_m[i]));
                        if (random_double_r(seed) < p)
                        {
                            int temp = sN_local[n][m];
                            sN_local[n][m] = sN_local[n2][m2];
                            sN_local[n2][m2] = temp;
                            sN_next[n][m] = sN_local[n][m];
                            sN_next[n2][m2] = sN_local[n2][m2];
                            
                            energia_actual += dE;
                        }
                    }
                    suma_E += energia_actual;
                    suma_E2 += (energia_actual * energia_actual);
                    numero_medidas++;
                }
                
                if (capture_evolution && f_evol)
                {
                    #pragma omp critical
                    {
                        write_lattice(f_evol, (const int (*)[N])sN_local);
                        fclose(f_evol);
                    }
                }

                double m_domain1 = 0.0;
                double m_domain2 = 0.0;

                for (int jj = 0; jj < L; ++jj)
                {
                    for (int kk = 0; kk < L; ++kk)
                    {
                        if (kk < L_split) m_domain1 += sN_local[jj][kk];
                        else m_domain2 += sN_local[jj][kk];
                    }
                }
                
                double particles_1 = L_split * L;
                double particles_2 = (L - L_split) * L;
                
                double mag_1 = (particles_1 > 0) ? (fabs(m_domain1) / particles_1) : 0.0;
                double mag_2 = (particles_2 > 0) ? (fabs(m_domain2) / particles_2) : 0.0;
                
                double mag_domain;
                if (particles_1 == 0) mag_domain = mag_2;
                else if (particles_2 == 0) mag_domain = mag_1;
                else mag_domain = (mag_1 + mag_2) / 2.0;
                
                double num_particulas = L * L;
                
                medidas_M[rep] = mag_domain; 
                medidas_M2[rep] = (mag_domain * mag_domain);
                medidas_E[rep] = (suma_E / numero_medidas) / num_particulas;
                medidas_E2[rep] = (suma_E2 / numero_medidas) / (num_particulas * num_particulas);
                
                if (rep == N_experimentos - 1)
                {
                    #pragma omp critical
                    {
                        export_density_profile(files_prof[s_idx], (const int (*)[N_MAG])sN_local, L, T_m[i]);
                    }
                }
            }

            double media_M, err_M, media_M2, err_M2, media_E, err_E, media_E2, err_E2;
            calcular_error_montecarlo(medidas_M, media_M, err_M);
            calcular_error_montecarlo(medidas_M2, media_M2, err_M2);
            calcular_error_montecarlo(medidas_E, media_E, err_E);
            calcular_error_montecarlo(medidas_E2, media_E2, err_E2);

            double num_particulas = L * L;
            double T = T_m[i];
            
            double c_N = (num_particulas / (T * T)) * (media_E2 - (media_E * media_E));
            double chi_N = (num_particulas / T) * (media_M2 - (media_M * media_M));

            fprintf(files[s_idx], "%lf\t%e\t%e\t%e\t%e\t%e\t%e\n", T_m[i], media_M, err_M, media_E, err_E, c_N, chi_N);
        }
    }

    fclose(mag32); fclose(mag64); fclose(mag128);
    fclose(prof32); fclose(prof64); fclose(prof128);
}


// ====================================================================================
// FUNCIÓN MAIN
// ====================================================================================
int main()
{
    auto inicio = chrono::high_resolution_clock::now();
    int semilla = 310809;
    srand(semilla);

    std::cout << "Iniciando simulacion del modelo de Ising (Dinamica de Kawasaki)..." << std::endl;

    std::cout << "\n[1/3] Realizando Tarea 1: Generando fotogramas de la red..." << std::endl;
    
    double T_low = 1.0;
    double T_high = 3.5;

    int mc_steps = TOTAL_TRIALS / (N * N);
    int extra_trials = TOTAL_TRIALS % (N * N);

    std::cout << "      -> Simulando con Condiciones Periodicas..." << std::endl;
    
    FILE *f_per_ord_low = fopen("ising_per_ord_low.dat", "w");
    FILE *f_per_des_low = fopen("ising_per_des_low.dat", "w");
    FILE *f_per_ord_high = fopen("ising_per_ord_high.dat", "w");
    FILE *f_per_des_high = fopen("ising_per_des_high.dat", "w");

    if (!f_per_ord_low || !f_per_des_low || !f_per_ord_high || !f_per_des_high) {
        std::cerr << "Error abriendo archivos para condiciones periodicas." << std::endl; return 1;
    }

    initialize_ordered(sord);
    initialize_disordered(sdesord);
    
    copy_lattice((const int (*)[128])sord, (int (*)[128])sord_next, N);
    copy_lattice((const int (*)[128])sdesord, (int (*)[128])sdesord_next, N);

    for (int step = 0; step < mc_steps; ++step) {
        metropolis_sweep(sord, sord_next, T_low);
        metropolis_sweep(sdesord, sdesord_next, T_low);

        if (step % std::max(1, mc_steps / 10) == 0) {
            write_lattice(f_per_ord_low, sord);
            write_lattice(f_per_des_low, sdesord);
        }
    }

    if (extra_trials > 0)
    {
        for (int trial = 0; trial < extra_trials; ++trial)
        {
            int n = rand() % N; int m = rand() % N;
            int n2 = n; int m2 = m;
            int vecino = rand() % 4;
            if (vecino == 0) n2 = (n + 1) % N; else if (vecino == 1) n2 = (n - 1 + N) % N;
            else if (vecino == 2) m2 = (m + 1) % N; else m2 = (m - 1 + N) % N;
            int dE = delta_energy_swap(sord, n, m, n2, m2);
            double p = std::min(1.0, exp(-dE / T_low));
            if (random_double() < p) { int temp = sord[n][m]; sord[n][m] = sord[n2][m2]; sord[n2][m2] = temp; }

            n = rand() % N; m = rand() % N; n2 = n; m2 = m;
            vecino = rand() % 4;
            if (vecino == 0) n2 = (n + 1) % N; else if (vecino == 1) n2 = (n - 1 + N) % N;
            else if (vecino == 2) m2 = (m + 1) % N; else m2 = (m - 1 + N) % N;
            dE = delta_energy_swap(sdesord, n, m, n2, m2);
            p = std::min(1.0, exp(-dE / T_low));
            if (random_double() < p) { int temp = sdesord[n][m]; sdesord[n][m] = sdesord[n2][m2]; sdesord[n2][m2] = temp; }
        }
    }

    write_lattice(f_per_ord_low, sord);
    write_lattice(f_per_des_low, sdesord);

    initialize_ordered(sord);
    initialize_disordered(sdesord);
    copy_lattice((const int (*)[128])sord, (int (*)[128])sord_next, N);
    copy_lattice((const int (*)[128])sdesord, (int (*)[128])sdesord_next, N);
    
    for (int step = 0; step < mc_steps; ++step) {
        metropolis_sweep(sord, sord_next, T_high);
        metropolis_sweep(sdesord, sdesord_next, T_high);
        if (step % std::max(1, mc_steps / 10) == 0) {
            write_lattice(f_per_ord_high, sord);
            write_lattice(f_per_des_high, sdesord);
        }
    }

    if (extra_trials > 0)
    {
        for (int trial = 0; trial < extra_trials; ++trial)
        {
            int n = rand() % N; int m = rand() % N; int n2 = n; int m2 = m;
            int vecino = rand() % 4;
            if (vecino == 0) n2 = (n + 1) % N; else if (vecino == 1) n2 = (n - 1 + N) % N;
            else if (vecino == 2) m2 = (m + 1) % N; else m2 = (m - 1 + N) % N;
            int dE = delta_energy_swap(sord, n, m, n2, m2);
            double p = std::min(1.0, exp(-dE / T_high));
            if (random_double() < p) { int temp = sord[n][m]; sord[n][m] = sord[n2][m2]; sord[n2][m2] = temp; }

            n = rand() % N; m = rand() % N; n2 = n; m2 = m;
            vecino = rand() % 4;
            if (vecino == 0) n2 = (n + 1) % N; else if (vecino == 1) n2 = (n - 1 + N) % N;
            else if (vecino == 2) m2 = (m + 1) % N; else m2 = (m - 1 + N) % N;
            dE = delta_energy_swap(sdesord, n, m, n2, m2);
            p = std::min(1.0, exp(-dE / T_high));
            if (random_double() < p) { int temp = sdesord[n][m]; sdesord[n][m] = sdesord[n2][m2]; sdesord[n2][m2] = temp; }
        }
    }

    write_lattice(f_per_ord_high, sord);
    write_lattice(f_per_des_high, sdesord);

    fclose(f_per_ord_low); fclose(f_per_des_low);
    fclose(f_per_ord_high); fclose(f_per_des_high);


    std::cout << "      -> Simulando con Bordes Fijos (eje X)..." << std::endl;

    FILE *f_bc_ord_low = fopen("ising_bc_ord_low.dat", "w");
    FILE *f_bc_des_low = fopen("ising_bc_des_low.dat", "w");
    FILE *f_bc_ord_high = fopen("ising_bc_ord_high.dat", "w");
    FILE *f_bc_des_high = fopen("ising_bc_des_high.dat", "w");

    if (!f_bc_ord_low || !f_bc_des_low || !f_bc_ord_high || !f_bc_des_high) {
        std::cerr << "Error abriendo archivos para bordes fijos." << std::endl; return 1;
    }

    initialize_ordered(sord);
    initialize_disordered_bc(sdesord);
    copy_lattice((const int (*)[128])sord, (int (*)[128])sord_next, N);
    copy_lattice((const int (*)[128])sdesord, (int (*)[128])sdesord_next, N);

    for (int step = 0; step < mc_steps; ++step) {
        metropolis_sweep_bc(sord, sord_next, T_low);
        metropolis_sweep_bc(sdesord, sdesord_next, T_low);
        if (step % std::max(1, mc_steps / 10) == 0) {
            write_lattice(f_bc_ord_low, sord);
            write_lattice(f_bc_des_low, sdesord);
        }
    }

    if (extra_trials > 0)
    {
        for (int trial = 0; trial < extra_trials; ++trial)
        {
            int n = 1 + rand() % (N - 2); int m = rand() % N; int n2, m2;
            select_neighbor_bc(n, m, N, n2, m2);
            int dE = delta_energy_swap(sord, n, m, n2, m2);
            double p = std::min(1.0, exp(-dE / T_low));
            if (random_double() < p) { int temp = sord[n][m]; sord[n][m] = sord[n2][m2]; sord[n2][m2] = temp; }

            n = 1 + rand() % (N - 2); m = rand() % N;
            select_neighbor_bc(n, m, N, n2, m2);
            dE = delta_energy_swap(sdesord, n, m, n2, m2);
            p = std::min(1.0, exp(-dE / T_low));
            if (random_double() < p) { int temp = sdesord[n][m]; sdesord[n][m] = sdesord[n2][m2]; sdesord[n2][m2] = temp; }
        }
    }

    write_lattice(f_bc_ord_low, sord);
    write_lattice(f_bc_des_low, sdesord);

    initialize_ordered(sord);
    initialize_disordered_bc(sdesord);
    copy_lattice((const int (*)[128])sord, (int (*)[128])sord_next, N);
    copy_lattice((const int (*)[128])sdesord, (int (*)[128])sdesord_next, N);

    for (int step = 0; step < mc_steps; ++step) {
        metropolis_sweep_bc(sord, sord_next, T_high);
        metropolis_sweep_bc(sdesord, sdesord_next, T_high);
        if (step % std::max(1, mc_steps / 10) == 0) {
            write_lattice(f_bc_ord_high, sord);
            write_lattice(f_bc_des_high, sdesord);
        }
    }

    if (extra_trials > 0)
    {
        for (int trial = 0; trial < extra_trials; ++trial)
        {
            int n = 1 + rand() % (N - 2); int m = rand() % N; int n2, m2;
            select_neighbor_bc(n, m, N, n2, m2);
            int dE = delta_energy_swap(sord, n, m, n2, m2);
            double p = std::min(1.0, exp(-dE / T_high));
            if (random_double() < p) { int temp = sord[n][m]; sord[n][m] = sord[n2][m2]; sord[n2][m2] = temp; }

            n = 1 + rand() % (N - 2); m = rand() % N;
            select_neighbor_bc(n, m, N, n2, m2);
            dE = delta_energy_swap(sdesord, n, m, n2, m2);
            p = std::min(1.0, exp(-dE / T_high));
            if (random_double() < p) { int temp = sdesord[n][m]; sdesord[n][m] = sdesord[n2][m2]; sdesord[n2][m2] = temp; }
        }
    }

    write_lattice(f_bc_ord_high, sord);
    write_lattice(f_bc_des_high, sdesord);

    fclose(f_bc_ord_low); fclose(f_bc_des_low);
    fclose(f_bc_ord_high); fclose(f_bc_des_high);


    std::cout << "\n[2/3] Realizando Tareas 2-7: Termodinamica (m0 = 0)..." << std::endl;
    calculate_magnetization();


    std::cout << "\n[3/3] Realizando Tarea 8: Termodinamica (m0 != 0)..." << std::endl;
    double m0_deseada_1 = 0.5; 
    double m0_deseada_2 = 0.8; 
    calculate_magnetization_m0(m0_deseada_1);
    calculate_magnetization_m0(m0_deseada_2);
    std::cout << "\nSimulacion finalizada con exito. Revisa los archivos de salida generados." << std::endl;
    
    auto fin = chrono::high_resolution_clock::now();
    chrono::duration<double, milli> tiempo_ejecucion = fin - inicio;
    cout << "El código multi-run completo tardó: " << tiempo_ejecucion.count() << " milisegundos." << endl;

    return 0;
}