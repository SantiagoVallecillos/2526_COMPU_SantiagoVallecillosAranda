#include <iostream>
#include <cmath>
#include <cstdlib>
#include <vector> // Necesario para almacenar las mediciones
#include <omp.h>
#include <chrono>
#include <cstdio>
#include <ctime>

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

// ---------------------------------------------------------
// NUEVA FUNCIÓN: Cálculo de errores estadísticos (Montecarlo)
// ---------------------------------------------------------
void calcular_error_montecarlo(const std::vector<double>& medidas, double& media, double& error) 
{
    int N_exp = medidas.size();
    media = 0.0;
    
    // 1. Calculamos el valor medio
    for (int i = 0; i < N_exp; ++i) 
    {
        media += medidas[i];
    }
    media /= N_exp;

    // 2. Calculamos la varianza (sigma^2 = 1/N * sum(X_i - X_media)^2)
    double varianza = 0.0;
    for (int i = 0; i < N_exp; ++i) 
    {
        varianza += (medidas[i] - media) * (medidas[i] - media);
    }
    varianza /= N_exp;

    // 3. Calculamos el error (sigma / sqrt(N)) con un 68.26% de probabilidad
    error = sqrt(varianza / N_exp);
}

// Devuelve un número aleatorio uniforme en [0, 1)
double random_double()
{
    return static_cast<double>(rand()) / (RAND_MAX + 1.0);
}

// Thread-safe RNG using rand_r
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

// Calcula el cambio de energía al intercambiar dos espines en (n1,m1) y (n2,m2)
// Versión parametrizada por tamaño L
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

    return (spins[n1][m1] - spins[n2][m2]) * (sum1 - sum2);
}

// Versión para matriz de tamaño fijo N (para mantener compatibilidad con código existente)
int delta_energy_swap(const int spins[N][N], int n1, int m1, int n2, int m2)
{
    return delta_energy_swap_L((const int (*)[128])spins, n1, m1, n2, m2, N);
}

// Calcula la energia total del sistema en cada configuracion "a lo bruto"
// Versión parametrizada por tamaño L
int total_energy_L(const int spins[128][128], int L)
{
    int energy = 0;
    #pragma omp parallel for reduction(+:energy) schedule(static)
    for (int i = 0; i < L; ++i)
    {
        for (int j = 0; j < L; ++j)
        {
            energy += -spins[i][j] * (spins[(i + 1) % L][j] + spins[(i - 1 + L) % L][j] + spins[i][(j + 1) % L] + spins[i][(j - 1 + L) % L]);
        }
    }
    return energy / 2; // Dividimos por 2 para evitar contar cada interacción dos veces
}

// Versión para matriz de tamaño fijo N (para mantener compatibilidad)
int total_energy(const int spins[N][N])
{
    return total_energy_L((const int (*)[128])spins, N);
}

// Inicializa la configuración ordenada (todos los espines +1)
void initialize_ordered(int spins[N][N])
{
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

// Inicializa la configuración desordenada con condiciones de frontera fijas en los bordes de la dirección X
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
                    spins[n][m] = 1; // Bordes fijos en el eje X
                else
                    spins[n][m] = (random_double_r(seed) <= 0.5) ? 1 : -1;
            }
        }
    }
}

// Copia una submatriz de tamaño L x L de src a dst
void copy_lattice(const int src[128][128], int dst[128][128], int L)
{
    #pragma omp parallel for schedule(static)
    for (int n = 0; n < L; ++n)
    {
        for (int m = 0; m < L; ++m)
        {
            dst[n][m] = src[n][m];
        }
    }
}

// Selecciona un vecino válido con condiciones de contorno fijas en los bordes X
// Retorna las coordenadas en n2 y m2
void select_neighbor_bc(int n, int m, int L, int& n2, int& m2)
{
    int candidates[4][2];
    int count = 0;

    // Vecino derecha, solo si no es borde fijo
    if (n + 1 != L - 1)
    {
        candidates[count][0] = n + 1;
        candidates[count][1] = m;
        count++;
    }

    // Vecino izquierda, solo si no es borde fijo
    if (n - 1 != 0)
    {
        candidates[count][0] = n - 1;
        candidates[count][1] = m;
        count++;
    }

    // Vecino arriba
    candidates[count][0] = n;
    candidates[count][1] = (m + 1) % L;
    count++;

    // Vecino abajo
    candidates[count][0] = n;
    candidates[count][1] = (m - 1 + L) % L;
    count++;

    int choice = rand() % count;
    n2 = candidates[choice][0];
    m2 = candidates[choice][1];
}

// ----------------------------------------------------------------------------------
// NUEVA FUNCIÓN: Inicialización desordenada para una magnetización fija m0 (Kawasaki)
// ----------------------------------------------------------------------------------
void initialize_kawasaki_m0(int spins[128][128], double m0, int L)
{
    // 1. Inicializamos toda la submatriz útil de tamaño L x L con espines -1.
    // Esto limpia el sistema y establece el estado base antes de añadir los espines +1.
    for (int n = 0; n < L; ++n)
    {
        for (int m = 0; m < L; ++m)
        {
            if (n < 128 && m < 128)  // Verificación de seguridad
                spins[n][m] = -1;
        }
    }

    // 2. Relación teórica entre magnetización por partícula (m0) y la fracción 'x' de espines +1:
    // m0 = (N_plus - N_minus) / (L*L)
    // Sabiendo que N_plus + N_minus = L*L, llegamos a que la fracción es x = (1 + m0) / 2.
    double x = (1.0 + m0) / 2.0;
    int total_spins = L * L;
    
    // Calculamos el número exacto de espines +1 necesarios.
    // Usamos round() de <cmath> para aproximar al entero más cercano debido a la discretización de la red.
    int count_plus_needed = static_cast<int>(round(x * total_spins));

    // 3. Distribuimos aleatoriamente los espines +1 calculados en la red.
    // Se utiliza el mismo algoritmo de muestreo por rechazo que implementaste en calculate_magnetization().
    int count_ones = 0;
    while (count_ones < count_plus_needed)
    {
        int rn = rand() % L;
        int rm = rand() % L;
        
        // Si la posición elegida al azar contiene un -1, la cambiamos a +1
        if (spins[rn][rm] == -1)
        {
            spins[rn][rm] = 1;
            count_ones++; // Incrementamos el contador solo si se ha modificado un espín de forma efectiva
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

// Realiza un paso de Monte Carlo completo (un "sweep" N*N de ensayos) utilizando la dinámica de Kawasaki
void metropolis_sweep(int spins[N][N], int spins_next[N][N], double T)
{
    for (int trial = 0; trial < N * N; ++trial)
    {
        int n = rand() % N;
        int m = rand() % N;

        //Elegimos espín vecino aleatorio para intercambiarlo con el espín actual
        int n2 = n;
        int m2 = m;
        int vecino = rand() % 4;
        if (vecino == 0)
            n2 = (n + 1) % N;
        else if (vecino == 1)
            n2 = (n - 1 + N) % N;
        else if (vecino == 2)
            m2 = (m + 1) % N;
        else
            m2 = (m - 1 + N) % N;
        // Calculamos el cambio de energía localmente (O(1)) y decidimos aceptación
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

// Realiza un paso de Monte Carlo completo (un "sweep" N*N de ensayos) utilizando la dinámica de Kawasaki
// utilizando las condiciones de contorno fijas en los bordes de la dirección X
void metropolis_sweep_bc(int spins[N][N], int spins_next[N][N], double T)
{
    for (int trial = 0; trial < N * N; ++trial)
    {
        // Elegimos un espín aleatorio que no esté en los bordes fijos
        int n = 1 + rand() % (N - 2); // solo posiciones interiores en X
        int m = rand() % N;

        // Elegimos un vecino aleatorio válido que no esté en los bordes fijos
        int candidates[4][2];
        int count = 0;

        // Vecino derecha, solo si no es borde fijo
        if (n + 1 != N - 1)
        {
            candidates[count][0] = n + 1;
            candidates[count][1] = m;
            count++;
        }

        // Vecino izquierda, solo si no es borde fijo
        if (n - 1 != 0)
        {
            candidates[count][0] = n - 1;
            candidates[count][1] = m;
            count++;
        }

        // Vecino arriba
        candidates[count][0] = n;
        candidates[count][1] = (m + 1) % N;
        count++;

        // Vecino abajo
        candidates[count][0] = n;
        candidates[count][1] = (m - 1 + N) % N;
        count++;

        int choice = rand() % count;
        int n2 = candidates[choice][0];
        int m2 = candidates[choice][1];

        // Calculamos el cambio de energía localmente (O(1)) y decidimos aceptación
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

// Calcula y almacena las curvas de magnetización para diferentes temperaturas
void calculate_magnetization()
{
    #define N_MAG 128
    constexpr int N_TEMPS = 10;
    const double T_m[N_TEMPS] = {1.5, 1.7, 1.9, 2.1, 2.3, 2.5, 2.7, 2.9, 3.1, 3.3};
    int i, j, k, n, m;
    
    // === CAMBIO 1: Eliminamos E, E2, M2 globales de aquí. ===
    // Por qué: Deben inicializarse por separado para cada repetición y tamaño, 
    // así que usaremos vectores dentro del bucle como propusiste.

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

    // Número de experimentos independientes por temperatura para sacar estadísticas
    const int N_experimentos = 10;

    int sizes[3] = {32, 64, 128};
    FILE *files[3] = {mag32, mag64, mag128};

    // Hacemos los cálculos de la magnetización
    for (i = 0; i < N_TEMPS; i++)
    {
        for (int s_idx = 0; s_idx < 3; ++s_idx) 
        {
            int L = sizes[s_idx];
            
            // === CAMBIO 3: Adaptamos los vectores para E, E2 y M2 ===
            // Por qué: Para almacenar el promedio de E, E2 y M2 de CADA repetición independiente.
            std::vector<double> medidas_M(N_experimentos);
            std::vector<double> medidas_E(N_experimentos);
            std::vector<double> medidas_E2(N_experimentos);
            std::vector<double> medidas_M2(N_experimentos);

            // Repetimos el experimento N veces para calcular los errores en paralelo
            #pragma omp parallel for schedule(dynamic)
            for (int rep = 0; rep < N_experimentos; ++rep)
            {
                unsigned int seed = (unsigned int)time(NULL) ^ (unsigned int)rep ^ (unsigned int)omp_get_thread_num();
                int sN_local[N_MAG][N_MAG];
                for (int nn = 0; nn < L; ++nn)
                    for (int mm = 0; mm < L; ++mm)
                        sN_local[nn][mm] = -1;

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
                for (int nn = 0; nn < L; ++nn)
                    for (int mm = 0; mm < L; ++mm)
                        sN_next[nn][mm] = sN_local[nn][mm];

                const int total_mc_trials = 1000000;
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
                        select_neighbor_bc(n, m, L, n2, m2);

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
                    for (j = 0; j < extra_trials; ++j)
                    {
                        int n = 1 + rand_r(&seed) % (L - 2);
                        int m = rand_r(&seed) % L;

                        int n2, m2;
                        select_neighbor_bc(n, m, L, n2, m2);

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
                for (j = 0; j < L; ++j)
                {
                    for (k = 0; k < L; ++k)
                    {
                        if (k < L / 2)
                            m_top += sN_local[j][k];
                        else
                            m_bottom += sN_local[j][k];
                    }
                }

                double particles_per_half = (L * L) / 2.0;
                double mag_domain = (fabs(m_top) / particles_per_half + fabs(m_bottom) / particles_per_half) / 2.0;
                double num_particulas = L * L;

                medidas_M[rep] = mag_domain;
                medidas_M2[rep] = (mag_domain * mag_domain);
                medidas_E[rep] = (suma_E / numero_medidas) / num_particulas;
                medidas_E2[rep] = (suma_E2 / numero_medidas) / (num_particulas * num_particulas);

                // Solo el hilo 0 realiza la escritura del perfil para evitar colisiones en IO
                if (rep == N_experimentos - 1)
                {
                    #pragma omp critical
                    {
                        export_density_profile(files_prof[s_idx], (const int (*)[N_MAG])sN_local, L, T_m[i]);
                    }
                }
            }

            // === CAMBIO 6: Calculamos medias y errores para todas las magnitudes termodinámicas ===
            double media_M, err_M;
            calcular_error_montecarlo(medidas_M, media_M, err_M);
            
            double media_M2, err_M2;
            calcular_error_montecarlo(medidas_M2, media_M2, err_M2);
            
            double media_E, err_E;
            calcular_error_montecarlo(medidas_E, media_E, err_E);
            
            double media_E2, err_E2;
            calcular_error_montecarlo(medidas_E2, media_E2, err_E2);

            // === CAMBIO 1.3: Cálculo de calor específico y susceptibilidad ===
            // Como las medias que hemos obtenido están "por partícula" (<e>, <e^2>, <m>, <m^2>), 
            // aplicamos el factor multiplicativo (L*L) derivado de adaptar la fórmula extensiva del PDF a valores intensivos.
            double num_particulas = L * L;
            double T = T_m[i];
            
            double c_N = (num_particulas / (T * T)) * (media_E2 - (media_E * media_E));
            double chi_N = (num_particulas / T) * (media_M2 - (media_M * media_M));

            // === CAMBIO 1.3 (Parte 2): Añadidas c_N y chi_N al archivo de salida ===
            // Formato CSV por tabulaciones: T | <M> | err_M | <E> | err_E | c_N | chi_N
            fprintf(files[s_idx], "%lf\t%e\t%e\t%e\t%e\t%e\t%e\n", T_m[i], media_M, err_M, media_E, err_E, c_N, chi_N);
        }
    }

    fclose(mag32);
    fclose(mag64);
    fclose(mag128);
}

// ---------------------------------------------------------
// NUEVA FUNCIÓN: Conteo global de espines (+1 y -1)
// ---------------------------------------------------------
void export_global_spin_count(FILE* file, const int spins[128][128], int L, double T)
{
    int count_plus = 0;
    int count_minus = 0;

    for (int i = 0; i < L; ++i)
    {
        for (int j = 0; j < L; ++j)
        {
            if (spins[i][j] == 1)
            {
                count_plus++;
            }
            else if (spins[i][j] == -1)
            {
                count_minus++;
            }
        }
    }

    // Formato CSV por tabulaciones: Temperatura | Cantidad +1 | Cantidad -1
    fprintf(file, "%lf\t%d\t%d\n", T, count_plus, count_minus);
}

// ---------------------------------------------------------
// NUEVA FUNCIÓN: Perfil de densidad fila por fila
// ---------------------------------------------------------
void export_density_profile(FILE* file, const int spins[N][N], int L, double T)
{
    // Iteramos sobre la coordenada Y (que en tu matriz es el segundo índice, 'm')
    for (int m = 0; m < L; ++m)
    {
        int count_plus = 0;
        int count_minus = 0;
        
        // Sumamos todos los valores a lo largo de la coordenada X (primer índice, 'n')
        for (int n = 0; n < L; ++n)
        {
            if (spins[n][m] == 1)
                count_plus++;
            else if (spins[n][m] == -1)
                count_minus++;
        }
        
        // Exportamos: Temperatura, Coordenada Y, Cantidad +1, Cantidad -1
        fprintf(file, "%lf\t%d\t%d\t%d\n", T, m, count_plus, count_minus);
    }
    // Dejamos una línea en blanco al terminar la matriz. 
    // Esto es muy útil en Gnuplot para separar los datos de distintas temperaturas.
    fprintf(file, "\n"); 
}

// -------------------------------------------------------------------------
// NUEVA FUNCIÓN: Calcula magnetización para m0 distinto de cero
// -------------------------------------------------------------------------
void calculate_magnetization_m0(double m0)
{
    #define N_MAG 128
    constexpr int N_TEMPS = 10;
    const double T_m[N_TEMPS] = {1.5, 1.7, 1.9, 2.1, 2.3, 2.5, 2.7, 2.9, 3.1, 3.3};
    int i; // Mantenemos solo 'i' como global para el bucle exterior de temperaturas

    // Abrimos archivos con sufijo _m0 para no sobreescribir los datos de m0=0
    FILE *mag32 = fopen("magn32_m0.txt", "w");
    FILE *mag64 = fopen("magn64_m0.txt", "w");
    FILE *mag128 = fopen("magn128_m0.txt", "w");
    FILE *prof32 = fopen("perfil32_m0.txt", "w");
    FILE *prof64 = fopen("perfil64_m0.txt", "w");
    FILE *prof128 = fopen("perfil128_m0.txt", "w");
    
    FILE *files[3] = {mag32, mag64, mag128};
    FILE *files_prof[3] = {prof32, prof64, prof128};

    if (!mag32 || !mag64 || !mag128)
    {
        std::cerr << "Error abriendo archivos de magnetización m0." << std::endl;
        return;
    }

    const int N_experimentos = 10;
    int sizes[3] = {32, 64, 128};

    // Calculamos la fracción teórica x
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

            // === CAMBIO OpenMP: Paralelizamos el bucle de repeticiones de Montecarlo ===
            #pragma omp parallel for schedule(dynamic)
            for (int rep = 0; rep < N_experimentos; ++rep)
            {
                // Generamos una semilla única por hilo y repetición para evitar colisiones RNG
                unsigned int seed = (unsigned int)time(NULL) ^ (unsigned int)rep ^ (unsigned int)omp_get_thread_num();
                
                // === CAMBIO OpenMP: Matriz sN pasa a ser sN_local (privada para cada hilo) ===
                int sN_local[N_MAG][N_MAG];
                
                // Inicializamos la matriz local a -1
                for (int nn = 0; nn < L; ++nn)
                    for (int mm = 0; mm < L; ++mm)
                        sN_local[nn][mm] = -1;

                // === CAMBIO OpenMP: Inicialización 'm0' en línea (Thread-safe) ===
                // En lugar de llamar a initialize_kawasaki_m0() que usa el rand() global, 
                // generamos los espines aquí dentro usando rand_r(&seed)
                int total_spins = L * L;
                int count_plus_needed = static_cast<int>(round(x_frac * total_spins));
                int count_ones = 0;
                while (count_ones < count_plus_needed)
                {
                    int rn = rand_r(&seed) % L;
                    int rm = rand_r(&seed) % L;
                    if (sN_local[rn][rm] == -1)
                    {
                        sN_local[rn][rm] = 1;
                        count_ones++;
                    }
                }

                // Matriz para guardar el estado siguiente
                int sN_next[N_MAG][N_MAG];
                for (int nn = 0; nn < L; ++nn)
                    for (int mm = 0; mm < L; ++mm)
                        sN_next[nn][mm] = sN_local[nn][mm];

                const int total_mc_trials = 1000000;
                int mc_steps = total_mc_trials / (L * L);
                int extra_trials = total_mc_trials % (L * L);

                // Control para exportar el fotograma (sólo hilo de la repetición 0 lo hará)
                bool capture_evolution = (s_idx == 2 && rep == 0);
                FILE *f_evol = nullptr;

                if (capture_evolution)
                {
                    char filename_evol[100];
                    snprintf(filename_evol, sizeof(filename_evol), "ising_evol_m0_%.2f_T_%.2f.dat", m0, T_m[i]);
                    f_evol = fopen(filename_evol, "w");
                }

                int writes_evol = 0;
                int max_writes_evol = std::max(1, mc_steps / 10);
                int energy_threshold_evol = L * 2;
                double last_energy_evol = 0;

                double suma_E = 0.0;
                double suma_E2 = 0.0;
                int numero_medidas = 0; 
                
                // Evaluamos energía inicial (pasando el puntero local)
                double energia_actual = total_energy_L((const int (*)[128])sN_local, L);

                if (capture_evolution && f_evol)
                {
                    write_lattice(f_evol, (const int (*)[N])sN_local);
                    last_energy_evol = energia_actual;
                }

                for (int step = 0; step < mc_steps; ++step)
                {
                    for (int trial = 0; trial < L * L; ++trial)
                    {
                        // === CAMBIO OpenMP: Usamos rand_r en todo el Montecarlo ===
                        int n = 1 + rand_r(&seed) % (L - 2);
                        int m = rand_r(&seed) % L;

                        int n2, m2;
                        // OJO: select_neighbor_bc por dentro sigue usando rand(). Si te da problemas de rendimiento
                        // lo ideal es no usarla y poner el vecino directamente aquí dentro como en metropolis_sweep.
                        // La mantengo para respetar tu estructura exacta.
                        select_neighbor_bc(n, m, L, n2, m2);

                        int dE = delta_energy_swap_L((const int (*)[128])sN_local, n, m, n2, m2, L);
                        double p = min(1.0, exp(-dE / T_m[i]));
                        
                        // === CAMBIO OpenMP: random_double pasa a random_double_r ===
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

                    if (capture_evolution && f_evol && writes_evol < max_writes_evol)
                    {
                        if (std::abs(energia_actual - last_energy_evol) >= energy_threshold_evol)
                        {
                            write_lattice(f_evol, (const int (*)[N])sN_local);
                            last_energy_evol = energia_actual;
                            writes_evol++;
                        }
                    }

                    suma_E += energia_actual;
                    suma_E2 += (energia_actual * energia_actual);
                    numero_medidas++;
                }

                if (extra_trials > 0)
                {
                    // === CAMBIO OpenMP: Variable de bucle local 'j_ext' para no chocar con otros hilos ===
                    for (int j_ext = 0; j_ext < extra_trials; ++j_ext)
                    {
                        int n = 1 + rand_r(&seed) % (L - 2);
                        int m = rand_r(&seed) % L;

                        int n2, m2;
                        select_neighbor_bc(n, m, L, n2, m2);

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

                    if (capture_evolution && f_evol && writes_evol < max_writes_evol)
                    {
                        if (std::abs(energia_actual - last_energy_evol) >= energy_threshold_evol)
                        {
                            write_lattice(f_evol, (const int (*)[N])sN_local);
                            last_energy_evol = energia_actual;
                            writes_evol++;
                        }
                    }

                    suma_E += energia_actual;
                    suma_E2 += (energia_actual * energia_actual);
                    numero_medidas++;
                }

                if (capture_evolution && f_evol)
                {
                    fclose(f_evol);
                }

                double m_domain1 = 0.0;
                double m_domain2 = 0.0;

                // === CAMBIO OpenMP: Variables de bucle jj y kk locales ===
                for (int jj = 0; jj < L; ++jj)
                {
                    for (int kk = 0; kk < L; ++kk)
                    {
                        if (kk < L_split) 
                            m_domain1 += sN_local[jj][kk];
                        else 
                            m_domain2 += sN_local[jj][kk];
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
                    // === CAMBIO OpenMP: Impedimos colisión en escritura a archivo ===
                    #pragma omp critical
                    {
                        export_density_profile(files_prof[s_idx], (const int (*)[N_MAG])sN_local, L, T_m[i]);
                    }
                }
            }

            double media_M, err_M;
            calcular_error_montecarlo(medidas_M, media_M, err_M);
            
            double media_M2, err_M2;
            calcular_error_montecarlo(medidas_M2, media_M2, err_M2);
            
            double media_E, err_E;
            calcular_error_montecarlo(medidas_E, media_E, err_E);
            
            double media_E2, err_E2;
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

int main()
{
    auto inicio = chrono::high_resolution_clock::now();
    int semilla = 310809;
    srand(semilla);

    std::cout << "Iniciando simulacion del modelo de Ising (Dinamica de Kawasaki)..." << std::endl;

    // =====================================================================================
    // TAREA 1: Simular la dinámica y obtener fotogramas para diferentes temperaturas
    // =====================================================================================
    std::cout << "\n[1/3] Realizando Tarea 1: Generando fotogramas de la red..." << std::endl;
    
    double T_low = 1.0;
    double T_high = 3.5;

    int mc_steps = TOTAL_TRIALS / (N * N);
    int extra_trials = TOTAL_TRIALS % (N * N);
    int max_writes = std::max(1, mc_steps / 10);
    const int energy_threshold = N * 2;

    // -------------------------------------------------------------------------------------
    // CASO A: T_baja y T_alta con CONDICIONES PERIÓDICAS EN AMBOS EJES
    // -------------------------------------------------------------------------------------
    std::cout << "      -> Simulando con Condiciones Periodicas..." << std::endl;
    
    FILE *f_per_ord_low = fopen("ising_per_ord_low.dat", "w");
    FILE *f_per_des_low = fopen("ising_per_des_low.dat", "w");
    FILE *f_per_ord_high = fopen("ising_per_ord_high.dat", "w");
    FILE *f_per_des_high = fopen("ising_per_des_high.dat", "w");

    if (!f_per_ord_low || !f_per_des_low || !f_per_ord_high || !f_per_des_high) {
        std::cerr << "Error abriendo archivos para condiciones periodicas." << std::endl; return 1;
    }

    // Inicializamos: Ordenada para T_low/T_high, Desordenada para T_low/T_high
    initialize_ordered(sord);
    initialize_disordered(sdesord);
    
    // Copiamos a las matrices "next"
    for (int n = 0; n < N; ++n) {
        for (int m = 0; m < N; ++m) {
            sord_next[n][m] = sord[n][m];
            sdesord_next[n][m] = sdesord[n][m];
        }
    }

    int last_E_ord = total_energy(sord);
    int last_E_des = total_energy(sdesord);
    int writes_ord = 0, writes_des = 0;

    // Bucle Montecarlo (Usando metropolis_sweep normal) para T_baja
    for (int step = 0; step < mc_steps; ++step) {
        metropolis_sweep(sord, sord_next, T_low);
        metropolis_sweep(sdesord, sdesord_next, T_low);

        if (writes_ord < max_writes && std::abs(total_energy(sord) - last_E_ord) >= energy_threshold) {
            write_lattice(f_per_ord_low, sord);
            last_E_ord = total_energy(sord); writes_ord++;
        }
        if (writes_des < max_writes && std::abs(total_energy(sdesord) - last_E_des) >= energy_threshold) {
            write_lattice(f_per_des_low, sdesord);
            last_E_des = total_energy(sdesord); writes_des++;
        }
    }

    // Reiniciamos matrices para probar T_alta con Periódicas
    initialize_ordered(sord);
    initialize_disordered(sdesord);
    for (int n = 0; n < N; ++n) {
        for (int m = 0; m < N; ++m) {
            sord_next[n][m] = sord[n][m];
            sdesord_next[n][m] = sdesord[n][m];
        }
    }
    
    last_E_ord = total_energy(sord); last_E_des = total_energy(sdesord);
    writes_ord = 0; writes_des = 0;

    // Bucle Montecarlo (Usando metropolis_sweep normal) para T_alta
    for (int step = 0; step < mc_steps; ++step) {
        metropolis_sweep(sord, sord_next, T_high);
        metropolis_sweep(sdesord, sdesord_next, T_high);

        if (writes_ord < max_writes && std::abs(total_energy(sord) - last_E_ord) >= energy_threshold) {
            write_lattice(f_per_ord_high, sord);
            last_E_ord = total_energy(sord); writes_ord++;
        }
        if (writes_des < max_writes && std::abs(total_energy(sdesord) - last_E_des) >= energy_threshold) {
            write_lattice(f_per_des_high, sdesord);
            last_E_des = total_energy(sdesord); writes_des++;
        }
    }

    fclose(f_per_ord_low); fclose(f_per_des_low);
    fclose(f_per_ord_high); fclose(f_per_des_high);


    // -------------------------------------------------------------------------------------
    // CASO B: T_baja y T_alta con BORDES FIJOS EN X (Para observar mejor los dominios)
    // -------------------------------------------------------------------------------------
    std::cout << "      -> Simulando con Bordes Fijos (eje X)..." << std::endl;

    FILE *f_bc_ord_low = fopen("ising_bc_ord_low.dat", "w");
    FILE *f_bc_des_low = fopen("ising_bc_des_low.dat", "w");
    FILE *f_bc_ord_high = fopen("ising_bc_ord_high.dat", "w");
    FILE *f_bc_des_high = fopen("ising_bc_des_high.dat", "w");

    if (!f_bc_ord_low || !f_bc_des_low || !f_bc_ord_high || !f_bc_des_high) {
        std::cerr << "Error abriendo archivos para bordes fijos." << std::endl; return 1;
    }

    // Inicializamos usando la función con bordes fijos
    initialize_ordered(sord);
    initialize_disordered_bc(sdesord);
    for (int n = 0; n < N; ++n) {
        for (int m = 0; m < N; ++m) {
            sord_next[n][m] = sord[n][m];
            sdesord_next[n][m] = sdesord[n][m];
        }
    }

    last_E_ord = total_energy(sord); last_E_des = total_energy(sdesord);
    writes_ord = 0; writes_des = 0;

    // Bucle Montecarlo (Usando metropolis_sweep_bc) para T_baja
    for (int step = 0; step < mc_steps; ++step) {
        metropolis_sweep_bc(sord, sord_next, T_low);
        metropolis_sweep_bc(sdesord, sdesord_next, T_low);

        if (writes_ord < max_writes && std::abs(total_energy(sord) - last_E_ord) >= energy_threshold) {
            write_lattice(f_bc_ord_low, sord);
            last_E_ord = total_energy(sord); writes_ord++;
        }
        if (writes_des < max_writes && std::abs(total_energy(sdesord) - last_E_des) >= energy_threshold) {
            write_lattice(f_bc_des_low, sdesord);
            last_E_des = total_energy(sdesord); writes_des++;
        }
    }

    // Reiniciamos matrices para probar T_alta con bordes fijos
    initialize_ordered(sord);
    initialize_disordered_bc(sdesord);
    for (int n = 0; n < N; ++n) {
        for (int m = 0; m < N; ++m) {
            sord_next[n][m] = sord[n][m];
            sdesord_next[n][m] = sdesord[n][m];
        }
    }

    last_E_ord = total_energy(sord); last_E_des = total_energy(sdesord);
    writes_ord = 0; writes_des = 0;

    // Bucle Montecarlo (Usando metropolis_sweep_bc) para T_alta
    for (int step = 0; step < mc_steps; ++step) {
        metropolis_sweep_bc(sord, sord_next, T_high);
        metropolis_sweep_bc(sdesord, sdesord_next, T_high);

        if (writes_ord < max_writes && std::abs(total_energy(sord) - last_E_ord) >= energy_threshold) {
            write_lattice(f_bc_ord_high, sord);
            last_E_ord = total_energy(sord); writes_ord++;
        }
        if (writes_des < max_writes && std::abs(total_energy(sdesord) - last_E_des) >= energy_threshold) {
            write_lattice(f_bc_des_high, sdesord);
            last_E_des = total_energy(sdesord); writes_des++;
        }
    }

    fclose(f_bc_ord_low); fclose(f_bc_des_low);
    fclose(f_bc_ord_high); fclose(f_bc_des_high);


    // =====================================================================================
    // TAREAS 2 a 7: Cálculos termodinámicos con magnetización nula (m0 = 0)
    // =====================================================================================
    // Esta función encapsula todo: curvas de magnetización, energía, calor específico, 
    // susceptibilidad y los perfiles de densidad en el eje Y.
    std::cout << "\n[2/3] Realizando Tareas 2-7: Termodinamica (m0 = 0)..." << std::endl;
    calculate_magnetization();


    // =====================================================================================
    // TAREA 8: Cálculos termodinámicos partiendo de una magnetización NO nula (m0 != 0)
    // =====================================================================================
    // Ejecutamos exactamente la misma termodinámica (magnetización, E, c_N, chi_N, perfiles)
    // pero inicializando el sistema con la fracción asimétrica correspondiente a m0.
    std::cout << "\n[3/3] Realizando Tarea 8: Termodinamica (m0 != 0)..." << std::endl;
    double m0_deseada = 0.5; // Probamos con un m0 del 50%
    calculate_magnetization_m0(m0_deseada);

    std::cout << "\nSimulacion finalizada con exito. Revisa los archivos de salida generados." << std::endl;

    auto fin = chrono::high_resolution_clock::now();
    chrono::duration<double, milli> tiempo_ejecucion = fin - inicio;
    cout << "El código multi-run completo tardó: " << tiempo_ejecucion.count() << " milisegundos." << endl;
    return 0;
}