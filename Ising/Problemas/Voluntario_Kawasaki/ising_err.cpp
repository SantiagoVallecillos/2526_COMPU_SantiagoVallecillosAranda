#include <iostream>
#include <cmath>
#include <cstdlib>
#include <vector> // Necesario para almacenar las mediciones

#define N 1000            // Tamaño de las matrices
#define TOTAL_TRIALS 100000000

static int sord[N][N];
static int sdesord[N][N];

// ---------------------------------------------------------
// NUEVA FUNCIÓN: Cálculo de errores estadísticos (Montecarlo)
// Basado en los apuntes proporcionados.
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

// Calcula el cambio de energía para el espín en la posición (n, m)
int delta_energy(const int spins[N][N], int n, int m)
{
    int up = spins[(n + 1) % N][m];
    int down = spins[(n - 1 + N) % N][m];
    int right = spins[n][(m + 1) % N];
    int left = spins[n][(m - 1 + N) % N];
    return 2 * spins[n][m] * (up + down + right + left);
}

// Inicializa la configuración ordenada (todos los espines +1)
void initialize_ordered(int spins[N][N])
{
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
    for (int n = 0; n < N; ++n)
    {
        for (int m = 0; m < N; ++m)
        {
            spins[n][m] = (random_double() <= 0.5) ? 1 : -1;
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

// Realiza un paso de Monte Carlo completo (un "sweep" N*N de ensayos)
void metropolis_sweep(int spins[N][N], double T)
{
    for (int trial = 0; trial < N * N; ++trial)
    {
        int n = rand() % N;
        int m = rand() % N;
        int dE = delta_energy(spins, n, m);
        double p = std::min(1.0, exp(-dE / T));
        if (random_double() < p)
        {
            spins[n][m] = -spins[n][m];
        }
    }
}

// Calcula y almacena las curvas de magnetización para diferentes temperaturas
void calculate_magnetization()
{
    #define N_MAG 64
    int sN[N_MAG][N_MAG];
    double T_m[10];
    double mN16, mN32, mN64;
    int i, j, k, n, m, randn, randm;
    double deltaE, p, dseta;

    FILE *mag16 = fopen("magn16.txt", "w");
    FILE *mag32 = fopen("magn32.txt", "w");
    FILE *mag64 = fopen("magn64.txt", "w");

    if (!mag16 || !mag32 || !mag64)
    {
        std::cerr << "Error abriendo archivos de magnetización." << std::endl;
        return;
    }

    // Inicializamos las temperaturas
    T_m[0] = 1.5;
    T_m[1] = 1.7;
    T_m[2] = 1.9;
    T_m[3] = 2.1;
    T_m[4] = 2.3;
    T_m[5] = 2.5;
    T_m[6] = 2.7;
    T_m[7] = 2.9;
    T_m[8] = 3.1;
    T_m[9] = 3.3;

    // Número de experimentos independientes por temperatura para sacar estadísticas
    int N_experimentos = 10; 

    // Hacemos los cálculos de la magnetización
    for (i = 0; i < 10; i++)
    {
        std::vector<double> medidas16(N_experimentos);
        std::vector<double> medidas32(N_experimentos);
        std::vector<double> medidas64(N_experimentos);

        // Repetimos el experimento N veces para calcular los errores
            for (int rep = 0; rep < N_experimentos; ++rep)
            {
                mN16 = 0;
                mN32 = 0;
                mN64 = 0;

            // Generamos una configuración inicial de espines ordenada
            for (n = 0; n < N_MAG; n++)
            {
                for (m = 0; m < N_MAG; m++)
                {
                    sN[n][m] = 1;
                }
            }

            // Realizamos los pasos de Monte Carlo
            for (j = 0; j <= 1000000; j++)
            {
                randn = rand() % N_MAG;
                randm = rand() % N_MAG;

                // Calcular el cambio de energía deltaE con condiciones de frontera periódicas
                int up = (randn == 0) ? N_MAG - 1 : randn - 1;
                int down = (randn == N_MAG - 1) ? 0 : randn + 1;
                int left = (randm == 0) ? N_MAG - 1 : randm - 1;
                int right = (randm == N_MAG - 1) ? 0 : randm + 1;
                deltaE = 2 * (sN[randn][randm]) * (sN[down][randm] + sN[up][randm] + sN[randn][right] + sN[randn][left]);

                // Calcular probabilidad de aceptación p = min(1, exp(-deltaE/T))
                p = std::min(1.0, exp((-deltaE) / T_m[i]));

                // Generar número aleatorio y decidir si aceptar el flip
                dseta = random_double();
                if (dseta < p)
                {
                    sN[randn][randm] = -sN[randn][randm];
                }
            }

            // Para N = 16
            for (j = 0; j < 16; j++)
            {
                for (k = 0; k < 16; k++)
                {
                    mN16 += sN[j][k];
                }
            }
            mN16 = (fabs(mN16) / (16 * 16)) / 10000;
            medidas16[rep] = mN16; // Guardamos la medición

            // Para N = 32
            for (j = 0; j < 32; j++)
            {
                for (k = 0; k < 32; k++)
                {
                    mN32 += sN[j][k];
                }
            }
            mN32 = (fabs(mN32) / (32 * 32)) / 10000;
            medidas32[rep] = mN32; // Guardamos la medición

            // Para N = 64
            for (j = 0; j < 64; j++)
            {
                for (k = 0; k < 64; k++)
                {
                    mN64 += sN[j][k];
                }
            }
            mN64 = (fabs(mN64) / (64 * 64)) / 10000;
            medidas64[rep] = mN64; // Guardamos la medición
        }

        // Llamamos a nuestra nueva función para calcular las medias y los errores
        double media16, err16, media32, err32, media64, err64;
        calcular_error_montecarlo(medidas16, media16, err16);
        calcular_error_montecarlo(medidas32, media32, err32);
        calcular_error_montecarlo(medidas64, media64, err64);

        // Imprimimos Temperatura, Media de Magnetización, y Error Estadístico
        fprintf(mag16, "%lf\t%e\t%e\n", T_m[i], media16, err16);
        fprintf(mag32, "%lf\t%e\t%e\n", T_m[i], media32, err32);
        fprintf(mag64, "%lf\t%e\t%e\n", T_m[i], media64, err64);
    }

    fclose(mag16);
    fclose(mag32);
    fclose(mag64);
}

int main()
{
    int semilla = 310809;
    srand(semilla);

    double T_low = 1.0;
    double T_high = 3.5;

    FILE *ford1 = fopen("ising_1_data.dat", "w");
    FILE *ford2 = fopen("ising_2_data.dat", "w");
    FILE *fdesord1 = fopen("ising_desord1_data.dat", "w");
    FILE *fdesord2 = fopen("ising_desord2_data.dat", "w");

    if (!ford1 || !ford2 || !fdesord1 || !fdesord2)
    {
        std::cerr << "Error abriendo archivos de salida." << std::endl;
        return 1;
    }

    int mc_steps = TOTAL_TRIALS / (N * N);
    int extra_trials = TOTAL_TRIALS % (N * N);

    // Primera fase: T baja
    initialize_ordered(sord);
    initialize_disordered(sdesord);
    write_lattice(ford1, sord);
    write_lattice(fdesord1, sdesord);

    for (int step = 0; step < mc_steps; ++step) 
    {
        metropolis_sweep(sord, T_low);
        metropolis_sweep(sdesord, T_low);
        write_lattice(ford1, sord);
        write_lattice(fdesord1, sdesord);
    }

    if (extra_trials > 0)
    {
        for (int trial = 0; trial < extra_trials; ++trial)
        {
            int n = rand() % N;
            int m = rand() % N;
            int dE = delta_energy(sord, n, m);
            double p = std::min(1.0, exp(-dE / T_low));
            if (random_double() < p)
                sord[n][m] = -sord[n][m];

            n = rand() % N;
            m = rand() % N;
            dE = delta_energy(sdesord, n, m);
            p = std::min(1.0, exp(-dE / T_low));
            if (random_double() < p)
                sdesord[n][m] = -sdesord[n][m];
        }
        write_lattice(ford1, sord);
        write_lattice(fdesord1, sdesord);
    }

    // Segunda fase: T alta
    initialize_ordered(sord);
    initialize_disordered(sdesord);
    write_lattice(ford2, sord);
    write_lattice(fdesord2, sdesord);

    for (int step = 0; step < mc_steps; ++step)
    {
        metropolis_sweep(sord, T_high);
        metropolis_sweep(sdesord, T_high);
        write_lattice(ford2, sord);
        write_lattice(fdesord2, sdesord);
    }

    if (extra_trials > 0)
    {
        for (int trial = 0; trial < extra_trials; ++trial)
        {
            int n = rand() % N;
            int m = rand() % N;
            int dE = delta_energy(sord, n, m);
            double p = std::min(1.0, exp(-dE / T_high));
            if (random_double() < p)
                sord[n][m] = -sord[n][m];

            n = rand() % N;
            m = rand() % N;
            dE = delta_energy(sdesord, n, m);
            p = std::min(1.0, exp(-dE / T_high));
            if (random_double() < p)
                sdesord[n][m] = -sdesord[n][m];
        }
        write_lattice(ford2, sord);
        write_lattice(fdesord2, sdesord);
    }

    fclose(ford1);
    fclose(fdesord1);
    fclose(ford2);
    fclose(fdesord2);

    // Calcular curvas de magnetización
    calculate_magnetization();

    return 0;
}