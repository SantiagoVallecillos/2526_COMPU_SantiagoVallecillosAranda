#include <iostream>
#include <cmath>
#include <cstdlib>

#define N 25            // Tamaño de las matrices
#define TOTAL_TRIALS 10000000

static int sord[N][N];
static int sdesord[N][N];

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

    return 0;
}
