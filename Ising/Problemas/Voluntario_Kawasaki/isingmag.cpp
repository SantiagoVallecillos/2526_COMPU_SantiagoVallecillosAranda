#include <iostream>
#include <cmath>
#include <cstdlib>
#define N 20            //Tamaño de las matrices
double min(double x);

//Creamos una función mínimo entre 1 y x
double minimo(double x)
{
    double min;
    if(1 < x){min = 1;}
    else if(x <= 1){min = x;}

    return min;
}

int main()
{
    //Inicializamos para generar números aleatorios
    int semilla = 310809;
    srand(semilla);

    //Inicializamos variables
    int i, j, k, n, m, sord[N][N], sdesord[N][N], sN[64][64], randn, randm;
    double dseta, P_S, T[2], T_m[10], E_S, p, random_val, deltaE, mN16, mN32, mN64;

    FILE *ford1;
    ford1 = fopen("ising_data.dat", "w");    //Fichero para los espines de organización ordenada 1 con T baja

    FILE *ford2;
    ford2 = fopen("ising_2_data.dat", "w");    //Fichero para los espines de organización ordenada 1 con T alta

    FILE *fdesord1;
    fdesord1 = fopen("ising_desord1_data.dat", "w");    //Fichero para los espines de organización desordenada con T baja

    FILE *fdesord2;
    fdesord2 = fopen("ising_desord2_data.dat", "w");    //Fichero para los espines de organización desordenada con T alta

    FILE *mag16;
    mag16 = fopen("magn16.txt", "w"); //Fichero para la magnetización con N = 16

    FILE *mag32;
    mag32 = fopen("magn32.txt", "w"); //Fichero para la magnetización con N = 32

    FILE *mag64;
    mag64 = fopen("magn64.txt", "w"); //Fichero para la magnetización con N = 64

    T[0] = 1;
    T[1] = 3.5;
    
    //Generamos unas configuraciónes iniciales de espines
    for(n = 0; n < N; n++)
    {
        for(m = 0; m < N ; m++)
        {
            sord[n][m] = 1;

            random_val = (double)rand() / (RAND_MAX + 1.0);
            if(random_val <= 0.5){sdesord[n][m] = 1;}
            else if(random_val > 0.5){sdesord[n][m] = -1;} 
        }
        //Lo metemos todo en los ficheros 
        for(m = 0; m < N; m++)
        {
            if(m < N-1) fprintf(ford1, "%i,", sord[n][m]);
            else fprintf(ford1, "%i", sord[n][m]);
            if(m < N-1) fprintf(fdesord1, "%i,", sdesord[n][m]);
            else fprintf(fdesord1, "%i", sdesord[n][m]);
        }
        fprintf(ford1, "\n");
        fprintf(fdesord1, "\n");
    }
    fprintf(ford1, "\n\n");
    fprintf(fdesord1, "\n\n");

    //Hacemos el bucle para ver si cambian de espin
    for(i = 0; i <= 1000000; i++)
    {
        randn = rand() % N;
        randm = rand() % N;
        //printf("%i\t%i\n", randn, randm);

        // Aplicación del algoritmo de Metropolis para configuración ordenada
        // Calcular el cambio de energía deltaE
        if(randn == 0)
        {
            if(randm == 0){deltaE = 2*(sord[randn][randm])*(sord[randn+1][randm] + sord[N-1][randm] + sord[randn][randm+1] + sord[randn][0]);}
            else if(randm == N-1){deltaE = 2*(sord[randn][randm])*(sord[randn+1][randm] + sord[N-1][randm] + sord[randn][0] + sord[randn][randm-1]);}
            else{deltaE = 2*(sord[randn][randm])*(sord[randn+1][randm] + sord[N-1][randm] + sord[randn][randm+1] + sord[randn][randm-1]);}
        }
        else if(randn == N-1)
        {
            if(randm == 0){deltaE = 2*(sord[randn][randm])*(sord[0][randm] + sord[randn-1][randm] + sord[randn][randm+1] + sord[randn][N-1]);}
            else if(randm == N-1){deltaE = 2*(sord[randn][randm])*(sord[0][randm] + sord[randn-1][randm] + sord[randn][0] + sord[randn][randm-1]);}
            else{deltaE = 2*(sord[randn][randm])*(sord[0][randm] + sord[randn-1][randm] + sord[randn][randm+1] + sord[randn][randm-1]);}
        }
        else
        {
            if(randm == 0){deltaE = 2*(sord[randn][randm])*(sord[randn+1][randm] + sord[randn-1][randm] + sord[randn][randm+1] + sord[randn][N-1]);}
            else if(randm == N-1){deltaE = 2*(sord[randn][randm])*(sord[randn+1][randm] + sord[randn-1][randm] + sord[randn][0] + sord[randn][randm-1]);}
            else{deltaE = 2*(sord[randn][randm])*(sord[randn+1][randm] + sord[randn-1][randm] + sord[randn][randm+1] + sord[randn][randm-1]);}
        }
        //printf("%i\t%i\t%lf\n", randn, randm, deltaE);
        // Calcular probabilidad de aceptación p = min(1, exp(-deltaE/T))
        p = minimo(exp((-deltaE)/T[0]));

        // Generar número aleatorio y decidir si aceptar el flip
        dseta = (double)rand() / (RAND_MAX + 1.0);
        if(dseta < p){sord[randn][randm] = -sord[randn][randm];}
        else{continue;}

        if(i%(N*N) == 0)
        {
            //Lo volvemos a meter todo en el fichero
            for(n = 0; n < N; n++)
            { 
                for(m = 0; m < N; m++)
                {
                    if(m < N-1) fprintf(ford1, "%i,", sord[n][m]);
                    else fprintf(ford1, "%i", sord[n][m]);
                }
                fprintf(ford1, "\n");
            }
            fprintf(ford1, "\n");
        }

        //Para configuración desordenada
        // Aplicación del algoritmo de Metropolis para configuración desordenada
        // Calcular el cambio de energía deltaE
        if(randn == 0)
        {
            if(randm == 0){deltaE = 2*(sdesord[randn][randm])*(sdesord[randn+1][randm] + sdesord[N-1][randm] + sdesord[randn][randm+1] + sdesord[randn][N-1]);}
            else if(randm == N-1){deltaE = 2*(sdesord[randn][randm])*(sdesord[randn+1][randm] + sdesord[N-1][randm] + sdesord[randn][0] + sdesord[randn][randm-1]);}
            else{deltaE = 2*(sdesord[randn][randm])*(sdesord[randn+1][randm] + sdesord[N-1][randm] + sdesord[randn][randm+1] + sdesord[randn][randm-1]);}
        }
        else if(randn == N-1)
        {
            if(randm == 0){deltaE = 2*(sdesord[randn][randm])*(sdesord[0][randm] + sdesord[randn-1][randm] + sdesord[randn][randm+1] + sdesord[randn][N-1]);}
            else if(randm == N-1){deltaE = 2*(sdesord[randn][randm])*(sdesord[0][randm] + sdesord[randn-1][randm] + sdesord[randn][0] + sdesord[randn][randm-1]);}
            else{deltaE = 2*(sdesord[randn][randm])*(sdesord[0][randm] + sdesord[randn-1][randm] + sdesord[randn][randm+1] + sdesord[randn][randm-1]);}
        }
        else
        {
            if(randm == 0){deltaE = 2*(sdesord[randn][randm])*(sdesord[randn+1][randm] + sdesord[randn-1][randm] + sdesord[randn][randm+1] + sdesord[randn][N-1]);}
            else if(randm == N-1){deltaE = 2*(sdesord[randn][randm])*(sdesord[randn+1][randm] + sdesord[randn-1][randm] + sdesord[randn][0] + sdesord[randn][randm-1]);}
            else{deltaE = 2*(sdesord[randn][randm])*(sdesord[randn+1][randm] + sdesord[randn-1][randm] + sdesord[randn][randm+1] + sdesord[randn][randm-1]);}
        }

        // Calcular probabilidad de aceptación p = min(1, exp(-deltaE/T))
        p = minimo(exp((-deltaE)/T[0]));

        // Generar número aleatorio y decidir si aceptar el flip
        dseta = (double)rand() / (RAND_MAX + 1.0);
        if(dseta < p){sdesord[randn][randm] = -sdesord[randn][randm];}
        else{continue;}

        if(i%(N*N) == 0)
        {
            //Lo volvemos a meter todo en el fichero
            for(n = 0; n < N; n++)
            { 
                for(m = 0; m < N; m++)
                {
                    if(m < N-1) fprintf(fdesord1, "%i,", sdesord[n][m]);
                    else fprintf(fdesord1, "%i", sdesord[n][m]);
                }
                fprintf(fdesord1, "\n");
            }
            fprintf(fdesord1, "\n");
        }
    }

    //Regeneramos unas configuraciónes iniciales de espines
    for(n = 0; n < N; n++)
    {
        for(m = 0; m < N ; m++)
        {
            sord[n][m] = 1;

            random_val = (double)rand() / (RAND_MAX + 1.0);
            if(random_val <= 0.5){sdesord[n][m] = 1;}
            else if(random_val > 0.5){sdesord[n][m] = -1;} 
        }
        //Lo metemos todo en los ficheros 
        for(m = 0; m < N; m++)
        {
            if(m < N-1) fprintf(ford2, "%i,", sord[n][m]);
            else fprintf(ford2, "%i", sord[n][m]);
            if(m < N-1) fprintf(fdesord2, "%i,", sdesord[n][m]);
            else fprintf(fdesord2, "%i", sdesord[n][m]);
        }
        fprintf(ford2, "\n");
        fprintf(fdesord2, "\n");
    }
    fprintf(ford2, "\n");
    fprintf(fdesord2, "\n");

    for(i = 0; i <= 1000000; i++)
    {
        randn = rand() % N;
        randm = rand() % N;

        // Aplicación del algoritmo de Metropolis para configuración ordenada a T alta
        // Calcular el cambio de energía deltaE
        if(randn == 0)
        {
            if(randm == 0){deltaE = 2*(sord[randn][randm])*(sord[randn+1][randm] + sord[N-1][randm] + sord[randn][randm+1] + sord[randn][N-1]);}
            else if(randm == N-1){deltaE = 2*(sord[randn][randm])*(sord[randn+1][randm] + sord[N-1][randm] + sord[randn][0] + sord[randn][randm-1]);}
            else{deltaE = 2*(sord[randn][randm])*(sord[randn+1][randm] + sord[N-1][randm] + sord[randn][randm+1] + sord[randn][randm-1]);}
        }
        else if(randn == N-1)
        {
            if(randm == 0){deltaE = 2*(sord[randn][randm])*(sord[0][randm] + sord[randn-1][randm] + sord[randn][randm+1] + sord[randn][N-1]);}
            else if(randm == N-1){deltaE = 2*(sord[randn][randm])*(sord[0][randm] + sord[randn-1][randm] + sord[randn][0] + sord[randn][randm-1]);}
            else{deltaE = 2*(sord[randn][randm])*(sord[0][randm] + sord[randn-1][randm] + sord[randn][randm+1] + sord[randn][randm-1]);}
        }
        else
        {
            if(randm == 0){deltaE = 2*(sord[randn][randm])*(sord[randn+1][randm] + sord[randn-1][randm] + sord[randn][randm+1] + sord[randn][N-1]);}
            else if(randm == N-1){deltaE = 2*(sord[randn][randm])*(sord[randn+1][randm] + sord[randn-1][randm] + sord[randn][0] + sord[randn][randm-1]);}
            else{deltaE = 2*(sord[randn][randm])*(sord[randn+1][randm] + sord[randn-1][randm] + sord[randn][randm+1] + sord[randn][randm-1]);}
        }

        // Calcular probabilidad de aceptación p = min(1, exp(-deltaE/T))
        p = minimo(exp((-deltaE)/T[1]));

        // Generar número aleatorio y decidir si aceptar el flip
        dseta = (double)rand() / (RAND_MAX + 1.0);
        if(dseta < p){sord[randn][randm] = -sord[randn][randm];}
        else{continue;}

        if(i%(N*N) == 0)
        {
            //Lo volvemos a meter todo en el fichero
            for(n = 0; n < N; n++)
            { 
                for(m = 0; m < N; m++)
                {
                    if(m < N-1) fprintf(ford2, "%i,", sord[n][m]);
                    else fprintf(ford2, "%i", sord[n][m]);
                }
                fprintf(ford2, "\n");
            }
            fprintf(ford2, "\n");
        }

        //Para configuración desordenada
        // Aplicación del algoritmo de Metropolis para configuración desordenada a T alta
        // Calcular el cambio de energía deltaE
        if(randn == 0)
        {
            if(randm == 0){deltaE = 2*(sdesord[randn][randm])*(sdesord[randn+1][randm] + sdesord[N-1][randm] + sdesord[randn][randm+1] + sdesord[randn][N-1]);}
            else if(randm == N-1){deltaE = 2*(sdesord[randn][randm])*(sdesord[randn+1][randm] + sdesord[N-1][randm] + sdesord[randn][0] + sdesord[randn][randm-1]);}
            else{deltaE = 2*(sdesord[randn][randm])*(sdesord[randn+1][randm] + sdesord[N-1][randm] + sdesord[randn][randm+1] + sdesord[randn][randm-1]);}
        }
        else if(randn == N-1)
        {
            if(randm == 0){deltaE = 2*(sdesord[randn][randm])*(sdesord[0][randm] + sdesord[randn-1][randm] + sdesord[randn][randm+1] + sdesord[randn][N-1]);}
            else if(randm == N-1){deltaE = 2*(sdesord[randn][randm])*(sdesord[0][randm] + sdesord[randn-1][randm] + sdesord[randn][0] + sdesord[randn][randm-1]);}
            else{deltaE = 2*(sdesord[randn][randm])*(sdesord[0][randm] + sdesord[randn-1][randm] + sdesord[randn][randm+1] + sdesord[randn][randm-1]);}
        }
        else
        {
            if(randm == 0){deltaE = 2*(sdesord[randn][randm])*(sdesord[randn+1][randm] + sdesord[randn-1][randm] + sdesord[randn][randm+1] + sdesord[randn][N-1]);}
            else if(randm == N-1){deltaE = 2*(sdesord[randn][randm])*(sdesord[randn+1][randm] + sdesord[randn-1][randm] + sdesord[randn][0] + sdesord[randn][randm-1]);}
            else{deltaE = 2*(sdesord[randn][randm])*(sdesord[randn+1][randm] + sdesord[randn-1][randm] + sdesord[randn][randm+1] + sdesord[randn][randm-1]);}
        }
        

        // Calcular probabilidad de aceptación p = min(1, exp(-deltaE/T))
        p = minimo(exp((-deltaE)/T[1]));

        // Generar número aleatorio y decidir si aceptar el flip
        dseta = (double)rand() / (RAND_MAX + 1.0);
        if(dseta < p){sdesord[randn][randm] = -sdesord[randn][randm];}
        else{continue;}

        if(i%(N*N) == 0)
        {
            //Lo volvemos a meter todo en el fichero
            for(n = 0; n < N; n++)
            { 
                for(m = 0; m < N; m++)
                {
                    if(m < N-1) fprintf(fdesord2, "%i,", sdesord[n][m]);
                    else fprintf(fdesord2, "%i", sdesord[n][m]);
                }
                fprintf(fdesord2, "\n");
            }
            fprintf(fdesord2, "\n");
        }
    }

    //Inicializamos las temperaturas
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

    //Hacemos los cálculos de la magnetización
    for(i = 0; i < 10; i++)
    {
        mN16 = 0;
        mN32 = 0;
        mN64 = 0;
        //Generamos una configuración inicial de espines ordenada
        for(n = 0; n < 64; n++)
        {
            for(m = 0; m < 64; m++)
            {
                sN[n][m] = 1;
            }
        }
        //Hacemos el bucle para ver si cambia de espin
        for(j = 0; j <= 1000000; j++)
        {
            randn = rand() % 64;
            randm = rand() % 64;

            // Aplicación del algoritmo de Metropolis para cálculo de magnetización
            // Calcular el cambio de energía deltaE con condiciones de frontera periódicas
            int up = (randn == 0) ? 63 : randn - 1;
            int down = (randn == 63) ? 0 : randn + 1;
            int left = (randm == 0) ? 63 : randm - 1;
            int right = (randm == 63) ? 0 : randm + 1;
            deltaE = 2 * (sN[randn][randm]) * (sN[down][randm] + sN[up][randm] + sN[randn][right] + sN[randn][left]);

            // Calcular probabilidad de aceptación p = min(1, exp(-deltaE/T))
            p = minimo(exp((-deltaE)/T_m[i]));

            // Generar número aleatorio y decidir si aceptar el flip
            dseta = (double)rand() / (RAND_MAX + 1.0);
            if(dseta < p){sN[randn][randm] = -sN[randn][randm];}
            else{continue;}
        }

        //Para N = 16
        for(j = 0; j < 16; j++)
        {
            for(k = 0; k < 16; k++)
            {
                mN16 += sN[j][k];
            }
        }
        mN16 = ((fabs(mN16))/(16*16))/(10000);
        fprintf(mag16, "%lf\t%e\n", T_m[i], mN16);

        //Para N = 32
        for(j = 0; j < 32; j++)
        {
            for(k = 0; k < 32; k++)
            {
                mN32 += sN[j][k];
            }
        }
        mN32 = ((fabs(mN32))/(32*32))/(10000);
        fprintf(mag32, "%lf\t%e\n", T_m[i], mN32);

        //Para N = 64
        for(j = 0; j < 64; j++)
        {
            for(k = 0; k < 64; k++)
            {
                mN64 += sN[j][k];
            }
        }
        mN64 = ((fabs(mN64))/(64*64))/(10000);
        fprintf(mag64, "%lf\t%e\n", T_m[i], mN64);
    }

    fclose(ford1);
    fclose(fdesord1);
    fclose(ford2);
    fclose(fdesord2);
    fclose(mag16);
    fclose(mag32);
    fclose(mag64);

    return 0;
}