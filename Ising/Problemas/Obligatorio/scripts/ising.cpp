#include <iostream>
#include <cmath>
#include <cstdlib>
#define N 1000            //Tamaño de las matrices
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
    int i, n, m, sord[N][N], sdesord[N][N], randn, randm;
    double dseta, T[2], p, random_val, deltaE;

    FILE *ford1;
    ford1 = fopen("ising_1_data.dat", "w");    //Fichero para los espines de organización ordenada 1 con T baja

    FILE *ford2;
    ford2 = fopen("ising_2_data.dat", "w");    //Fichero para los espines de organización ordenada 1 con T alta

    FILE *fdesord1;
    fdesord1 = fopen("ising_desord1_data.dat", "w");    //Fichero para los espines de organización desordenada con T baja

    FILE *fdesord2;
    fdesord2 = fopen("ising_desord2_data.dat", "w");    //Fichero para los espines de organización desordenada con T alta

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
    for(i = 0; i <= 10000000; i++)
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

    for(i = 0; i <= 10000000; i++)
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

    fclose(ford1);
    fclose(fdesord1);
    fclose(ford2);
    fclose(fdesord2);

    return 0;
}