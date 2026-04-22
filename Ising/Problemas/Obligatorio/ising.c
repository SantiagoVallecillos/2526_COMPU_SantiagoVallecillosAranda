#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>    //Para generar números aleatorios
#define e 2.71828
#define N 20            //Tamaño de las matrices
double min(double x);

gsl_rng *tau;

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
    extern gsl_rng *tau;
    int semilla = 310809;
    tau = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(tau, semilla);

    //Inicializamos variables
    int i, j, k, n, m, sord[N][N], sdesord[N][N], sN[64][64], randn, randm;
    double dseta, P_S, T[2], T_m[10], E_S, p, rand, deltaE, mN16, mN32, mN64;

    FILE *ford1;
    ford1 = fopen("ising_1.txt", "w");    //Fichero para los espines de organización ordenada 1 con T baja

    FILE *ford2;
    ford2 = fopen("ising_2.txt", "w");    //Fichero para los espines de organización ordenada 1 con T alta

    FILE *fdesord1;
    fdesord1 = fopen("ising_desord1.txt", "w");    //Fichero para los espines de organización desordenada con T baja

    FILE *fdesord2;
    fdesord2 = fopen("ising_desord2.txt", "w");    //Fichero para los espines de organización desordenada con T alta

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

            rand = gsl_rng_uniform(tau);
            if(rand <= 0.5){sdesord[n][m] = 1;}
            else if(rand > 0.5){sdesord[n][m] = -1;} 
        }
        //Lo metemos todo en los ficheros 
        for(m = 0; m < N; m++)
        {
            fprintf(ford1, "%i\t%i\t%i\t", n, m, sord[n][m]);
            fprintf(fdesord1, "%i\t%i\t%i\t", n, m, sdesord[n][m]);
        }
    }
    fprintf(ford1, "\n");
    fprintf(fdesord1, "\n");

    //Hacemos el bucle para ver si cambian de espin
    for(i = 0; i <= 100000000; i++)
    {
        randn = gsl_rng_uniform_int(tau, N-1);
        randm = gsl_rng_uniform_int(tau, N-1);
        //printf("%i\t%i\n", randn, randm);

        //Para configuración ordenada 1
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
        p = minimo(exp((-deltaE)/T[0]));

        dseta = gsl_rng_uniform(tau);
        if(dseta < p){sord[randn][randm] = -sord[randn][randm];}
        else{continue;}

        if(i%(N*N) == 0)
        {
            //Lo volvemos a meter todo en el fichero
            for(n = 0; n < N; n++)
            { 
                for(m = 0; m < N; m++)
                {
                    fprintf(ford1, "%i\t%i\t%i\t", n, m, sord[n][m]);
                }
            }
            fprintf(ford1, "\n");
        }

        //Para configuración desordenada
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
            else if(randm == N-1){deltaE = 2*(sdesord[randn][randm])*(sdesord[randn+1][randm] + sdesord[randn-1][randm] + sdesord[randn][1] + sdesord[randn][randm-1]);}
            else{deltaE = 2*(sdesord[randn][randm])*(sdesord[randn+1][randm] + sdesord[randn-1][randm] + sdesord[randn][randm+1] + sdesord[randn][randm-1]);}
        }

        p = minimo(exp((-deltaE)/T[0]));

        dseta = gsl_rng_uniform(tau);
        if(dseta < p){sdesord[randn][randm] = -sdesord[randn][randm];}
        else{continue;}

        if(i%(N*N) == 0)
        {
            //Lo volvemos a meter todo en el fichero
            for(n = 0; n < N; n++)
            { 
                for(m = 0; m < N; m++)
                {
                    fprintf(fdesord1, "%i\t%i\t%i\t", n, m, sdesord[n][m]);
                }
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

            rand = gsl_rng_uniform(tau);
            if(rand <= 0.5){sdesord[n][m] = 1;}
            else if(rand > 0.5){sdesord[n][m] = -1;} 
        }
        //Lo metemos todo en los ficheros 
        for(m = 0; m < N; m++)
        {
            fprintf(ford2, "%i\t%i\t%i\t", n, m, sord[n][m]);
            fprintf(fdesord2, "%i\t%i\t%i\t", n, m, sdesord[n][m]);
        }
    }
    fprintf(ford2, "\n");
    fprintf(fdesord2, "\n");

    for(i = 0; i <= 100000000; i++)
    {
        randn = gsl_rng_uniform_int(tau, N-1);
        randm = gsl_rng_uniform_int(tau, N-1);

        //Para configuración ordenada 1
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

        p = minimo(exp((-deltaE)/T[1]));

        dseta = gsl_rng_uniform(tau);
        if(dseta < p){sord[randn][randm] = -sord[randn][randm];}
        else{continue;}

        if(i%(N*N) == 0)
        {
            //Lo volvemos a meter todo en el fichero
            for(n = 0; n < N; n++)
            { 
                for(m = 0; m < N; m++)
                {
                    fprintf(ford2, "%i\t%i\t%i\t", n, m, sord[n][m]);
                }
            }
            fprintf(ford2, "\n");
        }

        //Para configuración desordenada
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
        

        p = minimo(exp((-deltaE)/T[1]));

        dseta = gsl_rng_uniform(tau);
        if(dseta < p){sdesord[randn][randm] = -sdesord[randn][randm];}
        else{continue;}

        if(i%(N*N) == 0)
        {
            //Lo volvemos a meter todo en el fichero
            for(n = 0; n < N; n++)
            { 
                for(m = 0; m < N; m++)
                {
                    fprintf(fdesord2, "%i\t%i\t%i\t", n, m, sdesord[n][m]);
                }
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
        for(j = 0; j <= 100000; j++)
        {
            randn = gsl_rng_uniform_int(tau, 64);
            randm = gsl_rng_uniform_int(tau, 64);

            deltaE = 2*(sN[randn][randm])*(sN[randn+1][randm] + sN[randn-1][randm] + sN[randn][randm+1] + sN[randn][randm-1]);

            p = minimo(exp((-deltaE)/T_m[i]));

            dseta = gsl_rng_uniform(tau);
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