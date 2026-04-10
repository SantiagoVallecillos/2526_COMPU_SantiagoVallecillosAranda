/***********ALGORITMO DE VERLET PARA SIMULACIÓN DEL SISTEMA SOLAR***********/

/*
Este programa tiene como objetivo ejecutar el algoritmo de Verlet para simular el movimiento de los planetas en el sistema solar. 
En primer lugar, se tomarán los datos de un fichero donde se colocarán las condiciones iniciales y datos sobre los planetas.
Los datos requeridos son la distancia entre el sol y el planeta, su masa y su velocidad inicial.
Una vez recopilados estos datos, a continuación se hace un reescalado para no trabajar con una combinación de números grandes
y pequeños, lo que podría generar problemas de precisión en los cálculos.
A continuación, se implementará el algoritmo de Verlet para calcular la posición y velocidad de cada planeta en cada paso de tiempo. 
Luego, se deshace el reescalado de los datos obtenidos.
Finalmente, se almacenarán en un nuevo fichero para, posteriormente, poder animar la trayectoria de los planetas utilizando Python.
*/


#include <iostream>
#include <math.h>
#include <fstream>
#include <string>
#include <vector>

using namespace std;

// Dejo un hueco para definir funciones
//Función con el agoritmo de Verlet
void Verlet(double &t, double h, int N, double x0[][2], double v0[][2], double x[][4][2], double v[][4][2], double a[][4][2], double masa[]);
//Función para reescalar los datos de entrada
void reescalar(double x0[][2], double v0[][2], double masa[]);
//Función para deshacer el reescalado de los datos obtenidos
void deshacer_reescalado(double x[][4][2], double v[][4][2], double a[][4][2], double masa[]);
//Función para leer el fichero de datos de entrada
void leer_datos(ifstream data, double x0[][2], double v0[][2], double masa[]);

//Función para escribir en el fichero de salida
void escribir_datos(ofstream trayectoriasverlet, double x[][4][2], double v[][4][2], double a[][4][2], double t[], int N);


int main() {

    //Añado los ficheros
    ifstream data("condiciones_iniciales.txt");
    ofstream trayectoriasverlet("posiciones_planetas.dat");

    //Añado las variables del tiempo y el paso
    double t;
    double h;
    //Añado el número de pasos que queremos realizar
    int N;
    N=1000;

    //Añado el tiempo inicial y el paso
    t=0.0;
    h=0.01;

    //Dado el tiempo inicial, el paso y el número de pasos, genero un vector de tiempos para cada paso.
    double tiempo[N];
    for (int i=0; i<N; i++){
        tiempo[i]=i*h;
    };

    //Creo las matrices con las condiciones iniciales.
    double x0[4][2];    //Posiciones iniciales. x0[numero de planetas][coordenada]
    double v0[4][2];    //Velocidades iniciales. v0[numero de planetas][coordenada]

    //Creo matrices tridimensionales que almacenarán las posiciones, velocidades y aceleraciones de los planetas en cada paso de tiempo.
    double x[N][4][2];  //Posiciones. x[numero de pasos][numero de planetas][coordenada]
    double v[N][4][2];  //Velocidades. v[numero de pasos][numero de planetas][coordenada]
    double a[N][4][2];  //Aceleraciones. a[numero de pasos][numero de planetas][coordenada]

    //Necesito también un vector con las masas de cada planeta.
    double masa[4];

    return 0;
}

//Aquí escribo el código de las funciones.
void Verlet(double &t, double h, int N, double x0[][2], double v0[][2], double x[][4][2], double v[][4][2], double a[][4][2], double masa[]) {

    //Esta función ignora los reescalados, por lo que se asume que los datos de entrada ya han sido reescalados
    //y se deshará el reescalamiento a posteriori.

    //Para implementar el algoritmo, necesito la constante de gravitación universal, que es la siguiente:
    double G=6.67430e-11;
    //También necesitaré un vector auxiliar de velocidad angular.
    double w[4][2]; //Velocidad angular de cada planeta en cada coordenada.
    //Me aseguro de que t=0 al inicio del algoritmo.
    t=0.0;

    //En primer lugar, se añaden los datos iniciales a las matrices de posiciones y velocidades.
    for (int i=0; i<4; i++){
        for (int j=0; j<2; j++){
            x[0][i][j]=x0[i][j];
            v[0][i][j]=v0[i][j];
        }
    };

    //A continuación comienzo con el algoritmo de Verlet.
    for (int n=0; n<N-1;n++){         //Itero para cada paso de tiempo
        for (int i=0; i<4; i++){    //Itero para cada planeta
            for(int j=0; j<2; j++){ //Itero para cada componente de posición, velocidad y aceleración.
                //En primer lugar, calculo la aceleración de cada planeta en el paso t.
                a[n][i][j]=0.0; //Inicializo la aceleración a cero.
                a[n+1][i][j]=0.0; //Inicializo la aceleración en el paso t+h a cero.
                for (int k = 0; k < 4; k++) {
                    if (k == i) continue; // Saltamos el propio planeta
    
                    // Calculamos la distancia (mejor fuera del bucle 'j' para ahorrar CPU)
                    double dx = x[n][i][0] - x[n][k][0];
                    double dy = x[n][i][1] - x[n][k][1];
                    double dist3 = pow(dx*dx + dy*dy, 1.5);
    
                    // Acumulamos la aceleración con +=
                    a[n][i][j] += -G * masa[k] * (x[n][i][j] - x[n][k][j]) / dist3;
}
                //Una vez tengo la aceleración, calculo la posición en el paso t+h, utilizando la fórmula de Verlet.
                x[n+1][i][j]=x[n][i][j]+v[n][i][j]*h+0.5*a[n][i][j]*pow(h,2);
                //También calculo una velocidad angular que necesitaré más adelante.
                w[i][j]=v[n][i][j]+0.5*a[n][i][j]*h;
                //Calculo la aceleración en el paso t+h, utilizando la posición que acabo de calcular.
                for (int k = 0; k < 4; k++) {
                    if (k == i) continue; // Saltamos el propio planeta
    
                    // Calculamos la distancia (mejor fuera del bucle 'j' para ahorrar CPU)
                    double dx = x[n+1][i][0] - x[n+1][k][0];
                    double dy = x[n+1][i][1] - x[n+1][k][1];
                    double dist3 = pow(dx*dx + dy*dy, 1.5);
    
                    // Acumulamos la aceleración con +=
                    a[n+1][i][j] += -G * masa[k] * (x[n+1][i][j] - x[n+1][k][j]) / dist3;
                };
                //Calculo la velocidad en el paso t+h, utilizando la fórmula de Verlet.
                v[n+1][i][j]=w[i][j]+0.5*a[n+1][i][j]*h;
                //Por último, actualizo el tiempo.
                t+=h;
            }
        };

    };
    //Con esto ya tengo actualizados los vectores de datos de posición, velocidad y aceleración de cada planeta en cada paso de tiempo.
    return;
    
};

void leer_datos(ifstream data, double x0[][2], double v0[][2], double masa[]) {
    //Esta función lee los datos de entrada del fichero "condiciones_iniciales.txt" y los almacena en las matrices correspondientes.
    //El formato del fichero es el siguiente:
    //Planeta 1: masa, posición x, posición y, velocidad x, velocidad y
    //Planeta 2: masa, posición x, posición y, velocidad x, velocidad y
    //Planeta 3: masa, posición x, posición y, velocidad x, velocidad y
    //Planeta 4: masa, posición x, posición y, velocidad x, velocidad y

    for (int i=0; i<4; i++){
        data >> masa[i] >> x0[i][0] >> x0[i][1] >> v0[i][0] >> v0[i][1];
    };

    return;
};
//EDITAR MÁS ADELANTE ESCRIBIR DATOS
void escribir_datos(ofstream trayectoriasverlet, double x[][4][2], double v[][4][2], double a[][4][2], double t[], int N) {
    //Esta función escribe los datos de posición, velocidad y aceleración de cada planeta en cada paso de tiempo en el fichero "posiciones_planetas.dat".
    //El formato del fichero es el siguiente:
    //Paso de tiempo 1: tiempo, posición x planeta 1, posición y planeta 1, velocidad x planeta 1, velocidad y planeta 1, aceleración x planeta 1, aceleración y planeta 1, posición x planeta 2, posición y planeta 2, velocidad x planeta 2, velocidad y planeta 2, aceleración x planeta 2, aceleración y planeta 2, posición x planeta 3, posición y planeta 3, velocidad x planeta 3, velocidad y planeta 3, aceleración x planeta 3, aceleración y planeta 3, posición x planeta 4, posición y planeta 4, velocidad x planeta 4, velocidad y planeta 4, aceleración x planeta 4, aceleración y planeta 4
    //Paso de tiempo 2: tiempo, posición x planeta 1, posición y planeta 1, velocidad x planeta 1, velocidad y planeta 1, aceleración x planeta 1, aceleración y planeta 1, posición x planeta 2, posición y planeta 2, velocidad x planeta 2, velocidad y planeta 2, aceleración x planeta 2, aceleración y planeta 2, posición x planeta 3, posición y planeta 3, velocidad x planeta 3, velocidad y planeta 3, aceleración x planeta 3, aceleración y planeta 3, posición x planeta 4, posición y planeta 4, velocidad x planeta 4, velocidad y planeta 4, aceleración x planeta 4, aceleración y planeta 4
    //...
    //Paso de tiempo N: tiempo, posición x planeta 1, posición y planeta 1, velocidad x planeta 1, velocidad y planeta 1, aceleración x planeta 1, aceleración y planeta 1, posición x planeta 2, posición y planeta 2, velocidad x planeta 2, velocidad y planeta 2, aceleración x planeta 2, aceleración y planoeta

    for (int n=0; n<N; n++){
        trayectoriasverlet << t[n] << " ";
        for (int i=0; i<4; i++){
            for (int j=0; j<2; j++){
                trayectoriasverlet << x[n][i][j] << " " << v[n][i][j] << " " << a[n][i][j] << " ";
            }
        };
        trayectoriasverlet << endl;
    };
    return;
};

