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
#include <math.h> // esto es C
#include <fstream>
#include <string>
#include <vector>

using namespace std;

// Dejo un hueco para definir funciones
//Función con el agoritmo de Verlet
void Verlet(double &t, double h, int N, double x0[][2], double v0[][2], double x[][4][2], double v[][4][2], double a[][4][2], double masa[]);
//Función para reescalar los datos de entrada
void reescalar(double &h, double x0[][2], double v0[][2], double masa[]);
//Función para deshacer el reescalado de los datos obtenidos
void deshacer_reescalado(double x[][4][2], double v[][4][2], double a[][4][2], double masa[], int N);
//Función para leer el fichero de datos de entrada
void leer_datos(ifstream &data, double x0[][2], double v0[][2], double masa[]);

//Función para escribir en el fichero de salida
void escribir_datos_pos(ofstream &trayectoriasverlet, double x[][4][2], int N);


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


    //Creo las matrices con las condiciones iniciales.
    double x0[4][2];    //Posiciones iniciales. x0[numero de planetas][coordenada]
    double v0[4][2];    //Velocidades iniciales. v0[numero de planetas][coordenada]

    //Creo matrices tridimensionales que almacenarán las posiciones, velocidades y aceleraciones de los planetas en cada paso de tiempo.
    double x[N][4][2];  //Posiciones. x[numero de pasos][numero de planetas][coordenada]
    double v[N][4][2];  //Velocidades. v[numero de pasos][numero de planetas][coordenada]
    double a[N][4][2];  //Aceleraciones. a[numero de pasos][numero de planetas][coordenada]

    //Necesito también un vector con las masas de cada planeta.
    double masa[4];

    //Leo los datos de entrada del fichero y los almaceno en las matrices correspondientes.
    leer_datos(data, x0, v0, masa);

    //Reescalo los datos de entrada para no trabajar con números muy grandes o muy pequeños.
    reescalar(h, x0, v0, masa);

    //Dado el tiempo inicial, el paso y el número de pasos, genero un vector de tiempos para cada paso.
    //Este vector de tiempo está afectado por el reescalado, que se aplicó al paso de tiempo. Más adelante se deshará el reescalado de este vector de tiempo.
    double tiempo[N];
    for (int i=0; i<N; i++){
        tiempo[i]=i*h;
    };

    //Aplico el algoritmo de Verlet para calcular la posición, velocidad y aceleración de cada planeta en cada paso de tiempo.
    Verlet(t, h, N, x0, v0, x, v, a, masa);

    //Deshago el reescalado de los datos obtenidos.
    deshacer_reescalado(x, v, a, masa, N);

    //Escribo los datos de posición de cada planeta en cada paso de tiempo en el fichero de salida.
    escribir_datos_pos(trayectoriasverlet, x, N);

    //Una vez hecho esto, la simulación de las trayectorias está completa.

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

void leer_datos(ifstream &data, double x0[][2], double v0[][2], double masa[]) {
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

void escribir_datos_pos(ofstream &trayectoriasverlet, double x[][4][2], int N) {
    //Esta función escribe los datos de posición de cada planeta en cada paso de tiempo en el fichero "posiciones_planetas.dat".
    //El formato del fichero es el siguiente:
    //Paso de tiempo 1: posición x del planeta 1, posición y del planeta 1 [salto de línea] posición x del planeta 2, posición y del planeta 2, ...
    int i;
    int j;

    for (i=0; i<N; i++)
    {
        for(j=0; j<4; j++){
            trayectoriasverlet << x[i][j][0] << ", " << x[i][j][1] << endl;
        }
    };

    return;
};

void reescalar(double &h, double x0[][2], double v0[][2], double masa[]){
    double c;
    double G;
    double M;

    c=1.496e11;
    G=6.67430e-11;
    M=2e30;
    //Reescalo todas las posiciones, velocidades y masas.
    for (int i=0; i<4; i++){
        x0[i][0]=x0[i][0]/c;
        x0[i][1]=x0[i][1]/c;
        v0[i][0]=v0[i][0]*pow(c,1.5)/(c*pow(G*M,0.5));
        v0[i][1]=v0[i][1]*pow(c,1.5)/(c*pow(G*M,0.5));
        masa[i]=masa[i]/M;
    };
    h=h*pow(G*M/pow(c,3),0.5);

};

void deshacer_reescalado(double x[][4][2], double v[][4][2], double a[][4][2], double masa[], int N){
    double c;
    double G;
    double M;

    c=1.496e11;
    G=6.67430e-11;
    M=2e30;
    //Deshago el reescalado de todas las posiciones, velocidades, aceleraciones y masas.
    for (int n=0; n<N; n++){
        for (int i=0; i<4; i++){
            for (int j=0; j<2; j++){
                x[n][i][j]=x[n][i][j]*c;
                v[n][i][j]=v[n][i][j]*c*pow(G*M/pow(c,3),-0.5);
                a[n][i][j]=a[n][i][j]*G*M/pow(c,2);
            };
            masa[i]=masa[i]*M;
        }
    };
    return;
};
