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

using namespace std;

// Dejo un hueco para definir variables


// Dejo un hueco para definir funciones


int main() {

    //Añado los ficheros
    ifstream data("condiciones_iniciales.txt");
    ofstream trayectoriasverlet("posiciones_planetas.dat");

    return 0;
}


