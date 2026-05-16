#include <iostream>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <complex>
#include <fstream>

#define PI 3.14159265358979323846

using namespace std;

const complex<double> I(0.0, 1.0);

vector<double> InicializarPotencial(int N, int lambda, double k0);
vector<complex<double>> InicializarFuncionDeOnda(int N, double k0);
vector<complex<double>> CalcularAlpha(int N, double s, vector<double> V);
void EvolucionTemporal(int N, double s, vector<double> V, vector<complex<double>> phi, vector<complex<double>> alpha, ofstream &file_evolucion, ofstream &file_norma);


// Función principal.
int main()
{
    int N, nciclos, lambda; //Puntos de red, velocidad inicial y fuerza del potencial
    double k0, s; 

    ofstream file_norma("norma.dat");
    ofstream file_evolucion("evolucion.dat");

    //Inicializo los valores del momento angular reescalado y el paso de tiempo reescalado.
    k0=2*PI*nciclos/N;
    s=1/(4*k0*k0);

    vector<double> V(N+1);

    V=InicializarPotencial(N, lambda, k0); //Inicializo el potencial

    vector<complex<double>> phi(N+1); //Vector de la función de onda

    phi=InicializarFuncionDeOnda(N, k0); //Inicializo la función de onda

    vector<complex<double>> alpha; //Coeficiente alpha

    alpha=CalcularAlpha(N, s, V); //Calculo el coeficiente alpha



    file_norma.close();
    file_evolucion.close();
    return 0;
}

//Funciones utilizadas en el programa

vector<double> InicializarPotencial(int N, int lambda, double k0){
    vector<double> V(N+1);
    for(int i=0; i<=N; i++){
        if(i>=2*N/5 && i<=3*N/5)
            V[i]=lambda*k0*k0;
        else
            V[i]=0;
    }
    return V;
}

vector<complex<double>> InicializarFuncionDeOnda(int N, double k0){
    //Fuerzo las condiciones de contorno
    vector<complex<double>> phi(N+1);
    
    phi[0]=0;
    phi[N]=0;
    
    for(int j=1; j<N; j++){
        //Calculamos la amplitud gaussiana
        double termino_espacial = (4.0 * j - N) / (double)N; 
        double amplitud_gaussiana = exp(-8.0 * termino_espacial * termino_espacial);

        //Calculamos la fase compleja
        complex<double> fase_compleja = exp(I * k0 * (double)j);

        //Asignamos el valor de la función de onda en el punto j
        phi[j] = amplitud_gaussiana * fase_compleja;
    }
    
    return phi;
}   

vector<complex<double>> CalcularAlpha(int N, double s, vector<double> V){
    //Creamos el vector alpha y lo inicializamos con ceros
    vector<complex<double>> alpha(N+1, 0.0);
    //Creo la constante de la diagonal de A
    complex<double> A0;

    //Inicializo un bucle hacia atrás

    for (int j = N - 1; j >= 1; j--) {
        // Calculamos A_j^0 para este punto específico 
        complex<double> A0 = -2.0 + (2.0 * I / s) - V[j];
    
        // Calculamos alpha_{j-1}
        alpha[j - 1] = -1.0 / (A0 + alpha[j]);
    }

    
    return alpha;
}

void EvolucionTemporal(int N, int nciclos, double s, vector<double> V, vector<complex<double>> phi, vector<complex<double>> alpha, ofstream &file_evolucion, ofstream &file_norma){
    vector<complex<double>> beta(N+1);
    vector<complex<double>> chi(N+1);
    vector<complex<double>> b(N+1, 0.0);
    beta[N]=0.0;
    for(int n=1; n<=nciclos; n++){
        for(int j=N-1; j>=0; j--){
            // Calculamos b_j para este tiempo específico 
            b[j+1] = 4.0 * I * phi[j+1] / s;
            //Como A⁺_j = A⁻_j =1 para todo j, gamm[j+1]=alpha[j]
            beta[j] = alpha[j] * (b[j+1] - beta[j + 1]);
        }

        chi[0]=0.0;
        for(int j=1; j<=N; j++){
            chi[j]=alpha[j-1] * chi[j-1] + beta[j-1]; //Calculo la chi
            phi[j]=chi[j]-phi[j]; //Actualizamos la función de onda para el siguiente paso temporal
        }

        double normaphi=0.0;
        //Calculo el modulo de phi al cuadrado
        for(int j=0; j<=N; j++){
            normaphi+=norm(phi[j]);
        }

        
        
    }


    return;
}