#include <iostream>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <complex>
#include <fstream>
#include <string>
 
#define PI 3.14159265358979323846

using namespace std;

const complex<double> I(0.0, 1.0);

typedef complex<double> dcomplex; //Alivio para la notación

double norma(const vector<complex<double>>& v);
vector<double> potencial(int N, double lambda, double k0);
vector<complex<double>> calcular_phi0(int N, double k0);
void simular_caso(int N, double lambda, string nombre_datos, string nombre_norma);



int main()
{
    //Estudio de N (Resolución)
    simular_caso(700, 0.7, "datos_N700_L0.7.dat", "norma_N700_L0.7.dat");
    simular_caso(3000, 0.7, "datos_N3000_L0.7.dat", "norma_N3000_L0.7.dat");

    //Estudio de Lambda (Efecto túnel y reflexión)
    simular_caso(1000, 0.2, "datos_N1000_L0.2.dat", "norma_N1000_L0.2.dat");
    simular_caso(1000, 15.0, "datos_N1000_L15.dat", "norma_N1000_L15.dat");

    //Caso intermedio
    simular_caso(1000, 1.0, "datos_N1000_L1.0.dat", "norma_N1000_L1.0.dat");
    
    return 0;
}

//Código basado en el de mi compañero Alberto Jiménez, con algunas modificaciones de notación y optimizacion.

//Función para calcular la norma de un vector de números complejos
double norma(const vector<complex<double>>& v) {
    double sum = 0.0;
    for (const auto& elem : v) {
        sum += norm(elem); // norm(elem) devuelve el módulo al cuadrado de elem
    }
    return sqrt(sum);
}

//Función para calcular el potencial Vj
vector<double> potencial(int N, double lambda, double k0){
    //La barrera está definida entre 2N/5 y 3N/5. En el resto de puntos se anula
    vector<double> Vj(N, 0.0);
    
    for(int i=0; i<N; i++){
        if(i>=2*N/5 && i<=3*N/5){
            Vj[i]=lambda;
        } else {
            Vj[i]=0.0;
        }
    }
    return Vj;
}

//Función para inicializar la función de onda
vector<complex<double>> calcular_phi0(int N, double k0){
    vector<complex<double>> phi0(N);
    
    for(int j=0; j<N; j++){
        double exponente = -8.0 * (4.0*j-N)*(4.0*j-N) / (N*N);
        complex<double> fase(0,k0*j);
        phi0[j] = exp(exponente) * exp(fase);
    }
    return phi0;
}

//Función que realiza la simulación para los casos indicados.
void simular_caso(int N, double lambda, string nombre_datos, string nombre_norma){
    cout << "Simulando caso N=" << N << ", lambda=" << lambda << "..." << endl;

    //Parámetros de la simulación. Se pueden modificar para experimentar con diferentes condiciones iniciales y potenciales.
    const int n_ciclos = N/8; //Define cuán rápido se mueve al inicio la función de onda. N/8 nos debe dar una velocidad razonable.

    //El valor lambda representa la altura de la barrera de potencial y afecta a la facilidad de atravesarla.
    const int pasos_tiempo = 2000;  //Número total de pasos de tiempo a simular.
    const int n_D = 10;             //Guardamos datos cada n_D pasos para no saturar el archivo de salida.

    //Variables derivadas (cálculos previos)
    double k0_tilde = 2.0*PI*n_ciclos/N;
    double s_tilde = 1/(4*k0_tilde*k0_tilde);   //Este es el valor que dan los apuntes. Alberto tiene 0.5. Si explota el programa, lo cambio a lo que tiene Alberto.

    dcomplex i_unidad(0.0, 1.0);

    //Vectores que utilizaremos

    vector<dcomplex> phi(N); //Función de onda inicial
    vector<double> V(N+1); //Potencial Vj

    //Los siguientes vectores se utilizan en el método de Crank_Nicolson. Alpha y gamma se calculan hacia atrás
    //y beta es e coeficiente dinámico que se actualiza en cada paso de tiempo.

    vector<dcomplex> alpha(N);  //Coeficientes alpha estáticos
    vector<dcomplex> beta(N);   //Coeficientes beta dinámicos
    vector<dcomplex> gamma(N+1);//Coeficientes gamma auxiliares
    vector<dcomplex> chi(N+1);  //Vector chi intermedio para el paso temporal.

    //Archivos de salida
    ofstream file_data(nombre_datos); //Módulo de función de onda: x, probabilidad. Cada n_D pasos de tiempo, salto de línea doble que es un nuevo frame temporal
    ofstream file_norma(nombre_norma); //Norma total de la función de onda cada n_D pasos de tiempo. Debería ser constante.

    //Inicialización de los vectores del potencial y la función de onda.
    phi = calcular_phi0(N, k0_tilde);
    V = potencial(N, lambda, k0_tilde);

    //Normalización a 1
    double norma_inicial = norma(phi);
    for(int j=0; j<N; j++){
        phi[j] /= norma_inicial;
    }

    //Cálculo de los coeficientes alpha y gamma estáticos. Se pueden precomputar porque no dependen del tiempo.
    alpha[N-1] = 0.0; //Condición de frontera
    for(int j=N-1;j>0;j--){
        dcomplex A0=dcomplex(-2.0,0.0)+(2.0*i_unidad/s_tilde)-V[j];
        gamma[j]=1.0/(A0+alpha[j]);
        alpha[j-1]=-gamma[j];
    }

    //Bucle principal para la simulación temporal
    for(int n=0; n<pasos_tiempo; n++){
        //Cálculo de los coeficientes beta hacia atrás
        beta[N-1]=0.0; //Condición de frontera
        for(int j=N-1; j>0; j--){
            dcomplex b_j=(4.0*i_unidad/s_tilde)*phi[j];
            beta[j-1]=gamma[j]*(b_j-beta[j]);
        }

        //Cálculo del vector chi hacia adelante
        chi[0]=0.0; //Condición de frontera
        chi[N]=0.0; //Condición de frontera
        for(int j=0; j<N; j++){
            chi[j+1]=alpha[j]*chi[j]+beta[j];
        }

        //Actualización de la función de onda y cálculo de la norma
        for(int j=0; j<=N; j++){
            phi[j]=chi[j]-phi[j];
        }

        //Guardado de datos. Para la norma total se hace cada n_D pasos.
        if(n%n_D==0){
            double norma_actual = norma(phi);
            file_norma << n << ", " << norma_actual << endl;

            for(int j=0; j<=N; j++){
                file_data << j << ", " << norm(phi[j]) << endl;
            }

            file_data << endl << endl; //Salto de línea para indicar nuevo frame.
        }
        file_data.close();
        file_norma.close();

    }

}