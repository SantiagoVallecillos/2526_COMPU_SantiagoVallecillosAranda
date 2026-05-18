#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

using namespace std;

// Constantes físicas (SI)
const double G = 6.67e-11;
const double MT = 5.9736e24;
const double ML = 0.07349e24;
const double dTL = 3.844e8;
const double omega = 2.6617e-6;
const double RT = 6.37816e6;
const double RL = 1.7374e6;

// Parámetros reescalados (Diapositiva 10-11)
const double delta = G * MT / pow(dTL, 3);
const double mu = ML / MT;

void f(double t, const double y[], double dydt[]);
void paso_RK4(double t, double h, double y[], int n, int N);
void integracion_adaptativa(double &t, double &h, double y[], int N);
double calcular_H_prima(double t, double y[]);

int main() {
    
    double t=0.0;
    double h=30.0; //Paso inicial sugerido

    double r0=RT;
    double v0=11.2e3; //Velocidad de escape terrestre
    double phi0=0.0; //Ángulo inicial
    double pr0=0.0; //Momento radial inicial
    double pphi0=v0*r0; //Momento angular inicial

    // 1. Inicialización del vector de estado con las variables reescaladas (adimensionales)
    // Según las ecuaciones, definimos r_tilde y los momentos divididos por las escalas correspondientes.
    double y[4];
    y[0] = r0 / dTL;              // r tilde (distancia normalizada respecto a dTL)
    y[1] = phi0;                  // phi (ángulo)
    y[2] = pr0 / dTL;             // pr tilde (momento radial normalizado)
    y[3] = pphi0 / (dTL * dTL);   // pphi tilde (momento angular normalizado)

    // 2. Apertura del archivo de salida para registrar la trayectoria
    ofstream archivo("trayectoria.dat");
    if (!archivo.is_open()) {
        cerr << "Error al abrir el archivo para guardar los datos." << endl;
        return 1;
    }

    // Número de pasos internos para la función paso_RK4 (frecuencia de evaluación)
    int N = 1; 

    // Guardamos las condiciones iniciales en el archivo
    // Formato: t | x_nave | y_nave | x_luna | y_luna | H_prima | h
    archivo << t << "\t" 
            << y[0] * dTL * cos(y[1]) << "\t" << y[0] * dTL * sin(y[1]) << "\t"
            << dTL * cos(omega * t) << "\t" << dTL * sin(omega * t) << "\t"
            << calcular_H_prima(t, y) << "\t" << h << "\n";

    // 3. Bucle principal de la simulación con control de paso adaptativo
    bool continuar = true;
    while (continuar) {
        // Llamada a la rutina de integración adaptativa que actualiza t, h e y[]
        integracion_adaptativa(t, h, y, N);

        // Calculamos las coordenadas actuales para verificar las condiciones de parada
        double r_tilde = y[0];
        double phi = y[1];
        // Distancia reescalada de la nave a la Luna (r')
        double r_prime = sqrt(r_tilde * r_tilde + 1.0 - 2.0 * r_tilde * cos(phi - omega * t));

        // Guardamos los datos de la iteración actual en el archivo
        archivo << t << "\t" 
                << r_tilde * dTL * cos(phi) << "\t" << r_tilde * dTL * sin(phi) << "\t" 
                << dTL * cos(omega * t) << "\t" << dTL * sin(omega * t) << "\t" 
                << calcular_H_prima(t, y) << "\t" << h << "\n";

        // Condición de parada A: Colisión o reentrada en la Tierra 
        // Se añade (t > 100.0) para evitar que el programa termine en el instante inicial t=0
        if (t > 100.0 && r_tilde * dTL < RT) {
            cout << "La nave ha regresado o colisionado con la Tierra a los " << t << " segundos." << endl;
            continuar = false;
        }

        // Condición de parada B: Colisión con la superficie de la Luna
        if (r_prime * dTL < RL) {
            cout << "¡La nave ha colisionado con la superficie lunar a los " << t << " segundos!" << endl;
            continuar = false;
        }

        // Condición de parada C: Tiempo máximo de seguridad (ej. ~10 días) para evitar bucles infinitos
        if (t > 864000.0) {
            cout << "Se ha alcanzado el tiempo límite de la simulación sin colisiones." << endl;
            continuar = false;
        }
    }

    // Cerrar el flujo del archivo
    archivo.close();
    cout << "Simulación finalizada de manera exitosa. Resultados guardados en 'trayectoria.dat'." << endl;

    return 0;
}

void f(double t, const double y[], double dydt[]) {
    double r_tilde = y[0];
    double phi = y[1];
    double pr_tilde = y[2];
    double pphi_tilde = y[3];

    // Distancia nave-Luna reescalada (r')
    double r_prime = sqrt(r_tilde*r_tilde + 1.0 - 2.0*r_tilde*cos(phi - omega*t));

    // Ecuaciones de movimiento 
    dydt[0] = pr_tilde;
    dydt[1] = pphi_tilde / (r_tilde * r_tilde);
    
    double common_term = delta * mu / pow(r_prime, 3);
    
    dydt[2] = (pow(pphi_tilde, 2) / pow(r_tilde, 3)) - 
              delta * (1.0/(r_tilde*r_tilde) + (mu/pow(r_prime,3)) * (r_tilde - cos(phi - omega*t)));
              
    dydt[3] = - (delta * mu * r_tilde / pow(r_prime, 3)) * sin(phi - omega*t);
}

void paso_RK4(double t, double h, double y[], int n, int N) {
    double k1[4], k2[4], k3[4], k4[4], aux[4];

    // Realizamos N pasos RK4 para no tener que llamar a la función constantemente.

    for (int step = 0; step < N; step++) {
        f(t, y, k1);
        for (int i = 0; i < n; i++) aux[i] = y[i] + 0.5 * h * k1[i];

        f(t + 0.5 * h, aux, k2);
        for (int i = 0; i < n; i++) aux[i] = y[i] + 0.5 * h * k2[i];

        f(t + 0.5 * h, aux, k3);
        for (int i = 0; i < n; i++) aux[i] = y[i] + h * k3[i];

        f(t + h, aux, k4);

        for (int i = 0; i < n; i++) {
            y[i] += (h / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
        }

        t += h;
    }
}

void integracion_adaptativa(double &t, double &h, double y[], int N) {
    double y_h[4], y_h2[4];
    double epsilon_max = 1e-12; // Tolerancia sugerida

    // 1. Copiar estado actual
    for(int i=0; i<4; i++) y_h[i] = y_h2[i] = y[i];

    // 2. Calcular un paso de tamaño h
    paso_RK4(t, h, y_h, 4, N);

    // 3. Calcular dos pasos de tamaño h/2
    paso_RK4(t, h/2.0, y_h2, 4, N);
    paso_RK4(t + h/2.0, h/2.0, y_h2, 4, N);

    // 4. Estimar error (Ecuación diapositiva 10)
    double error = 0;
    for(int i=0; i<4; i++) {
        double e_i = (16.0/15.0) * fabs(y_h2[i] - y_h[i]);
        if(e_i > error) error = e_i;
    }

    // 5. Ajustar h según el algoritmo de la Diapositiva 9
    double s = pow(error / epsilon_max, 0.2);
    
    if (s > 2.0) {
        h = h / 2.0; // Error muy grande, reducir paso y repetir
    } else {
        t += h;
        for(int i=0; i<4; i++) y[i] = y_h2[i]; // Aceptar paso
        if (s < 0.1) h = 2.0 * h; // Error muy pequeño, aumentar paso
    }
}

double calcular_H_prima(double t, double y[]) {
    double r = y[0] * dTL; 
    double phi = y[1];
    double pr = y[2] * MT * dTL; // (Asumiendo masa nave m=1 para el chequeo)
    double pphi = y[3] * MT * dTL * dTL;
    
    double rL = sqrt(r*r + dTL*dTL - 2*r*dTL*cos(phi - omega*t));
    
    double H = (pr*pr)/2.0 + (pphi*pphi)/(2.0*r*r) - G*MT/r - G*ML/rL;
    return H - omega * pphi;
}