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
const double omega = 2.6617e-6; // Velocidad angular de la Luna en torno a la Tierra

// Parámetros reescalados
const double delta = G * MT / pow(dTL, 3);
const double mu = ML / MT;

// Declaración de funciones existentes
void f(double t, const double y[], double dydt[]);
void paso_RK4(double t, double h, double y[], int n, int N);
void integracion_adaptativa(double &t, double &h, double y[], int N);

// NUEVA FUNCIÓN: Calcula la constante de movimiento H' = H - omega * p_phi (adimensional)
double calcular_H_prima(const double y[], double t) {
    double r_tilde = y[0];
    double phi = y[1];
    double pr_tilde = y[2];
    double pphi_tilde = y[3];
    
    double r_prime = sqrt(r_tilde * r_tilde + 1.0 - 2.0 * r_tilde * cos(phi - omega * t));
    
    // Energía total menos el término de rotación
    double H_prima = 0.5 * (pr_tilde * pr_tilde + (pphi_tilde * pphi_tilde) / (r_tilde * r_tilde)) 
                     - delta * (1.0 / r_tilde + mu / r_prime) 
                     - omega * pphi_tilde;
    return H_prima;
}

int main() {
// Condiciones iniciales (NUEVAS: Órbita terrestre a 500 km)
    double altitud = 500e3; // 500 km sobre la superficie
    double r_real = 6.37816e6 + altitud;
    double r0 = r_real / dTL; // Distancia inicial normalizada
    
    // Velocidad orbital (menor a la de escape para quedar atrapado)
    double v0 = 8000.0; // m/s
    double phi0 = 0.0; 
    
    // Momentos adimensionales corregidos físicamente
    double pr0_tilde = 0.0; 
    double pphi0_tilde = r0 * (v0 / dTL); 

    // Guardamos el estado inicial para ambas simulaciones
    double y_inicial[4] = {r0, phi0, pr0_tilde, pphi0_tilde};
    
    // Calculamos el valor analítico constante de H' en t=0
    double H_prima_0 = calcular_H_prima(y_inicial, 0.0);

    // =======================================================
    // PARTE 1: SIMULACIÓN CON PASO FIJO
    // =======================================================
    double y_fijo[4] = {y_inicial[0], y_inicial[1], y_inicial[2], y_inicial[3]};
    double t_fijo = 0.0;
    double h_fijo = 30.0; // 1 minuto (sugerido por el problema)
    int iteraciones_fijo = 0;
    
    ofstream archivo_fijo("trayectoria_fija.dat");
    if (!archivo_fijo.is_open()) return 1;

    cout << "Iniciando simulacion con paso FIJO (h = " << h_fijo << " s)..." << endl;
    
    bool continuar_fijo = true;
    while (continuar_fijo) {
        paso_RK4(t_fijo, h_fijo, y_fijo, 4, 1);
        t_fijo += h_fijo; // Necesario sumar manualmente porque paso_RK4 toma 't' por valor
        iteraciones_fijo++;

        double r_tilde = y_fijo[0];
        double phi = y_fijo[1];
        double r_prime = sqrt(r_tilde * r_tilde + 1.0 - 2.0 * r_tilde * cos(phi - omega * t_fijo));

        // Guardamos puntos cada 100 iteraciones para no generar archivos inmensos
        if (iteraciones_fijo % 100 == 0) {
            double x_nave = r_tilde * cos(phi);
            double y_nave = r_tilde * sin(phi);
            double x_luna = cos(omega * t_fijo);
            double y_luna = sin(omega * t_fijo);
            archivo_fijo << x_nave << "," << y_nave << "\n" << x_luna << "," << y_luna << "\n\n";
        }

        // Condiciones de parada originales
        if (t_fijo > 100.0 && r_tilde * dTL < 6.37816e6) continuar_fijo = false;
        if (r_prime * dTL < 1.7374e6) continuar_fijo = false;
        if (t_fijo > 864000.0) continuar_fijo = false;
    }
    archivo_fijo.close();
    
    // Calculamos el error acumulado de la constante H'
    double error_H_fijo = fabs(calcular_H_prima(y_fijo, t_fijo) - H_prima_0);


    // =======================================================
    // PARTE 2: SIMULACIÓN CON PASO ADAPTATIVO
    // =======================================================
    double y_adap[4] = {y_inicial[0], y_inicial[1], y_inicial[2], y_inicial[3]};
    double t_adap = 0.0;
    double h_adap = 30.0; // Paso inicial sugerido
    int iteraciones_adap = 0;
    
    ofstream archivo_adap("trayectoria_adap.dat");
    if (!archivo_adap.is_open()) return 1;

    cout << "Iniciando simulacion con paso ADAPTATIVO..." << endl;
    
    bool continuar_adap = true;
    while (continuar_adap) {
        integracion_adaptativa(t_adap, h_adap, y_adap, 1);
        iteraciones_adap++; // La función 'integracion_adaptativa' sí actualiza 't_adap' por referencia

        double r_tilde = y_adap[0];
        double phi = y_adap[1];
        double r_prime = sqrt(r_tilde * r_tilde + 1.0 - 2.0 * r_tilde * cos(phi - omega * t_adap));

        if (iteraciones_adap % 100 == 0) {
            double x_nave = r_tilde * cos(phi);
            double y_nave = r_tilde * sin(phi);
            double x_luna = cos(omega * t_adap);
            double y_luna = sin(omega * t_adap);
            archivo_adap << x_nave << "," << y_nave << "\n" << x_luna << "," << y_luna << "\n\n";
        }

        if (t_adap > 100.0 && r_tilde * dTL < 6.37816e6) continuar_adap = false;
        if (r_prime * dTL < 1.7374e6) continuar_adap = false;
        if (t_adap > 864000.0) continuar_adap = false;
    }
    archivo_adap.close();
    
    double error_H_adap = fabs(calcular_H_prima(y_adap, t_adap) - H_prima_0);


    // =======================================================
    // RESULTADOS FINALES Y COMPARACIÓN OBLIGATORIA
    // =======================================================
    cout << "\n==============================================" << endl;
    cout << "  COMPARACION: H FIJA vs H ADAPTATIVA" << endl;
    cout << "==============================================" << endl;
    cout << "METODO FIJO:" << endl;
    cout << "  - Iteraciones (esfuerzo): " << iteraciones_fijo << endl;
    cout << "  - Error acumulado en H':  " << error_H_fijo << endl;
    cout << "\nMETODO ADAPTATIVO:" << endl;
    cout << "  - Iteraciones (esfuerzo): " << iteraciones_adap << endl;
    cout << "  - Error acumulado en H':  " << error_H_adap << endl;
    cout << "==============================================" << endl;
    cout << "Archivos 'trayectoria_fija.dat' y 'trayectoria_adap.dat' generados con exito." << endl;

    return 0;
}

// [MANTENER AQUÍ TUS FUNCIONES ORIGINALES: f, paso_RK4 e integracion_adaptativa]
void f(double t, const double y[], double dydt[]) {
    double r_tilde = y[0];
    double phi = y[1];
    double pr_tilde = y[2];
    double pphi_tilde = y[3];

    double r_prime = sqrt(r_tilde * r_tilde + 1.0 - 2.0 * r_tilde * cos(phi - omega * t));

    dydt[0] = pr_tilde;
    dydt[1] = pphi_tilde / (r_tilde * r_tilde);

    dydt[2] = (pphi_tilde * pphi_tilde) / (r_tilde * r_tilde * r_tilde)
              - delta * (1.0 / (r_tilde * r_tilde)
              + (mu / (r_prime * r_prime * r_prime)) * (r_tilde - cos(phi - omega * t)));

    dydt[3] = - (delta * mu * r_tilde / (r_prime * r_prime * r_prime)) * sin(phi - omega * t);
}

void paso_RK4(double t, double h, double y[], int n, int N) {
    double k1[4], k2[4], k3[4], k4[4], aux[4];

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
    double epsilon_max = 1e-12;

    for (int i = 0; i < 4; i++) {
        y_h[i] = y_h2[i] = y[i];
    }

    paso_RK4(t, h, y_h, 4, N);
    paso_RK4(t, h / 2.0, y_h2, 4, N);
    paso_RK4(t + h / 2.0, h / 2.0, y_h2, 4, N);

    double error = 0.0;
    for (int i = 0; i < 4; i++) {
        double e_i = (16.0 / 15.0) * fabs(y_h2[i] - y_h[i]);
        if (e_i > error) error = e_i;
    }

    double s = pow(error / epsilon_max, 0.2);
    if (s > 2.0) {
        h = h / 2.0;
    } else {
        t += h;
        for (int i = 0; i < 4; i++) y[i] = y_h2[i];
        if (s < 0.1) h = 2.0 * h;
    }
}