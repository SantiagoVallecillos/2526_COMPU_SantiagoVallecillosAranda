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

void f(double t, const double y[], double dydt[]);
void paso_RK4(double t, double h, double y[], int n, int N);
void integracion_adaptativa(double &t, double &h, double y[], int N);

int main() {
    double t = 0.0;
    double h = 30.0; // Paso inicial sugerido

    double r0 = 6.37816e6 / dTL; // Distancia inicial normalizada desde la Tierra
    double v0 = 11.2e3; // Velocidad de escape terrestre
    double phi0 = 0.0; // Ángulo inicial
    double pr0 = 0.0; // Momento radial inicial
    double pphi0 = v0 * r0 * dTL / dTL; // Momento angular normalizado

    double y[4];
    y[0] = r0;              // r tilde (distancia normalizada respecto a dTL)
    y[1] = phi0;            // phi (ángulo alrededor de la Tierra)
    y[2] = pr0 / dTL;       // pr tilde (momento radial normalizado)
    y[3] = pphi0 / (dTL * dTL); // pphi tilde (momento angular normalizado)

    ofstream archivo("trayectoria.dat");
    if (!archivo.is_open()) {
        cerr << "Error al abrir el archivo para guardar los datos." << endl;
        return 1;
    }

    int N = 1;
    double x_nave = y[0] * cos(y[1]);
    double y_nave = y[0] * sin(y[1]);
    double x_luna = cos(omega * t);
    double y_luna = sin(omega * t);

    archivo << x_nave << "," << y_nave << "\n"
            << x_luna << "," << y_luna << "\n\n";

    bool continuar = true;
    while (continuar) {
        integracion_adaptativa(t, h, y, N);

        double r_tilde = y[0];
        double phi = y[1];
        double r_prime = sqrt(r_tilde * r_tilde + 1.0 - 2.0 * r_tilde * cos(phi - omega * t));

        double x_nave = r_tilde * cos(phi);
        double y_nave = r_tilde * sin(phi);
        double x_luna = cos(omega * t);
        double y_luna = sin(omega * t);

        archivo << x_nave << "," << y_nave << "\n"
                << x_luna << "," << y_luna << "\n\n";

        if (t > 100.0 && r_tilde * dTL < 6.37816e6) {
            cout << "La nave ha regresado o colisionado con la Tierra a los " << t << " segundos." << endl;
            continuar = false;
        }

        if (r_prime * dTL < 1.7374e6) {
            cout << "¡La nave ha colisionado con la superficie lunar a los " << t << " segundos!" << endl;
            continuar = false;
        }

        if (t > 864000.0) {
            cout << "Se ha alcanzado el tiempo límite de la simulación sin colisiones." << endl;
            continuar = false;
        }
    }

    archivo.close();
    cout << "Simulación finalizada de manera exitosa. Resultados guardados en 'trayectoria.dat'." << endl;
    return 0;
}

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
