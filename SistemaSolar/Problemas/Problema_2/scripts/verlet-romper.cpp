#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

using namespace std;

constexpr size_t kNumPlanetas = 9;
constexpr size_t kNumCoordenadas = 2;
constexpr double kG = 6.67430e-11;
constexpr double kUnidadDistancia = 1.496e11;
constexpr double kUnidadMasa = 2.0e30;

using Vector2D = array<double, kNumCoordenadas>;
using PlanetArray = array<Vector2D, kNumPlanetas>;
using Trajectory = vector<PlanetArray>;
using EnergyArray = vector<array<double, kNumPlanetas>>;
using MassArray = array<double, kNumPlanetas>;

void leer_datos(ifstream& data, PlanetArray& x0, PlanetArray& v0, MassArray& masa);
void reescalar(double& h, PlanetArray& x0, PlanetArray& v0, MassArray& masa);
void deshacer_reescalado(Trajectory& x, Trajectory& v, Trajectory& a, MassArray& masa);
void Verlet(double& t, double h, int N, const PlanetArray& x0, const PlanetArray& v0, Trajectory& x, Trajectory& v, Trajectory& a, const MassArray& masa);
PlanetArray calcular_aceleraciones(const PlanetArray& posiciones, const MassArray& masa);
void escribir_datos(ofstream& out, const Trajectory& x);
void escribir_datos_energia(ofstream& out, const EnergyArray& data);
void escribir_datos_periodo(ofstream& out, const array<double, kNumPlanetas>& periodos);
void invariantes(const Trajectory& x, const Trajectory& v, const Trajectory& a, const MassArray& masa, EnergyArray& E, Trajectory& L, Trajectory& p, EnergyArray& mod_p);
void periodos(const EnergyArray& E, const MassArray& masa, array<double, kNumPlanetas>& periodos);
void convertir_periodo_a_dias(array<double, kNumPlanetas>& periodos);
int main() {
    ifstream data("condiciones_iniciales.txt");
    if (!data) {
        cerr << "Error: no se pudo abrir condiciones_iniciales.txt\n";
        return 1;
    }

    ofstream trayectorias("posiciones_planetas.dat");
    if (!trayectorias) {
        cerr << "Error: no se pudo crear posiciones_planetas.dat\n";
        return 1;
    }

    ofstream velocidades("velocidades_planetas.dat");
    if (!velocidades) {
        cerr << "Error: no se pudo crear velocidades_planetas.dat\n";
        return 1;
    }

    ofstream aceleraciones("aceleraciones_planetas.dat");
    if (!aceleraciones) {
        cerr << "Error: no se pudo crear aceleraciones_planetas.dat\n";
        return 1;
    }

    ofstream momento_angular("momento_angular.dat");
    if (!momento_angular) {
        cerr << "Error: no se pudo crear momento_angular.dat\n";
        return 1;
    }

    ofstream energia("energia.dat");
    if (!energia) {
        cerr << "Error: no se pudo crear energia.dat\n";
        return 1;
    }

    ofstream momento_lineal("momento_lineal.dat");
    if (!momento_lineal) {
        cerr << "Error: no se pudo crear momento_lineal.dat\n";
        return 1;
    }

    ofstream periodo_file("periodos.dat");
    if (!periodo_file) {
        cerr << "Error: no se pudo crear periodos.dat\n";
        return 1;
    }

    constexpr int N = 1000;
    double t = 0.0;
    double h = 3600;

    PlanetArray x0{};
    PlanetArray v0{};
    MassArray masa{};

    leer_datos(data, x0, v0, masa);
    reescalar(h, x0, v0, masa);

    vector<double> tiempo(N);
    for (int i = 0; i < N; ++i) {
        tiempo[i] = i * h;
    }

    Trajectory x(N, PlanetArray{});
    Trajectory v(N, PlanetArray{});
    Trajectory a(N, PlanetArray{});
    EnergyArray E(N);
    Trajectory L(N, PlanetArray{});
    Trajectory p(N, PlanetArray{});
    EnergyArray mod_p(N);
    array<double, kNumPlanetas> periodo{};

    Verlet(t, h, N, x0, v0, x, v, a, masa);
    escribir_datos(trayectorias, x);
    escribir_datos(velocidades, v);
    escribir_datos(aceleraciones, a);

    deshacer_reescalado(x, v, a, masa);
    invariantes(x, v, a, masa, E, L, p, mod_p);
    periodos(E, masa, periodo);
    convertir_periodo_a_dias(periodo);

    escribir_datos(momento_angular, L);
    escribir_datos_energia(momento_lineal, mod_p);
    escribir_datos_energia(energia, E);
    escribir_datos_periodo(periodo_file, periodo);

    return 0;
}

void leer_datos(ifstream& data, PlanetArray& x0, PlanetArray& v0, MassArray& masa) {
    for (size_t i = 0; i < kNumPlanetas; ++i) {
        data >> masa[i] >> x0[i][0] >> x0[i][1] >> v0[i][0] >> v0[i][1];
    }
}

void reescalar(double& h, PlanetArray& x0, PlanetArray& v0, MassArray& masa) {
    for (size_t i = 0; i < kNumPlanetas; ++i) {
        x0[i][0] /= kUnidadDistancia;
        x0[i][1] /= kUnidadDistancia;
        v0[i][0] *= pow(kUnidadDistancia, 1.5) / (kUnidadDistancia * sqrt(kG * kUnidadMasa));
        v0[i][1] *= pow(kUnidadDistancia, 1.5) / (kUnidadDistancia * sqrt(kG * kUnidadMasa));
        masa[i] /= kUnidadMasa;
    }

    h *= sqrt(kG * kUnidadMasa / (kUnidadDistancia * kUnidadDistancia * kUnidadDistancia));
}

void deshacer_reescalado(Trajectory& x, Trajectory& v, Trajectory& a, MassArray& masa) {
    const double factorVel = 1.0 / (pow(kUnidadDistancia, 1.5) / (kUnidadDistancia * sqrt(kG * kUnidadMasa)));
    const double factorAcel = kG * kUnidadMasa / (kUnidadDistancia * kUnidadDistancia);

    for (auto& paso : x) {
        for (auto& planeta : paso) {
            for (auto& componente : planeta) {
                componente *= kUnidadDistancia;
            }
        }
    }

    for (auto& paso : v) {
        for (auto& planeta : paso) {
            for (auto& componente : planeta) {
                componente *= factorVel;
            }
        }
    }

    for (auto& paso : a) {
        for (auto& planeta : paso) {
            for (auto& componente : planeta) {
                componente *= factorAcel;
            }
        }
    }

    for (auto& masa_planeta : masa) {
        masa_planeta *= kUnidadMasa;
    }
}

PlanetArray calcular_aceleraciones(const PlanetArray& posiciones, const MassArray& masa) {
    PlanetArray aceleraciones{};
    
    // Masa del sol en unidades escaladas
    constexpr double masa_sol = 1.989e30 / 2.0e30;

    for (size_t i = 0; i < kNumPlanetas; ++i) {
        aceleraciones[i] = {0.0, 0.0};
        
        // Contribución del sol (ubicado en el origen)
        const double dx_sol = posiciones[i][0];
        const double dy_sol = posiciones[i][1];
        const double dist_sol_cubica = pow(dx_sol * dx_sol + dy_sol * dy_sol, 1.5);
        
        aceleraciones[i][0] += -masa_sol * dx_sol / dist_sol_cubica;
        aceleraciones[i][1] += -masa_sol * dy_sol / dist_sol_cubica;
        
        // Contribución de los otros planetas
        for (size_t k = 0; k < kNumPlanetas; ++k) {
            if (k == i) {
                continue;
            }

            const double dx = posiciones[i][0] - posiciones[k][0];
            const double dy = posiciones[i][1] - posiciones[k][1];
            const double distanciaCubica = pow(dx * dx + dy * dy, 1.5);

            aceleraciones[i][0] += -masa[k] * dx / distanciaCubica;
            aceleraciones[i][1] += -masa[k] * dy / distanciaCubica;
        }
    }

    return aceleraciones;
}

void Verlet(double& t, double h, int N, const PlanetArray& x0, const PlanetArray& v0, Trajectory& x, Trajectory& v, Trajectory& a, const MassArray& masa) {
    t = 0.0;
    x[0] = x0;
    v[0] = v0;
    a[0] = calcular_aceleraciones(x[0], masa);

    for (int n = 0; n + 1 < N; ++n) {
        for (size_t i = 0; i < kNumPlanetas; ++i) {
            for (size_t j = 0; j < kNumCoordenadas; ++j) {
                x[n + 1][i][j] = x[n][i][j] + v[n][i][j] * h + 0.5 * a[n][i][j] * h * h;
            }
        }

        a[n + 1] = calcular_aceleraciones(x[n + 1], masa);

        for (size_t i = 0; i < kNumPlanetas; ++i) {
            for (size_t j = 0; j < kNumCoordenadas; ++j) {
                v[n + 1][i][j] = v[n][i][j] + 0.5 * (a[n][i][j] + a[n + 1][i][j]) * h;
            }
        }

        t += h;
    }
}

void escribir_datos(ofstream& out, const Trajectory& x) {
    for (const auto& paso : x) {
        for (const auto& planeta : paso) {
            out << planeta[0] << ", " << planeta[1] << '\n';
        }
        out << '\n';
    }
}

void escribir_datos_energia(ofstream& out, const EnergyArray& data) {
    for (const auto& paso : data) {
        for (size_t i = 0; i < kNumPlanetas; ++i) {
            out << paso[i] << '\n';
        }
        out << '\n';
    }
}

void escribir_datos_periodo(ofstream& out, const array<double, kNumPlanetas>& periodos) {
    for (size_t i = 0; i < kNumPlanetas; ++i) {
        out << periodos[i] << '\n';
    }
}

void invariantes(const Trajectory& x, const Trajectory& v, const Trajectory& a, const MassArray& masa, EnergyArray& E, Trajectory& L, Trajectory& p, EnergyArray& mod_p) {
    for (size_t n = 0; n < x.size(); ++n) {
        for (size_t i = 0; i < kNumPlanetas; ++i) {
            const double energia_cinetica = 0.5 * masa[i] * (v[n][i][0] * v[n][i][0] + v[n][i][1] * v[n][i][1]);
            // ✅ Inicializamos la energía potencial con la del Sol
            double masa_sol = 2e30; // Masa del sol
            double dist_sol = sqrt(x[n][i][0] * x[n][i][0] + x[n][i][1] * x[n][i][1]);
            double energia_potencial = -kG * masa_sol * masa[i] / dist_sol;

            for (size_t k = 0; k < kNumPlanetas; ++k) {
                if (k == i) {
                    continue;
                }

                const double dx = x[n][i][0] - x[n][k][0];
                const double dy = x[n][i][1] - x[n][k][1];
                const double distancia = sqrt(dx * dx + dy * dy);
                energia_potencial += -kG * masa[i] * masa[k] / distancia;
            }

            E[n][i] = energia_cinetica + energia_potencial;
            L[n][i][0] = masa[i] * (x[n][i][1] * v[n][i][0] - x[n][i][0] * v[n][i][1]);
            L[n][i][1] = masa[i] * (x[n][i][0] * v[n][i][1] - x[n][i][1] * v[n][i][0]);
            p[n][i][0] = masa[i] * v[n][i][0];
            p[n][i][1] = masa[i] * v[n][i][1];
            mod_p[n][i] = sqrt(p[n][i][0] * p[n][i][0] + p[n][i][1] * p[n][i][1]);
        }
    }
}

void periodos(const EnergyArray& E, const MassArray& masa, array<double, kNumPlanetas>& periodos) {
    array<double, kNumPlanetas> energia_media{};

    for (const auto& paso : E) {
        for (size_t i = 0; i < kNumPlanetas; ++i) {
            energia_media[i] += paso[i];
        }
    }

    for (size_t i = 0; i < kNumPlanetas; ++i) {
        energia_media[i] /= static_cast<double>(E.size());
        const double semieje_mayor = -kG * kUnidadMasa * masa[i] / (2.0 * energia_media[i]);
        periodos[i] = 2.0 * M_PI * pow(semieje_mayor, 1.5) / sqrt(kG * kUnidadMasa);
    }
}

void convertir_periodo_a_dias(array<double, kNumPlanetas>& periodos) {
    constexpr double segundos_por_dia = 86400.0;
    for (auto& periodo : periodos) {
        periodo /= segundos_por_dia;
    }
}