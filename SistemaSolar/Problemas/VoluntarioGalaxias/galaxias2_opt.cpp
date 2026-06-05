#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include <chrono>
#include <random>
#include <omp.h>

using namespace std;

constexpr size_t kNumPlanetas = 500;
constexpr size_t kNumCoordenadas = 2;
constexpr double kG = 6.67430e-11;
constexpr double kUnidadDistancia = 1e21; // Radio vía láctea
constexpr double kUnidadMasa = 8.2e36;
constexpr double kRadioColision = 4.5e9; // Radio planetario del Sistema Solar.
constexpr double kRadioRegeneracion = 0.3;
constexpr double kRadioHorizonte = 0.1;
constexpr double kPi = 3.14159265358979323846;
constexpr double kRadioColisionEscalada = kRadioColision / kUnidadDistancia;

using Vector2D = array<double, kNumCoordenadas>;
using PlanetArray = array<Vector2D, kNumPlanetas>;
using Trajectory = vector<PlanetArray>;
using EnergyArray = vector<array<double, kNumPlanetas>>;
using MassArray = array<double, kNumPlanetas>;

void leer_datos(ifstream& data, PlanetArray& x0, PlanetArray& v0, MassArray& masa);
void reescalar(double& h, PlanetArray& x0, PlanetArray& v0, MassArray& masa);
void deshacer_reescalado(Trajectory& x, Trajectory& v, Trajectory& a, MassArray& masa);
bool esta_cerca_origen(const Vector2D& posicion);
void generar_en_orbita(Vector2D& x, Vector2D& v, std::mt19937_64& rng);
void regenerar_si_cerca_origen(PlanetArray& x, PlanetArray& v, std::mt19937_64& rng);
void calcular_siguiente_paso(const PlanetArray& x_curr, const PlanetArray& v_curr, const PlanetArray& a_curr, PlanetArray& x_next, PlanetArray& v_next, PlanetArray& a_next, const MassArray& masa, double h);
void Verlet(double& t, double h, int N, const PlanetArray& x0, const PlanetArray& v0, Trajectory& x, Trajectory& v, Trajectory& a, const MassArray& masa, std::mt19937_64& rng, vector<double>& tiempos, vector<double>& energias);
PlanetArray calcular_aceleraciones(const PlanetArray& posiciones, const MassArray& masa);
double energia_total(const PlanetArray& x, const PlanetArray& v, const MassArray& masa);
void escribir_datos(ofstream& out, const Trajectory& x);
void escribir_datos_energia(ofstream& out, const EnergyArray& data);
void escribir_datos_periodo(ofstream& out, const array<double, kNumPlanetas>& periodos);
void invariantes(const Trajectory& x, const Trajectory& v, const MassArray& masa, EnergyArray& E, Trajectory& L, Trajectory& p, EnergyArray& mod_p);
void periodos(const EnergyArray& E, const MassArray& masa, array<double, kNumPlanetas>& periodos);
void convertir_periodo_a_dias(array<double, kNumPlanetas>& periodos);
double choque_elastico(double v1, double v2, double m1, double m2);
bool haycolision(const PlanetArray& posiciones, size_t i, size_t k);
void calcular_distribucion_radial(const Trajectory& x, const MassArray& masa, const string& nombre_archivo);
void calcular_momento_inercia_medio(const Trajectory& x, const MassArray& masa, const string& nombre_archivo);
void calcular_flujo_masa_absorbido(const Trajectory& x, const MassArray& masa, double h_reducido);

int main() {
    auto inicio = chrono::high_resolution_clock::now();

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

    ofstream estado_estacionario("estado_estacionario.dat");
    if (!estado_estacionario) {
        cerr << "Error: no se pudo crear estado_estacionario.dat\n";
        return 1;
    }

    constexpr int N_final = 10000;
    double t = 0.0;
    double h = 3.156e13; // 1 millón de años en segundos

    PlanetArray x0{};
    PlanetArray v0{};
    MassArray masa{};

    leer_datos(data, x0, v0, masa);
    reescalar(h, x0, v0, masa);

    std::mt19937_64 rng(static_cast<unsigned long>(chrono::high_resolution_clock::now().time_since_epoch().count()));

    PlanetArray x_curr = x0;
    PlanetArray v_curr = v0;
    regenerar_si_cerca_origen(x_curr, v_curr, rng);
    PlanetArray a_curr = calcular_aceleraciones(x_curr, masa);
    PlanetArray x_next{};
    PlanetArray v_next{};
    PlanetArray a_next{};

    double energia_anterior = energia_total(x_curr, v_curr, masa);
    constexpr double kEnergiaUmbral = 1e-5;
    constexpr int kPasosEstablesRequeridos = 10;
    int pasos_estables = 0;
    int iteracion = 0;
    while (pasos_estables < kPasosEstablesRequeridos) {
        regenerar_si_cerca_origen(x_curr, v_curr, rng);
        calcular_siguiente_paso(x_curr, v_curr, a_curr, x_next, v_next, a_next, masa, h);

        double energia_actual = energia_total(x_next, v_next, masa);
        if (fabs(energia_actual - energia_anterior) <= kEnergiaUmbral) {
            pasos_estables++;
        } else {
            pasos_estables = 0;
        }

        energia_anterior = energia_actual;
        x_curr = x_next;
        v_curr = v_next;
        a_curr = a_next;
        t += h;
        iteracion++;
    }

    Trajectory x(N_final, PlanetArray{});
    Trajectory v(N_final, PlanetArray{});
    Trajectory a(N_final, PlanetArray{});
    EnergyArray E(N_final);
    Trajectory L(N_final, PlanetArray{});
    Trajectory p(N_final, PlanetArray{});
    EnergyArray mod_p(N_final);
    array<double, kNumPlanetas> periodo{};
    vector<double> tiempos;
    vector<double> energias;

    Verlet(t, h, N_final, x_curr, v_curr, x, v, a, masa, rng, tiempos, energias);
    escribir_datos(trayectorias, x);
    escribir_datos(velocidades, v);
    escribir_datos(aceleraciones, a);

    deshacer_reescalado(x, v, a, masa);
    invariantes(x, v, masa, E, L, p, mod_p);
    periodos(E, masa, periodo);
    convertir_periodo_a_dias(periodo);
    calcular_distribucion_radial(x, masa, "densidad_radial.dat");
    calcular_distribucion_radial(x, masa, "densidad_radial.dat");
    calcular_momento_inercia_medio(x, masa, "momento_inercia.dat");
    calcular_flujo_masa_absorbido(x, masa, h); // Le pasamos h para calcular el tiempo total

    const double factorEnergia = kG * kUnidadMasa * kUnidadMasa / kUnidadDistancia;
    const double factorTiempo = pow(kUnidadDistancia, 1.5) / sqrt(kG * kUnidadMasa);
    for (size_t j = 0; j < tiempos.size(); ++j) {
        double tiempo_real = tiempos[j] * factorTiempo;
        double energia_real = energias[j] * factorEnergia;
        estado_estacionario << tiempo_real << " " << energia_real << "\n";
    }

    auto fin = chrono::high_resolution_clock::now();
    chrono::duration<double, milli> tiempo_ejecucion = fin - inicio;
    cout << "El código tardó: " << tiempo_ejecucion.count() << " milisegundos." << endl;
    cout << "Iteraciones hasta estabilidad: " << iteracion << endl;

    return 0;
}

void leer_datos(ifstream& data, PlanetArray& x0, PlanetArray& v0, MassArray& masa) {
    for (size_t i = 0; i < kNumPlanetas; ++i) {
        data >> masa[i] >> x0[i][0] >> x0[i][1] >> v0[i][0] >> v0[i][1];
    }
}

void reescalar(double& h, PlanetArray& x0, PlanetArray& v0, MassArray& masa) {
#pragma omp parallel for
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

bool esta_cerca_origen(const Vector2D& posicion) {
    return posicion[0] * posicion[0] + posicion[1] * posicion[1] < kRadioRegeneracion * kRadioRegeneracion;
}

void generar_en_orbita(Vector2D& x, Vector2D& v, std::mt19937_64& rng) {
    std::uniform_real_distribution<double> dist_radio(kRadioHorizonte, 1.0);
    std::uniform_real_distribution<double> dist_angulo(0.0, 2.0 * kPi);

    const double radio = dist_radio(rng);
    const double theta = dist_angulo(rng);
    x[0] = radio * cos(theta);
    x[1] = radio * sin(theta);

    const double velocidad_orbital = sqrt(1.0 / radio);
    v[0] = -velocidad_orbital * sin(theta);
    v[1] = velocidad_orbital * cos(theta);
}

void regenerar_si_cerca_origen(PlanetArray& x, PlanetArray& v, std::mt19937_64& rng) {
    for (size_t i = 0; i < kNumPlanetas; ++i) {
        if (esta_cerca_origen(x[i])) {
            generar_en_orbita(x[i], v[i], rng);
        }
    }
}

void calcular_siguiente_paso(const PlanetArray& x_curr, const PlanetArray& v_curr, const PlanetArray& a_curr, PlanetArray& x_next, PlanetArray& v_next, PlanetArray& a_next, const MassArray& masa, double h) {
#pragma omp parallel for
    for (size_t i = 0; i < kNumPlanetas; ++i) {
        for (size_t j = 0; j < kNumCoordenadas; ++j) {
            x_next[i][j] = x_curr[i][j] + v_curr[i][j] * h + 0.5 * a_curr[i][j] * h * h;
        }
    }

    a_next = calcular_aceleraciones(x_next, masa);

#pragma omp parallel for
    for (size_t i = 0; i < kNumPlanetas; ++i) {
        for (size_t j = 0; j < kNumCoordenadas; ++j) {
            v_next[i][j] = v_curr[i][j] + 0.5 * (a_curr[i][j] + a_next[i][j]) * h;
        }
    }

    for (size_t i = 0; i < kNumPlanetas; ++i) {
        for (size_t k = i + 1; k < kNumPlanetas; ++k) {
            if (haycolision(x_next, i, k)) {
                double v1_final_x = choque_elastico(v_curr[i][0], v_curr[k][0], masa[i], masa[k]);
                double v1_final_y = choque_elastico(v_curr[i][1], v_curr[k][1], masa[i], masa[k]);
                double v2_final_x = choque_elastico(v_curr[k][0], v_curr[i][0], masa[k], masa[i]);
                double v2_final_y = choque_elastico(v_curr[k][1], v_curr[i][1], masa[k], masa[i]);

                v_next[i][0] = v1_final_x;
                v_next[i][1] = v1_final_y;
                v_next[k][0] = v2_final_x;
                v_next[k][1] = v2_final_y;
            }
        }
    }
}

PlanetArray calcular_aceleraciones(const PlanetArray& posiciones, const MassArray& masa) {
    PlanetArray aceleraciones{};

#pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < kNumPlanetas; ++i) {
        aceleraciones[i] = {0.0, 0.0};

        const double dx_sol = posiciones[i][0];
        const double dy_sol = posiciones[i][1];
        const double dist_sol2 = dx_sol * dx_sol + dy_sol * dy_sol;
        const double dist_sol_cubica = pow(dist_sol2, 1.5);

        if (dist_sol_cubica > 0.0) {
            aceleraciones[i][0] += -dx_sol / dist_sol_cubica;
            aceleraciones[i][1] += -dy_sol / dist_sol_cubica;
        }

        for (size_t k = 0; k < kNumPlanetas; ++k) {
            if (k == i) {
                continue;
            }

            const double dx = posiciones[i][0] - posiciones[k][0];
            const double dy = posiciones[i][1] - posiciones[k][1];
            const double dist2 = dx * dx + dy * dy;

            if (dist2 <= 0.0) {
                continue;
            }

            const double distanciaCubica = pow(dist2, 1.5);
            aceleraciones[i][0] += -masa[k] * dx / distanciaCubica;
            aceleraciones[i][1] += -masa[k] * dy / distanciaCubica;
        }
    }

    return aceleraciones;
}

void Verlet(double& t, double h, int N, const PlanetArray& x0, const PlanetArray& v0, Trajectory& x, Trajectory& v, Trajectory& a, const MassArray& masa, std::mt19937_64& rng, vector<double>& tiempos, vector<double>& energias) {
    PlanetArray x_curr = x0;
    PlanetArray v_curr = v0;
    regenerar_si_cerca_origen(x_curr, v_curr, rng);
    PlanetArray a_curr = calcular_aceleraciones(x_curr, masa);
    PlanetArray x_next{};
    PlanetArray v_next{};
    PlanetArray a_next{};

    t = 0.0;
    if (N > 0) {
        x[0] = x_curr;
        v[0] = v_curr;
        a[0] = a_curr;
        tiempos.push_back(t);
        energias.push_back(energia_total(x_curr, v_curr, masa));
    }

    for (int n = 0; n + 1 < N; ++n) {
        calcular_siguiente_paso(x_curr, v_curr, a_curr, x_next, v_next, a_next, masa, h);

        x[n + 1] = x_next;
        v[n + 1] = v_next;
        a[n + 1] = a_next;

        x_curr = x_next;
        v_curr = v_next;
        a_curr = a_next;
        t += h;
        
        tiempos.push_back(t);
        energias.push_back(energia_total(x_curr, v_curr, masa));
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

void invariantes(const Trajectory& x, const Trajectory& v, const MassArray& masa, EnergyArray& E, Trajectory& L, Trajectory& p, EnergyArray& mod_p) {
    size_t num_pasos = x.size();
#pragma omp parallel for schedule(dynamic)
    for (size_t n = 0; n < num_pasos; ++n) {
        for (size_t i = 0; i < kNumPlanetas; ++i) {
            const double energia_cinetica = 0.5 * masa[i] * (v[n][i][0] * v[n][i][0] + v[n][i][1] * v[n][i][1]);
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

#pragma omp parallel for
    for (size_t i = 0; i < kNumPlanetas; ++i) {
        energia_media[i] /= static_cast<double>(E.size());
        if (energia_media[i] >= 0.0) {
            periodos[i] = 0.0;
            continue;
        }
        const double semieje_mayor = -kG * kUnidadMasa * masa[i] / (2.0 * energia_media[i]);
        periodos[i] = 2.0 * kPi * pow(semieje_mayor, 1.5) / sqrt(kG * kUnidadMasa);
    }
}

void convertir_periodo_a_dias(array<double, kNumPlanetas>& periodos) {
    constexpr double segundos_por_dia = 86400.0;
#pragma omp parallel for
    for (size_t i = 0; i < kNumPlanetas; ++i) {
        periodos[i] /= segundos_por_dia;
    }
}

double energia_total(const PlanetArray& x, const PlanetArray& v, const MassArray& masa) {
    double energia = 0.0;

#pragma omp parallel for reduction(+:energia)
    for (size_t i = 0; i < kNumPlanetas; ++i) {
        energia += 0.5 * masa[i] * (v[i][0] * v[i][0] + v[i][1] * v[i][1]);
        double dist_sol = sqrt(x[i][0] * x[i][0] + x[i][1] * x[i][1]);
        if (dist_sol > 0.0) {
            energia += -masa[i] / dist_sol;
        }

        for (size_t k = i + 1; k < kNumPlanetas; ++k) {
            double dx = x[i][0] - x[k][0];
            double dy = x[i][1] - x[k][1];
            double distancia = sqrt(dx * dx + dy * dy);
            if (distancia > 0.0) {
                energia += -masa[i] * masa[k] / distancia;
            }
        }
    }

    return energia;
}

double choque_elastico(double v1, double v2, double m1, double m2){
    double v1_final = (v1 * (m1 - m2) + 2 * m2 * v2) / (m1 + m2);
    return v1_final;
}

bool haycolision(const PlanetArray& posiciones, size_t i, size_t k) {
    double dx = posiciones[i][0] - posiciones[k][0];
    double dy = posiciones[i][1] - posiciones[k][1];
    double distancia = sqrt(dx * dx + dy * dy);
    return distancia < kRadioColisionEscalada;
}

void calcular_distribucion_radial(const Trajectory& x, const MassArray& masa, const string& nombre_archivo) {
    // Definimos el número de anillos concéntricos (puedes ajustar este número)
    constexpr size_t kNumAnillos = 100; 
    
    // El radio máximo de la galaxia es el radio de la Vía Láctea en metros
    const double r_max = kUnidadDistancia; 
    const double delta_r = r_max / kNumAnillos;

    // Vector para acumular la masa en cada anillo
    vector<double> masa_acumulada(kNumAnillos, 0.0);
    size_t num_pasos = x.size();

    // 1. Recorrer todos los pasos temporales y todos los planetas
    for (size_t n = 0; n < num_pasos; ++n) {
        for (size_t i = 0; i < kNumPlanetas; ++i) {
            // Calcular la distancia del sistema solar al centro (agujero negro)
            double r = sqrt(x[n][i][0] * x[n][i][0] + x[n][i][1] * x[n][i][1]);
            
            // Determinar en qué anillo cae según su radio
            size_t indice_anillo = static_cast<size_t>(r / delta_r);
            
            // Asegurarnos de no salirnos del rango del vector por errores numéricos
            if (indice_anillo < kNumAnillos) {
                masa_acumulada[indice_anillo] += masa[i];
            }
        }
    }

    // 2. Abrir el archivo de salida para guardar los resultados
    ofstream out(nombre_archivo);
    if (!out) {
        cerr << "Error: no se pudo crear el archivo " << nombre_archivo << "\n";
        return;
    }

    // 3. Calcular la densidad media para cada anillo y escribirla
    for (size_t i = 0; i < kNumAnillos; ++i) {
        double r_interno = i * delta_r;
        double r_externo = (i + 1) * delta_r;
        
        // Área del anillo circular en 2D: A = pi * (r_out^2 - r_in^2)
        double area_anillo = kPi * (r_externo * r_externo - r_interno * r_interno);
        
        // Masa promedio que hubo en este anillo por cada instante de tiempo
        double masa_promedio = masa_acumulada[i] / static_cast<double>(num_pasos);
        
        // Densidad superficial media de masa (kg / m^2)
        double densidad = masa_promedio / area_anillo;
        
        // Radio medio del anillo (punto medio útil para representar en gráficas)
        double r_medio = (r_interno + r_externo) / 2.0;

        // Guardamos: [Radio medio] [Densidad de masa]
        out << r_medio << " " << densidad << "\n";
    }
    
    out.close();
}

void calcular_momento_inercia_medio(const Trajectory& x, const MassArray& masa, const string& nombre_archivo) {
    double inercia_total_acumulada = 0.0;
    size_t num_pasos = x.size();

    // 1. Calcular el momento de inercia en cada instante de tiempo
    for (size_t n = 0; n < num_pasos; ++n) {
        double inercia_paso = 0.0;
        for (size_t i = 0; i < kNumPlanetas; ++i) {
            // Distancia al cuadrado al centro (agujero negro)
            double r2 = x[n][i][0] * x[n][i][0] + x[n][i][1] * x[n][i][1];
            inercia_paso += masa[i] * r2;
        }
        inercia_total_acumulada += inercia_paso;
    }

    // 2. Calcular la media temporal
    double inercia_media = inercia_total_acumulada / static_cast<double>(num_pasos);

    // 3. Imprimir por consola y guardar en archivo
    cout << "Momento de inercia medio del sistema: " << inercia_media << " kg*m^2\n";

    ofstream out(nombre_archivo);
    if (out) {
        out << "Momento_de_inercia_medio(kg*m^2)\n";
        out << inercia_media << "\n";
        out.close();
    } else {
        cerr << "Error: no se pudo crear el archivo " << nombre_archivo << "\n";
    }
}

void calcular_flujo_masa_absorbido(const Trajectory& x, const MassArray& masa, double h_reducido) {
    double masa_absorbida_total = 0.0;
    size_t num_pasos = x.size();
    
    // Como la trayectoria ya está des-reescalada (en metros), adaptamos el radio de regeneración
    double radio_absorcion_metros = kRadioRegeneracion * kUnidadDistancia;
    double r_abs_cuadrado = radio_absorcion_metros * radio_absorcion_metros;

    // 1. Detectar cuántas veces un cuerpo cruza el horizonte hacia adentro
    for (size_t i = 0; i < kNumPlanetas; ++i) {
        // Estado inicial del planeta i en el primer paso
        bool estaba_dentro = (x[0][i][0] * x[0][i][0] + x[0][i][1] * x[0][i][1]) < r_abs_cuadrado;
        
        for (size_t n = 1; n < num_pasos; ++n) {
            double r2 = x[n][i][0] * x[n][i][0] + x[n][i][1] * x[n][i][1];
            bool esta_dentro = r2 < r_abs_cuadrado;
            
            // Si en el paso anterior estaba fuera y ahora está dentro, ha sido absorbido
            if (!estaba_dentro && esta_dentro) {
                masa_absorbida_total += masa[i];
            }
            
            // Actualizamos el estado para la siguiente iteración temporal
            estaba_dentro = esta_dentro;
        }
    }

    // 2. Calcular el tiempo total real transcurrido durante las mediciones
    // Hay que deshacer el reescalado del tiempo (h_reducido) para tener segundos
    const double factorTiempo = pow(kUnidadDistancia, 1.5) / sqrt(kG * kUnidadMasa);
    double tiempo_total_real_segundos = static_cast<double>(num_pasos) * h_reducido * factorTiempo;

    // 3. Calcular el flujo (kg / s)
    double flujo_medio = 0.0;
    if (tiempo_total_real_segundos > 0.0) {
        flujo_medio = masa_absorbida_total / tiempo_total_real_segundos;
    }

    // 4. Mostrar resultados por consola
    cout << "Masa total absorbida durante medicion: " << masa_absorbida_total << " kg\n";
    cout << "Tiempo total de medicion: " << tiempo_total_real_segundos << " s\n";
    cout << "Flujo medio de masa absorbido: " << flujo_medio << " kg/s\n";
}
