#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

constexpr std::size_t kNumPlanetas = 4;
constexpr std::size_t kNumCoordenadas = 2;
constexpr double kG = 6.67430e-11;

using Vector2D = std::array<double, kNumCoordenadas>;
using PlanetArray = std::array<Vector2D, kNumPlanetas>;
using Trajectory = std::vector<PlanetArray>;
using MassArray = std::array<double, kNumPlanetas>;

void leer_datos(std::ifstream& data, PlanetArray& x0, PlanetArray& v0, MassArray& masa);
void reescalar(double& h, PlanetArray& x0, PlanetArray& v0, MassArray& masa);
void deshacer_reescalado(Trajectory& x, Trajectory& v, Trajectory& a, MassArray& masa);
void escribir_datos_pos(std::ofstream& out, const Trajectory& x);
PlanetArray calcular_aceleraciones(const PlanetArray& posiciones, const MassArray& masa);

int main() {
    std::ifstream data("condiciones_iniciales.txt");
    if (!data) {
        std::cerr << "Error: no se pudo abrir condiciones_iniciales.txt\n";
        return 1;
    }

    std::ofstream trayectorias("posiciones_planetas.dat");
    if (!trayectorias) {
        std::cerr << "Error: no se pudo crear posiciones_planetas.dat\n";
        return 1;
    }

    constexpr int N = 1000;
    double t = 0.0;
    double h = 0.01;

    PlanetArray x0{};
    PlanetArray v0{};
    MassArray masa{};

    leer_datos(data, x0, v0, masa);
    reescalar(h, x0, v0, masa);

    Trajectory x(N, PlanetArray{});
    Trajectory v(N, PlanetArray{});
    Trajectory a(N, PlanetArray{});

    x[0] = x0;
    v[0] = v0;
    a[0] = calcular_aceleraciones(x[0], masa);

    for (int n = 0; n + 1 < N; ++n) {
        for (std::size_t i = 0; i < kNumPlanetas; ++i) {
            for (std::size_t j = 0; j < kNumCoordenadas; ++j) {
                x[n + 1][i][j] = x[n][i][j] + v[n][i][j] * h + 0.5 * a[n][i][j] * h * h;
            }
        }

        a[n + 1] = calcular_aceleraciones(x[n + 1], masa);

        for (std::size_t i = 0; i < kNumPlanetas; ++i) {
            for (std::size_t j = 0; j < kNumCoordenadas; ++j) {
                const double w = v[n][i][j] + 0.5 * a[n][i][j] * h;
                v[n + 1][i][j] = w + 0.5 * a[n + 1][i][j] * h;
            }
        }

        t += h;
    }

    deshacer_reescalado(x, v, a, masa);
    escribir_datos_pos(trayectorias, x);

    return 0;
}

void leer_datos(std::ifstream& data, PlanetArray& x0, PlanetArray& v0, MassArray& masa) {
    for (std::size_t i = 0; i < kNumPlanetas; ++i) {
        data >> masa[i] >> x0[i][0] >> x0[i][1] >> v0[i][0] >> v0[i][1];
    }
}

void escribir_datos_pos(std::ofstream& out, const Trajectory& x) {
    for (std::size_t i = 0; i < x.size(); ++i) {
        for (std::size_t planeta = 0; planeta < kNumPlanetas; ++planeta) {
            out << x[i][planeta][0] << ", " << x[i][planeta][1] << '\n';
        };
        out << '\n';
    }
}

void reescalar(double& h, PlanetArray& x0, PlanetArray& v0, MassArray& masa) {
    constexpr double c = 1.496e11;
    constexpr double M = 2.0e30;

    for (std::size_t i = 0; i < kNumPlanetas; ++i) {
        x0[i][0] /= c;
        x0[i][1] /= c;
        v0[i][0] *= std::pow(c, 1.5) / (c * std::sqrt(kG * M));
        v0[i][1] *= std::pow(c, 1.5) / (c * std::sqrt(kG * M));
        masa[i] /= M;
    }

    h *= std::sqrt(kG * M / (c * c * c));
}

void deshacer_reescalado(Trajectory& x, Trajectory& v, Trajectory& a, MassArray& masa) {
    constexpr double c = 1.496e11;
    constexpr double M = 2.0e30;

    for (auto& paso : x) {
        for (auto& planeta : paso) {
            for (auto& componente : planeta) {
                componente *= c;
            }
        }
    }

    const double factorVel = c / std::sqrt(kG * M / (c * c * c));
    const double factorAcel = kG * M / (c * c);

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

    for (double& masa_planeta : masa) {
        masa_planeta *= M;
    }
}

PlanetArray calcular_aceleraciones(const PlanetArray& posiciones, const MassArray& masa) {
    PlanetArray aceleraciones{};

    for (std::size_t i = 0; i < kNumPlanetas; ++i) {
        aceleraciones[i] = {0.0, 0.0};
        for (std::size_t k = 0; k < kNumPlanetas; ++k) {
            if (k == i) {
                continue;
            }
            const double dx = posiciones[i][0] - posiciones[k][0];
            const double dy = posiciones[i][1] - posiciones[k][1];
            const double dist3 = std::pow(dx * dx + dy * dy, 1.5);
            aceleraciones[i][0] += -kG * masa[k] * dx / dist3;
            aceleraciones[i][1] += -kG * masa[k] * dy / dist3;
        }
    }

    return aceleraciones;
}
