#include <iostream>
#include <fstream>
#include <random>
#include <cmath>
#include <iomanip>
#include <string>

int main(int argc, char* argv[]) {
    // Defaults
    int N = 500;
    double baseMass = 1.989e30; // Masa base de los sistemas solares (masa del sistema solar)
    std::string outFile = "condiciones_iniciales.txt";
    double minR = 1.8e13;    
    double maxR = 1e21;  //diámetro de la Vía Láctea.  
    double centralMass = 8.2e36; // masa del agujero negro supermasivo del centro de nuestra galaxia.
    bool equalMass = false; // if true, all masses exactly baseMass

    // Simple argument parsing: allow passing N and baseMass and output filename
    if (argc > 1) N = std::stoi(argv[1]);
    if (argc > 2) baseMass = std::stod(argv[2]);
    if (argc > 3) outFile = argv[3];
    if (argc > 4) minR = std::stod(argv[4]);
    if (argc > 5) maxR = std::stod(argv[5]);
    if (argc > 6) centralMass = std::stod(argv[6]);
    if (argc > 7) equalMass = (std::string(argv[7]) == "1");

    std::ofstream ofs(outFile);
    if (!ofs) {
        std::cerr << "Error: no se pudo abrir " << outFile << " para escritura\n";
        return 1;
    }

    // Archivo adicional con las distancias al centro (una por línea)
    std::ofstream ofs_dist("distancias.dat");
    if (!ofs_dist) {
        std::cerr << "Error: no se pudo abrir distancias.dat para escritura\n";
        return 1;
    }

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> var(-0.05, 0.05); // ±5% variation

    const double G = 6.67430e-11;

    ofs << std::scientific << std::setprecision(3);

    for (int i = 0; i < N; ++i) {
        double r;
        if (N == 1) {
            r = 0.5 * (minR + maxR);
        } else {
            r = minR + (maxR - minR) * i / static_cast<double>(N - 1);
        }

        double mass = baseMass;
        if (!equalMass) mass = baseMass * (1.0 + var(gen));

        double theta = (2.0 * M_PI * i) / N;
        double x = r * std::cos(theta);
        double y = r * std::sin(theta);

        // approximate circular velocity around centralMass, tangential to the radius vector
        double v = std::sqrt(G * centralMass / r);
        double vx = -v * std::sin(theta);
        double vy = v * std::cos(theta);

        // Small random perturbations in both velocity components
        double delta_vx = v * 0.05 * var(gen);
        double delta_vy = v * 0.05 * var(gen);
        vx += delta_vx;
        vy += delta_vy;

        ofs << mass << '\t' << x << '\t' << y << '\t' << vx << '\t' << vy << '\n';

        // Escribimos la distancia radial al archivo de distancias
        ofs_dist << r << '\n';
    }

    ofs.close();
    ofs_dist.close();
    std::cout << "Archivo generado: " << outFile << " (" << N << " objetos)\n";
    return 0;
}
