#include <iostream>
#include <fstream>
#include <random>
#include <cmath>
#include <iomanip>
#include <string>

int main(int argc, char* argv[]) {
    // Defaults
    int N = 1000;
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
    
    // Distribuciones para posiciones completamente aleatorias
    std::uniform_real_distribution<double> dist_r(minR, maxR);                    // Radio aleatorio
    std::uniform_real_distribution<double> dist_theta(0.0, 2.0 * M_PI);          // Ángulo aleatorio
    std::uniform_real_distribution<double> dist_mass(-0.05, 0.05);               // Variación de masa
    
    // Distribuciones para velocidades aleatorias
    // La velocidad máxima se basa en la velocidad orbital a la distancia más cercana
    const double G = 6.67430e-11;
    double v_max = std::sqrt(G * centralMass / minR);
    std::uniform_real_distribution<double> dist_vx(-v_max, v_max);               // Velocidad x aleatoria
    std::uniform_real_distribution<double> dist_vy(-v_max, v_max);               // Velocidad y aleatoria

    ofs << std::scientific << std::setprecision(3);

    for (int i = 0; i < N; ++i) {
        // Generar posiciones completamente aleatorias
        double r = dist_r(gen);
        double theta = dist_theta(gen);
        double x = r * std::cos(theta);
        double y = r * std::sin(theta);

        // Generar masa con variación opcional
        double mass = baseMass;
        if (!equalMass) mass = baseMass * (1.0 + dist_mass(gen));

        // Generar velocidades completamente aleatorias
        double vx = dist_vx(gen);
        double vy = dist_vy(gen);

        ofs << mass << '\t' << x << '\t' << y << '\t' << vx << '\t' << vy << '\n';

        // Escribimos la distancia radial al archivo de distancias
        ofs_dist << r << '\n';
    }

    ofs.close();
    ofs_dist.close();
    std::cout << "Archivo generado: " << outFile << " (" << N << " objetos)\n";
    return 0;
}
