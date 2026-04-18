import matplotlib.pyplot as plt

# Read the data from energia.dat
data = []
with open('energia.dat', 'r') as f:
    for line in f:
        line = line.strip()
        if line:  # Skip empty lines
            data.append(float(line))

# Assuming 9 energies per iteration (8 planets + Pluto)
energies_per_iteration = 9
num_iterations = len(data) // energies_per_iteration

# List of planet names
planets = ['Mercurio', 'Venus', 'Tierra', 'Marte', 'Júpiter', 'Saturno', 'Urano', 'Neptuno', 'Plutón']

# Initialize lists for each planet's energies
energies = [[] for _ in range(energies_per_iteration)]

# Fill the energies lists
for i in range(num_iterations):
    for j in range(energies_per_iteration):
        energies[j].append(data[i * energies_per_iteration + j])

# Plot each planet's energy
iterations = list(range(1, num_iterations + 1))
for j in range(energies_per_iteration):
    plt.plot(iterations, energies[j], label=planets[j])

plt.xlabel('Número de iteraciones')
plt.ylabel('Energía (J)')
plt.title('Energía de cada planeta vs Iteraciones')
plt.legend()
plt.grid(True)
plt.savefig('energia_planetas.png', dpi=300, bbox_inches='tight')
plt.show()