import matplotlib.pyplot as plt

# Read the data from momento_lineal.dat
data = []
with open('momento_lineal.dat', 'r') as f:
    for line in f:
        line = line.strip()
        if line:  # Skip empty lines
            data.append(float(line))

# Assuming 9 momenta per iteration (8 planets + Pluto)
momenta_per_iteration = 9
num_iterations = len(data) // momenta_per_iteration

# List of planet names
planets = ['Mercurio', 'Venus', 'Tierra', 'Marte', 'Júpiter', 'Saturno', 'Urano', 'Neptuno', 'Plutón']

# Initialize lists for each planet's momenta
momenta = [[] for _ in range(momenta_per_iteration)]

# Fill the momenta lists
for i in range(num_iterations):
    for j in range(momenta_per_iteration):
        momenta[j].append(data[i * momenta_per_iteration + j])

# Plot each planet's linear momentum
iterations = list(range(1, num_iterations + 1))
for j in range(momenta_per_iteration):
    plt.plot(iterations, momenta[j], label=planets[j])

plt.xlabel('Número de iteraciones')
plt.ylabel('Módulo del momento lineal')
plt.title('Módulo del momento lineal de cada planeta vs Iteraciones')
plt.legend()
plt.grid(True)
plt.savefig('momento_lineal_planetas.png', dpi=300, bbox_inches='tight')
plt.show()