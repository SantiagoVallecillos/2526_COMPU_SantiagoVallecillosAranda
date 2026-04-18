import matplotlib.pyplot as plt

# Read the data from momento_angular.dat
data = []
with open('momento_angular.dat', 'r') as f:
    for line in f:
        line = line.strip()
        if line:
            parts = line.split(',')
            data.append((float(parts[0]), float(parts[1])))

# Assuming 9 momenta per iteration (8 planets + Pluto)
momenta_per_iteration = 9
num_iterations = len(data) // momenta_per_iteration

# List of planet names
planets = ['Mercurio', 'Venus', 'Tierra', 'Marte', 'Júpiter', 'Saturno', 'Urano', 'Neptuno', 'Plutón']

# Initialize lists for each planet's angular momentum components
angular_momenta_x = [[] for _ in range(momenta_per_iteration)]
angular_momenta_y = [[] for _ in range(momenta_per_iteration)]

# Fill the lists
for i in range(num_iterations):
    for j in range(momenta_per_iteration):
        angular_momenta_x[j].append(data[i * momenta_per_iteration + j][0])
        angular_momenta_y[j].append(data[i * momenta_per_iteration + j][1])

# Plot for x component
iterations = list(range(1, num_iterations + 1))
plt.figure(figsize=(10, 6))
for j in range(momenta_per_iteration):
    plt.plot(iterations, angular_momenta_x[j], label=planets[j])
plt.xlabel('Número de iteraciones')
plt.ylabel('Componente x del momento angular')
plt.title('Componente x del momento angular de cada planeta vs Iteraciones')
plt.legend()
plt.grid(True)
plt.savefig('momento_angular_x.png', dpi=300, bbox_inches='tight')
plt.close()

# Plot for y component
plt.figure(figsize=(10, 6))
for j in range(momenta_per_iteration):
    plt.plot(iterations, angular_momenta_y[j], label=planets[j])
plt.xlabel('Número de iteraciones')
plt.ylabel('Componente y del momento angular')
plt.title('Componente y del momento angular de cada planeta vs Iteraciones')
plt.legend()
plt.grid(True)
plt.savefig('momento_angular_y.png', dpi=300, bbox_inches='tight')
plt.close()