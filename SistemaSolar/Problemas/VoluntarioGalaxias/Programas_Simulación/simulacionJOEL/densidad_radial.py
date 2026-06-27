import numpy as np
import matplotlib.pyplot as plt

# 1. Cargar los datos del archivo
# Columna 0: Radio medio (m)
# Columna 1: Densidad media (kg/m^2)
# Columna 2: Error estadístico de la densidad (kg/m^2)
datos = np.loadtxt("densidad_radial.dat")
r_medio = datos[:, 0]
densidad = datos[:, 1]
error = datos[:, 2]

# 2. Configurar la figura
plt.figure(figsize=(10, 6))

# 3. Dibujar la densidad con barras de error clásicas
plt.errorbar(r_medio, densidad, yerr=error, 
             fmt='-o',          # Formato: línea sólida con marcadores de punto
             color='#1f77b4',   # Color de la línea y el punto
             ecolor='#333333',  # Color de las barras de error (gris oscuro/negro)
             elinewidth=1.2,    # Grosor de la barra de error
             capsize=3,         # Tamaño del remate horizontal (la "T") de la barra
             markersize=4,      # Tamaño del punto de los datos
             label='Densidad media con error (1$\sigma$)')

# 4. Configurar ejes y etiquetas con las unidades correctas
plt.title("Distribución Radial de Densidad del Sistema", fontsize=14, pad=15)
plt.xlabel("Posición radial r (m)", fontsize=12)
plt.ylabel("Densidad de masa superficial $\\rho$ (kg/m²)", fontsize=12)

# Formatear ejes en notación científica (distancias del orden de 1e21)
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

# 5. Detalles estéticos
plt.grid(True, linestyle='--', alpha=0.6)
plt.legend(loc='upper right', fontsize=11)
plt.tight_layout()

# 6. Guardar y mostrar
plt.savefig("densidad_radial_plot_barras.png", dpi=300)
print("¡Gráfica generada y guardada como 'densidad_radial_plot_barras.png'!")
plt.show()