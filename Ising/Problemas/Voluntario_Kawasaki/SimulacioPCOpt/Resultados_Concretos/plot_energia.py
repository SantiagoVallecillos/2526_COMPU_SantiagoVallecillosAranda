import numpy as np
import matplotlib.pyplot as plt
import os

# ==============================================================================
# CONFIGURACIÓN DEL SCRIPT
# ==============================================================================
# Nombres de los archivos a leer. Cámbialos por los que quieras analizar (ej: magn32.txt para m0=0)
archivos = ['magn32_m0_0.50.txt', 'magn64_m0_0.50.txt', 'magn128_m0_0.50.txt']
etiquetas = ['N = 32', 'N = 64', 'N = 128']
colores = ['blue', 'orange', 'green']

# Crear la figura (1 solo gráfico grande para ver bien las barras de error)
fig, ax = plt.subplots(figsize=(8, 6))
fig.suptitle('Energía Interna - Modelo de Ising (Dinámica de Kawasaki)', fontsize=16)

# ==============================================================================
# LECTURA Y GRAFICADO DE DATOS
# ==============================================================================
for archivo, etiqueta, color in zip(archivos, etiquetas, colores):
    if not os.path.exists(archivo):
        print(f"Advertencia: No se encontró el archivo {archivo}. Se omitirá.")
        continue

    # Cargar los datos. np.loadtxt separa automáticamente por tabulaciones o espacios
    # Columnas: 0:T | 1:<M> | 2:err_M | 3:<E> | 4:err_E | 5:c_N | 6:chi_N
    datos = np.loadtxt(archivo)
    
    T = datos[:, 0]
    E = datos[:, 3]       # Columna de la Energía media
    err_E = datos[:, 4]   # Columna del Error de la Energía media

    # Gráfica de Energía vs Temperatura con barras de error
    # Usamos errorbar en lugar de plot para visualizar la desviación estadística
    ax.errorbar(T, E, yerr=err_E, fmt='-o', markersize=5, capsize=4, 
                color=color, label=etiqueta, linewidth=1.5, alpha=0.8)

# ==============================================================================
# FORMATO FINAL Y GUARDADO
# ==============================================================================
ax.set_xlabel('Temperatura ($J/k_B$)', fontsize=12)
ax.set_ylabel('Energía media por partícula $\\langle e \\rangle$ ($J$)', fontsize=12)
ax.set_title('Evolución de la energía frente a la temperatura', fontsize=14)
ax.grid(True, linestyle='--', alpha=0.6)

# Añadir la leyenda
ax.legend(fontsize=11)

plt.tight_layout() # Ajusta los márgenes
plt.subplots_adjust(top=0.90) # Deja espacio para el título principal

# Guardar la gráfica en una imagen de alta calidad y mostrarla en pantalla
plt.savefig('energia_media_resultados_m0_0.50.png', dpi=300)
print("Gráfica guardada como 'energia_media_resultados_m0_0.50.png'")
plt.show()