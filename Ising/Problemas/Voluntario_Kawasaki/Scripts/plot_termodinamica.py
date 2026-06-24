import numpy as np
import matplotlib.pyplot as plt
import os

# ==============================================================================
# CONFIGURACIÓN DEL SCRIPT
# ==============================================================================
# Nombres de los archivos a leer. Puedes cambiarlos por los terminados en "_m0.txt"
archivos = ['magn32.txt', 'magn64.txt', 'magn128.txt']
etiquetas = ['N = 32', 'N = 64', 'N = 128']
colores = ['blue', 'orange', 'green']

# Crear la figura con 3 subgráficos (1 fila, 3 columnas)
fig, axs = plt.subplots(1, 3, figsize=(15, 5))
fig.suptitle('Análisis Termodinámico - Modelo de Ising (Dinámica de Kawasaki)', fontsize=16)

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
    M = datos[:, 1]
    c_N = datos[:, 5]
    chi_N = datos[:, 6]

    # 1. Gráfica de Magnetización vs Temperatura
    axs[0].plot(T, M, marker='o', markersize=4, linestyle='-', color=color, label=etiqueta)
    axs[0].set_xlabel('Temperatura ($J/k_B$)')
    axs[0].set_ylabel('Magnetización por partícula $\\langle m \\rangle$')
    axs[0].set_title('Curva de Magnetización')
    axs[0].grid(True, linestyle='--', alpha=0.6)

    # 2. Gráfica de Calor Específico vs Temperatura
    axs[1].plot(T, c_N, marker='^', markersize=4, linestyle='-', color=color, label=etiqueta)
    axs[1].set_xlabel('Temperatura ($J/k_B$)')
    axs[1].set_ylabel('Calor Específico $c_N$')
    axs[1].set_title('Calor Específico')
    axs[1].grid(True, linestyle='--', alpha=0.6)

    # 3. Gráfica de Susceptibilidad Magnética vs Temperatura
    axs[2].plot(T, chi_N, marker='s', markersize=4, linestyle='-', color=color, label=etiqueta)
    axs[2].set_xlabel('Temperatura ($J/k_B$)')
    axs[2].set_ylabel('Susceptibilidad $\\chi_N$')
    axs[2].set_title('Susceptibilidad Magnética')
    axs[2].grid(True, linestyle='--', alpha=0.6)

# ==============================================================================
# FORMATO FINAL Y GUARDADO
# ==============================================================================
# Poner la leyenda en el primer gráfico
axs[0].legend()

plt.tight_layout() # Ajusta los márgenes para que no se superpongan los textos
plt.subplots_adjust(top=0.88) # Deja espacio para el título principal

# Guardar la gráfica en una imagen de alta calidad y mostrarla en pantalla
plt.savefig('termodinamica_resultados.png', dpi=300)
print("Gráfica guardada como 'termodinamica_resultados.png'")
plt.show()