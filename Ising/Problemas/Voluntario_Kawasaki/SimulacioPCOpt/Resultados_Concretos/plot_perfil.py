import numpy as np
import matplotlib.pyplot as plt
import os

# ==============================================================================
# CONFIGURACIÓN DEL SCRIPT
# ==============================================================================
archivo_perfil = 'perfil128_m0_0.80.txt'
L = 128

# Seleccionamos temperaturas específicas para que el gráfico no esté saturado
temperaturas_a_graficar = [1.5, 2.3, 3.1] 

if not os.path.exists(archivo_perfil):
    print(f"Error: No se encontró el archivo {archivo_perfil}")
    exit()

datos = np.loadtxt(archivo_perfil)

# ==============================================================================
# PROCESAMIENTO Y GRAFICADO (Gráfico de Barras)
# ==============================================================================
plt.figure(figsize=(10, 6))

temperaturas_unicas = np.unique(datos[:, 0])
# Definimos un ancho para las barras y un desfase (offset) para que no se solapen
width = 0.25 
offsets = [-width, 0, width] 

for idx, T_target in enumerate(temperaturas_a_graficar):
    if T_target in temperaturas_unicas:
        datos_T = datos[datos[:, 0] == T_target]
        
        eje_Y = datos_T[:, 1]
        densidad = datos_T[:, 2] / L # Densidad = (nº espines +1) / L
        
        # Graficamos barras: 
        # eje_Y + offset para que las barras de distintas T queden una al lado de la otra
        plt.bar(eje_Y + offsets[idx], densidad, width=width, 
                label=f'T = {T_target} $J/k_B$', alpha=0.7)

# ==============================================================================
# FORMATO FINAL Y GUARDADO
# ==============================================================================
plt.title(f'Perfil de Densidad para N = {L}', fontsize=14)
plt.xlabel('Coordenada Y (Posición de la fila)', fontsize=12)
plt.ylabel('Densidad (Fracción de espines +1)', fontsize=12)

plt.ylim(-0.05, 1.05)
plt.grid(axis='y', linestyle='--', alpha=0.7)
plt.legend()

plt.tight_layout()
plt.savefig(f'perfil_densidad_{L}_.png', dpi=300)
print(f"Gráfica de barras guardada como 'perfil_densidad_barras_{L}_.png'")
plt.show()