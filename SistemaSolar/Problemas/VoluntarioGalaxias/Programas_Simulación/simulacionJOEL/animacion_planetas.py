# ================================================================================
# ANIMACION SISTEMA SOLAR (VERSIÓN ALTA EFICIENCIA)
# ================================================================================

from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.patches import Circle
import numpy as np

try:
    import cupy as cp
    gpu_available = True
except ImportError:
    cp = None
    gpu_available = False

# Parámetros
# ========================================
file_in = "posiciones_planetas.dat" 
file_out = "planetas" 
use_gpu = True  

# --- NUEVOS PARÁMETROS DE OPTIMIZACIÓN ---
frame_step = 100          # PROCESA 1 DE CADA X FOTOGRAMAS. Auméntalo si va lento (ej. 50 o 100).
max_trail_length = 50    # LONGITUD MÁXIMA DE LA ESTELA. Evita que la memoria colapse.
# -----------------------------------------

x_min = -1
x_max = 1
y_min = -1 
y_max = 1

interval = 20 # Tiempo entre fotogramas en milisegundos
show_trail = True 
trail_width = 1 
save_to_file = True 
dpi = 150 

planet_radius = 0.01 
planet_colors = ["tab:blue", "tab:orange", "tab:green", "tab:red",
                 "tab:purple", "tab:brown", "tab:pink", "tab:gray"]


# Lectura del fichero de datos
# ========================================
print("Leyendo fichero de datos...")
with open(file_in, "r") as f:
    data_str = f.read()

xp = cp if gpu_available and use_gpu else np
if gpu_available and use_gpu:
    def to_cpu(array): return xp.asnumpy(array)
else:
    def to_cpu(array): return array

frames_data = list()

# OPTIMIZACIÓN 1: Descartar bloques de texto antes de procesarlos numéricamente
bloques_texto = [b for b in data_str.split("\n\n") if b.strip()]
bloques_texto = bloques_texto[::frame_step] # Aplicar el salto de fotogramas

print(f"Procesando {len(bloques_texto)} fotogramas (con salto de {frame_step})...")

for frame_data_str in bloques_texto:
    frame_data = list()
    for planet_pos_str in frame_data_str.split("\n"):
        if planet_pos_str.strip():
            planet_pos = xp.fromstring(planet_pos_str, sep=",")
            if planet_pos.size > 0:
                frame_data.append(planet_pos)

    if frame_data:
        frames_data.append(xp.vstack(frame_data))

nplanets = int(frames_data[0].shape[0])

# OPTIMIZACIÓN 2: Precalcular el historial en una matriz 3D para acceso instantáneo
# Dimensión: (num_fotogramas, num_planetas, coordenadas_xy)
print("Construyendo matriz histórica 3D para acelerar el renderizado...")
historial = np.array([to_cpu(f) for f in frames_data])


# Creación de la animación/gráfico
# ========================================
print("Generando animación...")
fig, ax = plt.subplots()
ax.axis("equal")  
ax.set_xlim(x_min, x_max)
ax.set_ylim(y_min, y_max)
ax.set_facecolor('black') # Fondo negro opcional, queda más "espacial"

if not hasattr(planet_radius, "__iter__"):
    planet_radius = planet_radius*np.ones(nplanets)
else:
    if not nplanets == len(planet_radius):
        raise ValueError("El número de radios no coincide con el de planetas")

if not hasattr(planet_colors, "__iter__") or isinstance(planet_colors, str):
    planet_colors = [planet_colors]*nplanets
elif len(planet_colors) < nplanets:
    planet_colors = (planet_colors * ((nplanets // len(planet_colors)) + 1))[:nplanets]

planet_points = list()
planet_trails = list()

for j_planet, (radius, color) in enumerate(zip(planet_radius, planet_colors)):
    x, y = historial[0, j_planet, 0], historial[0, j_planet, 1]
    
    planet_point = Circle((x, y), radius, facecolor=color, edgecolor="none", zorder=3)
    ax.add_artist(planet_point)
    planet_points.append(planet_point)

    if show_trail:
        planet_trail, = ax.plot(x, y, "-", linewidth=trail_width, color=color, alpha=0.6, zorder=2)
        planet_trails.append(planet_trail)
 
def update(j_frame):
    for j_planet in range(nplanets):
        x = historial[j_frame, j_planet, 0]
        y = historial[j_frame, j_planet, 1]
        planet_points[j_planet].center = (x, y)

        if show_trail:
            # OPTIMIZACIÓN 3: Extraer la estela con slicing (super rápido) y con longitud límite
            inicio_estela = max(0, j_frame - max_trail_length)
            xs_new = historial[inicio_estela:j_frame+1, j_planet, 0]
            ys_new = historial[inicio_estela:j_frame+1, j_planet, 1]
            planet_trails[j_planet].set_data(xs_new, ys_new)

    return planet_points + planet_trails

def init_anim():
    if show_trail:
        for j_planet in range(nplanets):
            planet_trails[j_planet].set_data([], [])
    return planet_points + planet_trails

nframes = len(frames_data)

if nframes > 1:
    animation = FuncAnimation(
            fig, update, init_func=init_anim,
            frames=nframes, blit=True, interval=interval)

    if save_to_file:
        animation.save("{}.mp4".format(file_out), dpi=dpi)
        print(f"Vídeo guardado con éxito como {file_out}.mp4")
    else:
        plt.show()
else:
    if save_to_file:
        fig.savefig("{}.pdf".format(file_out))
        print("Imagen guardada con éxito.")
    else:
        plt.show()