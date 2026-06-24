#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
grafica_magnetizacion.py
Lee uno o varios ficheros de datos de magnetización (T, M, error)
y dibuja las curvas con barras de error.

Uso:
  python3 grafica_magnetizacion.py magn16.txt magn32.txt \
      --labels "N=16" "N=32" --output magnetizacion.png --show

Si no se proporcionan etiquetas, se usan los nombres de fichero.
"""

import argparse
import os
import numpy as np
import matplotlib.pyplot as plt
import itertools


def load_data(path):
    data = np.genfromtxt(path)
    if data.size == 0:
        raise ValueError(f"Fichero vacío o formato no válido: {path}")
    if data.ndim == 1:
        if data.shape[0] < 3:
            raise ValueError(f"Fichero con menos de 3 columnas: {path}")
        data = data.reshape(1, -1)
    if data.shape[1] < 3:
        raise ValueError(f"Fichero con menos de 3 columnas: {path}")
    return data[:, 0], data[:, 1], data[:, 2]


def plot_files(files, labels=None, output='magnetizacion.png', show=False):
    plt.figure(figsize=(7, 5))
    colors = itertools.cycle(plt.get_cmap('tab10').colors)

    for i, fp in enumerate(files):
        label = labels[i] if labels and i < len(labels) else os.path.basename(fp)
        try:
            T, M, err = load_data(fp)
        except Exception as e:
            print(f"No se pudo leer '{fp}': {e}")
            continue

        color = next(colors)
        plt.errorbar(T, M, yerr=err, fmt='o-', label=label, capsize=3, markersize=5, color=color)

    plt.xlabel('Temperatura (T)')
    plt.ylabel('Magnetización |M|')
    plt.title('Curva de magnetización con barras de error')
    plt.grid(alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(output, dpi=200)
    print(f"Gráfica guardada en: {output}")
    if show:
        plt.show()


def main():
    parser = argparse.ArgumentParser(description='Dibuja magnetización con barras de error')
    parser.add_argument('files', nargs='+', help='Ficheros de datos (T M error)')
    parser.add_argument('--labels', nargs='*', help='Etiquetas para cada fichero')
    parser.add_argument('--output', default='magnetizacion.png', help='Nombre del fichero de salida')
    parser.add_argument('--show', action='store_true', help='Mostrar la figura además de guardarla')
    args = parser.parse_args()

    plot_files(args.files, labels=args.labels, output=args.output, show=args.show)


if __name__ == '__main__':
    main()
