#!/usr/bin/env python3

import argparse
from pathlib import Path
import numpy as np
from matplotlib import pyplot as plt


def main():
    parser = argparse.ArgumentParser(
        description="Genera una gráfica de la energía total del sistema en función del tiempo desde estado_estacionario.dat"
    )
    parser.add_argument(
        "--input",
        default="estado_estacionario.dat",
        help="Fichero de entrada con tiempo y energía (por defecto: estado_estacionario.dat)",
    )
    parser.add_argument(
        "--output",
        default="energia_vs_tiempo.png",
        help="Nombre del fichero de salida de la gráfica (por defecto: energia_vs_tiempo.png)",
    )
    parser.add_argument(
        "--show",
        action="store_true",
        help="Mostrar la gráfica en pantalla además de guardarla",
    )
    args = parser.parse_args()

    data_path = Path(args.input)
    if not data_path.exists():
        raise SystemExit(f"Error: no existe el fichero de entrada '{args.input}'")

    data = np.loadtxt(data_path)
    if data.ndim == 1:
        data = data.reshape(1, -1)

    if data.shape[1] < 2:
        raise SystemExit("Error: el fichero de entrada debe tener al menos dos columnas: tiempo y energía")

    tiempos = data[:, 0]
    energias = data[:, 1]

    plt.figure(figsize=(10, 6))
    plt.plot(tiempos, energias, color="tab:blue", linewidth=1.5)
    plt.title("Energía total del sistema vs. tiempo")
    plt.xlabel("Tiempo [s]")
    plt.ylabel("Energía total [J]")
    plt.grid(True, linestyle="--", alpha=0.5)
    plt.tight_layout()
    plt.savefig(args.output, dpi=150)
    print(f"Gráfica guardada en '{args.output}'")

    if args.show:
        plt.show()


if __name__ == "__main__":
    main()
