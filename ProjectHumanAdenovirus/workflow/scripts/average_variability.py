import numpy as np
import sys

# Redirect standard output and error to the log file defined in the Snakemake rule
sys.stderr = sys.stdout = open(snakemake.log[0], "w")

WINDOW_SIZE = int(snakemake.params.window_size)
variability = np.loadtxt(snakemake.input.variability, delimiter=",", dtype=float)

average_variability = []

for i in range(0,len(variability),1):   # fix Minor: overlapping windows
    window_variability = variability[i:i+WINDOW_SIZE]
    average_variability.append(np.mean(window_variability))
    print(f"Window {i//WINDOW_SIZE + 1}: Average Variability = {average_variability[-1]}")

np.savetxt(snakemake.output.variability, average_variability, fmt="%.6f", delimiter="\t")

