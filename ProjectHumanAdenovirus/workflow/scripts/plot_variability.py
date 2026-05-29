import matplotlib.pyplot as plt
import numpy as np

WINDOW_SIZE = int(snakemake.params.window_size)
variability = np.loadtxt(snakemake.input.variability, delimiter=",", dtype=float)
# Create x-axis values based on the number of windows
x_values = np.arange(1, len(variability)*WINDOW_SIZE + 1, WINDOW_SIZE)

plt.figure(figsize=(10, 6))
plt.plot(x_values, variability, marker='o', linestyle='-', color='blue')
plt.title(f'Average Sequence Variability (Window Size = {WINDOW_SIZE})')
plt.xlabel('Genome Position (bp)')
plt.ylabel('Average Variability (Shannon Entropy)')
plt.grid()
plt.savefig(snakemake.output.variability_plot)
plt.close()
