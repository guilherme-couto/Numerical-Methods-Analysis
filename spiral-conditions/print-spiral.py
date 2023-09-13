import numpy as np
import matplotlib.pyplot as plt

with open('lastSpiralV.txt', 'r') as f:
    lines = f.readlines()
    file_matrix = []
    for line in lines:
        line = line.split()
        file_matrix.append([])
        for item in line:
            file_matrix[-1].append(float(item))
        

# Plot the matrix as a heatmap
plt.imshow(file_matrix, cmap='plasma', vmin=0, vmax=100)
plt.colorbar(label='V (mV)')
plt.title('Spiral')
plt.xticks([])
plt.yticks([])
plt.savefig('spiral.png')