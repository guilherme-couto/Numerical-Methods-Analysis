import numpy as np
import matplotlib.pyplot as plt

# Read csv file in format X,Y,Z,a,b,c,V where X,Y,Z are the coordinates of the point, a,b,c are the semi-axes of the ellipsoid and V is the potential at that point
# Store the data in a matrix where i and j represent X and Y coordinates and the value in the matrix is the potential at that point

with open('V_it_15000.txt', 'r') as f:
    lines = f.readlines()
    file_matrix = np.zeros((100, 100))
    for line in lines:
        line = line.split(',')
        x = int((float(line[0]) / 200))
        y = int((float(line[1]) / 200))
        V = float(line[6])
        file_matrix[x][y] = V

# Plot the matrix as a heatmap
plt.imshow(file_matrix, cmap='plasma', vmin=0, vmax=100)
plt.colorbar(label='V (mV)')
plt.title('V at t = 15000')
plt.xticks([])
plt.yticks([])
plt.savefig('V_it_15000.png')
plt.show()