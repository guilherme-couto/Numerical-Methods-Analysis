import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import imageio.v2
import math

def create_gif(num_threads, dt_ODE, dt_PDE, method, cell_model):

    times = []

    frames_file = f'./simulation-files/{cell_model}/{method}/frames-{num_threads}-{dt_ODE}-{dt_PDE}.txt'
    f = open(frames_file, 'r')
    lines = f.readlines()
    discretization_size = len(lines[1].split())

    for i in range(len(lines)):
        line = lines[i].split()
        lines[i] = [float(x) for x in line]
        if len(line) == 1:
            times.append(float(line[0]))

    totalframes = len(times)
    frame_rate = math.ceil(totalframes / 100)

    frames = []
    i = 1
    for n in range(totalframes):
        if n % frame_rate == 0:
            frame_name = f'frame-{n}.png'
            frames.append(frame_name)
            
            plt.imshow(lines[i:i+discretization_size], cmap='plasma', vmin=0, vmax=100)
            plt.colorbar(label='V (mV)')
            plt.title(f'{cell_model} {method} dtODE = {dt_ODE} dtPDE = {dt_PDE} t = {times[n]:.2f}')
            plt.xticks([])
            plt.yticks([])
            
            plt.savefig(frame_name)
            plt.close()

            i += (discretization_size+1)

    # Create gif directory
    if not os.path.exists(f'./gifs/{cell_model}/{method}'):
        os.makedirs(f'./gifs/{cell_model}/{method}')
        
    # Build gif
    with imageio.v2.get_writer(f'./gifs/{cell_model}/{method}/{num_threads}-{dt_ODE}-{dt_PDE}.gif', mode='I') as writer:
        for frame in frames:
            image = imageio.v2.imread(frame)
            writer.append_data(image)

    # Remove files
    for png in set(frames):
        os.remove(png)

def main():
    # Possible parameters
    cell_model = 'AFHN'
    numbers_threads = [4]
    dts_ODE = [0.02, 0.04, 0.06, 0.08, 0.1, 0.2]
    methods = ['ADI1', 'SSI-ADI', 'ADI1.5', 'FE']
    
    for number_threads in numbers_threads:
        for method in methods:
            for dt_ODE in dts_ODE:
                dts_PDE = [dt_ODE]
                if method == 'ADI1' or method == 'FE':
                    i = 2
                    while dt_ODE * i <= 0.2:
                        dts_PDE.append(dt_ODE * i)
                        i += 1
                for dt_PDE in dts_PDE:
                    create_gif(number_threads, f'{dt_ODE:.3f}', f'{dt_PDE:.3f}', method, cell_model)
                    print(f'Gif created for {number_threads} threads, {dt_ODE} dt_ODE, {dt_PDE} dt_PDE and {method} method.')

if __name__ == '__main__':
    main()
    