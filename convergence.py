import numpy as np
import matplotlib.pyplot as plt
import os
import math

# Simulation parameters (parameters.c)
# L = 2 cm
# T = 160 ms

# Parameters (old)
#dxs = ['0.010', '0.020', '0.025']
#dts_second = ['0.08000', '0.16000', '0.20000'] # A = 8 (dt = A dx)
#dts_first = ['0.03200', '0.12800', '0.20000'] # A = 320 (dt = A dx²)
#dts_second = ['0.010', '0.020', '0.025']
#dts_first = ['0.00010', '0.00040', '0.000625'] # A = 1 (dt = A dx²)

cell_models = ['AFHN']#, 'AFHN-Fibro']
numbers_threads = [6]
methods_first = ['OS-ADI', 'SSI-ADI', 'MOSI-ADI.3']
methods_second = ['OS-ADI', 'SSI-ADI', 'MOSI-ADI.3']

dxs = ['0.020', '0.025', '0.040', '0.050']
dts_first = ['0.00400', '0.00625', '0.01600', '0.02500']           # A = 10 (dt = A dx²)
dts_second = ['0.02000', '0.02500', '0.04000', '0.05000']           # A = 1 (dt = A dx)

dx_ref = 0.001
dt_ref_second = 0.001   # A = 1 (dt = A dx)
dt_ref_first = 0.00001  # A = 10 (dt = A dx²)
ref_methods = ['SSI-ADI']
ref = 0.001
spatial_ref = 0.2

print('For second order methods:')
for i in range(len(dxs)):
    print(f'dx = {dxs[i]}, dt = {dts_second[i]}')
print('\nFor first order methods:')
for i in range(len(dxs)):
    print(f'dx = {dxs[i]}, dt = {dts_first[i]}')

# Run convergence cases
# Compile
os.system(f'gcc -o main main.c -O3 -lm -lpthread -fopenmp')
for cell_model in cell_models:
    for num_threads in numbers_threads:
        for method in methods_first:
            for i in range(len(dxs)):
                spatial_rate = int(spatial_ref / float(dxs[i]))
                rate = int(float(dxs[i]) / ref)
                
                # Second order
                path = f'./simulation-files/{dxs[i]}/{cell_model}/{method}/last-{numbers_threads[-1]}-{(dts_second[i])}.txt'
                if not os.path.exists(path):
                    execution_line = f'./main {num_threads} {dxs[i]} {dts_second[i]} {method} {rate} {spatial_rate}'
                    print(f'Executing {execution_line}')
                    os.system(f'{execution_line}')
                    print(f'Simultation {execution_line} finished!\n')

                # First order
                path = f'./simulation-files/{dxs[i]}/{cell_model}/{method}/last-{numbers_threads[-1]}-{(float(dts_first[i])):.5f}.txt'
                if not os.path.exists(path):
                    execution_line = f'./main {num_threads} {dxs[i]} {dts_first[i]} {method} {rate} {spatial_rate}'
                    print(f'Executing {execution_line}')
                    os.system(f'{execution_line}')
                    print(f'Simultation {execution_line} finished!\n')

# Run reference case
print(f'Reference case {(float(dxs[0])/dx_ref):.2f} times smaller')
print(f'dx_ref = {dx_ref}, dt_ref_second = {dt_ref_second}')
print(f'dx_ref = {dx_ref}, dt_ref_first = {dt_ref_first}\n')

spatial_rate = int(spatial_ref / float(dx_ref))
rate = int(dx_ref / ref)

reference_path_second = f'./simulation-files/{(dx_ref):.3f}/{cell_models[0]}/{ref_methods[0]}/last-{numbers_threads[-1]}-{(dt_ref_second):.5f}.txt'
reference_path_first = f'./simulation-files/{(dx_ref):.3f}/{cell_models[0]}/{ref_methods[0]}/last-{numbers_threads[-1]}-{(dt_ref_first):.5f}.txt'

for ref_method in ref_methods:

    # For second order
    if not os.path.exists(reference_path_second):
        execution_line = f'./main {numbers_threads[-1]} {dx_ref} {dt_ref_second} {ref_method} {rate} {spatial_rate}'
        print(f'Executing Reference: {execution_line}')
        os.system(f'{execution_line}')
        print(f'Simultation Reference: {execution_line} finished!\n')

    # For first order
    if not os.path.exists(reference_path_first):
        execution_line = f'./main {numbers_threads[-1]} {dx_ref} {dt_ref_first} {ref_method} {rate} {spatial_rate}'
        print(f'Executing Reference: {execution_line}')
        os.system(f'{execution_line}')
        print(f'Simultation Reference: {execution_line} finished!\n')
      

# Get reference solution

reference_second = []

reference_first = []
with open(reference_path_second, 'r') as file:
    for line in file:
        line = line.split()
        for value in line:
            reference_second.append(float(value))
print(f'Reference second length = {len(reference_second)}')
reference_second = np.array(reference_second)
ref_norm_second = np.linalg.norm(reference_second)

with open(reference_path_first, 'r') as file:
    for line in file:
        line = line.split()
        for value in line:
            reference_first.append(float(value))
print(f'Reference first length = {len(reference_first)}')
reference_first = np.array(reference_first)
ref_norm_first = np.linalg.norm(reference_first)


# Calculate errors and plot
cases_first = {}
cases_second = {}
for method in methods_first:
    cases_first[method] = []
for method in methods_second:
    cases_second[method] = []

for method in methods_second:
    for i in range(len(dxs)):
        # Second order methods
        path = f'./simulation-files/{dxs[i]}/{cell_models[0]}/{method}/last-{numbers_threads[-1]}-{(dts_second[i])}.txt'
        case = []
        with open(path, 'r') as file:
            for line in file:
                line = line.split()
                for value in line:
                    case.append(float(value))
        case = np.array(case)
        error = np.linalg.norm(case - reference_second) / np.linalg.norm(reference_second)
        cases_second[method].append(error)
            
        # First order methods
        if dts_first[i] == '0.000625':
            dts_first[i] = '0.00062'
        if dts_first[i] == '0.003125':
            dts_first[i] = '0.00312'
        path = f'./simulation-files/{dxs[i]}/{cell_models[0]}/{method}/last-{numbers_threads[-1]}-{((dts_first[i]))}.txt'
        case = []
        with open(path, 'r') as file:
            for line in file:
                line = line.split()
                for value in line:
                    case.append(float(value))
        case = np.array(case)
        error = np.linalg.norm(case - reference_first) / np.linalg.norm(reference_first)
        cases_first[method].append(error)

# Plot log-log
dxs_float = [float(dx) for dx in dxs]
print('\ndxs_float =', dxs_float)
dts_second_float = [float(dt) for dt in dts_second]
plt.figure()
for method in methods_second:
    plt.loglog(dts_second_float, cases_second[method], label=method)
plt.legend()
plt.xlabel('dt')
plt.ylabel('Error')
plt.title(f'Convergence for second order methods')
plt.grid()
plt.savefig(f'./simulation-files/convergence_plot_2nd.png')

dts_first_float = [float(dt) for dt in dts_first]
plt.figure()
for method in methods_first:
    plt.loglog(dts_first_float, cases_first[method], label=method)
plt.legend()
plt.xlabel('dt')
plt.ylabel('Error')
plt.title(f'Convergence for first order methods')
plt.grid()
plt.savefig(f'./simulation-files/convergence_plot_1st.png')

# Print errors
slopes_file = open('slopes.txt', 'w') 
print('\nErrors for second order methods:')
slopes_file.write('\nErrors for second order methods:\n')
for method in methods_second:
    print(f'{method} = {cases_second[method]}')
    slopes_file.write(f'{method} = {cases_second[method]}\n')

print('\nErrors for first order methods:')
slopes_file.write('\nErrors for first order methods:\n')
for method in methods_first:
    print(f'{method} = {cases_first[method]}')
    slopes_file.write(f'{method} = {cases_first[method]}\n')

# Print slopes
print('\nSlopes for second order methods (last and first points):')
slopes_file.write('\nSlopes for second order methods (last and first points):\n')
for method in methods_second:
    print(f'{method} = {(math.log(cases_second[method][-1] / cases_second[method][0]) / math.log(dts_second_float[-1] / dts_second_float[0])):.3f}')
    slopes_file.write(f'{method} = {(math.log(cases_second[method][-1] / cases_second[method][0]) / math.log(dts_second_float[-1] / dts_second_float[0])):.3f}\n')

print('\nSlopes for first order methods (last and first points):')
slopes_file.write('\nSlopes for first order methods (last and first points):\n')
for method in methods_first:
    print(f'{method} = {(math.log(cases_first[method][-1] / cases_first[method][0]) / math.log(dts_first_float[-1] / dts_first_float[0])):.3f}')
    slopes_file.write(f'{method} = {(math.log(cases_first[method][-1] / cases_first[method][0]) / math.log(dts_first_float[-1] / dts_first_float[0])):.3f}\n')

    
   
# print('\nSlopes for second order methods (avg):')
# for method in methods_second:
# 	slope = []
# 	for i in range(1, len(cases_second[method])):
# 		slope.append(math.log(cases_second[method][i] / cases_second[method][i-1]) / math.log(dts_second_float[i] / dts_second_float[i-1]))
# 	slope = np.array(slope)
# 	print(f'{method} = {(np.mean(slope)):.3f}')

# print('\nSlopes for first order methods (avg):')
# for method in methods_first:
# 	slope = []
# 	for i in range(1, len(cases_first[method])):
# 		slope.append(math.log(cases_first[method][i] / cases_first[method][i-1]) / math.log(dts_first_float[i] / dts_first_float[i-1]))
# 	slope = np.array(slope)
# 	print(f'{method} = {(np.mean(slope)):.3f}')
