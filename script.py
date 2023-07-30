import os

# Possible parameters
numbers_threads = [4]
dts_ODE = [0.02, 0.04, 0.06, 0.08, 0.1, 0.2]
methods = ['ADI1', 'SSI-ADI', 'ADI1.5', 'FE']

for num_threads in numbers_threads:
    for method in methods:
        for dt_ODE in dts_ODE:
            dts_PDE = [dt_ODE]
            if method == 'ADI1' or method == 'FE':
                i = 2
                while dt_ODE * i <= 0.2:
                    dts_PDE.append(dt_ODE * i)
                    i += 1
            for dt_PDE in dts_PDE:
                print(f'Executing ./main {num_threads} {dt_ODE} {dt_PDE} {method}')
                os.system(f'./main {num_threads} {dt_ODE} {dt_PDE} {method}')
                print(f'Simultation {num_threads} {dt_ODE} {dt_PDE} {method} finished!\n')
                