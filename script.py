import os

from execution import *

for cell_model in cell_models:
    for num_threads in numbers_threads:
        for method in methods:
            for dt_ODE in dts_ODE:
                dts_PDE = [dt_ODE]
                if method == 'ADI1' or method == 'FE':
                    i = 2
                    while dt_ODE * i <= max_dt_PDE:
                        dts_PDE.append(dt_ODE * i)
                        i += 1
                for dt_PDE in dts_PDE:
                    execution_line = ''
                    if cell_model.find('Fibro') != -1:
                        execution_line = f'./main {num_threads} {dt_ODE} {dt_PDE} {method} 1'
                    else:
                        execution_line = f'./main {num_threads} {dt_ODE} {dt_PDE} {method} 0'
                    print(f'Executing {execution_line}\n')
                    os.system(f'{execution_line}')
                    print(f'Simultation {execution_line} finished!\n')
                    