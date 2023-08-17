from execution import *
import sys

# Check arguments
if len(sys.argv) == 1:
    for dx in dxs:
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
                                execution_line = f'./main {num_threads} {dt_ODE} {dt_PDE} {method} 1 0 0'
                            else:
                                execution_line = f'./main {num_threads} {dt_ODE} {dt_PDE} {method} 0 0 0'
                            print(f'Executing {execution_line}\n')
                            os.system(f'{execution_line}')
                            print(f'Simultation {execution_line} finished!\n')
                            
elif len(sys.argv) == 2:
    if sys.argv[1] == 'speedup':
        print('Executing speedup test\n')
        for dx in dxs:
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
                                    execution_line = f'./main {num_threads} {dt_ODE} {dt_PDE} {method} 1 0 1'
                                else:
                                    execution_line = f'./main {num_threads} {dt_ODE} {dt_PDE} {method} 0 0 1'
                                print(f'Executing {execution_line}\n')
                                os.system(f'{execution_line}')
                                print(f'Simultation {execution_line} finished!\n')
    else:
        print('Invalid argument\n')
        exit()   
                    