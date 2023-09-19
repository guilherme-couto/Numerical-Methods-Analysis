from execution import *

def read_file_matrix_to_vector(filename):
    f = open(filename, 'r')
    file_matrix = f.readlines()

    for i in range(len(file_matrix)):
        line = file_matrix[i].split()
        file_matrix[i] = [float(x) for x in line]

    vector = []
    for i in range(len(file_matrix)):
        for j in range(len(file_matrix[i])):
            vector.append(file_matrix[i][j])
    return np.array(vector)
    

def main():
    for dx in dxs:
        for cell_model in cell_models:
            # Reference solution
            # Check if reference file exists
            reference_filename = f'./simulation-files/{dx}/{cell_model}/{reference_method}/last-{numbers_threads[0]}-{reference_dt}-{reference_dt}.txt'
            
            if not os.path.exists(reference_filename):
                print(f'Reference file {reference_filename} does not exist. Creating it...')
                execution_line = f'./main {numbers_threads[0]} {reference_dt} {reference_dt} {reference_method} 1 0 0'
                print(f'Executing {execution_line}')
                os.system(execution_line)
                print(f'Simultation {execution_line} finished!\n')
            
            reference = read_file_matrix_to_vector(reference_filename)

            # Open file to write errors of each method configuration
            if not os.path.exists(f'./analysis/{dx}/{cell_model}'):
                os.makedirs(f'./analysis/{dx}/{cell_model}')
            f_errors = open(f'./analysis/{dx}/{cell_model}/errors.txt', 'w')

            for number_threads in numbers_threads:
                for method in methods:
                    for dt_ODE in dts_ODE:
                        dts_PDE = [dt_ODE]
                        if method == 'OS-ADI' or method == 'FE':
                            i = 2
                            while dt_ODE * i <= max_dt_PDE:
                                dts_PDE.append(dt_ODE * i)
                                i += 1
                        for dt_PDE in dts_PDE:
                            filename = f'./simulation-files/{dx}/{cell_model}/{method}/last-{number_threads}-{dt_ODE:.3f}-{dt_PDE:.3f}.txt'
                            case = read_file_matrix_to_vector(filename)

                            f_errors.write(f'\nError for {method} with dx = {dx}, dtODE = {dt_ODE:.3f} and dtPDE = {dt_PDE:.3f}:\n')

                            # If the simulation diverged
                            if np.any(np.isnan(case)):
                                f_errors.write('Simulation diverged\n')

                            else:
                                mae = np.mean(np.abs(case - reference))
                                mse = np.mean((case - reference)**2)
                                relative_error = np.linalg.norm(case - reference) / np.linalg.norm(reference)

                                f_errors.write(f'Mean absolute error: {mae:.4f}\n')
                                f_errors.write(f'Mean square error: {mse:.4f}\n')
                                f_errors.write(f'Relative error: {relative_error:.4f} ({100*relative_error:.4f} %)\n')
                                
                    f_errors.write('\n----------------------------------------------------------------------------\n')
                    f_errors.write('----------------------------------------------------------------------------\n\n')

if __name__ == '__main__':
    main()