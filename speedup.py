from execution import *
import sys

# File to save speedups
speedup_file = open('speedup_analysis.txt', 'w')

# Create a dict to save results for each delta_t
for dt in dts_ODE:
    speedup_file.write('\n------------------------------------------------\n')
    speedup_file.write(f'For dt = {dt}')
    speedup_file.write('\n------------------------------------------------\n')
    cases = {}
    for method in methods:
        speedup_file.write(f'Method {method}\n')
        cases[method] = {}
        avgs = []
        for num_threads in numbers_threads:
            speedup_file.write(f'{num_threads}: ')
            results = []
            cases[method][num_threads] = results
            #with open('teste.txt', 'r') as f:
            with open(f'./simulation-files/{dxs[-1]}/{cell_models[0]}/{method}/infos-{num_threads}-{(dt):.3f}-{(dt):.3f}.txt', 'r') as f:
                lines = f.readlines()
                for line in lines:
                    line = line.split()
                    cases[method][num_threads].append(float(line[-2]))
                    speedup_file.write(f'{float(line[-2])} ')
            results = cases[method][num_threads]
            results.remove(max(results))
            results.remove(min(results))
            results = np.array(results)
            avgs.append([num_threads, np.mean(results), np.std(results)])
            speedup_file.write(f' |  Mean: {(np.mean(results)):.3f} +- {(np.std(results)):.3f}')
            speedup_file.write(f'\n')
        speedup_file.write(f'\n')
        speedup_file.write(f'Speedups: ')
        for i in range(len(avgs)-1):
            speedup_file.write(f'({avgs[i][0]}) {(avgs[-1][1]/avgs[i][1]):.2f} | ')
        speedup_file.write(f'\n\n')

    

