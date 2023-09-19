import numpy as np
import matplotlib.pyplot as plt
import os
import math

# Possible parameters
dxs = ['0.005']
cell_models = ['MV-Fibro']
numbers_threads = [6]
dts_ODE = [0.01, 0.02, 0.04, 0.08, 0.16]
max_dt_PDE = 0.0
methods = ['OS-ADI', 'SSI-ADI', 'MOSI-ADI', 'MOSI-ADI.2', 'MOSI-ADI.3', 'FE']
repeats_to_speedup = 5
reference_method = 'SSI-ADI'
reference_dt = 0.001

# Single case
# cell_models = ['MV-Fibro']
# numbers_threads = [4]
# dxs = ['0.020']
# dts_ODE = [0.005]
# max_dt_PDE = 0.0
# methods = ['SSI-ADI']
# repeats_to_speedup = 5