import numpy as np
import matplotlib.pyplot as plt
import os
import math

# Possible parameters
dxs = ['0.005']
cell_models = ['AFHN-Fibro']
numbers_threads = [6]
dts_ODE = [0.01, 0.02, 0.04, 0.08, 0.16]
max_dt_PDE = 0.0
methods = ['ADI1', 'SSI-ADI', 'ADI1.5.3']
repeats_to_speedup = 5

# Single case
# cell_models = ['AFHN-Fibro']
# numbers_threads = [2]
# dxs = ['0.020']
# dts_ODE = [0.02]
# max_dt_PDE = 0.2
# methods = ['SSI-ADI']