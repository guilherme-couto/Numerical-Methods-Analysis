import numpy as np
import matplotlib.pyplot as plt
import os
import imageio.v2
import math

# Possible parameters
dxs = ['0.020']
cell_models = ['AFHN', 'AFHN-Fibro']
numbers_threads = [4]
dts_ODE = [0.02, 0.04, 0.06, 0.08, 0.1, 0.2]
max_dt_PDE = 0.2
methods = ['ADI1', 'SSI-ADI', 'ADI1.5', 'FE']

# Single case
# cell_models = ['AFHN']
# numbers_threads = [4]
# dts_ODE = [0.2]
# max_dt_PDE = 0.2
# methods = ['SSI-ADI']