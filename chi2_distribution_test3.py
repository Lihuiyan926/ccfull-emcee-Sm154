import csv
import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import platform
from iminuit import Minuit
from iminuit.util import propagate
from matplotlib import gridspec
# from iminuit.cost import LeastSquares
from scipy import interpolate #插值操作
import pandas as pd
import shutil
import uuid
import time
from joblib import Parallel, delayed
from scipy.stats import norm
import itertools
import cupy as cp
# data = pd.read_csv('expdata_jia.csv')
# Load the data from CSV file
data = pd.read_csv('filterdata07.csv')

Eeff = data['Sm154_Eeff'].values  # Convert to NumPy array
Ratio = data['Sm154_Ratio'].values  # Convert to NumPy array
Error = data['Sm154_Err'].values  # Convert to NumPy array

# Use NumPy arrays for CPU processing
expx = Eeff
expy = Ratio
expyerr = Error

# Get the length of expx
expn = len(expx)

fchain = 'Sm154-beta3-sup1.csv'
if os.path.exists(fchain):  # if existed
    os.remove(fchain)


def model(expx, par):
    beta_200, beta_400, beta_300 = par[0], par[1], par[2]

    # Create a temporary directory
    uuid_str = uuid.uuid4().hex
    tmp_file_name = f'tmpfile_{uuid_str}'
    fn = tmp_file_name

    if not os.path.exists(fn):
        os.makedirs(fn)

    # Change to the temporary directory
    os.chdir(fn)

    # Copy input files
    shutil.copyfile('../ccfull-sc.inp', './ccfull-sc.inp')
    shutil.copyfile('../TEST2.INP', './TEST2.INP')
    shutil.copyfile('../a', './a')

    # Modify the parameters and write them back to the file
    with open("ccfull-sc.inp", 'r') as file:
        lines = file.readlines()
        lines[2] = f'0.082,{beta_200},{beta_400},3\n'
        lines[3] = f'1.012,{beta_300},3,1\n'

    with open("ccfull-sc.inp", 'w') as file:
        file.writelines(lines)

    # Execute the program
    sysstr = platform.system()
    main = "ccfull-sc2.exe" if sysstr == "Windows" else "./a"

    if os.path.exists(main):
        os.system(f'chmod 777 {main}' if sysstr == "Linux" else '')
        os.system(main)

    fresult = "ANGULAR.DAT"

    if os.path.exists(fresult):
        print('ANGULAR finish')
        fcc_x, fcc_y = np.loadtxt(fresult, usecols=(0, 5), unpack=True)

        # Perform interpolation on the CPU
        finterpolate = interpolate.interp1d(fcc_x, fcc_y, kind='cubic')
        fcc_exp = finterpolate(expx)
    else:
        raise FileNotFoundError("Result file not found.")

    # Clean up the temporary directory
    os.chdir('..')
    shutil.rmtree(fn)

    return fcc_exp


# Define ranges for the parameters
beta_2_range = np.arange(0.237, 0.2751, 0.001)
beta_4_range = np.arange(0.040, 0.0601, 0.001)
beta_3_range = np.arange(0.131, 0.1371, 0.001)


def calculate_chisq(beta_2, beta_4, beta_3):
    par = [beta_2, beta_4, beta_3]
    ftheo = model(expx, par)

    # Calculate chi-square
    chisq = np.sum((ftheo - expy) ** 2 / expyerr ** 2)  # CPU calculation (no need for .get())
    return {'beta_2': beta_2, 'beta_4': beta_4, 'beta_3': beta_3, 'chi-square': chisq}


# Create combinations of parameters
param_combinations = itertools.product(beta_2_range.tolist(), beta_4_range.tolist(), beta_3_range.tolist())

# Perform parallel calculations
num_cores = -1  # Use all available CPU cores
results = Parallel(n_jobs=num_cores)(
    delayed(calculate_chisq)(beta_2, beta_4, beta_3) for beta_2, beta_4, beta_3 in param_combinations
)

# Write results to CSV
with open(fchain, 'a', newline='') as csvfile:  # Ensure it is replaced with the desired file name
    fieldnames = ['beta_2', 'beta_4', 'beta_3', 'chi-square']
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

    if csvfile.tell() == 0:  # If file is empty, write header
        writer.writeheader()

    # Write all results at once
    for result in results:
        writer.writerow(result)