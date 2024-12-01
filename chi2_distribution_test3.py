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
data = pd.read_csv('filterdata07.csv')

Eeff = data['Sm154_Eeff'].values  # Convert to NumPy array
Ratio = data['Sm154_Ratio'].values  # Convert to NumPy array
Error = data['Sm154_Err'].values  # Convert to NumPy array

# Convert NumPy arrays to CuPy arrays for GPU processing
expx = cp.asarray(Eeff)
expy = cp.asarray(Ratio)
expyerr = cp.asarray(Error)

# Get the length of expx
expn = len(expx)

fchain = 'Sm154-beta3-sup1.csv'
if os.path.exists(fchain):  # if existed
       os.remove(fchain)


def model(expx, par):
    beta_200, beta_400, beta_300 = par[0], par[1], par[2]

    # 创建临时目录
    uuid_str = uuid.uuid4().hex
    tmp_file_name = f'tmpfile_{uuid_str}'
    fn = tmp_file_name
    # print(fn)
    if not os.path.exists(fn):
        os.makedirs(fn)
    # os.mkdir('Test')
    # 切换到文件夹下
    os.chdir(fn)
    # 复制输入文件
    shutil.copyfile('../ccfull-sc.inp', './ccfull-sc.inp')
    shutil.copyfile('../TEST2.INP', './TEST2.INP')
    shutil.copyfile('../a', './a')


    # 修改参数并写回文件
    with open("ccfull-sc.inp", 'r') as file:
        lines = file.readlines()
        lines[2] = f'0.082,{beta_200},{beta_400},3\n'
        lines[3] = f'1.012,{beta_300},3,1\n'

    with open("ccfull-sc.inp", 'w') as file:
        file.writelines(lines)

    # 执行程序
    sysstr = platform.system()
    main = "ccfull-sc2.exe" if sysstr == "Windows" else "./a"

    if os.path.exists(main):
        os.system(f'chmod 777 {main}' if sysstr == "Linux" else '')
        os.system(main)

    fresult = "ANGULAR.DAT"

    if os.path.exists(fresult):
        print('ANGULAR finish')
        fcc_x, fcc_y = cp.loadtxt(fresult, usecols=(0, 5), unpack=True)

        # 在 GPU 上执行插值
        finterpolate = interpolate.interp1d(cp.asnumpy(fcc_x), cp.asnumpy(fcc_y), kind='cubic')
        fcc_exp = cp.asarray(finterpolate(cp.asnumpy(expx)))
    else:
        raise FileNotFoundError("Result file not found.")

    # 清理临时文件夹
    os.chdir('..')
    shutil.rmtree(fn)

    return fcc_exp


# Cut range
beta_2_range = cp.arange(0.237, 0.2751, 0.001)
beta_4_range = cp.arange(0.040, 0.0601, 0.001)
beta_3_range = cp.arange(0.131, 0.1371, 0.001)



def calculate_chisq(beta_2, beta_4, beta_3):
    par = [beta_2, beta_4, beta_3]
    ftheo = model(expx, par)

    chisq = cp.sum((ftheo - expy) ** 2 / expyerr ** 2).get()  # 从 GPU 转回 CPU 以返回结果
    return {'beta_2': beta_2, 'beta_4': beta_4, 'beta_3': beta_3, 'chi-square': chisq}


# Create combinations of parameters
param_combinations = itertools.product(beta_2_range.tolist(), beta_4_range.tolist(), beta_3_range.tolist())

# Perform parallel calculations
num_cores = -1  # Use all available CPU cores
results = Parallel(n_jobs=num_cores)(
    delayed(calculate_chisq)(beta_2, beta_4, beta_3) for beta_2, beta_4, beta_3 in param_combinations
)

# Write results to CSV
with open(fchain, 'a', newline='') as csvfile:  # 确保替换为你想要的文件名
    fieldnames = ['beta_2', 'beta_4', 'beta_3', 'chi-square']
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

    if csvfile.tell() == 0:  # If file is empty, write header
        writer.writeheader()
        # Write all results at once
    for result in results:
        writer.writerow(result)