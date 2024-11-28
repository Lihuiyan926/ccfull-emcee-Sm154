
import subprocess
import os
import numpy as np
import matplotlib
# matplotlib.use('Agg') #调用后端
import matplotlib.pyplot as plt
# import sys
import platform
from iminuit import Minuit
from iminuit.util import propagate
# from iminuit.cost import LeastSquares
from scipy import interpolate #插值操作
import pandas as pd
import uuid
import shutil
# import time
from joblib import Parallel, delayed
from scipy.optimize import differential_evolution

# from multiprocessing import cpu_count

#os.environ["OMP_NUM_THREADS"] = "1"
par=[0, 0,0] # 参数初始设置【

#add the experimental data
fit0 = 'fit1-Sm154.out'
if os.path.exists(fit0):  # 如果文件存在
    # 删除文件，可使用以下两种方法。
    os.remove(fit0)
    #os.unlink(path)

#add the experimental data
data = pd.read_csv('filterdata07.csv')
Eeff = data['Sm154_Eeff']
Ratio = data['Sm154_Ratio']
# Dqel = data['Sm152_Dqel']
# Derror = data['Sm152_Derror']
Error = data['Sm154_Err']



expx=Eeff
expy=Ratio
# expy=Dqel
expyerr=Error
expn=len(expx)


def model(expx, par):
    beta_200, beta_400,beta_300 = par[0], par[1],par[2]

    uuid_str = uuid.uuid4().hex
    tmp_file_name = 'tmpfile_%s' % uuid_str
    fn = tmp_file_name
    # print(fn)
    if not os.path.exists(fn):
        os.makedirs(fn)
    # os.mkdir('Test')
    # 切换到文件夹下
    os.chdir(fn)

    shutil.copyfile('../ccfull-sc.inp', './ccfull-sc.inp')
    shutil.copyfile('../TEST2.INP', './TEST2.INP')
    shutil.copyfile('../a', './a')

    with open("ccfull-sc.inp", "r") as file:
        lines = file.readlines()

        # 修改参数
    linereplace1 = 3
    linerplace2 = 4
    beta_20 = beta_200
    beta_40 = beta_400
    beta_30 = beta_300
    # 修改特定行的内容
    lines[linereplace1 - 1] = '0.082,' + str(beta_20) + ',' + str(beta_40) + ',5\n'
    lines[linerplace2 - 1] = '1.012,' + str(beta_30) + ',3,1\n'
    # 将修改后的内容写回文件
    with open("ccfull-sc.inp", "w") as file:
        file.writelines(lines)
    # print(theta)

    sysstr = platform.system()
    if (sysstr == "Windows"):
        main = "ccfull-sc2.exe"
        # main2 = "derivation.exe"
    elif (sysstr == "Linux"):
        main = "./a"
        # main2 = "./derivation"
        os.system('chmod 777 a')
        # os.system('chmod 777 derivation')

    if os.path.exists(main):
        os.system(main)

    # if os.path.exists(main) and os.path.exists(main2):
    #    os.system(main)
    #    os.system(main2)

    fresult = "ANGULAR.DAT"
    if os.path.exists(fresult):
        print('ANGULAR finish')
    # fresult = "db2.out"
    global fcc_x, fcc_y
    # AngDis 读取
    fcc_x, fcc_y = np.loadtxt(fresult, usecols=(0, 5), unpack=True)
    # Dqel 读取
    # fcc_x, fcc_y = np.genfromtxt(fresult, usecols=(0, 1), skip_header=1,skip_footer=1, unpack=True)

    finterpolate = interpolate.interp1d(fcc_x, fcc_y, kind='cubic')  # 三次样条插值
    fcc_exp = finterpolate(expx)


    # fresult2="db2.out"
    # # Dqel 读取
    # fcc_x1, fcc_y1 = np.genfromtxt(fresult2, usecols=(0, 1), skip_header=1,skip_footer=1, unpack=True)
    # finterpolate = interpolate.interp1d(fcc_x1, fcc_y1, kind='cubic')  # 三次样条插值
    # fcc_exp1 = finterpolate(expx)
    os.chdir('..')
    shutil.rmtree(fn)

    return fcc_exp

def sfactor_chi(beta_20,beta_40,beta_30):
    # to mimic the result of ccfull calculation
    par[0], par[1],par[2] = beta_20,beta_40,beta_30
    ftheo=model(expx,par)
    fcn =np.sum( (ftheo-expy)**2 /expyerr**2)
    tmp = np.hstack(([beta_20, beta_40,beta_30, fcn], ftheo.flatten()))
    fit = open(fit0, 'a')
    np.savetxt(fit, tmp.reshape(1,-1), fmt='%.5e',  header='', newline='\n', )
    # 计算一阶微分的均方差
    # diff_loss = np.mean((fdqel - expy1) ** 2)
    fit.close()
    return fcn

def pre_chi(par):
    # to mimic the result of ccfull calculation
    ftheo =model(expx,par)
    fcn =np.sum( (ftheo-expy)**2 /expyerr**2)
    print(fcn)
    # 计算一阶微分的均方差
    # diff_loss = np.mean((fdqel - expy1) ** 2)
    return fcn
# 定义拟合参数的边界
bounds = [(0.0, 0.4), (0, 0.2),(0.0,0.2)]  # 对应 beta_20 和 beta_40 的边界

# 使用全局优化算法获取最优的初始参数
initial_params = differential_evolution(pre_chi, bounds,workers=10).x

print("DE finished")
m = Minuit(sfactor_chi, beta_20=initial_params[0], beta_40=initial_params[1],beta_30=initial_params[2])
m.limits['beta_20'] =(0.1,0.4)
m.limits['beta_40'] =(0,0.2)
m.limits['beta_30'] = (0.0,0.2)
m.migrad()
# m.minos()
print("Minuit finished")

# m.draw_mnmatrix()
# plt.show()
# print(m.values)
# print(m.errors)
# print(m.fval)
# print(fit_info)
# print(m.fval / len(x1))
# beta_20,beta_40= m.values
# par[0], par[1] = beta_20,beta_40
# fcc_exp = model(expx,par)

print("Global best solution:", initial_params)
print(m.values)
print(m.errors)
print(m.covariance)
print(m.fval)
# run error propagation
# x1= np.arange(50, 72, 1.)
# y, ycov = propagate(lambda p: model(x1, p), m.values, m.covariance)
# print(y)
# print(ycov)
# plot everything
#plt.errorbar(expx, d, derr, fmt="o", label="data", zorder=0)