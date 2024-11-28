import numpy as np
import emcee
import shutil
import os
import platform
import uuid
import sys
import pandas as pd
import corner
from matplotlib import pyplot as plt
from scipy import interpolate
from multiprocessing import Pool



# 定义模型函数
def model(expx, par):
    beta2, beta3, beta4 = par[0], par[1], par[2]

    # 创建临时目录
    uuid_str = uuid.uuid4().hex
    tmp_file_name = f'tmpfile_{uuid_str}'
    fn = tmp_file_name
    if not os.path.exists(fn):
        os.makedirs(fn)

    os.chdir(fn)
    shutil.copyfile('../ccfull-sc.inp', './ccfull-sc.inp')
    shutil.copyfile('../TEST2.INP', './TEST2.INP')
    shutil.copyfile('../a', './a')

    # 修改参数并写回文件
    with open("ccfull-sc.inp", 'r') as file:
        lines = file.readlines()
        lines[2] = f'0.082,{beta2},{beta4},5\n'
        lines[3] = f'1.012,{beta3},3,1\n'

    with open("ccfull-sc.inp", 'w') as file:
        file.writelines(lines)

    # 执行模型
    sysstr = platform.system()
    main = "ccfull-sc2.exe" if sysstr == "Windows" else "./a"
    if os.path.exists(main):
        os.system(f'chmod 777 {main}' if sysstr == "Linux" else '')
        os.system(main)

    # 加载结果
    fresult = "ANGULAR.DAT"
    if os.path.exists(fresult):
        print('ANGULAR finish')
        fcc_x, fcc_y = np.loadtxt(fresult, usecols=(0, 5), unpack=True)

        # 使用CPU进行插值
        finterpolate = interpolate.interp1d(fcc_x, fcc_y, kind='cubic')
        fcc_exp = finterpolate(expx)
    else:
        raise FileNotFoundError("Result file not found.")

    # 清理临时文件
    os.chdir('..')
    shutil.rmtree(fn)

    return fcc_exp


# 定义对数似然函数
def log_likelihood(par, expx, expy, expyerr):
    model_values = model(expx, par)

    # 计算残差（模型 - 数据）以及加权的残差平方和
    residual = model_values - expy
    chi2 = np.sum((residual / expyerr) ** 2)

    # 似然是负的卡方值
    return -0.5 * chi2


# 定义对数先验函数（这里是均匀先验）
def log_prior(par):
    beta2, beta3, beta4 = par
    # 使用指定范围内的均匀先验
    if 0.0< beta2 < 0.4 and 0.0 < beta3 < 0.2 and 0.0 < beta4 < 0.2:
        return 0.0  # log(1)，表示平坦先验
    return -np.inf  # log(0)，表示超出先验范围


# 定义对数后验函数（先验 * 似然）
def log_posterior(par, expx, expy, expyerr):
    lp = log_prior(par)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(par, expx, expy, expyerr)


# 定义单个 MCMC 步骤的并行计算函数
def parallel_log_posterior(params, expx, expy, expyerr):
    return log_posterior(params, expx, expy, expyerr)


# 设置emcee采样器并行化
def run_mcmc(expx, expy, expyerr, nwalkers=128, nsteps=500, burn_in_steps=50):
    # 初始参数猜测（可以调整）
    initial_pos = [np.array([0.26, 0.07, 0.06]) + 1e-3 * np.random.randn(3) for i in range(nwalkers)]

    # 设置MCMC采样器
    ndim = 3  # 参数数量
    ncores = os.cpu_count()  # 获取可用的CPU核心数
    with Pool(processes=ncores) as pool:
        # 创建并设置移动策略为 StretchMove，步长设置为 0.001
        sampler = emcee.EnsembleSampler(
            nwalkers, ndim, log_posterior, args=(expx, expy, expyerr), pool=pool
        )

        # 运行MCMC链
        sampler.run_mcmc(initial_pos, burn_in_steps, progress=True)

        # 丢弃burn-in阶段的样本，保留后续的样本
        sampler.reset()  # 重置以开始正式的采样
        sampler.run_mcmc(initial_pos, nsteps, progress=True)

    # 获取采样结果
    samples = sampler.get_chain(flat=True)

    # 将样本保存到txt文件中
    np.savetxt('mcmc_samples_Sm154.txt', samples, fmt='%f', delimiter=',', header="beta2, beta3, beta4")

    return samples


# 加载数据
data = pd.read_csv('filterdata07.csv')
Eeff = data['Sm154_Eeff'].values  # 转换为NumPy数组
Ratio = data['Sm154_Ratio'].values  # 转换为NumPy数组
Error = data['Sm154_Err'].values  # 转换为NumPy数组

# 运行MCMC以拟合模型参数
samples = run_mcmc(Eeff, Ratio, Error)

# 打印最优拟合参数
best_params = np.median(samples, axis=0)
print(f"Best-fit parameters: beta2={best_params[0]}, beta3={best_params[1]}, beta4={best_params[2]}")


# 绘制corner图，显示95%置信区间
fig = corner.corner(samples,
                    labels=["beta2", "beta3", "beta4"],
                    # truths=best_params,
                    quantiles=[0.025,0.5, 0.975],  # 95%置信区间
                    show_titles=True,
                    title_fmt='.3f',
                    title_args={"fontsize": 12},
                    fill_contours=True,
                    )  # 设置标题字体大小
# 显示图形
plt.savefig('Sm154-emcee.png')
plt.savefig('Sm154-emcee.eps')
plt.show()

# 提取每个参数的样本
beta_2_samples = samples[:, 0]
beta_3_samples = samples[:, 1]
beta_4_samples = samples[:, 2]

# 计算皮尔逊相关系数矩阵
correlation_matrix = np.corrcoef([beta_2_samples, beta_3_samples, beta_4_samples])

# 打印相关系数矩阵
print("beta_2, beta_3, beta_4 之间的相关系数矩阵：")
print(correlation_matrix)

# 将相关系数矩阵写入txt文件
with open("correlation_matrix_Sm154.txt", "w") as f:
    f.write("beta_2, beta_3, beta_4 之间的相关系数矩阵：\n")
    np.savetxt(f, correlation_matrix, fmt="%.4f", delimiter="\t")