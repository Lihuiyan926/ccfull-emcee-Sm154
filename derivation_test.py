import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import math

# Ap=16
# AT=186

# theta_lab=175
# theta_l=math.radians(theta_lab)
# theta=theta_l+math.asin(Ap*math.sin(theta_l)/AT)

fresult = "ANGULAR.DAT"
Eeff, R=np.loadtxt(fresult, usecols=(0, 5), unpack=True)

# data = pd.read_csv('data186.csv')
# Elab = data['W186_Elab']
# R = data['W186_Ratio']
# Error=data['W186_Err']
#
# Ecm = AT/(Ap+AT)*Elab
#
# Eeff=2*Ecm*math.sin(theta/2)/(1+math.sin(theta/2))

n=len(Eeff)

Stp=2
accur=0.05

with open('derivation_8.dat', 'w') as f:
    f.write("Eeff   Derivative \n")  # 添加列名
    for i in range(0, n):
        l = 0

        for j in range(i, n):
            diff = abs(Eeff[j] - Eeff[i])

            if (diff <= (Stp + accur) and diff >= (Stp - accur)):
                l = j

            if (l != 0  and i != l ):
                e1 = Eeff[i]
                e3 = Eeff[l]
                e2=(e3+e1)/2
                s1 = R[i]
                s3 = R[l]
                # err1 = Error[i]
                # err3 = Error[l]

                # barrier distribution
                # three point formula
                deriv = (s3 - s1) / (e3 - e1)
                # Calculate deriv error
                # deriv_error = math.sqrt(((1 / (e3 - e1)) * err3) ** 2  + ((-1 / (e3 - e1)) * err1) ** 2)

                f.write('{:.2f} {:.4E} \n'.format(e2, -deriv))


                l = 0