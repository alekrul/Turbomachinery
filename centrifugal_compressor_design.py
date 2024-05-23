from math import *
import numpy as np
from slip_factor import slipFactor

#input
psi = 1.04 #power input factor
N = 290 #rev/s
d_imp = 0.5 #m
d_et = 0.3 #m eye tip diameter
d_er = 0.15 #m eye root diameter
m_dot = 9 #kg/s
T_t1 = 295 #K
p_t1 = 1.1 #bar
n_iso = 0.78 #isoentropic efficiency

R = 287
cp = 1005
gamma = 1.4

U = pi*d_imp*N #impeller tip speed

#find V1 = Va1 (Vt1 = 0)

A1 = pi*(d_et**2 - d_er**2)/4

rho1 = p_t1*100000/R*T_t1
Va1_ini = m_dot/(rho1*A1) #primeira aproximação de Va1

max_iter = 100
er = 0.02
yy = gamma/(gamma - 1)

for i in range(max_iter):
    T_1 = T_t1 - (Va1_ini**2)/(2*cp)
    p_1 = p_t1/((T_t1/T_1)**yy)
    rho1 = p_1*100000/(R*T_1)

    #checando
    Va1 = m_dot/(rho1*A1)

    error = abs(Va1 - Va1_ini)/Va1_ini

    Va1_ini = Va1

    if error <= er:
        break

U1_et = pi*d_et*N
U1_er = pi*d_er*N

beta1_r = atan(Va1/U1_er)
beta1_t = atan(Va1/U1_et)

print(degrees(beta1_r))
print(degrees(beta1_t))

