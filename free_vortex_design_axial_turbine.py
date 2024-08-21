# Axial turbine 1D design using free vortex design
#Author: Alexandre Mendonça Krul (krul@ita.br)

from math import *
import numpy as np
import axial_turbine_1D_design as turbine
import matplotlib.pyplot as plt

#rr = rm - (h/2) and rt = rm + (h/2)

#inicio:

rr2 = turbine.rm - turbine.h2/2
rt2 = turbine.rm + turbine.h2/2

rr3 = turbine.rm - turbine.h3/2
rt3 = turbine.rm + turbine.h3/2

rm_rr_2 = turbine.rm/rr2
rm_rt_2 = turbine.rm/rt2

rm_rr_3 = turbine.rm/rr3
rm_rt_3 = turbine.rm/rt3

r2 = np.linspace(rr2, rt2, num=20)
r3 = np.linspace(rr3, rt3, num=20)

flux_coeff = turbine.flux_coeff

tg_alpha2 = np.zeros(20)
tg_alpha3 = np.zeros(20)
tg_beta2 = np.zeros(20)
tg_beta3 = np.zeros(20)
alpha2 = np.zeros(20)
alpha3 = np.zeros(20)
beta2 = np.zeros(20)
beta3 = np.zeros(20)

for i in range(len(r2)):

    tg_alpha2[i] = (turbine.rm/r2[i])*tan(turbine.alpha2)
    tg_alpha3[i] = (turbine.rm/r3[i])*tan(turbine.alpha3)
    tg_beta2[i] = tg_alpha2[i] - ((r2[i]/turbine.rm)*(turbine.Um/turbine.Va2))
    tg_beta3[i] = tg_alpha3[i] + ((r3[i]/turbine.rm)*(turbine.Um/turbine.Va2))
    alpha2[i] = degrees(atan(tg_alpha2[i]))
    alpha3[i] = degrees(atan(tg_alpha3[i]))
    beta2[i] = degrees(atan(tg_beta2[i]))
    beta3[i] = degrees(atan(tg_beta3[i]))



# Plotting all functions on the same graph
plt.plot(r2, alpha2, label=r'$\alpha_2$')
plt.plot(r2, beta2, label=r'$\beta_2$')
plt.plot(r2, alpha3, label=r'$\alpha_3$')
plt.plot(r2, beta3, label=r'$\beta_3$')

# Labeling the graph
plt.xlabel('r')
plt.ylabel('degrees')
plt.title('Plot of blade angles vs r')
plt.legend()
plt.grid(True)

# Show the plot
plt.show()