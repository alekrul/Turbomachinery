# Axial turbine 1D design
#Author: Alexandre Mendonça Krul (krul@ita.br)

from math import *
import numpy as np

#data:
m_dot = 20 #[kg/s]
p_in = 4 #[bar]
T01 = 1100 #[K]
p_ratio = 1.873 #[adim]
temp_drop = 145 #[K]

#assumptions
loss_coeff = 0.05
eff = 0.90

#constants
cp = 1148
gamma = 1.333
R = 0.287

KK = gamma/(gamma-1)

#constraints
N = 250 #[rev/s] --- compressor constraint
Um = 340 #[m/s] --- blade tip stress constraint

#DESIGN

#assuming (a) Va2 = Va3 and (b) V1 = V3 and (c) single stage -> alpha1 = 0
alpha1 = 0


flux_coeff = 0.8 #assumption
temp_drop_coeff = 2*cp*temp_drop/(Um**2)

convergence = 0.10
max_it = 10
desired_reaction = 0.5
for i in range(max_it):
    alpha3 = alpha1
    tg_beta3 = tan(alpha3) + 1/flux_coeff

    reaction = (tg_beta3*(2*flux_coeff) - (0.5*temp_drop_coeff))/2

    print(reaction)

    if((desired_reaction-reaction)>convergence):

        alpha1 = alpha1+radians(5)
    else:
        break
print("final reaction:",reaction)
print("swirl angle:", degrees(alpha1))

tg_beta2 = (1/(2*flux_coeff))*((temp_drop_coeff/2)-(2*reaction))
tg_alpha2 = tg_beta2 + 1/flux_coeff
beta2 = atan(tg_beta2)
alpha2 = atan(tg_alpha2)

print(degrees(beta2))
print(degrees(alpha2))

#with data above the velocities diagram can be sketched
#now we will estimate density at statios 1, 2 and 3 and blade height h and tip/root radius

#starting at station 2
Va2 = Um*flux_coeff
V2 = Va2/cos(alpha2)

delta_T02_T2 = V2**2/(2*cp)
T02 = T01
T2 = T02-delta_T02_T2

delta_T2_T2prime = loss_coeff*delta_T02_T2
T2_prime = T2-delta_T2_T2prime

p01_p2 = (T01/T2_prime)**KK
p2 = p_in/p01_p2

p01_pc = ((gamma+1)/2)**KK
print(p01_pc>p01_p2)

rho2 = p2*100/(R*T2)
A2 = m_dot/(rho2*Va2)
A2N = m_dot/(rho2*V2) #throat area of the Nozzle

print(rho2)
print(A2)
print(A2N)

#station 1
#Va1 = V1 = V3
Va1 = Va2/cos(alpha3)
V1 = Va1
delta_in_out = (V1**2)/(2*cp)

T1 = T01 - delta_in_out
p1_p01 = (T1/T01)**KK
p1 = p_in*p1_p01

rho1 = (100*p1)/(R*T1)
A1 = m_dot/(rho1*Va1)
print(rho1)
print(A1)


#section 3

T03 = T01-temp_drop
T3 = T03-delta_in_out
p3 = (p_in/p_ratio)*((T3/T03)**KK)
rho3 = (100*p3)/(R*T3)
A3 = m_dot/(rho3*Va2) #Va2 = Va3
print(rho3)
print(A3)

#estimate blade height and annulus adius ation at 1, 2 and 3
rm = Um/(2*pi*N)

h1 = A1*N/Um
h2 = A2*N/Um
h3 = A3*N/Um

rt_rr_1 = (rm+(h1/2))/(rm-(h1/2))
rt_rr_2 = (rm+(h2/2))/(rm-(h2/2))
rt_rr_3 = (rm+(h3/2))/(rm-(h3/2))

print("heights: ",h1," ",h2," ", h3)
print("rt_rr: ",rt_rr_1," ",rt_rr_2," ", rt_rr_3)
