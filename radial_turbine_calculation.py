# Exemplo de cálculo de turbina radial dada em sala de aula
#Disciplina ME211 - Turbomáquinas prof. Dr. Takachi
#Autor: Alexandre Mendonça Krul (krul@ita.br)

from math import *
import numpy as np

#dados:

Nv = 37
tte = 1.5e-3 #m
r0 = 0.181
r1 = 0.168
r2 = 0.166
r3t = 0.132
r3h = 0.022
b = 0.016
N = 48000 #rpm
T_t0 = 1280 #K
p_t0 = 8.6 #Pa
alpha_1 = radians(70) #radians
M_cri1 = 1
M_cri3 = 0.47
T_t3 = 800 #K
eff = 0.88
gamma = 1.33
cp = 1158
R = 287

omega = 48000*2*pi/60 #rad/s

print("------- PART A -------")

GAMA = (2*pi*r1 - Nv*(tte/cos(alpha_1)))/(2*pi*r1)

yy = (2*gamma)/(gamma+1)
xx = (gamma-1)/(gamma+1)
zz = 1/(gamma-1)

A1_act = GAMA*2*pi*r1*b*cos(alpha_1)

print("Area: ", A1_act," mˆ2")

T_t1 = T_t0
p_t1 = p_t0

aux1 = (p_t1*100000*A1_act)/sqrt(T_t1)
aux2 = sqrt(yy/R)*M_cri1
aux3 = (1 - (xx*(M_cri1**2)))**zz
m_dot = aux1*aux2*aux3
print("m_dot: ", m_dot," kg/s")

print("------- PART B -------")


V_cr1 = sqrt(R*yy*T_t1)
V1 = M_cri1*V_cr1
print("V1: ", V1," m/s")

V_t1 = V1*sin(alpha_1)
V_r1 = V1*cos(alpha_1)

V_r1_affected = GAMA*V_r1

M_cri1_affected = sqrt(V_t1**2 + V_r1_affected**2)/V1

print("M_cri1': ", M_cri1_affected)
print("Change in M_cri1': ", round(100*(M_cri1-M_cri1_affected)/M_cri1),"%")

print("------- PART C -------")

#p_t3 = p_t0*((((T_t3/T_t1)-1)/eff)+1)**(gamma/(gamma-1))
p_t3 = (((((T_t3/T_t0)-1)/eff)+1)**(gamma/(gamma-1)))*p_t0
p_t3_corrected = p_t3/GAMA
print("p_t3': ", p_t3_corrected," bar")

rho_t3 = p_t3_corrected*100000/(R*T_t3)
rho_3 = rho_t3*(1-xx*(M_cri3**2))**zz
V_cr3 = sqrt(R*yy*T_t3)

V3 = V_cr3*M_cri3

V_z3 = (m_dot/rho_3)/(pi*(r3t**2 - r3h**2))
print("V3: ",V3," m/s")

alpha_3 = -acos(V_z3/V3)
print("Alpha_3:",degrees(alpha_3)," °")

r3 = (r3h+r3t)/2

W_t3 = V3*sin(alpha_3)-(omega*r3)

W3 = sqrt(W_t3**2 + V_z3**2)
print("W3: ", W3," m/s")

betha_3 = atan(W_t3/V_z3)
print("Betha_3:", degrees(betha_3)," °")

print("------- PART D -------")

T_tr3 = T_t3 - ((V3**2 - W3**2)/(2*cp))
W_cr3 = sqrt(R*yy*T_tr3)


M_crr3 = W3/W_cr3
if M_crr3>=1:

    W3 = W_cr3
    print("Choked, W3: ", W3," m/s")

print("M_crr3: ", M_crr3)

print("------- PART E -------")

alpha_2 = alpha_1
M_cri2 = M_cri1_affected
V_cri2 = V_cr1
T_t2 = T_t1

V2 = M_cri2*V_cri2
W_r2 = V2*cos(alpha_2)
V_r2 = W_r2
V_t2 = V2*sin(alpha_2)
W_t2 = V_t2 - (omega*r2)
W2 = sqrt(W_t2**2 + W_r2**2)
T_tr2 = T_t2 - ((V2**2 - W2**2)/(2*cp))
deltaTr_rotor = T_tr3 - T_tr2

print("V2: ", V2," m/s")
print("W2: ", W2," m/s")
print("Delta relative T in rotor: ", deltaTr_rotor," K")

print("------- PART F -------")

deltah_id = cp*T_tr3*(1 - (T_t3/T_tr3))/eff
Ns = omega*sqrt(m_dot/rho_3)/deltah_id**(3/4)
print("Specific speed:", Ns)


print("------- PART G -------")
U2 = omega*r2
U3 = omega*r3

reaction = ((W3**2 - W2**2)+(U2**2 - U3**2))/((W3**2 - W2**2)+(U2**2 - U3**2)+(V2**2 - V3**2))
print("Stage Reaction: ", reaction*100, "%")

print("------- PART H -------")

work_coeff = cp*(T_t0-T_t3)/(U2**2)
print("Stage work coefficient: ",work_coeff)



