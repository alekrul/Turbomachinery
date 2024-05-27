# Exemplo de cálculo de turbina radial dada em sala de aula
#Disciplina ME211 - Turbomáquinas prof. Dr. Takachi
#Autor: Alexandre Mendonça Krul (krul@ita.br)

from math import *
import numpy as np

#Given Values:

Nv = 37 # Number of vanes
tte = 1.5e-3 # Trailing edge thickness in meters
r0 = 0.181 #m
r1 = 0.168 # Radius at stator exit in meters
r2 = 0.166 # Radius at rotor inlet in meters
r3t = 0.132 # Rotor exit radius in meters
r3h = 0.022 #m
b = 0.016 # Blade height in meters
N = 48000 # Rotational speed in rpm
T_t0 = 1280 # Total temperature at stator exit in K
p_t0 = 8.6 # Total pressure at stator exit in Pa in bars
alpha_1 = radians(70) # Swirl flow angle in radians
M_cri1 = 1 # Critical Mach number at stator exit
M_cri3 = 0.47 # Critical Mach number at rotor exit
T_t3 = 800 # Total temperature at rotor exit in K
eff = 0.88

#Constants
gamma = 1.33 # Specific heat ratio (Isentropic flow)
cp = 1158 # Specific heat at constant pressure
R = 287 # Specific gas constant for air in J/(kg*K)

omega = 48000*2*pi/60 #rad/s

print("------- PART A -------")

GAMA = (2*pi*r1 - Nv*(tte/cos(alpha_1)))/(2*pi*r1)

yy = (2*gamma)/(gamma+1)
xx = (gamma-1)/(gamma+1)
zz = 1/(gamma-1)

A1_act = GAMA*2*pi*r1*b*cos(alpha_1)

print(f"Area: {A1_act:.5f} mˆ2")

T_t1 = T_t0
p_t1 = p_t0

aux1 = (p_t1*100000*A1_act)/sqrt(T_t1)
aux2 = sqrt(yy/R)*M_cri1
aux3 = (1 - (xx*(M_cri1**2)))**zz
m_dot = aux1*aux2*aux3

print(f"a) Mass flow rate at the stator exit: {m_dot:.2f} kg/s")

print("------- PART B -------")


V_cr1 = sqrt(R*yy*T_t1)
V1 = M_cri1*V_cr1
print(f"V1: {V1:.2f} m/s")

V_t1 = V1*sin(alpha_1)
V_r1 = V1*cos(alpha_1)

#adjust radial component for expansion
V_r1_affected = GAMA*V_r1

M_cri1_affected = sqrt(V_t1**2 + V_r1_affected**2)/V1

print(f"b) New Mach number just downstream (Mcr1'): {M_cri1_affected:.2f}")

print(f"Change in M_cri1': {round(100*(M_cri1-M_cri1_affected)/M_cri1)}%")

print("------- PART C -------")

#p_t3 = p_t0*((((T_t3/T_t1)-1)/eff)+1)**(gamma/(gamma-1))
p_t3 = (((((T_t3/T_t0)-1)/eff)+1)**(gamma/(gamma-1)))*p_t0
p_t3_corrected = p_t3/GAMA
print(f"p_t3': {p_t3_corrected:.3f} bar")

rho_t3 = p_t3_corrected*100000/(R*T_t3)
rho_3 = rho_t3*(1-xx*(M_cri3**2))**zz
V_cr3 = sqrt(R*yy*T_t3)

V3 = V_cr3*M_cri3

V_z3 = (m_dot/rho_3)/(pi*(r3t**2 - r3h**2))
print(f"V3: {V3:.2f} m/s")

alpha_3 = -acos(V_z3/V3)
print(f"c) Rotor-exit absolute flow angle (α3m): {degrees(alpha_3):.2f} degrees")

r3 = (r3h+r3t)/2

W_t3 = V3*sin(alpha_3)-(omega*r3)

W3 = sqrt(W_t3**2 + V_z3**2)
print(f"W3: {W3:.2f} m/s")

betha_3 = atan(W_t3/V_z3)

print(f"   Rotor-exit relative flow angle (β3m): {degrees(betha_3):.2f} degrees")

print("------- PART D -------")

T_tr3 = T_t3 - ((V3**2 - W3**2)/(2*cp))
W_cr3 = sqrt(R*yy*T_tr3)


M_crr3 = W3/W_cr3

is_choked = M_crr3 >= 1

print(f"d) The rotor is choked: {is_choked}")

if M_crr3>=1:
    W3 = W_cr3

print(f"M_crr3: {M_crr3:.3f}")

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

print(f"V2: {V2:.2f} m/s")
print(f"W2: {W2:.2f} m/s")
print(f"e) Change in total relative temperature (ΔTtr): {deltaTr_rotor:.2f} K")


print("------- PART F -------")

deltah_id = cp*T_tr3*(1 - (T_t3/T_tr3))/eff
Ns = omega*sqrt(m_dot/rho_3)/deltah_id**(3/4)

print(f"f) Stage specific speed (Ns): {Ns:.2f} radians")


print("------- PART G -------")
U2 = omega*r2
U3 = omega*r3

reaction = ((W3**2 - W2**2)+(U2**2 - U3**2))/((W3**2 - W2**2)+(U2**2 - U3**2)+(V2**2 - V3**2))


print(f"g) Stage reaction (R): {reaction * 100:.2f}%")


print("------- PART H -------")

work_coeff = cp*(T_t0-T_t3)/(U2**2)


print(f"h) Work coefficient (ψ): {work_coeff:.3f}")


