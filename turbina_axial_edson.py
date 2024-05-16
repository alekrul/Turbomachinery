import numpy as np
from math import *

# %% Turbina Axial
# Dados
m = 8.2  # kg/s
p_tin = 10 * 10**5  # Pa
T_tin = 1380  # K
p_tin_p_tex = 9.3
N = 35000 * 2 * np.pi / 60  # rad/s
R = 287  # J/kgK
cp = 1148  # J/kgK
Um = 414  # m/s
gama = 1.333
neta_tt = 0.91
neta_stg1 = 0.55
neta_stg2 = 0.45

# %% 1 Estágio
T_tex = T_tin * (1 + neta_tt * ((1 / p_tin_p_tex) ** ((gama - 1) / gama) - 1))
p_tex = p_tin * (1 / p_tin_p_tex)
rho_tex = p_tex / (R * T_tex)
delta_h_id = cp * T_tin * (1 - (1 / p_tin_p_tex) ** ((gama - 1) / gama))
Ns = (N * np.sqrt(m / rho_tex)) / (delta_h_id**(3 / 4))

# %% 2 Estágios
# Estágio 1
delta_h_id_stg1 = neta_stg1 * delta_h_id
p_t_int = p_tin * (1 - delta_h_id_stg1 / (cp * T_tin))**(gama / (gama - 1))
T_t_int = T_tin * (1 - neta_tt * (1 - (p_t_int / p_tin)**((gama - 1) / gama)))
rho_t_int = p_t_int / (R * T_t_int)
N_stg1 = N * np.sqrt(m / rho_t_int) / (delta_h_id_stg1)**(3 / 4)

# Estágio 2
delta_h_id_stg2 = neta_stg2 * delta_h_id
p_t_ex = p_t_int * (1 - delta_h_id_stg2 / (cp * T_t_int))**(gama / (gama - 1))
T_t_ex = T_t_int * (1 - neta_tt * (1 - (p_t_ex / p_t_int)**((gama - 1) / gama)))
rho_t_ex = p_t_ex / (R * T_t_ex)
N_stg2 = N * np.sqrt(m / rho_t_ex) / (delta_h_id_stg2)**(3 / 4)

# %% 3 Estágios
# Pressao por estagio cte
p_stg = (p_tin_p_tex)**(1 / 3)

# Analise termodinamica da Turbina
T_t0 = T_tin
p_t0 = p_tin

T_t2 = T_t0 * (1 - neta_tt * (1 - (1 / p_stg)**((gama - 1) / gama)))
p_t2 = p_t0 * 1 / p_stg
rho_t2 = p_t2 / (R * T_t2)

T_t4 = T_t2 * (1 - neta_tt * (1 - (1 / p_stg)**((gama - 1) / gama)))
p_t4 = p_t2 * 1 / p_stg
rho_t4 = p_t4 / (R * T_t4)

T_t6 = T_t4 * (1 - neta_tt * (1 - (1 / p_stg)**((gama - 1) / gama)))
p_t6 = p_t4 * 1 / p_stg
rho_t6 = p_t6 / (R * T_t6)

# Trabalho por estágio
w_t1 = int(100 * (T_t0 - T_t2) / (T_t0 - T_t6))
w_t2 = int(100 * (T_t2 - T_t4) / (T_t0 - T_t6))
w_t3 = 100 - (w_t1 + w_t2)

# Calculo do raio médio
rm = Um / N

# Ângulos de escoamento
# Escolhendo
alpha2 = -18  # degrees
alpha2_rad = np.radians(alpha2)  # rad
Phi = 0.65

# Calculo de betha 2
betha2_rad = np.arctan(np.tan(alpha2_rad) - 1 / Phi)
betha2 = np.degrees(betha2_rad)

# Calculo do coeficiente carregamento da blade Psi (estágio 1)
Psi_stg1 = cp * (T_t0 - T_t2) / (Um**2)

# Calculo do Grau de Reacao
R1 = -Phi * np.tan(betha2_rad) - 0.5 * Psi_stg1

# Calculo do angulo de entrada
betha1_rad = np.arctan(1 / (2 * Phi) * (Psi_stg1 - 2 * R1))
betha1 = np.degrees(betha1_rad)

alpha1_rad = 1 / Phi - np.tan(betha1_rad)
alpha1 = np.degrees(alpha1_rad)

# Estagio 2
# rm = cte; Phi = cte
# Escolhendo
alpha4 = -15  # degrees
alpha4_rad = np.radians(alpha4)  # rad

# Calculo de betha 4
betha4_rad = np.arctan(np.tan(alpha4_rad) - 1 / Phi)
betha4 = np.degrees(betha4_rad)

# Calculo do coeficiente carregamento da blade Psi (estágio 2)
Psi_stg2 = cp * (T_t2 - T_t4) / (Um**2)

# Calculo do Grau de Reacao
R2 = -Phi * np.tan(betha4_rad) - 0.5 * Psi_stg2

# Calculo do angulo de entrada
betha3_rad = np.arctan(1 / (2 * Phi) * (Psi_stg2 - 2 * R2))
betha3 = np.degrees(betha3_rad)

# Estagio 3
# rm = cte; Phi = cte
# Escolhendo
alpha6 = 0  # degrees
alpha6_rad = np.radians(alpha6)  # rad

# Calculo de betha 6
betha6_rad = np.arctan(np.tan(alpha6_rad) - 1 / Phi)
betha6 = np.degrees(betha6_rad)

# Calculo do coeficiente carregamento da blade Psi (estágio 3)
Psi_stg3 = cp * (T_t4 - T_t6) / (Um**2)

# Calculo do Grau de Reacao
R3 = -Phi * np.tan(betha6_rad) - 0.5 * Psi_stg3

# Calculo do angulo de entrada
betha5_rad = np.arctan(1 / (2 * Phi) * (Psi_stg3 - 2 * R3))
betha5 = np.degrees(betha5_rad)

# %% Calculo da velocidade estagio 1 - entrada
Ca_1 = Phi * Um
Wx_1 = Ca_1 * np.tan(betha1_rad)
Cx_1 = Um + Wx_1
C1 = np.sqrt(Ca_1**2 + Cx_1**2)
T_t1 = T_t0
# velocidade critica
V_cr1 = np.sqrt(R * T_t1 * 2 * gama / (gama + 1))
# Calculo da velocidade estagio 1 - saida
W_2 = Ca_1 / np.cos(betha2_rad)
C2 = Ca_1 / np.cos(alpha2_rad)
# Temperatura critica no rotor
T_tr2 = T_t1 - (W_2**2 + C2**2) / (2 * cp)
# velocidade critica no rotor
W_cr2 = np.sqrt(R * T_tr2 * 2 * gama / (gama + 1))
# Mach critico
M_cr1 = C1 / V_cr1
Mr_cr2 = W_2 / W_cr2

# %% Calculo da velocidade estagio 2 - entrada
Ca_3 = Phi * Um
Wx_3 = Ca_3 * np.tan(betha3_rad)
Cx_3 = Um + Wx_3
C3 = np.sqrt(Ca_3**2 + Cx_3**2)
T_t3 = T_t2
# velocidade critica
V_cr3 = np.sqrt(R * T_t3 * 2 * gama / (gama + 1))
# Calculo da velocidade estagio 2 - saida
W_4 = Ca_3 / np.cos(betha4_rad)
C4 = Ca_3 / np.cos(alpha4_rad)
# Temperatura critica no rotor
T_tr4 = T_t3 - (W_4**2 + C4**2) / (2 * cp)
# velocidade critica no rotor
W_cr4 = np.sqrt(R * T_tr4 * 2 * gama / (gama + 1))
# Mach critico
M_cr3 = C3 / V_cr3
Mr_cr4 = W_4 / W_cr4

# %% Calculo da velocidade estagio 3 - entrada
Ca_5 = Phi * Um
Wx_5 = Ca_5 * np.tan(betha5_rad)
Cx_5 = Um + Wx_5
C5 = np.sqrt(Ca_5**2 + Cx_5**2)
T_t5 = T_t4
# velocidade critica
V_cr5 = np.sqrt(R * T_t5 * 2 * gama / (gama + 1))
# Calculo da velocidade estagio 3 - saida
W_6 = Ca_5 / np.cos(betha6_rad)
C6 = Ca_5 / np.cos(alpha6_rad)
# Temperatura critica no rotor
T_tr6 = T_t5 - (W_6**2 + C6**2) / (2 * cp)
# velocidade critica no rotor
W_cr6 = np.sqrt(R * T_tr6 * 2 * gama / (gama + 1))
# Mach critico
M_cr5 = C5 / V_cr5
Mr_cr6 = W_6 / W_cr6

# %% Calculo do escoamento Meridional
# Estagio 1
perda = 5  # %
p_t1 = (1 - perda / 100) * p_t0
rho_t1 = p_t1 / (R * T_t1)
rho_1 = rho_t1 * (1 - (gama - 1) / (gama + 1) * M_cr1**2)**(1 / (gama - 1))
h1 = m / (rho_1 * 2 * np.pi * rm * Ca_1)

Ca_0 = Ca_1
V_cr0 = V_cr1
M_cr0 = Ca_0 / V_cr0
rho_t0 = p_t0 / (R * T_t0)
rho_0 = rho_t0 * (1 - (gama - 1) / (gama + 1) * M_cr0**2)**(1 / (gama - 1))
h0 = m / (rho_0 * 2 * np.pi * rm * Ca_0)

V_cr2 = np.sqrt(R * T_t2 * 2 * gama / (gama + 1))
M_cr2 = C2 / V_cr2
rho_t2 = p_t2 / (R * T_t2)
rho_2 = rho_t2 * (1 - (gama - 1) / (gama + 1) * M_cr2**2)**(1 / (gama - 1))
h2 = m / (rho_2 * 2 * np.pi * rm * Ca_1)

# Estagio 2
p_t3 = (1 - perda / 100) * p_t2
rho_t3 = p_t3 / (R * T_t3)
rho_3 = rho_t3 * (1 - (gama - 1) / (gama + 1) * M_cr3**2)**(1 / (gama - 1))
h3 = m / (rho_3 * 2 * np.pi * rm * Ca_3)

V_cr4 = np.sqrt(R * T_t4 * 2 * gama / (gama + 1))
M_cr4 = C4 / V_cr4
rho_t4 = p_t4 / (R * T_t4)
rho_4 = rho_t4 * (1 - (gama - 1) / (gama + 1) * M_cr4**2)**(1 / (gama - 1))
h4 = m / (rho_4 * 2 * np.pi * rm * Ca_3)

# Estagio 3
p_t5 = (1 - perda / 100) * p_t4
rho_t5 = p_t5 / (R * T_t5)
rho_5 = rho_t5 * (1 - (gama - 1) / (gama + 1) * M_cr5**2)**(1 / (gama - 1))
h5 = m / (rho_5 * 2 * np.pi * rm * Ca_5)

V_cr6 = np.sqrt(R * T_t6 * 2 * gama / (gama + 1))
M_cr6 = C6 / V_cr6
rho_t6 = p_t6 / (R * T_t6)
rho_6 = rho_t6 * (1 - (gama - 1) / (gama + 1) * M_cr6**2)**(1 / (gama - 1))
h6 = m / (rho_6 * 2 * np.pi * rm * Ca_5)

# Gap lengths (distribuição adotada)
Lambida = [1.5, 1.8, 2.1, 2.4, 2.7, 3]

# corda média
h = [h0, h1, h2, h3, h4, h5, h6]
hav = np.zeros(len(h) - 1)
for i in range(6):
    hav[i] = (h[i] + h[i + 1]) * 0.5

cz = hav / Lambida

# Distancia Axial
w = np.zeros(len(cz) - 1)
d = np.zeros(len(cz) - 1)
for i in range(5):
    w[i] = (hav[i] + hav[i + 1]) * 0.5
    if i % 2 != 0:
        d[i] = 0.25 * w[i]
    else:
        d[i] = 0.5 * w[i]

d = d * 100

# Velocidade especifica por estagio
delta_h_id_stg1_3 = (w_t1 / 100) * delta_h_id
N_stg1_3 = N * np.sqrt(m / rho_2) / (delta_h_id_stg1_3)**(3 / 4)

delta_h_id_stg2_3 = (w_t2 / 100) * delta_h_id
N_stg2_3 = N * np.sqrt(m / rho_4) / (delta_h_id_stg2_3)**(3 / 4)

delta_h_id_stg3_3 = (w_t3 / 100) * delta_h_id
N_stg3_3 = N * np.sqrt(m / rho_6) / (delta_h_id_stg3_3)**(3 / 4)

# Comprimento da Turbina
L = np.sum(d) + np.sum(cz)

# Angulo do canal axial
theta_rad = np.arctan(0.5 * (h4 - h0) / L)
theta = np.degrees(theta_rad)

print(f'Ns: {Ns}')
print(f'N_stg1: {N_stg1}')
print(f'N_stg2: {N_stg2}')
print(f'w_t1: {w_t1}')
print(f'w_t2: {w_t2}')
print(f'w_t3: {w_t3}')
print(f'rm: {rm}')
print(f'alpha2: {alpha2}')
print(f'betha2: {betha2}')
print(f'Psi_stg1: {Psi_stg1}')
print(f'R1: {R1}')
print(f'betha1: {betha1}')
print(f'alpha1: {alpha1}')
print(f'alpha4: {alpha4}')
print(f'betha4: {betha4}')
print(f'Psi_stg2: {Psi_stg2}')
print(f'R2: {R2}')
print(f'betha3: {betha3}')
print(f'alpha6: {alpha6}')
print(f'betha6: {betha6}')
print(f'Psi_stg3: {Psi_stg3}')
print(f'R3: {R3}')
print(f'betha5: {betha5}')
print(f'M_cr1: {M_cr1}')
print(f'Mr_cr2: {Mr_cr2}')
print(f'M_cr3: {M_cr3}')
print(f'Mr_cr4: {Mr_cr4}')
print(f'M_cr5: {M_cr5}')
print(f'Mr_cr6: {Mr_cr6}')
print(f'L: {L}')
print(f'theta: {degrees(theta)}')