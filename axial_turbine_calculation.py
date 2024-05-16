# Exemplo de cálculo de turbina axial dada em sala de aula
#Disciplina ME211 - Turbomáquinas prof. Dr. Takachi
#Autor: Alexandre Mendonça Krul (krul@ita.br)

from math import *
import numpy as np

#dados:
alpha_1 = 0 #[degrees]
m_dot = 8.2 #[kg/s]
p_in = 10 #[bar]
T_in = 1380 #[K]
p_ratio = 9.3 #[adim]
rpm = 35000 #[rpm]

#constraints
Um = 414 #[m/s]
max_stages = 3
aspect_ratio_3 = 3

#assumptions
eff = 0.91
gamma = 1.33
cp = 1158


yy = (gamma-1)/gamma
T_ex = T_in*(1-eff*(1-((1/p_ratio)**yy)))
p_ex = p_in/p_ratio

#assumindo que a densidade total na saida é próxima do valor da densidade estática na saída:
rho_ex = p_ex*100/(0.287*T_ex)

entalphy = cp*T_in*(1-((1/p_ratio)**yy))

N = rpm *2*pi/60
Ns = (N*sqrt(m_dot/rho_ex))/(entalphy**(3/4))
print("Considerando 2 estágios, com distribuição de trabalho 55:45")
#assumindo 2 estágios e com distribuição de trabalho 55:45
entalphy_s1 = 0.55*entalphy
entalphy_s2 = entalphy-entalphy_s1

#cálculo para o primeiro estágio
p_int = p_in*((1-entalphy_s1/(cp*T_in))**(1/yy))
T_int = T_in*(1-eff*(1-((p_int/p_in)**yy)))
rho_int = p_int*100/(0.287*T_int)

Ns_1 = (N*sqrt(m_dot/rho_int))/(entalphy_s1**(3/4))
print("Velocidade específica no primeiro estágio: ", Ns_1)
#cálculo para o segundo estágio
p_ext_o = p_int*((1-entalphy_s2/(cp*T_int))**(1/yy))
T_ext_o = T_int*(1-eff*(1-((p_ext_o/p_int)**yy)))
rho_ext_o = p_ext_o*100/(0.287*T_ext_o)

Ns_2 = (N*sqrt(m_dot/rho_ext_o))/(entalphy_s2**(3/4))
print("Velocidade específica no segundo estágio: ", Ns_2)
print("dado que Ns_1 é baixo, desconfigurando turbina axial, iremos assumir 3 estágios")

#dado que Ns_1 é baixo, desconfigurando turbina axial, iremos assumir 3 estágios
stages = max_stages
p_ratio_s = p_ratio**(1/stages)

T = np.zeros(stages+1)
T[0]=T_in
p = np.zeros(stages+1)
p[0] = p_in
Ns = np.zeros(stages)

for i in range(stages):
    p[i+1] = p[i]/p_ratio_s
    T[i+1] = T[i]*(1-eff*(1-((1/p_ratio_s)**yy)))
    
#usando a diferenca de temperaturas para calcular a distribuição de trabalho entre estágios
w_s1 = (T[1]-T[0])/(T[3]-T[0])
w_s2 = (T[2]-T[1])/(T[3]-T[0])
w_s3 = (T[3]-T[2])/(T[3]-T[0])

print("Distribuição de trabalho: ",w_s1,":",w_s2,":",w_s3)
#calculo do raio médio
rm = Um/N
print("Raio médio: ",rm)

#inicializando variaveis
beta_1 = np.zeros(stages)
beta_2 = np.zeros(stages)

#primeiro estágio

#considerando alfa_3 (saida do escoamento) igual alfa_2
alpha_2 = alpha_1
flux_coeff = 0.6 #escolhido
stage_work_coeff = (cp*(T[0]-T[1]))/(Um**2)

beta_2[0] = atan(tan(alpha_2) - (1/flux_coeff))
R = (-tan(beta_2[0])*2*flux_coeff - stage_work_coeff)/2
beta_1[0] = atan(1/(2*flux_coeff)*(stage_work_coeff-2*R))
print(R) 

#o valor do grau de reação é sub otimo, devemos chegar próximo de 50%
#usando alpha_2 = -18 graus e flow coeff 0.65
alpha_2 = radians(-18)
flux_coeff = 0.65

beta_2[0] = atan(tan(alpha_2) - (1/flux_coeff))
R = (-tan(beta_2[0])*2*flux_coeff - stage_work_coeff)/2
beta_1[0] = atan(1/(2*flux_coeff)*(stage_work_coeff-2*R))
print("Reação no primeiro estágio com alfa 2 = -18 graus: ", R) #bem próximo de 50%

#segundo estágio
#adotando alpha_4 = -15
alpha_4 = radians(-15)
stage_work_coeff = (cp*(T[1]-T[2]))/(Um**2)


beta_2[1] = atan(tan(alpha_4) - (1/flux_coeff))
R = (-tan(beta_2[1])*2*flux_coeff - stage_work_coeff)/2
beta_1[1] = atan(1/(2*flux_coeff)*(stage_work_coeff-2*R))
print("Reação no segundo estágio com alfa 4 = -15 graus: ", R)

#terceiro estágio
#adotando alpha_6 = 0
alpha_6 = 0
stage_work_coeff = (cp*(T[2]-T[3]))/(Um**2)

beta_2[2] = atan(tan(alpha_6) - (1/flux_coeff))
R = (-tan(beta_2[2])*2*flux_coeff - stage_work_coeff)/2
beta_1[2] = atan(1/(2*flux_coeff)*(stage_work_coeff-2*R))
print("Reação no terceiro estágio com alfa 6 = 0 graus: ", R)
Va = flux_coeff*Um

#Verificação do número de Mach

print("Verificação do número de Mach em cada estágio")
Vcr1 = np.zeros(stages)
W1 = np.zeros(stages)
V1 = np.zeros(stages)
Vt1 = np.zeros(stages)
W2 = np.zeros(stages)
M1 = np.zeros(stages)
M2 = np.zeros(stages)
W2 = np.zeros(stages)
V2 = np.zeros(stages)
Ttr2 = np.zeros(stages)
Wcr2 = np.zeros(stages)

for i in range(stages):
    print("AVALIAÇÃO DE M PARA ESTÁGIO ", i+1)
    Vcr1[i] = sqrt(((2*gamma/(gamma+1)))*287*T[i])
    W1[i] = Va/cos(beta_1[i])
    Vt1[i] = W1[0]*sin(beta_1[i]) + Um
    V1[i] = sqrt(Vt1[i]**2 + Va**2)
    M1[i] = V1[i]/Vcr1[i]
    print("M1 ", M1[i])
    W2[i] = Va/cos(beta_2[i])
    V2[i] = Va/cos(alpha_2)

    Ttr2[i] = T[i] - ((V2[i]**2 + W2[i]**2)/(2*cp))
    print(Ttr2[i])
    Wcr2[i] = sqrt(((2*gamma/(gamma+1)))*287*Ttr2[i])
    M2[i] = W2[i]/Wcr2[i]
    print("M2 ", M2[i])
    

#Cálculo do escoamento meridional
rho1t = np.zeros(stages)
rho2t = np.zeros(stages)
rho1 = np.zeros(stages)
rho2 = np.zeros(stages)
h1 = np.zeros(stages)
h2 = np.zeros(stages)
pressure_loss = 0.05

Vcr0 = sqrt(((2*gamma/(gamma+1)))*287*T[0])
M0 = Va/Vcr0
rho0 = p_in*100/(0.287*T_in)
h0 = m_dot/(rho0*2*pi*rm*Va)

print("Calculo da altura de cada pá por estágio")
for i in range(stages):
    print("Estágio ",i+1) 
    p[i] = (1-pressure_loss)*p[i]
    rho1t[i] = (p[i]*100)/(0.287*T[i])
    xx = (gamma-1)/(gamma+1)
    zz = 1/(gamma-1)
    rho1[i] = rho1t[i]*((1-(xx*M1[i]**2))**(zz))
    h1[i] = m_dot/(rho1[i]*2*pi*rm*Va)
    print("h1:",h1[i])
    rho2t[i] = (p[i+1]*100)/(0.287*T[i+1])
    rho2[i] = rho2t[i]*((1-(xx*M2[i]**2))**(zz))
    h2[i] = m_dot/(rho2[i]*2*pi*rm*Va)
    print("h2:",h2[i])

#distancias axiais

lambda_s = [1.5,2.1,2.7]
lambda_r = [1.8,2.4,3.0]

Cz_s = h1/lambda_s
Cz_r = h2/lambda_r


z1 = ((h1[0]+h2[0])/2)/4
z2 = ((h2[0]+h1[1])/2)/2
z3 = ((h1[1]+h2[1])/2)/4
z4 = ((h2[1]+h1[2])/2)/2
z5 = ((h1[2]+h2[2])/2)/4

L = z1+z2+z3+z4+z5+sum(Cz_r)+sum(Cz_s)

print("Comprimento axial:",L)

#cálculo da velocidade específica
Ns_1 = (N*sqrt(m_dot/rho2[0]))/((w_s1*entalphy)**(3/4))
Ns_2 = (N*sqrt(m_dot/rho2[1]))/((w_s2*entalphy)**(3/4))
Ns_3 = (N*sqrt(m_dot/rho2[2]))/((w_s3*entalphy)**(3/4))

print("Velocidades específicas:",Ns_1," ,",Ns_2," ,",Ns_3)

#construção do canal axial
theta = atan(0.5*(h2[1]-h0)/L) #<15 degrees
print("Angulo do canal axial: ", degrees(theta),"<15 graus")