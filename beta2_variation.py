from turbopython import *
from math import *
import numpy as np
import matplotlib.pyplot as plt

r2 = 150/1000 #m
b2 = 6.25/1000 #m
Q = 30*(1/1000) #from L/s to m3/s
omega = 750*2*pi*(1/60) #rpm to rad/s

beta2 = [pi/6, pi/4, pi/3, pi/6+pi/4, pi/2, pi/4+pi/3, 2*pi/3, 3*pi/4, 5*pi/6]

Vn2 = flowtoVel(Q,r2,b2)

U2 = TipSpeed(omega,r2)

W_dot = np.zeros(len(beta2))
V2abs = np.zeros(len(beta2))
W2abs = np.zeros(len(beta2))

for i in range(len(beta2)):
    print("------------------------------")
    print("ANGLE Beta2:",degrees(beta2[i]))
    W2abs[i] = Vn2/sin(beta2[i])
    W2 = Wvetorial(W2abs[i],beta2[i])
    print(W2)
    print("W2: ",W2)
    V2,V2abs[i] = Vabs(U2,W2)
    print("V2: ",V2)
    W_dot[i] = power(U2,[0,0],[0,0],V2,Q*999)
    print("W_dot: ",W_dot[i])

beta2_degrees = [30,45,60,75,90,105,120,135,150]
plt.plot(beta2_degrees,W_dot)
plt.xlabel("β")
plt.ylabel(r'$\dot{W}$')
plt.title("Variação de "+r'$\dot{W}$'+" a partir da mudança de β")
plt.show()


fig, ax1 = plt.subplots()

color = 'tab:blue'
ax1.set_xlabel('β')
ax1.set_ylabel('Absolute Velocity', color=color)
ax1.plot(beta2_degrees, V2abs, color=color)
ax1.tick_params(axis='y', labelcolor=color)

ax2 = ax1.twinx()  
color = 'tab:red'
ax2.set_ylabel('Relative Velocity', color=color)  
ax2.plot(beta2_degrees, W2abs, color=color)
ax2.tick_params(axis='y', labelcolor=color)

# Exibindo o gráfico
fig.tight_layout()  
plt.show()