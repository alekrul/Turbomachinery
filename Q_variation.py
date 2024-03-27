from turbopython import *
from math import *
import numpy as np
import matplotlib.pyplot as plt

r2 = 150/1000 #m
b2 = 6.25/1000 #m
Q = [0.10, 0.50, 0.90, 1.30, 1.70, 2.10, 2.50, 2.90, 3.30]
omega = 750*2*pi*(1/60) #rpm to rad/s

beta2 = 40

U2 = TipSpeed(omega,r2)
W_dot = np.zeros(len(Q))

for i in range(len(Q)):
    Vn2 = flowtoVel(Q[i],r2,b2)
    print("------------------------------")
    print("Q:",Q[i])
    W2abs = Vn2/sin(beta2)
    W2 = Wvetorial(W2abs,beta2)
    V2,V2abs = Vabs(U2,W2)
    W_dot[i] = power(U2,[0,0],[0,0],V2,Q[i]*999)

plt.plot(Q,W_dot)
plt.xlabel("Q")
plt.ylabel(r'$\dot{W}$')
plt.title("Variação de "+r'$\dot{W}$'+" a partir da mudança de Q")
plt.show()



