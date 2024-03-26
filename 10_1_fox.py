from turbopython import *
from math import *

#Given a centrifugal pump dimensions below 
r1 = 175/1000 #m
r2 = 500/1000 #m
b1 = 50/1000 #m
b2 = 30/1000 #m
beta1 = radians(65) #radians
beta2 = radians(70) #radians

#and operation:
omega = 750*2*pi*(1/60) #rpm to rad/s

Q = 0.75 #m3/s

#FIND W_dot

#W_dot = mass_flow(U2*Vt2 - U1*Vt1), so we need to find each variable

Vn1 = flowtoVel(Q,r1,b1)
Vn2 = flowtoVel(Q,r2,b2)


#W*sin(beta) = Vn so 
W1abs = Vn1/sin(beta1)
W2abs = Vn2/sin(beta2)


W1 = Wvetorial(W1abs,beta1)
W2 = Wvetorial(W2abs,beta2)

#find U = omega*r
U1 = TipSpeed(omega,r1)
U2 = TipSpeed(omega,r2)

#with U and W we can find V
V1vector, V1abs = Vabs(U1,W1)
V2vector, V2abs = Vabs(U2,W2)

#now we have Vt, U and mass_flow = Q*rho
mass_flow = Q*999

#find Work
W_dot = power(U2,U1,V1vector,V2vector,mass_flow)
H = W_dot/9.81
print("W dot:", W_dot)
print("_____________________")
print("V1:", V1vector)
print("V2:", V2vector)
print("W1:", W1)
print("W2:", W2)
print("U1:", U1)
print("U2:", U2)

#evaluate slip factor if there are 12 blades:

print("Slip factor:")
sf = slipFactor(12,U2,V2vector,W2,r1,r2,beta2,method="Balje")
print("Balje:",sf)
sf = slipFactor(12,U2,V2vector,W2,r1,r2,beta2,method="Saravanamuttoo")
print("Saravanamuttoo:",sf)
sf = slipFactor(12,U2,V2vector,W2,r1,r2,beta2,method="Stodola")
print("Stodola:",sf)
sf = slipFactor(12,U2,V2vector,W2,r1,r2,beta2,method="Stanitz")
print("Stanitz:",sf)