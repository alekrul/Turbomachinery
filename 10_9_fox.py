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

#find Q that makes Vt1 = 0, then finds W_dot

U1 = TipSpeed(omega,r1)
# U = Vt + W*cos(beta), but Vt=0 so W = U/cos(beta)
W1abs = U1[0]/cos(beta1)
W1vector = Wvetorial(W1abs,beta1)

V1n = W1vector[1]
Q = Volumeflow(V1n,r1,b1)
print("Q:",Q)

#now to find W_dot we need to find Vt1 and U2
U2 = TipSpeed(omega,r2)

V2n = flowtoVel(Q,r2,b2)
W2abs = V2n/sin(beta2)

W2vector = Wvetorial(W2abs,beta2)

V2,v2absvalue = Vabs(U2,W2vector)
V1 = [0,V1n]
mass_flow = 999*Q
W_dot = power(U2,U1,V1,V2,mass_flow)

print("W_dot:",W_dot)
