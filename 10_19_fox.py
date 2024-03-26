from turbopython import *
from math import *

#Given a centrifugal pump dimensions below 
r1 = 75/1000 #m
r2 = 150/1000 #m
b1 = 7.5/1000 #m
b2 = 6.25/1000 #m
beta1 = radians(25) #radians
beta2 = radians(40) #radians

#and operation:
Q = 30*(1/1000) #from L/s to m3/s

#find: velocities diagram at inlet, Project velocity for Vt1 = 0, velocities diagram at outlet, alfa2 and work_dot

V1n = flowtoVel(Q,r1,b1)
W1abs = V1n/sin(beta1)

W1 = Wvetorial(W1abs,beta1)

# U = Vt + W*cos(beta), but Vt=0 so U = W*cos(beta)
U1 = [-W1[0],0]
V1 = [0,V1n]

V2n = flowtoVel(Q,r2,b2)
W2abs = V2n/sin(beta2)
W2 = Wvetorial(W2abs,beta2)

omega = U1[0]/r1
print("Omega:", omega)

U2 = TipSpeed(omega,r2)

V2, V2abs = Vabs(U2,W2)

alfa_2 = flowAngle(V2)

print("alfa 2:",alfa_2)

mass_flow = 999*Q
W_dot = power(U2,U1,V1,V2,mass_flow)
print("W_dot:",W_dot)