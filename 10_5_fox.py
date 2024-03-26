from turbopython import *
from math import *

beta1 = radians(40)
beta2 = radians(60)
r1 = 380/1000
r2 = 1140/1000
b1 = 120/1000
b2 = 80/1000
omega = 575 * 2*3.1415*(1/60) #rpm to rad/s
Q = 18000/3600 #m3/h to m3/s

Vn1 = flowtoVel(Q,r1,b1)
Vn2 = flowtoVel(Q,r2,b2)
W1abs = Vn1/sin(beta1)
W2abs = Vn2/sin(beta2)
W1_vector = Wvetorial(W1abs,beta1)
W2_vector = Wvetorial(W2abs,beta2)

U1 = TipSpeed(omega,r1)
U2 = TipSpeed(omega,r2)

V1vector, V1  = Vabs(U1,W1_vector)
V2vector, V2 = Vabs(U2,W2_vector)

mass_flow = Q*997

T = torque(r2,r1,V1vector,V2vector,mass_flow)
Work = T*omega

print("Vn1", Vn1)
print("Vn2",Vn2)
print("W1", W1_vector)
print("W2",W2_vector)
print("U1", U1)
print("U2",U2)
print("V1", V1)
print("V1 vector", V1vector)
print("V2 vector",V2vector)
print("Torque", T)
print("Work",Work)