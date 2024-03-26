#!/usr/bin/env python
# coding: utf-8


from math import *


def TipSpeed(omega,r):
    return [omega*r,0]

def slipFactor(n,U = [1,1],V = [0,0],W=[0,0],r1 = 1,r2 = 0, beta = pi/6,method="Balje"):
    Sf = 1
    if method == "Saravanamuttoo":
        #from saravanamutto book
        Sf = 1 - (0.63*pi/n)
    elif method == "Stodola":
        #Stodola method
        Sf = 1 - ((pi*cos(beta))/n*(1-(V[1]/U[0])*tan(beta)))
    elif method == "Stanitz":
        #Stanitz method
        A = (W[1]/U[0])*(1/tan(beta))
        B = 0.63*pi/n
        Sf = 1 - B*(1-(1/A))
    elif method == "Balje":
        Sf = 1/(1+(6.2/(n*((r2/r1)**(2/3)))))
    else:
        print("Couldn't find method")
    return Sf
    

def power(U2,U1,V1,V2,mass_flow):
    #considering steady flow, no friction, 1D flow in the inlet and outlet, pressure effects negligible
    return mass_flow*(U2[0]*V2[0] - U1[0]*V1[0])

def powerWithSlipFactor(Sf,U):
    return Sf*(U**2)

def Wvetorial (W, beta):
    return [-W*cos(beta),W*sin(beta)]

def Vabs(U, W):
    V = [U[0]+W[0],U[1]+W[1]]
    return V, sqrt(V[0]**2 + V[1]**2)

def flowtoVel(Q, r, b):
    #continuity equation considering steady flow, uniform 1D flow in the inlet and outlet, incompressible flow
    return Q/(2*pi*r*b)

def torque(r2,r1,V1,V2,mass_flow):
    #considering steady flow, no friction, 1D flow in the inlet and outlet, pressure effects negligible
    return mass_flow*(r2*V2[0] - r1*V1[0])

def Volumeflow(Vn,r,b):
    return Vn*(2*pi*r*b)

def flowAngle(V):
    return degrees(atan(V[0]/V[1]))