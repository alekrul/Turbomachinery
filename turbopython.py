#!/usr/bin/env python
# coding: utf-8


from math import *


def TipSpeed(omega,r):
    return [omega*r,0]

def slipFactor(n):
    #from saravanamutto book -> only for gas turbines
    return (1 - (0.63*pi/n))

def power(U2,U1,V1,V2,mass_flow):
    #considering steady flow, no friction, 1D flow in the inlet and outlet, pressure effects negligible
    return mass_flow*(U2[0]*V2[0] - U1[0]*V1[0])

def powerWithSlipFactor(n,U,mass_flow):
    #from saravanamutto book -> only for gas turbines
    return slipFactor(n)*(U**2)

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