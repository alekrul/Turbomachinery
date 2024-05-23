from math import *
import numpy as np

def slipFactor(Z,U ,V = [0,0],r2 = 0, beta = pi/6,gamma = pi/2,t = 0.0015, turning_rate = 1,drho_dm = 0,rho2 = 1,b=1):
    Sf = 1
    
    s = 2*pi*r2/Z
    F = 1 - (sin(pi/Z)*sin((pi/Z) + beta)*cos(beta)*sin(gamma)) - t/(s*cos(beta))
    delta_sigma_radial = F*pi*cos(beta)*sin(gamma)/Z
    phi = V[1]/U
    delta_sigma_turn = F*s*phi*turning_rate/(4*cos(beta))
    delta_sigma_passage = -F*phi*s*sin(beta)*drho_dm/(4*rho2*b)
    Sf = 1 - delta_sigma_radial - delta_sigma_turn - delta_sigma_passage
    
    return Sf