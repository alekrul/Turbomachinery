from math import *
import numpy as np

def slipFactor(Z,U = [1,1],V = [0,0],W=[0,0],r1 = 1,r2 = 0, beta = pi/6,gamma = pi/2,t = 0.0015, turning_rate = 1,drho_dm = 0,rho2 = 1,b=1, method="Stodola"):
    Sf = 1
    if method == "Saravanamuttoo":
        #from saravanamutto book
        Sf = 1 - (0.63*pi/Z)
    elif method == "Stodola":
        #Stodola method
        Sf = 1 - ((pi*cos(beta))/Z*(1-(V[1]/U[0])*tan(beta)))
    elif method == "Stanitz":
        #Stanitz method
        A = (W[1]/U[0])*(1/tan(beta))
        B = 0.63*pi/Z
        Sf = 1 - B*(1-(1/A))
    elif method == "Balje":
        Sf = 1/(1+(6.2/(Z*((r2/r1)**(2/3)))))
    elif method == "qiu":
        s = 2*pi*r2/Z
        F = 1 - (sin(pi/Z)*sin((pi/Z) + beta)*cos(beta)*sin(gamma)) - t/(s*cos(beta))
        delta_sigma_radial = F*pi*cos(beta)*sin(gamma)/Z
        phi = V[1]/U[0]
        delta_sigma_turn = F*s*phi*turning_rate/(4*cos(beta))
        delta_sigma_passage = -F*phi*s*sin(beta)*drho_dm/(4*rho2*b)
        Sf = 1 - delta_sigma_radial - delta_sigma_turn - delta_sigma_passage
    else:
        print("Couldn't find method")
    return Sf

def slipFactor(Z,U = [1,1],V = [0,0],W=[0,0],r1 = 1,r2 = 0, beta = pi/6,gamma = pi/2,t = 0.0015, turning_rate = 1,drho_dm = 0,rho2 = 1,b=1, method="Stodola"):
    Sf = 1
    if method == "Saravanamuttoo":
        #from saravanamutto book
        Sf = 1 - (0.63*pi/Z)
    elif method == "Stodola":
        #Stodola method
        Sf = 1 - ((pi*cos(beta))/Z*(1-(V[1]/U[0])*tan(beta)))
    elif method == "Stanitz":
        #Stanitz method
        A = (W[1]/U[0])*(1/tan(beta))
        B = 0.63*pi/Z
        Sf = 1 - B*(1-(1/A))
    elif method == "Balje":
        Sf = 1/(1+(6.2/(Z*((r2/r1)**(2/3)))))
    elif method == "qiu":
        s = 2*pi*r2/Z
        F = 1 - (sin(pi/Z)*sin((pi/Z) + beta)*cos(beta)*sin(gamma)) - t/(s*cos(beta))
        delta_sigma_radial = F*pi*cos(beta)*sin(gamma)/Z
        phi = V[1]/U[0]
        delta_sigma_turn = F*s*phi*turning_rate/(4*cos(beta))
        delta_sigma_passage = -F*phi*s*sin(beta)*drho_dm/(4*rho2*b)
        Sf = 1 - delta_sigma_radial - delta_sigma_turn - delta_sigma_passage
    
    return Sf