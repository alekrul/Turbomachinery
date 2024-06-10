# Radial turbine calculation
#ME211 - Turbomáquinas prof. Dr. Takachi
#Autor: Alexandre KRUL (krul@ita.br) and Elias AHTITICH (Elias.AHTITICH@student.isae-supaero.fr)

from math import *
import numpy as np
from slip_factor import slipFactor

#Given Values:
psi = 1.04 #power input factor
N = 290 #rotational speed [rev/s]
D_imp = 0.5 #overall impeller diameter [m]
D_eye_tip = 0.3 #eye tip diameter [m]
D_eye_root = 0.15 #eye root diameter [m]
m_dot = 9 #air mass flow [kg/s]
T_01 = 295 #inlet stagnation temperatura [K]
p_01 = 1.1 #inlet stagnation pressure [bar]
eta = 0.78 #isoentropic efficiency
w = 0.05 #radial width of vaneless space
r_d_mean = 0.33 #aproximate mean radius of diffuser throat
depth = 0.0176 #depth of diffuser passages
Z_diff = 12 #number of diffuser vanes

#constants
cp = 1005 # Specific heat at constant pressure
R = 287 # Specific gas constant for air in J/(kg*K)
gama = 1.4

#Impeller parameters assumed based on Eckert rotors database
Z = 20
t = 0.003 #assuming 3mm blade tickness based on eckardt rotors (rotor A, B and O has 3mm) - a further investigations based on total weight and stress can be done here
meridional_angle = pi/2 #meridional inclination angle (radial)


#requirements:
#(a) to determine the pressure ratio of the compressor and the power required to drive it 
# assuming that the velocity of the air at the inlet is axial; 
#(b) to calculate the inlet angle of the impeller vanes at the root and tip radii of the eye, 
# assuming that the axial inlet velocity is constant across the eye annulus; 
# and (c) to estimate the axial depth of the impeller channels at the periphery of the impeller.

#Procedure:

#we will first assume Slip Factor as 0.9 for initial value
slip_factor_init = 0.9 #initial value for slip factor

U = pi*N*D_imp #rotational speed at tip [m/s]

max_iterations = 10
err = 1e-3

slip_factor = slip_factor_init

for i in range(max_iterations):
    Delta_T = psi*slip_factor*(U**2)/cp
    annulus_area = pi*((D_eye_tip**2) - (D_eye_root**2))/4

    #estimate Va1 (inlet axial velocity)

    rho1 = p_01*100000/(R*T_01)
    Va1 = m_dot/(rho1*annulus_area)
    for j in range(max_iterations):
    
        T_1 = T_01 - ((Va1**2)/(2*cp))

        p_1 = p_01/((T_01/T_1)**(gama/(gama-1)))

        rho1 = p_1*100000/(R*T_1)

        Va1_new = m_dot/(rho1*annulus_area)
        v_change = abs(Va1-Va1_new)/Va1
        Va1 = Va1_new
        if v_change <= err:
            break

    V1 = [0,Va1]

    U_eye_tip = pi*D_eye_tip*N
    U_eye_root = pi*D_eye_root*N

    beta_eye_root = atan(Va1/U_eye_root)
    beta_eye_tip = atan(Va1/U_eye_tip)

    # (c)
    Vr2 = Va1 #choosing radial velocity at outlet == axial velocity at inlet
    Vt2 = slip_factor*U #whirl velocity (tangencial)

    V2 = sqrt(Vr2**2 + Vt2**2)

    #assuming half the total loss occurs in the impeller
    impeller_eta = 1 - (0.5*(1-eta))


    p_02 = p_01*((1+(impeller_eta*Delta_T/T_01))**(gama/(gama-1)))

    T_02 = Delta_T + T_01

    T_2 = T_02 - ((V2**2)/(2*cp))

    p_2 = p_02*((T_2/T_02)**(gama/(gama-1)))

    rho2 = p_2*100000/(R*T_2)

    A_out = m_dot/(rho2*Vr2)

    d = A_out/(pi*D_imp)

    Wr = Vr2
    Wt = U - Vt2
    W = sqrt(Wr**2 + Wt**2)

    beta_tip = atan(Wt/Wr)

    turning_rate = (beta_tip-beta_eye_tip)/(D_imp-D_eye_tip)
    

    
    slip_factor_new = slipFactor(Z,[U,0],[Vt2,Vr2],[Wt,Wr],D_eye_tip/2,D_imp/2,beta_tip,meridional_angle,t,turning_rate,0,rho2,d,"qiu")
    change = abs(slip_factor_new-slip_factor)/slip_factor
    slip_factor=slip_factor_new

    if change <= err:
        break


p_ratio = (1 + ((eta*Delta_T)/T_01))**(gama/(gama-1))

Power = m_dot*cp*Delta_T

#Diffuser calculation
r_imp = D_imp/2
r_d = r_imp+w #diffuser vane leading edge radius

Vtd = Vt2*(r_imp/r_d)

#find radial component with Vr2 as first assumption

Vrd = Vr2
A_diffuser = 2*pi*r_d*depth

for j in range(max_iterations):
    Vd = sqrt(Vrd**2 + Vtd**2)
    T_d = T_02 - ((Vd**2)/(2*cp))

    p_d_p_02 = (T_d/T_02)**(gama/(gama-1))

    p_d = p_01*p_d_p_02*((1+(impeller_eta*Delta_T/T_01))**(gama/(gama-1)))

    rhod = p_d*100000/(R*T_d)

    Vrd_new = m_dot/(rhod*A_diffuser)
    Vrd_change = abs((Vrd-Vrd_new)/Vrd)
    Vrd = Vrd_new
    if  Vrd_change <= err:   
        break


beta_d = atan(Vrd/Vtd)

#calculate throat width

Vtd_throat = Vt2*(r_imp/r_d_mean)

#find radial component with Vrd as first assumption
Vrd_throat = Vrd
A_diffuser_t = 2*pi*r_d_mean*depth

for j in range(max_iterations):
    Vd_throat = sqrt(Vrd_throat**2 + Vtd_throat**2)
    T_d = T_02 - ((Vd_throat**2)/(2*cp))

    p_d_p_02 = (T_d/T_02)**(gama/(gama-1))

    p_d = p_01*p_d_p_02*((1+(impeller_eta*Delta_T/T_01))**(gama/(gama-1)))

    rhod = p_d*100000/(R*T_d)

    Vrd_t_new = m_dot/(rhod*A_diffuser_t)
    Vrd_t_change = abs((Vrd_throat-Vrd_t_new)/Vrd_throat)
    Vrd_throat = Vrd_t_new
    if  Vrd_t_change <= err:   
        break

beta_t = atan(Vrd_throat/Vtd_throat)

throat_area = A_diffuser_t*sin(beta_t)
width_throat = throat_area/(Z_diff*depth)


#RESULTS WITH NO INLET GUIDE VANES
print("RESULTS ASSUMING INLET AIR VELOCITY IS AXIAL")

print(f"\u0394T:{Delta_T:.1f} K")
print(f"Pressure ratio:{p_ratio:.3f}")
print(f"Power required: {(Power/1000):.1f} kW")
print(f"Inlet angle at impeller eye tip β1t:{degrees(beta_eye_tip):.2f} degrees")
print(f"Inlet angle at impeller eye root β1r:{degrees(beta_eye_root):.2f} degrees")
print(f"Depth of impeller channel: {d:.5f} m")
print(f"Slip factor calculated \u03C3: {slip_factor:.5f}")
print(f"Diffuser vanes leading edge angles β1d: {degrees(beta_d):.2f} degrees")
print(f"Width of the diffuser throat: {width_throat:.5f} m")
#--------------------------------------------------------------

#Compressibility Effects

W1 = sqrt(Va1**2 + U_eye_tip**2)
a1 = sqrt(gama*R*T_1)

M_eye_tip = W1/a1 #mach number at tip of the eye
print(f"Relative Mach Number at eye tip: {M_eye_tip:.3f}")

print("Too high, we will assume a Inlet guide vane to make a pre-swirl of 30 degrees")
swirl_angle = radians(30)

#re-calculating everything to check changes in blade parameters, work, geometry

slip_factor_init = slip_factor #considering previous slip factor calculation as initial value

slip_factor = slip_factor_init

annulus_area = pi*((D_eye_tip**2) - (D_eye_root**2))/4
for i in range(max_iterations):

    #estimate Va1 (inlet axial velocity)

    rho1 = p_01*100000/(R*T_01)
    Va1 = m_dot/(rho1*annulus_area)
    
    for j in range(max_iterations):
        V1 = Va1/cos(swirl_angle)
        T_1 = T_01 - ((V1**2)/(2*cp))

        p_1 = p_01/((T_01/T_1)**(gama/(gama-1)))

        rho1 = p_1*100000/(R*T_1)

        Va1_new = m_dot/(rho1*annulus_area)
        v_change = abs(Va1-Va1_new)/Va1
        Va1 = Va1_new
        if v_change <= err:
            break

    U_eye_tip = pi*D_eye_tip*N
    U_eye_root = pi*D_eye_root*N
    U_eye = (U_eye_root+U_eye_tip)/2

    Vt1 = Va1*tan(swirl_angle)

    Wt1_eye_tip = U_eye_tip - Vt1
    Wt1_eye_root = U_eye_root - Vt1

    beta_eye_root = atan(Va1/Wt1_eye_root)
    beta_eye_tip = atan(Va1/Wt1_eye_tip)

    Vr2 = Va1 #choosing radial velocity at outlet == axial velocity at inlet
    Vt2 = slip_factor*U #whirl velocity (tangencial)

    V2 = sqrt(Vr2**2 + Vt2**2)

    #assuming half the total loss occurs in the impeller
    impeller_eta = 1 - (0.5*(1-eta))

    Delta_T = psi*(slip_factor*(U**2) - Vt1*U_eye)/cp

    p_02 = p_01*((1+(impeller_eta*Delta_T/T_01))**(gama/(gama-1)))

    T_02 = Delta_T + T_01

    T_2 = T_02 - ((V2**2)/(2*cp))

    p_2 = p_02*((T_2/T_02)**(gama/(gama-1)))

    rho2 = p_2*100000/(R*T_2)

    A_out = m_dot/(rho2*Vr2)

    d = A_out/(pi*D_imp)

    Wr = Vr2
    Wt = U - Vt2
    W = sqrt(Wr**2 + Wt**2)

    beta_tip = atan(Wt/Wr)

    turning_rate = (beta_tip-beta_eye_tip)/(D_imp-D_eye_tip)
    
    slip_factor_new = slipFactor(Z,[U,0],[Vt2,Vr2],[Wt,Wr],D_eye_tip/2,D_imp/2,beta_tip,meridional_angle,t,turning_rate,0,rho2,d,"qiu")
    change = abs(slip_factor_new-slip_factor)/slip_factor
    slip_factor=slip_factor_new

    if change <= err:
        break


p_ratio = (1 + ((eta*Delta_T)/T_01))**(gama/(gama-1))

Power = m_dot*cp*Delta_T

#Diffuser calculation
r_imp = D_imp/2
r_d = r_imp+w #diffuser vane leading edge radius

Vtd = Vt2*(r_imp/r_d)

#find radial component with Vr2 as first assumption

Vrd = Vr2
A_diffuser = 2*pi*r_d*depth

for j in range(max_iterations):
    Vd = sqrt(Vrd**2 + Vtd**2)
    T_d = T_02 - ((Vd**2)/(2*cp))

    p_d_p_02 = (T_d/T_02)**(gama/(gama-1))

    p_d = p_01*p_d_p_02*((1+(impeller_eta*Delta_T/T_01))**(gama/(gama-1)))

    rhod = p_d*100000/(R*T_d)

    Vrd_new = m_dot/(rhod*A_diffuser)
    Vrd_change = abs((Vrd-Vrd_new)/Vrd)
    Vrd = Vrd_new
    if  Vrd_change <= err:   
        break


beta_d = atan(Vrd/Vtd)

#calculate throat width

Vtd_throat = Vt2*(r_imp/r_d_mean)

#find radial component with Vrd as first assumption
Vrd_throat = Vrd
A_diffuser_t = 2*pi*r_d_mean*depth

for j in range(max_iterations):
    Vd_throat = sqrt(Vrd_throat**2 + Vtd_throat**2)
    T_d = T_02 - ((Vd_throat**2)/(2*cp))

    p_d_p_02 = (T_d/T_02)**(gama/(gama-1))

    p_d = p_01*p_d_p_02*((1+(impeller_eta*Delta_T/T_01))**(gama/(gama-1)))

    rhod = p_d*100000/(R*T_d)

    Vrd_t_new = m_dot/(rhod*A_diffuser_t)
    Vrd_t_change = abs((Vrd_throat-Vrd_t_new)/Vrd_throat)
    Vrd_throat = Vrd_t_new
    if  Vrd_t_change <= err:   
        break

beta_t = atan(Vrd_throat/Vtd_throat)

throat_area = A_diffuser_t*sin(beta_t)
width_throat = throat_area/(Z_diff*depth)

#RESULTS WITH INLET GUIDE VANES (PRE-SWIRL OF 30 DEGREES)
print("RESULTS WITH INLET GUIDE VANES (PRE-SWIRL OF 30 DEGREES)")

print(f"\u0394T:{Delta_T:.1f} K")
print(f"Pressure ratio:{p_ratio:.3f}")
print(f"Power required: {(Power/1000):.1f} kW")
print(f"Inlet angle at impeller eye tip β1t:{degrees(beta_eye_tip):.2f} degrees")
print(f"Inlet angle at impeller eye root β1r:{degrees(beta_eye_root):.2f} degrees")
print(f"Depth of impeller channel: {d:.5f} m")
print(f"Slip factor calculated \u03C3: {slip_factor:.5f}")
print(f"Diffuser vanes leading edge angles β1d: {degrees(beta_d):.2f} degrees")
print(f"Width of the diffuser throat: {width_throat:.5f} m")