from thermo import *
#https://pypi.org/project/thermo/#description

water = Chemical('water')
water.calculate(T=298, P=101325)
print("Water Rho at T=298K and P=101325",water.rho)

air = Mixture('air',T=1000)
print("Air Cp at 1000K:",air.Cp)
