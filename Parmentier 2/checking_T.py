import numpy as np


# Hd209458b
Ab = 7.1010310640033283E-002 # Bond albedo
f = 0.25                     # averaged profile
sigma_b = 5.670374419E-8     # Stefan-Boltzmann constant
T_int = 280.8020223944342    # internal temperature
T_star = 6026.3549454750     # stellar temperature
R_star = 1.1999759942400 * (695700E3) # stellar radius [m]
D = R_star * 8.81            # orbital distance [m]
mu = 1/np.sqrt(3)            # cosine of the zenith angle


# A. Check the Equilibrium Temperature (T_eq0 直接從母恆星射過來的) 
# Equation (1) in Parmentier 2015
T_eq0 = T_star * np.sqrt(R_star/(2*D))
print("(1)T_eq0  = {:.3f} K ".format(T_eq0))

T_eq0 = 1448

# Equation (2) in Parmentier 2015
T_irr0 = (4*T_eq0**4)**(1/4)
print("(2)T_irr0 = {:.3f} K ".format(T_irr0))
# Equation (3) in Parmentier 2015
T_mu0 = (mu*T_irr0**4)**(1/4)
print("(3)T_mu0  = {:.3f} K ".format(T_mu0))
# Equation (4) in Parmentier 2015
T_mu = (mu*(T_irr0**4)*(1-Ab))**(1/4)
print("(4)T_mu   = {:.3f} K ".format(T_mu))
# 4-1
T_irr = ((1-Ab)*T_irr0**4)**(1/4)
print("(4-1)T_irr= {:.3f} K ".format(T_irr))
# Equation (5) in Parmentier 2015
T_mu2 = ((T_eq0**4)*(1-Ab)*(4*f))**(1/4)
print("(5)T_mu2  = {:.3f} K ".format(T_mu2))

T_mu3 = (mu*(T_irr)**4)**(1/4)
print("( )T_mu3  = {:.3f} K ".format(T_mu3))

# Equation (6) in Parmentier 2015
T_eff = (T_mu**4 + T_int**4 )**(1/4)
print("(6)T_eff  = {:.3f} K ".format(T_eff))


T_irr_gu = (1-Ab)*T_star * np.sqrt(R_star/(D))
print("T_irr_gu = {:.3f} K ".format(T_irr_gu))