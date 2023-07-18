import numpy as np

# WASP-19 b
Q_int = 0.00009297509039270726E-4
Mp = 1.133*1.89813E27 #kg
Rp = 1.386*69911E3 #m

# HD 209458 b
# Q_int = 0.0003053355867099633E-4
# Mp = 0.69*1.89813E27 #kg
# Rp = 1.359*69911E3 #m

sigma_b = 5.670374419E-8 #W m^-2 K^-4

result = (Q_int*Mp/(4*np.pi*Rp**2*sigma_b))**(1/4)
print(result)
