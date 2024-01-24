import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import dedalus.public as d3
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.interpolate import interp1d
import logging


sigma_B = 5.67037442E-8 #kg s^-3 K^-4

#-------Checking--------
CHECK_BASIS = False
#-----------------------
Nx = 16
Ny = 16


coords = d3.CartesianCoordinates('phi', 'ln_sigma')
dist = d3.Distributor(coords, dtype=np.complex128)
#phi_basis = d3.ComplexFourier(coords['phi'], size=Nx, bounds=(-0.5*np.pi, 0.5*np.pi)) # phi can be 0 to 2pi
phi_basis = d3.Chebyshev(coords['phi'], size=Nx, bounds=(0, 0.5*np.pi)) # phi can be 0 to 2pi
#ln_sigma_basis = d3.ComplexFourier(coords['ln_sigma'], size=Ny, bounds=(np.log(0.5E-5), 0)) # log(sigma) can be  np.log(0.5E-5) to 0 -> sigma = 0.5E-5 to 1
ln_sigma_basis = d3.Chebyshev(coords['ln_sigma'], size=Ny, bounds=(np.log(0.5E-5), 0)) # log(sigma) can be  np.log(0.5E-5) to 0 -> sigma = 0.5E-5 to 1

if CHECK_BASIS:
    plt.figure(figsize=(6, 1.5))
    plt.plot(ln_sigma_basis.global_grid().ravel(),0*ln_sigma_basis.global_grid().ravel()+1,'o',markersize=2,label='Sigma')
    plt.legend()
    plt.ylim(0.5,1.5)
    plt.gca().yaxis.set_ticks([])
    plt.tight_layout()
    plt.savefig('basis_distribution_sigma.png',dpi=300)

    plt.figure(figsize=(6, 1.5))
    plt.plot(phi_basis.global_grid().ravel(),0*phi_basis.global_grid().ravel()+1,'o',markersize=2,label='Phi')
    plt.legend()
    plt.ylim(0.5,1.5)
    plt.gca().yaxis.set_ticks([])
    plt.tight_layout()
    plt.savefig('basis_distribution_phi.png',dpi=300)


phi,ln_sigma = dist.local_grids(phi_basis,ln_sigma_basis) # sigma = np.exp(ln_sigma)

# Set the field (See whether the parameter is (phi,ln(sigma)))
u = dist.Field(name='u', bases=(phi_basis, ln_sigma_basis), dtype=np.complex128) # u is function of phi,ln(sigma)
v = dist.Field(name='v', bases=(phi_basis, ln_sigma_basis), dtype=np.complex128) # v is function of phi,ln(sigma)
P = dist.Field(name='P', bases=(phi_basis, ln_sigma_basis), dtype=np.complex128) # P is function of phi,ln(sigma)


# u0 = dist.Field(name='u0', bases=(phi_basis),dtype=np.complex128)
epsilon = dist.Field(name='epsilon', bases=(phi_basis,ln_sigma_basis),dtype=np.complex128)
Kappa = dist.Field(name='kappa',bases=(ln_sigma_basis),dtype=np.complex128)
GAMMA0 = dist.Field(name='Gamma0',bases=(ln_sigma_basis),dtype=np.complex128)

# set the constant's field
Rp = dist.Field(name='Rp',dtype=np.complex128)
Omega0 = dist.Field(name='Omega0',dtype=np.complex128)
OMEGA = dist.Field(name='OMEGA',dtype=np.complex128)
R = dist.Field(name='R',dtype=np.complex128)
Gamma0 = dist.Field(name='Gamma0', bases=(ln_sigma_basis),dtype=np.complex128)
cp = dist.Field(name='cp',dtype=np.complex128)

# Define fields for sine and cosine of phi
sin_phi = dist.Field(name='sin_phi', bases=(phi_basis),dtype=np.complex128)
cos_phi = dist.Field(name='cos_phi', bases=(phi_basis),dtype=np.complex128)

# Initialize these fields
sin_phi['g'] = np.sin(phi)
sin_phi['g'].imag = 0
cos_phi['g'] = np.cos(phi)
print(cos_phi['g'])
cos_phi['g'].imag=0
# Define sigma field and set its values
sigma = dist.Field(name='sigma', bases=(phi_basis, ln_sigma_basis),dtype=np.complex128)
sigma['g'] = np.exp(ln_sigma)
sigma['g'].imag=0
# Asign the value

phi_grid = phi_basis.global_grid().ravel()
np.savetxt('phi.csv', phi_grid, delimiter=',')
ln_sigma_grid = ln_sigma_basis.global_grid().ravel()
# Save to CSV file
np.savetxt('ln_sigma_grid.csv', ln_sigma_grid, delimiter=',')

## import the data
data = pd.read_csv('./PTprofile(20bar)(Hd209458b)(1206).csv',sep=',')
# print(data.columns)
Betav1 = data['Betav1'][0]
Betav2 = data['Betav2'][0]
Betav3 = data['Betav3'][0]
Gv1 = data['Gv1'][0]
Gv2 = data['Gv2'][0]
Gv3 = data['Gv3'][0]
sigma_ = data['sigma'] 
kappa = data['Kappa(m^2/kg)']
gamma0_ = data['GAMMA_0']
tau = data['tau']
T_irr = data['T_irr(K)'][0]
T = data['T(K)']

def interpolating(old_data_grid,old_data,new_data_grid):
    old_data_grid = np.array(old_data_grid)
    old_data = np.array(old_data)
    f = interp1d(old_data_grid,old_data,kind='slinear', fill_value='extrapolate')
    # New set of grid-values : new_data_grid
    new_data = f(new_data_grid)
    return new_data

def plot_interpolation(new_data,old_data,name):
    plt.figure(figsize=(5,4))
    plt.plot(np.log(sigma_),old_data,'o',color='blue',label='Original data',markersize=0.75)
    plt.plot(np.array(ln_sigma),np.array(new_data),'x',color='red',markersize=0.75)
    plt.yscale('log')
    plt.xlabel('ln(sigma)')
    plt.title(f"Interpolating {name}")
    plt.ylabel(name)
    plt.legend()
    plt.savefig(f'1003_check_interpolation_{name}.png',dpi=300)


# interpolating kappa,gamma0_
new_kappa = interpolating(np.log(sigma_),kappa,ln_sigma)
new_gamma0 = interpolating(np.log(sigma_),gamma0_,ln_sigma)
new_tau = interpolating(np.log(sigma_),tau,ln_sigma)
new_T = interpolating(np.log(sigma_),T,ln_sigma)

GAMMA0['g']=new_gamma0
GAMMA0['g'].imag=0
Kappa['g']=new_kappa
Kappa['g'].imag=0
# Assign the value :

# Hd209458b
G = 6.67E-11
M_star = 1.06917519587500 *(1988500E24) #kg
R_star = 1.1999759942400 * (695700E3) #m
T_star = 6026.3549454750 #K
Teq0 = 1448 #K 
Tirr0 = (4*Teq0**4)**(0.25)
mu=1/np.sqrt(3.0) 
Ab = 7.1010310640033283E-002 #from Parmentier
Tmu = (mu*(Tirr0**4)*(1-Ab))**(0.25)
Omega = np.sqrt(G*M_star/(R_star*((T_star/Tmu)**2))**3)


Omega0['g'] = 0.000000001
Omega0['g'].imag=0
OMEGA['g'] = Omega
OMEGA['g'].imag= 0
Gamma0['g'] = new_gamma0
Gamma0['g'].imag=0
cp['g']= 12330.08474
cp['g'].imag=0
Rp['g'] = 15.6 * (6371000)
Rp['g'].imag =0 
R['g'] = data['R(J/(K*g))'][0] * 1000 # SI : J/K/kg
R['g'].imag = 0
print("Rp:"+str(Rp['g']))
print("R:"+str(R['g']))
print("Omega0 : "+str(Omega0['g']))
print("Omega : "+str(OMEGA['g']))
#exit()
# Create a 2D meshgrid for phi and ln_sigma
PHI, LN_SIGMA = np.meshgrid(phi, ln_sigma)

new_tau_reshaped = new_tau.reshape(Ny, 1)
new_kappa_reshaped = new_kappa.reshape(Ny, 1)


epsilon['g'] = ((np.pi*2*Betav1*Gv1*new_kappa_reshaped*sigma_B*T_irr**4*np.exp(-np.sqrt(3)*Gv1*new_tau_reshaped)*np.cos(PHI))\
                +(np.pi*2*Betav2*Gv2*new_kappa_reshaped*sigma_B*T_irr**4*np.exp(-np.sqrt(3)*Gv2*new_tau_reshaped)*np.cos(PHI))\
                +(np.pi*2*Betav3*Gv3*new_kappa_reshaped*sigma_B*T_irr**4*np.exp(-np.sqrt(3)*Gv3*new_tau_reshaped)*np.cos(PHI))).T
# print(epsilon['g'])
epsilon['g'].imag=0
print(epsilon['g'])
np.savetxt('epsilon.csv', epsilon['g'], delimiter=',')
np.savetxt('tau.csv', new_tau_reshaped, delimiter=',')
np.savetxt('Gamma0.csv', Gamma0['g'], delimiter=',')
# 4.

# tau1 = dist.Field(name='tau1')
# tau2 = dist.Field(name='tau2')
# tau3 = dist.Field(name='tau3')
#tau_1 = dist.Field(name='tau_1')
tau_1 = dist.Field(name='tau_1', bases=phi_basis)
#tau_2 = dist.Field(name='tau_2')
tau_2 = dist.Field(name='tau_2', bases=phi_basis)
#tau_3 = dist.Field(name='tau_3')
tau_3 = dist.Field(name='tau_3', bases=ln_sigma_basis)
#tau_4 = dist.Field(name='tau_4')
tau_4 = dist.Field(name='tau_4', bases=ln_sigma_basis)

#substitution
dln_sigma = lambda A: d3.Differentiate(A, coords['ln_sigma'])
lift_ln_sigma_basis = ln_sigma_basis.derivative_basis(2)
lift_ln_sigma = lambda A, n: d3.Lift(A, lift_phi_basis, n)

#lift_basis_y = ybasis.derivative_basis(2)
#lift_x = lambda A, n: d3.Lift(A, lift_basis_y, n)

dphi = lambda B:d3.Differentiate(B, coords['phi'])
lift_phi_basis =phi_basis.derivative_basis(2)
lift_phi = lambda B, n : d3.Lift(B, lift_ln_sigma_basis, n)
print("lift define done!")
#problem = d3.LBVP([u,v,P,tau_1,tau_3,tau_4], namespace=locals())
problem = d3.LBVP([u,v,P,tau_1,tau_2,tau_3,tau_4], namespace=locals())
# Add Equation 1:
problem.add_equation("1j*Omega0*u - (2*OMEGA+Omega0)*sin_phi*v + 1j*P/cos_phi = 0")
#problem.add_equation("cos_phi*1j*Omega0*u - cos_phi*(2+Omega0)*sin_phi*v + 1j*P = 0")
print("Eq. 1 succeed")
# Add Equation 2:
problem.add_equation("1j*Omega0*v + (2*OMEGA+Omega0)*sin_phi*u + dphi(P)/Rp+lift_phi(tau_1,-1) = 0")
#problem.add_equation("1j*Omega0*v + (2*OMEGA+Omega0)*sin_phi*u + dphi(P)/Rp+lift_phi(tau_1,-1) = 0")
#problem.add_equation("1j*Omega0*v + (2+Omega0)*sin_phi*u - dphi(P) = 0")
print("Eq. 2 succeed")

# Add Equation 3:
#problem.add_equation("1j*Omega0*dln_sigma(1/Gamma0*dln_sigma(P)) -sigma*(dphi(v) - v*sin_phi/cos_phi + 1j*u/cos_phi) = -1/cp*dln_sigma(epsilon/Gamma0)")
problem.add_equation("Rp*1j*Omega0*dln_sigma(1/(Gamma0)*dln_sigma(P)) -sigma*(dphi(v) - v*sin_phi/cos_phi + 1j*u/cos_phi)+lift_phi(tau_2,-1)+lift_ln_sigma(tau_3,-1)+lift_ln_sigma(tau_4,-2) = -Rp/cp*dln_sigma(epsilon/Gamma0)")
#problem.add_equation("1j*Omega0*dln_sigma(1/(Gamma0*R)*dln_sigma(P)) -sigma*(dphi(v) - v*sin_phi/cos_phi + 1j*u/cos_phi)+lift_phi(tau_2,-1)+lift_ln_sigma(tau_3,-1)+lift_ln_sigma(tau_4,-2) = -1/cp*dln_sigma(epsilon/Gamma0)")
#problem.add_equation("cos_phi*1j*Omega0*dln_sigma(1/Gamma0*dln_sigma(P)) -sigma*(cos_phi*dphi(v) - v*sin_phi + 1j*u)+lift_phi(tau_2,-1)+lift_ln_sigma(tau_3,-1)+lift_ln_sigma(tau_4,-2) = -cos_phi/cp*dln_sigma(epsilon/Gamma0)")

print("Add equation done!")



# Add boundary condition 1:
problem.add_equation("v(phi=0)=0")
#problem.add_equation("u(phi=np.pi/2)=0")
# Add boundary condition 2:
problem.add_equation("P(phi=np.pi/2)=0")
# Add boundary condition 3:
ln_sigma_left = np.log(0.5E-5)
ln_sigma_right = 0
problem.add_equation("dln_sigma(P)(ln_sigma=np.log(0.5E-5)) = 0")
# Add boundary condition 4:
GAMMA0_bc = dist.Field(name='GAMMA0_bc',dtype=np.complex128)
GAMMA0_bc['g'] = new_gamma0[0][-1]
GAMMA0_bc['g'].imag=0
print("new_gamma0 = ",new_gamma0[0][-1])
T0_bc = dist.Field(name='T0_bc',dtype=np.complex128)
T0_bc['g'] = new_T[0][-1]
T0_bc['g'].imag=0
print("T0_bc = ", new_T[0][-1])
#exit()

print(new_gamma0[0][-1]/new_T[0][-1])
#problem.add_equation("dln_sigma(P)(ln_sigma=0)+GAMMA0_bc/T0_bc*P=0")
problem.add_equation("dln_sigma(P)(ln_sigma=0)=0")


print('add equation...successful!')
#exit()
# # Print data types of fields
# fields = [u, v, P, epsilon, Kappa, GAMMA0, Rp, Omega0, R, Gamma0, cp, sin_phi, cos_phi, sigma]
# for field in fields:
#     print(f"Data type of {field.name}: {field['g'].dtype}")

# # Print data types of arrays
# arrays = [PHI, LN_SIGMA]
# for array in arrays:
#     print(f"Data type of {array}: {array.dtype}")
# exit()

print("Solving..")
# Set up root logger
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)  # Set to DEBUG for verbose output

# You can also set up a handler for logging to a file or console with a specific format
handler = logging.StreamHandler()
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)
solver = problem.build_solver()
solver.solve()
logger.setLevel(logging.WARNING)  # Change to WARNING or INFO after debugging
print("Solving...Done")
# Gather global data
phi = phi_basis.global_grid()
ln_sigma = ln_sigma_basis.global_grid()
ug = u.allgather_data('g')
vg = v.allgather_data('g')
Pg = P.allgather_data('g')
#print('--------------------')
#print(ug)
#print('--------------------')
#print(vg)
#print('--------------------')
#print(Pg)
#exit()
print("Saving Pg....")
np.savetxt('Pg.csv', Pg, delimiter=',')
print("Saving Pg... Done")
np.savetxt('ug.csv', ug, delimiter=',')
np.savetxt('vg.csv', vg, delimiter=',')
array_str = np.array2string(Pg)
exit()

#with open('Pg_result.txt','w') as f:
#    f.write(array_str)


if dist.comm.rank == 0:
#    plt.figure(figsize=(6, 4))
    count = 0
    for name in [ug,vg,Pg]:
        names = ['ug','vg','Pg']
        fig, axs = plt.subplots(1, 2, figsize=(10, 5))
        magnitude = np.abs(name)
        phase = np.angle(name)
        c1 = axs[0].imshow(magnitude, cmap='viridis')
        c2 = axs[1].imshow(phase, cmap='hsv')
   
        fig.colorbar(c1, ax=axs[0], label='Magnitude')
        fig.colorbar(c2, ax=axs[1], label='Phase')

        axs[0].set_title('Magnitude')
        axs[1].set_title('Phase')
        plt.savefig(f'{names[count]}.png')
        

        x = np.arange(name.shape[1])
        y = np.arange(name.shape[0])
        X, Y = np.meshgrid(x, y)
        Z = np.abs(name)

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.plot_surface(X, Y, Z, cmap='viridis')

        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Magnitude')
        plt.savefig(f'{names[count]}_3D.png')
        count += 1
#    plt.pcolormesh(phi.ravel(), ln_sigma.ravel(), ug.T, cmap='viridis', shading='gouraud', rasterized=True)
#    plt.gca().set_aspect('equal')
#    plt.xlabel('phi')
#    plt.ylabel('ln_sigma')
#    plt.title("TEST")
#    plt.tight_layout()
#    plt.colorbar()
#    plt.savefig('test_1103.png', dpi=200)
#    plt.show()


