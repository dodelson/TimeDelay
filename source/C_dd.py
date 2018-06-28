import numpy as np

import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import scipy.special as special
from scipy import interpolate

filename = ''
font = {#'family' : 'normal',
                #'weight' : 'bold',
                'size'  : 18}

plt.rc('font', **font)

omega_m = 0.3

def D_A(z): #angular diameter distance, in h^-1*Mpc
    ans= integrate.quad(lambda x:1.0/
                    np.sqrt(omega_m*(1.+x)*(1.+x)*(1.+x)+1-omega_m),0.0,z)[0]
    return 3000.0 * ans/(1.0+z)

def chi(z): #comiving distance, assuming flat universe, in h^-1*Mpc
    return D_A(z) * (1.0+z)

z = np.concatenate((np.loadtxt('OUTPUT1/matter_power_lin/z.txt'),
                    np.loadtxt('OUTPUT2/matter_power_lin/z.txt'),
                    np.loadtxt('OUTPUT3/matter_power_lin/z.txt')))
nz = len(z)
k_lin = np.loadtxt('OUTPUT1/matter_power_lin/k_h.txt')
k_nl = np.loadtxt('OUTPUT1/matter_power_nl/k_h.txt')

p_lin = np.loadtxt('OUTPUT1/matter_power_lin/p_k.txt')[0]
p_nl = np.loadtxt('OUTPUT1/matter_power_nl/p_k.txt')[0]
nk_nl = len(k_nl)
nk_lin = len(k_lin)

#growth function
D_delta = np.concatenate((np.loadtxt('OUTPUT1/growth_parameters/d_z.txt'),
                          np.loadtxt('OUTPUT2/growth_parameters/d_z.txt'),
                          np.loadtxt('OUTPUT3/growth_parameters/d_z.txt')))
D_phi = D_delta + D_delta*z

#comiving distance in redshift bins
chi_R = np.zeros(nz)
for i in range(nz):
    chi_R[i] = chi(z[i])
chi_Rmax = chi_R[nz-1]

#Interpolation on D_phi
interp_D = interpolate.interp1d(z, D_phi)
#Interpolation on redshift
interp_z = interpolate.interp1d(chi_R,z)
#First integration
def integrand(chi,k,l):
    return special.spherical_jn(l, (chi_Rmax - chi)*k)\
                *interp_D(interp_z(chi))

def interp1d(lnk,karray,funcarray):
    interp_eg = interpolate.interp1d(karray,karray*funcarray)
    return interp_eg(np.exp(lnk))

func_k = np.zeros(nk_nl)
C_l = np.zeros(160)
for l in range(1,161):
    for i, kloop in enumerate(k_nl):
        int_0 = integrate.quad(integrand, 0.0, chi_Rmax,args= (kloop,l),epsrel=1e-6)[0]
        func_k[i] = int_0*int_0*p_nl[i]/(k_nl[i]*k_nl[i])
   
    k_min = np.log(k_nl.min())
    k_max = np.log(k_nl.max())
    term = integrate.quad(interp1d, k_min, k_max,args=(k_nl,func_k), epsrel = 1e-6)[0]
    units = l*(l+1)*1e-12*term*omega_m*omega_m/(9.*np.pi*np.pi*chi_Rmax*chi_Rmax)
    print(l,units)
    C_l[l-1] = units

X = range(1,161)
Y = C_l

# plt.figure(figsize=(25, 16))
plt.figure(figsize=(14, 6))
plt.loglog(X,Y, label = 'Time Delay Spectrum')
plt.xlabel('$l$')
plt.ylabel('$l(l+1)C_{l}^{dd}/2\pi)$')
plt.legend()
plt.savefig('test.png')
plt.show()
