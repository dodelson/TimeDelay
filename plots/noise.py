import numpy as np

import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import scipy.special as special
from scipy import interpolate
from sympy.physics.wigner import wigner_3j

filename = ''
font = {#'family' : 'normal',
                #'weight' : 'bold',
                'size'  : 18}

plt.rc('font', **font)

omega_m = 0.3
h=0.7
def D_A(z): #angular diameter distance, in Mpc
    ans= integrate.quad(lambda x:1.0/
                    np.sqrt(omega_m*(1.+x)*(1.+x)*(1.+x)+1-omega_m),0.0,z)[0]
    return 3000.0 * ans/((1.0+z)*h)

def chi(z): #comiving distance, assuming flat universe, in Mpc
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

#matter power spectrum at last scattering surface
p_lin_R = np.loadtxt('OUTPUT3/matter_power_lin/p_k.txt')[499]
p_nl_R = np.loadtxt('OUTPUT3/matter_power_nl/p_k.txt')[499]

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
# def integrand(chi,k,l):
#     return special.spherical_jn(l, (chi_Rmax - chi)*k)\
#                 *interp_D(interp_z(chi))

monopole = np.loadtxt('../source/data/monopole.txt')
Psi_N = np.loadtxt('../source/data/Psi_N.txt')
klin_cmb = np.loadtxt('../source/data/klin_cmb.txt')
nklin_cmb = len(klin_cmb)
#monopole + Psi_N spectrum
P_0 = (monopole+Psi_N)*(monopole+Psi_N)

def interp1d(lnk,karray,funcarray):
    interp_eg = interpolate.interp1d(karray,karray*funcarray)
    return interp_eg(np.exp(lnk))

# func_k = np.zeros(nk_nl)
# C_l = np.zeros(160)
# for l in range(1,161):
#     for i, kloop in enumerate(k_nl):
#         int_0 = integrate.quad(integrand, 0.0, chi_Rmax,args= (kloop,l),epsrel=1e-6)[0]
#         func_k[i] = int_0*int_0*p_nl[i]/(k_nl[i]*k_nl[i])
   
#     k_min = np.log(k_nl.min())
#     k_max = np.log(k_nl.max())
#     term = integrate.quad(interp1d, k_min, k_max,args=(k_nl,func_k), epsrel = 1e-6)[0]
#     units = l*(l+1)*1e-12*term*omega_m*omega_m/(9.*np.pi*np.pi*chi_Rmax*chi_Rmax)
#     print(l,units)
#     C_l[l-1] = units

C_l1 = np.zeros(2601) #eq 25 in the paper
func_k1 = np.zeros(nklin_cmb)

for l in range(0,2601):
    print(l)
    for i,kloop in enumerate(klin_cmb):
        func_k1[i] = kloop*kloop*special.spherical_jn(l,kloop*chi_Rmax)*P_0[i]*\
                    (special.spherical_jn(l,kloop*chi_Rmax)*l - special.spherical_jn(l+1,kloop*chi_Rmax)*kloop*chi_Rmax)\
                    /(8.*np.pi*np.pi*np.pi)
        # func_k1[i] = p_nl_R[i]*kloop*kloop*special.spherical_jn(l,kloop*chi_Rmax)*\
        #             (special.spherical_jn(l,kloop*chi_Rmax)*l - special.spherical_jn(l+1,kloop*chi_Rmax)*kloop*chi_Rmax)\
        #             /(np.power(2.*np.pi,3.))
    k_min = np.log(klin_cmb.min())
    k_max = np.log(klin_cmb.max())
    term = integrate.quad(interp1d, k_min, k_max,args=(klin_cmb,func_k1), epsrel = 1e-6)[0]
    C_l1[l] = term

print(chi_Rmax)

