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

C_l1 = np.zeros(161) #eq 25 in the paper
func_k1 = np.zeros(nk_nl)

for l in range(0,161):
    for i,kloop in enumerate(k_nl):
        func_k1[i] = omega_m*omega_m*(1.+z[nz-1])*(1.+z[nz-1])*p_nl_R[i]*special.spherical_jn(l,kloop*chi_Rmax)*\
                    (special.spherical_jn(l,kloop*chi_Rmax)*l - special.spherical_jn(l+1,kloop*chi_Rmax)*kloop*chi_Rmax)\
                    /(16.*np.power(4*np.pi,3)*9e12*kloop*kloop)
        # func_k1[i] = p_nl_R[i]*kloop*kloop*special.spherical_jn(l,kloop*chi_Rmax)*\
        #             (special.spherical_jn(l,kloop*chi_Rmax)*l - special.spherical_jn(l+1,kloop*chi_Rmax)*kloop*chi_Rmax)\
        #             /(np.power(2.*np.pi,3.))
    k_min = np.log(k_nl.min())
    k_max = np.log(k_nl.max())
    term = integrate.quad(interp1d, k_min, k_max,args=(k_nl,func_k1), epsrel = 1e-6)[0]
    C_l1[l] = term


def f_expect(l1,l2,L):
    return np.sqrt((2*l1+1)*(2*l2+1)/(4*np.pi*(2*L+1)))*wigner_3j(l2,l1,L,0,0,0)*(C_l1[l1]+np.power(-1,l1+l2+L)*C_l1[l2])

def g_weight(l1,l2,L):
    return 1

C_CMB = np.arange(200)
def d_variance(L,M):
    sum3 = 0.
    for l1 in range(161):
        print(l1)
        sum2 = 0.
        for l2 in range(161):
            sum1 = 0.
            for m1 in range(-l1,l1+1):
                sum0 = 0.
                if M==0:
                    for m2 in range(-l2,l2+1):
                        sum0 = sum0 + g_weight(l1,l1,L)*np.conj(g_weight(l2,l2,L))*wigner_3j(l1,l1,L,m1,-m1,0)*np.conj(wigner_3j(l2,l2,L,m2,-m2,0))
                    sum0 = sum0 + g_weight(l1,l2,L)*wigner_3j(l2,l1,L,m1,-m1,0)*np.conj(wigner_3j(l2,l1,L,m1,-m1,0))*(np.conj(g_weight(l1,l2,L))+np.power(-1,l1+l2+L)*np.conj(g_weight(l2,l1,L)))
                else:
                    sum0 = sum0 + g_weight(l1,l2,L)*wigner_3j(l2,l1,L,m1-M,-m1,M)*np.conj(wigner_3j(l2,l1,L,m1-M,-m1,M))*(np.conj(g_weight(l1,l2,L))+np.power(-1,l1+l2+L)*np.conj(g_weight(l2,l1,L)))
                sum1 = sum1 + sum0
            sum2 = sum2 + sum1*C_CMB[l2]
        sum3 = sum3 + sum2*C_CMB[l1]
    return sum3

print(d_variance(8,0))

