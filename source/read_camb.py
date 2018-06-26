import sys, platform, os
from matplotlib import pyplot as plt
import numpy as np
from IPython.display import display
import six
import warnings
import scipy.integrate as integrate
warnings.filterwarnings('ignore', category=DeprecationWarning, module='.*/IPython/.*')

sys.path.insert(0,r'c:\work\dist\git\camb\pycamb')
import camb
from camb.symbolic import *
sympy.init_printing()

monopole_source, ISW, doppler, quadrupole_source = get_scalar_temperature_sources()

pars = camb.set_params(H0=70.0, ombh2=0.01809, omch2=0.12400, As=2.1e-9, ns=0.96, tau=0.08)
from matplotlib import rcParams
rcParams.update( {'axes.labelsize': 14,
							'font.size': 14,
							'legend.fontsize': 14,
							'xtick.labelsize': 13,
							'ytick.labelsize': 13})
cl_label= r'$\ell(\ell+1)C_\ell/2\pi\quad [\mu {\rm K}^2]$'

#Example of plotting the time evolution of Newtonian gauge variables and the monopole sources 

# data= camb.get_background(pars)
# conformal_times = np.linspace(1, 800, 300)
# ks = [0.01,0.05]
# Delta_g_N = make_frame_invariant(Delta_g, 'Newtonian')
# display('Temperature monopole source in general', monopole_source)
# display('Temperature monopole source in Newtonian gauge', newtonian_gauge(monopole_source))

# ev = data.get_time_evolution(ks, conformal_times, [Psi_N, monopole_source])
# _, axs= plt.subplots(1,2, figsize=(14,6))
# for i, ax in enumerate(axs):
#     ax.plot(conformal_times,ev[i,:, 0])
#     ax.plot(conformal_times,ev[i,:, 1]*1000)

#     ax.set_title('$k= %s$'%ks[i])
#     ax.set_xlabel(r'$\eta/\rm{Mpc}$')
#     ax.set_xlim(conformal_times[0], conformal_times[-1])
# plt.legend([r'$\Psi_N$',
#             r'Monopole (x1000)'])
omega_m = 0.3
h=0.7
def D_A(z): #angular diameter distance, in h^-1*Mpc
		ans= integrate.quad(lambda x:1.0/
										np.sqrt(omega_m*(1.+x)*(1.+x)*(1.+x)+1-omega_m),0.0,z)[0]
		return 3000.0 * ans/(1.0+z)

def chi(z): #comiving distance, assuming flat universe, in Mpc, h= 0.7
		return D_A(z) * (1.0+z) / h

z_cmb = np.array([1100.0])

data= camb.get_background(pars)
ks = np.logspace(-5.,np.log10(2.),500)
Delta_g_N = make_frame_invariant(Delta_g, 'Newtonian')
display('Temperature monopole source in general', monopole_source)
display('Temperature monopole source in Newtonian gauge', newtonian_gauge(monopole_source))

ev = data.get_redshift_evolution(ks, z_cmb, [Psi_N, monopole_source])

plt.figure(figsize=(14,6))
plt.semilogx(ks,ev[:,0, 0])
plt.semilogx(ks,ev[:,0, 1])
plt.semilogx(ks,ev[:,0, 0]+ev[:,0, 1])

plt.title(r'monopole and $\Psi_N$ at recombination')
plt.xlabel(r'$k/\rm{Mpc^{-1}}$')

plt.legend([r'$\Psi_N$',r'Monopole',r'$\Psi_N+Monopole$'])
plt.show()

Psi_N = ev[np.arange(len(ks)),0,0]
monopole = ev[np.arange(len(ks)),0,0]
np.savetxt('Psi_N.txt',Psi_N,header='Psi_N')
np.savetxt('monopole.txt',monopole,header='monopole')




