from numpy import sqrt,linspace,exp,arange,zeros
from scipy import integrate
from scipy.interpolate import interp1d



omegam=0.31
littleh=0.68
omegar=2.47e-5*(1.+0.68)/littleh**2
infty = 10000.

def hubinv(z):
	tmp = sqrt(omegam*(1.+z)**3+(1.-omegam)+omegar*(1.+z)**4)
	hubinv = 3000./tmp
	return hubinv
	
def chi(z): 
	# input z and omegam and compute comoving distance, angular diameter distance, growth function, and H^{-1}(z) [in h Mpc^{-1}]
	tmp, err = integrate.quad(hubinv, 0., z)
	return tmp

def dist(z):
    # input z and omegam and compute comoving distance, angular diameter distance, growth function, and H^{-1}(z) [in h Mpc^{-1}]
    tmp, err = integrate.quad(hubinv, z, infty)
    return tmp

def growthint(z):
	growthint = (1.+z)*(hubinv(z)/3000.)**3
	return growthint
	
def growth(z):
	tmp, errg = integrate.quad(growthint,z,100.)
	growth = 2.5*omegam*tmp/(hubinv(z)/3000.)
	return growth

z=.83
print chi(z)
zl=.3
zs=1.
print chi(zs)/(1.+zs),chi(zl)/(1.+zl)

#z = arange(0.5,.72,.1)
#da = zeros(z.size)
#for i in range(z.size):
#	da[i] = chi(z[i])/(1.+z[i])
#	print z[i],chi(z[i]),da[i]/littleh,hubinv(z[i])

#z_to_da = interp1d(z, da, kind='linear')

#print chi(z) # chi
#print chi(z)/(1.+z) # d_A
#print hubinv(z)
#print growth(z)
