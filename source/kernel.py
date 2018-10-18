import numpy as np
import matplotlib.pyplot as plt
import distance as d

z=np.arange(0.1,10.,.1)
n=np.size(z)
chistar=d.chi(1189.)
chia=np.zeros(n)
k1=np.zeros(n)
for i,za in enumerate(z):
	chia[i]=d.chi(za)
	k1[i]=d.hubinv(za)*(chistar-chia[i])/(chistar**2)
plt.plot(z,(1+z)*d.hubinv(z)/chistar,'b',label='Time Delay')
plt.plot(z,(1+z)*k1,'--r',label='Deflection Angle')
#plt.xscale('log')
plt.axis([0.1,10,0,.4])
plt.xlabel('$z$',fontsize=18)
plt.xscale('log')
plt.legend()
plt.ylabel('$K(z_*,z)$',fontsize=18)
plt.savefig('../plots/kernel.pdf')

