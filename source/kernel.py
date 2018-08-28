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
plt.plot(z,d.hubinv(z)/chistar,'b',label='Time Delay')
plt.plot(z,k1,'--r',label='Deflection Angle')
#plt.xscale('log')
plt.xlabel('Redshift $z$',fontsize=18)
plt.legend()
plt.ylabel('Kernel (Unnormalized)',fontsize=18)
plt.savefig('../plots/kernel.pdf')

