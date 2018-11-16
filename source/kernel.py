import numpy as np
import matplotlib.pyplot as plt
import distance as d

z=np.arange(0.01,10.,.01)
n=np.size(z)
chistar=d.chi(1189.)
chia=np.zeros(n)
k1=np.zeros(n)
for i,za in enumerate(z):
	chia[i]=d.chi(za)
	k1[i]=d.hubinv(za)*(chistar-chia[i])/(chistar*chia[i])
plt.plot(z,d.hubinv(z)/chistar,'b',label='Time Delay')
plt.plot(z,k1,'--r',label='Deflection Angle')
#plt.xscale('log')
plt.axis([0,10,0.003,100])
plt.xlabel('$z$',fontsize=18)
plt.yscale('log')
plt.legend()
plt.ylabel('$K(z)$',fontsize=18)
plt.savefig('../plots/kernel.pdf')

