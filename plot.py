import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

rc('text', usetex=True)

data1 = np.loadtxt('out/euler.csv', delimiter=',', skiprows=1)
data2 = np.loadtxt('out/pre-cor.csv', delimiter=',', skiprows=1)
t = data1[:,0]
ke1 = data1[:,1]
pe1 = data1[:,2]
ke2 = data2[:,1]
pe2 = data2[:,2]
KE = ke1[0]
PE = pe1[0]

fig,(ax1,ax2) = plt.subplots(2,1,sharex=True)
ax1.plot(t,ke1/KE,'-b',label='Euler')
ax2.plot(t,ke2/KE,'-r',label='Predictor-Corrector')
ax2.set_xlabel('time $t$')
ax1.set_title('Error comparison for iterative schemes: KE($t$)/KE($t=0$)')
ax1.legend(loc='upper right')
ax2.legend(loc='upper right')
fig.tight_layout()
fig.savefig('out/ke.png')

fig,(ax1,ax2) = plt.subplots(2,1,sharex=True)
ax1.plot(t,pe1/PE,'-b',label='Euler')
ax2.plot(t,pe2/PE,'-r',label='Predictor-Corrector')
ax2.set_xlabel('time $t$')
ax1.set_title('Error comparison for iterative schemes: PE($t$)/PE($t=0$)')
ax1.legend(loc='upper right')
ax2.legend(loc='upper right')
fig.tight_layout()
fig.savefig('out/pe.png')
