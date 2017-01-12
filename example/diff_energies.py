import numpy as np 
import matplotlib.pyplot as plt


data = np.loadtxt('energies.dat')

index_occ = data[:,1]==2
Eocc = data[index_occ,2]
Eocc = np.unique(Eocc[::2])

index_virt = data[:,1]==0
Evirt = data[index_virt,2]
Evirt = np.unique(Evirt[::2])


xocc = 0.5*(Eocc + np.roll(Eocc,1))[1:]
deocc = (Eocc - np.roll(Eocc,1))[1:]

xvirt = 0.5*(Evirt + np.roll(Evirt,1))[1:]
devirt = (Evirt - np.roll(Evirt,1))[1:]

#plt.plot(xocc,deocc)
#plt.plot(xvirt,devirt)
#plt.show()


Wph = 0.017
Wph_max = 2*Wph

Xocc = np.linspace(min(Eocc)-1,max(Eocc)+1,500)
dos_occ = np.zeros(500)
for i in range(len(deocc)):
	if deocc[i] < Wph_max:
		dos_occ += Wph**2/((Xocc-xocc[i])**2+Wph**2)

Wph = 0.017
Xvirt = np.linspace(min(Evirt)-1,max(Evirt)+1,500)
dos_virt = np.zeros(500)
for i in range(len(devirt)):
	if devirt[i] < Wph_max:
		dos_virt += Wph**2/((Xvirt-xvirt[i])**2+Wph**2)
	

plt.plot(Xocc,dos_occ)
plt.plot(Xvirt,dos_virt)
plt.show()

