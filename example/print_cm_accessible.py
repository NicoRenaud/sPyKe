import numpy as np
import matplotlib.pyplot as plt

# read the data
data = np.loadtxt('energies.dat')

# find the energies of the occpuied orbitals
index_occ = data[:,1]==2
Eocc = data[index_occ,2]
nocc = len(Eocc)
print ' %d occupied orbitals found' %nocc

# find the energies of the virtual orbitals
index_virt = data[:,1]==0
Evirt = data[index_virt,2]
nVirt = len(Evirt)
print ' %d virtual orbitals found' %nVirt



# gap of the system 
Egap = Evirt[0]-Eocc[-1]
print 'Fondamental gap of the system %1.3f eV' %Egap



# HOMO/LUMO energy
EH = Eocc[-1]
EL = Evirt[0]

# width of the Lorentzian
Wph = 0.001



################################################
# transition from an excited electron
################################################

data = []
k=0
for ivirt in range(1,nVirt):

	# initial energy of the excited electrons
	E0 = Evirt[ivirt]
	Eex = E0-EH

	# for all the possible relaxation
	for ivirt2 in range(ivirt-1):

		# final energy of the first electrons
		E1 = Evirt[ivirt2]

		# for all the Impact ionization excitation 
		for ivirt3 in range(ivirt-1):
			E2 = Evirt[ivirt3]

			DE_loss = E0-E1
			DE_gained = EH-E2

			# delta energy of the process
			DE = DE_loss + DE_gained

			# compute a lorentzian of width equal to the phoonon energy
			L = Wph**2/(DE**2+Wph**2)
			
			# store the results
			if L>1E-6:
				data.append([Eex,L,nocc+ivirt,nocc+ivirt2,nocc-1,nocc+ivirt3])
				k+=1

maxEx = Evirt[-1]-EH
data = np.array(data)
np.savetxt('data_ex.dat',data)
X = np.linspace(0,maxEx*1.25,500)
DOS = np.zeros(500)
for i in range(len(data)):
	DOS += data[i,1]*0.010**2/((X-data[i,0])**2+0.010**2)


data_save = np.array([X/Egap,DOS])
plt.plot(X/Egap,DOS)
plt.xlabel('E_{gamma}/E_{gap}')
np.savetxt('cm_acceissble_elec.dat',data_save.T)

################################################
# transition from an excited hole
################################################
data = []
k=0

for iocc in range(nocc-2,-1,-1):

	# initial energy of the hole
	E0 = Eocc[iocc]
	Eex = EL-E0

	# for ll the possible relaxation
	for iocc2 in range(nocc-1,iocc,-1):

		# energy of the relaxed state
		E1 = Eocc[iocc2]

		for iocc3 in range(nocc-1,iocc,-1):

			E2 = Eocc[iocc3]

			DE_loss = E1-E0
			DE_gained = E2-EL

			DE = DE_loss + DE_gained

			# compute a lorentzian of width equal to the phoonon energy
			L = Wph**2/(DE**2+Wph**2)
			
			# store the results
			if L>1E-6:
				data.append([Eex,L,iocc,iocc2,iocc3,nocc+1])
				k+=1

maxEx = EL-Eocc[0]
data = np.array(data)
np.savetxt('data_hole.dat',data)
X = np.linspace(0,maxEx*1.25,500)
DOS = np.zeros(500)
for i in range(len(data)):
	DOS += data[i,1]*0.010**2/((X-data[i,0])**2+0.010**2)


data_save = np.array([X/Egap,DOS])
np.savetxt('cm_acceissble_hole.dat',data_save.T)
plt.plot(X/Egap,DOS)
plt.xlabel('E_{gamma}/E_{gap}')
plt.show()


