import numpy as np
from myclass import *
from scipy.constants import *
import matplotlib.pyplot as plt 
import sys
epsilon = 22.258
epsilon_ = 2
hbar = 6.582119 * 10**-16 #ev.s

###################################################
##				READ THE DATA
###################################################

# read the energies
data = np.loadtxt('energies.dat')
data_el,data_hole,dos_hole,dos_el,Eocc,Evirt, = determine_cm_accessibe_v2(data)
Egap = Evirt[0]-Eocc[-1]

print len(data_el)
print len(data_hole)

plt.semilogy(dos_el[0,:]/Egap,dos_el[1,:])
plt.semilogy(dos_hole[0,:]/Egap,dos_hole[1,:])
plt.show()

data_save_el = np.array([dos_el[0,:]/Egap,dos_el[1,:]])
data_save_hole = np.array([dos_hole[0,:]/Egap,dos_hole[1,:]])
np.savetxt('dos_el.dat',data_save_el.T)
np.savetxt('dos_hl.dat',data_save_hole.T)



# read the orbitals
eigenvectors = np.loadtxt('eigenvectors.dat')

# read the system
system = atomic_system('Pbse.in')


###################################################
##				Intra atomic coulomb integrals
##		COMPUTED WITH EPS = 1
###################################################
se_U_ss = 24.074874  
se_U_sp = 14.307149 
se_U_pp = 12.881552 
 

pb_U_ss = 9.642478 
pb_U_sp = 7.085438
pb_U_pp = 6.648948


###################################################
##		Form the itra atomic sub matrices
###################################################
v_diag_se = np.array([
[se_U_ss, se_U_sp, se_U_sp, se_U_sp],
[se_U_sp,se_U_pp, se_U_pp, se_U_pp],
[se_U_sp,se_U_pp, se_U_pp, se_U_pp],
[se_U_sp,se_U_pp, se_U_pp, se_U_pp] ])

v_diag_pb = np.array([
[pb_U_ss, pb_U_sp, pb_U_sp, pb_U_sp],
[pb_U_sp,pb_U_pp, pb_U_pp, pb_U_pp],
[pb_U_sp,pb_U_pp, pb_U_pp, pb_U_pp],
[pb_U_sp,pb_U_pp, pb_U_pp, pb_U_pp] ])

###################################################
##				COMPUTE THE V MATRIX
###################################################
norb = 4*system.nAtom
V = np.zeros((norb,norb))

# for all the atoms
for iA in range(system.nAtom):

	# atom type
	at_typ = system.atomType[iA]
	index_iA = range(4*iA,4*(iA+1))

	# for all tge atoms
	for jA in range(iA,system.nAtom):

		index_jA = range(4*jA,4*(jA+1))

		# if we are on the same atom
		# we have to put the intra atomic Coulomb integrals
		if iA == jA:

			if at_typ == 'Se':
				V[np.ix_(index_iA,index_iA)] = v_diag_se/epsilon_

			elif at_typ == 'Pb':
				V[np.ix_(index_iA,index_iA)] = v_diag_pb/epsilon_

		else:

			rij = np.sqrt(np.sum((system.xyz[iA]-system.xyz[jA])**2))/0.529 # in bohr
			vij = (1./rij)*27.211	# eV
			V[np.ix_(index_iA,index_jA)] = vij*np.ones((4,4))/epsilon
			V[np.ix_(index_jA,index_iA)] = vij*np.ones((4,4))/epsilon
print V

g2 = (0.005)**2

##########################################
##	Compute the CM coupling
##########################################
nT_EL = len(data_el)
Wph = 0.017
maxEx = Evirt[-1]-Eocc[-1]
X_el = dos_el[0,:]
DOS_el = np.zeros(len(dos_el[0,:]))
data_CM_el = []

for iT in range(nT_EL):

	print '===== %d/%d' %(iT,nT_EL)
	
	# index of the initia/final states
	i1 = int(data_el[iT,1])
	i2 = int(data_el[iT,3])
	f1 = int(data_el[iT,2])
	f2 = int(data_el[iT,4])
	print '\t Eex = %1.3f eV' %(data_el[iT,0])
	print '\t I1 = %d I2 = %d' %(i1,i2)
	print '\t F1 = %d F2 = %d' %(f1,f2)
	Eex = data_el[iT,0]
	de = data_el[iT,-1]
	Vcm = compute_cm_coupling_baer(i1,i2,f1,f2,eigenvectors,V,system.nAtom)
	data_CM_el.append([Eex,Vcm])
	DOS_el += 2*np.pi/hbar*Vcm*Wph**2/((X_el-Eex)**2+Wph**2)
	print '\t Vcm = %1.3e delta = %1.3e eV' %(Vcm,de)
data_CM_el = np.array(data_CM_el)

##########################################
##	Compute the CM coupling from home
##########################################
nT_HL = len(data_hole)
Wph = 0.017
maxEx = Evirt[0]-Eocc[0]
X_hl = dos_hole[0,:]
DOS_hl = np.zeros(len(dos_hole[0,:]))
data_CM_hl = []
for iT in range(nT_HL):

	print '===== %d/%d' %(iT,nT_HL)
	
	# index of the initia/final states
	i1 = int(data_hole[iT,1])
	i2 = int(data_hole[iT,3])
	f1 = int(data_hole[iT,2])
	f2 = int(data_hole[iT,4])
	print '\t Eex = %1.3f eV' %(data_hole[iT,0])
	print '\t I1 = %d I2 = %d' %(i1,i2)
	print '\t F1 = %d F2 = %d' %(f1,f2)

	Eex = data_hole[iT,0]
	de = data_hole[iT,-1]
	Vcm = compute_cm_coupling_baer(i1,i2,f1,f2,eigenvectors,V,system.nAtom)
	data_CM_hl.append([Eex,Vcm])
	DOS_hl += 2*np.pi/hbar*Vcm*Wph**2/((X_hl-Eex)**2+Wph**2)
	print '\t Vcm = %1.3e delta = %1.3e eV' %(Vcm,de)

data_CM_hl = np.array(data_CM_hl)

#np.savetxt('data_cm_cplg_el.dat',data_CM_el)
#np.savetxt('data_cm_cplg_hl.dat',data_CM_hl)


#plt.plot(data_CM_el[:,0]/Egap,data_CM_el[:,1],'o')
#plt.plot(data_CM_hl[:,0]/Egap,data_CM_hl[:,1],'o')
#plt.show()


plt.plot(X_hl/Egap,DOS_hl*dos_hole[1,:])
plt.plot(X_el/Egap,DOS_el*dos_el[1,:])
plt.show()

data_save = np.array([X_hl/Egap,DOS_hl*dos_hole[1,:]])
np.savetxt('data_cm_cplg_hl.gnuplot',data_save.T)

data_save = np.array([X_el/Egap,DOS_el*dos_el[1,:]])
np.savetxt('data_cm_cplg_el.gnuplot',data_save.T)
