from numpy import *
##############################################
# class of the repetition vector
# and dimensionality of the system
##############################################
class atomic_system:

	
	# define the molecule from the input files
	# file_name contaons the xyz

	def __init__(self,file_name):

		####################################################
		##
		## 	Read the XYZ Information
		##
		####################################################

		vect_raw = []

		#open the input file, read the info, close the file
		inputFile1=open(file_name,'r')
		inputFileLines=inputFile1.readlines()
		inputFile1.close()

		#remove any newlines at the end of the input file
		lenInput=len(inputFileLines)
		while inputFileLines[-1]=='\n':
			del inputFileLines[-1]
			lenInput=lenInput-1

		#make a matrix of the atom information (this happens quite a bit in this code - get used to it)
		for x in range(lenInput):
			if inputFileLines[x].split("\n")[0]== "ATOMS":
				while inputFileLines[x+1].split("\n")[0] != "END":
					vect_raw.append(inputFileLines[x+1].split())
					x+=1

		# number of atoms
		self.nAtom = int(len(vect_raw))

		# type of atoms
		self.atomType = chararray(self.nAtom,itemsize=2)
		for i in range(self.nAtom):
			self.atomType[i] = vect_raw[i][0]

		# atom positions	
		self.xyz = zeros((self.nAtom,3))	
		for i in range(self.nAtom):
			pos = [ float(vect_raw[i][1]) , float(vect_raw[i][2]),float(vect_raw[i][3])  ]
			self.xyz[i] = pos


##############################################
# determine the states accessible for CM
##############################################
import numpy as np
def determine_cm_accessibe(data):


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

	data_el = []
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
					data_el.append([Eex,L,nocc+ivirt,nocc+ivirt2,nocc-1,nocc+ivirt3,DE])
					k+=1

	maxEx = Evirt[-1]-EH
	X_el = np.linspace(0,maxEx*1.25,500)
	DOS_el = np.zeros(500)
	data_el = np.array(data_el)
	for i in range(len(data_el)):
		DOS_el += data_el[i,1]*0.010/((X_el-data_el[i,0])**2+0.010**2)
	dos_el = np.array([X_el/Egap,DOS_el])


	################################################
	# transition from an excited hole
	################################################
	data_hole = []
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
					data_hole.append([Eex,L,iocc,iocc2,iocc3,nocc+1,DE])
					k+=1

	maxEx = EL-Eocc[0]
	X_hole = np.linspace(0,maxEx*1.25,500)
	DOS_hole = np.zeros(500)
	data_hole = np.array(data_hole)
	for i in range(len(data_hole)):
		DOS_hole += data_hole[i,1]*0.010/((X_hole-data_hole[i,0])**2+0.010**2)
	dos_hole = np.array([X_hole/Egap,DOS_hole])

	return data_el,data_hole,dos_hole,dos_el


##############################################
# determine the states accessible for CM
##############################################
import numpy as np
def determine_cm_accessibe_v2(data):


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
	Wph = 0.020

	Wph_max = 2*Wph

	THRESHOLD = 0.005

	################################################
	# transition from an excited electron
	################################################
	maxEx = Evirt[-1]-EH
	X_el = np.linspace(0,maxEx*1.25,500)
	DOS_el = np.zeros(500)
	data_el = []

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
				if np.abs(DE) < THRESHOLD:
					data_el.append([Eex,nocc+ivirt,nocc+ivirt2,nocc-1,nocc+ivirt3,DE])
					DOS_el += Wph/((X_el-Eex)**2+Wph**2)
					k+=1


	data_el = np.array(data_el)
	dos_el = np.array([X_el,DOS_el])


	################################################
	# transition from an excited hole
	################################################
	data_hole = []
	maxEx = EL-Eocc[0]
	X_hole = np.linspace(0,maxEx*1.25,500)
	DOS_hole = np.zeros(500)
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
				if np.abs(DE) < THRESHOLD:
					DOS_hole += Wph/((X_hole-Eex)**2+Wph**2)			
					data_hole.append([Eex,iocc,iocc2,iocc3,nocc+1,DE])
					k+=1

	
	data_hole = np.array(data_hole)
	dos_hole = np.array([X_hole,DOS_hole])

	return data_el,data_hole,dos_hole,dos_el,Eocc,Evirt


########################################
##	Compute the CM coupling 
########################################

def compute_cm_coupling(I1,I2,F1,F2,U,V,nat):

	C1 = U[:,I1]*U[:,I2]
	C2 = U[:,F1]*U[:,F2]
	CM = np.outer(C1,C2)
	CM = np.sum(CM*V)

	return CM

def compute_cm_coupling_baer(I1,I2,F1,F2,U,V,nat):

	W2 = 0

	# first term
	C1 = U[:,I1]*U[:,I2]
	C2 = U[:,F1]*U[:,F2]
	V1 = np.outer(C1,C2)
	V1 = np.sum(V1*V)
	

	# second term
	C1 = U[:,I1]*U[:,F1]
	C2 = U[:,I2]*U[:,F2]
	V2 = np.outer(C1,C2)
	V2 = np.sum(V2*V)

	W2 = 2*(V1-V2)**2+V1**2+V2**2

	return W2










