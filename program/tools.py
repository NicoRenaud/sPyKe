import numpy as np
from string import rjust
from string import ljust
import commands
import os
import sys
import copy
import hkl_module as hkl
import ctypes

################################################
#
#	 Count the number of electrons
#
################################################
def count_electron(cell,hkl_param):

		param_raw = []
		atom_type = []

		#open the input file, read the info, close the file
		inputFile1=open(hkl_param,'r')
		inputFileLines=inputFile1.readlines()
		inputFile1.close()

		#make a matrix of the atom information
		lenInput=len(inputFileLines)
		while inputFileLines[-1]=='\n':
			del inputFileLines[-1]
			lenInput=lenInput-1
		
		for x in range(lenInput):
					param_raw.append(inputFileLines[x].split())
					atom_type.append([param_raw[x][0] , param_raw[x][1]])

		nelec = 0
		for i in range(cell.nAtom):
			for j in range(len(atom_type)):
				if cell.atomType[i] == atom_type[j][0]:
					nelec += int(atom_type[j][1])

		return nelec

################################################
#
#	 Block diagonalize  Hcell Hint
#											Hint  Hcell
################################################
def block_diag(hcell,scell,hint,sint):

	
	# diagonalize the cell hamiltonian
	wcell,ucell = np.linalg.eigh(hcell)
	norb = len(hcell)
	
	# index of the two cell in the total matrices
	ind_cell1 = range(norb)
	ind_cell2 = range(norb,2*norb)
	
	# form the oration matrix
	Urot = np.zeros((2*norb,2*norb))
	Urot[np.ix_(ind_cell1,ind_cell1)] = ucell
	Urot[np.ix_(ind_cell2,ind_cell2)] = ucell



	# form the total Hamiltonian
	Htot = np.zeros((2*norb,2*norb))
	Htot[np.ix_(ind_cell1,ind_cell1)] = hcell
	Htot[np.ix_(ind_cell2,ind_cell2)] = hcell
	Htot[np.ix_(ind_cell1,ind_cell2)] = hint
	Htot[np.ix_(ind_cell2,ind_cell1)] = hint.T

	# form the total overlap matrnp.ix
	Stot = np.zeros((2*norb,2*norb))
	Stot[np.ix_(ind_cell1,ind_cell1)] = scell
	Stot[np.ix_(ind_cell2,ind_cell2)] = scell
	Stot[np.ix_(ind_cell1,ind_cell2)] = sint
	Stot[np.ix_(ind_cell2,ind_cell1)] = sint.T

	# rotate the matrices
	Htot = np.dot(Urot.T,np.dot(Htot,Urot))
	Stot = np.dot(Urot.T,np.dot(Stot,Urot))

	# extract the new matrices
	hcell = Htot[np.ix_(ind_cell1,ind_cell1)]
	hint  = Htot[np.ix_(ind_cell1,ind_cell2)]
	scell = Stot[np.ix_(ind_cell1,ind_cell1)]
	sint  = Stot[np.ix_(ind_cell1,ind_cell2)]


	return hcell,scell,hint,sint

################################################
#
#	Orthogonalize the basis
#
################################################


def ortho(H,S):

	n = S.shape[0]
	s, U = eigh(S)
	V = dot(dot(U, diag(1 / sqrt(s))), U.T)
	Ho = dot(dot(V, H), V)

	return Ho

################################################
#
#	Compute the electronic struture
#
################################################

def compute_h(system,hkl_param,keht,out):

	# name f the input file for huckel
	hkl_in = out+"/_hkl_.in"


	# create a input file for huckel
	export_huckel(hkl_in,system,hkl_param,keht)


	# determine the number of orbitals
	nb_orb = hkl.nbOrb(hkl_in)

	# create the matrix
	h = np.zeros(nb_orb*nb_orb)
	s = np.zeros(nb_orb*nb_orb)

	# create the matrices
	hkl.hklHam(h,s,nb_orb,hkl_in)

	# reshape
	h = h.reshape(nb_orb,nb_orb)
	s = s.reshape(nb_orb,nb_orb)

	# execute huckel
	#cmd = huckel_exe + ' ' + hkl_in + ' ' + out + '/'
	#os.system(cmd)

	# import the results in pyhton
	#h = iter_loadtxt(out+"/H.dat")
	#s = iter_loadtxt(out+"/S.dat")

	return h, s

################################################
#
#       load text file into numpy
#
################################################

def iter_loadtxt(filename, skiprows=0, dtype=float):
    def iter_func():
        with open(filename, 'r') as infile:
            for _ in range(skiprows):
                next(infile)
            for line in infile:
                line = line.rstrip().split()
                for item in line:
					yield dtype(item)
        iter_loadtxt.rowlength = len(line)

    data = np.fromiter(iter_func(), dtype=dtype)
    data = data.reshape((-1, iter_loadtxt.rowlength))
    return data

################################################
#
#	Export the husky file for the junction
#
################################################

def export_huckel(huckelFileName,system,hkl_param,keht):

	natom = len(system.xyz)

	fp = open(huckelFileName,'w')
	fp.write(" //////////////////////////////\n //   pos for Huckel \n //////////////////////////////\n")
	fp.write("nb_atom %d\n" %(natom) )
	fp.write("parameters %s\n" %(hkl_param))
	fp.write("Keht %1.3f\n" %(keht))

	for x in range(natom):
			fp.write("%s %s %s %s\n"  %(ljust(str(system.atomType[x]),3), rjust(str(system.xyz[x][0]),20), rjust(str(system.xyz[x][1]),20), rjust(str(system.xyz[x][2]),20)) )

	fp.write("\n//////////////////////////////\n")
	fp.close


################################################
#
#	Export the eigenvector files
#
################################################

def write_vect(data,iorb,Ef,K,fname):

	# number of K points
	nKpoint = len(data)
	
	#  total number of orbitals
	norb = len(data[0])-1
	
	# number of band to be pritned
	nBand = len(iorb)
	
	# energy scale
	nE = 250
	E = np.linspace(-15,0,nE)
	
	# sigma for the Gaussian
	sigma=0.1
	
	
	f = open(fname,'w')
	gmax = 0
	
	# for all the Kpoints
	for iK in range(nKpoint):
	
		# energy of the band
		e0 = data[iK][0] #-Ef
		
		# for all the energy in the scale
		for iE in range(nE):
			c = 0
			
			# sum up all the coefficients +1 because [0] is the energy of the band
			for ip in range(len(iorb)):
				c += data[iK][iorb[ip]+1]
			
			# compute the Gaussian
			g = c * np.exp( -(E[iE]-e0)**2/(sigma*c) )
			
			if g>gmax:
				gmax = g
			
			# print in the file
			f.write('%f %f %f\n' %(K[iK],E[iE],g))
		f.write('\n')
	f.close()
	return gmax



################################################
#
#	Get the atomic orbital
#
################################################
def get_atomicOrbital(filename,atoms):

	f = open(filename,'r')
	data = f.readlines()
	f.close

	ao = []
	for i in range(len(data)):
		data[i] = data[i].split()

	index_mat = 0
	# for all the atoms in the cell
	for i in range(len(atoms)):

		# current atom type
		at = atoms[i]

		# check what type it is
		for j in range(len(data)):

			if len(data[j])>0:
				if at == data[j][0]:

					orb = []

					# number of orbtital for each
					nS = int(data[j][2])
					nP = int(data[j][3])
					nD = int(data[j][4])

					if nS != 0:

						orb.append(['s',index_mat])
						index_mat += 1

					if nP != 0:

						orb.append(['px',index_mat])
						index_mat += 1

						orb.append(['py',index_mat])
						index_mat += 1

						orb.append(['pz',index_mat])
						index_mat += 1

					if nD != 0:

						for iorb in range(5):
							oname = 'd%d' %(iorb+1)
							orb.append([oname,index_mat])
							index_mat += 1
		ao.append(orb)
	return(ao)	

################################################
#
#	Apply the spin orbit coupling
#
################################################
def spinOrbit(hcell,ao,spin_orbit,atomType):

	n = len(hcell)
	I = np.eye(2)
	hcell_so = np.kron(I,np.matrix.copy(hcell)).astype('complex')
	

	for iat in range(len(ao)):

		for jat in range(len(spin_orbit)):
			if atomType[iat] == spin_orbit[jat][0]:
				so = float(spin_orbit[jat][1])

		for iorb in range(len(ao[iat])):

			if ao[iat][iorb][0]=='px':

				# index of the orbitals concerned
				index_px_p = ao[iat][iorb][1]
				index_py_p = ao[iat][iorb][1]+1
				index_pz_p = ao[iat][iorb][1]+2

				index_px_m = ao[iat][iorb][1]+n
				index_py_m = ao[iat][iorb][1]+1+n
				index_pz_m = ao[iat][iorb][1]+2+n

				# < p_x +| Hso| p_z - >
				hcell_so[index_px_p,index_pz_m] += so
				hcell_so[index_pz_m,index_px_p] += so

				# < p_x -| Hso| p_z + >
				hcell_so[index_px_m,index_pz_p] -= so
				hcell_so[index_pz_p,index_px_m] -= so

				# <p_x +| Hso| p_y +>
				hcell_so[index_px_p,index_py_p] -= 1.0j*so
				hcell_so[index_py_p,index_px_p] += 1.0j*so

				# <p_x -| Hso| p_y ->
				hcell_so[index_px_m,index_py_m] += 1.0j*so
				hcell_so[index_py_m,index_px_m] -= 1.0j*so


				# < p_y +| Hso| p_z - >
				hcell_so[index_py_p,index_pz_m] -= 1.0j*so
				hcell_so[index_pz_m,index_py_p] += 1.0j*so

				# < p_y +- Hso| p_z + >
				hcell_so[index_py_m,index_pz_p] -= 1.0j*so
				hcell_so[index_pz_p,index_py_m] += 1.0j*so


	return(hcell_so)


################################################
#
#	print a mopac file for MO
#
################################################
def print_mo_mopac(wH,uH,s,indexHomo,system,filename,PARAM,SO,nb_print_orb):

	natom = system.nAtom
	natom_print = natom

	if SO:
		natom_print*=2

	#nb_orb = len(wH)
	nb_orb = len(uH)
	count = 0

	# comute the inverse of S
	invS = np.linalg.inv(s)


	# open the output file
	f = open(filename,'w')

	# header
	f.write("        %d MOPAC-Graphical data Version 2012.13.084W\n" %natom_print)

	# print the atoms
	for i in range(natom):
		at = find_index(system.atomType[i])
		f.write('%4d %*f%*f%*f  0.0000\n' %(at,12,system.xyz[i][0],12,system.xyz[i][1],12,system.xyz[i][2]))

	#if we have spin orbit we double the atoms
	if SO:
		for i in range(natom):
			at = find_index(system.atomType[i])
			f.write('%4d %*f%*f%*f  0.0000\n' %(at,12,system.xyz[i][0],12,system.xyz[i][1],12,system.xyz[i][2]))



	# print the slater exponents
	for i in range(natom):
		sc1,sc2,sc3 = get_slater_coeff(system.atomType[i],PARAM)
		f.write("  %1.7f  %1.7f  %1.7f\n" %(sc1,sc2,sc3))

	#if we have spin orbit we double the coeff
	if SO:
		# print the slater exponents
		for i in range(natom):
			sc1,sc2,sc3 = get_slater_coeff(system.atomType[i],PARAM)
			f.write("  %1.7f  %1.7f  %1.7f\n" %(sc1,sc2,sc3))

	# determine what to print
	index_print_orb = []
	occ = []

	for i in range(nb_print_orb):
		index_print_orb.append(indexHomo-nb_print_orb/2+1+i)
		if i<nb_print_orb/2:
			occ.append(2)
		else:
			occ.append(0)

	# print the orbitals
	for iprint in range(nb_print_orb):

		iorb = index_print_orb[iprint]
		f.write(" ORBITAL %3d  A  %2.5g\n" %(occ[iprint],wH[iorb]-wH[indexHomo]))

		for jcomp in range(nb_orb):
			f.write("% 1.8E" %(uH[jcomp,iorb]).real)
			count += 1
			if count == 5:
				f.write("\n")
				count  = 0
		if count>0:
				f.write("\n")
				count  = 0

	# print the inverse matrix
	count = 0
	f.write("INVERSE_MATRIX[%dx%d]=\n"%(nb_orb,nb_orb))
	for i in range(nb_orb):
		for j in range(i+1):
			f.write("% 1.8E" %(invS[j,i]))
			count+=1
			if count == 5:
				f.write("\n")
				count  = 0
		if count>0:
				f.write("\n")
				count  = 0
	f.write(" Keywords: SYMMETRY GRAPHF")
	f.close()


###########################################################
#
# print the energy level
#
###########################################################
def print_level(w,index_homo,outdir):

	filename = outdir+'/level_occ.dat'
	f = open(filename,'w')
	for i in range(index_homo):
		w0 = w[i]-w[index_homo]
		f.write("%f %f \n %f %f\n\n" %(-0.5,w0,0.5,w0))
	f.close()

	filename = outdir+'/level_virt.dat'
	f = open(filename,'w')
	for i in range(index_homo+1,len(w)):
		w0 = w[i]-w[index_homo]
		f.write("%f %f \n %f %f\n\n" %(-0.5,w0,0.5,w0))
	f.close()




###########################################################
#
# print the DOS
#
###########################################################
def print_dos(w,index_homo,filename,width,emin,emax,nE):

	f = open(filename,'w')
	E = np.linspace(emin,emax,nE)
	for iE in range(len(E)):
		g = 0
		e = E[iE]
		for iW in range(len(w)):
			w0 = w[iW]-w[index_homo]
			g += np.exp( -(e-w0)**2/width**2 )
		f.write("%f %f\n" %(e,g))
	f.close()

###########################################################
#
# find the index that corresponds to a given atom name
#
###########################################################
def find_index(name):
	#list = ["H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Ti", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr"]
	# replaced Fe by Pb otherwise Jmol can't print MO on Pb ....
	list = ["H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Pb", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Ti", "Fe", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr"]
	index= [i for i, x in enumerate(list) if x == name]
	return index[0]+1



###########################################################
#
# find the slater expoants of a given atom
#
###########################################################
def get_slater_coeff(name,param):

	f =open(param,'r')
	data = f.readlines()
	f.close

	for i in range(len(data)):
		data[i] = data[i].split()
		if data[i][0] == name:
			sc1 = data[i][7]
			sc2 = data[i][13]
			sc3 = data[i][19]
			break
	return float(sc1),float(sc2),float(sc3)



###########################################################
#
# Get the number of atomic orbitals per atoms
#
###########################################################
def get_nbAO(atomType,param):

	nbAo = []

	f =open(param,'r')
	data = f.readlines()
	f.close

	for i in range(len(data)):
		data[i] = data[i].split()

	for i in range(len(atomType)):

		for j in range(len(data)):
			if data[j][0] == atomType[i]:
				if int(data[j][2]) != 0:
					nbAo.append([i])
				if int(data[j][3]) != 0:
					nbAo.append([i,i,i])
				if int(data[j][4]) != 0:
					print 'd orbitals detected '
					nbAo.append([i,i,i,i,i])
				break;
	nbAo = sum(nbAo,[])
	return nbAo


###########################################################
#
# Compute the optical properties
#
###########################################################
def compute_optical(w,u,s,indexHomo,pos,nbAO,ntrans,filename,_SO_,_ZDO_):

	if _ZDO_ :
		print "\n \t\t\t       Neglect of the overalp during calculation",
		sys.stdout.flush()

	f = open(filename,'w')

	# if not SO everything normal
	size = len(u)
	MULT = 1

	# if SO we consider only the spin up MO
	if _SO_:
		size = len(u)/2
		MULT = 2

	# constant for conversion of dipoles
	l0 = 1.085*0.00001/0.00012

	# create the output
	fTrans = []
	mu = []

	# for all the occ MO we consider
	for i in range(ntrans):

		# we skip the spin down MO --> 2*i
		index_occ = indexHomo - MULT*i

		for j in range(ntrans):

			# we also skip the spin down MO --> 2*j
			index_virt = indexHomo + 1 + MULT*j

			# frequency
			freq = w[index_virt].real-w[index_occ].real

			dx = 0
			dy = 0
			dz = 0

			# here we neglect the overlap
			if _ZDO_:
				for k1 in range(size):
					dx += u[k1,index_occ].real*u[k1,index_virt].real*pos[nbAO[k1]][0]
					dy += u[k1,index_occ].real*u[k1,index_virt].real*pos[nbAO[k1]][1]
					dz += u[k1,index_occ].real*u[k1,index_virt].real*pos[nbAO[k1]][2]

			else:
				for k1 in range(size):
					for k2 in range(size):
						dx += u[k1,index_occ].real*u[k2,index_virt].real*pos[nbAO[k2]][0]*s[k1,k2]
						dy += u[k1,index_occ].real*u[k2,index_virt].real*pos[nbAO[k2]][1]*s[k1,k2]
						dz += u[k1,index_occ].real*u[k2,index_virt].real*pos[nbAO[k2]][2]*s[k1,k2]

			DX = dx*dx*2*l0/freq
			DY = dy*dy*2*l0/freq
			DZ = dz*dz*2*l0/freq

			DTOT = DX+DY+DZ
			f.write("%02d %02d % f % f % f % f % f\n" %(i,j,freq,DX,DY,DZ,DTOT))
			#print("%02d %02d % f % f % f % f % f" %(i,j,freq,DX,DY,DZ,DTOT))

			fTrans.append(freq)
			mu.append(DTOT)

	f.close()
	return fTrans,mu

###########################################################
#
# print the DOS
#
###########################################################
def print_abs(w,mu,width,filename):
	f = open(filename,'w')

	emin = w[0]*0.9
	emax = w[-1]*1.1
	nE = 1000


	E = np.linspace(emin,emax,nE)
	for iE in range(len(E)):
		g = 0
		e = E[iE]
		enm = 1240./e
		for iW in range(len(w)):
			g += mu[iW]*np.exp( -(e-w[iW])**2/width**2 )
		f.write("%f %f %f\n" %(e,enm,g))
	f.close()







		
