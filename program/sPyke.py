import numpy as np
import scipy.linalg as scila
import os
import sys
import argparse
from myclass import *
from tools import *
import datetime as dt
import time
import hkl_module as hkl
import spec_module as spec


#sys.stdout = open('sPyKe.out', 'w')


def main(argv):


  	#huckel path laptop
	#huckel_exe='/Users/nicolasrenaud/Documents/PROJECTS/PbSe_Frank/sPyKe/program/huckel/huckel'
		
	# which methof we want to diagonalize the matrix
	__diag__ = 'lapack'

	print "\n\n"
	print "%s ====================================================" %dt.datetime.now()
	print "%s ==         sPyKe :                                ==" %dt.datetime.now()
	print "%s ==         Extended Huckel                        ==" %dt.datetime.now()
	print "%s ==         Calculation in python                  ==" %dt.datetime.now()
	print "%s ==================================================== \n\n" %dt.datetime.now()

	##########################
	#
	# 	Parse the arguments
	#
	##########################
	
	parser = argparse.ArgumentParser(description='Compute the band structure of periodic crystals')
	
	#required (positional) arguments
	parser.add_argument('input', help = 'input file for the calculation')
	parser.add_argument('hkl_param', help = 'EHMO parameter for the calculations')
		
	#optional arguments
	parser.add_argument('-odir', '--output_directory', default='./',
											help = 'output directory where the file will be stored')

	parser.add_argument('-keht', '--KEHT', default=1.75,
											help = 'Keht constant for the calculation of the couplings',type=float)


	########## paramter to compute the MO
	parser.add_argument('-cMO', '--compute_MO', default=100,
								help = 'number of MO to compute',type=float)


	######### parameter for the dos calculation
	parser.add_argument('-dw', '--dos_width', default=0.1,
							help = 'peak width for the dos',type=float)

	parser.add_argument('-dMin', '--dos_min', default=-2,
								help = 'minimum energy for the DOS calculation',type=float)

	parser.add_argument('-dMax', '--dos_max', default=3,
								help = 'maximum energy for the DOS calculation',type=float)

	parser.add_argument('-nDos', '--nb_dos', default=251,
								help = 'number of point in the dos',type=int)


	########### parameter to print the molecular orbital
	parser.add_argument('-pMO', '--print_MO', default=10,
								help = 'number of MO to print in mo.dat',type=int)

	######### parameter for the optical properties calculations
	parser.add_argument('-nTrans', '--nb_opt_trans', default=10,
								help = 'number of MO above/below the HOMO for the transition',type=int)

	parser.add_argument('-zdo', '--Zero_Diff_overlap', default=1,
								help = 'neglect (1) or not (0) the overlap during absorption calculation (default 1)',type=int)

	parser.add_argument('-opw', '--opticalPeak_width', default=0.05,
								help = 'width of the peak in the absorption spectra',type=float)

	# done
	args=parser.parse_args()
	out_dir = args.output_directory
	dos_width = args.dos_width
	dos_min = args.dos_min
	dos_max = args.dos_max
	dos_nE = args.nb_dos
	nMO_print = args.print_MO
	nOptTrans = args.nb_opt_trans
	_ZDO_ = args.Zero_Diff_overlap
	optPeak_width = args.opticalPeak_width
	nMO_compute = args.compute_MO

	# check
	even_test = (nMO_compute/2 - int(nMO_compute/2))==0
	if not even_test:
		print "odd number of MO to compute"
		nMO_compute = nMO_compute+1

	##########################
	#
	# 	Read the data
	#
	##########################
	
	print "%s === Read the input file" %dt.datetime.now(),
	sys.stdout.flush()

	## read the position of the unit cell
	system = atomic_system(args.input)

	# get the type of ao and their indexes
	ao = get_atomicOrbital(args.hkl_param,system.atomType)

	## determine the number of electrons
	nelec = count_electron(system,args.hkl_param)

	## read the spin orbit information 
	_SO_, spinorbit = read_spinOrbit(args.input)
	
	## index of homo/lumo
	index_homo = nelec/2-1
	if _SO_:
		index_homo=nelec-1

	print "\r\t\t\t\t\t\t\t\t\t\t done"

	# compute the electronic structure
	# extended huckel done with huckel code
  	print "%s === Compute the Hamiltonian with EHT" %dt.datetime.now(),
	sys.stdout.flush()

	h, s = compute_h(system,args.hkl_param,args.KEHT,out_dir)
	print "\r\t\t\t\t\t\t\t\t\t\t done"

	# size pf the cell
	norb = len(h)
	print "\t\t\t     System contains:"
  	print "\t\t\t\t --  %d atoms \n\t\t\t\t --  %d orbitals \n\t\t\t\t --  %d electrons" %(system.nAtom,norb,nelec)

	###############################################
	#
	#		Add the spin orbit interactions
	#
	###############################################
	if(_SO_):
		print "%s === Add the spin orbit coupling to the Hamiltonian" %dt.datetime.now(),
		sys.stdout.flush()

		h = spinOrbit(h,ao,spinorbit,system.atomType)
		s = np.kron(np.eye(2),s)
		norb *= 2
		nMO_print *= 2
		nMO_compute *= 2

	  	print "\r\t\t\t\t\t\t\t\t\t\t done"
	print "\t\t\t     System contains:"
  	print "\t\t\t\t --  %d atoms \n\t\t\t\t --  %d orbitals \n\t\t\t\t --  %d electrons" %(system.nAtom,norb,nelec)



	###############################################
	#
	#		Diagonalize the Hamiltonian
	#
	###############################################

	# range of eigenvalues to compute
	index_mo = [index_homo-nMO_compute/2+1,index_homo+nMO_compute/2]

	# check
	if norb < nMO_compute:
		print "\n%s === Warning : Too many orbitals to compute (%d / %d)" %(dt.datetime.now(),nMO_compute,norb)
		print "%s === Warning : All orbitals will be computed\n " %(dt.datetime.now())
		nMO_compute = norb
		index_mo = [0,norb-1]

	start = time.time()
	print "%s === Compute %d eigenvalues/vectors (%d --> %d)" %(dt.datetime.now(),nMO_compute,index_mo[0],index_mo[1]),
	sys.stdout.flush()


	# diagonalize the system using scipy
	if __diag__ == "scipy":

		# diagonalize
		w2, u2 = scila.eigh(a=h,eigvals=index_mo,b=s)

		#sort the eigenvalues
		ind_sort = np.argsort(w2)
		w2 = np.sort(w2)

		#deal with the eigenvectors
		u2 = u2[:,ind_sort]

	# diagonalize the system using cblas
	if __diag__ == 'lapack':

		# prepare for diag
		w2 = np.zeros(norb)
		u2 = np.zeros((norb,norb))

		# waring h and s are replaced by temp var
		# so we need to copy the overlap
		scopy = np.copy(s)

		# diagonalize
		spec.spec_pencil_zinger(w2, u2, h, scopy, norb)

		# transpose the eigen vectors
		u2=u2.T

	end = time.time()
	print "\r\t\t\t\t\t\t\t\t\t\t done in %1.6f s" %(end-start)


	# Normalize the eigenvectors
	for i in range(len(u2[0])):
		u2[:,i] /= np.linalg.norm(u2[:,i])

	index_homo_reduced = nMO_compute/2-1

	print "\t\t\t     Energy Information:"
	print "\t\t\t\t --  Energy HOMO (orb #%d): %f eV" %(index_homo,w2[index_homo_reduced])
	print "\t\t\t\t --  Energy LUMO (orb #%d): %f eV" %(index_homo+1,w2[index_homo_reduced+1])
	print "\t\t\t\t --  Energy Gap %f eV" %(w2[index_homo_reduced+1]-w2[index_homo_reduced])


	###############################################
	#
	#		Print the level/dos/mo
	#
	###############################################

	# print the energy levels
	print "%s === Compute the DOS" %dt.datetime.now(),
	sys.stdout.flush()
	print_level(w2,index_homo_reduced,out_dir)
	print_dos(w2,index_homo_reduced,out_dir+'/dos.dat',dos_width,dos_min,dos_max,dos_nE)
	print "\r\t\t\t\t\t\t\t\t\t\t done"


	if nMO_print > nMO_compute:
		print "\n%s === Warning : Too many orbitals to be printed (%d / %d)" %(dt.datetime.now(),nMO_print,nMO_compute)
		print "%s === Warning : All orbitals will be printed\n " %(dt.datetime.now())
		nMO_print = nMO_compute


	# print the mo
	print "%s === Compute the molecular orbitals" %dt.datetime.now(),
	sys.stdout.flush()
	print_mo_mopac(w2,u2,s,index_homo_reduced,system,out_dir+'/mo.dat',args.hkl_param,_SO_,nMO_print)
	print "\r\t\t\t\t\t\t\t\t\t\t done"



	###############################################
	#
	#		Compute the optical properties
	#
	###############################################

	if nOptTrans>nMO_compute:
		print "\n%s === Warning : Too many transitions demanded (%d / %d)" %(dt.datetime.now(),nOptTrans,nMO_compute)
		print "%s === Warning : All available transitions will be computed\n " %(dt.datetime.now())

	print "%s === Compute the optical absorption spectra" %dt.datetime.now(),
	sys.stdout.flush()
	nbAO = get_nbAO(system.atomType,args.hkl_param)
	if _SO_:
		nbAO = nbAO+nbAO
	ftrans,mu = compute_optical(w2,u2,s,index_homo_reduced,system.xyz,nbAO,nOptTrans,out_dir+'/optical.dat',_SO_,_ZDO_)
	print_abs(ftrans,mu,optPeak_width,out_dir+'/abs.dat')
	sys.stdout.write('\r\t\t\t\t\t\t\t\t\t\t done\n')

if __name__=='__main__':
	main(sys.argv[1:])




