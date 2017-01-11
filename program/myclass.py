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
# read the spin orbit infrmation
##############################################
def read_spinOrbit(filename):

		spin_orbit = []
		SO_FLAG = 0
		#open the input file, read the info, close the file
		inputFile1=open(filename,'r')
		inputFileLines=inputFile1.readlines()
		inputFile1.close()

		#remove any newlines at the end of the input file
		lenInput=len(inputFileLines)
		while inputFileLines[-1]=='\n':
			del inputFileLines[-1]
			lenInput=lenInput-1

		#make a matrix of the atom information 
		for x in range(lenInput):
			if inputFileLines[x].split("\n")[0]== "SPIN_ORBIT":
				SO_FLAG = 1
				while inputFileLines[x+1].split("\n")[0] != "END":
					spin_orbit.append(inputFileLines[x+1].split())
					x+=1
		return SO_FLAG, spin_orbit










