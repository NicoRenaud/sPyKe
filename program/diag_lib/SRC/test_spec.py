import spec_module as spec
import numpy as np
import ctypes
import scipy.linalg as scila

# show what's inside the module
dir(spec)

# dim
n=2

# create a random H
h = np.random.rand(n,n)
h = h+h.T

h = np.zeros((2,2))
h[0,1] = 1
h[1,0] = 1

# create a S 
s = np.eye(n)
s[0,1] = 0.1
s[1,0] = 0.1

print s
# crete void containers
eig_val = np.zeros(n)
eig_vect = np.zeros((n,n))

# diagonalize

spec.spec_pencil_zinger(eig_val, eig_vect, h, s, n)
eig_vect = eig_vect.T

print "Eigen value with C"
print eig_val

print "Eigen vectors with C"
print eig_vect


print s

# diagonalize with scipy
w2, u2 = scila.eigh(a=h,eigvals=[0,1],b=s)

print "Eigen value with scipy"
print w2

print "Eigen vector with scipy"
print u2