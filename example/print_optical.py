import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('optical.dat')

k = 0
M = np.zeros((10,10))
for i in range(10):
	for j in range(10):
		M[i,j] = data[k,6]
		k+= 1
print M
plt.imshow(M)
plt.show()