import os.path as path
import numpy as np
from loadbinary import importbinary_old

H = importbinary_old('smallBinary.obj')

dim = 1024

filename = path.join('c:\\temp\\','newfile.dat')

Hreal = H[:,0]
Himag = H[:,1]

H = Hreal + 1j*Himag

H = H.reshape(dim,dim)

fp = np.memmap(filename, dtype='complex', mode='w+', shape=(dim,dim))

fp[:] = H[:]

del H

print np.dot(fp,fp)
