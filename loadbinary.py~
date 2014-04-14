import numpy as np
import time


def importbinary_old(path):
    start = time.time()
    fd = open(path, 'rb')

    try:
        dim = np.fromfile(fd, dtype=np.dtype('i'), count=1)
        H = np.fromfile(fd, dtype=np.dtype('d'),
                count=2*dim*dim).reshape((2, dim*dim)).T
        res_states = np.fromfile(fd, dtype=np.dtype('i'), count=2)
        check = np.fromfile(fd, dtype=np.dtype('l'), count=1)

        print check
    finally:
        fd.close()
    return H

'''
def importbinary(path):
    start = time.time()
    fd = open(path, 'rb')

    try:
        dim = np.fromfile(fd, dtype=np.dtype('i'), count=1)
        for i in range(dim):

        H = np.fromfile(fd, dtype=np.dtype('d'),
                count=2*dim*dim).reshape((2, dim*dim)).T

        res_states = np.fromfile(fd, dtype=np.dtype('i'), count=2)
        check = np.fromfile(fd, dtype=np.dtype('l'), count=1)

        print check
    finally:
        fd.close()
'''

if __name__ == '__main__':
    print 10*'-'
    importbinary_old('smallBinary.obj')
