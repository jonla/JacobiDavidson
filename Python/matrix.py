from __future__ import division
import numpy as np
# from nhqm.bases.momentum import MomentumBasis, gauss_contour, triangle_contour
# from nhqm.problems import Helium5


def realsymmetric(k, D):
    ''' Creates a k by k sized real symmetric matrix
    with D non-zero elements on each row'''
    A = np.zeros((k, k))

    for i in xrange(k):
        for j in xrange(i + 1):
            if np.abs(i - j) < D:
                A[i, j] = 2 * np.random.rand(1) - 1
                A[j, i] = A[i, j]
    return A


def complexhermitian(k, D):
    ''' Creates a k by k sized symmetric complex matrix
    with D non-zero elements on each row'''
    A = np.zeros((k, k), dtype=complex)

    for i in xrange(k):
        A[i, i] = 2 * np.random.rand(1) - 1
        for j in xrange(i):
            if np.abs(i - j) < D:
                A[i, j:j + 1] = (2 * np.random.rand(1) - 1
                        + 2 * np.random.rand(1) * 1j + 1j)
                A[j, i] = np.conj(A[i, j])
    return A


def complexsymmetric(k, D):
    ''' Creates a k by k sized symmetric complex matrix
    with D non-zero elements on each row'''
    A = np.zeros((k, k), dtype=complex)

    for i in xrange(k):
        for j in xrange(i + 1):
            if np.abs(i - j) < D:
                A[i, j:j + 1] = (2 * np.random.rand(1) - 1
                        + 2 * np.random.rand(1) * 1j + 1j)
                A[j, i] = A[i, j]
    return A


def matrixgen(k, D):
    ''' Creates a k by k sized symmetric matrix
    with D non-zero elements on each row. This
    old method is left here just to not break stuff
    '''
    A = np.zeros((k, k))

    for i in xrange(k):
        for j in xrange(i + 1):
            if np.abs(i - j) < D:
                A[i, j] = np.random.rand(1)
                A[j, i] = A[i, j]
    return A


def helium5berggrenmatrix(x_peak=0.17, y_peak=-0.07, basis_state_count=30,
        k_max=30):
    problem = Helium5()

    # Bases and contours
    berg_contour = triangle_contour(x_peak, y_peak, k_max, basis_state_count, 5)
    berg = MomentumBasis(berg_contour)

    return berg.hamiltonian(problem, quantum_numbers)
