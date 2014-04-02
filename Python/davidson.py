import numpy as np
from numpy import dot
# import matplotlib.pyplot as plt
# from largest import largestEig
from numpy.linalg import norm
from matrix import * # noqa
import time


def davidsolver(A, guess, iterations, N_gmres, eps):
    '''
    This function is meant to compute the eigenpair with eigenvector
    overlapping the inital guess.

    TODO:
        * Modify the correction equation to use a suitable
        solver (GMRES) with a preconditioner.
    '''
    # Timing stuff
    startsolver = time.time()
    timeCEQ = 0

    V, M, theta, r, iterations = davidinit(A, guess)

    # theta1 = theta
    n = len(A)
    u = V
    for m in range(iterations):

        ## CORRECTION EQUATION
        startCEQ = time.time()
        # t = solvecorrectioneq(A, u, theta, r, n)
        K = np.diag(np.diag(A)) - theta*np.identity(n)
        t = gmres(A, K, u, r, theta, N_gmres)
        timeCEQ = timeCEQ + time.time() - startCEQ

        ## MOD. GRAM-SCHMIDT
        vplus = modgramshmidt(t, V)

        ## CONSTRUCT M, EXPAND V
        M = np.vstack([np.hstack([M, dot(V.T, dot(A, vplus))]),
                       np.hstack([dot(vplus.T, dot(A, V)),
                                  dot(vplus.T, dot(A, vplus))])])

        V = np.hstack((V, vplus))

        ## SELECT EIGENPAIR OF M, CALCULATE u AND r
        evals, evecs = np.linalg.eig(M)
        # thetai = abs(evals - theta1).argmin()
        # thetai = abs(evals).argmax()
        thetai = abs(evecs[0, :]).argmax()

        theta = evals[thetai]
        s = evecs[:, [thetai]]
        u = dot(V, s)
        r = dot(A, u) - theta * u
        f = norm(r)
        print "Iteration:", m, " Theta:", theta, " norm(r):", f
        if f < eps:
            print "Total time:", (time.time() - startsolver)
            print "time on CEQ:", (timeCEQ / (time.time() - startsolver))
            return theta, u
    print "Total time:", (time.time() - startsolver)
    print "time on CEQ:", (timeCEQ / (time.time() - startsolver))
    return theta, u


def davidinit(A, guess):
    iterations = len(guess)
    V = guess / norm(guess)
    theta = dot(V.T, dot(A, V))
    M = theta
    r = dot(A, V) - theta * V
    return V, M, theta, r, iterations


def solvecorrectioneq(A, u, theta, r, n):
    '''
    This is a very costly solution
    and should be replaced when we study larger
    matrices. The method we should use relies on finding a suitable
    preconditioner to (A-theta*I) and then the iterativ method GMRES.
    '''
    # Davidson's suggested correction equation
    # t = np.linalg.solve(np.diag(A) * np.eye(n) - theta * np.eye(n), -r)

    # The JD correction equation
    K = dot(np.eye(n) - np.outer(u, u), dot(A - theta * np.eye(n),
                                            np.eye(n) - np.outer(u, u)))
    K = np.vstack([K, 100 * u.T])   # lstsq weight 100 for u * t = 0
    t = np.linalg.lstsq(K, np.vstack([-r, 0]))[0]
    return t


def modgramshmidt(tin, V, kappah=0.25):
    t = tin

    #if len(V.shape) == 1 or len(V[0,:]) == 1:
     #   t = t - dot(t, V) * V

    if len(V[1]) == 1:
        t = t - dot(t.T, V) * V

    else:
        for j in range(len(V.T)):
            t = t - dot(t.T, V[:, [j]]) * V[:, [j]]
        if norm(t) / norm(tin) < kappah:
            for j in range(len(V.T)):
                t = t - dot(t.T, V[:, [j]]) * V[:, [j]]
    return t / norm(t)


def gmres(A, K, u, r, theta, iter):

    # Parameters needed
    I = np.identity(len(u))  # Hoping this is saved as sparse?
    uhat = np.linalg.solve(K, u)
    mu = dot(u.T, uhat)
    rhat = np.linalg.solve(K, r)
    rtilde = rhat - dot(u.T, rhat) / mu * uhat

    # Initiate GMRES
    beta = norm(rtilde)
    V = -rtilde / beta
    W = np.zeros([len(u), 0])

    for i in range(iter):
        # print 'GMRES iteration:', i
        y = dot((A - theta * I), V[:, i])
        yhat = np.linalg.solve(K, y)
        # W[:, i] = yhat - dot(u.T, yhat) / mu * uhat
        W = np.append(W, yhat - dot(u.T, yhat) / mu * uhat, 1)
        H = np.zeros([i + 2, i + 1], dtype=complex)

        for j in range(i):
            H[j, i] = dot(W[:, i], V[:, j])
            W[:, j] = W[:, j] - H[j, i] * V[:, j]

        H[i + 1, i] = norm(W[:, i])
        b = np.append(np.array([[beta]]), np.zeros([i + 1, 1]), 0)
        y = np.linalg.lstsq(H, b)[0]
        V = np.append(V, W[:, [i]] / H[i + 1, i], 1)

    return dot(V[:, :-1], y)


def davidsontest(k=100, TOL=1.e-3, N=45, hermitian=False, real=False,
                 target=18, guessoffset=0.15):

    print 50 * '-'
    # k = 100                     # Matrix size
    # TOL = 1.e-3                 # Margin of error
    D = 1000                     # Diagonal shape
    N_gmres = 5
    # N = 45                      # Iterations

    if hermitian:
        A = realsymmetric(k, D) if real else complexhermitian(k, d)
    else:
        A = complexsymmetric(k, D)

    eig, vec = np.linalg.eig(A)
    Eig = np.sort(eig)
    # target = 18  # = eig.argmax()
    eigtarget = eig[target]
    vtarget = vec[:, [target]]

    # guess = np.random.rand(k, 1)
    # guess = np.ones((k, 1))

    guess = vtarget + guessoffset * np.ones((k, 1))
    guess = guess / norm(guess)
    theta1 = dot(guess.T, dot(A, guess))
    print "Matrix size:", k
    print "Target eigenvalue:", eigtarget
    print "Guess'*Vtarget:", dot(guess.T, vtarget)
    print "Theta 1:", theta1

    theta, u = davidsolver(A, guess, N, N_gmres, TOL)
    neari = abs(eig - theta1).argmin()
    neareig = eig[neari]
    nearvec = vec[:, [neari]]
    print "RESULTING EIGENVALUE:", theta
    print "Nearest eigenvalue to Theta 1:", neareig
    print "Guess*nearest:", dot(nearvec.T, guess)
    print "Computed smallest and largest using eig:"
    print Eig[-1], ", ", Eig[0]


def run_small():

    # Running parameters
    TOL = 1.e-10
    N = 10
    N_gmres = 5
    Path = '/Users/jonathanlarsson/shule/kandidat/singleparticle/H_small.dat'
    H = np.load(Path)
    H = H[:, 0] + 1j * H[:, 1]
    dim = np.sqrt(len(H))
    H = H.reshape(dim, dim)

    double_reso = 160 * 128 + 120  # 0-index: 129 resp 120
    guess = np.zeros([dim, 1])
    guess[double_reso] = 1 + 0.0001j  # v'*A*v - theta fix

    theta, u = davidsolver(A, guess, N, N_gmres, TOL)


if __name__ == "__main__":
    davidsontest()
