

import numpy as np
from scipy.linalg import solve
from pseudospectral import diffmat

def bvp(x,p,q,f,a,b):
    n = len(x)
    D = diffmat(x)
    P = np.diag(p)
    Q = np.diag(q)
    y = np.zeros(n)
    A = -np.dot(np.dot(D,P),D)+Q
    A[-1,:] = D[-1,:] # Modify for Neumann condition
    g = f
    y[0] = a
    y[-1] = b
    y[1:] = solve(A[1:,1:],g[1:]-a*A[1:,0])
    return y

if __name__ == '__main__':
    from matplotlib.pyplot import plot, savefig
    import os
    x,u = modal_poisson1d(50,lambda x: np.exp(x))
    plot(x,u)    
    filename = "poisson1d.png"
    print "Generating output file " + filename

    savefig(filename)