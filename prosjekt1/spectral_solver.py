import numpy as np
import matplotlib.pyplot as plt

def lagrange (x_values, indeks, x):
    teller = 1
    nevner = 1
    for i in range (len(x_values)):
        if i != indeks:
            teller = teller*(x - x_values[i])
            nevner = nevner*(x_values[indeks] - x_values[i])
    return teller / nevner

def dd_lagrange (x_values, indeks, x):
    teller = 0
    nevner = 1
    for i in range(len(x_values)):
        if i != indeks:
            nevner *= (x_values[indeks] - x_values[i])
    for i in range(len(x_values)):
        if i != indeks:
            for j in range(len(x_values)):
                if j != indeks and j != i:
                    ledd = 1
                    for l in range(len(x_values)):
                        if l != indeks and l != i and l != j:
                            ledd *= (x - x_values[l])
                    teller += ledd
    return teller / nevner
    
def spectral_laplace(x_values, dd_math_function, sigma, ua, ub):
    """Takes in interpolation points x_points, double derivative and start end points to
    numerically solve the Poisson-problem with the two Dirichlet-values
    """
    B = []
    for x in x_values:
        B += [-dd_math_function(x, sigma)]
    B[0] = ua
    B[len(x_values) - 1] = ub
    #B ferdig
    A=[]
    for i in range (len(x_values)):
        a = []
        for j in range (len(x_values)):
            if i == 0 or i == len(x_values) - 1:
                a.append(lagrange(x_values, j, x_values[i]))
            else:
                a.append(dd_lagrange(x_values, j, x_values[i]))
        A.append(a)
    #A ferdig
    return np.linalg.solve(A, B)
    
def main():

    ua = -1
    ub = 1
    x_start = 0
    x_end = 10
    u_avg = 3.0
    N = 35
    upper = 2000.
    lower = 0.
    
    def f(x, sigma):
        mu = 2.0
        return np.exp(-(x-mu)**2/(sigma**2))
    
    x_values = (x_end + x_start)/2 - (x_end - x_start)/2*np.cos(np.arange(N)*np.pi/(N-1))
    print(x_values)
    
    diff_u = 1
        
    while np.abs(diff_u) > 0.000001:
        sigma = (lower + upper)/2
        print(sigma)
        u_values = spectral_laplace(x_values, f, sigma, ua, ub)
        
        approx_integral = 0
        for j in range(len(u_values)-1):
            approx_integral += (u_values[j] + u_values[j+1])/2*(x_values[j+1]-x_values[j])
        
        diff_u = u_avg - approx_integral/(x_end - x_start)
        print(diff_u)
        if diff_u > 0:
            lower = sigma
        else:
            upper = sigma
            
        plt.axis([0, 10, -1, 5])
        plt.plot(x_values, u_values)
        plt.show(block=True)

if __name__ == "__main__":
    main()
