import math
import numpy as np
import scipy
import matplotlib.pyplot as plt

def spectral_laplace_lhs(x_points):
    """
    Takes in an array of x grid points and generates a fitting diffmatrix
    """
    x_length = len(x_points)
    matrix_a = np.zeros((x_length, x_length))
    for i in range(x_length-1):
        matrix_a[i, i] = -2
        matrix_a[i+1, i] = 1
        matrix_a[i, i+1] = 1
    matrix_a[x_length-1, x_length-1] = -2
    print np.matrix(matrix_a)
    return matrix_a

def spectral_laplace_rhs(x_points, math_function, start, end):
    """
    Returns the right hand side of the eqation AU = B which needs
    to be solved in a Dirichlet boundary problem.
    Params:
     x - array of x grid points
     f - mathematical function
     ua - initial boundary value
     ub - concluding boundary value
   """

    x_length = len(x_points)
    start_x = x_points[0]
    end_x = x_points[x_length-1]
    diff_x = (end_x-start_x)/x_length
    matrix_b = np.zeros(x_length)
    for i in range(1, x_length-2):
        matrix_b[i] = math_function(x_points[i])*diff_x**2
    matrix_b[0] = math_function(x_points[0])*diff_x**2-start
    matrix_b[-1] = math_function(x_points[-1])*diff_x**2-end
    print matrix_b
    return matrix_b

def spectral_laplace(x_points, math_function, start_value, end_value):
    """
    Fetches both LHS and RHS of the equation UA = B and solves for matrix U
    """
    matrix_a = spectral_laplace_lhs(x_points)
    matrix_b = spectral_laplace_rhs(x_points, math_function, start_value, end_value)
    return np.linalg.solve(matrix_a, matrix_b)

def calculate_laplace_polynomial(x_uniform, math_f, start_value, end_value):
    """
    NOT YET FINISHED
    """
    u_matrix = spectral_laplace(x_uniform, math_f, start_value, end_value)
    


def main():
    """
    main
    """

    def math_f(x_value):
        """
        math_f
        """
        return math.exp(x_value)*math.cos(8*math.pi*x_value)

    start_value = 0
    end_value = 1
    x_uniform = np.linspace(start_value, end_value, 10)
    print u_matrix
    plt.plot(x_uniform, u_matrix)
    plt.show()

if __name__ == "__main__":
    main()
