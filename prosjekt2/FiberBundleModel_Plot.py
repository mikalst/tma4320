# -*- coding: utf-8 -*-

import FiberBundleModel_Constants as FBM_c
import FiberBundleModel_ModelParameters as FBM_MP
import numpy as np
import matplotlib.pyplot as plt
from decimal import Decimal

# 'Import' proper ALPHA and BETA from model constants file
BETA = FBM_MP.BETA
ALPHA = FBM_MP.LAMBDA * 2.125E6 * FBM_MP.LENGTH_BETWEEN**3 / (0.3E9 * np.pi * FBM_MP.RADIUS**4 / 4)
ALPHA = 1E4
BETA = 1E4


def main():
    assert(FBM_MP.SYSTEM_SIZE // 4)
    plot_numerical_analytical()


# noinspection PyPep8Naming
def plot_numerical_analytical():
    """The main function
    """
    N = FBM_MP.SYSTEM_SIZE  # Number of equations in system

    A = makeAMatrix(N, BETA)
    u = makeRightHandSide(N)

    # Solve the system to get the coefficients for the polynomial
    x = np.linalg.solve(A, u)
    x_analytical = [ALPHA / 12, - ALPHA / 24, 0, - ALPHA / BETA]
    diff_x = x_analytical - x[0:4]

    print('Solution x =', x[0:4])
    print("Difference A - N = ", diff_x[0:4])

    print('a= ', "%.5f" % x[0])
    print('b= ', "%.5f" % x[1])
    print('c= ', "%.5f" % x[2])
    print('d= ', "%.5f" % x[3])

    # Plot the solution
    print("Using alpha = %.5f" % ALPHA)
    plotSolution(x, ALPHA)


def makeAMatrix(systemSize, beta):
    """Creates matrix A
    """
    A = np.zeros((systemSize, systemSize))

    # Creates the k matrix and k-1 matrix
    k_matrix = 4 * [[]]
    k_matrix[0] = [6, 0, 0, beta]
    k_matrix[1] = [0, 2, 0, 0]
    k_matrix[2] = [0, 0, 1, 0]
    k_matrix[3] = [0, 0, 0, 1]

    k_1_matrix = 4 * [[]]
    k_1_matrix[0] = [-6, 0, 0, 0]
    k_1_matrix[1] = [-6, -2, 0, 0]
    k_1_matrix[2] = [-3, -2, -1, 0]
    k_1_matrix[3] = [-1, -1, -1, -1]

    for k in range(0, systemSize, 4):

        for j in range(0, systemSize, 4):
            if k == j:  # legger inn k-1 matrisen
                for m in range(4):
                    for n in range(4):
                        A[k + m][j + n] = k_1_matrix[m][n]
            elif k == (j + 4) and j + 4 <= systemSize:  # legger inn k-matrisen
                for m in range(4):
                    for n in range(4):
                        A[k + m][j + n] = k_matrix[m][n]

                        # Om vi går ut av matrisen vil vi legge k i høyre hjørne
            for m in range(4):
                for n in range(4):
                    A[m][n + systemSize - 4] = k_matrix[m][n]

    print(A)

    return A


def makeRightHandSide(systemSize):
    """Makes right hand side of Ax = u
    """

    u = np.zeros(systemSize)

    for i in range(int(systemSize/4)):
        u[4*i:4*i+4] = [-1 * ALPHA, -0.5 * ALPHA, -1 / 6 * ALPHA, -1 / 24 * ALPHA]

    return u


def plotSolution(solutionVector, alpha):
    """
    :param solutionVector: x vector contain a_i, b_i, c_i and d_i
    :param alpha: alpha value to used in computation
    :return:
    """
    # The plotting is fine, now I should focus on creating working matrices for arbitrary N

    # Initialize empty figures, with handles, names, sizes, background color and edge color
    fig_solution = plt.figure("Solution", figsize=FBM_c.PLT_SIZE_WIDE, facecolor=FBM_c.PLT_BG_COLOR,
                              edgecolor=FBM_c.PLT_E_COLOR)
    fig_solution.suptitle("Løsning av kabelligningen")
    # Initialize the individual subfigures
    ax_y = fig_solution.add_subplot(411)  # Subfigure for the solution, preparing a 4 row
    # ax_y.set_title(r"$N = 16, \alpha = 1, \beta = 1$")  # α = 1
    ax_y.set_title(r"$Kroker$" + " " + r"$N = 4$" + ", " + r"$\alpha = %.2E $" % Decimal(ALPHA) + ", " r"$\beta = %.2E $" % Decimal(BETA))  # α = λE0l**3/B
    # 1 column canvas, placing this axis in the first location/row
    ax_dy = fig_solution.add_subplot(412, sharex=ax_y)  # Subfigure for the first derivative,
    #  with exactly the same x-axis as ax_y, in the second location/row
    ax_ddy = fig_solution.add_subplot(413, sharex=ax_y)  # Subfig for the second derivative,
    # with exactly the same x-axis as ax_y
    ax_dddy = fig_solution.add_subplot(414, sharex=ax_y)  # Subfig for the third derivative,
    # with exactly the same x-axis as ax_y
    # ax_y.set_xlim([0, 4])  # Set the x-axis limits
    # (and thus those of all the other subplots in the figure)

    # Set the x-label only on the bottom figure
    ax_dddy.set_xlabel("$\\xi$", fontsize=FBM_c.PLT_FNT_SIZE)

    # Set all the y-axis labels for the function and its derivatives
    ax_y.set_ylabel("$\eta$", fontsize=FBM_c.PLT_FNT_SIZE)
    ax_dy.set_ylabel("$d\eta/d\\xi$", fontsize=FBM_c.PLT_FNT_SIZE)
    ax_ddy.set_ylabel("$dˆ2\eta/d\\xiˆ2$", fontsize=FBM_c.PLT_FNT_SIZE)
    ax_dddy.set_ylabel("$dˆ3\eta/d\\xiˆ3$", fontsize=FBM_c.PLT_FNT_SIZE)

    # Solution domain between springs, just as in eq. (11)
    xi_minus_k = np.linspace(0, 1.0, 200)

    # Solution domain total
    total_xi_minus_k = np.linspace(0, 1.0 * FBM_MP.SYSTEM_SIZE / 4, 200 * FBM_MP.SYSTEM_SIZE / 4)

    # Prepare empty vectors for the values
    total_eta = []
    total_eta_exact = []

    dtotal_eta = []
    dtotal_eta_exact = []

    ddtotal_eta = []
    ddtotal_eta_exact = []

    dddtotal_eta = []
    dddtotal_eta_exact = []

    # Calculate solutions
    for k in range(int(FBM_MP.SYSTEM_SIZE / 4)):

        # Solution polynomial, eq. (11)
        eta = -(alpha / 24) * xi_minus_k ** 4 + solutionVector[4 * k] * xi_minus_k ** 3 \
              + solutionVector[4 * k + 1] * xi_minus_k ** 2 \
              + solutionVector[4 * k + 2] * xi_minus_k \
              + solutionVector[4 * k + 3]

        eta_exact = -(alpha / 24) * xi_minus_k ** 4 + alpha / 12 * xi_minus_k ** 3 \
              - alpha / 24 * xi_minus_k ** 2 \
              - alpha / BETA

        deta = -(alpha / 6) * xi_minus_k ** 3 + 3 * solutionVector[4 * k] * xi_minus_k ** 2 \
               + 2 * solutionVector[4 * k + 1] * xi_minus_k \
               + solutionVector[4 * k + 2]

        deta_exact = -(alpha / 6) * xi_minus_k ** 3 + 3 * alpha / 12 * xi_minus_k ** 2 \
               + 2 * - alpha / 24 * xi_minus_k

        ddeta = -(alpha / 2) * xi_minus_k ** 2 + 6 * solutionVector[4 * k] * xi_minus_k \
                + 2 * solutionVector[4 * k + 1]

        ddeta_exact = -(alpha / 2) * xi_minus_k ** 2 + 6 * alpha / 12 * xi_minus_k \
                + 2 * - alpha / 24

        dddeta = -alpha * xi_minus_k + 6 * solutionVector[4 * k]

        dddeta_exact = -alpha * xi_minus_k + 6 * alpha / 12

        # Append all the values to the value lists
        total_eta.extend(eta)
        total_eta_exact.extend(eta_exact)
        dtotal_eta.extend(deta)
        dtotal_eta_exact.extend(deta_exact)
        ddtotal_eta.extend(ddeta)
        ddtotal_eta_exact.extend(ddeta_exact)
        dddtotal_eta.extend(dddeta)
        dddtotal_eta_exact.extend(dddeta_exact)

        # Plot the cable location in both the spring and the cable plots
        #  (it's clear when noticing the difference between ax_y and ax_eta

        # Plotting the derivatives
        """
        ax_dy.plot(xi_minus_k, deta)
        ax_ddy.plot(xi_minus_k, ddeta)
        ax_dddy.plot(xi_minus_k, dddeta)
        """

    # Plotting eta and its derivatives, both numerical and analytical
    ax_y.plot(total_xi_minus_k, total_eta, lw=2.0)
    ax_y.plot(total_xi_minus_k, total_eta_exact, lw=2.0)
    ax_y.legend(["Numerisk", "Analytisk"], framealpha=0.8, loc=4)
    ax_dy.plot(total_xi_minus_k, dtotal_eta)
    ax_dy.plot(total_xi_minus_k, dtotal_eta_exact)
    # ax_dy.legend(["$d\eta/dx$", "$d\eta/dx_{exact}$"], framealpha=0.5)
    ax_ddy.plot(total_xi_minus_k, ddtotal_eta)
    ax_ddy.plot(total_xi_minus_k, ddtotal_eta_exact)
    # ax_ddy.legend(["$d^2\eta/dx^2$", "$d^2\eta/dx^2_{exact}$"], framealpha=0.5)
    ax_dddy.plot(total_xi_minus_k, dddtotal_eta)
    ax_dddy.plot(total_xi_minus_k, dddtotal_eta_exact)
    # ax_dddy.legend(["$d^3\eta/dx^3$", "$d^3\eta/dx^3_{exact}$"], framealpha=0.5)

    # Make the figure appear on-screen:
    plt.show()

if __name__ == "__main__":
    main()