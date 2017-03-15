"""Model of Springs and Cable
"""
# -*- coding: utf-8 -*-

import FiberBundleModel_Constants as FBM_c
import FiberBundleModel_ModelParameters as FBM_MP
import numpy as np
import matplotlib.pyplot as plt
from scipy import sparse
from scipy.sparse.linalg import spsolve
from operator import itemgetter
from decimal import Decimal

THRESHOLD = FBM_MP.THRESHOLD
BETAS = [1E-8, 1E-1, 1E8]
N = FBM_MP.SYSTEM_SIZE
N_SPRINGS = int(N / 4)
L_SYSTEM = 5e-9
NUMBER_OF_SIMULATIONS = FBM_MP.SIMULATIONS


def main():
    # Assert that each spring has 4 corresponding equations
    assert(N//4)

    secession_several()


def plot_sequence():
    return;


def secession_several():
    """Calculates the alpha values of several beta
    :return:
    """

    beta_and_corresponding_alpha_values = []

    # Iterate over the 3 beta values, store the alpha sequence of each
    for beta_value in BETAS:

        # Add initial alpha sequence
        alpha_sequence = np.array(secession(beta_value))

        # Calculate additional alpha sequences to calculate the average
        for simulation in range(1, NUMBER_OF_SIMULATIONS):

            # Print progress
            print(simulation/NUMBER_OF_SIMULATIONS * 100)
            alpha_sequence += np.array(secession(beta_value))

        # Calculate the average alpha sequence for the current beta
        alpha_sequence_average = [element / NUMBER_OF_SIMULATIONS for element in alpha_sequence]

        # Append beta and its corresponding average alpha sequence
        beta_and_corresponding_alpha_values.append([beta_value, alpha_sequence_average])

    # Plot alpha sequences for beta 0, 1 and 2
    plotAlpha([item[0] for item in beta_and_corresponding_alpha_values],
              [item[1] for item in beta_and_corresponding_alpha_values],
              N_SPRINGS)


# noinspection PyPep8Naming
def secession(beta):
    """The ripping function
    """
    # Convert from x to eta
    thresholdMax = 0.1e-9 / L_SYSTEM

    # Array for storing beta values
    beta_values = [beta for each in range(N_SPRINGS)]

    # Check for when to stop iteration
    beta_values_final = [0 for each in range(N_SPRINGS)]

    # Create random threshold values between 0 and thresholdMax at which point the springs break
    threshold_values = np.random.uniform(0, thresholdMax, N)
    only_one_string_left = False
    alpha_values = []
    alpha_values_simple = []
    alpha_one = 1

    while not only_one_string_left:

        # print("Beta values: ", beta_values)

        A = makeAMatrix(N, beta_values)
        # A = makeSparseAMatrix(N, beta_values)
        u = makeRHS(N, alpha_one)

        # Solve the system to get the coefficients for the polynomial
        x = np.linalg.solve(A, u)
        # x = sparse.linalg.spsolve(A.tocsc(), u)
        # print(x)

        # Find the index of spring that is going to rip first
        r_k_of_all_springs = [-beta_values[k]*x[4*k+3]/threshold_values[k] for k in range(N_SPRINGS)]
        r_spring_max_value = max(r_k_of_all_springs)
        r_spring_max_index = max(enumerate(r_k_of_all_springs), key=itemgetter(1))[0]

        # Calculate the resulting alpha and store it
        r_spring_max_alpha = beta_values[r_spring_max_index] / r_spring_max_value

        # print(r_spring_max_alpha)
        alpha_values.append([r_spring_max_index, r_spring_max_alpha])
        alpha_values_simple.append(r_spring_max_alpha)

        # Set necessary decrements to prepare for next loop
        beta_values[r_spring_max_index] = 0

        if beta_values == beta_values_final:
            only_one_string_left = True

    # print(alpha_values)

    return alpha_values_simple


def makeAMatrix(systemSize, beta_values):
    """Creates matrix A
    """
    A = np.zeros((systemSize, systemSize))
    n_springs = int(systemSize / 4)

    k_matrices_by_spring = []

    zero_four_by_four = [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]]

    total_matrix = []

    for spring in range(n_springs):

        k_matrix = 4 * [[]]
        k_matrix[0] = [6, 0, 0, beta_values[spring]]
        k_matrix[1] = [0, 2, 0, 0]
        k_matrix[2] = [0, 0, 1, 0]
        k_matrix[3] = [0, 0, 0, 1]

        k_1_matrix = 4 * [[]]
        k_1_matrix[0] = [-6, 0, 0, 0]
        k_1_matrix[1] = [-6, -2, 0, 0]
        k_1_matrix[2] = [-3, -2, -1, 0]
        k_1_matrix[3] = [-1, -1, -1, -1]

        k_matrices_by_spring.append([k_matrix, k_1_matrix])

    total_matrix.append([k_matrices_by_spring[0][0]] + (n_springs-2)*[zero_four_by_four] + [k_matrices_by_spring[0][1]])

    for spring in range(1, n_springs):
        total_matrix.append(
            (spring-1) * [zero_four_by_four] + [k_matrices_by_spring[spring][1]] + [k_matrices_by_spring[spring][0]]
            + (n_springs - spring - 1) * [zero_four_by_four])

    A = np.bmat(total_matrix)

    return A


def makeRHS(systemSize, alpha):
    """Makes right hand side of Ax = u
    """

    u = np.zeros(systemSize)

    for i in range(int(systemSize/4)):
        u[4*i:4*i+4] = [-1*alpha, -0.5*alpha, -1/6*alpha, -1/24*alpha]

    return u


def makeSparseAMatrix(systemSize, beta_values):
    """Creates sparse matrix A
    """
    A_sparseMatrix = sparse.dok_matrix((systemSize, systemSize))

    for row in range(0, systemSize, 4):
        A_sparseMatrix[row, row] = -6
        A_sparseMatrix[row+1, row+1] = -2
        A_sparseMatrix[row+2, row+2] = -1
        A_sparseMatrix[row+3, row+3] = -1

        A_sparseMatrix[row+1, row] = -6

        A_sparseMatrix[row+2, row] = -3
        A_sparseMatrix[row+2, row+1] = -2

        A_sparseMatrix[row+3, row] = -1
        A_sparseMatrix[row+3, row+1] = -1
        A_sparseMatrix[row+3, row+2] = -1

    for row in range(4, systemSize, 4):
        A_sparseMatrix[row, row-4] = 6
        A_sparseMatrix[row+1, row-3] = 2
        A_sparseMatrix[row+2, row-2] = 1
        A_sparseMatrix[row+3, row-1] = 1

        A_sparseMatrix[row, row-1] = beta_values[int(row/4)]

    # Element in top right corner
    A_sparseMatrix[0, systemSize-4] = 6
    A_sparseMatrix[1, systemSize-3] = 2
    A_sparseMatrix[2, systemSize-2] = 1
    A_sparseMatrix[3, systemSize-1] = 1

    A_sparseMatrix[0, systemSize-1] = beta_values[0]

    return A_sparseMatrix


def makeSparseRHS(systemSize, alpha):
    """Makes right hand side of Ax = u
    """

    u = np.zeros(systemSize)

    for i in range(int(systemSize/4)):
        u[4*i:4*i+4] = [-1*alpha, -0.5*alpha, -1/6*alpha, -1/24*alpha]

    return u


def plotAlpha(beta_values, alpha_lists, n_total_springs):
    """Plot alpha curves of several beta

    :param beta_values: list of ints
    :param alpha_lists: list of lists of ints
    :param n_total_springs: int
    :return: void
    """
    # Initialize empty figures, with handles, names, sizes, background color and edge color
    fig_solution = plt.figure("Solution",
                              figsize=FBM_c.PLT_SIZE_WIDE,
                              facecolor=FBM_c.PLT_BG_COLOR,
                              edgecolor=FBM_c.PLT_E_COLOR)
    fig_solution.suptitle("Løsrivelse av kabelligningen")

    # Initialize the individual subfigures
    ax_y = fig_solution.add_subplot(311)  # Subfigure for the solution, preparing a 4 row

    # ax_y.set_title(r"$N = 16, \alpha = 1, \beta = 1$")  # α = 1
    ax_y.set_title(r"$Kroker$" + " " + r"$N = %s $" % N_SPRINGS + ",   "
                   + r"$Simuleringer$" + " " + r"$M = %s $" % NUMBER_OF_SIMULATIONS)
    # 1 column canvas, placing this axis in the first location/row
    ax_dy = fig_solution.add_subplot(312, sharex=ax_y)  # Subfigure for the first derivative

    #  with exactly the same x-axis as ax_y, in the second location/row
    ax_ddy = fig_solution.add_subplot(313, sharex=ax_y)  # Subfig for the second derivative

    # Set the x-label only on the bottom figure
    ax_ddy.set_xlabel(r"$\frac{n}{N}$", fontsize=FBM_c.PLT_FNT_SIZE)

    # Set all the y-axis labels for the function and its derivatives
    ax_y.set_ylabel(r"$\alpha$", fontsize=FBM_c.PLT_FNT_SIZE)
    ax_dy.set_ylabel(r"$\alpha$", fontsize=FBM_c.PLT_FNT_SIZE)
    ax_ddy.set_ylabel(r"$\alpha$", fontsize=FBM_c.PLT_FNT_SIZE)

    # Initialize x values
    x_values = [n/n_total_springs for n in range(len(alpha_lists[0]))]

    # Plot beta_x for x in [0, 1, 2] and the corresponding alpha sequence
    ax_y.plot(x_values, alpha_lists[0])
    legend_one = r"$\beta = $" + '%.2E' % Decimal(beta_values[0])
    ax_y.legend([legend_one])
    ax_dy.plot(x_values, alpha_lists[1])
    legend_one = r"$\beta = $" + '%.2E' % Decimal(beta_values[1])
    ax_dy.legend([legend_one])
    ax_ddy.plot(x_values, alpha_lists[2])
    legend_one = r"$\beta = $" + '%.2E' % Decimal(beta_values[2])
    ax_ddy.legend([legend_one])

    # Show figure
    plt.show()

if __name__ == "__main__":
    main()
