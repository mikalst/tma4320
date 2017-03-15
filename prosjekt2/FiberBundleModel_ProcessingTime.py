from FiberBundleModel_Ripping import makeSparseAMatrix
from FiberBundleModel_Ripping import makeRHS
from FiberBundleModel_Ripping import makeAMatrix
import FiberBundleModel_Constants as FBM_c
import numpy as np
from timeit import Timer
import matplotlib.pyplot as plt
from scipy import sparse
from scipy.sparse.linalg import spsolve
from operator import itemgetter
from decimal import Decimal

BETA = 1E-1
L_SYSTEM = 5e-9


def secession(N_SPRINGS):
    """The ripping function
    """
    N = 4*N_SPRINGS
    beta = 0.1
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
        A = makeAMatrix(N, beta_values)
        # A = makeSparseAMatrix(N, beta_values)
        u = makeRHS(N, alpha_one)

        # Solve the system to get the coefficients for the polynomial
        x = np.linalg.solve(A, u)
        # x = sparse.linalg.spsolve(A.tocsc(), u)

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

    # print(alpha_values_simple)

    return alpha_values_simple


def secession_sparse(N_SPRINGS):
    """The ripping function
    """
    N = 4*N_SPRINGS
    beta = 0.1
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
        # A = makeAMatrix(N, beta_values)
        A = makeSparseAMatrix(N, beta_values)
        u = makeRHS(N, alpha_one)

        # Solve the system to get the coefficients for the polynomial
        # x = np.linalg.solve(A, u)
        x = sparse.linalg.spsolve(A.tocsc(), u)

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

    # print(alpha_values_simple)

    return alpha_values_simple


def main():
    # Assert that each spring has 4 corresponding equations
    normal_times = []
    sparse_times = []
    N = [16 + 4*x for x in range(0, 47)]
    N_SPRINGS = [int(n/4) for n in N]
    print(N)
    REPETITIONS = 5
    for n in N:
        print(n)
        t_1 = Timer(lambda: secession(n))
        t_2 = Timer(lambda: secession_sparse(n))
        normal_times.append(t_1.timeit(number=REPETITIONS))
        sparse_times.append(t_2.timeit(number=REPETITIONS))

# Initialize empty figures, with handles, names, sizes, background color and edge color
    fig_solution = plt.figure("Solution", figsize=FBM_c.PLT_SIZE_WIDE, facecolor=FBM_c.PLT_BG_COLOR,
                              edgecolor=FBM_c.PLT_E_COLOR)
    fig_solution.suptitle("Kjøretid løsrivelse, Dense vs Sparse")

    # Initialize the individual subfigures
    ax_y = fig_solution.add_subplot(111)  # Subfigure for the solution, preparing a 4 row

    # ax_y.set_title(r"$N = 16, \alpha = 1, \beta = 1$")  # α = 1
    ax_y.set_title(r"Kjøretid", fontsize=FBM_c.PLT_FNT_SIZE)
    # 1 column canvas, placing this axis in the first location/row
    # Set the x-label only on the bottom figure
    ax_y.set_xlabel(r"$Kroker$" + " " + r"$N$", fontsize=FBM_c.PLT_FNT_SIZE)

    # Set all the y-axis labels for the function and its derivatives
    ax_y.set_ylabel(r"$tid$" + " " + r"$[s]$", fontsize=FBM_c.PLT_FNT_SIZE)

    # Plot beta_x for x in [0, 1, 2] and the corresponding alpha sequence
    ax_y.plot(N_SPRINGS, normal_times)
    ax_y.plot(N_SPRINGS, sparse_times)
    legend_one = [r"$Dense$", r"$Sparse$"]
    ax_y.legend(legend_one, fontsize=FBM_c.PLT_FNT_SIZE)

    # Show figure
    plt.grid(True)
    plt.show()


if __name__ == "__main__":
    main()