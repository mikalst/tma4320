"""Basically just a list of model constants kept neatly in its own file
"""

import numpy as np

SYSTEM_SIZE = 16
SIMULATIONS = 50
# ALPHA = 1
LENGTH_BETWEEN = 5E-9
LENGTH_SPRING = 5E-10
K_CONST = 1E-1
RADIUS = 1E-9
THRESHOLD = 6E-10
LAMBDA = 2*1.60E-19/3.4E-10
E_CONST = 3E8
I = np.pi * RADIUS**4 / 4
BETA = K_CONST * LENGTH_BETWEEN **3 / (E_CONST * I)


def main():
    print(BETA)
    return

if __name__ == "__main__":
    main()