# -*- coding: utf-8 -*-
"""
Created on Sat Mar  4 17:06:33 2017

@author: eriko
"""

import numpy as np
import math

# values:

N = 4
systemSize = 4 * N


# make matrices

def B(i):
    B = i
    return B


k_matrisen = 4 * [[]]
k_matrisen[0] = [6, 0, 0, 1]
k_matrisen[1] = [0, 2, 0, 0]
k_matrisen[2] = [0, 0, 1, 0]
k_matrisen[3] = [0, 0, 0, 1]

k_1_matrisen = 4 * [[]]
k_1_matrisen[0] = [-6, 0, 0, 0]
k_1_matrisen[1] = [-6, -2, 0, 0]
k_1_matrisen[2] = [-3, -2, -1, 0]
k_1_matrisen[3] = [-1, -1, -1, -1]


def makeAMatrix(systemSize):
    A = np.zeros((systemSize, systemSize))
    for k in range(0, systemSize, 4):
        for j in range(0, systemSize, 4):
            if k == j:  # legger inn k-1 matrisen
                for m in range(4):
                    for n in range(4):
                        A[k + m][j + n] = k_1_matrisen[m][n]
                print(A)
            elif k == (j + 4) and j + 4 <= systemSize:  # legger inn k-matrisen
                for m in range(4):
                    for n in range(4):
                        A[k + m][j + n] = k_matrisen[m][n]



                        # Om vi går ut av matrisen vil vi legge k i høyre hjørne
            for m in range(4):
                for n in range(4):
                    A[m][n + systemSize - 4] = k_matrisen[m][n]
    return A


u = [1, 0.5, 1 / 6, 1 / 24, 1, 0.5, 1 / 6, 1 / 24, 1, 0.5, 1 / 6, 1 / 24, 1, 0.5, 1 / 6, 1 / 24]

# def main():
A = makeAMatrix(systemSize)

print('\n Matrisen A ser slik ut: \n')
for x in range(16):
    for y in range(16):
        print(int(A[x][y]), '\t', end='')
    print('\n')

x = np.linalg.solve(A, u)

print('Løsningen er x= ', x[0:4])

print('a= ', "%.5f" % x[0])
print('b= ', "%.5f" % x[1])
print('c= ', "%.5f" % x[2])
print('d= ', "%.5f" % x[3])
