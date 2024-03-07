import numpy as np

x = np.zeros(len(A))
y = np.zeros(len(A))
A = np.array([[1, 1, 0, 3],
              [2, 1, -1, 1],
              [3, -1, -1, 2],
              [-1, 2, 3, -1]
              ])
print(A.shape)
L = np.zeros(A.shape)
for i in range(len(L)):
    for j in range(len(L)):
        if i == j:
            L[i, j] = 1
print(L)
b = np.array([4, 1, -3, 4])
n = len(A)
for k in range(0, n - 1):
    for i in range(k + 1, n):
        factor = A[i, k] / A[k, k]
        for j in range(k, n):
            A[i, j] = A[i, j] - A[k, j] * factor
        b[i] = b[i] - b[k] * factor
        L[i, k] = factor
print(A)
print(b)
print(L)

# _____ forward substitution_____________

y[0] =

def GaussElimination_LU(A_matrix, b_matrix):
    return


# import numpy as np
# from scipy.linalg import lu, solve
#
# # Given matrix A and vector b
# A = np.array([
#     [1, 1, 0, 3],
#     [2, 1, -1, 1],
#     [3, -1, -1, 2],
#     [-1, 2, 3, -1]
# ])
#
# b = np.array([4, 1, -3, 4])
#
# # Perform LU decomposition
# P, L, U = lu(A)
#
# # Solve Ly = b for y
# y = solve(L, b)
#
# # Solve Ux = y for x
# x = solve(U, y)
#
# print("Solution x:", x)

# # ------------------------------------------------------------------
# # LU factorization of a matrix A
# # -------------------------------------------------------------------
# import sys
# import os
# import numpy as np
# import scipy as sp
#
#
# def LU_factor(A):
#     '''The routine gives, A = LU (no pivoting)
#     -- A is the rectangular matrix of size (M,N)
#     -- L is lower traiangular matrix of size (M,M)
#        with unit diagonal elements
#     -- U is upper traiangular matrix of size (M,N)
#     '''
#     M, N = A.shape
#     L = np.identity(M)
#     U = np.zeros((M, N))
#     if M != N:
#         print("Matrix should be square")
#         return None
#
#     # The routine will replace upper triangular part of A with U
#
#     # Step-1
#     U[0, 0] = A[0, 0]
#     if np.isclose(U[0, 0], 0.0, atol=1.0E-12):
#         print("Singlar matrix: factorization not possible")
#         return None
#
#     # Step-2: 1st row U & 1st column of L
#     for j in range(1, N):
#         U[0, j] = A[0, j]
#     for i in range(1, M):
#         L[i, 0] = A[i, 0] / U[0, 0]
#
#     # Step-3: Other rows & columns
#     for i in range(1, M - 1):
#
#         # -------------------------------------
#         sum = 0.0
#         for k in range(0, i):
#             sum += L[i, k] * U[k, i]
#         U[i, i] = A[i, i] - sum
#         if np.isclose(U[i, i], 0.0, atol=1.0E-12):
#             print("Singlar matrix: factorization not possible")
#             return None
#
#         # j-th row of U
#         for j in range(i + 1, N):
#             sum = 0.0
#             for k in range(0, i):
#                 sum += L[i, k] * U[k, j]
#             U[i, j] = A[i, j] - sum
#
#         # j-th col of L
#         for j in range(i + 1, M):
#             sum = 0.0
#             for k in range(0, i):
#                 sum += L[j, k] * U[k, i]
#             L[j, i] = (A[j, i] - sum) / U[i, i]
#
#     # Step-4:
#     sum = 0.0
#     for k in range(0, M - 1):
#         sum += L[M - 1, k] * U[k, M - 1]
#     U[M - 1, M - 1] = A[M - 1, M - 1] - sum
#
#     # LU
#     return L, U
#
#
# A = np.array([[1, 1, 0, 3], [2, 1, -1, 1], [3, -1, -1, 2], [-1, 2, 3, -1]])
#
# #P, L, U = sp.linalg.lu(A)
# # sols = GaussElim_Pivot(A)
# # print("PLU = \n", P @ L @ U)
# #print("L = \n", L)
# #print("U = \n", U)
#
# L, U = LU_factor(A)
# print("\nL = \n", L)
# print("\nU = \n", U)
# print("\nLU = \n", L@U)
