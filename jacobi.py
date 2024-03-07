import numpy.linalg as alg
import numpy as np


# Jacobi fixed point iterative method to solve system of linear equations
def jacobi_iterative(A_matrix, b_matrix, max_iter: int = 100, tol: float = 0.0001):
    x_guess = np.zeros(len(A_matrix))  # random guess initialization
    x = np.zeros(len(A_matrix))  # copy of x_guess to store the next i+1 guess by update
    for iteration in range(max_iter):  # loop through maximum iterations
        for i in range(len(A_matrix)):  # i = 0, 1, 2,...n- dimension of system of linear equation
            dum_sum = 0  # initialize summation
            for j in range(len(A_matrix)):  # j = 0, 1, 2,...n- dimension of system of linear equation
                if j != i:  # taking non-diagonal elements
                    dum_sum += A_matrix[i][j] * x_guess[j]  # summation term
            x[i] = (- dum_sum + b_matrix[i]) / A_matrix[i][i]  # calculating next guess
        if abs(alg.norm(x - x_guess)) < tol:  # checking if relative error is exceeding tolerance value
            print(f"The solution at {iteration} iteration\nx: {x}")
            raise Warning('The procedure is successful')
    # maximum iterations reached and returns solution as x matrix
    return print(f"Maximum iterations reached!\nThe solution at {iteration} iteration\nx: {x}")


# matrix of system of linear equations
A = np.array([[1, 1, 0, 3],
              [2, 1, -1, 1],
              [3, -1, -1, 2],
              [-1, 2, 3, -1]
              ])
# constants
b = np.array([4, 1, -3, 4])

jacobi_iterative(A, b, max_iter=1000, tol=0.001)