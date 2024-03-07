import numpy as np
from scipy.linalg import norm



# matrix of system of linear equations
A = np.array([[10, -1, 2, 0.],
              [-1, 11, -1, 3],
              [2, -1, 10, -1],
              [0., 3, -1, 8]])

# constants
b = np.array([6, 25, -11, 15])


n = len(A)
x_guess = np.zeros(len(A))
x = np.zeros(len(A))
for iteration in range(1000):
    for i in range(n):
        sum_x = 0
        for j in range(0, i -1):
            sum_x = sum_x + A[i, j] * x[j]
        sum_xo = 0
        for j in range(i+1, n):
            sum_xo = sum_xo + A[i, j] * x_guess[j]

        x[i] = (-sum_x - sum_xo + b[i]) / A[i, i]


    if abs(norm(x - x_guess)) < 1.0e-3:
        print(f"The solution at {iteration} iteration\nx: {x}")
        raise Warning('The procedure is successful')

# maximum iterations reached and returns solution as x matrix
    print(f"Maximum iterations reached!\nThe solution at {iteration} iteration\nx: {x}")


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
            #print(f"The solution at {iteration} iteration\nx: {x}")
            raise Warning('The procedure is successful')
    # maximum iterations reached and returns solution as x matrix
    return print(f"Maximum iterations reached!\nThe solution at {iteration} iteration\nx: {x}")


# matrix of system of linear equations
A = np.array([[10, -1, 2, 0.],
              [-1, 11, -1, 3],
              [2, -1, 10, -1],
              [0., 3, -1, 8]])

# constants
b = np.array([6, 25, -11, 15])

jacobi_iterative(A, b, max_iter=1000, tol=0.001)