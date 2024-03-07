import numpy as np

A = np.array([[3, -2, 5, 0],
              [1, 1, 2, 1],
              [2, 7, 6, 5],
            [4, 5, 8, 1]], float)
b = np.array([2, 5, 7, 4], float)


n = len(b)
x = np.zeros(n)
# ____________Partial pivoting_______________
# for k in range(0, n - 1):
#     max_ind = k
#     for i in range(k+1, n):
#         if A[k, k] < A[i, k]:
#             max_ind = i
#
#     A[[k, max_ind]] = A[[max_ind, k]]
#     b[[k, max_ind]] = b[[max_ind, k]]

# ________________Complete pivoting___________________
# for k in range(0, n-1):
#     max_row_index, max_col_index = np.unravel_index(np.argmax(A[k:, k:]), A[k:, k:].shape)
#     A[[k, max_row_index + k]] = A[[max_row_index + k, k ]]
#     b[[k, max_row_index + k]] = b[[max_row_index + k, k]]
#     A[:, [k, max_col_index + k]] = A[:, [max_col_index + k, k]]
# print(A)

# _________Row echelon_________________
for k in range(0, n - 1):
    for i in range(k+1, n):
        if A[k,k] == 0:
            if A[i, k] != 0:
                A[[k, i]] = A[[i,k]]
                b[[k, i]] = b[[i, k]]
        if A[i, k] == 0:
            pass
        factor = A[i, k] / A[k, k]
        for j in range(k, n):
            A[i,j] = A[i, j] - A[k, j] * factor
        b[i]= b[i] - b[k] * factor

print(A)
print(b)
# ___________Back Substitution________________-
x[n-1] = b[n-1]/ A[n-1, n-1]
for i in range(n-2, -1, -1):
    sum = 0
    for j in range(i+1, n):
        sum = sum + A[i, j] * x[j]
    x[i] = (b[i] - sum) / A[i, i]

print(x)
#_____________Using package______________
print(np.linalg.solve(A, b))