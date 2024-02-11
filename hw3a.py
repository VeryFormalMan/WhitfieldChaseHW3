# Chase Whitfield
# MAE 3403
# 02/12/2024
# Homework 3 Part A

import math

def transpose(matrix):
    """
    Transpose a matrix.
    """
    return [[matrix[j][i] for j in range(len(matrix))] for i in range(len(matrix[0]))]

def determinant(matrix):
    """

    :param matrix: Square matrix which determinant needs to be calculated
    :return:determinant(float): The determinant of matrix input.

    """
    if len(matrix) == 1:
        return matrix[0][0]
    elif len(matrix) == 2:
        return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]
    else:
        det = 0
        for i in range(len(matrix)):
            minor = [row[:i] + row[i+1:] for row in matrix[1:]]
            det += (-1) ** i * matrix[0][i] * determinant(minor)
        return det

def is_positive_definite(matrix):
    """
    Check if matrix is positive
    :param matrix: The square matrix to check for positive definiteness
    :return: bool: True if the matrix is positive definite, false otherwise.
    """
    n = len(matrix)
    for i in range(1, n + 1):
        submatrix = [row[:i] for row in matrix[:i]]
        if determinant(submatrix) <= 0:
            return False
    return True

# Used GPT to help create positive definite matrix

def is_symmetric(matrix):
    """
    Check if matrix is symmetric
    :param matrix: List to check for symmetry.
    :return: bool: True if matrix is symmetric, False otherwise.
    """
    n = len(matrix)
    for i in range(n):
        for j in range(i+1, n):
            if matrix[i][j] != matrix[j][i]:
                return False
    return True

# Rest of the code remains the same...



def multiply_matrices(matrix1, matrix2):
    """
    Multiply two matrices.
    """
    result = []
    for i in range(len(matrix1)):
        row = []
        for j in range(len(matrix2[0])):
            total = 0
            for k in range(len(matrix2)):
                total += matrix1[i][k] * matrix2[k][j]
            row.append(total)
        result.append(row)
    return result

# was having issues with defining multiply_matrices function, used GPT


def cholesky_decomposition(matrix):
    """
    Perform Cholesky decomposition
    :param matrix: Symmetric positive definite matrix to be decomposed
    :return: Lower triangular matrix L such that LL^T = matrix.
    """
    n = len(matrix)
    L = [[0.0] * n for _ in range(n)]

    for i in range(n):
        for j in range(i+1):
            if i == j:
                L[i][j] = (matrix[i][i] - sum(L[i][k]**2 for k in range(j))) ** 0.5
            else:
                L[i][j] = (matrix[i][j] - sum(L[i][k]*L[j][k] for k in range(j))) / L[j][j]

    return L

def solve_with_cholesky(matrix, b):
    """
    Solve equation using Cholesky method.
    :param matrix: The symmetric positive definite matrix
    :param b: The constant vector.
    :return: Solution vector
    """
    L = cholesky_decomposition(matrix)

    # Forward substitution
    n = len(matrix)
    y = [0.0] * n
    for i in range(n):
        y[i] = (b[i] - sum(L[i][j] * y[j] for j in range(i))) / L[i][i]

    # Backward substitution
    x = [0.0] * n
    for i in range(n-1, -1, -1):
        x[i] = (y[i] - sum(L[j][i] * x[j] for j in range(i+1, n))) / L[i][i]

    return x
# Used GPT to help with solve_with_cholesky functions

def lu_decomposition(matrix):
    """
    Perform LU decomposition using Doolittle.
    :param matrix: Matrix to be decomposed.
    :return: Lower triangular matrix L and upper triangle matrix so that LU = matrix.
    """
    n = len(matrix)
    L = [[0.0] * n for _ in range(n)]
    U = [[0.0] * n for _ in range(n)]

    for i in range(n):
        L[i][i] = 1.0

    for i in range(n):
        for j in range(i, n):
            U[i][j] = matrix[i][j] - sum(L[i][k] * U[k][j] for k in range(i))
        for j in range(i+1, n):
            L[j][i] = (matrix[j][i] - sum(L[j][k] * U[k][i] for k in range(i))) / U[i][i]

    return L, U

def solve_with_doolittle(L, U, b):
    """
    Solve using Doolittle.
    :param L: Lower triangular matrix.
    :param U: Upper triangular matrix.
    :param b: The constant vector.
    :return: Solution vector.
    """
    n = len(L)
    y = [0.0] * n
    x = [0.0] * n

    # Forward substitution (Ly = b)
    for i in range(n):
        y[i] = b[i] - sum(L[i][j] * y[j] for j in range(i))

    # Backward substitution (Ux = y)
    for i in range(n-1, -1, -1):
        x[i] = (y[i] - sum(U[i][j] * x[j] for j in range(i+1, n))) / U[i][i]

    return x
# Was having issues with LU decomposition from other file, used GPT to help.

def solve_matrix_equation(matrix, b):
    """
    Solve using Cholesky, if not use Doolittle.
    :param matrix: Coefficient matrix of system of equations.
    :param b: Constant vector.
    :return: Solution vector.
    """
    if is_symmetric(matrix) and is_positive_definite(matrix):
        print("Using Cholesky Method")
        return solve_with_cholesky(matrix, b)
    else:
        print("Using Doolittle Method")
        L, U = lu_decomposition(matrix)
        return solve_with_doolittle(L, U, b)

# Solve Matrices
A1 = [[1, -1, 3, 2], [-1, 5, -5, -2], [3, -5, 19, 3], [2, -2, 3, 21]]
b1 = [15, -35, 94, 1]

A2 = [[4, 2, 4, 0], [2, 2, 3, 2], [4, 3, 6, 3], [0, 2, 3, 9]]
b2 = [20, 36, 60, 122]

solution1 = solve_matrix_equation(A1, b1)
solution2 = solve_matrix_equation(A2, b2)

print("Solution for Problem 1:", solution1)
print("Solution for Problem 2:", solution2)


# Used ChatGPT to help with docstring applications.
