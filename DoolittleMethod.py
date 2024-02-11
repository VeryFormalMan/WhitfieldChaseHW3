#region imports
from copy import deepcopy as dcpy
from math import cos,pi
import Gauss_Seidel as GS
#endregion

#region Functions
def LUFactorization(A):
    """
    This is the Lower-Upper factorization part of Doolittle's method.
    :param A: a nxn matrix
    :return: a tuple with (L, U)
    """
    n = len(A)
    U = [[0] * n for _ in range(n)]
    L = [[0] * n for _ in range(n)]

    for i in range(n):
        # Check for zero values on the diagonal of U
        if A[i][i] == 0:
            raise ValueError("Matrix is not suitable for LU factorization (zero on diagonal of U)")

        # Calculate the diagonal elements of U
        U[i][i] = A[i][i]

        # Calculate the elements of L and U
        for k in range(i + 1, n):
            L[k][i] = A[k][i] / U[i][i]
            U[i][k] = A[i][k]
        for j in range(i + 1, n):
            for k in range(i + 1, n):
                A[j][k] = A[j][k] - L[j][i] * U[i][k]

    # Set diagonal elements of L to 1
    for i in range(n):
        L[i][i] = 1

    return L, U


def BackSolve(A,b,UT=True):
    """
    This is a backsolving algorithm for a matrix and b vector where A is triangular
    :param A: A triangularized matrix (Upper or Lower)
    :param b: the right hand side of a matrix equation Ax=b
    :param UT: boolean of upper triangular (True) or lower triangular (False)
    :return: the solution vector x, from Ax=b
    """
    nRows=len(b)
    nCols=nRows
    x=[0]*nRows
    if UT:
        for nR in range(nRows-1,-1,-1):
            s=0
            for nC in range(nR+1,nRows):
                s+=A[nR][nC]*x[nC]
            x[nR]=1/A[nR][nR]*(b[nR]-s)
    else:
        for nR in range(nRows):
            s=0
            for nC in range(nR):
                s+=A[nR][nC]*x[nC]
            x[nR]=1/A[nR][nR]*(b[nR]-s)
    B = GS.checkMatrixSoln(A, x, False)
    return x

def Doolittle(Aaug):
    """

    :param Aaug: the augmented matrix
    :return: the solution vector x
    """
    A,b=GS.separateAugmented(Aaug)
    L,U=LUFactorization(A)
    B=GS.matrixMult(L,U)
    y=BackSolve(L,b, UT=False)
    x=BackSolve(U,y, UT=True)
    return x

def main():
    A=[[3, 5, 2],[0,8,2],[6,2,8]]
    L,U=LUFactorization(A)
    print("L:")
    for r in L:
        print(r)

    print("\nU:")
    for r in U:
        print(r)

    aug=[[3,9,6,4.6],[18,48,39,27.2], [9,-27,42,9]]
    aug = [[3, 1, -1, 2],
          [1, 4, 1, 12],
          [2, 1, 2, 10]]
    x=Doolittle(aug)
    x=[round(y,3) for y in x]
    print("x: ", x)
    y=GS.GaussSeidel(aug,[0,0,0])
    y=[round(z,3) for z in y]
    b=GS.checkMatrixSoln(aug,y)
    print("b: ",b)
#endregion

if __name__ == "__main__":
    main()