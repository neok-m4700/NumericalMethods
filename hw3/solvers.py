import numpy as np
import scipy as sp
import scipy.sparse as sparse
import scipy.linalg

class Iterative(object):
    def __init__(self, A):
        self.iter = 1000
        self.epsilon = 1e-6
        self.Pinv = self.get_Pinv(A)
        self.t1 = sparse.identity(A.shape[0]) - self.Pinv * A

    def get_Pinv(self, A):
        pass

    def solve(self, b):
        xnew = xold = b
        t2 = self.Pinv * b
        for i in range(self.iter):
            xnew = self.t1 * xold + t2
            xdiff = (xold - xnew)
            err = abs(xdiff.mean())
            if err < self.epsilon:
                break;
            xold = xnew
        return xnew;


class Jacobi(Iterative):

    def get_Pinv(self, A):
        return sparse.identity(A.shape[0])/A[0, 0]


class Gauss(Iterative):

    def get_Pinv(self, A):
        P = sparse.tril(A)
        return sparse.linalg.inv(P.tocsc())


class SOR(Iterative):

    def get_Pinv(self, A):
        w = 1.9
        D = sparse.identity(A.shape[0]) * A[0, 0]
        L = sparse.tril(A)
        P = D = w * L
        return sparse.linalg.inv(P.tocsc())
