import numpy as np
from rref import rref
from pprint import pprint as pp

def inv(mat):
    augmat = np.hstack([mat, np.eye(len(mat))])
    rref(augmat)
    return augmat[:,len(mat):]

def main():
    A = np.random.rand(4,4)
    print "Original Matrix"
    pp(A)

    imat = inv(A)

    print
    print "Inverted Matrix"
    pp(A)

    print
    isInv = np.isclose(
            (np.matrix(A) * np.matrix(imat)),
            np.eye(len(A)))
    print "Is Inverse: %s" % isInv.all()


if __name__ == "__main__":
    main()
