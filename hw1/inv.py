# This program requires the numpy package
# It calls rref.py as a subroutine
# Execute: python inv.py


import numpy as np
from rref import rref
from pprint import pprint as pp

def inv(mat):
    #augment an identity matrix to the input matrix
    #The rref the augmented matrix
    #returns the modified identity matrix
    augmat = np.hstack([mat, np.eye(len(mat))])
    rref(augmat)
    return augmat[:,len(mat):]

def main():
    # Generate and invert a random matrix
    A = np.random.rand(4,4)
    print "Original Matrix"
    pp(A)

    imat = inv(A)

    print
    print "Inverted Matrix"
    pp(imat)

    print

    #Check that A*imatA == Identity Matrix

    isInv = np.isclose(
            (np.matrix(A) * np.matrix(imat)),
            np.eye(len(A)))
    print "Is Inverse: %s" % isInv.all()


if __name__ == "__main__":
    main()
