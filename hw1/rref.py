#This program needs the numpy python package
import numpy as np
from pprint import pprint as pp

def rref(mat):
    #reduced echelon form
    ref(mat)
    #row reduced echelon form
    backsub(mat)

def ref(mat):
    #row echelon form
    n = len(mat)
    #for each row
    for i in range(n):
        submat = mat[i:]
        #find non zero column
        cols = (submat != 0).any(0).nonzero()[0];
        #return if no nonzero columns
        if len(cols) == 0:
            return 
        col = cols[0]
        #find first nonzero element in column
        row = i + (submat[:,col] != 0).nonzero()[0][0]
        #pivot = (row,col)
        #exchange first row with pivot row
        temp = mat[i]
        mat[i] = mat[row]
        mat[row] = temp
        #divide first row by pivot
        mat[i,:] = mat[i,:]/mat[i][col]
        #subtract scaled first row from all subsequent row
        for j in range(i+1,n):
            mat[j,:] = mat[j,:] - mat[i,:]*(mat[j][col])

def backsub(mat):
    #take a row echelon form and remove non zeros on pivot columns
    n = len(mat)
    for i in range(n-1,0,-1):
        #skip row if all zeros
        if (mat[i] == 0).all():
            continue
        col = mat[i].nonzero()[0][0]
        for j in range(i-1, -1, -1):
            mat[j] = mat[j] - mat[i]*mat[j][col]


def main():
    A = np.matrix( """[1 1 0 0 0 0 -1 -8;
                       0 0 1 1 0 0 -1 -5;
                       0 0 0 0 1 1 -1 -1;
                       0 0 1 0 1 0 -1 -8;
                       1 0 0 0 0 0 -1 -6;
                       0 1 0 1 0 1 -1  0;
                       0 0 0 0 0 1 -1 -13;
                       0 1 0 0 1 0 -1 -5;
                       1 1 1 1 1 1  0 31]""" )
    A = np.array(A)
    print "Initial Matrix"
    pp(A)
    rref(A)
    print
    print "rref'd Matrix"
    pp(A)


if __name__=="__main__":
    main()

