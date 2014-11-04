"""
Lattice.py
This program has a function stringBitHailtonian(L) which takes in the
side legnth of a square toric lattice and outputs the Hamiltonian of the
lattice in a binary string bit 2n by n matrix representation of the form
    [ 0 A ]
    [ B 0 ]
where A represents the block matrix due to the vertices, B represents
the block matrix due to the placets and 2n = 2L^2 is the number of qubits
on the lattice's edges.

The way each vertex lattice is numbered starts from 0 to n-1 and goes row
by row. The plaquettes are numbered in a similar manner and the qubits are
also numbered row by row, but from 0 to 2n-1.

The tensor product of Pauli and Identity matrices for each vertex and
placet is kind of effectively what the output matrix represents.
It should hold some information about the commutation relations of the
placet and vertex operators. The Hamiltonian takes the form
   $$H = -\sum B_p - \sum A_v.$$
where B_p and A_v are the placet and vertex operators respectively.

The main() function loops through L = 2 to 8 and prints out the matrix
each for each L. It's just here to test the stringBitHamiltonian(L)
function.
"""

import sys
import numpy as np

def main():
    for L in xrange(2,9):
        print('String bit matrix representation of a '+str(L)+' by '+str(L)+' toric lattice')
        show(stringBitHamiltonian(L))

def stringBitHamiltonian(L):
    n = L*L
    A = np.zeros((2*n,n),dtype=bool)
    B = np.zeros((2*n,n),dtype=bool)
    for x in xrange(L):
        for y in xrange(L):
            i = x + L*y
            A[i,i] = True
            A[n+i,i] = True
            A[(x-1)%L+L*y,i] = True
            A[n+x+(y-1)%2*L,i] = True
            B[i,i] = True
            B[n+i,i] = True
            B[(x+1)%L+L*y,i] = True
            B[n+x+(y+1)%2*L,i] = True
    z = np.zeros((2*n,n),dtype=bool)
    return np.bmat([[z,A],[B,z]])

def show(M):
    (n,m) = np.shape(M)
    for row in xrange(n):
        s = ''
        for col in xrange(m):
            if M[row,col]:
                s += '1'
            else:
                s += ' '
        print(s)

if __name__ == '__main__':
    main()
