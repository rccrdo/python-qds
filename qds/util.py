# Copyright (c) 2009 Riccardo Lucchese, riccardo.lucchese at gmail.com
#
# This software is provided 'as-is', without any express or implied
# warranty. In no event will the authors be held liable for any damages
# arising from the use of this software.
#
# Permission is granted to anyone to use this software for any purpose,
# including commercial applications, and to alter it and redistribute it
# freely, subject to the following restrictions:
#
#    1. The origin of this software must not be misrepresented; you must not
#    claim that you wrote the original software. If you use this software
#    in a product, an acknowledgment in the product documentation would be
#    appreciated but is not required.
#
#    2. Altered source versions must be plainly marked as such, and must not be
#    misrepresented as being the original software.
#
#    3. This notice may not be removed or altered from any source
#    distribution.

import numpy
import hs

def comm(A,B):
    """
    commutator [A,B]
    """
    return numpy.dot(A,B) - numpy.dot(B,A)

def acomm(A,B):
    """
    anti-commutator {A,B}
    """
    return numpy.dot(A,B) + numpy.dot(B,A)

def hermitian_subspace_basis(n):
    """
    returns a basis set for the real subspace of Hermitian n*n matrices
    """
    sqrt2over2 = numpy.sqrt(2.)/2.
    assert(n)
    basis = []
    for row in range(0,n):
        for col in range(row,n):
            if row == col:
                mat = numpy.zeros((n,n), dtype=complex)
                mat[row][col] = 1.
                basis.append(mat)
            else:
                mat = numpy.zeros((n,n), dtype=complex)
                mat[col][row] = sqrt2over2
                mat[row][col] = sqrt2over2
                basis.append(mat)                
                mat = numpy.zeros((n,n), dtype=complex)
                mat[row][col] = numpy.complex(0.,sqrt2over2)
                mat[col][row] = numpy.complex(0.,-sqrt2over2)
                basis.append(mat)

    # unit vectors ?
#    for M in basis:
#        for F in basis:
#           print hs.dot(F,M)
        # assert hs.dot(M,M) == 1.
              
    return basis            


def hermitian_subspace_generator(generator_func, n):
    """
    return the generator acting on the real Hermitian subspace of n*n complex matrices
	generator_func  function accepting a n*n matrix and returning the
	                time derivative for the given input
	n               system's order
    """
    assert generator_func
    assert n>0

    basis = hermitian_subspace_basis(n)
    _map = numpy.zeros((n*n,n*n))
    row = 0
    col = 0
    for row in range(0, n*n):
        for col in range(0, n*n):
            _map[row][col] = hs.dot(basis[col], generator_func(basis[row]))

    return _map

