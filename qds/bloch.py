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

import hs
import pauli 
import numpy

def vector(evo):
    assert len(evo)
    (time,dop) = evo[0]
    assert dop.shape[0] == dop.shape[1] == 2

    # bloch vector coordinates
    vector = None
    for i in range(0,len(evo)):
        (time,dop) = evo[i]
        assert (dop == numpy.conjugate(numpy.transpose(dop))).all()
        x = hs.dot(pauli.X,dop)
        y = hs.dot(pauli.Y,dop)
        z = hs.dot(pauli.Z,dop)
        entry = numpy.array([x,y,z])
        if vector is None:
            vector = entry
        else:
            vector = numpy.vstack([vector,entry])
            
    return vector


def bloch_sphere(osr_generators):
    assert len(osr_generators)
    for E in osr_generators:
        assert E.shape[0] == E.shape[1] == 2

    import pauli 

    alpha_i = []
    a_i_k = []
    for E in osr_generators:
        # coordinates in the basis I, X, Y e Z
        alpha_i.append(operator.hs_dot(pauli.I,E)/2.)
        a_i_k.append([operator.hs_dot(pauli.X,E)/2., \
                      operator.hs_dot(pauli.Y,E)/2., \
                      operator.hs_dot(pauli.Z,E)/2. ])

    print alpha_i
    print a_i_k
        
    # build M and c
    M = numpy.zeros((3,3))
    c = numpy.zeros((3,1))
    for j in [0,1,2]:
        for k in [0,1,2]:
            for l in range(0,len(osr_generators)):
                M[j][k] = M[j][k] + a_i_k[l][j]*numpy.conj(a_i_k[l][k]) + \
                                    numpy.conj(a_i_k[l][j])*a_i_k[l][k]
                if j == k:
                    abs_alpha = numpy.abs(alpha_i[l]) 
                    M[j][k] = M[j][k] + abs_alpha*abs_alpha
                    for p in [0,1,2]:
                        M[j][k] = M[j][k] - a_i_k[l][p]*numpy.conj(a_i_k[l][p])

                for p in [0,1,2]:
                    e_j_k_p = 0.
                    if (j,k,p) in [(0,1,2),(1,2,0),(2,0,1)]:
                        e_j_k_p = 1.
                    elif (j,k,p) in [(2,1,0),(1,0,2),(0,2,1)]:
                        e_j_k_p = -1.
                    M[j][k] = M[j][k] + e_j_k_p*complex(0, alpha_i[l]*numpy.conj(a_i_k[l][p]) \
                                                            -numpy.conj(alpha_i[l])*a_i_k[l][p])
              
    for k in [0,1,2]:
        for l in range(0,len(osr_generators)):
            for j in [0,1,2]:
                for p in [0,1,2]:
                    e_j_p_k = 0.
                    if (j,p,k) in [(0,1,2),(1,2,0),(2,0,1)]:
                        e_j_p_k = 1.
                    elif (j,p,k) in [(2,1,0),(1,0,2),(0,2,1)]:
                        e_j_p_k = -1.

                    c[k] = c[k] + complex(0,2)*e_j_p_k*a_i_k[l][j]*numpy.conj(a_i_k[l][p])
    print M
    print c
