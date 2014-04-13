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

def dot(A,B):
    """
    Hilbert-Schmidt inner product: <A,B> = tr(adj(A)B)
    """
    assert A.shape[0] == B.shape[0] \
           and A.shape[1] == B.shape[1]

    return numpy.trace(numpy.dot(A.conj().transpose(), B))

def norm(A):
    """
    Norm induced by the Hilbert-Schmidt inner product
    """
    return numpy.sqrt(dot(A,A))

