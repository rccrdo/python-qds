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
import lme

def corresponding_lme_operators(H, M, F):
    """
    converts the H, M, F operators used in the FME to the H,L
    operators to be used in the LME
    H      system's Hamiltonian
    M      measurement operator
    F      feedback operator

    returns a tuple (H_lme, L_lme) with the Hamiltonian and the noise
    operator to be used to integrate the associated LME
    """
    # check matrices orders
    assert H.ndim == M.ndim == F.ndim == 2
    assert H.shape[0] == H.shape[1] == M.shape[0] == M.shape[1] == F.shape[0] == F.shape[1]

    # check matrices props
    assert (H == H.conj().transpose()).all()

    FM = numpy.dot(F,M)
    Madj_F = numpy.transpose(numpy.conjugate(FM))
    H_lme = H + 0.5*(FM + Madj_F)
    L_lme = M - numpy.complex(0., -1.)*F
    return (H_lme, L_lme)


def integrate(dop_0, H, M, F, tstep, tf, integrator='euler'):
    """
    integrate the Wiseman-Wilburn Markovian Feedback Master equation
    	dop_0	      system's initial state as a density operator
	H	      system's Hamiltonian
	M             measurement operator
        F             feedback operator
    	tstep	      simulation step time
    	tf	      simulation finish time
        integrator    'euler' or 'rk4' for runge-kutta' 4th order method
    """
    # matrices checks are done in corresponding_lme_operators()
    # and in lme.integrate()
    
    (H_lme, L_lme) = corresponding_lme_operators(H, M, F)
    print  (H_lme, L_lme)
    return lme.integrate(dop_0,  H_lme, (L_lme), tstep, tf)

