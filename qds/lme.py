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

import time
import numpy
import operator
import util
import hs

from signal import Signal

def hamiltonian_dt(dop, H):
    """
    Hamiltonian part of the LME time derivative
    	dop	system's state as a density operator
	H       system's Hamiltonian
    """
    return numpy.complex(0.,-1.)*util.comm(H, dop)


def lindbladian_dt(dop, Lk):
    """
    lindbladian part of the LME time derivative
    	dop	system's state as a density operator
	Lk      the sequence of Lindblad operators
    """
    dop_dt = numpy.zeros(dop.shape)
    for L in Lk:
        L_adj = L.conj().transpose()
        L_adj_L = numpy.dot(L_adj, L)
        dop_dt = dop_dt + numpy.dot(L,numpy.dot(dop,L_adj)) \
                        -0.5*util.acomm(L_adj_L,dop)
    return dop_dt

def _dt(dop, H, Lk, dt_func_data=None, integrator_time=None):
    """
    LME time derivative
        dop              system's state as a density operator
	H                system's Hamiltonian
	Lk	         the sequence of Lindblad operators
        dt_func_data     unused here
        integrator_time  unused here
    """
    return hamiltonian_dt(dop, H) + lindbladian_dt(dop, Lk)


def _Delta_euler(dt_func, dt_func_data, integrator_time, dop, H, Lk, tstep):
    """
    state 'increment' using Euler's method
        dt_func         function returning the time derivative
        dt_func_data    custom object available used by custom dt_func
        integrator_time current time as reported by the integrator [s]
    	dop             system's state as a density operator
	H               system's Hamiltonian
	Lk	        the sequence of Lindblad operators
        tstep           time length of the integration step [s]
    """
    return tstep*dt_func(dop, H, Lk, dt_func_data, integrator_time)


def _Delta_rk4(dt_func, dt_func_data, integrator_time, dop, H, Lk, tstep):
    """
    state 'increment' using Runge-Kutta's 4th order method
        dt_func         function returning the time derivative
        dt_func_data    custom object available used by custom dt_func
        integrator_time current time as reported by the integrator [s]
    	dop             system's state as a density operator
	H               system's Hamiltonian
	Lk	        the sequence of Lindblad operators
        tstep           time length of the integration step [s]
    """
    k1 = dt_func(dop, H, Lk, dt_func_data, integrator_time)
    k2 = dt_func(dop + 0.5*tstep*k1, H, Lk, dt_func_data, integrator_time)
    k3 = dt_func(dop + 0.5*tstep*k2, H, Lk, dt_func_data, integrator_time)
    k4 = dt_func(dop + tstep*k3, H, Lk, dt_func_data, integrator_time)
    return tstep*(k1/6. + k2/3. + k3/3. + k4/6.)


def integrate(dop_0, H, Lk, tstep, tf, integrator='euler', dt_func=None, dt_func_data=None):
    """
    integrate the Lindblad Master Equation
    	dop_0	      system's initial state as a density operator
	H	      system's Hamiltonian
	Lk	      the sequence of Lindblad operators
    	tstep	      simulation step time
    	tf	      simulation finish time
        integrator    'euler' or 'rk4' for runge-kutta' 4th order method
        dt_func       a function returning the LME time derivative; if defined
                      it will be used in place
                      of the dafult one. Use it to hook
                      different control algorithms into the integrator
                      (ie.  a controller modulating some Hamiltonian
                      terms). ! You can still use the
                      {hamiltonian,lindbladian}_dt functions defined in
                      this module if you need them.
        dt_func_data  accessory data that wil be passed to dt_func
    """
    # check matrices orders
    assert dop_0.ndim == H.ndim == 2
    assert dop_0.shape[0] == dop_0.shape[1] == H.shape[0] == H.shape[1]
    for L in Lk:
        assert L.ndim == dop_0.ndim
        assert L.shape[0] == L.shape[1] == dop_0.shape[0]

    # check matrices props
    assert (dop_0 == dop_0.conj().transpose()).all()
    assert numpy.abs(1. - numpy.trace(dop_0)) < 0.000000001
    assert (H == H.conj().transpose()).all()

    # check timing params
    assert (tstep > 0) and (tf >= tstep)

    # select integrator
    assert integrator in ('euler', 'rk4')
    Delta_func = ({'euler' : _Delta_euler, 'rk4' : _Delta_rk4})[integrator]

    # default to _dt for calculating LME's time derivative
    if not dt_func:
        dt_func = _dt

    # init simulation loop
    max_delta_dop = 0.
    algo_start_time = time.time()
    algo_last_time = algo_start_time

    evo = Signal('Density op. evolution')
    integrator_time = 0.
    evo.append(integrator_time, dop_0)
    integrator_time = integrator_time + tstep 
    dop = dop_0

    while integrator_time <= tf:
        # integrator_time may be needed in state tracking controller
        delta_dop = Delta_func(dt_func, dt_func_data, integrator_time, dop, H, Lk, tstep)
        new_max_delta_dop = numpy.max(numpy.abs(delta_dop))
        if new_max_delta_dop > max_delta_dop:
            max_delta_dop = new_max_delta_dop
        dop = dop + delta_dop
        #TODO: should print out a measure of the 'drift' from 'hermitianicity'
        #assert (dop == dop.conj().transpose()).all()
        #assert numpy.trace(dop) == 1
        evo.append(integrator_time,dop)
        integrator_time = integrator_time + tstep 
 
        if (time.time() - algo_last_time) > 20.:
            algo_last_time = time.time()
            progress = integrator_time/float(tf)*100.
            ETA_min = (algo_last_time - algo_start_time) \
                            /integrator_time*(tf - integrator_time)/60.
            print 'in lme.sim: %d%%, ETA: %.1fmin' % (progress, ETA_min)

    if max_delta_dop > 0.001:
        print ' ! warning in lme.sim, max_delta_dop: ', max_delta_dop

    return evo

