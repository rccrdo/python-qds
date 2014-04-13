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
from signal import Signal


def _gen_exp(start_value, end_value, lambda_exp, tf, tstep):
    """
    Returns a trajectory of the form (end_value - start_value)*(1 - e^(lambda_exp*t)) + start_value
    """
    exp_evo = []
    exp_evo.append(start_value)
    exp_time = tstep

    while exp_time < tf:
            value = (end_value - start_value)*(1. -numpy.exp(lambda_exp*exp_time)) + start_value
            exp_evo.append(value)
            exp_time = exp_time + tstep

#    print '\n****', start_value, end_value, p_0, tf, tstep
#    print exp_evo
    return exp_evo

def _stream_columns(data_name, columns):
    assert data_name
    buffer = data_name + ' = ['
    line_start_space = ';\n' + ' '*len(buffer)

    # case where columns is just one list
    try:
        col_size = len(columns[0])
        cols = columns
    except:
        cols = [columns]
        col_size = len(cols[0])

    i = 0
    _col_distance = 22
    while (i < col_size):
        line = ''
        if i:
            line = line_start_space
        for col in cols:

            value = col[i]
            value_re = numpy.real(value)
            value_im = numpy.imag(value)

            if value_im:
                if value_re:
                    value_str = '%.6f%+.6fi' % (value_re, value_im)
                else:
                    value_str = '%.6fi' % value_im
            else:
                value_str = '%.6f' % value_re

            value_str_len = len(value_str)
#            print value_str, len(value_str)
            assert (value_str_len < _col_distance)

            line += ' '*(_col_distance - value_str_len) + value_str

        i = i +1
        buffer += line   

    buffer += ' ];'
    return buffer


def stream_trajectories_with_bound(signal, bound_exp):
    assert len(signal)

    # only plot trajectories for the upper triangular part
    timeline = signal.timeline()
    tstep = timeline[1] - timeline[0]
    tf = timeline[len(timeline) -1] + tstep

    trajs = signal.upper_triang_trajectories()
    buffer = _stream_columns('dop_data', trajs)
    exp_traj = _gen_exp(0., 1., bound_exp, tf, tstep)
    buffer += '\n\n'
    buffer += _stream_columns('bound_data', exp_traj)

    return buffer


_COL_DISTANCE = 16

def _value_str(value):
    value_re = numpy.real(value)
    value_im = numpy.imag(value)

    if value_im:
        if value_re:
            value_str = '%.4f%+.4fi' % (value_re, value_im)
        else:
            value_str = '%.4fi' % value_im
    else:
        value_str = '%.4f' % value_re

    return value_str

def stream_numpy_matrix(name, matrix):
    assert name

    buffer = name + ' = ['
    line_start_space = ';\n' + ' '*len(buffer)

    # handle vectors
    if len(matrix.shape)== 1:
#        print "doing vector: %s\n" % name
        for row in range(0, matrix.shape[0]):
            line = ''
            if row:
                line = line_start_space
            value_str = _value_str(matrix[row])
            value_str_len = len(value_str)
            #            print value_str, len(value_str)
            assert (value_str_len < _COL_DISTANCE)
            line += ' '*(_COL_DISTANCE - value_str_len) + value_str
            buffer += line   
    else:
#        print "doing matrix: %s\n" % name
        for row in range(0, matrix.shape[0]):
            line = ''
            if row:
                line = line_start_space

            for col in range(0, matrix.shape[1]):
                value_str = _value_str(matrix[row][col])

                value_str_len = len(value_str)
                #            print value_str, len(value_str)
                assert (value_str_len < _COL_DISTANCE)
                line += ' '*(_COL_DISTANCE - value_str_len) + value_str
            buffer += line   

    buffer += ' ];\n'
    return buffer
