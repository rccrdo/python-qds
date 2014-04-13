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
from pylab import *

from signal import Signal

def trajectories(signal, of=None):
    assert len(signal)

    # only plot trajectories for the upper triangular part
    timeline = signal.timeline()
    trajs = signal.upper_triang_trajectories()

    figure(figsize=(20,7))      

#    title('traiettorie dei posti in \rho')     
#    ylabel('stato')
#    xlabel('tempo [s]')
    n = signal.shape()[0]
    num_plots = len(trajs) 
    num_rows = int(numpy.ceil(float(num_plots)/n))
    row = col = 1
    for i in range(0, num_plots):
        if col != col%(n+1):
            row = row +1
            col = row

        subplot(num_rows, n, i+1)
#        real_traj =  numpy.real(trajctories[i])
#        print numpy.real(trajectories[i])
        plot(timeline, numpy.real(trajs[i]), 'k')
        hold(True)
 
        # special case text and plot switch for elements on the diag.
        head = '('
        tail = ')'
        if row != col:
            plot(timeline, numpy.imag(trajs[i]), 'k--')
        else:
            head = '<'
            tail = '>'

        title(''.join([head,str(row),',',str(col),tail]))
        ylim((0,1))
        grid(True)
        col = col + 1

    draw()
    if not of:
        show()
    else:
        savefig(of)

    close()


def _gen_exp(start_value, end_value, p_0, tf, tstep):
    exp_evo = []
    exp_evo.append(start_value)
    exp_time = tstep

    while exp_time < tf:
            value = (end_value - start_value)*(1. -numpy.exp(p_0*exp_time)) + start_value
            exp_evo.append(value)
            exp_time = exp_time + tstep

#    print '\n****', start_value, end_value, p_0, tf, tstep
#    print exp_evo
    return exp_evo


def trajectories_with_bound(signal, z_0, of=None):
    assert len(signal)

    # only plot trajectories for the upper triangular part
    timeline = signal.timeline()
    tstep = timeline[1] - timeline[0]
    tf = timeline[len(timeline) -1] + tstep

    trajs = signal.upper_triang_trajectories()

    figure(figsize=(20,7))      

#    title('traiettorie dei posti in \rho')     
#    ylabel('stato')
#    xlabel('tempo [s]')
    n = signal.shape()[0]
    num_plots = len(trajs) 
    num_rows = int(numpy.ceil(float(num_plots)/n))
    row = col = 1
    for i in range(0, num_plots):
        if col != col%(n+1):
            row = row +1
            col = row

        subplot(num_rows, n, i+1)
#        real_traj =  numpy.real(trajctories[i])
#        print numpy.real(trajectories[i])
        plot(timeline, numpy.real(trajs[i]), 'k')
        hold(True)
 
        # special case text and plot switch for elements on the diag.
        head = '('
        tail = ')'
        if row != col:
            plot(timeline, numpy.imag(trajs[i]), 'k--')
        else:
            head = '<'
            tail = '>'
            # plot bounds too 
            if row==1:
                exp_bound = _gen_exp(trajs[i][0], 1., z_0, tf, tstep)
                plot(timeline, exp_bound, 'r')

        title(''.join([head,str(row),',',str(col),tail]))
        ylim((0,1))
        grid(True)
        col = col + 1

    draw()
    if not of:
        show()
    else:
        savefig(of)

    close()


def mag_bars(signal, step=1):
    tmp_file_suffix = '/tmp/qds_fig'
    file_cleanup_list = []

    from mpl_toolkits.mplot3d import Axes3D
    from matplotlib.collections import PolyCollection
    from matplotlib.colors import colorConverter
    import pylab
    import random
    import numpy as np

    (time, dop) = signal[0]
    xpos, ypos = np.meshgrid(range(0,dop.shape[0]), range(0,dop.shape[1]))
    xpos = xpos.flatten()
    ypos = ypos.flatten()
    zpos = [0] * len(xpos)
    dx = [0.5] * len(xpos)
    dy = [0.5] * len(xpos)

#    pylab.show()
    for i in range(0, len(signal), step):
        print 'frame ' + str(i) +' of ' + str(len(signal))
        tmp_filename = ''.join([tmp_file_suffix,str(i), '.png'])

        (time, dop) = signal[i]
        assert dop.shape[0] == dop.shape[1]
        fig = figure(figsize=(2, 2)) 
        ax = Axes3D(fig)
        dz = numpy.abs(dop.flatten())
        ax.bar3d(xpos, ypos, zpos, dx, dy, dz, color='b')
        ax.set_zlim3d([0., 1.])
        pylab.draw()

        savefig(tmp_filename, format='png')
        file_cleanup_list.append(tmp_filename)
        close()

def traj(signal):
    assert len(signal)
    (time, dop) = signal[0]

    if dop.shape[1] == 1:
        dop = numpy.reshape(dop, (numpy.sqrt(dop.shape[0]), numpy.sqrt(dop.shape[0])), order='F')

    colors = []
    num_colors = dop.shape[0]*(dop.shape[0] +1)/2
    #print num_colors
    for n in range(0, num_colors):
        colors.append( (numpy.random.uniform(low=0., high=1.), 
                        numpy.random.uniform(low=0., high=1.), 
                        numpy.random.uniform(low=0., high=1.)) ) 

    figure(figsize=(4, 4))                 
    axis([-2, 2, -2, 2])
    hold(True)
    grid(True)
    for row in range(0, dop.shape[0]):
        for col in range(0, dop.shape[1]):
            if col >= row:
                num_color = row*dop.shape[0] + col  - row*(row +1)/2
            else:
                continue            
                num_color = col*dop.shape[1] + row  - col*(col +1)/2

            print 'traj  [' + str(row) +', ' + str(col) + ']'
            traj_re = []
            traj_im = []
            for i in range(0, len(signal)):
                (time, dop) = signal[i]
                if dop.shape[1] == 1:
                    dop = numpy.reshape(dop, (numpy.sqrt(dop.shape[0]), numpy.sqrt(dop.shape[0])), order='F')
                val = dop[row][col]
                traj_re.append(numpy.real(val))
                traj_im.append(numpy.imag(val))

            plot(traj_re, traj_im)

    show()


def bloch_vector(signal):
    assert len(signal)
    (time,dop) = signal[0]
    assert dop.shape[0] == dop.shape[1] == 2

    import pauli 
    import hs

    x = []
    y = []
    z = []
    for i in range(0,len(signal)):
        (time,dop) = signal[i]
        # coordinates in the basis I, X, Y e Z
        i_tmp = hs.dot(pauli.I,dop)
        x_tmp = hs.dot(pauli.X,dop)
        y_tmp = hs.dot(pauli.Y,dop)
        z_tmp = hs.dot(pauli.Z,dop)
#        print i_tmp, x_tmp, y_tmp, z_tmp
#        assert numpy.imag(x_tmp) == numpy.imag(y_tmp) == numpy.imag(z_tmp) == 0.
#        print 'dop',  dop
#        print 'sum',  (pauli.I*i_tmp + x_tmp*pauli.X+y_tmp*pauli.Y+z_tmp*pauli.Z)/2.
#        assert (pauli.I*i_tmp + x_tmp*pauli.X+y_tmp*pauli.Y+z_tmp*pauli.Z == dop).all()
        x.append(numpy.real(x_tmp))
        y.append(numpy.real(y_tmp))
        z.append(numpy.real(z_tmp))

#    x = numpy.array(x)   
#    y = numpy.array(y)   
#    z = numpy.array(z)   

#    print x
#    print y
#    print z
#   print x.shape
#   print y.shape
#   print z.shape

    from mpl_toolkits.mplot3d import Axes3D
    import pylab

    mpl.rcParams['legend.fontsize'] = 10
    
    fig = pylab.figure()
    ax = Axes3D(fig)

    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)

    sx = 1 * np.outer(np.cos(u), np.sin(v))
    sy = 1 * np.outer(np.sin(u), np.sin(v))
    sz = 1 * np.outer(np.ones(np.size(u)), np.cos(v))
    ax.plot_surface(sx, sy, sz,  rstride=25, cstride=25, color='b', alpha=0.05)
    
    pylab.hold(True)


    ax.plot(x, y, z, label='parametric curve')
    ax.legend()

    pylab.show()
    
def bloch_vector2(signal):
    assert len(signal)
    (time,dop) = signal[0]
    assert dop.shape[0] == dop.shape[1] == 2

    import pauli 
    import hs
    import bloch
    vector = bloch.vector(signal)

    from mpl_toolkits.mplot3d import Axes3D
    import pylab

    mpl.rcParams['legend.fontsize'] = 10
    
    fig = pylab.figure()
    ax = Axes3D(fig)

    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)

    sx = 1 * np.outer(np.cos(u), np.sin(v))
    sy = 1 * np.outer(np.sin(u), np.sin(v))
    sz = 1 * np.outer(np.ones(np.size(u)), np.cos(v))
    ax.plot_surface(sx, sy, sz,  rstride=25, cstride=25, color='b', alpha=0.05)
    
    pylab.hold(True)


    ax.plot(vector[:,0], vector[:,1], vector[:,2], label='bloch vector')
    ax.legend()

    pylab.show()
