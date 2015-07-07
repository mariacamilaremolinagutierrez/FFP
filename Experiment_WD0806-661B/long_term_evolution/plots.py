import os, math
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from sympy.solvers import solve
from sympy import Symbol

from amuse.support import io
from amuse.units import units
from amuse.datamodel import Particles
from amuse.ext.orbital_elements import orbital_elements_from_binary

def stop_code():
    print '\nSTOP'
    import sys
    sys.exit()

def solve_for_x(m0, m_bp, m_ffp):

    coefficients = [m_bp+m_ffp, 2*m_bp+3*m_ffp, m_bp+3*m_ffp, -(m_bp+3*m0), -(3*m0+2*m_bp), -(m0+m_bp)]
    solutions = np.roots(coefficients)

    return abs(solutions[math.ceil(len(solutions)/2.0)-1])

def is_hill_stable(m_values, a_values, e_values):

    m0 = m_values[0]
    m_ffp = m_values[1]
    m_bp = m_values[2]

    M = m0 + m_bp + m_ffp
    mu = m0 + m_bp

    a_1 = a_values[2]
    a_2 =a_values[1]

    e_1 = e_values[2]
    e_2 = e_values[1]

    x = solve_for_x(m0, m_bp, m_ffp)

    f_x = m0*m_bp + (m0*m_ffp)/(1+x) + (m_bp*m_ffp)/(x)
    g_x = m_bp*m_ffp + m0*m_ffp*(1+x)**2 + m0*m_bp*(x**2)

    A = -(f_x**2)*g_x/(m_ffp**3 * mu**3 * (1-e_2**2))

    if ((1-e_1**2)/(1-e_2**2) < 0.0):
        return False
    if ((a_1*m_ffp*mu)/(a_2*m0*m_bp) < 0.0):
        return False

    beta = (m0*m_bp/m_ffp)**(3.0/2.0)*(M/(mu**4))**(1.0/2.0)*((1-e_1**2)/(1-e_2**2))**(1.0/2.0)
    y = ((a_1*m_ffp*mu)/(a_2*m0*m_bp))**(1.0/2.0)

    equation = (1+y**2)*(beta**2*y**2 + 2*beta*y + 1) - A*y**2

    if(equation >= 0.0):
        return True
    else:
        return False

def plot_trajectory(x, y, filename):

    f = plt.figure(figsize=(70,30))

    x_star = x[:,0]-x[:,0]
    x_ffp = x[:,1]-x[:,0]

    y_star = y[:,0]-y[:,0]
    y_ffp = y[:,1]-y[:,0]

    x_planet = x[:,2]-x[:,0]
    y_planet = y[:,2]-y[:,0]

    plt.plot(x_star,y_star,c='y',label='Star')
    plt.scatter(x_star[0],y_star[0],c='black',marker='*', lw = 0)
    plt.scatter(x_star[-1],y_star[-1],c='y',marker='*', lw = 0)

    plt.plot(x_ffp,y_ffp,c='c',label='FFP', lw = 2)
    plt.scatter(x_ffp[0],y_ffp[0],c='black', lw = 0)
    plt.scatter(x_ffp[-1],y_ffp[-1],c='c', lw = 0)

    plt.plot(x_planet,y_planet,c='m',label='BP',alpha=0.5)
    plt.scatter(x_planet[0],y_planet[0],c='black', lw = 0)
    plt.scatter(x_planet[-1],y_planet[-1],c='m', lw = 0)

    plt.axhline(y=0, xmin=-80, xmax=10, c='black', linestyle='--')
    plt.axvline(x=0, ymin=-5, ymax=2, c='black', linestyle='--')

    plt.title('Trajectory FFP', fontsize=40)
    plt.axes().set_aspect('equal', 'datalim')
    plt.xlabel('$x$ (AU)', fontsize=40)
    plt.ylabel('$y$ (AU)', fontsize=40)
    plt.legend(fontsize=40)
    plt.savefig('./trajectory_'+filename+'.png')

    plt.close()

def plot_orbital_elements(times, eccs, smas, filename):

    f = plt.figure(figsize=(35,15))

    ecc_starbp_ffp = eccs[:,1]
    ecc_star_bp = eccs[:,2]

    sma_starbp_ffp = smas[:,1]
    sma_star_bp = smas[:,2]

    #eccentricity
    subplot = f.add_subplot(1,2,1)

    subplot.plot(times,ecc_starbp_ffp,c='red',label='FFP and Star+BP')
    subplot.plot(times,ecc_star_bp,c='green',label='BP and Star')

    subplot.set_title('Eccentricity', fontsize=20)
    subplot.set_xlabel('$t$ (yr)', fontsize=20)
    subplot.set_ylabel('$e$', fontsize=20)
    subplot.set_ylim(-0.5,1.5)
    subplot.legend(fontsize=20)

    #semimajoraxis
    subplot = f.add_subplot(1,2,2)

    sma_starbp_ffp_log = [math.log10(s) for s in sma_starbp_ffp]
    sma_star_bp_log = [math.log10(s) for s in sma_star_bp]

    subplot.plot(times,sma_starbp_ffp_log,c='red',label='FFP and Star+BP')
    subplot.plot(times,sma_star_bp_log,c='green',label='BP and Star')

    subplot.set_title('Semimajor Axis', fontsize=20)
    subplot.set_xlabel('$t$ (yr)', fontsize=20)
    subplot.set_ylabel('$a$ (AU)', fontsize=20)
    subplot.legend(fontsize=20)

    plt.savefig('./orbital_elements_'+filename+'.png')
    plt.close()

def make_individual_plots(filenames, ms_bp, as_bp, bs_ffp, phis_bp):

    num_files = len(filenames)

    for i in range(num_files):

        filename = filenames[i]
        m_bp = ms_bp[i]
        a_bp = as_bp[i]
        b_ffp = bs_ffp[i]
        phi_bp = phis_bp[i]

        #Bodies Order: star - ffp - bp
        bodies = io.read_set_from_file('./'+filename+'.hdf5', 'hdf5')

        #To append body history
        eccs = []
        smas = []
        xs = []
        ys = []
        times = []

        for data in bodies.history:

            #Order: star - ffp - bp
            e_values = data.eccentricity
            a_values = data.semimajoraxis.value_in(units.AU)
            x_values = data.x.value_in(units.AU)
            y_values = data.y.value_in(units.AU)
            t_value = data.time.value_in(units.yr)[0]

            eccs.append(e_values)
            smas.append(a_values)
            xs.append(x_values)
            ys.append(y_values)
            times.append(t_value)

        plot_trajectory(np.array(xs), np.array(ys), filename)
        plot_orbital_elements(np.array(times), np.array(eccs), np.array(smas), filename)

        print is_hill_stable([607,7.5,m_bp], a_values, e_values)

if __name__ in ('__main__', '__plot__'):

    #This one is at -300 and seem to be going further
    m2 = 1.188889
    a2 = 1.000000
    b2 = -8.023452
    p2 = 275.241762
    f2 = 'm1.188889e+00_a1.000000e+00_b-8.023452e+00_p2.752418e+02'

    #This one starts down but is similar to the second
    m3 = 1.188889
    a3 = 1.000000
    b3 = 11.232833
    p3 = 229.458778
    f3 = 'm1.188889e+00_a1.000000e+00_b1.123283e+01_p2.294588e+02'

    filenames = [f2,f3]
    ms_bp = [m2,m3]
    as_bp = [a2,a3]
    bs_ffp = [b2,b3]
    phis_bp = [p2,p3]

    #Make individual plots
    make_individual_plots(filenames, ms_bp, as_bp, bs_ffp, phis_bp)
