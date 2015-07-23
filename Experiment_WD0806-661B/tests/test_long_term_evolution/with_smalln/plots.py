import os, math
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from sympy.solvers import solve
from sympy import Symbol

from amuse.support import io
from amuse.units import units,constants
from amuse.datamodel import Particles

def stop_code():
    print '\nSTOP'
    import sys
    sys.exit()

def plot_trajectory(x, y, filename, x_start_longterm, y_start_longterm):

    f = plt.figure(figsize=(70,30))

    x_star = x[:,0]-x[:,0]
    x_ffp = x[:,1]-x[:,0]
    x_planet = x[:,2]-x[:,0]

    y_star = y[:,0]-y[:,0]
    y_ffp = y[:,1]-y[:,0]
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

    plt.scatter(x_start_longterm[1],y_start_longterm[1], c='c', lw = 0, marker='s')
    plt.scatter(x_start_longterm[2],y_start_longterm[2], c='m', lw = 0, marker='s')

    plt.axhline(y=0, xmin=-80, xmax=10, c='black', linestyle='--')
    plt.axvline(x=0, ymin=-5, ymax=2, c='black', linestyle='--')

    plt.title('Trajectory FFP', fontsize=40)
    plt.axes().set_aspect('equal', 'datalim')
    plt.xlabel('$x$ (AU)', fontsize=40)
    plt.ylabel('$y$ (AU)', fontsize=40)
    plt.legend(fontsize=40)
    plt.savefig('./trajectory_'+filename+'.png')

    plt.close()

def plot_orbital_elements(times_array, e_array, a_array, i_array, l_array, g_array, fout):

    #BP
    figure = plt.figure(figsize=(25,20))
    N_subplots = 5

    plot_e = figure.add_subplot(N_subplots,1,1)
    plot_a = figure.add_subplot(N_subplots,1,2)
    plot_i = figure.add_subplot(N_subplots,1,3)
    plot_l = figure.add_subplot(N_subplots,1,4)
    plot_g = figure.add_subplot(N_subplots,1,5)

    plot_e.plot(times_array,e_array[:,0],c='c',lw=2)
    a_array_subset = a_array[:,0]
    log10_a_array = [math.log10(aa) for aa in a_array_subset]
    plot_a.plot(times_array,log10_a_array,c='c',lw=2)
    plot_i.plot(times_array,i_array[:,0]*180.0/np.pi,c='c',lw=2)
    plot_l.plot(times_array,l_array[:,0]*180.0/np.pi,c='c',lw=2)
    plot_g.plot(times_array,g_array[:,0]*180.0/np.pi,c='c',lw=2)

    t_max_yr = max(times_array)
    plot_e.set_xlim(0,t_max_yr)
    plot_a.set_xlim(0,t_max_yr)
    plot_i.set_xlim(0,t_max_yr)
    plot_l.set_xlim(0,t_max_yr)
    plot_g.set_xlim(0,t_max_yr)

    plot_e.set_xlabel('$t/\mathrm{yr}$')
    plot_a.set_xlabel('$t/\mathrm{yr}$')
    plot_i.set_xlabel('$t/\mathrm{yr}$')
    plot_l.set_xlabel('$t/\mathrm{yr}$')
    plot_g.set_xlabel('$t/\mathrm{yr}$')

    plot_e.set_ylabel('$e_\mathrm{BP}$')
    plot_a.set_ylabel('$\log{(a_\mathrm{BP})}$')
    plot_i.set_ylabel('$i_\mathrm{BP} ({}^\circ)$')
    plot_l.set_ylabel('$l_\mathrm{BP} ({}^\circ)$')
    plot_g.set_ylabel('$g_\mathrm{BP} ({}^\circ)$')

    figure.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)

    plt.savefig('./orbital_elements_bp_'+fout+'.png')
    plt.close()

    #FFP
    figure = plt.figure(figsize=(25,20))
    N_subplots = 5

    plot_e = figure.add_subplot(N_subplots,1,1)
    plot_a = figure.add_subplot(N_subplots,1,2)
    plot_i = figure.add_subplot(N_subplots,1,3)
    plot_l = figure.add_subplot(N_subplots,1,4)
    plot_g = figure.add_subplot(N_subplots,1,5)

    plot_e.plot(times_array,e_array[:,1],c='c',lw=2)
    a_array_subset = a_array[:,1]
    log10_a_array = [math.log10(aa) for aa in a_array_subset]
    plot_a.plot(times_array,log10_a_array,c='c',lw=2)
    plot_i.plot(times_array,i_array[:,1]*180.0/np.pi,c='c',lw=2)
    plot_l.plot(times_array,l_array[:,1]*180.0/np.pi,c='c',lw=2)
    plot_g.plot(times_array,g_array[:,1]*180.0/np.pi,c='c',lw=2)

    t_max_yr = max(times_array)
    plot_e.set_xlim(0,t_max_yr)
    plot_a.set_xlim(0,t_max_yr)
    plot_i.set_xlim(0,t_max_yr)
    plot_l.set_xlim(0,t_max_yr)
    plot_g.set_xlim(0,t_max_yr)

    plot_e.set_xlabel('$t/\mathrm{yr}$')
    plot_a.set_xlabel('$t/\mathrm{yr}$')
    plot_i.set_xlabel('$t/\mathrm{yr}$')
    plot_l.set_xlabel('$t/\mathrm{yr}$')
    plot_g.set_xlabel('$t/\mathrm{yr}$')

    plot_e.set_ylabel('$e_\mathrm{FFP}$')
    plot_a.set_ylabel('$\log{(a_\mathrm{FFP})} (\mathrm{AU})$')
    plot_i.set_ylabel('$i_\mathrm{FFP} ({}^\circ)$')
    plot_l.set_ylabel('$l_\mathrm{FFP} ({}^\circ)$')
    plot_g.set_ylabel('$g_\mathrm{FFP} ({}^\circ)$')

    figure.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)

    plt.savefig('./orbital_elements_ffp_'+fout+'.png')

# def plot_orbital_elements(times, eccs, smas, filename):
#
#     f = plt.figure(figsize=(35,15))
#
#     ecc_starbp_ffp = eccs[:,1]
#     ecc_star_bp = eccs[:,2]
#
#     sma_starbp_ffp = smas[:,1]
#     sma_star_bp = smas[:,2]
#
#     #eccentricity
#     subplot = f.add_subplot(1,2,1)
#
#     subplot.plot(times,ecc_starbp_ffp,c='red',label='FFP and Star+BP')
#     subplot.plot(times,ecc_star_bp,c='green',label='BP and Star')
#
#     subplot.set_title('Eccentricity', fontsize=20)
#     subplot.set_xlabel('$t$ (yr)', fontsize=20)
#     subplot.set_ylabel('$e$', fontsize=20)
#     subplot.set_ylim(-0.5,1.5)
#     subplot.legend(fontsize=20)
#
#     #semimajoraxis
#     subplot = f.add_subplot(1,2,2)
#
#     sma_starbp_ffp_log = []
#
#     i=0
#     for s in sma_starbp_ffp:
#         if(i!=0):
#             try:
#                 sma_starbp_ffp_log.append(math.log10(s))
#             except:
#                 # print i
#                 sma_starbp_ffp_log.append(sma_starbp_ffp_log[-1])
#         i += 1
#
#     sma_star_bp_log = []
#
#     i=0
#     for s in sma_star_bp:
#         if(i!=0):
#             try:
#                 sma_star_bp_log.append(math.log10(s))
#             except:
#                 # print i
#                 sma_star_bp_log.append(sma_star_bp_log[-1])
#         i += 1
#
#     subplot.plot(times[1:],sma_starbp_ffp_log,c='red',label='FFP and Star+BP')
#     subplot.plot(times[1:],sma_star_bp_log,c='green',label='BP and Star')
#
#     subplot.set_title('Semimajor Axis', fontsize=20)
#     subplot.set_xlabel('$t$ (yr)', fontsize=20)
#     subplot.set_ylabel('$\log{(a)}$', fontsize=20)
#     subplot.legend(fontsize=20)
#
#     plt.savefig('./orbital_elements_'+filename+'.png')
#     plt.close()

def my_orbital_elements_from_binary(mass1, mass2, binary):

    position = binary[1].position-binary[0].position
    velocity = binary[1].velocity-binary[0].velocity
    total_mass = mass1 + mass2

    specific_energy = (1.0/2.0)*velocity.lengths_squared() - constants.G*total_mass/position.lengths()
    specific_angular_momentum = position.cross(velocity)
    specific_angular_momentum_norm = specific_angular_momentum.lengths()
    specific_angular_momentum_unit=specific_angular_momentum/specific_angular_momentum_norm

    semimajor_axis = -constants.G*total_mass/(2.0*specific_energy)

    eccentricity_argument = 2.0*specific_angular_momentum_norm**2*specific_energy/(constants.G**2*total_mass**2)
    if (eccentricity_argument <= -1): eccentricity = 0.0
    else: eccentricity = np.sqrt(1.0 + eccentricity_argument)

    ### Orbital inclination ###
    inclination = np.degrees(np.arccos(specific_angular_momentum.z/specific_angular_momentum_norm))

    ### Longitude of ascending nodes, with reference direction along x-axis ###
    z_vector = [0.,0.,1.] | units.none
    ascending_node_vector = z_vector.cross(specific_angular_momentum)
    if ascending_node_vector.lengths().number==0:
      ascending_node_vector_unit= np.array([1.,0.,0.])
    else:
      ascending_node_vector_unit = ascending_node_vector/ascending_node_vector.lengths()
    long_asc_node=np.degrees(np.arctan2(ascending_node_vector_unit[1],ascending_node_vector_unit[0]))

    ### Argument of periapsis and true anomaly, using eccentricity a.k.a. Laplace-Runge-Lenz vector ###
    mu = constants.G*total_mass ### Argument of pericenter ###
    position_unit = position/position.lengths()
    e_vector = ( (1.0/mu)*velocity.cross(specific_angular_momentum) - position_unit ) | units.none
    if (e_vector.lengths() == 0.0): ### Argument of pericenter and true anomaly cannot be determined for e = 0, in this case return 1.0 for the cosines ###
        cos_arg_per = 1.0
        arg_per=0.
        cos_true_anomaly = 1.0
        true_anomaly=0.
    else:
        e_vector_unit = e_vector/e_vector.lengths()

        cos_arg_per = np.dot(e_vector_unit,ascending_node_vector_unit)
        #cos_arg_per = e_vector_unit.dot(ascending_node_vector_unit)
        e_cross_an=np.cross(e_vector_unit,ascending_node_vector_unit)
        ss=-np.sign(np.dot(specific_angular_momentum_unit,e_cross_an))
        #ss=-np.sign(specific_angular_momentum_unit.dot(e_cross_an))
        sin_arg_per = ss*(e_cross_an**2).sum()**0.5
        arg_per=np.degrees(np.arctan2(sin_arg_per,cos_arg_per))


        cos_true_anomaly = np.dot(e_vector_unit,position_unit)
        #cos_true_anomaly = e_vector_unit.dot(position_unit)
        e_cross_pos=np.cross(e_vector_unit,position_unit)
        ss=np.sign(np.dot(specific_angular_momentum_unit,e_cross_pos))
        #ss=np.sign(specific_angular_momentum_unit.dot(e_cross_pos))
        sin_true_anomaly = ss*(e_cross_pos**2).sum()**0.5
        true_anomaly=np.degrees(np.arctan2(sin_true_anomaly,cos_true_anomaly))

    return mass1, mass2, semimajor_axis, eccentricity, true_anomaly, inclination, long_asc_node, arg_per

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
        incs = []
        lans = []
        aps = []

        xs = []
        ys = []
        times = []

        snapshot = 1

        x_start_longterm = [0.0,0.0,0.0]
        y_start_longterm = [0.0,0.0,0.0]

        for data in bodies.history:

            #Order: star - ffp - bp
            e_values = data.eccentricity
            a_values = data.semimajoraxis.value_in(units.AU)
            i_values = data.inclination
            l_values = data.longitude_of_ascending_node
            g_values = data.argument_of_periapsis

            x_values = data.x.value_in(units.AU)
            y_values = data.y.value_in(units.AU)
            t_value = data.time.value_in(units.yr)[0]

            eccs.append(e_values)
            smas.append(a_values)
            incs.append(i_values)
            lans.append(l_values)
            aps.append(g_values)

            xs.append(x_values)
            ys.append(y_values)
            times.append(t_value)

            if(snapshot==600):
                x_start_longterm = x_values
                y_start_longterm = y_values

            snapshot += 1

        plot_trajectory(np.array(xs), np.array(ys), filename, x_start_longterm, y_start_longterm)
        plot_orbital_elements(np.array(times), np.array(eccs), np.array(smas), np.array(incs), np.array(lans), np.array(aps), filename)

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
