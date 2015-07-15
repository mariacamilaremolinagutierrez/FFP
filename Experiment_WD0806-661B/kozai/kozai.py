from amuse.community.seculartriple.interface import SecularTriple
from amuse.datamodel import Particles
from amuse.units import units,constants,quantities
from amuse.units.optparse import OptionParser
from amuse.ext.orbital_elements import new_binary_from_orbital_elements

import matplotlib.pyplot as plt
import numpy as np
import math

def make_triple_system(m_star=0.58|units.MSun, m_bp=0.1|units.MJupiter, m_ffp=7.5|units.MJupiter, a_star_bp=1.0|units.AU, a_star_ffp=1540.5|units.AU, e_star_bp=0.0, e_star_ffp=0.998, om_star_bp=0.0, om_star_ffp=0.0, inc=0.0):

    binaries = Particles(2)
    inner_stars = Particles(2)
    outer_stars = Particles(2)

    binaries[0].child1 = inner_stars[0]
    binaries[0].child2 = inner_stars[1]
    binaries[1].child1 = outer_stars[0]
    binaries[1].child2 = outer_stars[1]

    binaries[0].child1.mass = m_star
    binaries[0].child2.mass = m_bp
    binaries[1].child1.mass = m_ffp

    binaries[0].semimajor_axis = a_star_bp
    binaries[1].semimajor_axis = a_star_ffp

    binaries[0].eccentricity = e_star_bp
    binaries[1].eccentricity = e_star_ffp
    binaries[0].argument_of_pericenter = om_star_bp*np.pi/180.0
    binaries[1].argument_of_pericenter = om_star_ffp*np.pi/180.0
    binaries[0].longitude_of_ascending_node = 0.0
    binaries[1].longitude_of_ascending_node = 0.0

    binaries[0].inclination = inc*np.pi/180.0 ### this is the inclination between the orbital planes of binaries[0] and binaries[1] ###
    binaries[1].inclination = 0.0

    return binaries

def give_position(semimajor_axis, eccentricity, inclination, argument_of_periapsis):

    inclination = np.radians(inclination)
    argument_of_periapsis = np.radians(argument_of_periapsis)

    #Supposing this:
    #longitude_of_the_ascending_node = 0.0
    #true_anomaly = 0.0

    cos_inclination = np.cos(inclination)
    sin_inclination = np.sin(inclination)

    cos_arg_per = np.cos(argument_of_periapsis)
    sin_arg_per = np.sin(argument_of_periapsis)

    ### alpha is a unit vector directed along the line of node ###
    alphax = cos_arg_per
    alphay = sin_arg_per*cos_inclination
    alphaz = sin_arg_per*sin_inclination

    alpha = [alphax,alphay,alphaz]

    ### Relative position and velocity ###
    separation = semimajor_axis*(1.0 - eccentricity**2)/(1.0 + eccentricity) # Compute the relative separation
    position_vector = separation*alpha

    return position_vector[0], position_vector[1], position_vector[2]

def evolve_triple_system(binaries,end_time, output_time_step, fout):

    code = SecularTriple()
    code.binaries.add_particles(binaries)
    code.parameters.equations_of_motion_specification = 0
    code.parameters.f_quad = 1.0
    code.parameters.f_oct = 1.0
    code.parameters.f_mass_transfer = 0.0 #leave off
    code.parameters.f_1PN_in = 0.0 #general relativity specifications
    code.parameters.f_1PN_out = 0.0 #general relativity specifications
    code.parameters.f_25PN_in = 0.0 #general relativity specifications
    code.parameters.f_25PN_out = 0.0 #general relativity specifications

    ### quantities later used for plotting ###
    x = quantities.AdaptingVectorQuantity()
    y = quantities.AdaptingVectorQuantity()
    z = quantities.AdaptingVectorQuantity()

    times_array = []
    e_array = []
    a_array = []
    g_array = []
    i_array = []

    print '\nBINARIES AT START:'
    print code.binaries

    time = 0.0 | units.yr
    while (time < end_time):

        code.evolve_model(time)
        time += output_time_step

        ecc_star_bp = code.binaries[0].eccentricity
        ecc_star_ffp = code.binaries[1].eccentricity
        sma_star_bp = code.binaries[0].semimajor_axis
        sma_star_ffp = code.binaries[1].semimajor_axis
        ap_star_bp = code.binaries[0].argument_of_pericenter
        ap_star_ffp = code.binaries[1].argument_of_pericenter
        inc_star_bp = code.binaries[0].inclination
        inc_star_ffp = code.binaries[1].inclination

        # print ecc_star_bp, ecc_star_ffp, sma_star_bp, sma_star_ffp, ap_star_bp, ap_star_ffp, inc_star_bp, inc_star_ffp

        times_array.append(time.value_in(units.yr))
        e_array.append([ecc_star_bp, ecc_star_ffp])
        a_array.append([sma_star_bp.value_in(units.AU), sma_star_ffp.value_in(units.AU)])
        g_array.append([ap_star_bp, ap_star_ffp])
        i_array.append([inc_star_bp, inc_star_ffp])

        x_star_bp,y_star_bp,z_star_bp = give_position(sma_star_bp, ecc_star_bp, inc_star_bp, ap_star_bp)
        x_star_ffp,y_star_ffp,z_star_ffp = give_position(sma_star_ffp, ecc_star_ffp, inc_star_ffp, ap_star_ffp)

        xs = [0.0, x_star_bp.value_in(units.AU), x_star_ffp.value_in(units.AU)] |units.AU
        ys = [0.0, y_star_bp.value_in(units.AU), y_star_ffp.value_in(units.AU)] |units.AU
        zs = [0.0, z_star_bp.value_in(units.AU), z_star_ffp.value_in(units.AU)] |units.AU

        x.append(xs)
        y.append(ys)
        z.append(zs)

    print '\nBINARIES AT END:'
    print code.binaries

    # m_star = 0.58 | units.MSun
    # m_bp = 1.188889 | units.MJupiter
    # m_ffp = 7.5|units.MJupiter
    #
    # binnn = new_binary_from_orbital_elements(m_star, m_bp, sma_star_bp, ecc_star_bp, 0.0, inc_star_bp, 0.0, ap_star_bp, G=units.constants.G)
    # binnn2 = new_binary_from_orbital_elements(m_star, m_ffp, sma_star_ffp, ecc_star_ffp, 0.0, inc_star_ffp, 0.0, ap_star_ffp, G=units.constants.G)
    # print '\nX:'
    # print xs
    # print (binnn[1].x-binnn[0].x).as_quantity_in(units.AU)
    # print (binnn2[1].x-binnn2[0].x).as_quantity_in(units.AU)
    # print '\nY:'
    # print ys
    # print (binnn[1].y-binnn[0].y).as_quantity_in(units.AU)
    # print (binnn2[1].y-binnn2[0].y).as_quantity_in(units.AU)
    # print '\nZ:'
    # print zs
    # print (binnn[1].z-binnn[0].y).as_quantity_in(units.AU)
    # print (binnn2[1].z-binnn2[0].y).as_quantity_in(units.AU)

    # fig=plt.figure()
    # plt.plot(times_array.value_in(units.yr), g_in_array)
    # plt.plot(times_array.value_in(units.yr), g1_in_array)
    # plt.show()
    # plt.close()

    times_array = np.array(times_array)
    e_array = np.array(e_array)
    a_array = np.array(a_array)
    g_array = np.array(g_array)
    i_array = np.array(i_array)

    return times_array,e_array,a_array,g_array,i_array,x,y,z

def plot_trajectory(x, y, fout):

    f = plt.figure(figsize=(70,30))

    x_star = (x[:,0]-x[:,0]).value_in(units.AU)
    x_planet = (x[:,1]-x[:,0]).value_in(units.AU)
    x_ffp = (x[:,2]-x[:,0]).value_in(units.AU)

    y_star = (y[:,0]-y[:,0]).value_in(units.AU)
    y_planet = (y[:,1]-y[:,0]).value_in(units.AU)
    y_ffp = (y[:,2]-y[:,0]).value_in(units.AU)

    plt.plot(x_star,y_star,c='y',label='Star')
    plt.scatter(x_star[0],y_star[0],c='black',marker='*', lw = 0)
    plt.scatter(x_star[-1],y_star[-1],c='y',marker='*', lw = 0)

    plt.plot(x_ffp,y_ffp,c='c',label='FFP', lw = 2)
    plt.scatter(x_ffp[0],y_ffp[0],c='black', lw = 0)
    plt.scatter(x_ffp[-1],y_ffp[-1],c='c', lw = 0)

    plt.plot(x_planet,y_planet,c='m',label='BP',alpha=0.5)
    plt.scatter(x_planet[0],y_planet[0],c='black', lw = 0)
    plt.scatter(x_planet[-1],y_planet[-1],c='m', lw = 0)

    plt.axhline(y=0, xmin=-80, xmax=10, c='black', linestyle='--', alpha=0.5)
    plt.axvline(x=0, ymin=-5, ymax=2, c='black', linestyle='--', alpha=0.5)

    plt.title('Trajectory FFP', fontsize=40)
    plt.axes().set_aspect('equal', 'datalim')
    plt.xlabel('$x$ (AU)', fontsize=40)
    plt.ylabel('$y$ (AU)', fontsize=40)
    plt.legend(fontsize=40)
    plt.savefig('./trajectory_'+fout+'.png')

    plt.close()

def plot_orbital_elements(times_array, e_array, a_array, g_array, i_array, fout):

    #BP
    figure = plt.figure(figsize=(25,15))
    N_subplots = 4

    plot_e = figure.add_subplot(N_subplots,1,1)
    plot_i = figure.add_subplot(N_subplots,1,2)
    plot_g = figure.add_subplot(N_subplots,1,3)
    plot_a = figure.add_subplot(N_subplots,1,4)

    for i,ti in enumerate(times_array[:]):
        times_array_yr = times_array[i]
        plot_e.scatter(times_array_yr,e_array[i,0],s=1, lw=0)
        plot_i.scatter(times_array_yr,i_array[i,0]*180.0/np.pi,s=1, lw=0)
        plot_g.scatter(times_array_yr,g_array[i,0],s=1, lw=0)
        plot_a.scatter(times_array_yr,math.log10(a_array[i,0]),s=1, lw=0)

    t_max_yr = max(times_array)
    plot_e.set_xlim(0,t_max_yr)
    plot_i.set_xlim(0,t_max_yr)
    plot_g.set_xlim(0,t_max_yr)
    plot_a.set_xlim(0,t_max_yr)

    plot_e.set_xlabel('$t/\mathrm{yr}$')
    plot_i.set_xlabel('$t/\mathrm{yr}$')
    plot_g.set_xlabel('$t/\mathrm{yr}$')
    plot_a.set_xlabel('$t/\mathrm{yr}$')

    plot_e.set_ylabel('$e_\mathrm{BP}$')
    plot_i.set_ylabel('$i_\mathrm{BP} ({}^\circ)$')
    plot_g.set_ylabel('$g_\mathrm{BP} ({}^\circ)$')
    plot_a.set_ylabel('$\log{(a_\mathrm{BP})} (\mathrm{AU})$')
    figure.subplots_adjust(left=0.2, right=0.85, top=0.8, bottom=0.15)

    plt.savefig('./orbital_elements_bp_'+fout+'.png')
    plt.close()

    #FFP
    figure = plt.figure(figsize=(25,15))
    N_subplots = 4

    plot_e = figure.add_subplot(N_subplots,1,1)
    plot_i = figure.add_subplot(N_subplots,1,2)
    plot_g = figure.add_subplot(N_subplots,1,3)
    plot_a = figure.add_subplot(N_subplots,1,4)

    for i,ti in enumerate(times_array[:]):
        times_array_yr = times_array[i]
        plot_e.scatter(times_array_yr,e_array[i,1],s=1, lw=0)
        plot_i.scatter(times_array_yr,i_array[i,1]*180.0/np.pi,s=1, lw=0)
        plot_g.scatter(times_array_yr,g_array[i,1],s=1, lw=0)
        plot_a.scatter(times_array_yr,math.log10(a_array[i,1]),s=1, lw=0)

    t_max_yr = max(times_array)
    plot_e.set_xlim(0,t_max_yr)
    plot_i.set_xlim(0,t_max_yr)
    plot_g.set_xlim(0,t_max_yr)
    plot_a.set_xlim(0,t_max_yr)

    plot_e.set_xlabel('$t/\mathrm{yr}$')
    plot_i.set_xlabel('$t/\mathrm{yr}$')
    plot_g.set_xlabel('$t/\mathrm{yr}$')
    plot_a.set_xlabel('$t/\mathrm{yr}$')

    plot_e.set_ylabel('$e_\mathrm{FFP}$')
    plot_i.set_ylabel('$i_\mathrm{FFP} ({}^\circ)$')
    plot_g.set_ylabel('$g_\mathrm{FFP} ({}^\circ)$')
    plot_a.set_ylabel('$\log{(a_\mathrm{BP})} (\mathrm{AU})$')
    figure.subplots_adjust(left=0.2, right=0.85, top=0.8, bottom=0.15)

    plt.savefig('./orbital_elements_ffp_'+fout+'.png')

def run(m_bp=0.1, e_star_bp=0.0, e_star_ffp=0.998, sma_star_bp=1.0, sma_star_ffp=1540.5, om_star_bp=0.0, om_star_ffp=0.0, incl=0.0, dt_snapshots=0.000065, endtime=0.065, fout='kozai_out.dat'):

    ### set parameters
    m_star = 0.58|units.MSun
    m_bp = m_bp|units.MJupiter
    m_ffp = 7.5|units.MJupiter
    a_star_bp = sma_star_bp|units.AU
    a_star_ffp = sma_star_ffp|units.AU

    ### initial conditions for binary
    binaries = make_triple_system(m_star,m_bp,m_ffp,a_star_bp,a_star_ffp,e_star_bp,e_star_ffp,om_star_bp,om_star_ffp,incl)

    ### solve equations of motions for the system in the octuple approximation
    output_time_step = dt_snapshots | units.yr
    end_time = endtime | units.yr

    times_array, e_array, a_array, g_array, i_array, x, y, z = evolve_triple_system(binaries, end_time, output_time_step, fout)

    # plot_trajectory(x, y, fout)
    plot_orbital_elements(times_array, e_array, a_array, g_array, i_array, fout)

if __name__ in ('__main__','__plot__'):

    # #This one starts down, is similar to next one
    m_bp = 1.188889
    e_star_bp = 0.002831325835564035
    e_star_ffp = 0.9982296497127007
    sma_star_bp = 1.0000772686971517
    sma_star_ffp = 2318.9515580915295
    om_star_bp = -13.1687421576
    om_star_ffp = -12.3821332867
    incl = 0.0
    f = 'm1.188889e+00_a1.000000e+00_b1.123283e+01_p2.294588e+02'

    run(m_bp, e_star_bp, e_star_ffp, sma_star_bp, sma_star_ffp, om_star_bp, om_star_ffp, incl, dt_snapshots=6.5, endtime=65000.0, fout=f)

    #This one is at -300 and seem to be going further
    m_bp = 1.188889
    e_star_ffp = 0.9979511755073169
    e_star_bp = 0.004490745565236455
    sma_star_ffp = 1411.8796326270449
    sma_star_bp = 0.9992253527307187
    om_star_ffp = -8.07219911366
    om_star_bp = -47.4605952228
    # om_star_ffp = 0.0
    # om_star_bp = -8.07219911366+47.4605952228
    incl = 180.0
    f = 'm1.188889e+00_a1.000000e+00_b-8.023452e+00_p2.752418e+02'

    run(m_bp, e_star_bp, e_star_ffp, sma_star_bp, sma_star_ffp, om_star_bp, om_star_ffp, incl, dt_snapshots=6.5, endtime=65000.0, fout=f)
