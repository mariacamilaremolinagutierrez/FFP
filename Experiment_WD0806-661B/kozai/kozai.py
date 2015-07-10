from amuse.community.seculartriple.interface import SecularTriple
from amuse.datamodel import Particles
from amuse.units import units,constants,quantities
from amuse.units.optparse import OptionParser
from amuse.ext.orbital_elements import new_binary_from_orbital_elements

from matplotlib import pyplot

import numpy

def make_triple_system(m1=0.58|units.MSun, m2=7.5|units.MJupiter, m3=0.1|units.MJupiter, a1=1540.5|units.AU, a2=1.0|units.AU, e1=0.998, e2=0.0, inc=0.0, om1=0.0, om2=0.0):

    binaries = Particles(2)
    inner_stars = Particles(2)
    outer_stars = Particles(2)

    binaries[0].child1 = inner_stars[0]
    binaries[0].child2 = inner_stars[1]
    binaries[1].child1 = outer_stars[0]
    binaries[1].child2 = outer_stars[1]

    binaries[0].child1.mass = m1
    binaries[0].child2.mass = m2
    binaries[1].child1.mass = m3

    binaries[0].semimajor_axis = a1
    binaries[1].semimajor_axis = a2

    binaries[0].eccentricity = e1
    binaries[1].eccentricity = e2
    binaries[0].argument_of_pericenter = om1*numpy.pi/180.0
    binaries[1].argument_of_pericenter = om2*numpy.pi/180.0
    binaries[0].longitude_of_ascending_node = 0.0
    binaries[1].longitude_of_ascending_node = 0.0

    binaries[0].inclination = inc*numpy.pi/180.0 ### this is the inclination between the orbital planes of binaries[0] and binaries[1] ###
    binaries[1].inclination = 0.0

    # binaries[0].child1.envelope_mass = 0.0 | units.MSun
    # binaries[0].child2.envelope_mass = 0.0 | units.MSun
    # binaries[1].child1.envelope_mass = 0.0 | units.MSun
    # binaries[0].child1.radius = 1.0 | units.RSun
    # binaries[0].child2.radius = 1.0 | units.RSun
    # binaries[1].child1.radius = 1.0 | units.RSun
    # binaries[0].child1.luminosity = 1.0 | units.LSun
    # binaries[0].child2.luminosity = 1.0 | units.LSun
    # binaries[1].child1.luminosity = 1.0 | units.LSun
    # binaries[1].child1.envelope_radius = 0.0 | units.RSun
    # binaries[0].child1.envelope_radius = 0.0 | units.RSun
    # binaries[0].child2.envelope_radius = 0.0 | units.RSun
    # binaries[0].child1.apsidal_motion_constant = 0.1
    # binaries[0].child2.apsidal_motion_constant = 0.1
    # binaries[1].child1.apsidal_motion_constant = 0.1
    # binaries[0].child1.gyration_radius = 0.1
    # binaries[0].child2.gyration_radius = 0.1
    # binaries[1].child1.gyration_radius = 0.1
    # binaries[0].child1.stellar_type = 1
    # binaries[0].child2.stellar_type = 1
    # binaries[1].child1.stellar_type = 1
    # binaries[0].child1.spin_angular_frequency = 1.0e3 | 1.0/units.yr
    # binaries[0].child2.spin_angular_frequency = 1.0e3 | 1.0/units.yr
    # binaries[1].child1.spin_angular_frequency = 1.0e3 | 1.0/units.yr

    # binaries[0].child1.mass_transfer_rate = -1.0e-6 | units.MSun / units.yr
    # binaries[0].mass_transfer_accretion_parameter = 1.0
    # binaries[1].mass_transfer_accretion_parameter = 1.0
    # binaries[0].mass_transfer_angular_momentum_loss_parameter = 0.0
    # binaries[1].mass_transfer_angular_momentum_loss_parameter = 0.0

    return binaries

def give_position(semimajor_axis, eccentricity, inclination, argument_of_periapsis):

    inclination = numpy.radians(inclination)
    argument_of_periapsis = numpy.radians(argument_of_periapsis)

    #Supposing this:
    #longitude_of_the_ascending_node = 0.0
    #true_anomaly = 0.0

    cos_inclination = numpy.cos(inclination)
    sin_inclination = numpy.sin(inclination)

    cos_arg_per = numpy.cos(argument_of_periapsis)
    sin_arg_per = numpy.sin(argument_of_periapsis)

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
    code.parameters.f_mass_transfer = 0.0
    code.parameters.f_1PN_in = 0.0
    code.parameters.f_1PN_out = 0.0
    code.parameters.f_25PN_in = 0.0
    code.parameters.f_25PN_out = 0.0

    ### quantities later used for plotting ###
    times_array = quantities.AdaptingVectorQuantity()
    x = quantities.AdaptingVectorQuantity()
    y = quantities.AdaptingVectorQuantity()
    z = quantities.AdaptingVectorQuantity()
    a_in_array = quantities.AdaptingVectorQuantity()
    e_in_array = []
    i_tot_array = []
    g_in_array = []
    a1_in_array = quantities.AdaptingVectorQuantity()
    e1_in_array = []
    i1_tot_array = []
    g1_in_array = []

    mc = 1

    print code.binaries

    time = 0.0 | units.Myr
    while (time < end_time):

        code.evolve_model(time)
        time += output_time_step

        inc = code.binaries[0].inclination
        inc1 = code.binaries[1].inclination
        ecc = code.binaries[0].eccentricity
        ecc1 = code.binaries[1].eccentricity
        sma = code.binaries[0].semimajor_axis
        sma1 = code.binaries[1].semimajor_axis
        ap = code.binaries[0].argument_of_pericenter
        ap1 = code.binaries[1].argument_of_pericenter

        times_array.append(time)
        e_in_array.append(ecc)
        a_in_array.append(sma)
        g_in_array.append(ap)
        i_tot_array.append(inc)
        e1_in_array.append(ecc1)
        a1_in_array.append(sma1)
        g1_in_array.append(ap1)
        i1_tot_array.append(inc1)

        print code.binaries

        xx,yy,zz = give_position(sma, ecc, inc, ap)
        xx1,yy1,zz1 = give_position(sma1, ecc1, inc1, ap1)

        xs = [0.0, xx.value_in(units.AU), xx1.value_in(units.AU)] |units.AU
        ys = [0.0, yy.value_in(units.AU), yy1.value_in(units.AU)] |units.AU
        zs = [0.0, zz.value_in(units.AU), zz1.value_in(units.AU)] |units.AU

        # print xs.value_in(units.AU)
        # print ys.value_in(units.AU)
        # print zs.value_in(units.AU)

        x.append(xs)
        y.append(ys)
        z.append(zs)

        if (mc == 3):
            code.stop()
        mc += 1

    e_in_array = numpy.array(e_in_array)
    g_in_array = numpy.array(g_in_array)
    i_tot_array = numpy.array(i_tot_array)

    return times_array,a_in_array,e_in_array,g_in_array,i_tot_array,x,y,z

def plot_trajectory(x, y):

    f = pyplot.figure(figsize=(70,30))

    x_star = (x[:,0]-x[:,0]).value_in(units.AU)
    x_ffp = (x[:,1]-x[:,0]).value_in(units.AU)
    x_planet = (x[:,2]-x[:,0]).value_in(units.AU)

    y_star = (y[:,0]-y[:,0]).value_in(units.AU)
    y_ffp = (y[:,1]-y[:,0]).value_in(units.AU)
    y_planet = (y[:,2]-y[:,0]).value_in(units.AU)

    pyplot.plot(x_star,y_star,c='y',label='Star')
    pyplot.scatter(x_star[0],y_star[0],c='black',marker='*', lw = 0)
    pyplot.scatter(x_star[-1],y_star[-1],c='y',marker='*', lw = 0)

    pyplot.plot(x_ffp,y_ffp,c='c',label='FFP', lw = 2)
    pyplot.scatter(x_ffp[0],y_ffp[0],c='black', lw = 0)
    pyplot.scatter(x_ffp[-1],y_ffp[-1],c='c', lw = 0)

    pyplot.plot(x_planet,y_planet,c='m',label='BP',alpha=0.5)
    pyplot.scatter(x_planet[0],y_planet[0],c='black', lw = 0)
    pyplot.scatter(x_planet[-1],y_planet[-1],c='m', lw = 0)

    pyplot.axhline(y=0, xmin=-80, xmax=10, c='black', linestyle='--')
    pyplot.axvline(x=0, ymin=-5, ymax=2, c='black', linestyle='--')

    pyplot.title('Trajectory FFP', fontsize=40)
    pyplot.axes().set_aspect('equal', 'datalim')
    pyplot.xlabel('$x$ (AU)', fontsize=40)
    pyplot.ylabel('$y$ (AU)', fontsize=40)
    pyplot.legend(fontsize=40)
    pyplot.savefig('./kozai_trajectory.png')

    pyplot.close()

def get_omega_stat(omega):
  """
  statistics of the time evolution of omega
  """
  o_min = min(omega)
  o_max = max(omega)
  o_med = numpy.median(omega)
  o_mean = numpy.mean(omega)
  return o_min, o_max, o_med, o_mean

def run(m_bp=0.1, e_star_ffp=0.998, e_star_bp=0.0, sma_star_ffp=1540.5, sma_star_bp=1.0, incl=0.0, om1=0.0, om2=0.0, dt_snapshots=0.000065, endtime=0.065, fout='kozai_out.dat',mx=5.0):

    ### set parameters
    m1=0.58|units.MSun
    m2=7.5|units.MJupiter
    m3=m_bp|units.MJupiter
    e1=e_star_ffp
    e2=e_star_bp
    a1=sma_star_ffp|units.AU
    a2=sma_star_bp|units.AU

    ### initial conditions for binary
    binaries = make_triple_system(m1,m2,m3,a1,a2,e1,e2,incl,om1,om2)

    ### solve equations of motions for the system in the octuple approximation
    output_time_step = dt_snapshots | units.Myr
    end_time = endtime | units.Myr

    times_array_o,a_in_array_o,e_in_array_o,g_in_array_o,i_tot_array_o,x,y,z = evolve_triple_system(binaries, end_time, output_time_step, fout)

    plot_trajectory(x, y)


if __name__ in ('__main__','__plot__'):

    #This one starts down, is similar to next one
    m1 = 1.188889
    e_star_ffp1 = 0.9982296497127007
    e_star_bp1 = 0.002831325835564035
    sma_star_ffp1 = 2318.9515580915295
    sma_star_bp1 = 1.0000772686971517
    inc1 = 0.0
    om11 = -12.3821332867
    om21 = -13.1687421576
    f1 = 'm1.188889e+00_a1.000000e+00_b1.123283e+01_p2.294588e+02.dat'

    #This one is at -300 and seem to be going further
    m2 = 1.188889
    e_star_ffp2 = 0.9979511755073169
    e_star_bp2 = 0.004490745565236455
    sma_star_ffp2 = 1411.8796326270449
    sma_star_bp2 = 0.9992253527307187
    inc2 = 180.0
    om12 = -8.07219911366
    om22 = -47.4605952228
    f2 = 'm1.188889e+00_a1.000000e+00_b-8.023452e+00_p2.752418e+02.dat'

    run(m1, e_star_ffp1, e_star_bp1, sma_star_ffp1, sma_star_bp1, inc1, om11, om21, dt_snapshots=0.000065, endtime=0.065, fout=f1,mx=5.0)
    #run(m2, e_star_ffp2, e_star_bp2, sma_star_ffp2, sma_star_bp2, inc2, om12, om22, dt_snapshots=0.000065, endtime=0.0065, fout=f2,mx=5.0)
