import math, numpy, argparse, os

from matplotlib.pyplot import *
from sympy.solvers import solve
from sympy import Symbol

from amuse.io import write_set_to_file
from amuse.io import base
from amuse.units import units, constants, nbody_system
from amuse.units.quantities import AdaptingVectorQuantity
from amuse.datamodel import Particles, ParticlesSuperset
from amuse.community.smalln.interface import SmallN
from amuse.ext.orbital_elements import orbital_elements_from_binary, new_binary_from_orbital_elements

def initialize_code(bodies, code=SmallN, timestep_parameter=0.0169):
    """
    initialize gravity code for bodies
    """

    stars_gravity = code()
    stars_gravity.particles.add_particles(bodies)
    stars_gravity.commit_particles()
    stars_gravity.parameters.timestep_parameter = timestep_parameter  # default 0.0169 time

    return stars_gravity

def get_parabolic_velocity(m0, m_ffp, b_ffp, r_inf, m_bp, a_bp, phi_bp):

    cm_M = m0 + m_bp
    cm_x = m_bp*a_bp*math.cos(np.radians(phi_bp))/cm_M
    cm_y = m_bp*a_bp*math.sin(np.radians(phi_bp))/cm_M

    M = m0 + m_ffp + m_bp
    r = ((r_inf-cm_x)**2+(b_ffp-cm_y)**2).sqrt()

    parabolic_velocity_squared = 2.0*M*(1.0 | nbody_system.length**3 * nbody_system.time**-2 * nbody_system.mass**-1)/r

    return parabolic_velocity_squared.sqrt()

def get_bodies_in_orbit(m0, m_ffp, m_bp, a_bp, e_bp, phi_bp, b_ffp, r_inf):

    #Bodies
    bodies = Particles()

    ##Get BP in orbit
    #Binary
    star_planet = new_binary_from_orbital_elements(m0, m_bp, a_bp, e_bp, true_anomaly=phi_bp)
    #Planet attributes
    star_planet.eccentricity = e_bp
    star_planet.semimajoraxis = a_bp
    #Center on the star
    star_planet.position -= star_planet[0].position
    star_planet.velocity -= star_planet[0].velocity
    cm_p = star_planet.center_of_mass()
    cm_v = star_planet.center_of_mass_velocity()

    ##Get FFP in orbit
    #Particle set
    m0_ffp = Particles(2)
    #Zeros and parabolic velocity
    zero_p = 0.0 | nbody_system.length
    zero_v = 0.0 | nbody_system.speed
    parabolic_velocity = get_parabolic_velocity(m0, m_ffp, b_ffp, r_inf, m_bp, a_bp, phi_bp)
    #Central star
    m0_ffp[0].mass = m0
    m0_ffp[0].position = (zero_p,zero_p,zero_p)
    m0_ffp[0].velocity = (zero_v,zero_v,zero_v)
    #Free-floating planet
    m0_ffp[1].mass = m_ffp
    m0_ffp[1].position = (-r_inf-cm_p[0],-b_ffp-cm_p[1],zero_p)
    m0_ffp[1].velocity = (parabolic_velocity,zero_v,zero_v)
    #Orbital Elements
    star_planet_as_one = Particles(1)
    star_planet_as_one.mass = m0 + m_bp
    star_planet_as_one.position = cm_p
    star_planet_as_one.velocity = cm_v
    binary = [star_planet_as_one[0], m0_ffp[1]]
    m1, m2, sma, e, ta, i, lan, ap = orbital_elements_from_binary(binary)
    #For the star it sets the initial values of semimajoraxis and eccentricity of the ffp around star+bp
    m0_ffp.eccentricity = e
    m0_ffp.semimajoraxis = sma

    #Order: star, ffp, bp
    bodies.add_particle(m0_ffp[0])
    bodies.add_particle(m0_ffp[1])
    bodies.add_particle(star_planet[1])

    return bodies

def save_particles_to_file(bodies, bodies_to_save, bodies_filename,time,converter):
    #Add attributes that I'm interested in to the bodies_to_save
    bodies_to_save.position = converter.to_si(bodies.position).as_quantity_in(units.AU)
    bodies_to_save.velocity = converter.to_si(bodies.velocity).as_quantity_in(units.kms)
    bodies_to_save.semimajoraxis = converter.to_si(bodies.semimajoraxis).as_quantity_in(units.AU)
    bodies_to_save.eccentricity = bodies.eccentricity
    bodies_to_save.time = converter.to_si(time).as_quantity_in(units.yr)

    write_set_to_file(bodies_to_save, bodies_filename, "hdf5")

def solve_for_x(m0, m_bp, m_ffp):

    coefficients = [m_bp+m_ffp, 2*m_bp+3*m_ffp, m_bp+3*m_ffp, -(m_bp+3*m0), -(3*m0+2*m_bp), -(m0+m_bp)]
    solutions = numpy.roots(coefficients)

    return abs(solutions[math.ceil(len(solutions)/2.0)-1])

def is_hill_stable(m_values, a_values, e_values, converter):

    # m0 = converter.to_si(m_values[0]).value_in(units.MJupiter)
    # m_ffp = converter.to_si(m_values[1]).value_in(units.MJupiter)
    # m_bp = converter.to_si(m_values[2]).value_in(units.MJupiter)

    m0 = m_values[0].value_in(nbody_system.mass)
    m_ffp = m_values[1].value_in(nbody_system.mass)
    m_bp = m_values[2].value_in(nbody_system.mass)

    M = m0 + m_bp + m_ffp
    mu = m0 + m_bp

    # a_1 = converter.to_si(a_values[2]).value_in(units.AU)
    # a_2 = converter.to_si(a_values[1]).value_in(units.AU)

    a_1 = a_values[2].value_in(nbody_system.length)
    a_2 =a_values[1].value_in(nbody_system.length)

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

def evolve_gravity(bodies, number_of_planets, converter, t_end, n_steps, n_snapshots, bodies_filename):

    #Particles that will be saved in the file
    bodies_to_save = Particles(number_of_planets+2)

    #Positions and velocities centered on the center of mass
    bodies.move_to_center()

    time = 0. | nbody_system.time
    dt = t_end / float(n_steps)
    dt_snapshots = t_end / float(n_snapshots)

    gravity = initialize_code(bodies, timestep_parameter = dt.value_in(nbody_system.time))
    channel_from_gr_to_framework = gravity.particles.new_channel_to(bodies)

    x = AdaptingVectorQuantity()
    y = AdaptingVectorQuantity()
    times = AdaptingVectorQuantity()
    system_energies = AdaptingVectorQuantity()

    E_initial = gravity.kinetic_energy + gravity.potential_energy
    DeltaE_max = 0.0 | nbody_system.energy

    while time<=t_end:

        gravity.evolve_model(time)
        channel_from_gr_to_framework.copy()

        bodies.collection_attributes.timestamp = time

        gravity.particles.collection_attributes.timestamp = time

        x.append(bodies.x)
        y.append(bodies.y)
        times.append(time)

        E = gravity.kinetic_energy + gravity.potential_energy
        system_energies.append(E)

        DeltaE = abs(E-E_initial)
        if ( DeltaE > DeltaE_max ):
            DeltaE_max = DeltaE

        #Orbital Elements
        star = bodies[0]
        ffp = bodies[1]
        bp = bodies[2]

        #Star+BP and FFP
        cm_x = (star.mass*star.x + bp.mass*bp.x)/(star.mass+bp.mass)
        cm_y = (star.mass*star.y + bp.mass*bp.y)/(star.mass+bp.mass)
        cm_vx = (star.mass*star.vx + bp.mass*bp.vx)/(star.mass+bp.mass)
        cm_vy = (star.mass*star.vy + bp.mass*bp.vy)/(star.mass+bp.mass)

        star_bp = Particles(1)
        star_bp.mass = star.mass + bp.mass
        star_bp.position = [cm_x, cm_y, 0.0 | nbody_system.length]
        star_bp.velocity = [cm_vx, cm_vy, 0.0 | nbody_system.speed]

        binary = [star_bp[0], ffp]
        m1, m2, sma_starbp_ffp, e_starbp_ffp, ta, i, lan, ap = orbital_elements_from_binary(binary)
        bodies[1].eccentricity = e_starbp_ffp
        bodies[1].semimajoraxis = sma_starbp_ffp

        #Star and BP
        binary = [star, bp]
        m1, m2, sma_star_bp, e_star_bp, ta, i, lan, ap = orbital_elements_from_binary(binary)
        bodies[2].eccentricity = e_star_bp
        bodies[2].semimajoraxis = sma_star_bp

        save_particles_to_file(bodies, bodies_to_save, bodies_filename,time,converter)

        time += dt_snapshots

    max_energy_change = DeltaE_max/E_initial

    is_stable = is_hill_stable(bodies.mass, bodies.semimajoraxis, bodies.eccentricity, converter)

    gravity.stop()

    return x,y,times,system_energies,max_energy_change,is_stable

def plot_trajectory(x,y,number_of_planets, j):

    colors = ['magenta', 'green', 'DarkOrange', 'red']

    f=figure(figsize=(35,15))

    x_star = x[:,0].value_in(nbody_system.length)
    x_ffp = x[:,1].value_in(nbody_system.length)

    y_star = y[:,0].value_in(nbody_system.length)
    y_ffp = y[:,1].value_in(nbody_system.length)

    plot(x_star,y_star,'y',label='Star')
    scatter(x_star[0],y_star[0],c='black',marker='*')
    scatter(x_star[-1],y_star[-1],c='y',marker='*')

    plot(x_ffp,y_ffp,'c',label='FFP')
    scatter(x_ffp[0],y_ffp[0],c='black')
    scatter(x_ffp[-1],y_ffp[-1],c='c')

    for i in range(0, number_of_planets):

        x_planet = x[:,i+2].value_in(nbody_system.length)
        y_planet = y[:,i+2].value_in(nbody_system.length)

        color_planet = colors[i]

        plot(x_planet,y_planet,color=color_planet,label='BP',alpha=0.5)
        scatter(x_planet[0],y_planet[0],c='black')
        scatter(x_planet[-1],y_planet[-1],color=color_planet)

    axhline(y=0, xmin=-80, xmax=10, c='black', linestyle='--')
    axvline(x=0, ymin=-5, ymax=2, c='black', linestyle='--')

    title('Trajectory FFP (nbody units)')
    xlabel("$x$", fontsize=20)
    ylabel("$y$", fontsize=20)
    legend()

    savefig('trajectory'+str(j)+'.png')

    xlim(-3.1,1.1)
    ylim(-1,1)

    close()

def plot_energy_change(times, energies):

    times = times.value_in(nbody_system.time)
    energies = energies.value_in(nbody_system.energy)

    initial_energy = energies[0]
    energies = (energies-initial_energy)/initial_energy

    f=figure(figsize=(15,15))

    plot(times,numpy.log10(abs(energies)),color='black')

    axhline(y=0, xmin=0, xmax=times[-1], c='m', linestyle='--')

    title('Total Energy of the System (nbody units)')
    xlabel("$t$", fontsize=20)
    ylabel("$\log{(|\Delta E / E|)}$", fontsize=20)

    savefig('energy.png')
    close()

def convert_units(converter, t_end_p, m0_p, m_ffp_p, e_bp_p, m_bp_p, a_bp_p, b_ffp_p, phi_bp_p):

    #time of integration in yr
    t_end = converter.to_nbody(t_end_p | units.yr)
    #mass of the disk-central star in MSun
    m0 = converter.to_nbody(m0_p | units.MSun)
    #mass of the FFP MJupiter
    m_ffp = converter.to_nbody(m_ffp_p | units.MJupiter)
    #eccentricity of the bounded planets
    e_bp = e_bp_p
    #mass of the bounded planets in MJupiter
    m_bp = converter.to_nbody(m_bp_p | units.MJupiter)
    #semimajor axis of the bounded planets in AU
    a_bp = converter.to_nbody(a_bp_p | units.AU)
    #impact parameter of the FFP in AU
    b_ffp = converter.to_nbody(b_ffp_p | units.AU)
    #initial angle for of the bounded planet in degrees
    phi_bp = phi_bp_p

    return t_end, m0, m_ffp, e_bp, m_bp, a_bp, b_ffp, phi_bp

def stop_code():
    print '\nSTOP'
    import sys
    sys.exit()

def run_capture(t_end_p=650.0, m0_p=0.58, m_ffp_p=7.5, e_bp_p=0.0, m_bp_p=0.1, a_bp_p=1.0, b_ffp_p=1.0, phi_bp_p=0.0, n_steps=10000, path_filename='./particles/m/m_a_b_p_.hdf5', n_snapshots=500):
    """
    Units: t_end_p(yr), m0_p(MSun), m_ffp_p(MJupiter), e_bp_p(None), m_bp_p(MJupiter), a_bp_p(AU), b_ffp(AU), phi_p(degrees)
    """

    #Converter used in this program
    converter = nbody_system.nbody_to_si(1 | units.MSun,  5 | units.AU)

    #Conversion of units
    t_end, m0, m_ffp, e_bp, m_bp, a_bp, b_ffp, phi_bp = convert_units(converter, t_end_p, m0_p, m_ffp_p, e_bp_p, m_bp_p, a_bp_p, b_ffp_p, phi_bp_p)

    #Number of planets
    number_of_planets = 1
    #Initial distance to the planet (x-cordinate)
    r_inf = 40.0*a_bp

    #Particle superset: star, FFP, planets
    bodies = get_bodies_in_orbit(m0, m_ffp, m_bp, a_bp, e_bp, phi_bp, b_ffp, r_inf)

    #Evolve time
    x,y,times,energies,max_energy_change,is_stable = evolve_gravity(bodies, number_of_planets, converter, t_end, n_steps, n_snapshots, path_filename)

    #plot_trajectory(x,y,number_of_planets,phi_bp)
    #plot_energy_change(times, energies)

    return max_energy_change, is_stable

if __name__ in ('__main__', '__plot__'):

    run_capture()
