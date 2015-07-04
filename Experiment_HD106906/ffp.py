import math, numpy, argparse, os

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
    cm_x = m_bp*a_bp*math.cos(numpy.radians(phi_bp))/cm_M
    cm_y = m_bp*a_bp*math.sin(numpy.radians(phi_bp))/cm_M

    M = m0 + m_ffp + m_bp
    r = ((r_inf-cm_x)**2+(b_ffp-cm_y)**2).sqrt()

    parabolic_velocity_squared = 2.0*M*(1.0 | nbody_system.length**3 * nbody_system.time**-2 * nbody_system.mass**-1)/r

    return parabolic_velocity_squared.sqrt()

def get_bodies_in_orbit(m0, m_ffp, m_bp, a_bp, e_bp, phi_bp, lan_bp, b_ffp, r_inf):

    #Bodies
    bodies = Particles()

    ##Get BP in orbit
    #Binary
    star_planet = new_binary_from_orbital_elements(m0, m_bp, a_bp, e_bp, true_anomaly=phi_bp, inclination = 28, longitude_of_the_ascending_node = lan_bp)
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

def solve_for_x(m0, m_bp, m_ffp):

    coefficients = [m_bp+m_ffp, 2*m_bp+3*m_ffp, m_bp+3*m_ffp, -(m_bp+3*m0), -(3*m0+2*m_bp), -(m0+m_bp)]
    solutions = numpy.roots(coefficients)

    return abs(solutions[math.ceil(len(solutions)/2.0)-1])

def is_hill_stable(m_values, a_values, e_values, converter):

    m0 = m_values[0].value_in(nbody_system.mass)
    m_ffp = m_values[1].value_in(nbody_system.mass)
    m_bp = m_values[2].value_in(nbody_system.mass)

    M = m0 + m_bp + m_ffp
    mu = m0 + m_bp

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

def evolve_gravity(bodies, converter, t_end, n_steps, n_snapshots):

    #Positions and velocities centered on the center of mass
    bodies.move_to_center()

    time = 0. | nbody_system.time
    dt = t_end / float(n_steps)
    dt_snapshots = t_end / float(n_snapshots)

    gravity = initialize_code(bodies, timestep_parameter = dt.value_in(nbody_system.time))
    channel_from_gr_to_framework = gravity.particles.new_channel_to(bodies)

    E_initial = gravity.kinetic_energy + gravity.potential_energy
    DeltaE_max = 0.0 | nbody_system.energy

    orbital_elements = numpy.zeros((n_snapshots,4))

    snapshot = 0

    while time<=t_end:

        gravity.evolve_model(time)
        channel_from_gr_to_framework.copy()

        bodies.collection_attributes.timestamp = time
        gravity.particles.collection_attributes.timestamp = time

        E = gravity.kinetic_energy + gravity.potential_energy

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

        orbital_elements[snapshot,0] = e_starbp_ffp
        orbital_elements[snapshot,1] = e_star_bp
        orbital_elements[snapshot,2] = converter.to_si(sma_starbp_ffp).value_in(units.AU)
        orbital_elements[snapshot,3] = converter.to_si(sma_star_bp).value_in(units.AU)

        snapshot += 1
        time += dt_snapshots

    max_energy_change = DeltaE_max/E_initial

    is_stable = is_hill_stable(bodies.mass, bodies.semimajoraxis, bodies.eccentricity, converter)

    if(is_stable):
        numpy.save('./orbital_elements/orbital_elements',orbital_elements) #Change filenameeeeee

    gravity.stop()

    return max_energy_change,is_stable

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

def plot_trajectory(x,y,z):

    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D

    f = plt.figure(figsize=(15,15))

    ax = f.add_subplot(111, projection='3d')

    x_star = x[:,0].value_in(nbody_system.length)
    x_ffp = x[:,1].value_in(nbody_system.length)
    x_planet = x[:,2].value_in(nbody_system.length)

    y_star = y[:,0].value_in(nbody_system.length)
    y_ffp = y[:,1].value_in(nbody_system.length)
    y_planet = y[:,2].value_in(nbody_system.length)

    z_star = z[:,0].value_in(nbody_system.length)
    z_ffp = z[:,1].value_in(nbody_system.length)
    z_planet = z[:,2].value_in(nbody_system.length)

    ax.plot(x_star,y_star,z_star,'y',label='Star')
    ax.scatter(x_star[0],y_star[0],z_star[0],c='black',marker='*')
    ax.scatter(x_star[-1],y_star[-1],z_star[-1],c='y',marker='*')

    ax.plot(x_ffp,y_ffp,z_ffp,'c',label='FFP')
    ax.scatter(x_ffp[0],y_ffp[0],z_ffp[0],c='black')
    ax.scatter(x_ffp[-1],y_ffp[-1],z_ffp[-1],c='c')

    ax.plot(x_planet,y_planet,z_planet,color='magenta',label='BP',alpha=0.5)
    ax.scatter(x_planet[0],y_planet[0],z_planet[0],c='black')
    ax.scatter(x_planet[-1],y_planet[-1],z_planet[-1],color='magenta')

    ax.set_title('Trajectory FFP (nbody units)')
    ax.legend()
    ax.set_xlim(-40,4)
    ax.set_ylim(-22,22)
    ax.set_zlim(-22,22)

    plt.savefig('trajectory.png')
    plt.show()
    plt.close()

def run_capture(t_end_p=650.0, m0_p=1.5, m_ffp_p=11, e_bp_p=0.0, m_bp_p=0.1, a_bp_p=5.0, b_ffp_p=1.0, phi_bp_p=0.0, lan_bp=0.0, n_steps=12000, n_snapshots=600):
    """
    Units: t_end_p(yr), m0_p(MSun), m_ffp_p(MJupiter), e_bp_p(None), m_bp_p(MJupiter), a_bp_p(AU), b_ffp(AU), phi_p(degrees)
    """
    #Converter used in this program
    converter = nbody_system.nbody_to_si(1 | units.MSun,  5 | units.AU)

    #Conversion of units
    t_end, m0, m_ffp, e_bp, m_bp, a_bp, b_ffp, phi_bp = convert_units(converter, t_end_p, m0_p, m_ffp_p, e_bp_p, m_bp_p, a_bp_p, b_ffp_p, phi_bp_p)

    #Initial distance to the planet (x-cordinate)
    r_inf = 40.0*a_bp

    #Particle superset: star, FFP, planets
    bodies = get_bodies_in_orbit(m0, m_ffp, m_bp, a_bp, e_bp, phi_bp, lan_bp, b_ffp, r_inf)

    #Evolve time
    max_energy_change,is_stable = evolve_gravity(bodies, converter, t_end, n_steps, n_snapshots)

    return max_energy_change, is_stable

if __name__ in ('__main__', '__plot__'):

    print run_capture()
