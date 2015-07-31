import math
import numpy as np

from amuse.units import units, nbody_system
from amuse.datamodel import Particles
from amuse.community.smalln.interface import SmallN
from amuse.community.kepler.interface import Kepler
from amuse.ext.orbital_elements import new_binary_from_orbital_elements

def stop_code():
    print '\nSTOP'
    import sys
    sys.exit()

def initialize_code(bodies, code=SmallN, timestep_parameter=0.0169):
    """
    initialize gravity code for bodies
    """
    stars_gravity = code(channel_type="sockets")
    stars_gravity.particles.add_particles(bodies)
    stars_gravity.commit_particles()
    stars_gravity.parameters.timestep_parameter = timestep_parameter  # default 0.0169 time

    return stars_gravity

def get_parabolic_velocity(m0, m_bp, b_ffp, r_inf):

    M = m0 + m_bp
    r = (r_inf**2 + b_ffp**2).sqrt()

    parabolic_velocity_squared = 2.0*M*(1.0 | nbody_system.length**3 * nbody_system.time**-2 * nbody_system.mass**-1)/r

    return parabolic_velocity_squared.sqrt()

def my_orbital_elements_from_binary(binary, G=nbody_system.G):

    mass1=binary[0].mass
    mass2=binary[1].mass
    position = binary[1].position-binary[0].position
    velocity = binary[1].velocity-binary[0].velocity
    total_mass = mass1 + mass2

    specific_energy = (1.0/2.0)*velocity.lengths_squared() - G*total_mass/position.lengths()
    specific_angular_momentum = position.cross(velocity)
    specific_angular_momentum_norm = specific_angular_momentum.lengths()
    specific_angular_momentum_unit=specific_angular_momentum/specific_angular_momentum_norm

    ### Semimajor axis ###
    if (specific_energy == 0.0 | nbody_system.length**2 * nbody_system.time**-2):
        semimajor_axis = -float('Inf') | nbody_system.length
    else:
        semimajor_axis = -G*total_mass/(2.0*specific_energy)

    ### Eccentricity ###
    eccentricity_argument = 2.0*specific_angular_momentum_norm**2*specific_energy/(G**2*total_mass**2)
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
    mu = G*total_mass ### Argument of pericenter ###
    position_unit = position/position.lengths()
    e_vector = ( (1.0/mu)*velocity.cross(specific_angular_momentum) - position_unit ) | units.none
    if (e_vector.lengths() == 0.0): ### Argument of pericenter and true anomaly cannot be determined for e = 0, in this case return 1.0 for the cosines ###
        cos_arg_per = 1.0
        arg_per=0.
    else:
        e_vector_unit = e_vector/e_vector.lengths()

        cos_arg_per = np.dot(e_vector_unit,ascending_node_vector_unit)
        e_cross_an=np.cross(e_vector_unit,ascending_node_vector_unit)
        ss=-np.sign(np.dot(specific_angular_momentum_unit,e_cross_an))
        sin_arg_per = ss*(e_cross_an**2).sum()**0.5
        arg_per=np.degrees(np.arctan2(sin_arg_per,cos_arg_per))

    return semimajor_axis, eccentricity, inclination, long_asc_node, arg_per

def get_bodies_in_orbit(m0, m_ffp, m_bp, a_bp, e_bp, phi_bp, inc_bp, lan_bp, b_ffp, r_inf):

    #Bodies
    bodies = Particles()

    ##Get BP in orbit
    #Binary
    star_planet = new_binary_from_orbital_elements(m0, m_bp, a_bp, e_bp, true_anomaly=phi_bp, inclination = inc_bp, longitude_of_the_ascending_node = lan_bp)
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
    parabolic_velocity = get_parabolic_velocity(m0, m_bp, b_ffp, r_inf)
    #Central star
    m0_ffp[0].mass = m0
    m0_ffp[0].position = (zero_p,zero_p,zero_p)
    m0_ffp[0].velocity = (zero_v,zero_v,zero_v)
    #Free-floating planet
    m0_ffp[1].mass = m_ffp
    m0_ffp[1].position = (-r_inf+cm_p[0], b_ffp+cm_p[1], cm_p[2])
    m0_ffp[1].velocity = (parabolic_velocity+cm_v[0], cm_v[1], cm_v[2])

    #To find the orbital period of the BP
    G = (1.0 | nbody_system.length**3 * nbody_system.time**-2 * nbody_system.mass**-1)
    orbital_period_bp = 2*math.pi*((a_bp**3)/(G*m0)).sqrt()

    #To find the distance and time to periastron
    kep = Kepler()
    kep.initialize_code()

    star_planet_as_one = Particles(1)
    star_planet_as_one.mass = m0 + m_bp
    star_planet_as_one.position = cm_p
    star_planet_as_one.velocity = cm_v

    kepler_bodies = Particles()
    kepler_bodies.add_particle(star_planet_as_one[0])
    kepler_bodies.add_particle(m0_ffp[1])

    kep.initialize_from_particles(kepler_bodies)

    kep.advance_to_periastron()
    time_pericenter = kep.get_time()
    
#    print kepler_bodies
#    print kep.get_separation_vector()
#    print time_pericenter

    kep.stop()

    binary = [star_planet_as_one[0], m0_ffp[1]]
    sma, e, inclination, long_asc_node, arg_per = my_orbital_elements_from_binary(binary)
    m0_ffp.eccentricity = e
    m0_ffp.semimajoraxis = sma

    #Adding bodies. Order: star, ffp, bp
    bodies.add_particle(m0_ffp[0])
    bodies.add_particle(m0_ffp[1])
    bodies.add_particle(star_planet[1])

    return bodies, time_pericenter, orbital_period_bp

def solve_for_x(m0, m_bp, m_ffp):

    coefficients = [m_bp+m_ffp, 2*m_bp+3*m_ffp, m_bp+3*m_ffp, -(m_bp+3*m0), -(3*m0+2*m_bp), -(m0+m_bp)]
    solutions = np.roots(coefficients)

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

    if (e_1 > 1.0 or e_2 > 1.0 or a_1 < 0.0 or a_2 < 0.0):
        return False
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

def evolve_gravity(bodies, converter, t_end, dt_integration, n_snapshots):

    #Positions and velocities centered on the center of mass
    bodies.move_to_center()
    
    t = 0.0 | nbody_system.time
    dt_snapshots = t_end / float(n_snapshots)

    gravity = initialize_code(bodies, timestep_parameter = dt_integration.value_in(nbody_system.time))
    channel_from_gr_to_framework = gravity.particles.new_channel_to(bodies)

    E_initial = gravity.kinetic_energy + gravity.potential_energy
    DeltaE_max = 0.0 | nbody_system.energy
    
    if(t_end != 0.0 | nbody_system.time):
        
        while (t <= t_end):
            
            gravity.evolve_model(t)
            channel_from_gr_to_framework.copy()
    
            E = gravity.kinetic_energy + gravity.potential_energy
    
            DeltaE = abs(E-E_initial)
            if ( DeltaE > DeltaE_max ):
                DeltaE_max = DeltaE
    
            t += dt_snapshots
        
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
        sma_starbp_ffp, e_starbp_ffp, inc_starbp_ffp, lan_starbp_ffp, ap_starbp_ffp = my_orbital_elements_from_binary(binary)
        bodies[1].eccentricity = e_starbp_ffp
        bodies[1].semimajoraxis = sma_starbp_ffp
    
        #Star and BP
        binary = [star, bp]
        sma_star_bp, e_star_bp, inc_star_bp, lan_star_bp, ap_star_bp = my_orbital_elements_from_binary(binary)
        bodies[2].eccentricity = e_star_bp
        bodies[2].semimajoraxis = sma_star_bp
    
        is_stable = is_hill_stable(bodies.mass, bodies.semimajoraxis, bodies.eccentricity, converter)
    
        #Star and FFP
        binary = [star, ffp]
        sma_star_ffp, e_star_ffp, inc_star_ffp, lan_star_ffp, ap_star_ffp = my_orbital_elements_from_binary(binary)
    
        max_energy_change = DeltaE_max/E_initial
    
        gravity.stop()
    
        return max_energy_change, is_stable, e_star_ffp, e_star_bp, sma_star_ffp, sma_star_bp, inc_star_ffp, inc_star_bp, lan_star_ffp, lan_star_bp, ap_star_ffp, ap_star_bp
    
    else:
        return 1.0, False, 0.0, 0.0, 1.0 | nbody_system.length, 1.0 | nbody_system.length, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0

def convert_units(converter, m0_p, m_ffp_p, e_bp_p, m_bp_p, a_bp_p, b_ffp_p, phi_bp_p, inc_bp_p, lan_bp_p):

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
    #inclination in degrees
    inc_bp = inc_bp_p
    #longitude_of_the_ascending_node in degrees
    lan_bp = lan_bp_p

    return m0, m_ffp, e_bp, m_bp, a_bp, b_ffp, phi_bp, inc_bp, lan_bp

def run_capture(m0_p=0.58, m_ffp_p=7.5, e_bp_p=0.0, m_bp_p=0.1, a_bp_p=1.0, b_ffp_p=1.0, phi_bp_p=0.0, inc_bp_p=0.0, lan_bp_p=0.0, n_snapshots=600, n_r0_in_rinf=45.0):
    """
    Units: m0_p(MSun), m_ffp_p(MJupiter), e_bp_p(None), m_bp_p(MJupiter), a_bp_p(AU), b_ffp(AU), phi_p(degrees)
    """
    #Converter used in this program
    converter = nbody_system.nbody_to_si(1 | units.MSun,  5 | units.AU)

    #Conversion of units
    m0, m_ffp, e_bp, m_bp, a_bp, b_ffp, phi_bp, inc_bp, lan_bp = convert_units(converter, m0_p, m_ffp_p, e_bp_p, m_bp_p, a_bp_p, b_ffp_p, phi_bp_p, inc_bp_p, lan_bp_p)

    #Initial distance to the planet (x-cordinate)
    r_inf = n_r0_in_rinf*a_bp

    #Particles: star, FFP, planets
    bodies, time_pericenter, orbital_period_bp = get_bodies_in_orbit(m0, m_ffp, m_bp, a_bp, e_bp, phi_bp, inc_bp, lan_bp, b_ffp, r_inf)
    
    #Evolve time
    t_end = 5.0*time_pericenter
    dt_integration = orbital_period_bp/50.0
    max_energy_change, is_stable, e_star_ffp, e_star_bp, sma_star_ffp, sma_star_bp, inc_star_ffp, inc_star_bp, lan_star_ffp, lan_star_bp, ap_star_ffp, ap_star_bp = evolve_gravity(bodies, converter, t_end, dt_integration, n_snapshots)

    return converter.to_si(t_end).value_in(units.yr), max_energy_change, is_stable, e_star_ffp, e_star_bp, converter.to_si(sma_star_ffp).value_in(units.AU), converter.to_si(sma_star_bp).value_in(units.AU), inc_star_ffp, inc_star_bp, lan_star_ffp, lan_star_bp, ap_star_ffp, ap_star_bp

if __name__ in ('__main__', '__plot__'):
    
    m=0.1
    a=6.44444444444 
    b=-75.1528396747
    p=55.4986232567
    i=0.0	
    l=0.0
    
    print run_capture(m_bp_p=m, a_bp_p=a, b_ffp_p=b, phi_bp_p=p, inc_bp_p=i, lan_bp_p=l)
