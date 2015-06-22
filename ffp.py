import math, numpy, argparse, os

from matplotlib.pyplot import *

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

def get_planets(m0, a_planets, m_planets, e_planets, phi=0.):

    #Particle set with all planets
    planets = Particles()

    #Initialize star-planet orbit
    for i,a_i in enumerate(a_planets):
        #Binary
        star_planet = new_binary_from_orbital_elements(m0, m_planets[i], a_planets[i], e_planets[i], true_anomaly=phi)
        #Planets attributes
        star_planet.eccentricity = e_planets[i]
        star_planet.semimajoraxis = a_planets[i]
        #Center on the star
        star_planet.position -= star_planet[0].position
        star_planet.velocity -= star_planet[0].velocity
        planets.add_particle(star_planet[1])

    return planets

def get_parabolic_velocity(m0, m_ffp, b, r_inf, m_planets, a_planets, phi):

    M = m0 + m_ffp + m_planets.sum()

    cm_x = (np.array(m_planets)*np.array(a_planets*math.cos(np.radians(phi)))).sum()/M
    cm_y = (np.array(m_planets)*np.array(a_planets*math.sin(np.radians(phi)))).sum()/M

    r = ((r_inf-cm_x)**2+(b-cm_y)**2).sqrt()

    parabolic_velocity_squared = 2.0*(1.0 | nbody_system.length**3 * nbody_system.time**-2 * nbody_system.mass**-1)*M / r

    return parabolic_velocity_squared.sqrt()

def get_ffp_in_orbit(m0, m_ffp, b, r_inf, parabolic_velocity):

    m0_and_ffp_in_orbit = Particles(2)

    zero_p = 0.0 | nbody_system.length
    zero_v = 0.0 | nbody_system.speed

    #Central star
    m0_and_ffp_in_orbit[0].mass = m0
    m0_and_ffp_in_orbit[0].position = (zero_p,zero_p,zero_p)
    m0_and_ffp_in_orbit[0].velocity = (zero_v,zero_v,zero_v)

    #Free-floating planet
    m0_and_ffp_in_orbit[1].mass = m_ffp
    m0_and_ffp_in_orbit[1].position = (-r_inf,-b,zero_p)
    m0_and_ffp_in_orbit[1].velocity = (parabolic_velocity,zero_v,zero_v)

    m1, m2, sma, e, ta, i, lan, ap = orbital_elements_from_binary(m0_and_ffp_in_orbit)

    #For the star it sets the initial values of semimajoraxis and eccentricity of the ffp
    m0_and_ffp_in_orbit.eccentricity = e
    m0_and_ffp_in_orbit.semimajoraxis = sma

    return m0_and_ffp_in_orbit

def energies_binaries(bodies, indexA, indexB):
    """
    function to calculate energy of a binary (particle set with two particles)
    """

    labels = ['Star','FFP', 'BP1', 'BP2', 'BP3', 'BP4'] #... temporary

    particleA, particleB = bodies[indexA], bodies[indexB]

    m_A, m_B = particleA.mass, particleB.mass

    v_A = particleA.velocity.value_in(nbody_system.speed)
    vsquared_A = sum(v_A*v_A) | nbody_system.speed*nbody_system.speed
    v_B = particleB.velocity.value_in(nbody_system.speed)
    vsquared_B = sum(v_B*v_B) | nbody_system.speed*nbody_system.speed

    kinetic_energy = (m_A*vsquared_A + m_B*vsquared_B)/2.0

    #distances = (bodies.distances_squared(bodies[0])).sqrt()
    r_AB = (particleA.position - particleB.position).value_in(nbody_system.length)
    rmag_AB = math.sqrt(sum(r_AB*r_AB)) | nbody_system.length

    potential_energy = -m_A*m_B/rmag_AB*(1|nbody_system.length**3 * nbody_system.time**(-2) / nbody_system.mass)

    binary_energy = kinetic_energy+potential_energy

    return binary_energy

def save_particles_to_file(bodies, bodies_to_save, bodies_filename,time,converter):
    #Add attributes that I'm interested in to the bodies_to_save
    bodies_to_save.position = converter.to_si(bodies.position).as_quantity_in(units.AU)
    bodies_to_save.velocity = converter.to_si(bodies.velocity).as_quantity_in(units.kms)
    bodies_to_save.semimajoraxis = converter.to_si(bodies.semimajoraxis).as_quantity_in(units.AU)
    bodies_to_save.eccentricity = bodies.eccentricity
    bodies_to_save.time = converter.to_si(time).as_quantity_in(units.yr)

    write_set_to_file(bodies_to_save, bodies_filename, "hdf5")

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

        #Order: 0_ffp, 0_bp1, 0_bp2, 0_bp3,...
        for j in range(1,number_of_planets+2):
            binary = [bodies[0], bodies[j]]
            m1, m2, sma, e, ta, i, lan, ap = orbital_elements_from_binary(binary)
            bodies[j].eccentricity = e
            bodies[j].semimajoraxis = sma

        E = gravity.kinetic_energy + gravity.potential_energy
        system_energies.append(E)

        DeltaE = abs(E-E_initial)
        if ( DeltaE > DeltaE_max ):
            DeltaE_max = DeltaE

        save_particles_to_file(bodies, bodies_to_save, bodies_filename,time,converter)

        time += dt_snapshots

    max_energy_change = DeltaE_max/E_initial

    gravity.stop()

    return x,y,times,system_energies,max_energy_change

def plot_trajectory(x,y,number_of_planets):

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

    savefig('trajectory.png')

    xlim(-3.1,1.1)
    ylim(-1,1)

    savefig('trajectory_zoom.png')
    close()

def plot_energy(times, energies):

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

def convert_units(converter, t_end_p, m0_p, m_ffp_p, m_planets_p, a_planets_p, e_planets_p, n_steps_p, phi_p, b_p):

    #time of integration in yr
    t_end = converter.to_nbody(t_end_p | units.yr)
    #mass of the disk-central star in MSun
    m0 = converter.to_nbody(m0_p | units.MSun)
    #mass of the FFP MJupiter
    m_ffp = converter.to_nbody(m_ffp_p | units.MJupiter)
    #list of the masses of the bounded planets in MJupiter
    m_planets = converter.to_nbody(m_planets_p | units.MJupiter)
    #list of the semimajor axes of the bounded planets in AU
    a_planets = converter.to_nbody(a_planets_p | units.AU)
    #list of the eccentricities of the bounded planets
    e_planets = e_planets_p
    #number of steps for the integrator
    n_steps = n_steps_p
    #initial angle for of the bounded planet in degrees
    phi = phi_p
    #impact parameter of the FFP in AU
    b = converter.to_nbody(b_p | units.AU)

    return t_end, m0, m_ffp, m_planets, a_planets, e_planets, n_steps, phi, b

def stop_code():
    print '\nSTOP'
    import sys
    sys.exit()

def run_capture(t_end_p=650.0, m0_p=1.0, m_ffp_p=1.0, m_planets_p=[1.0], a_planets_p=[5.0], e_planets_p=[0.0], n_steps_p=10000, phi_p=0.0, b_p=5.0, iteration_number=1, n_snapshots_p=350):
    """
    Units: t_end(yr), m0(MSun), m_ffp(MJupiter), m_planets(MJupiter), a_planets(AU), e_planets(None), n_steps(None), phi(degrees), b(AU)
    """

    #Converter used in this program
    converter = nbody_system.nbody_to_si(1 | units.MSun,  5 | units.AU)

    #Conversion of units
    t_end, m0, m_ffp, m_planets, a_planets, e_planets, n_steps, phi, b = convert_units(converter, t_end_p, m0_p, m_ffp_p, m_planets_p, a_planets_p, e_planets_p, n_steps_p, phi_p, b_p)

    #Number of planets
    number_of_planets = len(e_planets)
    #Initial distance to the planet (x-cordinate)
    r_inf = 40.*max(a_planets)

    #Initialize planets
    planets = get_planets(m0,a_planets,m_planets,e_planets,phi)

    #Set the parabolic orbit of the ffp around the star
    star_and_ffp_in_orbit = get_ffp_in_orbit(m0, m_ffp, b, r_inf, get_parabolic_velocity(m0, m_ffp, b, r_inf, m_planets, a_planets, phi))

    #Particle superset: star, FFP, planets
    bodies = ParticlesSuperset([star_and_ffp_in_orbit, planets])

    #File to store results
    bodies_filename = './particles/'+str(iteration_number)+'.hdf5'

    #Evolve time
    x,y,times,energies,max_energy_change = evolve_gravity(bodies,number_of_planets,converter,t_end,n_steps,n_snapshots_p,bodies_filename)

    plot_trajectory(x,y,number_of_planets)
    plot_energy(times, energies)

    return max_energy_change

if __name__ in ('__main__', '__plot__'):

    run_capture()
