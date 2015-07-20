from amuse.community.seculartriple.interface import SecularTriple
from amuse.datamodel import Particles
from amuse.units import units

import numpy as np
import os, time

def stop_code():
    print '\nSTOP'
    import sys
    sys.exit()

def create_parameters_and_stables_files():

    masses_directories = os.listdir('./particles/')

    with open('./parameters.txt', 'w') as outfile:
        outfile.write('m_bp\ta_bp\tb_ffp\tphi_bp\te_star_ffp\te_star_bp\tsma_star_ffp\tsma_star_bp\tinc_star_ffp\tinc_star_bp\tlan_star_ffp\tlan_star_bp\tap_star_ffp\tap_star_bp\trun_time\tenergy_change\n')
        for i in range(len(masses_directories)):
            mass_dir = masses_directories[i]
            fname = './particles/'+mass_dir+'/parameters_'+mass_dir+'.txt'
            with open(fname) as infile:
                for line in infile:
                    outfile.write(line)
            infile.close()
    outfile.close()

    with open('./stables.txt', 'w') as outfile:
        outfile.write('m_bp\ta_bp\tb_ffp\tphi_bp\te_star_ffp\te_star_bp\tsma_star_ffp\tsma_star_bp\tinc_star_ffp\tinc_star_bp\tlan_star_ffp\tlan_star_bp\tap_star_ffp\tap_star_bp\trun_time\tenergy_change\n')
        for i in range(len(masses_directories)):
            mass_dir = masses_directories[i]
            fname = './particles/'+mass_dir+'/stables_'+mass_dir+'.txt'
            with open(fname) as infile:
                for line in infile:
                    outfile.write(line)
            infile.close()
    outfile.close()

def make_triple_system(m_star, m_bp, m_ffp, e_star_bp, e_star_ffp, a_star_bp, a_star_ffp, i_star_bp, i_star_ffp, lan_star_bp, lan_star_ffp, ap_star_bp, ap_star_ffp):

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

    binaries[0].eccentricity = e_star_bp
    binaries[1].eccentricity = e_star_ffp
    binaries[0].semimajor_axis = a_star_bp
    binaries[1].semimajor_axis = a_star_ffp
    binaries[0].inclination = 0.0
    binaries[1].inclination = i_star_ffp - i_star_bp #inclination between the orbital planes of binaries[0] and binaries[1]
    binaries[0].longitude_of_ascending_node = lan_star_bp
    binaries[1].longitude_of_ascending_node = lan_star_ffp
    binaries[0].argument_of_pericenter = ap_star_bp
    binaries[1].argument_of_pericenter = ap_star_ffp

    return binaries

def evolve_triple_system(binaries, end_time):

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

    # print '\nBINARIES AT START:'
    # print code.binaries

    code.evolve_model(end_time)

    ecc_star_bp = code.binaries[0].eccentricity
    ecc_star_ffp = code.binaries[1].eccentricity
    sma_star_bp = code.binaries[0].semimajor_axis
    sma_star_ffp = code.binaries[1].semimajor_axis
    inc_star_bp = code.binaries[0].inclination
    inc_star_ffp = code.binaries[1].inclination
    lan_star_bp = code.binaries[0].longitude_of_ascending_node
    lan_star_ffp = code.binaries[1].longitude_of_ascending_node
    ap_star_bp = code.binaries[0].argument_of_pericenter
    ap_star_ffp = code.binaries[1].argument_of_pericenter

    # print '\nBINARIES AT END:'
    # print code.binaries

    return ecc_star_bp, ecc_star_ffp, sma_star_bp, sma_star_ffp, inc_star_bp, inc_star_ffp, lan_star_bp, lan_star_ffp, ap_star_bp, ap_star_ffp

def run(endtime, m_star, m_bp, m_ffp, e_star_bp, e_star_ffp, a_star_bp, a_star_ffp, i_star_bp, i_star_ffp, lan_star_bp, lan_star_ffp, ap_star_bp, ap_star_ffp):

    #setting units
    end_time = endtime | units.yr
    m_star = m_star | units.MSun
    m_bp = m_bp | units.MJupiter
    m_ffp = m_ffp | units.MJupiter
    a_star_bp = a_star_bp | units.AU
    a_star_ffp = a_star_ffp | units.AU
    i_star_bp = np.deg2rad(i_star_bp)
    i_star_ffp = np.deg2rad(i_star_ffp)
    lan_star_bp = np.deg2rad(lan_star_bp)
    lan_star_ffp = np.deg2rad(lan_star_ffp)
    ap_star_bp = np.deg2rad(ap_star_bp)
    ap_star_ffp = np.deg2rad(ap_star_ffp)

    #initial conditions for binary
    binaries = make_triple_system(m_star, m_bp, m_ffp, e_star_bp, e_star_ffp, a_star_bp, a_star_ffp, i_star_bp, i_star_ffp, lan_star_bp, lan_star_ffp, ap_star_bp, ap_star_ffp)

    #solve equations of motions for the system in the octuple approximation
    ecc_star_bp, ecc_star_ffp, sma_star_bp, sma_star_ffp, inc_star_bp, inc_star_ffp, lan_star_bp, lan_star_ffp, ap_star_bp, ap_star_ffp = evolve_triple_system(binaries, end_time)

    return ecc_star_bp, ecc_star_ffp, sma_star_bp.value_in(units.AU), sma_star_ffp.value_in(units.AU), np.rad2deg(inc_star_bp), np.rad2deg(inc_star_ffp), np.rad2deg(lan_star_bp), np.rad2deg(lan_star_ffp), np.rad2deg(ap_star_bp), np.rad2deg(ap_star_ffp)

if __name__ in ('__main__','__plot__'):

    endtime = 10000000.0
    m_star = 0.58
    m_ffp = 7.5

    create_parameters_and_stables_files()

    # stop_code()

    secular_stables_file = open('./secular_stables.txt','w')
    stables_file = open('./stables.txt','r')

    #write the header
    line = stables_file.readline()
    secular_stables_file.write(line)

    cont = 1

    line = stables_file.readline()
    while(line != ''):

        start_time = time.time()

        parts = line.split('\t')

        m_bp = float(parts[0])
        a_bp = float(parts[1])
        b_ffp = float(parts[2])
        phi_bp = float(parts[3])

        e_star_ffp = float(parts[4])
        e_star_bp = float(parts[5])
        a_star_ffp = float(parts[6])
        a_star_bp = float(parts[7])
        i_star_ffp = float(parts[8])
        i_star_bp = float(parts[9])
        lan_star_ffp = float(parts[10])
        lan_star_bp = float(parts[11])
        ap_star_ffp = float(parts[12])
        ap_star_bp = float(parts[13])

        energy_change = float(parts[15].split('\t')[0])

        if (e_star_bp < 1.0 and e_star_ffp < 1.0):

            print cont

            e_sb, e_sf, a_sb, a_sf, i_sb, i_sf, lan_sb, lan_sf, ap_sb, ap_sf = run(endtime, m_star, m_bp, m_ffp, e_star_bp, e_star_ffp, a_star_bp, a_star_ffp, i_star_bp, i_star_ffp, lan_star_bp, lan_star_ffp, ap_star_bp, ap_star_ffp)

            run_time = float(parts[14]) + time.time() - start_time

            secular_stables_file.write(('%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%.4f\t%.4e\n')%(m_bp,a_bp,b_ffp,phi_bp,e_sf,e_sb,a_sf,a_sb,i_sf,i_sb,lan_sf,lan_sb,ap_sf,ap_sb,run_time,energy_change))

            if(cont==10):
                break

            cont += 1

        line = stables_file.readline()

    secular_stables_file.close()
    stables_file.close()
