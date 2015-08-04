from amuse.community.seculartriple.interface import SecularTriple
from amuse.datamodel import Particles
from amuse.units import units

import numpy as np
import os, time, sys

def stop_code():
    print '\nSTOP'
    import sys
    sys.exit()

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

def save_fractional_changes(e_star_bp, e_star_ffp, a_star_bp, a_star_ffp, i_star_bp, i_star_ffp, lan_star_bp, lan_star_ffp, ap_star_bp, ap_star_ffp, e_sb, e_sf, a_sb, a_sf, i_sb, i_sf, lan_sb, lan_sf, ap_sb, ap_sf, filename_fractional_changes):
    #first = initial, second = final
    if (e_star_ffp != 0.0 and e_sf != 0.0):
        e_star_ffp_change = (e_sf-e_star_ffp)/e_star_ffp
        e_sfs = ''
    else:
        e_star_ffp_change = (e_sf-e_star_ffp)
        e_sfs = 'D'

    if (e_star_bp != 0.0 and e_sb != 0.0):
        e_star_bp_change = (e_sb-e_star_bp)/e_star_bp
        e_sbs = ''
    else:
        e_star_bp_change = (e_sb-e_star_bp)
        e_sbs = 'D'

    a_star_ffp_change = (a_sf-a_star_ffp)/a_star_ffp
    a_star_bp_change = (a_sb-a_star_bp)/a_star_bp

    if (i_star_ffp != 0.0 and i_sf != 0.0):
        i_star_ffp_change = (i_sf-i_star_ffp)/i_star_ffp
        i_sfs = ''
    else:
        i_star_ffp_change = (i_sf-i_star_ffp)
        i_sfs = 'D'

    if (i_star_bp != 0.0 and i_sb != 0.0):
        i_star_bp_change = (i_sb-i_star_bp)/i_star_bp
        i_sbs = ''
    else:
        i_star_bp_change = (i_sb-i_star_bp)
        i_sbs = 'D'

    if (lan_star_ffp != 0.0 and lan_sf != 0.0):
        lan_star_ffp_change = (lan_sf-lan_star_ffp)/lan_star_ffp
        lan_sfs = ''
    else:
        lan_star_ffp_change = (lan_sf-lan_star_ffp)
        lan_sfs = 'D'

    if (lan_star_bp != 0.0 and lan_sb != 0.0):
        lan_star_bp_change = (lan_sb-lan_star_bp)/lan_star_bp
        lan_sbs = ''
    else:
        lan_star_bp_change = (lan_sb-lan_star_bp)
        lan_sbs = 'D'

    if (ap_star_ffp != 0.0 and ap_sf != 0.0):
        ap_star_ffp_change = (ap_sf-ap_star_ffp)/ap_star_ffp
        ap_sfs = ''
    else:
        ap_star_ffp_change = (ap_sf-ap_star_ffp)
        ap_sfs = 'D'

    if (ap_star_bp != 0.0 and ap_sb != 0.0):
        ap_star_bp_change = (ap_sb-ap_star_bp)/ap_star_bp
        ap_sbs = ''
    else:
        ap_star_bp_change = (ap_sb-ap_star_bp)
        ap_sbs = 'D'

    to_write = (e_sfs+'%f\t'+e_sbs+'%f\t%f\t%f\t'+i_sfs+'%f\t'+i_sbs+'%f\t'+lan_sfs+'%f\t'+lan_sbs+'%f\t'+ap_sfs+'%f\t'+ap_sbs+'%f')%(e_star_ffp_change,e_star_bp_change,a_star_ffp_change,a_star_bp_change,i_star_ffp_change,i_star_bp_change,lan_star_ffp_change,lan_star_bp_change,ap_star_ffp_change,ap_star_bp_change)
    command = 'echo "'+to_write+'" >> '+filename_fractional_changes
    os.system(command)

    
def evolve_secular_per_mass(mass_dir):
    
    m_star = 0.58
    m_ffp = 7.5
    
    filename_stables = './particles/'+mass_dir+'/stables_'+mass_dir+'.txt'
    filename_secular_stables = './particles/'+mass_dir+'/secular_stables_'+mass_dir+'.txt'
    filename_fractional_changes = './particles/'+mass_dir+'/fractional_changes_'+mass_dir+'.txt'
    filename_status = './particles/'+mass_dir+'/secular_status.txt'
    
    fs=open(filename_stables)    
    num_stables = sum(1 for line in fs)
    fs.close()
    
    stables_file = open(filename_stables,'r')    
    
    cont = 1    
    
    line = stables_file.readline()
    
    while(line != ''):

        start_time = time.time()

        parts = line.split('\t')

        m_bp = float(parts[0])
        a_bp = float(parts[1])
        b_ffp = float(parts[2])
        phi_bp = float(parts[3])
        inc_bp = float(parts[4])
        lan_bp = float(parts[5])
        
        e_star_ffp = float(parts[6])
        e_star_bp = float(parts[7])
        a_star_ffp = float(parts[8])
        a_star_bp = float(parts[9])
        i_star_ffp = float(parts[10])
        i_star_bp = float(parts[11])
        lan_star_ffp = float(parts[12])
        lan_star_bp = float(parts[13])
        ap_star_ffp = float(parts[14])
        ap_star_bp = float(parts[15])

        run_time_before = float(parts[16])
        energy_change = float(parts[17])
        integration_time = float(parts[18].split('\t')[0])

        endtime = integration_time*20.0
        new_integration_time = endtime + integration_time
        
        to_write = ('%d/%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f')%(cont,num_stables,endtime,m_bp,e_star_bp,e_star_ffp,a_star_bp,a_star_ffp,i_star_bp,i_star_ffp,lan_star_bp,lan_star_ffp,ap_star_bp,ap_star_ffp)
        command = 'echo "'+to_write+'" >> '+filename_status
        os.system(command)
        
        #Secular Evolve
        e_sb, e_sf, a_sb, a_sf, i_sb, i_sf, lan_sb, lan_sf, ap_sb, ap_sf = run(endtime, m_star, m_bp, m_ffp, e_star_bp, e_star_ffp, a_star_bp, a_star_ffp, i_star_bp, i_star_ffp, lan_star_bp, lan_star_ffp, ap_star_bp, ap_star_ffp)

        #Running time        
        run_time = run_time_before + time.time() - start_time
        
        #Save final elements in file after secular evolution        
        to_write = ('%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%.4f\t%.4e\t%f')%(m_bp,a_bp,b_ffp,phi_bp,inc_bp,lan_bp,e_sf,e_sb,a_sf,a_sb,i_sf,i_sb,lan_sf,lan_sb,ap_sf,ap_sb,run_time,energy_change,new_integration_time)
        command = 'echo "'+to_write+'" >> '+filename_secular_stables
        os.system(command)
        
        #Save fractional changes in orbital elements
        save_fractional_changes(e_star_bp, e_star_ffp, a_star_bp, a_star_ffp, i_star_bp, i_star_ffp, lan_star_bp, lan_star_ffp, ap_star_bp, ap_star_ffp, e_sb, e_sf, a_sb, a_sf, i_sb, i_sf, lan_sb, lan_sf, ap_sb, ap_sf, filename_fractional_changes)

        cont += 1

        line = stables_file.readline()

    stables_file.close()


if __name__ in ('__main__','__plot__'):

    #Obtain mass directory
    m_dir = sys.argv[1]
    
    #Evolve
    evolve_secular_per_mass(m_dir)


