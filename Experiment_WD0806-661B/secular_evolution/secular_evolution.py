from amuse.community.seculartriple.interface import SecularTriple
from amuse.datamodel import Particles
from amuse.units import units

import numpy as np
import os, time, sys, signal

def stop_code():
    print '\nSTOP'
    import sys
    sys.exit()
    
# Register an handler for the timeout
# From http://stackoverflow.com/questions/492519/timeout-on-a-python-function-call
def handler(signum, frame):
    raise Exception("timeout")

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

def evolve_secular_per_mass(mass_dir):
    
    m_star = 0.58
    m_ffp = 7.5
    
    filename_stables = './particles/'+mass_dir+'/stables_'+mass_dir+'.txt'
    filename_secular_stables = './particles/'+mass_dir+'/secular_stables_'+mass_dir+'.txt'
    filename_status = './particles/'+mass_dir+'/secular_status.txt'
    filename_stuck = './particles/'+mass_dir+'/stuck.txt'
    
    fs=open(filename_stables)    
    num_stables = sum(1 for line in fs)
    fs.close()
    
    #Counting the lines
    try:
        ss = open(filename_secular_stables,'r')
        num_runs = sum(1 for line in ss)
        ss.close()
    except:
        num_runs = 0
    
    stables_file = open(filename_stables,'r')    
    
    cont = 1    
    
    line = stables_file.readline()
    
    while(line != ''):
        
        #Check that you don't run the last parameters combination again, it starts on the last one (in case the program has to be restarted)        
        if(cont > num_runs+1):

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
    
            endtime = integration_time*100.0
            new_integration_time = endtime + integration_time
            
            to_write = ('%d/%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f')%(cont,num_stables,endtime,m_bp,e_star_bp,e_star_ffp,a_star_bp,a_star_ffp,i_star_bp,i_star_ffp,lan_star_bp,lan_star_ffp,ap_star_bp,ap_star_ffp)
            command = 'echo "'+to_write+'" >> '+filename_status
            os.system(command)
            
            #Start Handler
            #Register the signal function handler
            signal.signal(signal.SIGALRM, handler)
            #Define a timeout for your function
            signal.alarm(10)
            
            try:        
                #Secular Evolve
                e_sb, e_sf, a_sb, a_sf, i_sb, i_sf, lan_sb, lan_sf, ap_sb, ap_sf = run(endtime, m_star, m_bp, m_ffp, e_star_bp, e_star_ffp, a_star_bp, a_star_ffp, i_star_bp, i_star_ffp, lan_star_bp, lan_star_ffp, ap_star_bp, ap_star_ffp)
        
                #Running time        
                run_time = run_time_before + time.time() - start_time
                
                #Save final elements in file after secular evolution        
                to_write = ('%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%.4f\t%.4e\t%f')%(m_bp,a_bp,b_ffp,phi_bp,inc_bp,lan_bp,e_sf,e_sb,a_sf,a_sb,i_sf,i_sb,lan_sf,lan_sb,ap_sf,ap_sb,run_time,energy_change,new_integration_time)
                command = 'echo "'+to_write+'" >> '+filename_secular_stables
                os.system(command)
            
            #Raised when the program gets stucked -> inestable                            
            except Exception, exc: 
                #Save final elements in file after secular evolution        
                to_write = ('%d/%d\t%f\t%f\t%f\t%f\t%f\t%f')%(cont,num_stables,m_bp,a_bp,b_ffp,phi_bp,inc_bp,lan_bp)
                command = 'echo "'+to_write+'" >> '+filename_stuck
                os.system(command)
            
            #Cancel the timer if the function returned before timeout
            signal.alarm(0)

        cont += 1

        line = stables_file.readline()

    stables_file.close()


if __name__ in ('__main__','__plot__'):

    #Obtain mass directory
    m_dir = sys.argv[1]
    
    #Evolve
    evolve_secular_per_mass(m_dir)


