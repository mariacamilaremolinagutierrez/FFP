import numpy as np
import os, time, math, sys
import ffp

def stop_code():
    print '\nSTOP'
    import sys
    sys.exit()
    
def get_critical_velocity(m0, m_ffp, m_bp, r_inf):

    MJupiter_in_kg = 1.89813e27 #kg
    MSun_in_kg = 1.988435e30 #kg
    AU_in_m = 1.496e11 #m  

    G = 6.674e-11 #m**3/(kg*s**2)
    M = m0*MSun_in_kg+m_bp*MJupiter_in_kg #kg
    r = r_inf*AU_in_m #m

    critical_velocity = math.sqrt(2.0*M*G/r)

    return critical_velocity/1000.0 #km/s

def find_limit_b(r_0, number_r0_in_rinf, mass_ratio):
    return r_0*math.sqrt(( ((1+mass_ratio)**6)/(4.0*(mass_ratio**4)) + ((1+mass_ratio)**3)*math.sqrt( ((1+mass_ratio)**6)/(16.0*(mass_ratio**4)) + number_r0_in_rinf**2 )/(mass_ratio**2) )/2.0)

def squared_space(b_initial, b_final, n_bs):

    if (b_initial<0.0):
        
        #assuming b_initial = b_final
        ss = np.zeros(n_bs)
        ss_squared = np.linspace(-b_initial**2,b_final**2,n_bs)
        for i in range(len(ss_squared)):
            sssi = ss_squared[i]            
            if (sssi >= 0.0):
                ss[i] = math.sqrt(sssi)
            else:
                ss[i] = -math.sqrt(-sssi)
    else:
        ss = np.sqrt(np.linspace(b_initial**2,b_final**2,n_bs))
    return ss

def permute(m_ffp):
    
    #Fixed Values
    m0 = 1.0 #MSun
    m_bp = 1.0 #MJupiter
    e_bp = 0.0
    a_bp = 5.20336301 #AU

    n_snapshots = 600
    n_r0_in_rinf = 40.0
    b_limit_fraction = 1.0 #number from 0 to 1 that defines how much of the upper and lower limits do I want to take

    #Numbers of each parameter
    n_vs_ffp = 40
    n_bs_ffp = 4000
    #this should last 3.4 days
    
    total_permutations = n_vs_ffp*n_bs_ffp
    
    #Setting random seed    
    np.random.seed(12)

    #Variable parameters
    critical_velocity = get_critical_velocity(m0, m_ffp, m_bp, n_r0_in_rinf*a_bp) #km/s
    vs_ffp_log = 10**(np.linspace(np.log10(0.5*critical_velocity),np.log10(10.0*critical_velocity),n_vs_ffp)) #km/s

    #Filenames to write info
    m_ffp_filename = ('m%.4e')%(m_ffp)
    filename_parameters = './particles/'+m_ffp_filename+'/parameters_'+m_ffp_filename+'.txt'
    filename_stables = './particles/'+m_ffp_filename+'/stables_'+m_ffp_filename+'.txt'
    filename_status = './particles/'+m_ffp_filename+'/status_'+m_ffp_filename+'.txt'

    #Counting the lines
    try:
        par = open(filename_parameters,'r')
        num_runs = sum(1 for line in par)
        par.close()
    except:
        num_runs = 0
        os.makedirs('./particles/'+m_ffp_filename)

    i=1

    mass_ratio = m_ffp/m_bp

    for v_ffp in vs_ffp_log:

        b_limit = find_limit_b(a_bp, n_r0_in_rinf, mass_ratio)*b_limit_fraction
        bs_ffp = squared_space(0.0, b_limit, n_bs_ffp)

        for b_ffp in bs_ffp:
            
            phi_bp = np.random.random()*360.0            
            inc_bp = np.rad2deg(np.arccos(np.random.random()))
            lan_bp = np.random.random()*360.0
            ap_bp = np.random.random()*360.0

            #Check that you don't run the last parameters combination again, it starts on the last one (in case the program has to be restarted)
            if (i>num_runs):
                
                #status
                to_write = str(i)+'/'+str(int(total_permutations))+'\t'+str(v_ffp)+'\t'+str(m_ffp)+'\t'+str(b_ffp)+'\t'+str(phi_bp)+'\t'+str(inc_bp)+'\t'+str(lan_bp)+'\t'+str(ap_bp)
                res = os.system('echo "'+to_write+'" >> '+filename_status)

                #Time starts
                start_time = time.time()

                #Filename
                t_end, max_energy_change, is_stable, e_star_ffp, e_star_bp, sma_star_ffp, sma_star_bp, inc_star_ffp, inc_star_bp, lan_star_ffp, lan_star_bp, ap_star_ffp, ap_star_bp = ffp.run_capture(m0_p=m0,
                                                                                                                                                                                                       m_bp_p=m_bp,
                                                                                                                                                                                                       e_bp_p=e_bp,                                                
                                                                                                                                                                                                       a_bp_p=a_bp,
                                                                                                                                                                                                       m_ffp_p=m_ffp,
                                                                                                                                                                                                       v_ffp_p=v_ffp,
                                                                                                                                                                                                       b_ffp_p=b_ffp,
                                                                                                                                                                                                       phi_bp_p=phi_bp,
                                                                                                                                                                                                       inc_bp_p=inc_bp,
                                                                                                                                                                                                       lan_bp_p=lan_bp,
                                                                                                                                                                                                       ap_bp_p=ap_bp,  
                                                                                                                                                                                                       n_snapshots=n_snapshots,
                                                                                                                                                                                                       n_r0_in_rinf=n_r0_in_rinf)
                
                #Time stops
                duration = time.time()-start_time
                
                #Write in files
                line = ('%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%.4f\t%.4e\t%f')%(m_ffp,v_ffp,b_ffp,phi_bp,inc_bp,lan_bp,ap_bp,e_star_ffp,e_star_bp,sma_star_ffp,sma_star_bp,inc_star_ffp,inc_star_bp,lan_star_ffp,lan_star_bp,ap_star_ffp,ap_star_bp,duration,max_energy_change,t_end)
                command = 'echo "'+line+'" >> '+filename_parameters
                os.system(command)
                
                if(is_stable):
                    command = 'echo "'+line+'" >> '+filename_stables
                    os.system(command)

            #Advance counter
            i += 1

if __name__ in ('__main__', '__plot__'):

    #Obtain mass
    m_ffp = float(sys.argv[1])

    #Permute all of its a, b, phi
    permute(m_ffp)
