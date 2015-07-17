import numpy as np
import os, time, math, sys
import ffp

def find_limit_b(r_0, number_r0_in_rinf, mass_ratio):
    return r_0*math.sqrt(( ((1+mass_ratio)**6)/(4.0*(mass_ratio**4)) + ((1+mass_ratio)**3)*math.sqrt( ((1+mass_ratio)**6)/(16.0*(mass_ratio**4)) + number_r0_in_rinf**2 )/(mass_ratio**2) )/2.0)

def permute(m_bp):

    #Fixed Values
    m0 = 0.58 #MSun
    m_ffp = 7.5 #MJupiter
    e_bp = 0.0

    t_end = 650.0 #yr
    n_steps = 12000
    n_snapshots = int(n_steps/20.0)
    n_r0_in_rinf = 45.0
    b_limit_fraction = 0.3 #number from 0 to 1 that defines how much of the upper and lower limits do I want to take

    #Numbers of each parameter
    n_as_bp = 10
    n_bs_ffp = 10
    n_phis_bp = 100

    #Variable parameters
    as_bp = np.linspace(1.0,50.0,n_as_bp) #AU

    m_bp_filename = ('m%.4e')%(m_bp)

    #Counting the lines
    try:
        par = open('./particles/'+m_bp_filename+'/parameters_'+m_bp_filename+'.txt','r')
        num_runs = sum(1 for line in par)
        par.close()
        file_parameters = open('./particles/'+m_bp_filename+'/parameters_'+m_bp_filename+'.txt','a')
    except:
        num_runs = 0
        os.makedirs('./particles/'+m_bp_filename)
        file_parameters = open('./particles/'+m_bp_filename+'/parameters_'+m_bp_filename+'.txt','w')

    i=1

    #File for the stables
    file_stables = open('./particles/'+m_bp_filename+'/stables_'+m_bp_filename+'.txt','a')

    mass_ratio = m_ffp/m_bp

    for a_bp in as_bp:

        b_limit = find_limit_b(a_bp, n_r0_in_rinf, mass_ratio)*b_limit_fraction
        bs_ffp = np.linspace(-b_limit, b_limit, n_bs_ffp)

        for b_ffp in bs_ffp:

            np.random.seed(12)

            for j in range(n_phis_bp):

                phi_bp = np.random.random()*360.0 #degrees (this should be 1000 phis)

                #Check that you don't run the last parameters combination again, it starts on the last one (in case the program has to be restarted)
                if (i>num_runs):
                    #Time starts
                    start_time = time.time()

                    #Filename
                    max_energy_change, is_stable, e_star_ffp, e_star_bp, sma_star_ffp, sma_star_bp, inc_star_ffp, inc_star_bp, lan_star_ffp, lan_star_bp, ap_star_ffp, ap_star_bp = ffp.run_capture(t_end_p=t_end,
                                                                                                                                                                                                    m0_p=m0,
                                                                                                                                                                                                    m_ffp_p=m_ffp,
                                                                                                                                                                                                    e_bp_p=e_bp,
                                                                                                                                                                                                    m_bp_p=m_bp,
                                                                                                                                                                                                    a_bp_p=a_bp,
                                                                                                                                                                                                    b_ffp_p=b_ffp,
                                                                                                                                                                                                    phi_bp_p=phi_bp,
                                                                                                                                                                                                    n_steps=n_steps,
                                                                                                                                                                                                    n_snapshots=n_snapshots,
                                                                                                                                                                                                    n_r0_in_rinf=n_r0_in_rinf)

                    #Time stops
                    duration = time.time()-start_time

                    #Write in File
                    file_parameters.write(('%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%.4f\t%.4e\n')%(m_bp,a_bp,b_ffp,phi_bp,e_star_ffp,e_star_bp,sma_star_ffp,sma_star_bp,inc_star_ffp,inc_star_bp,lan_star_ffp,lan_star_bp,ap_star_ffp,ap_star_bp,duration,max_energy_change))

                    if(is_stable):
                        file_stables.write(('%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%.4f\t%.4e\n')%(m_bp,a_bp,b_ffp,phi_bp,e_star_ffp,e_star_bp,sma_star_ffp,sma_star_bp,inc_star_ffp,inc_star_bp,lan_star_ffp,lan_star_bp,ap_star_ffp,ap_star_bp,duration,max_energy_change))

                #Advance counter
                i += 1

    file_parameters.close()
    file_stables.close()

if __name__ in ('__main__', '__plot__'):

    #Obtain mass
    m_bp = float(sys.argv[1])

    #Permute all of its a, b, phi
    permute(m_bp)
