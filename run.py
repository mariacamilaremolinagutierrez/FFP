import numpy as np
import os, time, math
import ffp

def squared_space(b_initial, b_final, n_bs):
    if (b_initial<0.0): #assuming b_initial = b_final
        ss_half = np.sqrt(np.linspace(0.0,b_final**2,int(n_bs/2)))
        if(n_bs%2 == 0):
            ss = np.append(-1.0*ss_half[::-1],ss_half)
        else:
            ss = np.append(-1.0*ss_half[::-1],zeros(1),ss_half)
    else:
        ss = np.sqrt(np.linspace(b_initial**2,b_final**2,n_bs))
    return ss

def find_limit_b(mass_ratio):
    return math.sqrt(-(-mass_ratio**6-6*mass_ratio**5-15*mass_ratio**4-20*mass_ratio**3-15*mass_ratio**2-(mass_ratio+1)**3*math.sqrt(mass_ratio**6+6*mass_ratio**5+25615*mass_ratio**4+20*mass_ratio**3+15*mass_ratio**2+6*mass_ratio+1)-6*mass_ratio-1)/(8*mass_ratio**4))

def create_info_file():
    info = open('./info.txt','w')
    info.write('PARTICLES:\n')
    info.write('key\teccentricity\tsemimajoraxis\ttime\tvx\tvy\tvz\tx\ty\tz\n')
    info.write('-\tnone\tAU\tyr\tkms\tkms\tkms\tAU\tAU\tAU\n')
    #info.write('========\t===========\t===========\t===========\t===========\t===========\t===========\t===========\t===========\t===========\n')
    info.write('key_star\tinitial_eccentricity_ffp_star\tinitial_semimajoraxis_ffp_star\ttime\tvx_star\tvy_star\tvz_star\tx_star\ty_star\tz_star\n')
    info.write('key_ffp\teccentricity_ffp_star\tsemimajoraxis_ffp_star\ttime\tvx_ffp\tvy_ffp\tvz_ffp\tx_ffp\ty_ffp\tz_ffp\n')
    info.write('key_bp\teccentricity_bp_star\tsemimajoraxis_bp_star\ttime\tvx_bp\tvy_bp\tvz_bp\tx_bp\ty_bp\tz_bp\n')
    info.write('\nPARAMETERS:\n')
    info.write('iteration_number\tmass_ffp(MJupiter)\timpact_parameter(AU)\tphi(degrees)\tduration_run(s)\tenergy_conservation\n')
    info.close()

if __name__ in ('__main__', '__plot__'):

    #Create file with information of units and order of particle sets
    create_info_file()

    #Fixed Values
    t_end = 650.0 #yr
    m0 = 1.0 #MSun
    m_planets = [1.0] #MJupiter
    a_planets = [5.0] #AU
    e_planets = [0.0]
    n_steps = 10000
    n_snapshots = int(n_steps/20.0)

    #earth_mass = 0.0031452 #In MJupiter
    #ms_ffp = np.linspace(earth_mass, 10.0, 10) #MJupiter
    ms_ffp = [1.0] #MJupiter

    #np.random.seed(12)
    #phis = np.random.rand(10) #degrees (this should be 1000 phis)
    phis = [0.0] #degrees

    #File to follow this parameters
    file_parameters = open('./parameters.txt','w')

    i=1

    for m_ffp in ms_ffp:
        #mass_ratio = m_ffp/m_planets[-1]
        #b_limit = find_limit_b(mass_ratio)
        #bs = squared_space(-b_limit, b_limit, 10) #AU
        bs = [5.0] #AU
        for b in bs:
            for phi in phis:
                #Time starts
                start_time = time.time()

                max_energy_change = ffp.run_capture(t_end_p=t_end,
                                                    m0_p=m0,
                                                    m_ffp_p=m_ffp,
                                                    m_planets_p=m_planets,
                                                    a_planets_p=a_planets,
                                                    e_planets_p=e_planets,
                                                    n_steps_p=n_steps,
                                                    phi_p=phi,
                                                    b_p=b,
                                                    iteration_number=i,
                                                    n_snapshots_p=n_snapshots)

                #Time stops
                duration = time.time()-start_time

                #Write in File (iteration, mass of ffp, impact parameter, phi, duration)
                file_parameters.write(('%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.9f\n')%(i,m_ffp,b,phi,duration,max_energy_change))

                #Advance counter
                i += 1

    file_parameters.close()
