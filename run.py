import numpy as np
import os, time
import ffp

def squared_space(b_initial, b_final, n_bs):
    if (b_initial<0.0): #assuming b_initial = b_final
        ss_half = np.sqrt(np.linspace(0.0,b_final**2,int(n_bs/2)))
        if(n_bs%2 == 0):
            ss = np.append(-1.0*ss_half,ss_half)
        else:
            ss = np.append(-1.0*ss_half,zeros(1),ss_half)
    else:
        ss = np.sqrt(np.linspace(b_initial**2,b_final**2,n_bs))
    return ss

#Fixed Values
t_end = 650.0 #yr
m0 = 1.0 #MSun
m_planets = [1.0] #MJupiter
a_planets = [5.0] #AU
e_planets = [0.0]
n_steps = 10000

#3 parameters
#earth_mass = 0.0031452 #In MJupiter
#ms_ffp = np.linspace(earth_mass, 10.0, 10) #MJupiter
ms_ffp = [1.0] #MJupiter
#bs = squared_space(-7*5.0, 7*5.0, 10) #AU
bs = [5.0] #AU
#phis = np.random.rand(100) #degrees
phis = [0.0] #degrees

#Files
filename_status = 'status.txt'
filename_results = 'results.txt'

status=open(filename_status,'a')
results = open(filename_results,'a')

num_lines = sum(1 for line in open(filename_status))

if(num_lines == 0):
    status.write('iteration_number\tduration(seconds)\n')
    results.write('m_ffp\tb\tphi\te_0_ffp\ta_0_ffp\te_0_bp\ta_0_bp\te_ffp_bp\ta_ffp_bp\n')#still don't know what to print here
    num_lines += 1

#Counter
i = 1

for m_ffp in ms_ffp:
    for b in bs:
        for phi in phis:

            if (i>=num_lines):
                #Time starts
                start_time = time.time()

                e_0_ffp,a_0_ffp,e_0_bp,a_0_bp,e_ffp_bp,a_ffp_bp = ffp.run_capture(t_end_p=t_end,
                                                                                    m0_p=m0,
                                                                                    m_ffp_p=m_ffp,
                                                                                    m_planets_p=m_planets,
                                                                                    a_planets_p=a_planets,
                                                                                    e_planets_p=e_planets,
                                                                                    n_steps_p=n_steps,
                                                                                    phi_p=phi,
                                                                                    b_p=b)

                #Time stops
                duration = time.time()-start_time

                #Write status and results
                status.write(('%d\t%.4f\n')%(i,duration))
                results.write(('%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n')%(m_ffp,b,phi,e_0_ffp,a_0_ffp,e_0_bp,a_0_bp,e_ffp_bp,a_ffp_bp))

            #Advance counter
            i += 1

#Closing files
status.close()
results.close()
