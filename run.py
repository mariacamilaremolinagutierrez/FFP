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

# def find_upper_b(a): #a is lambda
#     return math.sqrt(-(-a**6-6*a**5-15*a**4-20*a**3-15*a**2-(a+1)**3*sqrt(a**6+6*a**5+25615*a**4+20*a**3+15*a**2+6*a+1)-6*a-1)/(8*a**4))

#3 parameters
#earth_mass = 0.0031452 #In MJupiter
#ms_ffp = np.linspace(earth_mass, 10.0, 10) #MJupiter
ms_ffp = [1.0] #MJupiter
bs = squared_space(-7*5.0, 7*5.0, 10) #AU
#bs = [5.0] #AU
phis = np.random.rand(10) #degrees (this should be 1000 phis)
#phis = [0.0] #degrees

#File to follow this parameters
file_parameters = open('/parameters.txt','w')

i=1

for m_ffp in ms_ffp:
    for b in bs:
        for phi in phis:
            #Time starts
            start_time = time.time()

            ffp.run_capture(t_end_p=t_end,
                            m0_p=m0,
                            m_ffp_p=m_ffp,
                            m_planets_p=m_planets,
                            a_planets_p=a_planets,
                            e_planets_p=e_planets,
                            n_steps_p=n_steps,
                            phi_p=phi,
                            b_p=b,
                            iteration_number=i)

            #Time stops
            duration = time.time()-start_time

            #Write in File (iteration, mass of ffp, impact parameter, phi, duration)
            file_parameters.write(('%d\t%.4f\t%.4f\t%.4f\t%.4f\n')%(i,m_ffp,b,phi,duration))

            #Advance counter
            i += 1

file_parameters.close()
