import numpy as np
import ffp

#Fixed Values
t_end = 650.0 #yr
m0 = 1.0 #MSun
m_planets = [1.0] #MJupiter
a_planets = [5.0] #AU
e_planets = [0.0]
n_steps = 10000

#3 parameters
ms_ffp = [1.0] #MJupiter
bs = [5.0] #AU
#phis = np.random.rand(100) #degrees
phis = [0.0] #degrees

for m_ffp in ms_ffp:
    for b in bs:
        for phi in phis:
            ffp.run_capture(t_end_p=t_end,
                            m0_p=m0,
                            m_ffp_p=m_ffp,
                            m_planets_p=m_planets,
                            a_planets_p=a_planets,
                            e_planets_p=e_planets,
                            n_steps_p=n_steps,
                            phi_p=phi,
                            b_p=b)
