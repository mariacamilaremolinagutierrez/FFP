import os
import numpy as np

if __name__ in ('__main__', '__plot__'):

    n_ms_bp = 10
    ms_bp = np.linspace(0.1,5.0,n_ms_bp) #MJupiter
    
    m_bp = ms_bp[6]
    
    os.system( 'amuse.sh run.py '+str(m_bp) )
