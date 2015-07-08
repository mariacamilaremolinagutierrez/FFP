import os
import numpy as np

def create_info_file():
    info = open('./info.txt','w')
    info.write('PARAMETERS (All of the ones run):\n')
    info.write('mass_bp(MJupiter)\tsemimajoraxis_bp(AU)\timpact_parameter_ffp(AU)\tphi_bp(degrees)\teccentricity_star_ffp\teccentricity_star_bp\tsemimajoraxis_star_ffp(AU)\tsemimajoraxis_star_bp(AU)\tduration_run(s)\tenergy_conservation\n')
    info.write('\nSTABLES (Only the Hill stables):\n')
    info.write('mass_bp(MJupiter)\tsemimajoraxis_bp(AU)\timpact_parameter_ffp(AU)\tphi_bp(degrees)\teccentricity_star_ffp\teccentricity_star_bp\tsemimajoraxis_star_ffp(AU)\tsemimajoraxis_star_bp(AU)\tduration_run(s)\tenergy_conservation\n')
    info.close()

if __name__ in ('__main__', '__plot__'):

    #Create file with information of units and order of particle sets
    create_info_file()

    n_ms_bp = 10
    ms_bp = np.linspace(0.1,5.0,n_ms_bp) #MJupiter

    for m_bp in ms_bp:
        os.system('(mpiexec amuse.sh run.py '+str(m_bp)+' >> log 2>&1 &)')
