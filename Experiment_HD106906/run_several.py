import os
import numpy as np

def create_info_file():
    info = open('./info.txt','w')
    info.write('PARTICLES:\n')
    info.write('key\teccentricity\tsemimajoraxis\ttime\tvx\tvy\tvz\tx\ty\tz\n')
    info.write('-\tnone\tAU\tyr\tkms\tkms\tkms\tAU\tAU\tAU\n')
    info.write('key_star\tinitial_eccentricity_ffp_star\tinitial_semimajoraxis_ffp_star\ttime\tvx_star\tvy_star\tvz_star\tx_star\ty_star\tz_star\n')
    info.write('key_ffp\teccentricity_ffp_star\tsemimajoraxis_ffp_star\ttime\tvx_ffp\tvy_ffp\tvz_ffp\tx_ffp\ty_ffp\tz_ffp\n')
    info.write('key_bp\teccentricity_bp_star\tsemimajoraxis_bp_star\ttime\tvx_bp\tvy_bp\tvz_bp\tx_bp\ty_bp\tz_bp\n')
    info.write('\nSTABLES:\n')
    info.write('folder\tfilename\tmass_bp(MJupiter)\tsemimajoraxis_bp(AU)\timpact_parameter_ffp(AU)\tphi_bp(degrees)\tduration_run(s)\tenergy_conservation\n')
    info.close()

if __name__ in ('__main__', '__plot__'):

    #Create file with information of units and order of particle sets
    create_info_file()

    n_ms_bp = 10
    ms_bp = np.linspace(0.1,5.0,n_ms_bp) #MJupiter

    for m_bp in ms_bp:
        os.system('(mpiexec amuse.sh run.py '+str(m_bp)+' >> log 2>&1 &)')
