import os
import numpy as np

#def create_info_file():
#    info = open('./info.txt','w')
#    info.write('PARAMETERS (All of the ones run):\n')
#    info.write('mass_bp(MJupiter)\tsemimajoraxis_bp(AU)\timpact_parameter_ffp(AU)\tphi_bp(degrees)\teccentricity_star_ffp\teccentricity_star_bp\tsemimajoraxis_star_ffp(AU)\tsemimajoraxis_star_bp(AU)\tinclination_star_ffp(degrees)\tinclination_star_bp(degrees)\tlongitude_of_ascending_node_star_ffp(degrees)\tlongitude_of_ascending_node_star_bp(degrees)\targument_of_pericenter_star_ffp(degrees)\targument_of_pericenter_star_bp(degrees)\tduration_run(s)\tenergy_conservation\ttime_integration(yr)\n')
#    info.write('\nSTABLES (Only the Hill stables):\n')
#    info.write('mass_bp(MJupiter)\tsemimajoraxis_bp(AU)\timpact_parameter_ffp(AU)\tphi_bp(degrees)\teccentricity_star_ffp\teccentricity_star_bp\tsemimajoraxis_star_ffp(AU)\tsemimajoraxis_star_bp(AU)\tinclination_star_ffp(degrees)\tinclination_star_bp(degrees)\tlongitude_of_ascending_node_star_ffp(degrees)\tlongitude_of_ascending_node_star_bp(degrees)\targument_of_pericenter_star_ffp(degrees)\targument_of_pericenter_star_bp(degrees)\tduration_run(s)\tenergy_conservation\ttime_integration(yr)\n')
#    info.close()

if __name__ in ('__main__', '__plot__'):

    #Create file with information of units and order of particle sets
#    create_info_file()

    n_ms_ffp = 1
    ms_ffp = np.linspace(1.0,10.0,n_ms_ffp) #MJupiter

    os.system('mkdir logs/')

    for m_ffp in ms_ffp:
        os.system('(amuse.sh run.py '+str(m_ffp)+' >> logs/log'+str(m_ffp)+'.txt 2>&1 &)')
