import os
import numpy as np

def create_info_file(info_filename):
    info = open(info_filename, 'w')
    info.write('PARAMETERS (All of the ones run):\n')
    info.write('mass_ffp(MJupiter)\tvelocity_ffp(AU)\timpact_parameter_ffp(AU)\ttrue_anomaly_bp(degrees)\tinclination_bp(degrees)\tlogitude_of_ascending_node_bp(degrees)\targument_of_periastron_bp(degrees)\teccentricity_star_ffp\teccentricity_star_bp\tsemimajoraxis_star_ffp(AU)\tsemimajoraxis_star_bp(AU)\tinclination_star_ffp(degrees)\tinclination_star_bp(degrees)\tlongitude_of_ascending_node_star_ffp(degrees)\tlongitude_of_ascending_node_star_bp(degrees)\targument_of_pericenter_star_ffp(degrees)\targument_of_pericenter_star_bp(degrees)\tduration_run(s)\tenergy_conservation\ttime_integration(yr)\n')
    info.write('\nSTABLES (Only the Hill stables):\n')
    info.write('mass_ffp(MJupiter)\tvelocity_ffp(AU)\timpact_parameter_ffp(AU)\ttrue_anomaly_bp(degrees)\tinclination_bp(degrees)\tlogitude_of_ascending_node_bp(degrees)\targument_of_periastron_bp(degrees)\teccentricity_star_ffp\teccentricity_star_bp\tsemimajoraxis_star_ffp(AU)\tsemimajoraxis_star_bp(AU)\tinclination_star_ffp(degrees)\tinclination_star_bp(degrees)\tlongitude_of_ascending_node_star_ffp(degrees)\tlongitude_of_ascending_node_star_bp(degrees)\targument_of_pericenter_star_ffp(degrees)\targument_of_pericenter_star_bp(degrees)\tduration_run(s)\tenergy_conservation\ttime_integration(yr)\n')
    info.close()

if __name__ in ('__main__', '__plot__'):

    #Create file with information of units and order of results    
    info_filename = './info.txt'
    if (os.path.isfile(info_filename) == False):
        create_info_file()

    MEarth_in_Mjupiter = 0.0031452 #MJupiter 
    
    ms_ffp_1 = np.linspace(MEarth_in_Mjupiter, 1.0, 6) #MJupiter
    ms_ffp_2 = np.linspace(1.0, 10.0, 5) #MJupiter
    
    ms_ffp = np.concatenate((ms_ffp_1[:5],ms_ffp_2))
    
    print ms_ffp
    
    os.system('mkdir logs/')

#    for m_ffp in ms_ffp:
#        os.system('(amuse.sh run.py '+str(m_ffp)+' >> logs/log'+str(m_ffp)+'.txt 2>&1 &)')

    m_ffp = 1.0
    os.system('(amuse.sh run.py '+str(m_ffp)+' >> logs/log'+str(m_ffp)+'.txt 2>&1 &)')