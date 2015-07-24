import matplotlib.pyplot as plt
import numpy as np
import os

os.system('mkdir plots_dt100/')

parameters = np.loadtxt('./particles/m2.5000e+00/parameters_m2.5000e+00.txt')

number_times_pericenter = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25]

def plot_result(index, name):

    fig = plt.figure()
    plt.plot(number_times_pericenter, parameters[:,index])
    plt.title(name)
    plt.xlabel('number_times_pericenter')
    plt.ylabel(name)
    plt.savefig('plots_dt100/'+name+'.png')
    plt.close()

indices_results = [6,7,8,9,10,11,12,13,14,15,16,17,18]
names_results = ['e_star_ffp','e_star_bp','sma_star_ffp','sma_star_bp','inc_star_ffp','inc_star_bp','lan_star_ffp','lan_star_bp','ap_star_ffp','ap_star_bp','duration','max_energy_change','t_end']

#(m_bp,a_bp,b_ffp,phi_bp,inc_bp,lan_bp,e_star_ffp,e_star_bp,sma_star_ffp,sma_star_bp,inc_star_ffp,inc_star_bp,lan_star_ffp,lan_star_bp,ap_star_ffp,ap_star_bp,duration,max_energy_change,t_end)

for i in range(len(indices_results)):
    plot_result(indices_results[i], names_results[i])
