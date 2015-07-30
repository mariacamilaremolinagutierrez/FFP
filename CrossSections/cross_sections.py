import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

rc('font',**{'size': 16})
rc('text', usetex=True)

def create_parameters_and_stables_files():

    os.system('mkdir results/')

    masses_directories = os.listdir('./particles/')

    with open('./results/parameters.txt', 'w') as outfile:
        outfile.write('m_ffp\tv_ffp\tb_ffp\tphi_bp\tinc_bp\tlan_bp\tap_bp\te_star_ffp\te_star_bp\tsma_star_ffp\tsma_star_bp\tinc_star_ffp\tinc_star_bp\tlan_star_ffp\tlan_star_bp\tap_star_ffp\tap_star_bp\trun_time\tenergy_change\tintegration_time\n')
        for i in range(len(masses_directories)):
            mass_dir = masses_directories[i]
            fname = './particles/'+mass_dir+'/parameters_'+mass_dir+'.txt'
            with open(fname) as infile:
                for line in infile:
                    outfile.write(line)
            infile.close()
    outfile.close()

    with open('./results/stables.txt', 'w') as outfile:
        outfile.write('m_ffp\tv_ffp\tb_ffp\tphi_bp\tinc_bp\tlan_bp\tap_bp\te_star_ffp\te_star_bp\tsma_star_ffp\tsma_star_bp\tinc_star_ffp\tinc_star_bp\tlan_star_ffp\tlan_star_bp\tap_star_ffp\tap_star_bp\trun_time\tenergy_change\tintegration_time\n')
        for i in range(len(masses_directories)):
            mass_dir = masses_directories[i]
            fname = './particles/'+mass_dir+'/stables_'+mass_dir+'.txt'
            with open(fname) as infile:
                for line in infile:
                    outfile.write(line)
            infile.close()
    outfile.close()

def read_parameters_df():

    df = pd.read_csv('./results/parameters.txt', sep='\t', dtype=np.float64)

    ms_ffp = df['m_ffp']
    vs_ffp = df['v_ffp']
    bs_ffp = df['b_ffp']

    return df, ms_ffp, vs_ffp, bs_ffp

def read_stables_df():

    df = pd.read_csv('./results/stables.txt', sep='\t', dtype=np.float64)

    ms_ffp = df['m_ffp']
    vs_ffp = df['v_ffp']
    bs_ffp = df['b_ffp']

    return df, ms_ffp, vs_ffp, bs_ffp
    
def create_plots_folders():
    os.system('mkdir plots/')
    
def calculate_cross_section(b_ffp_max, number_captures, number_total):
    
    cs = np.pi*b_ffp_max**2*number_captures/number_total
    
    return cs
    
def plot_cross_sections_vs_velocity(m_ffp, vs_ffp_captured, bs_ffp_captured, vs_ffp_total, bs_ffp_total):
    
    x_axis = []
    y_axis = []
    
    unique_velocities, number_unique_velocities = np.unique(vs_ffp_total, return_counts=True)
    
    for uniq_vel in unique_velocities:
        
        indices_captured = np.where(vs_ffp_captured == uniq_vel)
        indices_total = np.where(vs_ffp_total == uniq_vel)
        bs_ffp_captured_for_v = bs_ffp_captured[indices_captured]
        
        if (len(bs_ffp_captured_for_v)==0):
            b_ffp_captured_for_v_max = 0.0
        else:
            b_ffp_captured_for_v_max = max(bs_ffp_captured_for_v)
        
        cs = calculate_cross_section(b_ffp_captured_for_v_max, len(indices_captured[0]), len(indices_total[0]))
        
        if(cs != 0.0):
            log_cs = np.log10(cs)
        else:
            log_cs = 0.0
            
        x_axis.append(uniq_vel)
        y_axis.append(log_cs)
        
    plt.figure(figsize=(20,10))    
    
    plt.plot(x_axis, y_axis, c='black')
    plt.scatter(x_axis, y_axis, c='black',s=2)
    plt.axvline(unique_velocities[0]*2, 0.0, max(y_axis), c='black', linestyle='--', alpha=0.6)
    plt.xlim(0,unique_velocities[-1]*1.01)
    
    m_ffp_string = r"${0:.3f}$".format(m_ffp)
    
    plt.title('$\mathrm{Cross\quad Sections:}\quad m_{FFP}=$'+m_ffp_string+'$\mathrm{\quad (MJupiter)}$')    
    plt.xlabel('$v_{FFP}\quad \mathrm{(km/s)}$')    
    plt.ylabel('$\log(\sigma_{capt})$')    
    
    plt.savefig('./plots/cross_sections_'+str(m_ffp)+'.png')
    plt.close()

if __name__ in ('__main__', '__plot__'):
    
    #Parameters and Stables Files
    create_parameters_and_stables_files()

    #Read df
    df, ms_ffp, vs_ffp, bs_ffp = read_stables_df()
    df_par, ms_ffp_par, vs_ffp_par, bs_ffp_par = read_parameters_df()
    
    ms_ffp_par_uniq = np.unique(ms_ffp_par)
    
    print ms_ffp_par_uniq
    
    print 'Stables: ', len(ms_ffp), '/', len(ms_ffp_par)
    
    #Plots
    create_plots_folders()
    
    #Make statistics plots
    for m_ffp in ms_ffp_par_uniq:
        
        vs_ffp_captured = np.array(vs_ffp)[np.where(ms_ffp == m_ffp)]
        bs_ffp_captured = np.array(bs_ffp)[np.where(ms_ffp == m_ffp)]
        
        vs_ffp_total = np.array(vs_ffp_par)[np.where(ms_ffp_par == m_ffp)]
        bs_ffp_total = np.array(bs_ffp_par)[np.where(ms_ffp_par == m_ffp)]
        
        plot_cross_sections_vs_velocity(m_ffp, vs_ffp_captured, bs_ffp_captured, vs_ffp_total, bs_ffp_total)



