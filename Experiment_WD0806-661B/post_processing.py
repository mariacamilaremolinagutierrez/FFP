import os, math
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.pyplot as plt

from sympy.solvers import solve
from sympy import Symbol

from amuse.support import io
from amuse.units import units
from amuse.datamodel import Particles
from amuse.ext.orbital_elements import orbital_elements_from_binary

def stop_code():
    print '\nSTOP'
    import sys
    sys.exit()

def create_parameters_and_status_file():

    masses_directories = os.listdir('./particles/')

    with open('./parameters.txt', 'w') as outfile:
        outfile.write('m_bp\ta_bp\tb_ffp\tphi_bp\te_star_ffp\te_star_bp\tsma_star_ffp\tsma_star_bp\trun_time\tenergy_change\n')
        for i in range(len(masses_directories)):
            mass_dir = masses_directories[i]
            fname = './particles/'+mass_dir+'/parameters_'+mass_dir+'.txt'
            with open(fname) as infile:
                for line in infile:
                    outfile.write(line)
            infile.close()
    outfile.close()

    with open('./stables.txt', 'w') as outfile:
        outfile.write('m_bp\ta_bp\tb_ffp\tphi_bp\te_star_ffp\te_star_bp\tsma_star_ffp\tsma_star_bp\trun_time\tenergy_change\n')
        for i in range(len(masses_directories)):
            mass_dir = masses_directories[i]
            fname = './particles/'+mass_dir+'/stables_'+mass_dir+'.txt'
            with open(fname) as infile:
                for line in infile:
                    outfile.write(line)
            infile.close()
    outfile.close()

def read_parameters_df():

    df = pd.read_csv('./parameters.txt', sep='\t')

    ms_bp = df['m_bp']
    as_bp = df['a_bp']
    bs_ffp = df['b_ffp']
    phis_bp = df['phi_bp']

    return ms_bp, as_bp, bs_ffp, phis_bp

def read_stables_df():

    df = pd.read_csv('./stables.txt', sep='\t', dtype=np.float64)

    ms_bp = df['m_bp']
    as_bp = df['a_bp']
    bs_ffp = df['b_ffp']
    phis_bp = df['phi_bp']

    return df, ms_bp, as_bp, bs_ffp, phis_bp

def plot_trajectory(x, y, filename):

    f = plt.figure(figsize=(70,30))

    x_star = x[:,0]-x[:,0]
    x_ffp = x[:,1]-x[:,0]

    y_star = y[:,0]-y[:,0]
    y_ffp = y[:,1]-y[:,0]

    x_planet = x[:,2]-x[:,0]
    y_planet = y[:,2]-y[:,0]

    plt.plot(x_star,y_star,c='y',label='Star')
    plt.scatter(x_star[0],y_star[0],c='black',marker='*', lw = 0)
    plt.scatter(x_star[-1],y_star[-1],c='y',marker='*', lw = 0)

    plt.plot(x_ffp,y_ffp,c='c',label='FFP', lw = 2)
    plt.scatter(x_ffp[0],y_ffp[0],c='black', lw = 0)
    plt.scatter(x_ffp[-1],y_ffp[-1],c='c', lw = 0)

    plt.plot(x_planet,y_planet,c='m',label='BP',alpha=0.5)
    plt.scatter(x_planet[0],y_planet[0],c='black', lw = 0)
    plt.scatter(x_planet[-1],y_planet[-1],c='m', lw = 0)

    plt.axhline(y=0, xmin=-80, xmax=10, c='black', linestyle='--')
    plt.axvline(x=0, ymin=-5, ymax=2, c='black', linestyle='--')

    plt.title('Trajectory FFP', fontsize=40)
    plt.axes().set_aspect('equal', 'datalim')
    plt.xlabel('$x$ (AU)', fontsize=40)
    plt.ylabel('$y$ (AU)', fontsize=40)
    plt.legend(fontsize=40)
    plt.savefig('./plots/trajectories/'+filename+'.png')

    plt.close()

def plot_orbital_elements(times, eccs, smas, filename):

    f = plt.figure(figsize=(35,15))

    ecc_starbp_ffp = eccs[:,1]
    ecc_star_bp = eccs[:,2]

    sma_starbp_ffp = smas[:,1]
    sma_star_bp = smas[:,2]

    #eccentricity
    subplot = f.add_subplot(1,2,1)

    subplot.plot(times,ecc_starbp_ffp,c='red',label='FFP and Star+BP')
    subplot.plot(times,ecc_star_bp,c='green',label='BP and Star')

    subplot.set_title('Eccentricity', fontsize=20)
    subplot.set_xlabel('$t$ (yr)', fontsize=20)
    subplot.set_ylabel('$e$', fontsize=20)
    subplot.set_ylim(-0.5,1.5)
    subplot.legend(fontsize=20)

    #semimajoraxis
    subplot = f.add_subplot(1,2,2)

    subplot.plot(times,sma_starbp_ffp,c='red',label='FFP and Star+BP')
    subplot.plot(times,sma_star_bp,c='green',label='BP and Star')

    subplot.set_title('Semimajor Axis', fontsize=20)
    subplot.set_xlabel('$t$ (yr)', fontsize=20)
    subplot.set_ylabel('$a$ (AU)', fontsize=20)
    subplot.legend(fontsize=20)

    plt.savefig('./plots/orbital_elements/'+filename+'.png')
    plt.close()

def create_plots_folders():
    os.system('mkdir plots/')
    os.system('mkdir plots/statistics')

def uniq_list(seq):
    #http://www.peterbe.com/plog/uniqifiers-benchmark
    seen = set()
    seen_add = seen.add
    return [ x for x in seq if not (x in seen or seen_add(x))]

def plot_histograms(parameter, amount, df_name):

    f = plt.figure(figsize=(30,15))

    parameter = np.array(parameter)

    if (amount == 10):
        n,b,p = plt.hist(parameter, color = 'c')
    else:
        n,b,p = plt.hist(parameter, bins = int(amount/20.0), color = 'c')

    plt.title('Histogram of captures for each '+df_name, fontsize=20)
    plt.xlabel(df_name, fontsize=20)
    plt.ylabel('number of captures', fontsize=20)
    plt.ylim(0,max(n)*1.05)
    plt.savefig('./plots/statistics/number_captures_'+df_name+'.png')
    plt.close()

def plot_parameters(parameters_x_axis, parameters_y_axis, parameter_color, df_names, latex_names):

    #names: [x axis parameter, y axis parameter, color parameter, fixed parameter]

    f = plt.figure(figsize=(30,15))

    sc = plt.scatter(parameters_x_axis, parameters_y_axis, c=parameter_color, vmin=np.amin(parameter_color), vmax=np.amax(parameter_color), lw=0, s=10)
    cbar = plt.colorbar(sc, orientation="horizontal")

    plt.title('$\mathrm{CAPTURES}$', fontsize=20)
    plt.xlabel(latex_names[0], fontsize=30)
    plt.ylabel(latex_names[1], fontsize=30)
    cbar.set_label(latex_names[2], rotation=0, fontsize=30)

    plt.savefig('./plots/statistics/parameters_'+df_names[0]+'_'+df_names[1]+'_'+df_names[2]+'_'+df_names[3]+'.png')
    plt.close()

def plot_threesome(parameters_xaxis, parameters_color, parameters_marker, df_names, latex_names):
    #names: [x axis parameter, color parameter, marker parameter]

    num_lines = len(parameters_xaxis)
    list_to_uniq = []

    for i in range(num_lines):
        list_to_uniq.append(str(parameters_xaxis[i])+'\t'+str(parameters_color[i])+'\t'+str(parameters_marker[i]))

    unique_combinations, number_unique_combinations = np.unique(list_to_uniq, return_counts=True)

    xaxis_values = []
    colorbar_values = []
    marker_values = []

    for uc in unique_combinations:
        uc_parts = uc.split('\t')
        xaxis_values.append(float(uc_parts[0]))
        colorbar_values.append(float(uc_parts[1]))
        marker_values.append(float(uc_parts[2]))

    xaxis_values = np.array(xaxis_values)
    colorbar_values = np.array(colorbar_values)
    marker_values = np.array(marker_values)

    uniq_markers = np.array(uniq_list(marker_values))
    markers = [".",">","v","d","s","D","p","H","*","o"]
    marker_i = 0

    f = plt.figure(figsize=(40,20))

    m = cm.ScalarMappable(cmap=cm.jet)
    m.set_array(colorbar_values)

    for um in uniq_markers:

        par_indices = np.where(marker_values == um)
        sc = plt.scatter(xaxis_values[par_indices], number_unique_combinations[par_indices], c=colorbar_values[par_indices], lw=0, s=40, marker=markers[marker_i], label=latex_names[2]+'$ = $'+str(um))

        marker_i += 1

    cbar = plt.colorbar(m, orientation="horizontal")
    cbar.set_label(latex_names[1], fontsize=30)
    plt.title('$\mathrm{CAPTURES}$', fontsize=40)
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=20)
    plt.xlabel(latex_names[0], fontsize=30)
    plt.ylabel('$\mathrm{Captures}$', fontsize=30)
    plt.ylim(0, max(number_unique_combinations)*1.05)

    plt.savefig('./plots/statistics/threesome_'+df_names[0]+'_'+df_names[1]+'_'+df_names[2]+'.png')
    plt.close()


def make_statistical_plots(df, ms_bp, as_bp, bs_ffp, phis_bp):

    #Reading parameters dataframe
    ms_bp_par, as_bp_par, bs_ffp_par, phis_bp_par = read_parameters_df()

    #Getting all the ms, as, bs and phis that were combined
    ms_bp_par_uniq = uniq_list(ms_bp_par)
    as_bp_par_uniq = uniq_list(as_bp_par)
    bs_ffp_par_uniq = uniq_list(bs_ffp_par)
    phis_bp_par_uniq = uniq_list(phis_bp_par)

    #Numbers of each parameter
    n_ms_bp_par = len(ms_bp_par_uniq)
    n_as_bp_par = len(as_bp_par_uniq)
    n_bs_ffp_par = len(bs_ffp_par_uniq)
    n_phis_bp_par = len(phis_bp_par_uniq)

    print n_ms_bp_par, n_as_bp_par, n_bs_ffp_par, n_phis_bp_par

    #Number of captures
    # #m_bp
    # plot_histograms(ms_bp, n_ms_bp_par, 'm_bp')
    # #a_bp
    # plot_histograms(as_bp, n_as_bp_par, 'a_bp')
    # #b_ffp
    # plot_histograms(bs_ffp, n_bs_ffp_par, 'b_ffp')
    # #phi_bp
    # plot_histograms(phis_bp, n_phis_bp_par, 'phi_bp')

    #Plot parameters
    df_names = ['phi_bp', 'b_ffp', 'a_bp', 'm_bp']
    latex_names = ['$\phi_{BP}$', '$b_{FFP}$', '$a_{BP}$', '$m_{BP}$']
    plot_parameters(phis_bp, bs_ffp, as_bp, df_names, latex_names)

    df_names = ['phi_bp', 'a_bp', 'b_ffp', 'm_bp']
    latex_names = ['$\phi_{BP}$', '$a_{BP}$', '$b_{FFP}$', '$m_{BP}$']
    plot_parameters(phis_bp, as_bp, bs_ffp, df_names, latex_names)

    df_names = ['phi_bp', 'm_bp', 'a_bp', 'b_ffp']
    latex_names = ['$\phi_{BP}$', '$m_{BP}$', '$a_{BP}$', '$b_{FFP}$']
    plot_parameters(phis_bp, ms_bp, as_bp, df_names, latex_names)

    #Plot threesomes
    df_names = ['m_bp', 'b_ffp', 'a_bp']
    latex_names = ['$m_{BP} \quad \mathrm{(M_{Jupiter})}$', '$b_{FFP} \quad \mathrm{(AU)}$', '$a_{BP}$']
    plot_threesome(ms_bp, bs_ffp, as_bp, df_names, latex_names)

    df_names = ['b_ffp', 'm_bp', 'a_bp']
    latex_names = ['$b_{FFP} \quad \mathrm{(AU)}$', '$m_{BP} \quad \mathrm{(M_{Jupiter})}$', '$a_{BP}$']
    plot_threesome(bs_ffp, ms_bp, as_bp, df_names, latex_names)

    df_names = ['a_bp', 'b_ffp', 'm_bp']
    latex_names = ['$a_{BP} \quad \mathrm{(AU)}$', '$b_{FFP} \quad \mathrm{(AU)}$', '$m_{BP}$']
    plot_threesome(as_bp, bs_ffp, ms_bp, df_names, latex_names)

    df_names = ['b_ffp', 'a_bp', 'm_bp']
    latex_names = ['$b_{FFP} \quad \mathrm{(AU)}$', '$a_{BP} \quad \mathrm{(AU)}$', '$m_{BP}$']
    plot_threesome(bs_ffp, as_bp, ms_bp, df_names, latex_names)


if __name__ in ('__main__', '__plot__'):

    #Create file with all the runs made per line
    # create_parameters_and_status_file()

    #Read df
    df, ms_bp, as_bp, bs_ffp, phis_bp = read_stables_df()

    #Plottts
    # create_plots_folders()
    #Make statistics plots
    make_statistical_plots(df, ms_bp, as_bp, bs_ffp, phis_bp)
