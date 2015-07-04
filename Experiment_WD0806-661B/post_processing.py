import os, math
import pandas as pd
import numpy as np
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
        outfile.write('folder\tfilename\tm_bp\ta_bp\tb_ffp\tphi_bp\trun_time\tenergy_change\n')
        for i in range(len(masses_directories)):
            mass_dir = masses_directories[i]
            fname = './particles/'+mass_dir+'/parameters_'+mass_dir+'.txt'
            with open(fname) as infile:
                for line in infile:
                    outfile.write(mass_dir+'\t')
                    outfile.write(line)
            infile.close()
    outfile.close()

    with open('./stables.txt', 'w') as outfile:
        outfile.write('filename\tm_bp\ta_bp\tb_ffp\tphi_bp\trun_time\tenergy_change\n')
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

    filenames = df['filename']
    ms_bp = df['m_bp']
    as_bp = df['a_bp']
    bs_ffp = df['b_ffp']
    phis_bp = df['phi_bp']

    return ms_bp, as_bp, bs_ffp, phis_bp

def read_stables_df():

    df = pd.read_csv('./stables.txt', sep='\t')

    filenames = df['filename']
    ms_bp = df['m_bp']
    as_bp = df['a_bp']
    bs_ffp = df['b_ffp']
    phis_bp = df['phi_bp']

    return df, filenames, ms_bp, as_bp, bs_ffp, phis_bp

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

def make_individual_plots(filenames, ms_bp, as_bp, bs_ffp, phis_bp):

    num_files = len(filenames)

    for i in range(num_files):

        filename = filenames[i]
        m_bp = ms_bp[i]
        a_bp = as_bp[i]
        b_ffp = bs_ffp[i]
        phi_bp = phis_bp[i]

        #Bodies Order: star - ffp - bp
        bodies = io.read_set_from_file('./stables/'+filename, 'hdf5')

        #To append body history
        eccs = []
        smas = []
        xs = []
        ys = []
        times = []

        for data in bodies.history:

            #Order: star - ffp - bp
            e_values = data.eccentricity
            a_values = data.semimajoraxis.value_in(units.AU)
            x_values = data.x.value_in(units.AU)
            y_values = data.y.value_in(units.AU)
            t_value = data.time.value_in(units.yr)[0]

            eccs.append(e_values)
            smas.append(a_values)
            xs.append(x_values)
            ys.append(y_values)
            times.append(t_value)

        os.system('mkdir plots/trajectories')
        os.system('mkdir plots/orbital_elements')

        plot_trajectory(np.array(xs), np.array(ys), filename)
        plot_orbital_elements(np.array(times), np.array(eccs), np.array(smas), filename)

def uniq_list(seq):
    #http://www.peterbe.com/plog/uniqifiers-benchmark
    seen = set()
    seen_add = seen.add
    return [ x for x in seq if not (x in seen or seen_add(x))]

def plot_parameter_number_captures(parameter, n_total_parameters, df, df_name):

    f = plt.figure(figsize=(30,15))

    max_number = 0

    for par in parameter:

        subset = np.where(df[df_name] == par)
        number = len(subset[0])
        plt.scatter(float(par), number, c='r', lw=0, s=30)

        if(number > max_number):
            max_number = number

    plt.title('Number of '+df_name+' that turned a capture for each all of the ones tried', fontsize=20)
    plt.xlabel(df_name, fontsize=20)
    plt.ylabel('number of captures', fontsize=20)
    plt.ylim(-0.02*max_number,max_number*1.02)
    plt.savefig('./plots/statistics/number_captures_'+df_name+'.png')
    plt.close()

def plot_parameters(parameters_x_axis, parameters_y_axis, parameter_color, parameters_shape, df, df_names, latex_names):

    #names: [x axis parameter, y axis parameter, color parameter, shape parameter]

    markers = ["o",",",".","v","^","<",">","1","2","3","4","8","s","p","*","h","H","+","D","d","|","_"]

    f = plt.figure(figsize=(30,15))

    marker_i = 0

    for par_shape in parameters_shape:

        par_indices = np.where(df[df_names[3]] == par_shape)

        sc = plt.scatter(np.array(parameters_x_axis)[par_indices], np.array(parameters_y_axis)[par_indices], c=parameter_color, vmin=np.amin(parameter_color), vmax=np.amax(parameter_color), lw=0, s=60, marker = markers[marker_i], label = latex_names[3]+' = '+str(par_shape))
        cbar = plt.colorbar(sc)

        marker_i += 1

    plt.title('Captures', fontsize=20)
    plt.legend(fontsize=30)
    plt.xlabel(latex_names[0], fontsize=30)
    plt.ylabel(latex_names[1], fontsize=30)
    cbar.set_label(latex_names[2], rotation=0, fontsize=30)

    plt.savefig('./plots/statistics/parameters_'+df_names[0]+'_'+df_names[1]+'_'+df_names[2]+'_'+df_names[3]+'.png')
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
    #m_bp
    plot_parameter_number_captures(ms_bp_par_uniq, n_ms_bp_par, df, 'm_bp')
    #a_bp
    plot_parameter_number_captures(as_bp_par_uniq, n_as_bp_par, df, 'a_bp')
    #b_ffp
    plot_parameter_number_captures(bs_ffp_par_uniq, n_bs_ffp_par, df, 'b_ffp')
    #phi_bp
    plot_parameter_number_captures(phis_bp_par_uniq, n_phis_bp_par, df, 'phi_bp')

    # df_names = ['phi_bp', 'b_ffp', 'a_bp', 'm_bp']
    # latex_names = ['$\phi_{BP}$', '$b_{FFP}$', '$a_{BP}$', '$m_{BP}$']
    # plot_parameters(phis_bp, bs_ffp, as_bp, uniq_list(ms_bp), df, df_names, latex_names)
    #
    # df_names = ['phi_bp', 'a_bp', 'b_ffp', 'm_bp']
    # latex_names = ['$\phi_{BP}$', '$a_{BP}$', '$b_{FFP}$', '$m_{BP}$']
    # plot_parameters(phis_bp, as_bp, bs_ffp, uniq_list(ms_bp), df, df_names, latex_names)

    # df_names = ['phi_bp', 'm_bp', 'a_bp', 'b_ffp']
    # latex_names = ['$\phi_{BP}$', '$m_{BP}$', '$a_{BP}$', '$b_{FFP}$']
    # plot_parameters(phis_bp, ms_bp, as_bp, uniq_list(bs_ffp), df, df_names, latex_names)

if __name__ in ('__main__', '__plot__'):

    #Create file with all the runs made per line
    create_parameters_and_status_file()

    #Read df
    df, filenames, ms_bp, as_bp, bs_ffp, phis_bp = read_stables_df()

    #Plottts
    create_plots_folders()
    #Make individual plots
    #make_individual_plots(filenames, ms_bp, as_bp, bs_ffp, phis_bp)
    #Make statistics plots
    make_statistical_plots(df, ms_bp, as_bp, bs_ffp, phis_bp)
