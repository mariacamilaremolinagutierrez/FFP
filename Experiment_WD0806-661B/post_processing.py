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

def create_stables_folder():

    os.system('mkdir stables/')

    masses_directories = os.listdir('./particles/')

    for mass_dir in masses_directories:
        #Copy the file to other folder
        os.system('cp particles/'+mass_dir+'/*.hdf5 stables/')

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

    # print type(sma_star_bp)
    # print np.shape(sma_star_bp)
    # print sma_star_bp

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
    os.system('mkdir plots/trajectories')
    os.system('mkdir plots/orbital_elements')
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

        plot_trajectory(np.array(xs), np.array(ys), filename)
        plot_orbital_elements(np.array(times), np.array(eccs), np.array(smas), filename)

def uniq_list(seq):
    #http://www.peterbe.com/plog/uniqifiers-benchmark
    seen = set()
    seen_add = seen.add
    return [ x for x in seq if not (x in seen or seen_add(x))]

def plot_combination_percentage(first_iteration, second_iteration, third_iteration, parameter, n_total_parameters, df, df_names):

    combinations = []
    percentages = []

    combination = 1

    for fi in first_iteration:
        for si in second_iteration:
            for ti in third_iteration:

                subset = np.where((df[df_names[0]] == fi) & (df[df_names[1]] == si) & (df[df_names[2]] == ti))

                combinations.append(combination)
                percentages.append(len(subset[0])*100.0/float(n_total_parameters))

                combination += 1

    f = plt.figure(figsize=(30,15))

    plt.scatter(combinations,percentages,c='r')

    plt.title('Percentages of '+df_names[3]+' that turned a capture for each combination of the other three', fontsize=20)
    plt.xlabel('Combination Number', fontsize=20)
    plt.ylabel('Percentage of '+df_names[3], fontsize=20)
    plt.savefig('./plots/statistics/percentage_'+df_names[3]+'.png')

    plt.close()

def make_percentages_plots(df):

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

    #Combinations
    #m_bp
    df_names = ['a_bp','b_ffp','phi_bp','m_bp']
    plot_combination_percentage(as_bp_par_uniq, bs_ffp_par_uniq, phis_bp_par_uniq, ms_bp_par_uniq, n_ms_bp_par, df, df_names)
    #a_bp
    df_names = ['b_ffp','phi_bp','m_bp','a_bp']
    plot_combination_percentage(bs_ffp_par_uniq, phis_bp_par_uniq, ms_bp_par_uniq, as_bp_par_uniq, n_as_bp_par, df, df_names)
    #b_ffp
    df_names = ['phi_bp','m_bp','a_bp','b_ffp']
    plot_combination_percentage(phis_bp_par_uniq, ms_bp_par_uniq, as_bp_par_uniq, bs_ffp_par_uniq, n_bs_ffp_par, df, df_names)
    #phi_bp
    df_names = ['m_bp','a_bp','b_ffp','phi_bp']
    plot_combination_percentage(ms_bp_par_uniq, as_bp_par_uniq, bs_ffp_par_uniq, phis_bp_par_uniq, n_phis_bp_par, df, df_names)


if __name__ in ('__main__', '__plot__'):

    #Create file with all the runs made per line
    #create_parameters_and_status_file()

    #Move all the stables to another folder
    #create_stables_folder()

    #Read df
    df, filenames, ms_bp, as_bp, bs_ffp, phis_bp = read_stables_df()

    #Plottts
    #create_plots_folders()
    #Make individual plots
    #make_individual_plots(filenames, ms_bp, as_bp, bs_ffp, phis_bp)
    #Make statistics plots
    make_percentages_plots(df)
    
