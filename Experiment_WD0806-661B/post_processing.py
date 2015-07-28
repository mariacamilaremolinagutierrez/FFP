import os
import pandas as pd
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt

def stop_code():
    print '\nSTOP'
    import sys
    sys.exit()

def read_parameters_df():

    df = pd.read_csv('./results/stables.txt', sep='\t', dtype=np.float64)

    ms_bp = df['m_bp']
    as_bp = df['a_bp']
    bs_ffp = df['b_ffp']
    phis_bp = df['phi_bp']
    incs_bp = df['inc_bp']
    lans_bp = df['lan_bp']
    es_star_bp = df['e_star_bp']
    es_star_ffp = df['e_star_ffp']
    smas_star_ffp = df['sma_star_ffp']
    smas_star_bp = df['sma_star_bp']
    incs_star_ffp = df['inc_star_ffp']
    incs_star_bp = df['inc_star_bp']
    lans_star_ffp = df['lan_star_ffp']
    lans_star_bp = df['lan_star_bp']
    aps_star_ffp = df['ap_star_ffp']
    aps_star_bp = df['ap_star_bp']

    return ms_bp, as_bp, bs_ffp, phis_bp, incs_bp, lans_bp, es_star_bp, es_star_ffp, smas_star_ffp, smas_star_bp, incs_star_ffp, incs_star_bp, lans_star_ffp, lans_star_bp, aps_star_ffp, aps_star_bp

def read_stables_df():

    df = pd.read_csv('./results/stables.txt', sep='\t', dtype=np.float64)

    ms_bp = df['m_bp']
    as_bp = df['a_bp']
    bs_ffp = df['b_ffp']
    phis_bp = df['phi_bp']
    incs_bp = df['inc_bp']
    lans_bp = df['lan_bp']
    es_star_bp = df['e_star_bp']
    es_star_ffp = df['e_star_ffp']
    smas_star_ffp = df['sma_star_ffp']
    smas_star_bp = df['sma_star_bp']
    incs_star_ffp = df['inc_star_ffp']
    incs_star_bp = df['inc_star_bp']
    lans_star_ffp = df['lan_star_ffp']
    lans_star_bp = df['lan_star_bp']
    aps_star_ffp = df['ap_star_ffp']
    aps_star_bp = df['ap_star_bp']

    return df, ms_bp, as_bp, bs_ffp, phis_bp, incs_bp, lans_bp, es_star_bp, es_star_ffp, smas_star_ffp, smas_star_bp, incs_star_ffp, incs_star_bp, lans_star_ffp, lans_star_bp, aps_star_ffp, aps_star_bp

def create_parameters_and_stables_files():

    os.system('mkdir results/')

    masses_directories = os.listdir('./particles/')

    with open('./results/parameters.txt', 'w') as outfile:
        outfile.write('m_bp\ta_bp\tb_ffp\tphi_bp\tinc_bp\tlan_bp\te_star_ffp\te_star_bp\tsma_star_ffp\tsma_star_bp\tinc_star_ffp\tinc_star_bp\tlan_star_ffp\tlan_star_bp\tap_star_ffp\tap_star_bp\trun_time\tenergy_change\tintegration_time\n')
        for i in range(len(masses_directories)):
            mass_dir = masses_directories[i]
            fname = './particles/'+mass_dir+'/parameters_'+mass_dir+'.txt'
            with open(fname) as infile:
                for line in infile:
                    outfile.write(line)
            infile.close()
    outfile.close()

    with open('./results/stables.txt', 'w') as outfile:
        outfile.write('m_bp\ta_bp\tb_ffp\tphi_bp\tinc_bp\tlan_bp\te_star_ffp\te_star_bp\tsma_star_ffp\tsma_star_bp\tinc_star_ffp\tinc_star_bp\tlan_star_ffp\tlan_star_bp\tap_star_ffp\tap_star_bp\trun_time\tenergy_change\tintegration_time\n')
        for i in range(len(masses_directories)):
            mass_dir = masses_directories[i]
            fname = './particles/'+mass_dir+'/stables_'+mass_dir+'.txt'
            with open(fname) as infile:
                for line in infile:
                    outfile.write(line)
            infile.close()
    outfile.close()

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
        n,b,p = plt.hist(parameter, bins = 10, color = 'c')
    elif(amount == 100):
        n,b,p = plt.hist(parameter, bins = 30, color = 'c')
    else:
        n,b,p = plt.hist(parameter, bins = 30, color = 'c')

    plt.title('Histogram of captures for each '+df_name, fontsize=20)
    plt.xlabel(df_name, fontsize=20)
    plt.ylabel('number of captures', fontsize=20)
    plt.ylim(0,max(n)*1.05)
    plt.savefig('./plots/statistics/number_captures_'+df_name+'.png')
    plt.close()

def plot_histogram_phi(parameter, df_name):

    f = plt.figure(figsize=(30,15))

    parameter = np.array(parameter)

    #n,b,p = plt.hist(parameter, bins = int(len(parameter)/3.0), color = 'c')
    n,b,p = plt.hist(parameter, bins = 10, color = 'c')

    plt.title('Histogram of captures for each '+df_name, fontsize=20)
    plt.xlabel(df_name, fontsize=20)
    plt.ylabel('number of captures', fontsize=20)
    plt.xlim(0.0,360.0)
    # plt.ylim(0,20)
    plt.savefig('./plots/statistics/'+df_name+'.png')
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
    ms_bp_par, as_bp_par, bs_ffp_par, phis_bp_par, incs_bp_par, lans_bp_par, es_star_bp_par, es_star_ffp_par, smas_star_ffp_par, smas_star_bp_par, incs_star_ffp_par, incs_star_bp_par, lans_star_ffp_par, lans_star_bp_par, aps_star_ffp_par, aps_star_bp_par = read_parameters_df()

    #Getting all the ms, as, bs and phis that were combined
    ms_bp_par_uniq = uniq_list(ms_bp_par)
    as_bp_par_uniq = uniq_list(as_bp_par)
    bs_ffp_par_uniq = uniq_list(bs_ffp_par)
    phis_bp_par_uniq = uniq_list(phis_bp_par)
    incs_bp_par_uniq = uniq_list(incs_bp_par)
    lans_bp_par_uniq = uniq_list(lans_bp_par)

    #Numbers of each parameter
    n_ms_bp_par = len(ms_bp_par_uniq)
    n_as_bp_par = len(as_bp_par_uniq)
    n_bs_ffp_par = len(bs_ffp_par_uniq)
    n_phis_bp_par = len(phis_bp_par_uniq)
    n_incs_bp_par = len(incs_bp_par_uniq)
    n_lans_bp_par = len(lans_bp_par_uniq)

    print n_ms_bp_par, n_as_bp_par, n_bs_ffp_par, n_phis_bp_par, n_incs_bp_par, n_lans_bp_par

    #Number of captures
    #m_bp
    plot_histograms(ms_bp, n_ms_bp_par, 'm_bp')
    #a_bp
    plot_histograms(as_bp, n_as_bp_par, 'a_bp')
    #b_ffp
    bs_ffp_over_a_bp = np.array([])
    for a_bp_par in as_bp_par_uniq:
        to_add = np.array(np.array(bs_ffp)[ np.where( df['a_bp']==a_bp_par) ]) / a_bp_par
        bs_ffp_over_a_bp = np.concatenate((bs_ffp_over_a_bp, to_add))
    plot_histograms(bs_ffp_over_a_bp, n_bs_ffp_par, 'b_ffp_over_a_bp')
    #phi_bp
    plot_histograms(phis_bp, n_phis_bp_par, 'phi_bp')

    for aa in as_bp_par_uniq:
        a_indices = np.where(as_bp == aa)
        plot_histogram_phi(np.array(phis_bp)[a_indices], 'cuts_a/phi_bp_a_'+str(aa))

    for mm in ms_bp_par_uniq:
        m_indices = np.where(ms_bp == mm)
        plot_histogram_phi(np.array(phis_bp)[m_indices], 'cuts_m/phi_bp_m_'+str(mm))

    bss = [-2.54438100e+00, -2.29727000e+00,-2.17912800e+00,
            -2.12390100e+00,  -2.12128300e+00,  -2.10380800e+00, -2.09698600e+00,
            -2.09346900e+00,  -1.90848000e+00,  -1.72405100e+00,  -1.52662800e+00,
            -1.37836200e+00,  -1.30747700e+00,  -1.27434100e+00,  -1.27277000e+00,
            -1.26228500e+00,  -1.25819200e+00,  -1.25608100e+00,  -6.36160000e-01,
            -5.08876000e-01,  -4.59454000e-01,  -4.35826000e-01,  -4.24780000e-01,
            -4.24257000e-01,  -4.20762000e-01,  -4.19397000e-01,  -4.18694000e-01,
            4.18694000e-01,   4.19397000e-01,   4.20762000e-01,   4.24257000e-01,
            4.24780000e-01,   4.35826000e-01,   4.59454000e-01,   5.08876000e-01,
            6.36160000e-01,   1.25608100e+00,   1.25819200e+00,   1.26228500e+00,
            1.27277000e+00,   1.27434100e+00,   1.30747700e+00,   1.37836200e+00,
            1.52662800e+00,   1.72405100e+00,   1.90848000e+00,   2.09346900e+00,
            2.09698600e+00,   2.10380800e+00,   2.12128300e+00,   2.12390100e+00,
            2.17912800e+00,   2.29727000e+00,   2.54438100e+00]

    for bb in bss:
        b_indices = np.where(bs_ffp == bb)
        b_subset = np.array(phis_bp)[b_indices]
        if (len(b_subset)!=0):
            plot_histogram_phi(b_subset, 'cuts_b/phi_bp_b_'+str(bb))

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

    #Parameters and Stables Files
    create_parameters_and_stables_files()

    #Read df
    df, ms_bp, as_bp, bs_ffp, phis_bp, incs_bp, lans_bp, es_star_bp, es_star_ffp, smas_star_ffp, smas_star_bp, incs_star_ffp, incs_star_bp, lans_star_ffp, lans_star_bp, aps_star_ffp, aps_star_bp = read_stables_df()

    #Plottts
    # create_plots_folders()
    #Make statistics plots
    make_statistical_plots(df, ms_bp, as_bp, bs_ffp, phis_bp)
