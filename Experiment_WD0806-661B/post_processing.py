import os, math
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from amuse.support import io
from amuse.units import units


def create_parameters_file():

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

def read_df():
    df = pd.read_csv('./parameters.txt', sep='\t')

    folders = df['folder']
    filenames = df['filename']
    ms_bp = df['m_bp']
    as_bp = df['a_bp']
    bs_ffp = df['b_ffp']
    phis_bp = df['phi_bp']
    energy_changes = df['energy_change']

    return df, folders, filenames, ms_bp, as_bp, bs_ffp, phis_bp, energy_changes

def find_closest_approach(b_ffp, a_bp_initial):

    return (b_ffp**2)/math.sqrt(b_ffp**2 + 1600*a_bp_initial)


def plot_trajectory(x,y,number_of_planets, fname):

    colors = ['magenta', 'green', 'DarkOrange', 'red']

    f=plt.figure(figsize=(35,15))

    x_star = x[:,0]
    x_ffp = x[:,1]

    y_star = y[:,0]
    y_ffp = y[:,1]

    plt.plot(x_star,y_star,'y',label='Star')
    plt.scatter(x_star[0],y_star[0],c='black',marker='*')
    plt.scatter(x_star[-1],y_star[-1],c='y',marker='*')

    plt.plot(x_ffp,y_ffp,'c',label='FFP')
    plt.scatter(x_ffp[0],y_ffp[0],c='black')
    plt.scatter(x_ffp[-1],y_ffp[-1],c='c')

    for i in range(0, number_of_planets):

        x_planet = x[:,i+2]
        y_planet = y[:,i+2]

        color_planet = colors[i]

        plt.plot(x_planet,y_planet,color=color_planet,label='BP',alpha=0.5)
        plt.scatter(x_planet[0],y_planet[0],c='black')
        plt.scatter(x_planet[-1],y_planet[-1],color=color_planet)

    plt.axhline(y=0, xmin=-80, xmax=10, c='black', linestyle='--')
    plt.axvline(x=0, ymin=-5, ymax=2, c='black', linestyle='--')

    plt.title('Trajectory FFP (nbody units)')
    plt.xlabel("$x$", fontsize=20)
    plt.ylabel("$y$", fontsize=20)
    plt.legend()

    plt.savefig(fname)

    plt.xlim(-3.1,1.1)
    plt.ylim(-1,1)

    plt.savefig(fname+'_zoom.png')
    plt.close()


if __name__ in ('__main__', '__plot__'):

    #Create file with all the runs made per line
    create_parameters_file()

    #Masses
    m0 = 0.58 #MSun
    m_ffp = 7.5 #MJupiter

    df, folders, filenames, ms_bp, as_bp, bs_ffp, phis_bp, energy_changes = read_df()

    num_files = len(filenames)

    for i in range(num_files):

        folder = folders[i]
        filename = filenames[i]
        m_bp = ms_bp[i]
        b_ffp = bs_ffp[i]

        mass_ratio = m_ffp/m_bp

        bodies = io.read_set_from_file('./particles/'+folder+'/'+filename, 'hdf5')
        snapshot = 1

        fig = plt.figure()

        #0 if stable 1 if not
        stability = []
        times = []
        xs = []
        ys = []

        for data in bodies.history:

            #Order: star - ffp - bp
            e_values = data.eccentricity
            a_values = data.semimajoraxis.value_in(units.AU)
            vx_values = data.vx.value_in(units.kms)
            vy_values = data.vy.value_in(units.kms)
            x_values = data.x.value_in(units.AU)
            y_values = data.y.value_in(units.AU)

            t_value = data.time.value_in(units.yr)[0]

            if (snapshot == 1):
                closest_approach = find_closest_approach(b_ffp, a_values[2])

            #inequality: closest_approach >= a_values[2]*((1+mass_ratio)**(3.0/2.0)-(1-e_values[2]**2)**(1.0/2.0))**2.0/(2.0*mass_ratio**2.0)
            right_side = a_values[2]*((1+mass_ratio)**(3.0/2.0)-(1-e_values[2]**2)**(1.0/2.0))**2.0/(2.0*mass_ratio**2.0)

            if(closest_approach >= right_side):
                stability.append(0)
            else:
                stability.append(1)

            times.append(t_value)
            xs.append(x_values)
            ys.append(y_values)

            snapshot += 1

        plt.scatter(times, stability)
        plt.ylim(-0.5,1.5)
        plt.savefig('./tests/test_stability'+str(i)+'.png')
        plt.close()


        plot_trajectory(np.array(xs),np.array(ys),1,'./tests/trajectory'+str(i)+'.png')
