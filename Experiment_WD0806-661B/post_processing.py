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

def is_hill_stable(m_ffp, m_bp, closest_approach, a_values, e_values):

    mass_ratio = m_ffp/m_bp

    #inequality: closest_approach >= a_values[2]*((1+mass_ratio)**(3.0/2.0)-(1-e_values[2]**2)**(1.0/2.0))**2.0/(2.0*mass_ratio**2.0)
    right_side = a_values[2]*((1+mass_ratio)**(3.0/2.0)-(1-e_values[2]**2)**(1.0/2.0))**2.0/(2.0*mass_ratio**2.0)

    if(closest_approach >= right_side):
        return True
    else:
        return False

#def is_petrovich_stable(m0, m_ffp, m_bp, a_values, e_values):
def is_petrovich_stable(m0, m_ffp, m_bp, a_values, e_values, snapshot):

    masses = [m0, m_ffp, m_bp]
    inner_index = np.argmin(a_values[1:]) + 1
    outer_index = np.argmax(a_values[1:]) + 1

    a_in = a_values[inner_index]
    mu_in = masses[inner_index]/masses[0]
    e_in = e_values[inner_index]

    a_out = a_values[outer_index]
    mu_out = masses[outer_index]/masses[0]
    e_out = e_values[outer_index]

    #inequality: (a_out*(1-e_out))/(a_in*(1+e_in)) > 2.4*((max(mu_in, mu_out))**(1.0/3.0))*((a_out/a_in)**(1.0/2.0)) + 1.15
    left_side = (a_out*(1-e_out))/(a_in*(1+e_in))
    right_side = 2.4*((max(mu_in, mu_out))**(1.0/3.0))*((a_out/a_in)**(1.0/2.0))+1.15

    plt.scatter(snapshot,left_side,c='r',s=3, linewidth='0')
    plt.scatter(snapshot,right_side,c='g',s=3, linewidth='0')

    if(left_side > right_side):
        return True
    else:
        return False


def plot_trajectory(x,y,number_of_planets, fname):

    colors = ['magenta', 'green', 'DarkOrange', 'red']

    f=plt.figure(figsize=(70,30))

    x_star = x[:,0]
    x_ffp = x[:,1]

    y_star = y[:,0]
    y_ffp = y[:,1]

    plt.plot(x_star,y_star,'y',label='Star')
    plt.scatter(x_star[0],y_star[0],c='black',marker='*')
    plt.scatter(x_star[-1],y_star[-1],c='y',marker='*')

    plt.plot(x_ffp,y_ffp,'c',label='FFP', lw = 2)
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

    plt.title('Trajectory FFP (nbody units)', fontsize=40)
    plt.xlabel("$x$", fontsize=40)
    plt.ylabel("$y$", fontsize=40)
    plt.legend(fontsize=40)
    plt.savefig(fname)
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

        bodies = io.read_set_from_file('./particles/'+folder+'/'+filename, 'hdf5')
        snapshot = 1

        #0 if stable 1 if not
        hill_stability = []
        petrovich_stability = []

        times = []
        xs = []
        ys = []

        fig = plt.figure()

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

            if(is_hill_stable(m_ffp, m_bp, closest_approach, a_values, e_values)):
                hill_stability.append(0)
            else:
                hill_stability.append(1)

            #if(is_petrovich_stable(m0, m_ffp, m_bp, a_values, e_values)):
            if(is_petrovich_stable(m0, m_ffp, m_bp, a_values, e_values, snapshot)):
                petrovich_stability.append(0)
            else:
                petrovich_stability.append(1)

            times.append(t_value)
            xs.append(x_values)
            ys.append(y_values)

            snapshot += 1

        plt.savefig('./tests/petrovich'+str(i)+'.png')
        plt.close()

        if(sum(hill_stability) < snapshot-1):
            fig = plt.figure()
            plt.plot(times, hill_stability)
            plt.ylim(-0.5,1.5)
            plt.savefig('./tests/hill_stability'+str(i)+'.png')
            plt.close()

        if(sum(petrovich_stability) < snapshot-1):
            fig = plt.figure()
            plt.plot(times, petrovich_stability)
            plt.ylim(-0.5,1.5)
            plt.savefig('./tests/petrovich_stability'+str(i)+'.png')
            plt.close()

        plot_trajectory(np.array(xs),np.array(ys),1,'./trajectories/trajectory'+str(i)+'.png')
