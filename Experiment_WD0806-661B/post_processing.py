import os, math
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from sympy.solvers import solve
from sympy import Symbol

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

def read_parameters_df():

    df = pd.read_csv('./parameters.txt', sep='\t')

    folders = df['folder']
    filenames = df['filename']
    ms_bp = df['m_bp']
    as_bp = df['a_bp']
    bs_ffp = df['b_ffp']
    phis_bp = df['phi_bp']

    return df, folders, filenames, ms_bp, as_bp, bs_ffp, phis_bp

def read_stables_df():

    df = pd.read_csv('./stables.txt', sep='\t')

    filenames = df['filename']
    ms_bp = df['m_bp']
    as_bp = df['a_bp']
    bs_ffp = df['b_ffp']
    phis_bp = df['phi_bp']

    return df, filenames, ms_bp, as_bp, bs_ffp, phis_bp

def find_closest_approach(b_ffp, a_bp_initial):

    return (b_ffp**2)/math.sqrt(b_ffp**2 + 1600*a_bp_initial)

def solve_for_x(m0, m_bp, m_ffp):
    x = Symbol('x')
    solution = solve((m_bp+m_ffp)*x**5 + (2*m_bp+3*m_ffp)*x**4 + (m_bp+3*m_ffp)*x**3 - (m_bp+3*m0)*x**2 - (3*m0+2*m_bp)*x - (m0+m_bp), x)

    return solution

def is_hill_stable(m0, m_ffp, m_bp, a_values, e_values):

    M = m0 + m_bp + m_ffp
    mu = m0 + m_bp

    a_1 = a_values[2]
    a_2 = a_values[1]

    e_1 = e_values[2]
    e_2 = e_values[1]

    x = solve_for_x(m0, m_bp, m_ffp)[0]

    f_x = m0*m_bp + (m0*m_ffp)/(1+x) + (m_bp*m_ffp)/(x)
    g_x = m_bp*m_ffp + m0*m_ffp*(1+x)**2 + m0*m_bp*(x**2)

    A = -(f_x**2)*g_x/(m_ffp**3 * mu**3 * (1-e_2**2))
    beta = (m0*m_bp/m_ffp)**(3.0/2.0)*(M/(mu**4))**(1.0/2.0)*((1-e_1**2)/(1-e_2**2))**(1.0/2.0)
    y = ((a_1*m_ffp*mu)/(a_2*m0*m_bp))**(1.0/2.0)

    equation = (1+y**2)*(beta**2*y**2 + 2*beta*y + 1) - A*y**2

    if(equation >= 0.0):
        return True
    else:
        return False

def process_stable(stables_file, folder, filename, m_bp, a_bp, b_ffp, phi_bp):

    #Add its characteristics to another text
    stables_file.write(filename+'\t'+str(m_bp)+'\t'+str(a_bp)+'\t'+str(b_ffp)+'\t'+str(phi_bp)+'\n')

    #Copy the file to other folder
    os.system('cp particles/'+folder+'/'+filename+' stables/'+filename)

def extract_stables():

    #Create file with all the runs made per line
    create_parameters_file()

    #Masses
    m0 = 0.58 #MSun
    m_ffp = 7.5 #MJupiter

    #Number of snapshots
    n_snapshots = 500

    df, folders, filenames, ms_bp, as_bp, bs_ffp, phis_bp = read_parameters_df()

    num_files = len(filenames)

    #To record the stables
    stables_file = open('./stables.txt','w')
    stables_file.write('filename\tm_bp\ta_bp\tb_ffp\tphi_bp\n')

    for i in range(num_files):

        folder = folders[i]
        filename = filenames[i]
        m_bp = ms_bp[i]
        a_bp = as_bp[i]
        b_ffp = bs_ffp[i]
        phi_bp = phis_bp[i]

        bodies = io.read_set_from_file('./particles/'+folder+'/'+filename, 'hdf5')

        #Order: star - ffp - bp
        e_values = bodies.eccentricity
        a_values = bodies.semimajoraxis.value_in(units.AU)

        if(is_hill_stable(m0, m_ffp, m_bp, a_values, e_values)):
            process_stable(stables_file, folder, filename, m_bp, a_bp, b_ffp, phi_bp)

    stables_file.close()


if __name__ in ('__main__', '__plot__'):

    extract_stables()

    df, filenames, ms_bp, as_bp, bs_ffp, phis_bp = read_stables_df()

    num_files = len(filenames)

    for i in range(num_files):

        filename = filenames[i]
        m_bp = ms_bp[i]
        a_bp = as_bp[i]
        b_ffp = bs_ffp[i]
        phi_bp = phis_bp[i]

        bodies = io.read_set_from_file('./stables/'+filename, 'hdf5')

        for data in bodies.history:

            #Order: star - ffp - bp
            e_values = data.eccentricity
            a_values = data.semimajoraxis.value_in(units.AU)
            vx_values = data.vx.value_in(units.kms)
            vy_values = data.vy.value_in(units.kms)
            x_values = data.x.value_in(units.AU)
            y_values = data.y.value_in(units.AU)

            t_value = data.time.value_in(units.yr)[0]
