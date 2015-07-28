import os

if __name__ in ('__main__', '__plot__'):

    #Masses Directories
    masses_directories = os.listdir('./particles/')

    for m_dir in masses_directories:
        command = 'amuse.sh secular_evolution.py '+m_dir
        os.system(command)
